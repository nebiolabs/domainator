[index](README.md)

# Key algorithms

Notes on the non-obvious algorithms behind Domainator's network tools: how the maximum
spanning tree is built (in batch and in a single streaming pass), how the threshold
tables that every cut-by-score feature reads are derived from it, and how the
split-event series that drives the threshold sliders and split plots is scanned out of
it.

This is an internals document, aimed at people changing the code. For what these
features *do*, see the [Similarity Network Viewer guide](ssn_viewer.md), the
[agent workflow guide](agent_workflows.md), and the `.dsnv` section of
[File Formats](file_formats.md).

- [Why a maximum spanning tree](#why-a-maximum-spanning-tree)
- [Batch MST: Kruskal over a sorted edge table](#batch-mst-kruskal-over-a-sorted-edge-table)
- [Threshold tables](#threshold-tables)
- [Streaming MST ∪ kNN](#streaming-mst--knn)
- [Exact MST-kNN edge counts per threshold](#exact-mst-knn-edge-counts-per-threshold)
- [Split-event scanning](#split-event-scanning)
- [Replaying the merge order in the browser](#replaying-the-merge-order-in-the-browser)
- [Ports that must stay in step](#ports-that-must-stay-in-step)

## Why a maximum spanning tree

The tools in this family all start from an all-vs-all score matrix — `seq_dist.py`,
`compare_contigs.py`, `kmer_dist.py`, `structure_dist.py` — and all ask the same
question: at a given similarity cutoff, what are the connected components?

The naive answer re-thresholds the O(n²) matrix for every candidate cutoff. The
observation that makes everything else possible is that **the connected components of a
graph thresholded at any cutoff are exactly the connected components of its maximum
spanning forest thresholded at that cutoff**. Dropping an edge from a cycle never
disconnects anything, and Kruskal only ever drops the weakest edge of a cycle, so
component structure is preserved at every cutoff at once. A forest has fewer than `n`
edges, so once it is computed, every threshold question is answered against something
that scales with the node count instead of the edge count.

Domainator computes that forest once, in `MaxTree`
([ssn_edges.py:650](../src/domainator/ssn_edges.py#L650), re-exported from
`data_matrix`), and everything downstream — `build_ssn.py --mst`, the `matrix_report.py`
plots, the `.dsnv` bundle written by `build_ssn_viewer.py`, `ssn_navigator.py` — reads
that one structure.

Scores are similarities, so the interesting tree is the *maximum* spanning tree, not the
minimum. Where the code says "MST" it means maximum spanning tree throughout. It is
really a spanning *forest*: with `skip_zeros=True` (the default) zero-score pairs are not
edges at all, so a matrix with genuinely unrelated blocks yields one tree per connected
component, which is the desired answer rather than a degenerate case.

## Batch MST: Kruskal over a sorted edge table

`MaxTree.__init__` ([ssn_edges.py:708](../src/domainator/ssn_edges.py#L708)) takes either
a `DataMatrix` or an already-sorted `SortedUndirectedEdges` table. From a matrix it calls
`sorted_undirected_edges(skip_zeros=..., agg=max)`, which takes the strict lower triangle
and sorts it by descending score. `agg=max` is how asymmetry is resolved: an undirected
edge's weight is `max(M[i, j], M[j, i])`, which matters for inputs like DIAMOND bit
scores that are not exactly symmetric.

Two details of the sort are load-bearing:

- The dense fast path packs source, target and score into **one structured array** and
  sorts it in place ([data_matrix.py:928](../src/domainator/data_matrix.py#L928)), rather
  than `argsort`-ing a score array and applying the permutation. On an n² edge list the
  separate index array is the difference between fitting in memory and not.
- The sort is **stable** and the fast path is keyed on the builtin `max` object
  identity. Any other `agg` callable — including an equivalent lambda — falls back to the
  generic base-class path, which is correct but slower. That is intentional, not a bug.

The tree itself is textbook Kruskal with a union-find over `n` nodes, path-compressing on
`find` and linking roots without rank. Edges are scanned in descending weight; an edge
joining two distinct components is kept, everything else is discarded. Complexity is
dominated by the sort: **O(E log E)** time, **O(E)** memory for the edge table, where
E = n(n−1)/2 for a dense matrix.

The scan also records, per kept edge, the two counters that make the threshold tables
cheap later:

- `group_edges_above` — how many edges in the whole graph score strictly higher than this
  edge. Because the scan is descending, that is just the size of the prefix preceding the
  edge's **tie group**, which is recomputed whenever the score changes.
- `group_mst_above` — the same count restricted to MST edges.

A fourth column on `MaxTree.edges` holds the number of graph edges scanned since the
previous MST edge (the "gap"), which is what lets a caller reason about how much of the
matrix sat between two consecutive merges.

## Threshold tables

`_build_threshold_tables` ([ssn_edges.py:825](../src/domainator/ssn_edges.py#L825)) turns
those counters into the tables every consumer actually queries: `thresholds`,
`edges_above_threshold`, `mst_edges_above_threshold`, `edges_by_threshold`,
`cluster_count_by_threshold`, `cluster_count_by_edge_count`.

Three conventions are worth knowing before touching this code, because breaking any of
them silently desynchronizes the reports from what `build_ssn.py` would actually emit:

**Rows are keyed by distinct threshold, not by MST edge.** Many MST edges can share a
weight. Keying rows by edge would give one threshold several contradictory answers;
keying by distinct weight gives it exactly one. `thresholds` is therefore the deduplicated
descending weight list, and the per-edge counters are subsampled at each tie group's first
edge.

**A cut keeps scores *strictly greater* than the threshold.** That is the `--lb`
convention of `build_ssn.py`, so every row of these tables reproduces an actual
`--lb <threshold>` run. It also means a cut *at* the lowest MST weight still drops that
weight's own tie group — the weakest merge comes apart again. That is why, when every
score is positive, an extra `threshold = 0` row is appended for the complete graph (the
`--lb 0` default), and why the viewer's slider needs a floor stop (below).

**Cluster counts are derived, not counted.** `cluster_count = n_nodes - mst_edges_above`:
each MST edge kept by a cut merges exactly two components, so the component count falls by
one per kept MST edge from a start of `n_nodes` singletons. No second pass is needed.

`threshold_row_index(value)` is the lookup from a threshold value back to its row, used by
the slider stops so a stop can carry the exact edge counts for its cut.

## Streaming MST ∪ kNN

`--mst_knn k` keeps the union of the MST and an OR-symmetric k-nearest-neighbor graph:
the kNN part gives each node its local context, and the MST guarantees the result has the
same connected components as the full graph at every threshold. In batch form that is
`transform_matrix.apply_mst_knn_sparsification`, which needs the full matrix in memory.

`StreamingMstKnnAccumulator` ([ssn_edges.py:370](../src/domainator/ssn_edges.py#L370))
computes the same selection from a stream of edges fed in one at a time via `add_edge`,
never materializing the matrix. `compare_contigs.py` uses it to sparsify while it computes
chunks of the similarity matrix, discarding each chunk after feeding it in
([compare_contigs.py:290](../src/domainator/compare_contigs.py#L290)).

### The MST half is exact, via the matroid property

Spanning forests obey

```
MSF(E₁ ∪ E₂) = MSF( MSF(E₁) ∪ E₂ )
```

— an edge rejected as the weakest in a cycle stays rejected no matter what arrives later,
because any later edge can only add cycles. So the accumulator never needs more than the
running forest (≤ `n − 1` edges) plus a bounded buffer of unprocessed edges. When the
buffer fills (`mst_buffer_cap`, default `max(4n, 100_000)`), it recomputes
`MaxTree(running_forest ∪ buffer)` with the ordinary batch Kruskal and clears the buffer.

The result is **exact**: the final forest has the same weight as a single batch `MaxTree`
over all the edges. Ties among equal-weight edges may resolve to a different but equally
valid forest, which is the only way the two paths can differ.

### The kNN half is bounded, and exact in practice

Each node keeps a dict `neighbor -> symmetric_score`, where the stored score is the
running `max(M[i, j], M[j, i])` — both directions of a pair update the same slot, so an
entry can never be left holding a stale one-directional maximum. This matches how
`build_symmetric_neighbor_rankings` ranks neighbors in the batch path.

The dicts are the only unbounded state, so they are trimmed: when a node's dict exceeds
`2 * knn_soft_cap` entries it is cut back to the `knn_soft_cap` highest-scoring neighbors.
The cap is never allowed below `k` (default `max(4k, k + 16)`), so a true top-k neighbor
cannot be evicted by the trim itself. For symmetric input this is exact. For
slightly-asymmetric input it is exact as long as no node ever holds more than
`knn_soft_cap` neighbors above its true top-k — i.e. always, short of a pathological
high-degree hub.

Two shortcuts keep the common cases cheap: at `k = 0` (`--mst_knn 0`, MST only) the
adjacency is never allocated, since it would never be read back; with `include_mst=False`
(a pure `--knn k` graph, which does **not** preserve components) the MST buffer is never
allocated.

`finalize()` merges the two halves into one `{(min, max): score}` dict, taking the MST
edge's recorded weight as authoritative — adjacency trimming may since have dropped that
pair from both neighbor dicts, but the forest still needs the edge. `to_csr()` writes each
kept edge to both cells as the symmetric score it was selected on, so the streaming and
batch paths agree on values as well as on membership.

## Exact MST-kNN edge counts per threshold

`matrix_report.py` shows how many edges a `--lb t --mst_knn k` run would emit, for every
threshold in the table and every k in `[2, max_k]`. Computing that by re-running the
selection per (threshold, k) pair would be quadratic in the table size.

`mst_knn_edge_counts_by_threshold` ([ssn_edges.py:555](../src/domainator/ssn_edges.py#L555))
instead does one descending sweep. Every (node, rank) pair from the neighbor rankings
becomes an **activation event** carrying `(score, rank, source, target)`, packed into a
single structured array and sorted in place by score. Sweeping thresholds from high to low:

- Activations that pass the cut are consumed. Each undirected pair remembers the **minimum
  rank** at which either endpoint selected it, since a pair present at rank r is present
  for every k ≥ r. A pair whose minimum rank improves moves between rank buckets.
- MST edges the cut keeps are absorbed into a prefix set. An absorbed edge is removed from
  the rank buckets so the union is not double-counted.
- The answer for each k is then `|MST prefix| + cumsum(rank buckets)[k]` — one cumulative
  sum over a `max_k`-length array per threshold.

The invariant the loop depends on is that an MST edge is always *activated* before it is
*absorbed*, which holds because both loops advance under the same descending threshold.

## Split-event scanning

Everything above answers "what does this cut look like?". The split-event series answers
the more useful question: **which cuts are worth looking at?**

Read the MST bottom-up (adding edges in descending weight) and it is a sequence of merges.
Read the threshold slider top-down (raising the cutoff) and the same sequence runs
backwards as splits. The code computes merges and the UI presents splits; the two names
refer to the same events, which is why `merge_event_*` functions produce a payload field
called `split_events`.

### Pass 1: replay the merges

`component_size_summary_by_threshold`
([ssn_hierarchy.py:95](../src/domainator/ssn_hierarchy.py#L95)) replays the MST edges in
order through a union-find, emitting one summary row per edge (plus a row 0 at
`threshold = ∞`, all singletons). Each row records the threshold, the largest component
size, the mean size over non-singleton components, this edge's merge impact, and the
deltas from the previous row.

Merge impact is one of two metrics:

- `min_child` (default) — `min(left_size, right_size)`, the node count of the smaller
  piece. This is what "how many nodes get displaced by this split" means.
- `product` — `left_size * right_size`. Not a node count, which is why
  `merge_impact_axis_labels` exists: both the Plotly chart in `matrix_report` and the
  canvas chart in the viewer read their axis titles from it so neither can label a product
  as a count.

Two running statistics need care. The mean over non-singletons is maintained as a
running `(count, sum)` pair, adjusting for the two components leaving and one arriving.
The largest component is also a scalar running value. A merge only replaces two
components with their larger combined component, so component sizes never decrease and
`largest = max(largest, merged_size)` is exact. This avoids rescanning the components or
maintaining a heap. The Python and browser implementations both use this invariant.

### Pass 2: group by threshold

`threshold_merge_event_rows` ([ssn_hierarchy.py:191](../src/domainator/ssn_hierarchy.py#L191))
collapses the per-edge rows into one event row per **distinct threshold**, for the same
reason the threshold tables are keyed that way: a cut cannot separate edges that share a
weight.

Each row carries the summed impact over the tie group *and* the terms behind it —
`merge_size_counts` (a histogram of individual impacts), `largest_merge`, `merge_count`.
The sum alone cannot distinguish one large cluster splitting off from a swarm of tiny
ones, which is exactly the distinction a user is looking for. An edge whose endpoints were
already in one component has zero impact and is not counted as a merge.

The row's `edge_index` is `first_summary_row_idx - 2`: summary row 0 is the pre-merge
state and summary row i corresponds to `mst_edges[i-1]`, so that expression names the last
MST edge scoring *strictly above* this threshold — the edge index that reproduces the
`--lb threshold_to` cut exactly.

### Capping without blanking the axis

A large network has far more threshold groups than are worth plotting, so
`filter_merge_event_rows` ([ssn_hierarchy.py:341](../src/domainator/ssn_hierarchy.py#L341))
caps the series at `max_merge_events` (default 500) ranked by impact.

The cap is applied **where the series is displayed**, never where a file is written: the
viewer's "Split events" box, `ssn_navigator.py --max_merge_events`, `matrix_report.py
--max_merge_events`. A `.dsnv` bundle stores no series at all (see
[File Formats](file_formats.md#what-is-derived-not-stored)), so raising the cap re-runs
the selection rather than requiring a rebuild from the source matrix.

### Spending the cap on a window

`filter_merge_event_rows` takes an optional `window` — a `(low, high)` threshold range —
which narrows **the top-N pool only**. The three parts of the result are then:

1. the strongest events within the window, up to `max_merge_events`;
2. the strongest event in each 5% band of the **whole** range that nothing above already
   covers — computed over every row, not just the window, so a band outside it still gets
   its stop;
3. the row at `pinned_threshold`, if the first two missed it.

`window=None` collapses (1) to a whole-range ranking, which is exactly the behaviour this
function had before windows existed, and is what the two consumers without a viewport
(`matrix_report`, `ssn_navigator`) pass.

The viewer passes its chart's zoom window, and this is what makes a *small* cap resolve
more than a large one. Zoomed out, 50 events is the strongest 50 overall. Zoom into a
band and the same 50 slots are re-spent on that band, so merges far too small to rank
globally become visible — and keep becoming visible as the zoom goes further in. Hence
`DEFAULT_WINDOWED_MAX_MERGE_EVENTS = 50` for the viewer against `DEFAULT_MAX_MERGE_EVENTS
= 500` for the two tools that get only one shot at the whole range.

Part (3) is a correctness requirement, not a nicety. The slider's stops are the selected
events; if a zoom could drop the stop the thumb is sitting on, the thumb would resolve to
a neighbouring stop and **the displayed clustering would change as a side effect of
zooming**. The viewer captures `selectedThresholdValue()` before every re-selection and
pins it, so the cut in effect survives any window change.

Only the selection is redone on zoom. The union-find replay that produces the event rows
is window-independent and runs once per bundle — `deriveMergeSeries` in the browser, the
`ssn_bundle` functions in Python — while `selectMergeEvents` (a sort and a scan) is what
each zoom and pan re-runs.

Ranking by impact alone produces a broken chart. On a connected MST-kNN graph the weak
tail of the MST is individual outliers being attached to the giant component one or two
nodes at a time — the smallest impacts there are, so the cap drops every one of them —
while the axis and the slider still span the full threshold range, because the moving sum
and the floor stop are derived from the *unfiltered* data. The left half of the plot comes
out blank and the left half of the slider has no stops.

So after the top-N pass the axis is cut into `MERGE_EVENT_DENSITY_BINS` (20) equal bands
and the strongest event in each otherwise-empty band is added back. Because the ranked
list is strongest-first, the first row seen in an empty band is the best available
representative. Back-filled rows are additions, never replacements, so the output is
bounded at `max_merge_events + density_bins` and the top-N rows are untouched. This is why
the chart caption reports a count somewhere between N and N + 20.

### The moving sum

The per-threshold stems show the largest single split; the moving sum answers the
complementary question, "how much of the graph is coming apart around here in total?".

`merge_event_moving_sum` ([ssn_hierarchy.py:251](../src/domainator/ssn_hierarchy.py#L251))
evaluates a centred window of 5% of the threshold range at 800 grid points. A nested loop
would be O(n · grid); instead the impacts are sorted by threshold once, a prefix sum is
built, and each grid point's window becomes two `searchsorted` calls and a subtraction —
O(n log n + grid). `side="left"`/`side="right"` make the window inclusive at both ends,
matching `|t − g| ≤ W/2`.

Call it on the **unfiltered** rows. Filtering keeps the highest-impact rows, so a moving
sum taken afterwards undercounts exactly the small events this series exists to reveal.

### Slider stops and the floor

`threshold_slider_stops` ([ssn_hierarchy.py:434](../src/domainator/ssn_hierarchy.py#L434))
emits `∞`, one stop per event row, and then a **floor stop**. The floor is what makes the
fully merged network reachable at all: under the strictly-above `--lb` convention every
other stop excludes its own tie group, so the lowest event stop still splits the weakest
merge back apart.

The floor sits 1% of the MST weight range *below* the weakest merge rather than at 0
(`floor_threshold_value`, [ssn_hierarchy.py:413](../src/domainator/ssn_hierarchy.py#L413)).
A network whose scores run 350–650 would otherwise spend more than half its slider track
on empty space below the data, and negative scores are not cleared by 0 at all. It is
derived from the tree rather than from the event rows, which the cap may have thinned.

Stops optionally carry a `threshold_index` into the threshold tables, via the
`threshold_index_lookup` argument. Only `matrix_report` passes it, because only that
report carries those tables — the floor stop then reads the `threshold = 0` row, which is
the complete-graph cut it stands for and is exact whenever no graph edge scores at or
below the floor. A `.dsnv` bundle has no threshold tables, so its stops omit the field
rather than index something that is not there.

## Replaying the merge order in the browser

Both HTML front-ends need to answer "what are the clusters at this stop?" fast enough to
feel instant while a slider is dragged. They do it two different ways.

**`matrix_report.py` checkpoints a union-find.** Replaying MST edges from scratch per
stop is O(stops · edges). Instead `buildClusterCheckpoints` walks the stops once and
snapshots the union-find's `parent`/`size` arrays every `CLUSTER_CHECKPOINT_STRIDE` (50)
stops ([matrix_report.py:913](../src/domainator/matrix_report.py#L913)). A query copies
the nearest preceding checkpoint and replays at most 50 stops' worth of edges from there,
trading `stops / stride` snapshots of memory for bounded per-query work.

**The viewer cuts a precomputed dendrogram.** `build_mst_component_hierarchy`
([ssn_hierarchy.py:491](../src/domainator/ssn_hierarchy.py#L491)) turns the merge order
into a binary tree: `n` leaves plus one internal node per merge, each recording the
threshold at which it formed, its size, and — via a post-order pass that assigns
`leaf_start` / `leaf_count` into a shared `leaf_order` array — a contiguous slice naming
its members. Children are ordered by minimum leaf index so the layout is deterministic.

Cutting is then a descent from the roots: recurse into a cluster whose `threshold ≤ cut`,
otherwise keep it whole (`activeClustersAtThreshold`,
[ssn_viewer_html.py:7116](../src/domainator/ssn_viewer_html.py#L7116)). The `≤` is the
strictly-above `--lb` convention again — a component that merged exactly *at* the
threshold is split back apart. Membership needs no traversal at all, just the
`leaf_order` slice. This is what lets the `.dsnv` bundle scale with the node count rather
than with O(n²) edges. [A worked example on a small
forest](file_formats.md#hierarchy) walks all of that through three cuts.

The descent costs Θ(clusters returned), not Θ(nodes): each descent turns one frontier
entry into two, so `d` descents from `r` roots leave a frontier of `r + d`, and since the
frontier *is* the result, the traversal visits exactly `2k − r` nodes for `k` clusters. It
never walks below the cut, so depth does not enter. That ranges from Θ(components) at the
floor cut to Θ(n) at the `∞` cut — against a union-find replay, which costs Θ(edges above
the cut) however few clusters come out. That difference is why `matrix_report` needs its
checkpoint stride and the viewer does not.

Asking about **one node** rather than the whole partition is cheaper still: climb from
that node's leaf while the parent merge scores strictly above the cut
(`ssn_navigator._cluster_for_node`). Ancestors merged later and so score lower, so the
condition fails once and the walk stops — Θ(depth), typically ~20 on a 20,000-node
network. `--mode node` asks exactly this at every cut-point, so it walks rather than
partitioning: on that network the per-cut-point work is 0.6 ms instead of 296 ms.

## Ports that must stay in step

Most of this section's algorithms exist twice: once in Python, once in JavaScript embedded
in `ssn_viewer_html.py`. The JS side began as an extraction-only port and is now the
viewer's primary path — since the `.dsnv` bundle stores nothing derived, the JS *is* what
produces the split chart and the slider every time a bundle is opened. The Python side
serves `matrix_report`, `ssn_navigator` and `ssn_bundle`. The pairs are:

| Python | JavaScript |
| --- | --- |
| `build_mst_component_hierarchy` | `buildExtractionHierarchy` |
| `component_size_summary_by_threshold` + `threshold_merge_event_rows` | `mergeEventRows` |
| `filter_merge_event_rows` (window + pin) | `filterMergeEventRows` |
| `merge_event_rank_key` | `compareMergeEventRank` |
| `merge_event_density_bin` | `mergeEventDensityBin` |
| `merge_event_moving_sum` | `mergeEventMovingSum` |
| `threshold_slider_stops` | `buildSliderStops` |
| `floor_threshold_value` | `floorThresholdValue` |
| `format_threshold_value` | `formatThresholdValue` |
| `clusters_at_threshold` (`ssn_bundle`) | `activeClustersAtThreshold` |
| `merge_impact_axis_labels` | read from the payload, not reimplemented |

`ssn_bundle.merge_event_rows` / `merge_event_series` / `moving_sum` / `slider_stops` are
the Python entry points that compose the first group; `deriveMergeSeries` is the
JavaScript one. The constants they share (`MOVING_SUM_WINDOW_FRACTION`,
`MOVING_SUM_GRID_POINTS`, `DEFAULT_MAX_MERGE_EVENTS`, `MERGE_EVENT_DENSITY_BINS`) are
duplicated near the top of the viewer's script.

If you change one of these algorithms, change both sides. The parity is asserted directly:
`_assert_js_series_matches_python` in `test/test_ssn_viewer_browser.py` runs the Python
implementation over the same bundle the page has open and compares the series, the stops
and the moving sum field by field. `test/test_ssn_hierarchy.py` covers the Python side on
its own.
