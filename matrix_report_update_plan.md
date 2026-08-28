# Update plan: replace "Cluster Splits vs Threshold" in matrix_report

Target repo: `domainator` (`/home/sean/scripts/python/domainator`)

Replace the current split-events chart in generated `matrix_report.html` files with an
**overlapping-bars + 5% moving-sum** visualization.

---

## 1. Motivation

The current chart draws one stem per threshold at height `merge_impact`, which is the
**sum** of every merge occurring at that threshold. A single 88-node merge and
thirty-four merges totalling 36 nodes therefore render as comparable stems. Those are
opposite conclusions — *"one large cluster split off"* versus *"the graph is crumbling
into dust"* — and the reader cannot currently tell them apart.

The replacement shows both quantities against one x axis:

- **Bars + beads (left axis)** — the **largest single merge** at each threshold, with one
  bead per distinct merge size. Answers *"how big are the new clusters?"*
- **Moving-sum line (right axis)** — total nodes displaced within a sliding window of 5%
  of the threshold range. Answers *"how much is decomposing around here overall?"*

Read together: a tall bar with a single bead is one large clean split; short bars sitting
under a high line are wholesale fragmentation into small pieces.

## 2. Files to change

| File | Role |
|---|---|
| `src/domainator/ssn_hierarchy.py` | computes merge events — **data change** |
| `src/domainator/matrix_report.py` | emits the Plotly block — **view change** |

---

## 3. Part A — data layer (`ssn_hierarchy.py`)

### What already exists

`component_size_summary_by_threshold()` already computes, **per individual MST edge**,

```python
merge_impact = min(left_size, right_size)        # MERGE_IMPACT_MIN_CHILD (default)
merge_impact = left_size * right_size            # MERGE_IMPACT_PRODUCT
```

storing it in `summary[row_idx, COMPONENT_MERGE_IMPACT_COL]`. The per-merge detail this
plan needs is already computed and correct. It is discarded downstream.

### Where it is lost

`threshold_merge_event_rows()` walks tie groups of equal threshold and accumulates:

```python
merge_impact += float(component_summary[row_idx, COMPONENT_MERGE_IMPACT_COL])
```

**This `+=` is the entire information loss.** Keep the sum, but also keep the terms.

### The change

Collect the individual non-zero per-row values in each tie group (rows where the edge's
endpoints were already in one component contribute `0.0` and must be skipped), and emit
them as a **size → count histogram**. Keep `merge_impact` unchanged so existing consumers
keep working:

```python
"merge_impact":      float(merge_impact),      # unchanged: the SUM over the tie group
"merge_size_counts": {int(size): count, ...},  # NEW: how many merges of each size
"largest_merge":     float(max(sizes)) if sizes else 0.0,   # NEW
"merge_count":       len(sizes),               # NEW: total merges in the tie group
```

Emit `largest_merge` and `merge_count` as their own keys rather than deriving them in
JavaScript — if the histogram is ever truncated for payload size, both must have been
computed *before* truncation.

Both payload paths (`merge_event_series` and `split_events`) call
`filter_merge_event_rows(threshold_merge_event_rows(...))`, so the new keys propagate to
both with no further edits.

### Why a histogram rather than a list

Not for payload size — that argument does not survive gzip. Measured on the two reference
networks, JSON bytes before/after the report's own gzip step:

| | raw list | raw histogram | gzip list | gzip histogram |
|---|---|---|---|---|
| CBM45 (92 events, 284 merges) | 763 B | 956 B | **169 B** | **197 B (+17%)** |
| CBM45_subset (19 events, 29 merges) | 97 B | 171 B | **53 B** | **65 B (+23%)** |

The histogram is *larger* after compression, and it is ~30 bytes either way inside a 31 KB
report. Do not "optimize" this back to a list on size grounds.

The real reason is rendering. Beads at the same `(threshold, size)` are the same pixels, so
the expanded list draws marks that cannot be seen:

| metric | network | merges | distinct (threshold, size) | marks avoided |
|---|---|---|---|---|
| `min_child` (default) | CBM45 | 284 | 126 | **56%** |
| `min_child` | CBM45_subset | 29 | 22 | 24% |
| `product` | CBM45 | 284 | 228 | 20% |
| `product` | CBM45_subset | 29 | 25 | 14% |

At `t=172`, 34 markers collapse to 2. (Under `product` the gain is smaller — products vary
more, so there are fewer ties.)

This is **not a visual no-op**: rendering both ways differs by 757 pixels, up to 193/255 per
channel. Overplotted markers accumulate antialiasing, so a stack of 34 identical beads
renders as a slightly bolder, fatter dot than a single merge — an accidental encoding of
multiplicity that no legend explains. Deduplicating makes every bead weigh the same. If
multiplicity should be visible, encode it deliberately from the count (marker size, or
hover text), not through antialiasing.

Separately, and independent of structure: **emit `int`, not `float`.** `[1.0, 1.0, …]` →
`[1, 1, …]` cut CBM45's raw payload from 1331 B to 763 B (−43%).

### Invariants (assert these in a unit test)

- `sum(merge_size_counts.values()) == merge_count`
- `max(map(int, merge_size_counts)) == largest_merge`
- `sum(int(k) * v for k, v in merge_size_counts.items()) == merge_impact` — exact, by construction
- `merge_count <= summary_row_to - summary_row_from + 1`
- `largest_merge <= merge_impact`, with equality iff `merge_count == 1`

### Integrality guard

Histogram keys must be integral, or `int()` silently collides distinct sizes. Both current
metrics are safe (`min_child` = min of two ints; `product` = product of two ints). Assert
`float(size).is_integer()` and fall back to a plain list if a future metric breaks that.

---

## 4. Part B — the moving sum (and a trap)

### The trap

`filter_merge_event_rows()` caps at `DEFAULT_MAX_MERGE_EVENTS = 500`, keeping the top N rows
*ranked by `merge_impact`*. A moving sum computed in JavaScript from the **filtered** series
silently undercounts on any network with more than 500 threshold groups — and it undercounts
precisely the small events the line exists to reveal.

### The fix

Compute the moving sum **in Python from the unfiltered rows** and emit it as its own payload
key. This also avoids an O(n·m) window scan in the browser.

```python
MOVING_SUM_WINDOW_FRACTION = 0.05   # module-level constant, keep configurable
MOVING_SUM_GRID_POINTS = 800
```

Definition — centred, over the plotted x range:

```
lo, hi   = min / max threshold_value over the UNFILTERED event rows
W        = MOVING_SUM_WINDOW_FRACTION * (hi - lo)
grid     = MOVING_SUM_GRID_POINTS points evenly spaced over [lo, hi]
value(g) = sum of merge_impact for every unfiltered event with |threshold_value - g| <= W/2
```

Emit as:

```python
"merge_moving_sum": {"window": W, "x": [...], "y": [...]}
```

Guard `hi == lo` (single threshold) — emit an empty series.

Note it sums **`merge_impact`** (the per-threshold total), *not* `largest_merge`. That
contrast between the bars and the line is the entire point of the chart.

---

## 5. Part C — the chart (`matrix_report.py`)

Replace the two traces in `component_signal_plot_block` with three. Keep the existing
`chartLayout`, `chartConfig`, `SPLIT_PLOT_MARGIN_*` constants and the null-separated
polyline idiom already used for the stems.

### Trace 1 — stems

Same `mergeEventStemX` / `mergeEventStemY` construction as now, but push `largest_merge`
instead of `merge_impact`:

```js
mergeEventStemY.push(0, d.largest_merge, null);
```

`mode: 'lines'`, `line: {width: 1.6}`, `hoverinfo: 'skip'`, `legendgroup: 'splits'`,
`showlegend: false`.

### Trace 2 — beads

One marker per distinct merge size:

```js
Object.entries(d.merge_size_counts).forEach(([size, n]) => {
    beadX.push(d.threshold_value);
    beadY.push(Number(size));     // JSON object keys are strings — Number() is required
    beadN.push(n);                // for hover text, or marker.size
});
```

`mode: 'markers'`, `marker: {size: 5, line: {width: 0}}`, `legendgroup: 'splits'`,
`showlegend: true`, `name: 'Individual cluster split'`.

**`line: {width: 0}` is required.** An outlined marker looks fine in isolation, but on a
stem carrying several closely spaced beads the outlines merge into a solid band that erases
the stem behind them. Separate bead from stem by *shade*, not by outline.

### Trace 3 — moving sum

`yaxis: 'y2'`, `mode: 'lines'`, `line: {width: 2, shape: 'hv'}`.

Use `shape: 'hv'` — a moving sum over discrete events genuinely is a step function.
Smoothing it misrepresents where the events sit.

### Layout

```js
yaxis:  {title: 'Largest single split (nodes)'},
yaxis2: {title: 'Moving sum of split size (5% window)',
         overlaying: 'y', side: 'right', rangemode: 'tozero', showgrid: false},
margin: {..., r: <increase to fit the right-hand axis title>}
```

- **Change the left axis title.** It currently reads "Size of smallest new cluster", which
  describes the per-merge metric; the axis now shows the max over the tie group.
- `showgrid: false` on `y2` — two sets of gridlines over one plot area is unreadable.
- The legend already sits at `y: 1.08`; three entries still fit.

### Hover

- **Beads** — threshold, this merge size, and the count at that size, plus `merge_count` and
  `merge_impact` via `customdata`. E.g. *"33 merges of 1 node; 34 merges here, 36 nodes total"*.
- **Moving sum** — threshold, the sum, and the window width in score units.

---

## 6. Styling

Bead and stem must be **the same hue at different lightness** (they are one series); the
moving-sum line must be a **clearly different hue**. Values from the reference
implementation, checked with a colorblind-separation validator in both light and dark:

| Role | Light | Dark |
|---|---|---|
| Bead | `#1baf7a` | `#199e70` |
| Stem | `#80d2b4` | `#196349` |
| Moving sum | `#eda100` | `#c98500` |

Stem = bead blended 38% toward the page background. Substituting the report's existing
palette is fine as long as those two relationships hold.

---

## 7. Acceptance checks

Real values from `ssns/CBM45/CBM45.full.xgmml` (285 nodes, 24 573 edges) and
`ssns/CBM45_subset/CBM45_subset.full.xgmml` (30 nodes, 435 edges):

| Network | threshold | `merge_size_counts` | count | largest | impact |
|---|---|---|---|---|---|
| CBM45 | 62.00 | `{"88":1}` | 1 | 88 | 88 |
| CBM45 | 144.00 | `{"21":1,"14":1,"7":1,"6":1,"1":1}` | 5 | 21 | 49 |
| CBM45 | 172.00 | `{"3":1,"1":33}` | 34 | 3 | 36 |
| CBM45 | 174.00 | `{"1":19}` | 19 | 1 | 19 |
| CBM45_subset | 125.00 | `{"5":1,"1":1}` | 2 | 5 | 6 |
| CBM45_subset | 149.00 | `{"3":1,"2":1,"1":1}` | 3 | 3 | 6 |
| CBM45_subset | 166.00 | `{"1":5}` | 5 | 1 | 5 |

- CBM45: **48 of 92** events have `merge_count > 1`. CBM45_subset: **5 of 19**.
- Moving-sum peaks (CBM45): ≈188 at t≈58.5, ≈114 at t≈170.5, ≈96 at t≈148.0.
- Moving-sum peak (CBM45_subset): ≈13 at t≈128.1.

### The visual regression test

At `t=172` the bar must be height **3**, not 36, while the moving-sum line peaks near 114.
**If the bar is 36, the patch is still plotting the sum and the change has not taken.**

`{"3":1,"1":33}` at that threshold is the case the whole redesign exists for: one 3-node
cluster and thirty-three single nodes falling off at once — a bar of height 3 under a
moving-sum peak of ~114.

---

## 8. Provenance and caveats

- The `merge_size_counts` values above were derived by independently replaying union-find
  over **all** edges in descending score order, not over `tree.mst_edges`. They matched the
  report's `merge_impact` for **105 of 111** events across both networks; the six that
  differ do so by 1–3 nodes, from tie-breaking among equal-score edges. Every row in the
  acceptance table is taken from a matched event.
- Treat `merge_impact` as authoritative and rely on the exact
  `sum(size * count) == merge_impact` invariant, which holds by construction in the
  domainator implementation because both come from the same summary rows.
- Reference implementation: `compare_splits.py` in this directory, `--bead-style overlapping`.
  `reconstruct_merges()` is the working version of Part A and `moving_sum()` of Part B.
  Rendered output: `closeness_out/*.splits_vs_closeness_overlapping.{light,dark}.{png,svg}`.
