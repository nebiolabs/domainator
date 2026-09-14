[index](README.md)

# The Domainator Similarity Network Viewer (DSNV)

A guide to using the interactive viewer. For the on-disk `.dsnv` bundle schema see
[File Formats](file_formats.md#domainator-similarity-network-viewer-bundles-dsnv); for the
scriptable, read-only companion see `ssn_navigator.py` in the
[program table](../README.md).

- [What it is](#what-it-is)
- [Building a viewer](#building-a-viewer)
- [Opening a viewer](#opening-a-viewer)
- [The threshold: Cluster Splits vs Threshold](#the-threshold-cluster-splits-vs-threshold)
- [The network canvas](#the-network-canvas)
- [Selecting nodes](#selecting-nodes)
- [Selection presets](#selection-presets)
- [View Settings](#view-settings)
- [Colors](#colors)
- [The Node Metadata table](#the-node-metadata-table)
- [Editing metadata](#editing-metadata)
- [Select by value](#select-by-value)
- [Per-column charts and tables](#per-column-charts-and-tables)
- [Exports](#exports)
- [Saving a session](#saving-a-session)
- [Saving an extraction](#saving-an-extraction)
- [Keyboard and mouse reference](#keyboard-and-mouse-reference)
- [Notes and limits](#notes-and-limits)

## What it is

A sequence similarity network (SSN) is usually explored by picking a similarity cutoff,
thresholding an all-vs-all matrix, and looking at the connected components. Picking the
cutoff is the hard part, and re-thresholding a large matrix for each guess is slow.

The DSNV inverts that. `build_ssn_viewer.py` computes a **maximum spanning tree** of the
similarity matrix once and records the order in which clusters merge along it. The
viewer replays that merge order in the browser, so **moving the threshold slider is
instant** and works the same on a 200-node network and a 200,000-node one — the bundle
scales with the number of nodes, not with the O(n²) edges. Everything else in the app
hangs off that: cluster outlines, the metadata table, coloring, charts, and extractions
are all computed from the tree at the threshold you are currently looking at.

The whole viewer is **one self-contained HTML file**. It has no CDN links, no server, and
no network access: it runs from `file://`, from a USB stick, or from a shared drive, and
it works offline. That also means everything you do in it — edits, colors, selections —
lives in the browser tab until you explicitly save it back out.

## Building a viewer

```bash
# a matrix from any of the comparison tools (seq_dist.py, compare_contigs.py,
# kmer_dist.py, structure_dist.py, ...)
seq_dist.py -i proteins.fasta -r proteins.fasta -o scores.hdf5 --mode score

# a bundle plus a standalone page with the data already inside it
build_ssn_viewer.py -i scores.hdf5 \
    --html network.html --embed_data \
    --metadata families.tsv taxonomy.tsv \
    --color_by family --label_by locus_tag \
    --name "Nuclease candidates"
```

| Option | Why you want it |
| --- | --- |
| `--html PATH` | Write the viewer page. Without `--embed_data` it is an empty shell that opens bundles you pick from disk. |
| `--embed_data` | Inline the bundle into the page, so opening the file shows the network immediately. This is the shareable artifact — one file, no dependencies. |
| `-o PATH` | Write the `.dsnv` bundle on its own (for `ssn_navigator.py`, or to load into an existing shell page). |
| `--metadata A.tsv B.tsv` | Tab-separated tables keyed by node id in the first column. Several are merged left-to-right. Column types are inferred as text, integer or float. |
| `--color_by`, `--label_by` | The columns the viewer selects on open. Both are changeable in the UI. |
| `--categorical COL ...` | Numeric columns to color as discrete categories rather than a gradient — cluster numbers, plasmid groups, anything where the number is a name. Toggleable per column later. |
| `--subset`, `--subset_file` | Restrict to a list of node ids before building. This is the right way to narrow a network when the subset is not MST-connected (see [Saving an extraction](#saving-an-extraction)). |
| `--max_merge_events N` | How many of the strongest merge events reach the slider and the split plot (default 500, `0` for all). A few more are added on top so that no 5% band of the threshold axis is left empty, so the real count is between N and N + 20 — the chart's caption reports it. Raise it if the threshold you care about falls between two stops. |
| `--merge_impact_metric` | What the split plot's y-axis measures: `min_child` (default; the node count of the smaller piece) or `product` (the product of the two piece sizes). |

## Opening a viewer

- A page built with `--embed_data` loads its network on open. Nothing else to do.
- A page built without it, or one where you want a different network, has a file picker at
  the top. It accepts `.dsnv` bundles, gzipped or plain JSON, and saved sessions.
- A **saved session** is just a bundle with an extra `app_state` section, so you open it
  exactly the same way and the UI state comes back with it.

If a bundle was written by a newer version of Domainator than the page, the page loads it
anyway and says so in the status line — refusing a file it can very likely still draw is
the worse failure.

## The threshold: Cluster Splits vs Threshold

The top panel is the control that matters most.

Its four counters — **Clusters**, **Split Links** (MST links still joining two drawn
clusters), **Shown Nodes** and **Hidden Nodes** — describe what is on the canvas at the
current threshold and view settings.

The chart plots, against similarity threshold:

- **Largest single split** — the biggest thing that breaks apart at that threshold,
  measured by `--merge_impact_metric`.
- **Moving sum (5% window)** — the same quantity totalled over a sliding 5% window of the
  threshold range, which reveals *regions* of instability that individual spikes hide.
- **Current threshold** — where the slider is now.

Read it as a stability map. Flat stretches are thresholds where the clustering does not
change much, and are the defensible places to cut; tall spikes are thresholds where one
decision reshapes the whole picture.

The threshold axis is ticked at round numbers, and each tick label carries the precision
its step needs, so a label always names the value the tick is drawn at — on a chart whose
axis spans a few thousandths you get `0.867 0.868 0.869`, not five ticks all reading `0.87`.

**Hover any part of the chart for a readout**, and **click to set the threshold**.
Hovering a lollipop names the merge under the pointer — how many splits of that size
happened at that threshold, how many happened there in total, and the moving sum in the
window around it; hovering out on the moving-sum line reads that line's value at the
cursor. Clicking a lollipop jumps straight to its threshold; clicking anywhere else jumps
to the nearest stop. Either way your pan and zoom are kept, so you can click along the
chart and watch one region break up.

**Scroll to zoom, drag to pan, double-click to reset** — or use the **Reset zoom**
button beside the export buttons, which lights up only while the chart is zoomed. On a
network whose merges crowd into a narrow band of scores, the whole axis draws them as one
unreadable column; zooming in is how you read the band and click a particular event in it.
Shift-scroll (or a trackpad's horizontal scroll) pans without changing the zoom. The line
beside the buttons says which stretch of the axis is showing while you are in it.

Zooming changes only what the chart draws. The slider still reaches every stop, and the
dashed **Current threshold** line simply disappears while the threshold sits outside the
window rather than being pinned to an edge, where it would claim to be somewhere it is not.
Both chart exports follow the zoom, so a figure can be of one interesting band rather than
of the whole axis. The window is saved with a session.

**The stems are a capped selection, and the caption under the chart says so** — for
example *"515 of 161,763 merge events plotted — the strongest 500 by impact, plus 15 so
that every 5% of the axis with an event to show has one."* The axis and the moving-sum
line always span *every* event, so without that caption a stretch with no stems would be
indistinguishable from a stretch where nothing happens. Raise `--max_merge_events` (or set
it to `0`) to plot more of them.

Below it, the **slider** moves between *stops*. A stop is a threshold at which the graph
actually splits — nothing between two stops produces a different clustering, so the
slider has nowhere else meaningful to go. The lowest stop sits just below the weakest MST
edge, which is the fully merged network (the true connected components); the highest stop
is `∞`, where every node is alone.

**Jump to threshold** moves the slider by value:

- Type a number and press Enter. It snaps to the nearest stop, so any number is legal.
- **←** and **→** step to the previous and next stop. Every press lands on a threshold
  the view has not just been at, even where dozens of stops crowd into one position on the
  slider track. This is the way to walk a cluster's breakup one event at a time, and it
  keeps your current pan and zoom so you can watch one region rather than being re-fitted
  on every step.
- **↑** and **↓** in the field do the same as **→** and **←**. (The field is a plain text
  box rather than a number input on purpose: a spinner's fixed step of 1.0 is meaningless
  on an axis whose splits can be 0.001 apart or 30 apart.)
- An arrow greys out at the ends of the stop list.

## The network canvas

Nodes are laid out per cluster, and how depends on **Layout algorithm**:

| Layout | What it draws |
| --- | --- |
| **Tree** | A tidy tree per network component, with the MST links between clusters drawn as edges. Shows how clusters hang off one another. |
| **Force-directed** | Clusters positioned by a force simulation, with edges. Best sense of overall topology; the slowest, and the one computed off the main thread. |
| **Grid** | Clusters laid out in rows, no edges. Fast and dense. |
| **Packed** (default) | Clusters as tightly packed circles sized by node count, no edges. The default because it reads well at almost any node count. |
| **Treemap** | Every node as a fixed-size square on one global lattice, no edges. A node's position never changes with the threshold, so you can watch clusters form and dissolve around stationary nodes. |

Layouts with no edges are not a limitation so much as a choice: past a few thousand nodes
the edges are a solid mat and the cluster *shapes* carry all the information.

## View Settings

The sidebar holds everything that changes what is drawn rather than what is selected.

**Minimum cluster size** hides clusters smaller than N, which is how you get a
thousand-singleton network to render usefully. **Minimum cluster size trims leaf clusters
only** softens it: instead of dropping every small cluster, it repeatedly prunes small
clusters that hang off the graph by a single edge, so a small cluster *bridging* two large
ones survives. **Reduce subcluster elongation** adds a PCA axis to the within-cluster dot
layout, which stops big clusters from being drawn as long thin smears.

**Collapse long paths** is a sub-option of leaf trimming, and deals with what leaf trimming
leaves behind. The small clusters that survive it are the bridges, and a run of them joined
end to end is a spindle: it can stretch across the whole canvas while saying nothing except
"these two big clusters are related, weakly". With the option on, a below-minimum cluster
with exactly two edges is contracted away, and a whole run of them becomes a **single edge
carrying the run's weakest link** — the weakest link is all a path can actually support,
which is the same minimum-along-the-path rule the threshold slider itself applies. A
below-minimum cluster with three or more edges is a branch point, not a pass-through, and
stays.

A contracted edge is drawn **dashed**, on screen and in the SVG export, because it is not a
measured similarity between the two clusters it joins — it is the weakest step of a path
that is no longer drawn. Contracted clusters count towards the hidden-node total, and the
pill above the chart says how many paths went away.

| Toggle | Effect |
| --- | --- |
| **Show node count labels** | Label each cluster with how many nodes it holds |
| **Show edge score labels** | Label the inter-cluster edges with their MST weight (edge-drawing layouts only) |
| **Render cluster bounds** | Draw the cluster outlines |
| **Render nodes** | Draw the individual node dots — turn it off for an outline-only figure of a very large network |
| **Sort clusters by size** | Lay clusters out largest-first rather than in hierarchy order |

Labels appear when they fit the mark they label: a cluster's node count once its bubble
is wide enough for the number, an edge score once the link is longer than the badge. That
is a property of the mark's size on screen rather than of the zoom, so a network of a few
enormous clusters keeps its labels however far out you are, while a thousand small ones
drop the labels that have no room — and the SVG export shows exactly what the screen did.

**Focus largest cluster**, **Focus selection** (centers on the selection, zooming out only
if it does not already fit) and **Reset view** (re-fit everything) move the camera without
changing the selection.

## Selecting nodes

The selection is the viewer's central idea: almost everything else — the metadata table,
the charts, `Set column`, `Export table TSV`, `Save extraction`, `Focus selection` — acts
on the selected nodes, and on **all** nodes when nothing is selected.

On the canvas:

| Gesture | Effect |
| --- | --- |
| Wheel | Zoom |
| Drag the background | Pan |
| Shift-drag a box | Select the clusters in the box |
| Ctrl + Shift-drag a box | Select the individual nodes in the box |
| Alt + Shift-drag | Deselect instead of select (add Ctrl for individual nodes) |
| Ctrl-click a node | Toggle that one node |
| Click a cluster bubble | Toggle every node inside it |

From the table: click rows to stage them (Shift-click for a range), then **Select nodes**
to promote the staged rows into the graph selection. **Deselect rows** clears the staging
without touching the graph selection.

By query: see [Select by value](#select-by-value).

**Clear selection** empties it; **Focus selection** centers the view on it, zooming out
only if it does not already fit.

## Selection presets

Ten numbered slots, above the canvas, that remember a selection.

| Action | Effect |
| --- | --- |
| `Shift`+`0`–`9` | Store the current selection into that slot |
| `0`–`9` | Recall it (replacing the selection) |
| Click a slot | Recall it |
| Shift-click a slot | Add it to the current selection |
| Alt+Shift-click a slot | Subtract it from the current selection |
| Hover or focus a slot | Outline its nodes in gray without changing anything |

Presets are saved with a session. They hold node *ids*, not positions, so they survive a
bundle being rebuilt.

The digit shortcuts are suppressed while you are typing in a field and while any dialog is
open, so they never fire when you meant to type a number.

## Colors

**Color by** picks the column that colors nodes; **Label by** picks the column drawn as
node labels (or `node_id`). **Customize colors…** opens the picker for whichever column
`Color by` names.

### Categorical columns

Text columns, and numeric columns marked categorical, get a **discrete palette**: one
color per distinct value.

The default is **Domainator distinct (64)** — the same palette `get_palette()` hands out,
so a column's colors in the viewer match what `build_ssn.py` and friends would have
drawn for the same values. Colors are assigned in sorted value order and cycle if the
column has more than 64 values.

The picker lets you:

- Choose another named palette (Tableau 10, Okabe-Ito, several ColorBrewer sets). The note
  under the menu says how many colors the palette has and whether it will repeat.
- Click any swatch to set one value's color by hand. The menu then reads *Custom colors*.
- Set the **No-value color** used for empty cells.
- **Reset to defaults** — back to Domainator distinct.
- **Load color table…** / **Save color table** — a two-column `value<TAB>#RRGGBB` TSV,
  with an em-dash (`—`) row for the no-value color. Saving writes the effective color of
  every distinct value, so it round-trips, and it is the way to hold one coloring fixed
  across several networks.

Because colors are assigned by sorted position, editing a cell can shift the colors of
other values — adding a value that sorts first pushes everything along by one. Load a
color table if you need the assignment pinned.

### Numeric columns

Numeric columns get a **gradient**: an ordered list of `{value, color}` stops — at least
two, at most twelve. Coloring interpolates between the bracketing pair and holds the end
colors flat outside the ends, so two stops are a plain ramp, three reproduce a
low/mid/high midpoint, and more shape the ramp arbitrarily. **Add stop** inserts one; the
knobs under the histogram drag them; **Reset values to data range** re-spreads them across
the data. The two end stops cannot be removed.

The histogram behind the knobs bins the **whole column**, not the current selection, since
the knobs ride that axis.

**Treat values as categories (discrete colors)** flips a numeric column to a discrete
palette — right for integer cluster numbers, plasmid groups and the like. Each mode keeps
its own palette, so flipping back and forth does not discard the gradient you set up.

**Export legend SVG / PNG** writes a standalone legend for whatever the picker is showing.

## The Node Metadata table

The table shows the selected nodes, or every node when nothing is selected.

- **Click a column header's name** to sort by it; click again to reverse. **Reset table
  sort** returns to bundle order. **Nulls last / Nulls first** places empty cells.
- **Search node_id and metadata** filters rows by substring across every column.
- **Rows per page** and **Previous/Next page** page through the result.
- **Drag the divider** between headers to resize a column.
- **Copy** in a column header copies that column's values, for the rows the table is
  currently showing — filter applied, all pages, in display order. Its inverse is
  `Paste column…`.
- **Export table TSV** writes the displayed rows, with a leading `SSN_cluster` column
  holding each node's cluster number *at the current threshold* — so the file records which
  threshold you were looking at, not just the annotations.
- **Double-click a cell** to edit it. `node_id` is read-only.

Three scopes recur throughout the app and are worth keeping straight:

| Scope | What it means |
| --- | --- |
| **Selected nodes** | The graph selection. |
| **Displayed rows** | Selected nodes (or all nodes), with the search filter applied, ignoring pagination. This is what `Copy`, `Export table TSV` and every chart use. |
| **Staged table rows** | Rows clicked in the table, awaiting `Select nodes`. |

## Editing metadata

Metadata edited in the viewer is real data: `Save session` writes it into the bundle's
`metadata` table, so `ssn_navigator.py` and any other reader see the annotations without
knowing anything about the viewer. This is the intended way to record conclusions —
"these 40 nodes are the ones worth ordering" — against the network you drew them from.

Five panels sit above the table, one open at a time:

**Add column** — a new Text or Number column. Or **Fill with cluster numbers**, which
numbers every node by the cluster it is in *at the current threshold*, largest cluster
first. That is the same column `build_ssn.py --lb <threshold> --cluster` writes, under the
same default name `SSN_cluster`, so the two are interchangeable downstream. Cluster
numbers are labels rather than magnitudes, so the column is marked categorical
automatically.

**Set column** — write one value into a column for the Selected nodes, the Staged table
rows, the Rows on this page, or All nodes. A blank value clears the cells. **Paste
column…** takes a column of values, one per line, applied to the rows on the current page
in display order — the inverse of a header's `Copy`, and the way to round-trip a column
through a spreadsheet.

**Rename column** and **Delete column** do what they say. The delete menu labels each
column *(added here)* or *(from bundle)*, because deleting a bundle column discards data
the viewer cannot regenerate — reload the bundle to get it back.

**Select by value** — see below.

## Select by value

A metadata query that edits the selection, for when the nodes you want are defined by
their annotations rather than by where they sit on the canvas.

Pick a column (`node_id` is offered alongside the metadata columns), a comparison, and a
value:

| Comparison | Meaning |
| --- | --- |
| **contains** | Case-insensitive substring |
| **is exactly** | Case-insensitive whole-value match; on a numeric column, numeric equality too, so `5`, `5.0` and `05` all match a stored 5 |
| **greater than**, **less than** | Ordering comparison |
| **between** | Inclusive at both ends; a second field appears for the upper bound |
| **matches regex** | Case-insensitive, unanchored JavaScript regular expression |

Then one of three actions:

| Button | Searches | Effect |
| --- | --- | --- |
| **Add to selection** | Every node in the network | Adds the matches |
| **Remove from selection** | The current selection | Drops the matches |
| **Subset selection** | The current selection | Keeps only the matches |

Only *Add* can grow the selection. That asymmetry is what lets the buttons compose:
`Add to selection` on one query, then `Subset selection` on a second, is the intersection
of the two — so `family contains PF00001`, then `score greater than 60`, then
`Subset selection` again for a third condition. The two narrowing buttons are disabled
until something is selected.

Details worth knowing:

- **Blank cells never match anything**, including ordering comparisons. An empty cell is an
  absence, not a zero.
- **Ordering uses the same comparator as the table's column sort**, so *greater than* means
  what sorting that column already showed you. Numbers compare as numbers (so `less than
  10` includes 9, which text order would not), and text compares as text.
- **Comparisons run against the stored value, not the table's rendering of it.** A cell
  displayed as `1,234` is matched by `1234`; you never have to type a thousands separator.
- A malformed regular expression says so in the panel and leaves the selection alone.
- Matching ignores the table's search filter and pagination — *Add* really does mean every
  node in the network. The panel reports how many of how many matched.

## Per-column charts and tables

Every metadata column header except `node_id` carries a small **chart glyph**. Clicking it
opens a menu of the chart and table kinds that make sense for that column; picking one
opens a dialog. A dropdown in the dialog switches kinds without reopening the menu.

| Kind | Offered for |
| --- | --- |
| **Frequency table** — every distinct value, its count, its two shares, and a colored square | any column |
| **Bar chart** | any column |
| **Pie chart** | any column |
| **Summary statistics** — rows, how many have a value, how many do not, distinct values, most common; plus min / quartiles / median / max / mean / SD for numeric columns | any column |
| **Histogram** | integer and float columns |
| **Box plot** — median, quartiles, Tukey 1.5×IQR whiskers, outliers | integer and float columns |
| **Cumulative distribution (ECDF)** | integer and float columns |

The numeric kinds are offered for any numeric column even while it is being colored as
discrete categories — that toggle is a coloring choice, not a claim that the numbers are
not numbers.

**Every chart summarizes exactly the rows the header's `Copy` button would copy**:
selected nodes (or all nodes), search filter applied, all pages. The scope and the row
count are written into the chart's own subtitle, so an exported file says what it counted.
Charts refresh live as the selection, the filter, the palette, or a cell value changes.

**The frequency table measures each value twice.** `Percent` is its share of the rows
being charted — a third of this selection is `beta`. `% of all` is the other direction: the
share of *that value's own* network-wide population these rows caught — and those betas are
half of every beta in the network. `All nodes` beside it is that population, so the second
denominator is never a mystery. With nothing selected and no filter the two columns say the
same thing, because the scope is then the whole network. Both go into the TSV, for the bar
and pie charts too.

**Colors are the column's own.** A frequency table gets a swatch column; bars and slices
are filled from the column's palette — including a palette you have never looked at,
because a column that has never been the `Color by` column still has the colors it would
get if it were. Numeric marks use the whole-column range, so a histogram of the color
column matches the nodes on the canvas. Charts with many categories roll the tail into one
gray `Other (n values)` entry rather than drawing an unreadable number of marks — a pie
shows 12 slices, a bar chart 24, and the frequency table 500 rows. The note under the
chart says when that has happened, and **the TSV export always carries every value.**

Each dialog exports **SVG**, **PNG** (at the resolution in `Export view PNG`'s scale
menu), and **TSV** — the numbers behind the picture, including the frequency table's
colors, the histogram's bin edges, the box plot's five-number summary and the ECDF's
points.

## Exports

| Button | Writes |
| --- | --- |
| **Export view PNG** + scale | The cluster canvas as a raster image, at 1×–8× screen resolution |
| **Export view SVG** | The cluster canvas as vector art, for a figure |
| **Export chart PNG / SVG** | The split-event chart itself, at its current zoom, under the split chart |
| **Export table TSV** | The displayed rows, plus a cluster number per node at the current threshold |
| A column header's **Copy** | One column's displayed values, to the clipboard |
| **Export legend SVG / PNG** | A standalone legend for the color column |
| **Save color table** | `value<TAB>#RRGGBB` for every distinct value |
| A chart dialog's four buttons | That chart as SVG, PNG or TSV |
| **Save session…** | The bundle plus the UI state |
| **Save extraction…** | A new bundle over the selected nodes |

Every one of these goes through the browser's ordinary download flow, into whatever
download directory it is configured to use — a page has no other way to write a file.

Both PNG buttons follow the one **PNG resolution** selector in View Settings, and both
re-rasterize the chart at that density rather than upscaling a screenshot, so a 4× export
is genuinely four times the detail. The SVG exports are built from the same geometry the
canvas was painted from — the same axis box, the same tick values, the same marks — so the
file is what you were looking at, with real text you can restyle in Illustrator or Inkscape.

One thing an exported split chart leaves out: the dashed **Current threshold** line and its
orange dot. Those say where your slider is sitting, which is not something the chart
measures, so they would only be noise in a figure.

## Saving a session

**Save session…** writes a `.dsnv` file containing the bundle *and* an `app_state` section:
layout and view toggles, threshold, pan and zoom, `Color by`/`Label by`, custom palettes,
table sort/filter/paging/column widths, the selection, and all ten presets. Metadata edits
go into the bundle's own metadata table, not into `app_state`.

Open it like any other bundle and the state comes back. It is version-tolerant in both
directions: a session written by a different build of the viewer still opens, because every
saved field is optional and independently skippable — anything this build does not
recognize, or cannot apply, is skipped and reported in the load status rather than aborting
the load.

Rename the network (double-click the heading, or the pencil button) and the name follows it
into the file and the browser tab.

## Saving an extraction

**Save extraction…** writes a *new, smaller bundle* containing only the selected nodes,
with its hierarchy and merge series rebuilt from the induced MST edges. Use it to peel one
interesting cluster off a large network and hand it on as its own standalone viewer.

It has one precondition: **every original network component must contribute a single
MST-connected piece to the selection.** The reason is honest rather than technical. The MST
kept one path between any two nodes and threw the rest away, so a selection that omits the
nodes along that path leaves no evidence of how the surviving pieces relate — writing them
out as separate components would assert an absence of similarity the data does not support.
Nodes in different components of the original network are exempt, since they were already
unrelated.

If your subset is not MST-connected, subset the source matrix instead:
`build_ssn_viewer.py --subset` / `--subset_file`, which *measures* those relationships
rather than inferring them.

One field cannot be rebuilt and is written empty: `graph.edges_by_threshold` counts edges
of the full graph, which a bundle never carried.

## Keyboard and mouse reference

| Key | Where | Effect |
| --- | --- | --- |
| `0`–`9` | Anywhere | Recall selection preset |
| `Shift`+`0`–`9` | Anywhere | Store the selection into a preset |
| `↑` / `↓` | Jump to threshold field | Next / previous split stop |
| `Enter` | Jump to threshold field | Snap to the nearest stop to what you typed |
| `Enter` | Select by value fields | Add to selection |
| `Enter` | Add / Set / Rename column fields | Apply |
| `Escape` | Any dialog | Close it |
| `↑` / `↓` / `Home` / `End` | Chart glyph menu | Move between kinds |
| `Enter` | Chart glyph menu | Open the highlighted kind |
| Hover | Split chart | Readout for the merge or the moving-sum line under the pointer |
| Click | Split chart | Jump to that lollipop's threshold, or to the nearest stop |
| Scroll | Split chart | Zoom the threshold axis about the pointer |
| `Shift`+scroll | Split chart | Pan the threshold axis |
| Drag | Split chart | Pan the threshold axis (does not set the threshold) |
| Double-click | Split chart | Show the whole threshold range again (its clicks do not move the threshold) |

The digit shortcuts do not fire while a text field has focus or a dialog is open. The
canvas gestures in [Selecting nodes](#selecting-nodes) are the other half of this table.

## Notes and limits

- **Nothing is saved automatically.** Edits, colors, selections and presets live in the
  tab. Use `Save session…` before closing it.
- **The matrix must be a similarity matrix, not a distance matrix.** The spanning tree
  keeps the *highest*-weight edges, and the whole app reads "higher threshold" as "more
  similar, so finer clusters". Feed it a distance matrix and every reading inverts. Use one
  of the comparison tools' score modes, or convert first.
- **Zero-weight edges are dropped** when the tree is built, so a pair that scored exactly
  zero is treated as unrelated. Nodes with no non-zero edge at all end up as singletons.
- **`--max_merge_events` bounds the slider and the split chart's stems.** With the default
  500, a very large network's slider carries the 500 strongest merges plus up to 20 more
  for axis coverage — not every one. If the threshold you want falls between two stops,
  rebuild with a higher value or `0`. The caption under the chart always says how many of
  how many are drawn.
- **A connected network has a long straggler tail.** An MST must attach every outlier, so
  the weakest MST edges are single nodes joining the giant component. Those merges have the
  smallest impact there is, so they are the first thing the cap drops — which is why the
  left end of a large network's axis can look empty. It is not: it is a region where only
  1-3 nodes move at a time.
- **Node ids must match** between the matrix labels and the first column of each
  `--metadata` TSV. Unmatched metadata rows are dropped, and nodes with no metadata row get
  empty cells.
- **Metadata column types are text, integer or float.** Types are inferred at build time;
  a column added in the viewer is Text or Number.
- **`Label by` labels are capped, not zoom-gated.** At most 250 per-node labels are drawn
  at once, and only beside dots large enough to see, so they fill in as you zoom into a
  region rather than all at once.
- Large layouts (force-directed especially) are computed off the main thread and may take a
  moment after a threshold change; the canvas says so while it works.
- In the `Treemap` layout, `Minimum cluster size` suppresses small clusters' outlines rather
  than hiding their nodes — the shown/hidden node counts above the chart say which rule is
  in force.
- **`Collapse long paths` changes an edge's meaning, never the topology.** Two clusters are
  connected after a collapse exactly when they were connected before it; what changes is
  that the edge between them now reports a path minimum rather than one MST edge, which is
  why it is dashed. It does nothing while `Minimum cluster size` is 1, when nothing is below
  the minimum, and nothing in the `Treemap` layout, which never hides a cluster. In `Grid`
  and `Packed` the contracted clusters do disappear, but those layouts draw no links, so
  nothing shows what replaced them — the option is worth having in `Tree` and
  `Force-directed`.
