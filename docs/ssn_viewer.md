[index](README.md)

# Networks in the Domainator Viewer

A guide to using the interactive viewer with the similarity networks `domainator build_ssn`
writes. The Domainator Viewer also opens the trees `domainator build_tree` writes and the
projections `domainator build_projection` writes; this guide covers networks, and the table,
colors, presets and sessions work the same way for the other two. For the on-disk `.dvb`
bundle schema see [File Formats](file_formats.md#domainator-viewer-bundles-dvb); for the
scriptable, read-only companion see `domainator ssn_navigator` in the
[program table](https://github.com/nebiolabs/domainator#program-descriptions).

- [What it is](#what-it-is)
- [Building a viewer](#building-a-viewer)
- [Opening a viewer](#opening-a-viewer)
- [The threshold: Cluster Splits vs Threshold](#the-threshold-cluster-splits-vs-threshold)
- [The network canvas](#the-network-canvas)
- [Selecting nodes](#selecting-nodes)
- [Selection presets](#selection-presets)
- [View Settings](#view-settings)
- [Rooting trees](#rooting-trees)
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

The viewer inverts that. `domainator build_ssn --html_out` computes a **maximum spanning tree** of the
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
# a matrix from any of the comparison tools (seq_dist.py, domainator compare_contigs,
# domainator kmer_dist, domainator structure_dist, ...)
domainator seq_dist -i proteins.fasta -r proteins.fasta --hdf_out scores.hdf5 --mode score

# a standalone page with the network already inside it
domainator build_ssn -i scores.hdf5 \
    --html_out network.html \
    --metadata families.tsv taxonomy.tsv \
    --color_by family --label_by locus_tag \
    --name "Nuclease candidates"
```

| Option | Why you want it |
| --- | --- |
| `--html_out PATH` | Write the viewer page, with the network inside it, so opening the file shows the network immediately. This is the shareable artifact — one file, no dependencies. |
| `--dvb_out PATH` | Write the `.dvb` bundle on its own (for `domainator ssn_navigator`, to open in any viewer page, or to turn into a page later with `domainator viewer -i`). |
| `--threshold T` | Leave out edges scoring T or less, as for `--xgmml_out`, so the slider reaches down only to the weakest edge above T. The default, 0, keeps every edge. `--mst` and `--knn` change only the xgmml: the viewer always shows the network's maximum spanning tree. |
| `--metadata A.tsv B.tsv` | Tab-separated tables keyed by node id in the first column. Several are merged left-to-right. Column types are inferred as text, integer or float. |
| `--color_by`, `--label_by` | The columns the viewer selects on open. Both are changeable in the UI. |
| `--color_table PATH` | A two-column TSV (value, hex color; no header), the same file `build_ssn --color_table` reads, giving the `--color_by` column its starting colors. Values not in the table are gray, and a numeric column is colored as categories. Recorded in the bundle as `defaults.color_table`; "Reset to defaults" in the color picker returns to it. |
| `--categorical COL ...` | Numeric columns to color as discrete categories rather than a gradient — cluster numbers, plasmid groups, anything where the number is a name. Toggleable per column later. |
| `--subset`, `--subset_file` | Restrict to a list of node ids before building. This is the right way to narrow a network when the subset is not MST-connected (see [Saving an extraction](#saving-an-extraction)). |
| `--merge_impact_metric` | What the split plot's y-axis measures: `min_child` (default; the node count of the smaller piece) or `product` (the product of the two piece sizes). Recorded as the bundle's default. |
| `--name` | The network's name, shown in the heading and the browser tab. Defaults to the output's file name. |

`build_ssn` writes `--xgmml_out` for Cytoscape alongside these, from the same options, and
`domainator viewer --html_out viewer.html` writes an empty page that opens the bundles you
choose in it.

There is no `--max_merge_events`: how many split events reach the slider and the plot is
a **viewer control**, not something the bundle is built with. See
[the split chart](#the-threshold-cluster-splits-vs-threshold) below.

## Opening a viewer

- A page from `build_ssn --html_out` (or `domainator viewer -i`) loads its network on open.
  Nothing else to do.
- An empty page from `domainator viewer`, or one where you want a different network, has a
  file picker at the top. It accepts `.dvb` bundles, gzipped or plain JSON, saved sessions,
  and the `.dsnv` network bundles of earlier Domainator versions.
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
- **Moving sum (5% window)** — the same quantity totaled over a sliding 5% window of the
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

Zooming does not change the clustering, only what the chart draws — and how the slider
below it is laid out (see below). The dashed **Current threshold** line simply disappears
while the threshold sits outside the window rather than being pinned to an edge, where it
would claim to be somewhere it is not. Both chart exports follow the zoom, so a figure can
be of one interesting band rather than of the whole axis. The window is saved with a
session.

**The stems are a capped selection, and the caption under the chart says so** — for
example *"515 of 161,763 merge events plotted — the strongest 500 by impact, plus 15 so
that every 5% of the axis with an event to show has one."* The axis and the moving-sum
line always span *every* event, so without that caption a stretch with no stems would be
indistinguishable from a stretch where nothing happens.

**The `Split events` box beside that caption sets the cap** (default 50; `0` plots every
event). Changing it re-selects from the merge order the bundle already carries — no
rebuild, no source matrix, and the threshold you are looking at is kept. The cap is saved
with a session.

**The cap is spent on the part of the axis you are looking at, so zooming in shows more.**
The events plotted at any moment are:

- the strongest event in each 5% band of the **whole** range, so no stretch of the axis or
  of the slider is ever left without a stop;
- up to `Split events` of the strongest events **inside the current window**;
- the cut currently selected, which is pinned so that zooming can never move it.

Zoomed out, that second group is the strongest 50 overall — the same thing a fixed cap
would give you. Zoom into a band and those 50 slots are re-spent on the events in that
band, so merges too small to make the global list appear, and keep appearing as you go
further in. This is why the default is 50 rather than something larger: on a zoomable
chart a small number resolves *more*, because you can always go and get the detail where
you want it. `domainator matrix_report`'s split chart has no window to spend a cap on, so it keeps
a default of 500.

Below it, the **slider** moves between *stops*. A stop is a threshold at which the graph
actually splits — nothing between two stops produces a different clustering, so the
slider has nowhere else meaningful to go. The lowest stop sits just below the weakest MST
edge, which is the fully merged network (the true connected components); the highest stop
is `∞`, where every node is alone.

**The track is warped, not scaled.** A slider has about a thousand positions to offer and a
network can have five hundred stops, most of them bunched into a narrow band of scores — so
scaling stops onto the track by value alone would land a dozen of them on one position and
leave only the first of each dozen selectable. Instead every stop is given a position of its
own, and the track left over after that is spent saying where the stops are. Two consequences
worth knowing:

- **Every stop can be reached by dragging**, up to about 900 of them (raise `Split events`
  past that and the crowding returns — nothing can put a thousand stops on a thousand
  positions and still leave them apart on screen).
- **Zooming the chart magnifies the matching stretch of the track.** The zoomed band carries
  more track per unit of threshold than the rest, so the stops you are looking at spread out
  under the thumb. The threshold does not move when this happens; only the map under it does,
  so the thumb will visibly shift when you zoom.

Where stops are spread out to begin with, none of this is visible: the track is the plain
proportional scale it looks like.

**Jump to threshold** moves the slider by value:

- Type a number and press Enter. It snaps to the nearest stop, so any number is legal.
- **←** and **→** step to the previous and next stop. Every press lands on a threshold the
  view has not just been at. This is the way to walk a cluster's breakup one event at a
  time, and it keeps your current pan and zoom so you can watch one region rather than being
  re-fitted on every step.
- **↑** and **↓** in the field do the same as **→** and **←**. (The field is a plain text
  box rather than a number input on purpose: a spinner's fixed step of 1.0 is meaningless
  on an axis whose splits can be 0.001 apart or 30 apart.)
- An arrow grays out at the ends of the stop list.

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

### Nodes and Labels

**Nodes** configures the marks for individual sequences. **Color by** maps a metadata
property to the existing categorical palette or numeric gradient; **Shape by** maps a
property to categorical shapes. Shape values are always categories, even when numeric:
mapping shape to a score does not change its color gradient.

**Customize shapes...** shows each category's glyph and lets you change its assignment,
including the **No value** shape. Automatic assignments use sorted category order and
cycle through circle, square, diamond, triangle up, triangle down, and hexagon. The
dialog reports when shapes repeat; they are not unique identifiers for large category
sets. **Reset to automatic** removes overrides. With no shape mapping, dots stay circles
and treemap cells stay squares; with a mapping, the glyph fits inside the existing cell.

**Labels** chooses the **Text by** property: a metadata column, `node_id`, or None.
Label color, font, and other styling are not configurable here. Networks and projections
have one node configuration and one text-label configuration, retaining their existing
density-based label limits.

Trees additionally have **Add node**, **Add label**, **Remove**, **Move up**, and
**Move down** controls. Each leaf gets its configured markers side by side starting at
the actual branch tip, followed by its text labels. These are leaf annotations, not
metadata on internal nodes. Either list may be empty: removing annotations does not
remove branches, support values, or selection targets. **Show leaf labels** temporarily
hides the text without deleting its configuration.

The **first node configuration** supplies tree branch colors: a branch inherits a color
only when all leaves below it share that color. Reordering the first configuration
therefore changes branch colors. With no node configurations, branches are neutral.
Collapsed clades keep their structural wedge and, when text labels are enabled, their
clade-name/leaf-count summary rather than invented per-leaf annotations.

Dense trees show a compact strip for each mapped color when individual annotations are
too small to read. Long text may be shortened in the viewport; **Entire content**
image exports include all annotation text. Colors and shape assignments are shared per
property, so two node configurations using the same property use the same assignments.
Metadata edits update all affected configurations; column renames retain mappings,
and deletion clears the affected mappings with a status message.

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
that is no longer drawn. Contracted clusters count toward the hidden-node total, and the
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

**Zoom** shows the current zoom as a number (screen pixels per layout unit) and sets it when
you type one and press Enter, keeping the center of the view in place. Copy it from one
viewer into another to draw two networks, or two layouts of one network, at the same scale.

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

**Color by**, under **Nodes**, picks the column that colors a node configuration.
**Text by**, under **Labels**, picks its label text (or `node_id`).
**Customize colors...** opens the picker for that node configuration's color property,
including when it is a secondary tree marker. Palette edits, color-table import/export,
and legend exports all use that property, not necessarily the first marker's property.

### Categorical columns

Text columns, and numeric columns marked categorical, get a **discrete palette**: one
color per distinct value.

The default is **Domainator distinct (64)** — the same palette `get_palette()` hands out,
so a column's colors in the viewer match what `domainator build_ssn` and friends would have
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
`metadata` table, so `domainator ssn_navigator` and any other reader see the annotations without
knowing anything about the viewer. This is the intended way to record conclusions —
"these 40 nodes are the ones worth ordering" — against the network you drew them from.

Five panels sit above the table, one open at a time:

**Add column** — a new Text or Number column. Or **Fill with cluster numbers**, which
numbers every node by the cluster it is in *at the current threshold*, largest cluster
first. That is the same column `domainator build_ssn --threshold <threshold> --cluster` writes, under the
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

## Rooting trees

Tree mode has three rooting buttons below the tree:

- **Root at midpoint** places the root halfway along the longest leaf-to-leaf
  path, using the stored branch lengths even when the display is a cladogram.
  It requires at least two leaves, finite nonnegative lengths, and a positive
  diameter.
- **Root on selection** uses the selected group as the **outgroup**, placing the
  root halfway along its boundary branch. Select exactly one complete clade
  under the current root, or one leaf. Empty selections, the whole tree, and
  partial or nonmonophyletic groups cannot be used. Collapsing a clade does not
  change which leaves it contains.
- **Restore original root** restores the original topology, branch lengths,
  internal labels, and child order. Unlike **Reset view**, it changes the root
  rather than the zoom and scroll.

Rerooting preserves leaf-to-leaf distances, leaf metadata, colors, selections,
presets, and the alignment's sequence mapping. Branch support follows its
original split; internal labels stay with their vertices. Collapsed groups remain
collapsed if they are still the same clades, and otherwise expand. The status
line reports those expansions. Unary root stems are omitted when rerooting so
they cannot become extra leaves; root-only annotations remain recoverable in the
original-root snapshot rather than being reassigned.

The root-relative **Cut tree at** slider resets to its left endpoint (0, with the
cut line hidden), and vertical scroll returns to the top. Row height, alignment column width, and other display settings
are kept. Previously filled metadata columns are not rewritten; fill cluster
numbers again if you want a column for a new cut.

Newick export and **Save session…** use the current rooted tree. Sessions also
retain the original root, so **Restore original root** still works after reopening.

**Swap children** exchanges the top and bottom **child subtrees** of the selected
internal clade, leaving any middle children and the order inside each subtree
unchanged. Select exactly one complete internal clade; selecting the entire tree
swaps its root's children. Empty, partial, nonmonophyletic, and single-leaf
selections cannot be used. The same action again swaps the subtrees back.

Manual child order overrides ladderization at that node only; **Ladderize** remains
enabled and still applies inside its subtrees. The alignment rows move with their
leaves. Swaps do not change branch lengths, annotations, cluster assignments, or
the selection, and are included in Newick and image exports and saved sessions.
Rerooting retains manual order where the same clade and child subtrees survive;
orders that no longer apply are cleared and reported.

The **Cut tree at** slider shows its dashed marker only at interior positions,
both on screen and in PNG/SVG exports. Either endpoint hides the marker. The
**left endpoint (0) is the default** and puts all leaves in one cluster; the
right endpoint means **no cut**, with one cluster per leaf. Hiding the marker
does not change these clustering rules. Sessions preserve the slider position;
files without a saved position open at the left endpoint.

## Exports

| Button | Writes |
| --- | --- |
| **Export view PNG** + scale | The active network, tree/alignment, or projection as a raster image, at 1×–8× resolution |
| **Export view SVG** | The active view as vector art, for a figure |
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

**Image extent** controls both view image exports:

- **Current viewport** (the default) exports the portion currently visible.
- **Entire content** exports the complete current scene, including off-screen
  content, at the current zoom/detail level. The image grows instead of shrinking
  the scene to fit the viewport. In tree mode this includes every displayed row
  and every alignment column, even after scrolling horizontally or vertically.
  In projection mode it uses the currently chosen axis pair.

Entire-content export keeps collapsed clades collapsed and respects network
thresholds, filters, and display toggles. Existing label-density rules still
apply; it does not force every node to have a label. Selections are included,
but transient hover previews and drag boxes are not. Exporting does not move,
zoom, scroll, or resize the live view. Entire-content filenames end in
`_entire.svg` or `_entire.png` (with the network PNG's resolution suffix when
applicable); viewport filenames are unchanged. The extent choice is saved in
sessions, independently of chart and legend exports.

Large scenes can produce very large images. PNG export refuses dimensions over
32,767 pixels on either side or an area over 64 million pixels before allocating
the bitmap. Browsers with smaller limits can also refuse encoding; the status
line explains the failure rather than cropping or silently reducing resolution.
Choose a lower PNG resolution, zoom out, or use SVG, which has no bitmap limit.

Both PNG buttons follow the one **PNG resolution** selector in View Settings, and both
re-rasterize the chart at that density rather than upscaling a screenshot, so a 4× export
is genuinely four times the detail. The SVG exports are built from the same geometry the
canvas was painted from — the same axis box, the same tick values, the same marks — so the
file is what you were looking at, with real text you can restyle in Illustrator or Inkscape.

One thing an exported split chart leaves out: the dashed **Current threshold** line and its
orange dot. Those say where your slider is sitting, which is not something the chart
measures, so they would only be noise in a figure.

## Saving a session

**Save session…** writes a `.dvb` file containing the bundle *and* an `app_state` section:
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
`domainator build_ssn --subset` / `--subset_file`, which *measures* those relationships
rather than inferring them.

An extraction is a complete bundle: it carries the node ids, the induced MST edges, the
rebuilt hierarchy and the impact metric, and everything else is derived from those on
open, exactly as it is for the bundle it came from.

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
- **The split-event cap bounds the slider and the split chart's stems.** With the default
  50, a very large network's slider carries the 50 strongest merges *in the current
  window* plus up to 20 more for whole-axis coverage — not every one. If the threshold you
  want falls between two stops, zoom into that stretch (the cap re-spends itself there) or
  raise the `Split events` box; no rebuild is needed either way. The caption under the
  chart always says how many of how many are drawn.
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
