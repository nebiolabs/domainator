import math

import numpy as np


COMPONENT_THRESHOLD_COL = 0
COMPONENT_LARGEST_COL = 1
COMPONENT_AVG_NON_SINGLETON_COL = 2
COMPONENT_MERGE_IMPACT_COL = 3
COMPONENT_DELTA_LARGEST_COL = 4
COMPONENT_DELTA_AVG_NON_SINGLETON_COL = 5

MERGE_IMPACT_PRODUCT = "product"
MERGE_IMPACT_MIN_CHILD = "min_child"
MERGE_IMPACT_CHOICES = (MERGE_IMPACT_PRODUCT, MERGE_IMPACT_MIN_CHILD)
DEFAULT_MAX_MERGE_EVENTS = 500
# The viewer picks its events from the chart's current window rather than from the whole
# axis, so a much smaller number shows more: zoom in and the cap re-spends itself on the
# stretch being looked at. A chart with no viewport (matrix_report, ssn_navigator) has
# only one shot at the whole range and keeps the larger default above.
DEFAULT_WINDOWED_MAX_MERGE_EVENTS = 50
# How many equal bands the plotted threshold axis is cut into when back-filling the
# capped event series. 20 bands = 5% of the axis each, so the filtered series can
# exceed max_merge_events by at most 20 rows -- and every 5% of the axis that has an
# event to show gets one. See filter_merge_event_rows.
MERGE_EVENT_DENSITY_BINS = 20

# Width of the centred moving-sum window, as a fraction of the plotted threshold range,
# and how many points to sample it at. The moving sum answers "how much of the graph is
# decomposing around here overall?", as a counterweight to the per-threshold largest
# single merge that the split chart's stems show.
MOVING_SUM_WINDOW_FRACTION = 0.05
MOVING_SUM_GRID_POINTS = 800
# The window fraction as a whole-number percent, for axis titles and legends.
MOVING_SUM_WINDOW_PERCENT = int(round(MOVING_SUM_WINDOW_FRACTION * 100))


def format_threshold_value(threshold):
    if math.isinf(float(threshold)):
        return "∞"
    return f"{float(threshold):.2f}"


def format_merge_impact_metric(metric: str) -> str:
    if metric == MERGE_IMPACT_PRODUCT:
        return "product"
    if metric == MERGE_IMPACT_MIN_CHILD:
        return "min_child"
    raise ValueError(f"Unsupported merge impact metric: {metric}")


def merge_impact_axis_labels(metric: str) -> dict:
    """Axis titles and prose for the split chart, per what merge_impact measures.

    Under ``min_child`` an impact is a count of nodes; under ``product`` it is a product of
    two component sizes, which is not a node count and must not be labelled as one. Both
    the Plotly chart in matrix_report and the canvas chart in the Domainator Similarity
    Network Viewer read these so the two cannot drift apart.

    ``impact_amount`` is a one-placeholder template for naming an impact *value* in a
    sentence -- a bare ``{}`` rather than Plotly's or JavaScript's interpolation syntax,
    since both charts' hover readouts fill it in their own way. ``moving_sum_hover``
    introduces the moving-sum readout, which measures the same quantity.
    """
    if metric == MERGE_IMPACT_PRODUCT:
        return {
            "largest": "Largest single split (size product)",
            "moving_sum": f"Moving sum of split impact ({MOVING_SUM_WINDOW_PERCENT}% window)",
            "moving_sum_short": f"Moving sum ({MOVING_SUM_WINDOW_PERCENT}% window)",
            "impact_amount": "impact {}",
            "moving_sum_hover": "Split impact within",
        }
    if metric == MERGE_IMPACT_MIN_CHILD:
        return {
            "largest": "Largest single split (nodes)",
            "moving_sum": f"Moving sum of split size ({MOVING_SUM_WINDOW_PERCENT}% window)",
            "moving_sum_short": f"Moving sum ({MOVING_SUM_WINDOW_PERCENT}% window)",
            "impact_amount": "{} nodes",
            "moving_sum_hover": "Nodes displaced within",
        }
    raise ValueError(f"Unsupported merge impact metric: {metric}")


def _merge_size_key(size):
    """Histogram key for one individual merge size.

    Integral sizes key as ``int`` so the payload carries ``{"1": 33}`` rather than
    ``{"1.0": 33}``, which measurably shrinks the embedded JSON. Both current metrics are
    integral (``min_child`` is a min of two component sizes, ``product`` their product), but
    a non-integral size from some future metric keys by its full float value rather than
    being truncated by ``int()``, which would silently collide distinct sizes.
    """
    value = float(size)
    if value.is_integer():
        return int(value)
    return value


def component_size_summary_by_threshold(tree, merge_impact_metric=MERGE_IMPACT_MIN_CHILD):
    if tree.n_nodes == 0:
        return np.zeros((0, 6), dtype=float)

    if merge_impact_metric not in MERGE_IMPACT_CHOICES:
        raise ValueError(f"merge_impact_metric must be one of {sorted(MERGE_IMPACT_CHOICES)}")

    parent = np.arange(tree.n_nodes, dtype=int)
    component_sizes = np.ones(tree.n_nodes, dtype=int)

    largest_cluster = 1
    non_singleton_count = 0
    non_singleton_sum = 0
    summary = np.zeros((len(tree.mst_edges) + 1, 6), dtype=float)

    def find(node_idx):
        root = node_idx
        while parent[root] != root:
            root = parent[root]
        while parent[node_idx] != node_idx:
            next_idx = parent[node_idx]
            parent[node_idx] = root
            node_idx = next_idx
        return root

    def record(row_idx, threshold, merge_impact=0.0):
        nonlocal non_singleton_count, non_singleton_sum
        avg_non_singleton = float(non_singleton_sum) / non_singleton_count if non_singleton_count > 0 else 0.0

        summary[row_idx, COMPONENT_THRESHOLD_COL] = threshold
        summary[row_idx, COMPONENT_LARGEST_COL] = largest_cluster
        summary[row_idx, COMPONENT_AVG_NON_SINGLETON_COL] = avg_non_singleton
        summary[row_idx, COMPONENT_MERGE_IMPACT_COL] = float(merge_impact)
        if row_idx > 0:
            summary[row_idx, COMPONENT_DELTA_LARGEST_COL] = largest_cluster - summary[row_idx - 1, COMPONENT_LARGEST_COL]
            summary[row_idx, COMPONENT_DELTA_AVG_NON_SINGLETON_COL] = avg_non_singleton - summary[row_idx - 1, COMPONENT_AVG_NON_SINGLETON_COL]

    record(0, float("inf"))

    for row_idx, (source_idx, target_idx, threshold) in enumerate(tree.mst_edges, start=1):
        left_root = find(source_idx)
        right_root = find(target_idx)
        merge_impact = 0.0
        if left_root != right_root:
            left_size = int(component_sizes[left_root])
            right_size = int(component_sizes[right_root])
            if merge_impact_metric == MERGE_IMPACT_PRODUCT:
                merge_impact = float(left_size * right_size)
            else:
                merge_impact = float(min(left_size, right_size))

            if left_size > 1:
                non_singleton_count -= 1
                non_singleton_sum -= left_size
            if right_size > 1:
                non_singleton_count -= 1
                non_singleton_sum -= right_size

            merged_size = left_size + right_size
            parent[left_root] = right_root
            component_sizes[right_root] = merged_size

            largest_cluster = max(largest_cluster, merged_size)
            non_singleton_count += 1
            non_singleton_sum += merged_size

        record(row_idx, float(threshold), merge_impact=merge_impact)

    return summary


def threshold_merge_event_rows(component_summary):
    if component_summary is None or len(component_summary) < 2:
        return []

    event_rows = []
    previous_row = component_summary[0]
    row_idx = 1

    while row_idx < len(component_summary):
        first_summary_row_idx = row_idx
        threshold_value = float(component_summary[row_idx, COMPONENT_THRESHOLD_COL])
        merge_impact = 0.0
        # The individual per-edge impacts behind that sum. Summing alone cannot distinguish
        # one large cluster splitting off from a swarm of tiny ones, so keep the terms too.
        merge_size_counts = {}
        merge_count = 0
        largest_merge = 0.0
        last_row = component_summary[row_idx]

        while row_idx < len(component_summary) and float(component_summary[row_idx, COMPONENT_THRESHOLD_COL]) == threshold_value:
            row_impact = float(component_summary[row_idx, COMPONENT_MERGE_IMPACT_COL])
            merge_impact += row_impact
            if row_impact > 0.0:
                # A zero impact means the edge's endpoints were already in one component,
                # so it joined nothing and is not a merge.
                size_key = _merge_size_key(row_impact)
                merge_size_counts[size_key] = merge_size_counts.get(size_key, 0) + 1
                merge_count += 1
                largest_merge = max(largest_merge, row_impact)
            last_row = component_summary[row_idx]
            row_idx += 1

        # The MST edge index of the last edge scoring strictly above this threshold, so a
        # stop labelled `threshold_to` reproduces the `--lb threshold_to` cut exactly (the
        # tie group at that threshold is excluded, as `--lb` excludes it).
        edge_index = first_summary_row_idx - 2

        event_rows.append({
            "edge_index": int(edge_index),
            "summary_row_from": int(first_summary_row_idx),
            "summary_row_to": int(row_idx - 1),
            "threshold_from_value": float(previous_row[COMPONENT_THRESHOLD_COL]),
            "threshold_from": format_threshold_value(previous_row[COMPONENT_THRESHOLD_COL]),
            "threshold_to": format_threshold_value(last_row[COMPONENT_THRESHOLD_COL]),
            "threshold_value": float(last_row[COMPONENT_THRESHOLD_COL]),
            "merge_impact": float(merge_impact),
            # merge_impact is the sum over the tie group; these three describe its terms.
            # largest_merge and merge_count are emitted rather than derived downstream so
            # they stay correct even if merge_size_counts is ever truncated for payload size.
            "merge_size_counts": merge_size_counts,
            "largest_merge": float(largest_merge),
            "merge_count": int(merge_count),
            "delta_largest": float(abs(last_row[COMPONENT_LARGEST_COL] - previous_row[COMPONENT_LARGEST_COL])),
            "delta_avg_non_singleton": float(abs(last_row[COMPONENT_AVG_NON_SINGLETON_COL] - previous_row[COMPONENT_AVG_NON_SINGLETON_COL])),
        })
        previous_row = last_row

    return event_rows


def merge_event_moving_sum(event_rows,
                           window_fraction=MOVING_SUM_WINDOW_FRACTION,
                           grid_points=MOVING_SUM_GRID_POINTS):
    """Centred moving sum of ``merge_impact`` over the plotted threshold range.

    Call this with the rows straight from :func:`threshold_merge_event_rows`, *before*
    :func:`filter_merge_event_rows` caps them. Filtering keeps the top rows ranked by
    ``merge_impact``, so a moving sum taken afterwards undercounts on any network with
    more than ``max_merge_events`` threshold groups -- and it undercounts precisely the
    small events this series exists to reveal.

    Sums ``merge_impact`` (the per-threshold total), not ``largest_merge``: the contrast
    between the two is what the split chart is for. Returns
    ``{"window": W, "x": [...], "y": [...]}``, with empty series when there is no finite
    threshold range to slide a window over.
    """
    empty = {"window": 0.0, "x": [], "y": []}
    if not event_rows or grid_points < 2:
        return empty

    thresholds = np.array([float(row["threshold_value"]) for row in event_rows], dtype=float)
    impacts = np.array([float(row["merge_impact"]) for row in event_rows], dtype=float)
    finite = np.isfinite(thresholds) & np.isfinite(impacts)
    if not finite.any():
        return empty
    thresholds = thresholds[finite]
    impacts = impacts[finite]

    lo = float(thresholds.min())
    hi = float(thresholds.max())
    if not math.isfinite(lo) or not math.isfinite(hi) or hi == lo:
        return empty

    window = float(window_fraction) * (hi - lo)
    half_window = window / 2.0

    # searchsorted against a prefix sum turns the window scan into O(n log n + grid)
    # instead of the O(n * grid) a nested loop would cost.
    order = np.argsort(thresholds, kind="stable")
    sorted_thresholds = thresholds[order]
    cumulative = np.concatenate(([0.0], np.cumsum(impacts[order])))

    grid = np.linspace(lo, hi, int(grid_points))
    # 'left'/'right' make the window inclusive at both ends, matching |t - g| <= W/2.
    starts = np.searchsorted(sorted_thresholds, grid - half_window, side="left")
    ends = np.searchsorted(sorted_thresholds, grid + half_window, side="right")
    totals = cumulative[ends] - cumulative[starts]

    return {
        "window": window,
        "x": [float(value) for value in grid],
        "y": [int(value) if float(value).is_integer() else float(value) for value in totals],
    }


def summarize_merge_events(component_summary, max_items=5):
    if component_summary is None or len(component_summary) < 2:
        return []

    merge_rows = threshold_merge_event_rows(component_summary)
    ranked_rows = sorted(merge_rows, key=lambda row: (-row["merge_impact"], -row["delta_largest"], -row["delta_avg_non_singleton"]))
    return ranked_rows[:max_items]


def merge_event_rank_key(row):
    """Strongest first. Shared by the cap and the density back-fill so both agree.

    Ported to JavaScript in `ssn_viewer_html.py` (`compareExtractionMergeEventRank`)
    for viewer-built extractions; the two must stay in step or an extraction's plot
    will not match the plot it was extracted from.
    """
    return (
        -row["merge_impact"],
        -row["delta_largest"],
        -row["delta_avg_non_singleton"],
        row["edge_index"],
    )


def merge_event_density_bin(threshold_value, lo, hi, density_bins):
    """Which band of the plotted axis a threshold falls in, or -1 if the axis is a point.

    `hi` itself lands in the last band rather than one past the end.
    """
    if density_bins < 1 or not math.isfinite(lo) or not math.isfinite(hi) or hi <= lo:
        return -1
    position = int(((float(threshold_value) - lo) / (hi - lo)) * density_bins)
    return max(0, min(density_bins - 1, position))


def filter_merge_event_rows(event_rows, max_merge_events=DEFAULT_MAX_MERGE_EVENTS,
                            density_bins=MERGE_EVENT_DENSITY_BINS,
                            window=None, pinned_threshold=None):
    """The strongest `max_merge_events` events, plus enough to keep the axis covered.

    Ranking by impact alone leaves long stretches of the threshold axis with nothing
    plotted on them. On a connected MST-kNN graph the weak tail of MST edges is
    individual outliers being attached to the giant component, one or two nodes at a
    time; their impact is the smallest there is, so the cap drops every one of them --
    while the chart's x-axis and the threshold slider still span the real threshold
    range, because the moving sum and the floor stop are both derived from the
    *unfiltered* data. The result is a plot whose left half is blank and a slider whose
    left half has no stops on it.

    So after the top-N pass the axis is cut into `density_bins` equal bands and the
    strongest event in each otherwise-empty band is added back. That bounds the output
    at `max_merge_events + density_bins` rows, guarantees no band is silently empty
    while it still has an event to offer, and leaves the top-N rows themselves
    untouched -- a back-filled row is an addition, never a replacement.

    `max_merge_events=0` means no cap, and no back-fill either: nothing was dropped.

    `window` is an optional ``(low, high)`` threshold range -- the part of the axis a
    zoomable chart is currently showing. It narrows only the **top-N pool**: the cap
    then picks the strongest events *within the window*, so zooming in reveals the small
    events a whole-range ranking buries, while the band back-fill still runs over the
    whole axis so every stretch of the slider keeps a stop. ``None`` is the whole range,
    which is what a chart without a viewport (``matrix_report``, ``ssn_navigator``) wants
    and is exactly the behaviour this function had before windows existed.

    `pinned_threshold` is a threshold whose event must survive the cut whatever the
    window is. The viewer passes the cut currently in effect: if the selection could
    drop it, zooming would silently re-cluster the network under the user.
    """
    if max_merge_events is None:
        max_merge_events = DEFAULT_MAX_MERGE_EVENTS
    if max_merge_events < 0:
        raise ValueError("max_merge_events must be >= 0")
    if density_bins < 0:
        raise ValueError("density_bins must be >= 0")
    if max_merge_events == 0 or len(event_rows) <= max_merge_events:
        return list(event_rows)

    ranked_rows = sorted(event_rows, key=merge_event_rank_key)

    if window is None:
        pool = ranked_rows
    else:
        low, high = (float(window[0]), float(window[1]))
        if low > high:
            low, high = high, low
        pool = [row for row in ranked_rows if low <= float(row["threshold_value"]) <= high]

    filtered_rows = pool[:max_merge_events]
    selected_ids = {id(row) for row in filtered_rows}

    # The band edges come from the whole series, because that is what the axis spans --
    # a window narrows what the cap ranks, never what the axis has to cover.
    thresholds = [value for value in (float(row["threshold_value"]) for row in event_rows)
                  if math.isfinite(value)]
    if thresholds and density_bins > 0:
        lo = min(thresholds)
        hi = max(thresholds)
        covered = {
            merge_event_density_bin(row["threshold_value"], lo, hi, density_bins)
            for row in filtered_rows
        }
        # ranked_rows is strongest-first, so the first row seen in an empty band is
        # the strongest one available to represent it. Candidates come from the whole
        # series, not the window: a band outside the window still needs its stop.
        for row in ranked_rows:
            if id(row) in selected_ids:
                continue
            bin_index = merge_event_density_bin(row["threshold_value"], lo, hi, density_bins)
            if bin_index < 0 or bin_index in covered:
                continue
            covered.add(bin_index)
            filtered_rows.append(row)
            selected_ids.add(id(row))

    if pinned_threshold is not None and math.isfinite(float(pinned_threshold)):
        pinned = float(pinned_threshold)
        if not any(float(row["threshold_value"]) == pinned for row in filtered_rows):
            for row in ranked_rows:
                if float(row["threshold_value"]) == pinned:
                    filtered_rows.append(row)
                    break

    filtered_rows.sort(key=lambda row: row["edge_index"])
    return filtered_rows


def merge_event_table_rows(component_summary, max_items=25):
    rows = []
    strongest_rows = summarize_merge_events(component_summary, max_items=max_items)
    strongest_rows.sort(key=lambda row: (row["threshold_from_value"], row["threshold_value"]), reverse=True)
    for merge_row in strongest_rows:
        rows.append({
            "threshold_from": merge_row["threshold_from"],
            "threshold_to": merge_row["threshold_to"],
            "merge_impact": int(round(merge_row["merge_impact"])),
        })
    return rows


FLOOR_THRESHOLD_SPAN_FRACTION = 0.01


def floor_threshold_value(lowest_merge_threshold: float, highest_merge_threshold: float = None) -> float:
    """A cut strictly below the weakest merge, so every edge is kept.

    Offset below the weakest merge by 1% of the merge-weight range, rather than
    dropped to ``0``. A network whose scores run 350-650 would otherwise spend more
    than half its slider track on an empty stretch below the data, and scores that
    go negative are not cleared by ``0`` at all.

    ``highest_merge_threshold`` is optional only so a caller with a single weight in
    hand still gets a usable floor; pass both whenever the range is known.
    """
    lowest = float(lowest_merge_threshold)
    highest = lowest if highest_merge_threshold is None else float(highest_merge_threshold)
    span = highest - lowest
    if span > 0:
        return lowest - (FLOOR_THRESHOLD_SPAN_FRACTION * span)
    # Every merge sits at one weight, so there is no range to take a fraction of.
    step = abs(lowest) * FLOOR_THRESHOLD_SPAN_FRACTION
    return lowest - (step if step > 0 else 1.0)


def threshold_slider_stops(merge_event_rows, tree=None, threshold_index_lookup=None):
    """Slider stops for a merge-event series: ∞, one per event, then a floor.

    Shared by ``matrix_report``, the ``.dsnv`` viewer (via its JavaScript port) and
    ``ssn_bundle`` so their sliders offer the same cuts. The trailing floor stop is
    what makes the fully merged network reachable at all: every other stop excludes
    its own tie group under the strictly-above ``--lb`` convention, so the lowest
    event stop still splits the weakest merge back apart. It is derived from the MST
    rather than from ``merge_event_rows``, which a cap may have thinned.

    ``threshold_index_lookup`` is an optional ``threshold -> row index`` callable
    (in practice :meth:`MaxTree.threshold_row_index`). When given, every stop also
    carries a ``threshold_index`` pointing into the tree's threshold tables, which is
    how ``matrix_report`` reports each cut's edge count. A ``.dsnv`` bundle carries no
    threshold tables, so the viewer's stops omit the field entirely rather than
    carrying an index into something that is not there.
    """
    stops = [{
        "edge_index": -1,
        "threshold_label": "∞",
        "threshold_value": None,
    }]
    if threshold_index_lookup is not None:
        stops[0]["threshold_index"] = -1
    for merge_row in merge_event_rows:
        stop = {
            "edge_index": int(merge_row["edge_index"]),
            "threshold_label": merge_row["threshold_to"],
            "threshold_value": float(merge_row["threshold_value"]),
        }
        if threshold_index_lookup is not None:
            # Row in the tree's threshold tables for this cut. Those are keyed by
            # distinct threshold, not by MST edge, so the lookup is separate.
            stop["threshold_index"] = threshold_index_lookup(merge_row["threshold_value"])
        stops.append(stop)

    mst_edges = [] if tree is None else list(tree.mst_edges)
    if len(mst_edges) > 0:
        # mst_edges are weight-descending.
        floor_value = floor_threshold_value(float(mst_edges[-1][2]), float(mst_edges[0][2]))
        floor_stop = {
            # Every MST edge scores strictly above this cut.
            "edge_index": len(mst_edges) - 1,
            "threshold_label": format_threshold_value(floor_value),
            "threshold_value": floor_value,
        }
        if threshold_index_lookup is not None:
            # This stop stands for the complete-graph cut, so it reads its edge counts
            # from that row of the threshold tables. Its own value sits just below the
            # weakest merge rather than at 0 purely so the slider track stays usable,
            # and no table row is keyed there. The counts are exact whenever no graph
            # edge scores at or below the floor, which is the usual case.
            floor_stop["threshold_index"] = threshold_index_lookup(0.0)
        stops.append(floor_stop)
    return stops


def build_mst_component_hierarchy(tree):
    parent = np.arange(tree.n_nodes, dtype=int)
    hierarchy_nodes = []
    component_id_by_root = {}
    component_size_by_id = {}
    component_min_leaf_by_id = {}

    for node_index in range(tree.n_nodes):
        hierarchy_nodes.append({
            "id": node_index,
            "kind": "leaf",
            "node_index": node_index,
            "size": 1,
            "parent": None,
        })
        component_id_by_root[node_index] = node_index
        component_size_by_id[node_index] = 1
        component_min_leaf_by_id[node_index] = node_index

    def find(node_index):
        root = node_index
        while parent[root] != root:
            root = parent[root]
        while parent[node_index] != node_index:
            next_index = parent[node_index]
            parent[node_index] = root
            node_index = next_index
        return root

    next_component_id = tree.n_nodes
    for source_idx, target_idx, threshold in tree.mst_edges:
        left_root = find(source_idx)
        right_root = find(target_idx)
        if left_root == right_root:
            continue

        left_component_id = component_id_by_root[left_root]
        right_component_id = component_id_by_root[right_root]
        if component_min_leaf_by_id[left_component_id] > component_min_leaf_by_id[right_component_id]:
            left_component_id, right_component_id = right_component_id, left_component_id
        merged_size = component_size_by_id[left_component_id] + component_size_by_id[right_component_id]
        component_id = next_component_id
        next_component_id += 1

        hierarchy_nodes[left_component_id]["parent"] = component_id
        hierarchy_nodes[right_component_id]["parent"] = component_id
        hierarchy_nodes.append({
            "id": component_id,
            "kind": "cluster",
            "left": left_component_id,
            "right": right_component_id,
            "threshold": float(threshold),
            "size": int(merged_size),
            "parent": None,
        })
        component_size_by_id[component_id] = int(merged_size)
        component_min_leaf_by_id[component_id] = min(
            component_min_leaf_by_id[left_component_id],
            component_min_leaf_by_id[right_component_id],
        )

        parent[left_root] = right_root
        component_id_by_root.pop(left_root)
        component_id_by_root[right_root] = component_id

    roots = sorted(component_id_by_root.values(), key=lambda component_id: (-hierarchy_nodes[component_id]["size"], component_id))

    leaf_order = []
    stack = []
    for root_id in reversed(roots):
        stack.append((root_id, False))

    while stack:
        component_id, visited = stack.pop()
        node = hierarchy_nodes[component_id]
        if node["kind"] == "leaf":
            node["leaf_start"] = len(leaf_order)
            node["leaf_count"] = 1
            leaf_order.append(node["node_index"])
            continue

        if visited:
            left_node = hierarchy_nodes[node["left"]]
            right_node = hierarchy_nodes[node["right"]]
            node["leaf_start"] = left_node["leaf_start"]
            node["leaf_count"] = left_node["leaf_count"] + right_node["leaf_count"]
            continue

        stack.append((component_id, True))
        stack.append((node["right"], False))
        stack.append((node["left"], False))

    return {
        "nodes": hierarchy_nodes,
        "roots": roots,
        "leaf_order": leaf_order,
    }