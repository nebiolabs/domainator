import heapq
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
    """Axis titles for the split chart, which depend on what merge_impact measures.

    Under ``min_child`` an impact is a count of nodes; under ``product`` it is a product of
    two component sizes, which is not a node count and must not be labelled as one. Both
    the Plotly chart in matrix_report and the canvas chart in the SSN viewer read these so
    the two cannot drift apart.
    """
    if metric == MERGE_IMPACT_PRODUCT:
        return {
            "largest": "Largest single split (size product)",
            "moving_sum": f"Moving sum of split impact ({MOVING_SUM_WINDOW_PERCENT}% window)",
            "moving_sum_short": f"Moving sum ({MOVING_SUM_WINDOW_PERCENT}% window)",
        }
    if metric == MERGE_IMPACT_MIN_CHILD:
        return {
            "largest": "Largest single split (nodes)",
            "moving_sum": f"Moving sum of split size ({MOVING_SUM_WINDOW_PERCENT}% window)",
            "moving_sum_short": f"Moving sum ({MOVING_SUM_WINDOW_PERCENT}% window)",
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

    # Max-heap with lazy deletion tracks the largest component size in O(log n)
    # per merge instead of the O(n) cost of maintaining a sorted list.
    heap = [-1] * tree.n_nodes  # negated for max-heap; all components start at size 1
    heapq.heapify(heap)
    heap_deleted = {}  # size -> pending deletion count

    def heap_remove(size):
        heap_deleted[size] = heap_deleted.get(size, 0) + 1

    def heap_push(size):
        heapq.heappush(heap, -size)

    def heap_max():
        while heap:
            s = -heap[0]
            if heap_deleted.get(s, 0) > 0:
                heap_deleted[s] -= 1
                heapq.heappop(heap)
            else:
                return float(s)
        return 0.0

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
        largest_cluster = heap_max()
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

            heap_remove(left_size)
            heap_remove(right_size)

            if left_size > 1:
                non_singleton_count -= 1
                non_singleton_sum -= left_size
            if right_size > 1:
                non_singleton_count -= 1
                non_singleton_sum -= right_size

            merged_size = left_size + right_size
            parent[left_root] = right_root
            component_sizes[right_root] = merged_size

            heap_push(merged_size)
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


def filter_merge_event_rows(event_rows, max_merge_events=DEFAULT_MAX_MERGE_EVENTS):
    if max_merge_events is None:
        max_merge_events = DEFAULT_MAX_MERGE_EVENTS
    if max_merge_events < 0:
        raise ValueError("max_merge_events must be >= 0")
    if max_merge_events == 0 or len(event_rows) <= max_merge_events:
        return list(event_rows)

    ranked_rows = sorted(
        event_rows,
        key=lambda row: (-row["merge_impact"], -row["delta_largest"], -row["delta_avg_non_singleton"], row["edge_index"]),
    )
    filtered_rows = ranked_rows[:max_merge_events]
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