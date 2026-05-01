import math
from bisect import bisect_left, insort

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


def _remove_sorted_size(sorted_sizes, size):
    position = bisect_left(sorted_sizes, size)
    if position >= len(sorted_sizes) or sorted_sizes[position] != size:
        raise ValueError(f"size {size} missing from sorted component sizes")
    sorted_sizes.pop(position)


def component_size_summary_by_threshold(tree, merge_impact_metric=MERGE_IMPACT_MIN_CHILD):
    if tree.n_nodes == 0:
        return np.zeros((0, 6), dtype=float)

    if merge_impact_metric not in MERGE_IMPACT_CHOICES:
        raise ValueError(f"merge_impact_metric must be one of {sorted(MERGE_IMPACT_CHOICES)}")

    parent = np.arange(tree.n_nodes, dtype=int)
    component_sizes = np.ones(tree.n_nodes, dtype=int)
    sorted_sizes = [1] * tree.n_nodes
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
        largest_cluster = float(sorted_sizes[-1]) if len(sorted_sizes) > 0 else 0.0
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

            _remove_sorted_size(sorted_sizes, left_size)
            _remove_sorted_size(sorted_sizes, right_size)

            if left_size > 1:
                non_singleton_count -= 1
                non_singleton_sum -= left_size
            if right_size > 1:
                non_singleton_count -= 1
                non_singleton_sum -= right_size

            merged_size = left_size + right_size
            parent[left_root] = right_root
            component_sizes[right_root] = merged_size

            insort(sorted_sizes, merged_size)
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
        last_row = component_summary[row_idx]

        while row_idx < len(component_summary) and float(component_summary[row_idx, COMPONENT_THRESHOLD_COL]) == threshold_value:
            merge_impact += float(component_summary[row_idx, COMPONENT_MERGE_IMPACT_COL])
            last_row = component_summary[row_idx]
            row_idx += 1

        edge_index = row_idx - 2

        event_rows.append({
            "edge_index": int(edge_index),
            "summary_row_from": int(first_summary_row_idx),
            "summary_row_to": int(row_idx - 1),
            "threshold_from_value": float(previous_row[COMPONENT_THRESHOLD_COL]),
            "threshold_from": format_threshold_value(previous_row[COMPONENT_THRESHOLD_COL]),
            "threshold_to": format_threshold_value(last_row[COMPONENT_THRESHOLD_COL]),
            "threshold_value": float(last_row[COMPONENT_THRESHOLD_COL]),
            "merge_impact": float(merge_impact),
            "delta_largest": float(abs(last_row[COMPONENT_LARGEST_COL] - previous_row[COMPONENT_LARGEST_COL])),
            "delta_avg_non_singleton": float(abs(last_row[COMPONENT_AVG_NON_SINGLETON_COL] - previous_row[COMPONENT_AVG_NON_SINGLETON_COL])),
        })
        previous_row = last_row

    return event_rows


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