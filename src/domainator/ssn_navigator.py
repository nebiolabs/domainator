"""Navigate a build_ssn_viewer ``.dsnv`` bundle and emit compact JSON summaries.

This is a read-only, agent-friendly companion to the interactive HTML SSN
viewer. Instead of rendering the whole network, it answers targeted questions
about the sequence-similarity network at a chosen similarity threshold and
returns compact JSON: which clusters exist, what a cluster contains, and how a
cluster's metadata is distributed. Higher thresholds give finer clusters (more
splitting), matching the viewer's slider: ``--threshold inf`` is the "∞" stop at
its fine end, where every node is its own cluster, and omitting ``--threshold``
uses the floor at the coarse end, one cluster per connected component.

The cut-points offered by ``thresholds`` and ``node`` are derived from the bundle's
merge order, not read out of it, so ``--max_merge_events`` controls how many are offered
(``0`` for one per split event) whatever the bundle was built with.

Modes:
  overview    network size, metadata columns, component count, defaults.
  thresholds  suggested similarity cut-points with the resulting cluster counts.
  clusters    list the clusters active at --threshold (with sizes).
  cluster     detail one cluster (--id or --node): size, members, metadata.
  node        a node's metadata and the cluster it falls in at each cut-point.
"""

import json
import math
import sys

from jsonargparse import ActionConfigFile, ArgumentParser

from domainator import __version__, RawAndDefaultsFormatter
from domainator.ssn_bundle import (
    clusters_at_threshold,
    coarsest_threshold,
    component_members,
    load_bundle,
    merge_event_rows,
    merge_event_series,
    node_index_by_name,
    slider_stops,
    summarize_cluster_metadata,
)
from domainator.ssn_hierarchy import DEFAULT_MAX_MERGE_EVENTS

MODES = ("overview", "thresholds", "clusters", "cluster", "node")


def _metadata_row(bundle: dict, node_index: int) -> dict:
    columns = bundle["metadata"].get("columns", [])
    rows = bundle["metadata"].get("rows", [])
    if node_index >= len(rows):
        return {}
    row = rows[node_index]
    return {column["name"]: row[idx] for idx, column in enumerate(columns)}


def _cluster_records(hierarchy, active_ids, min_size=1):
    nodes = hierarchy["nodes"]
    records = [
        {"id": component_id, "size": nodes[component_id]["size"]}
        for component_id in active_ids
        if nodes[component_id]["size"] >= min_size
    ]
    records.sort(key=lambda record: (-record["size"], record["id"]))
    return records


def _merge_event_coverage(event_rows, selected_rows, max_merge_events):
    """How much of the network's split-event series this run is reporting.

    The cut-points offered are one per *selected* split event, and ``--max_merge_events``
    keeps only the strongest (plus a back-fill so no stretch of the threshold range goes
    unrepresented). A caller that treats the offered cut-points as the complete set of
    interesting thresholds is wrong on any large network, so say what was left out.

    The series itself is derived from the merge order rather than read from the bundle,
    so the total is exact on every bundle version -- including pre-v6 files, whose stored
    series was already capped when it was written.
    """
    return {
        "merge_events": len(selected_rows),
        "merge_event_total": len(event_rows),
        "max_merge_events": int(max_merge_events),
    }


def _threshold_summary(bundle, stops):
    hierarchy = bundle["graph"]["hierarchy"]
    summary = []
    for stop in stops:
        threshold_value = stop["threshold_value"]
        active = clusters_at_threshold(hierarchy, threshold_value)
        sizes = [hierarchy["nodes"][cid]["size"] for cid in active]
        non_singleton = sum(1 for size in sizes if size > 1)
        summary.append({
            "threshold_value": threshold_value,
            "threshold_label": stop["threshold_label"],
            "clusters": len(active),
            "non_singleton_clusters": non_singleton,
            "largest_cluster": max(sizes) if sizes else 0,
        })
    return summary


def _resolve_threshold(threshold, hierarchy):
    """Resolve ``--threshold`` to a numeric cut, following the viewer's slider.

    ``inf``/``∞`` is the slider's ∞ stop -- every merge cut, every node its own
    cluster. Omitting ``--threshold`` picks the other end instead, the floor
    below every merge, since one cluster per connected component is the more
    useful thing to see without being asked for a threshold.
    """
    if threshold is None:
        return coarsest_threshold(hierarchy)
    if isinstance(threshold, str):
        if threshold.lower() in ("inf", "infinity", "∞"):
            return math.inf
        threshold = float(threshold)
    return float(threshold)


def _threshold_out(threshold):
    """The threshold as it should appear in output JSON.

    ``null`` for the ∞ cut, which is how the bundle's own slider_stops spell it,
    and which keeps the output strict JSON -- ``json.dumps(inf)`` emits the
    non-standard ``Infinity``.
    """
    return None if math.isinf(threshold) else threshold


def _cluster_for_node(hierarchy, node_index, threshold):
    """The cluster containing ``node_index`` at ``threshold``, by walking up from its leaf.

    Asking which cluster holds one node does not need the whole partition: climb from the
    node's leaf while the parent merge scores strictly above the cut, and the last node
    reached is the cluster :func:`clusters_at_threshold` would have emitted for it. A
    cluster's ancestors merged later and therefore score lower, so the condition fails
    once and the walk stops -- no need to look at any other branch of the tree.

    This is O(depth) rather than O(clusters), which matters because ``--mode node``
    repeats the question at every cut-point. The ``>`` is the same strictly-above rule
    the descent uses, so a merge sitting exactly on the cut is not taken.

    ``threshold`` is ``None`` for the ∞ cut, where every node stands alone.
    """
    nodes = hierarchy["nodes"]
    if threshold is None:
        return node_index
    cut = float(threshold)
    current = node_index
    while True:
        parent = nodes[current]["parent"]
        if parent is None or not (float(nodes[parent]["threshold"]) > cut):
            return current
        current = parent


def ssn_navigator(bundle, mode, threshold=None, cluster_id=None, node=None,
                  min_size=1, top_n=50, members=False,
                  max_merge_events=DEFAULT_MAX_MERGE_EVENTS):
    """Answer a single navigation query against a loaded ``.dsnv`` bundle dict.

    Returns a JSON-serializable dict. See the module docstring for modes.
    """
    graph = bundle["graph"]
    hierarchy = graph["hierarchy"]
    node_names = graph["nodes"]
    threshold = _resolve_threshold(threshold, hierarchy)

    # Replaying the merge order costs a pass over the MST edges, so only the modes that
    # report cut-points pay for it. Cached so a mode needing both the rows and the stops
    # replays once.
    derived = {}

    def event_rows():
        if "rows" not in derived:
            derived["rows"] = merge_event_rows(bundle)
        return derived["rows"]

    def selected_rows():
        if "selected" not in derived:
            derived["selected"] = merge_event_series(
                bundle, max_merge_events=max_merge_events, event_rows=event_rows()
            )
        return derived["selected"]

    def stops():
        if "stops" not in derived:
            derived["stops"] = slider_stops(
                bundle, max_merge_events=max_merge_events, event_rows=event_rows()
            )
        return derived["stops"]

    if mode == "overview":
        return {
            "name": bundle.get("name"),
            "domainator_version": bundle.get("domainator_version"),
            "nodes": len(node_names),
            "mst_edges": len(graph.get("mst_edges", [])),
            "connected_components": len(hierarchy["roots"]),
            "metadata_columns": bundle["metadata"].get("columns", []),
            "defaults": bundle.get("defaults", {}),
            "merge_impact_metric": graph.get("merge_impact_metric"),
            **_merge_event_coverage(event_rows(), selected_rows(), max_merge_events),
        }

    if mode == "thresholds":
        # The cut-points are one per selected merge event, so they inherit the cap. Say
        # so: picking a cut-point from a silently truncated list is the failure this
        # guards against. Raise --max_merge_events (0 for all of them) to see more.
        return {
            "thresholds": _threshold_summary(bundle, stops()),
            **_merge_event_coverage(event_rows(), selected_rows(), max_merge_events),
        }

    if mode == "clusters":
        active = clusters_at_threshold(hierarchy, threshold)
        records = _cluster_records(hierarchy, active, min_size=min_size)
        singletons = sum(1 for r in records if r["size"] == 1)
        result = {
            "threshold": _threshold_out(threshold),
            "clusters": len(records),
            "singletons": singletons,
            "shown": min(len(records), top_n),
            "cluster_list": records[:top_n],
        }
        if len(records) > top_n:
            result["truncated"] = len(records) - top_n
        return result

    if mode == "cluster":
        active = clusters_at_threshold(hierarchy, threshold)
        if node is not None:
            index_by_name = node_index_by_name(bundle)
            if node not in index_by_name:
                raise ValueError(f"Node '{node}' not found in the bundle.")
            cluster_id = _cluster_for_node(hierarchy, index_by_name[node], threshold)
            if cluster_id is None:
                raise ValueError(f"No active cluster found for node '{node}'.")
        elif cluster_id is not None:
            if cluster_id not in active:
                raise ValueError(
                    f"Cluster id {cluster_id} is not an active cluster at "
                    f"threshold {threshold}; use --mode clusters to list valid ids."
                )
        else:
            raise ValueError("cluster mode requires either --id or --node.")

        member_indices = component_members(hierarchy, cluster_id)
        result = {
            "threshold": _threshold_out(threshold),
            "cluster_id": cluster_id,
            "size": len(member_indices),
            "metadata_summary": summarize_cluster_metadata(
                member_indices, bundle["metadata"], top_n=top_n
            ),
        }
        if members:
            result["members"] = [node_names[i] for i in member_indices]
        return result

    if mode == "node":
        if node is None:
            raise ValueError("node mode requires --node.")
        index_by_name = node_index_by_name(bundle)
        if node not in index_by_name:
            raise ValueError(f"Node '{node}' not found in the bundle.")
        node_index = index_by_name[node]
        memberships = []
        for stop in stops():
            stop_threshold = stop["threshold_value"]
            # One walk up from this node's leaf per cut-point, rather than rebuilding the
            # whole partition at each one only to look up a single member.
            cid = _cluster_for_node(hierarchy, node_index, stop_threshold)
            memberships.append({
                "threshold_value": stop_threshold,
                "threshold_label": stop["threshold_label"],
                "cluster_id": cid,
                "cluster_size": hierarchy["nodes"][cid]["size"] if cid is not None else None,
            })
        return {
            "node": node,
            "node_index": node_index,
            "metadata": _metadata_row(bundle, node_index),
            "memberships": memberships,
        }

    raise ValueError(f"Unknown mode '{mode}'. Choose one of {MODES}.")


def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__,
                            formatter_class=RawAndDefaultsFormatter)
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="A Domainator Similarity Network Viewer bundle (.dsnv) produced by build_ssn_viewer.py.")
    parser.add_argument("-o", "--output", type=str, default=None,
                        help="Write JSON to this file. Default: stdout.")
    parser.add_argument("--mode", type=str, default="overview", choices=MODES,
                        help="Which navigation query to run.")
    parser.add_argument("--threshold", type=str, default=None,
                        help="Similarity threshold to cut the network at (higher = finer "
                             "clusters), matching the viewer's slider. Pass 'inf' for the "
                             "slider's fine end, where every node is its own cluster. Omit "
                             "it for the coarse end: one cluster per connected component.")
    parser.add_argument("--id", type=int, default=None, dest="cluster_id",
                        help="Cluster id to describe (see --mode clusters). For --mode cluster.")
    parser.add_argument("--node", type=str, default=None,
                        help="Node id/name to describe or locate. For --mode cluster/node.")
    parser.add_argument("--min_size", type=int, default=1,
                        help="Only report clusters with at least this many members (--mode clusters).")
    parser.add_argument("--top_n", type=int, default=50,
                        help="Bound list output: at most this many clusters, category values, etc.")
    parser.add_argument("--members", action="store_true", default=False,
                        help="Include the full member id list in --mode cluster output.")
    parser.add_argument("--max_merge_events", type=int, default=DEFAULT_MAX_MERGE_EVENTS,
                        help="How many split events the --mode thresholds/node cut-points are drawn from: the strongest this many by impact, plus a few so that no 5%% band of the threshold range with an event to offer is left without one. 0 offers a cut-point at every split event. The series is derived from the bundle's merge order, so any value works on any bundle.")
    parser.add_argument("--config", action=ActionConfigFile)
    params = parser.parse_args(argv)

    bundle = load_bundle(params.input)
    result = ssn_navigator(
        bundle,
        mode=params.mode,
        threshold=params.threshold,
        cluster_id=params.cluster_id,
        node=params.node,
        min_size=params.min_size,
        top_n=params.top_n,
        members=params.members,
        max_merge_events=params.max_merge_events,
    )

    out_handle = open(params.output, "w") if params.output is not None else sys.stdout
    try:
        json.dump(result, out_handle, separators=(",", ":"))
        out_handle.write("\n")
    finally:
        if out_handle is not sys.stdout:
            out_handle.close()


def _entrypoint():
    main(sys.argv[1:])


if __name__ == "__main__":
    _entrypoint()
