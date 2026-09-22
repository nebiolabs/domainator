"""Read-only helpers for navigating a build_ssn_viewer ``.dsnv`` bundle.

The bundle is a gzip-compressed JSON document (see ``build_ssn_viewer.py``) that
stores an agglomerative MST merge tree plus a positional metadata table. These
helpers let non-browser callers (notably ``ssn_navigator.py``) reconstruct the
clusters that are active at a chosen similarity threshold and summarize the
metadata of their members, mirroring the logic of the JavaScript viewer's
``activeClustersAtThreshold`` / ``componentMembers`` functions.

From v6 the bundle carries only what cannot be derived from the merge order, so the
split-event series, its moving sum and the threshold slider's stops are computed here
on demand -- :func:`merge_event_rows`, :func:`merge_event_series`, :func:`moving_sum`
and :func:`slider_stops` -- rather than read out of the file. The viewer does the same
thing in JavaScript. A caller therefore picks its own ``max_merge_events`` at the point
of display instead of inheriting one chosen when the file was written.

A v4 bundle may also carry a top-level ``app_state`` section, written by the HTML
viewer's "Save session" button. It holds viewer UI state only and is ignored
here; the rest of such a file is an ordinary bundle, so a saved session's
``metadata`` simply reflects whatever was edited in the viewer.
"""

import gzip
import json
from os import PathLike
from typing import Dict, List, Optional, Sequence, Union

import numpy as np

from domainator.ssn_hierarchy import (
    DEFAULT_MAX_MERGE_EVENTS,
    MERGE_IMPACT_MIN_CHILD,
    component_size_summary_by_threshold,
    filter_merge_event_rows,
    floor_threshold_value,
    merge_event_moving_sum,
    threshold_merge_event_rows,
    threshold_slider_stops,
)


SSN_VIEWER_BUNDLE_FORMAT = "domainator_ssn_viewer_bundle"
# v3: per-event merge_size_counts/largest_merge/merge_count + graph.merge_moving_sum
# v4: optional top-level "app_state" section written by the HTML viewer's "Save
#     session" button (viewer UI state + selection presets). Purely additive:
#     build_ssn_viewer.py never writes app_state, and a reader that ignores the
#     section can treat a v4 file exactly like a v3 file.
# v5: graph.merge_event_total and graph.max_merge_events, so a reader can say how
#     many of the network's split events the capped graph.merge_event_series is
#     showing. Purely additive: a reader that ignores both keys treats a v5 file
#     exactly like a v4 file.
# v6: the derived threshold data is gone. graph.cluster_count_by_threshold and
#     graph.edges_by_threshold described the *full* input graph, which a bundle never
#     carried and no viewer ever read. graph.merge_event_series, merge_event_total,
#     max_merge_events, merge_moving_sum and slider_stops (with its threshold_index)
#     are all a pure function of (nodes, mst_edges, merge_impact_metric), so readers
#     derive them instead -- see merge_event_rows / slider_stops below, and their
#     JavaScript counterparts in ssn_viewer_html.py. This removes the build-time
#     --max_merge_events cap: the selection is now made where it is displayed.
#     A REMOVAL, not an addition, so a v3/v4/v5 reader would find keys missing from a
#     v6 file. This build reads all of them, because everything v6 needs has been in
#     the format since v3.
SSN_VIEWER_BUNDLE_VERSION = 6
# Versions this build can read. Kept as a tuple rather than an equality check so
# that additive revisions do not strand previously written bundles. Pre-v6 files
# carry the derived keys as well; they are ignored in favor of recomputing, so that
# a capped series written into an old file cannot silently limit what is shown now.
SUPPORTED_SSN_VIEWER_BUNDLE_VERSIONS = (3, 4, 5, 6)


def load_bundle(path: Union[str, PathLike]) -> dict:
    """Load and validate a ``.dsnv`` bundle (gzip JSON, or plain JSON as a fallback)."""
    with open(path, "rb") as handle:
        raw = handle.read()
    if raw[:2] == b"\x1f\x8b":  # gzip magic number
        raw = gzip.decompress(raw)
    bundle = json.loads(raw.decode("utf-8"))
    fmt = bundle.get("format")
    if fmt != SSN_VIEWER_BUNDLE_FORMAT:
        raise ValueError(
            f"'{path}' is not a Domainator Similarity Network Viewer bundle (format='{fmt}', "
            f"expected '{SSN_VIEWER_BUNDLE_FORMAT}')."
        )
    version = bundle.get("version")
    if version not in SUPPORTED_SSN_VIEWER_BUNDLE_VERSIONS:
        supported = ", ".join(str(v) for v in SUPPORTED_SSN_VIEWER_BUNDLE_VERSIONS)
        raise ValueError(
            f"Unsupported Domainator Similarity Network Viewer bundle version {version}; "
            f"this build understands version(s) {supported}."
        )
    return bundle


class _BundleTree:
    """The slice of the :class:`~domainator.ssn_edges.MaxTree` surface that the
    ``ssn_hierarchy`` functions actually touch, backed by a bundle's ``graph``.

    They need only ``n_nodes`` and ``mst_edges``; everything else on a real ``MaxTree``
    is derived from the input matrix, which a bundle does not carry.
    """

    __slots__ = ("n_nodes", "mst_edges")

    def __init__(self, n_nodes: int, mst_edges: List[Sequence]):
        self.n_nodes = n_nodes
        self.mst_edges = mst_edges


def bundle_tree(bundle: dict) -> _BundleTree:
    """Adapt a loaded bundle to what the ``ssn_hierarchy`` functions expect."""
    graph = bundle["graph"]
    return _BundleTree(len(graph["nodes"]), graph["mst_edges"])


def merge_event_rows(bundle: dict, merge_impact_metric: Optional[str] = None) -> List[dict]:
    """The bundle's **complete, uncapped** split-event series.

    Recomputed from ``graph.mst_edges`` rather than read from the file. Before v6 a
    capped copy was stored under ``graph.merge_event_series``; deriving it here means
    a caller is never silently limited to whatever cap the file happened to be built
    with, and pre-v6 bundles get the full series too.

    The metric defaults to the one the bundle records. Pass ``merge_impact_metric``
    to recompute under the other one -- the series is a pure function of the merge
    order and the metric, so nothing about the file constrains the choice.
    """
    if merge_impact_metric is None:
        merge_impact_metric = bundle["graph"].get("merge_impact_metric", MERGE_IMPACT_MIN_CHILD)
    summary = component_size_summary_by_threshold(
        bundle_tree(bundle), merge_impact_metric=merge_impact_metric
    )
    return threshold_merge_event_rows(summary)


def merge_event_series(bundle: dict, max_merge_events: int = DEFAULT_MAX_MERGE_EVENTS,
                       merge_impact_metric: Optional[str] = None,
                       event_rows: Optional[List[dict]] = None) -> List[dict]:
    """The capped selection of split events, as the viewer's split chart plots it.

    ``max_merge_events=0`` means no cap. Pass ``event_rows`` from
    :func:`merge_event_rows` to avoid replaying the merge order twice.
    """
    if event_rows is None:
        event_rows = merge_event_rows(bundle, merge_impact_metric=merge_impact_metric)
    return filter_merge_event_rows(event_rows, max_merge_events=max_merge_events)


def moving_sum(bundle: dict, merge_impact_metric: Optional[str] = None,
               event_rows: Optional[List[dict]] = None) -> dict:
    """The centered moving sum of split impact, over the **uncapped** rows."""
    if event_rows is None:
        event_rows = merge_event_rows(bundle, merge_impact_metric=merge_impact_metric)
    return merge_event_moving_sum(event_rows)


def slider_stops(bundle: dict, max_merge_events: int = DEFAULT_MAX_MERGE_EVENTS,
                 merge_impact_metric: Optional[str] = None,
                 event_rows: Optional[List[dict]] = None) -> List[dict]:
    """The threshold cuts a slider over this bundle should offer.

    ``∞``, one stop per selected split event, then the floor stop that reaches the
    fully merged network. Stops carry no ``threshold_index``: that indexed threshold
    tables which described the full input graph, and no bundle ever carried those.
    """
    if event_rows is None:
        event_rows = merge_event_rows(bundle, merge_impact_metric=merge_impact_metric)
    selected = filter_merge_event_rows(event_rows, max_merge_events=max_merge_events)
    return threshold_slider_stops(selected, tree=bundle_tree(bundle))


def coarsest_threshold(hierarchy: dict) -> float:
    """The cut below which nothing splits: one cluster per connected component."""
    merge_thresholds = [
        node["threshold"] for node in hierarchy["nodes"] if node["kind"] == "cluster"
    ]
    if len(merge_thresholds) == 0:
        return 0.0
    return floor_threshold_value(min(merge_thresholds), max(merge_thresholds))


def clusters_at_threshold(hierarchy: dict, threshold: Optional[float]) -> List[int]:
    """Return the component ids of the clusters active at ``threshold``.

    Python port of the viewer's ``activeClustersAtThreshold``: descend into a
    cluster node when it merged at or below the cut, otherwise emit the whole
    component. A merge at *exactly* the cut is split back apart, because
    ``build_ssn --lb T`` keeps only edges scoring strictly above ``T`` -- the
    comparison here has to be ``<=`` to match that and the viewer.
    ``threshold=None`` is the bundle's "∞" slider stop (``threshold_value: null``)
    and, like the viewer, means +infinity: every merge is cut, so every node is
    its own cluster. Use :func:`floor_threshold_value` for the opposite end.
    Higher thresholds yield finer partitions (more splitting), matching
    similarity-score semantics.
    """
    nodes = hierarchy["nodes"]
    # null is how slider_stops spells the ∞ stop; the viewer maps it to Infinity.
    cut = float("inf") if threshold is None else float(threshold)
    active: List[int] = []
    stack = list(reversed(hierarchy["roots"]))
    while stack:
        component_id = stack.pop()
        component = nodes[component_id]
        if component["kind"] == "leaf":
            active.append(component_id)
            continue
        if component["threshold"] <= cut:
            stack.append(component["right"])
            stack.append(component["left"])
            continue
        active.append(component_id)
    return active


def component_members(hierarchy: dict, component_id: int) -> List[int]:
    """Return the node indices belonging to ``component_id`` (a leaf_order slice)."""
    component = hierarchy["nodes"][component_id]
    start = component["leaf_start"]
    count = component["leaf_count"]
    return list(hierarchy["leaf_order"][start:start + count])


def node_index_by_name(bundle: dict) -> Dict[str, int]:
    """Map node id string -> position in ``graph.nodes`` (== metadata row index)."""
    return {name: idx for idx, name in enumerate(bundle["graph"]["nodes"])}


def summarize_cluster_metadata(
    member_indices: List[int],
    metadata: dict,
    top_n: int = 50,
) -> List[dict]:
    """Summarize the metadata distribution across a set of member node indices.

    For string columns: the ``top_n`` most common values with counts, plus the
    number of distinct values and missing entries. For numeric columns:
    count/missing and min/max/mean/quartiles. Reads the positional
    ``metadata['columns']`` / ``metadata['rows']`` table.
    """
    columns = metadata.get("columns", [])
    rows = metadata.get("rows", [])
    summaries = []
    for col_idx, column in enumerate(columns):
        values = [rows[i][col_idx] for i in member_indices]
        present = [v for v in values if v is not None]
        n_missing = len(values) - len(present)
        summary = {
            "name": column["name"],
            "type": column["type"],
            "count": len(present),
            "missing": n_missing,
        }
        if column["type"] in ("int", "float") and present:
            arr = np.asarray(present, dtype=float)
            summary.update({
                "min": float(np.min(arr)),
                "max": float(np.max(arr)),
                "mean": float(np.mean(arr)),
                "p25": float(np.percentile(arr, 25)),
                "median": float(np.percentile(arr, 50)),
                "p75": float(np.percentile(arr, 75)),
            })
        else:
            counts: Dict[str, int] = {}
            for value in present:
                key = str(value)
                counts[key] = counts.get(key, 0) + 1
            ranked = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
            summary["unique"] = len(counts)
            summary["top"] = [
                {"value": value, "count": count}
                for value, count in ranked[:top_n]
            ]
            if len(ranked) > top_n:
                summary["truncated"] = len(ranked) - top_n
        summaries.append(summary)
    return summaries
