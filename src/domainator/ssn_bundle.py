"""Read-only helpers for navigating a build_ssn_viewer ``.ssnv`` bundle.

The bundle is a gzip-compressed JSON document (see ``build_ssn_viewer.py``) that
stores an agglomerative MST merge tree plus a positional metadata table. These
helpers let non-browser callers (notably ``ssn_navigator.py``) reconstruct the
clusters that are active at a chosen similarity threshold and summarize the
metadata of their members, mirroring the logic of the JavaScript viewer's
``activeClustersAtThreshold`` / ``componentMembers`` functions.

A v4 bundle may also carry a top-level ``app_state`` section, written by the HTML
viewer's "Save session" button. It holds viewer UI state only and is ignored
here; the rest of such a file is an ordinary bundle, so a saved session's
``metadata`` simply reflects whatever was edited in the viewer.
"""

import gzip
import json
from os import PathLike
from typing import Dict, List, Optional, Union

import numpy as np

from domainator.ssn_hierarchy import floor_threshold_value


SSN_VIEWER_BUNDLE_FORMAT = "domainator_ssn_viewer_bundle"
# v3: per-event merge_size_counts/largest_merge/merge_count + graph.merge_moving_sum
# v4: optional top-level "app_state" section written by the HTML viewer's "Save
#     session" button (viewer UI state + selection presets). Purely additive:
#     build_ssn_viewer.py never writes app_state, and a reader that ignores the
#     section can treat a v4 file exactly like a v3 file.
SSN_VIEWER_BUNDLE_VERSION = 4
# Versions this build can read. Kept as a tuple rather than an equality check so
# that additive revisions do not strand previously written bundles.
SUPPORTED_SSN_VIEWER_BUNDLE_VERSIONS = (3, 4)


def load_bundle(path: Union[str, PathLike]) -> dict:
    """Load and validate a ``.ssnv`` bundle (gzip JSON, or plain JSON as a fallback)."""
    with open(path, "rb") as handle:
        raw = handle.read()
    if raw[:2] == b"\x1f\x8b":  # gzip magic number
        raw = gzip.decompress(raw)
    bundle = json.loads(raw.decode("utf-8"))
    fmt = bundle.get("format")
    if fmt != SSN_VIEWER_BUNDLE_FORMAT:
        raise ValueError(
            f"'{path}' is not an SSN viewer bundle (format='{fmt}', "
            f"expected '{SSN_VIEWER_BUNDLE_FORMAT}')."
        )
    version = bundle.get("version")
    if version not in SUPPORTED_SSN_VIEWER_BUNDLE_VERSIONS:
        supported = ", ".join(str(v) for v in SUPPORTED_SSN_VIEWER_BUNDLE_VERSIONS)
        raise ValueError(
            f"Unsupported SSN viewer bundle version {version}; "
            f"this build understands version(s) {supported}."
        )
    return bundle


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
