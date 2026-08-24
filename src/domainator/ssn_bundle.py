"""Read-only helpers for navigating a build_ssn_viewer ``.ssnv`` bundle.

The bundle is a gzip-compressed JSON document (see ``build_ssn_viewer.py``) that
stores an agglomerative MST merge tree plus a positional metadata table. These
helpers let non-browser callers (notably ``ssn_navigator.py``) reconstruct the
clusters that are active at a chosen similarity threshold and summarize the
metadata of their members, mirroring the logic of the JavaScript viewer's
``activeClustersAtThreshold`` / ``componentMembers`` functions.
"""

import gzip
import json
from os import PathLike
from typing import Dict, List, Optional, Union

import numpy as np

from domainator.build_ssn_viewer import (
    SSN_VIEWER_BUNDLE_FORMAT,
    SSN_VIEWER_BUNDLE_VERSION,
)


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
    if version != SSN_VIEWER_BUNDLE_VERSION:
        raise ValueError(
            f"Unsupported SSN viewer bundle version {version}; "
            f"this build understands version {SSN_VIEWER_BUNDLE_VERSION}."
        )
    return bundle


def clusters_at_threshold(hierarchy: dict, threshold: Optional[float]) -> List[int]:
    """Return the component ids of the clusters active at ``threshold``.

    Python port of the viewer's ``activeClustersAtThreshold``: descend into a
    cluster node only when it merged at a weaker (lower) threshold than the cut,
    otherwise emit the whole component. ``threshold=None`` yields the coarsest
    partition (one cluster per connected component). Higher thresholds yield
    finer partitions (more splitting), matching similarity-score semantics.
    """
    nodes = hierarchy["nodes"]
    # None (the "∞" slider stop) coerces to 0 to match the JS `< null` behavior.
    cut = 0.0 if threshold is None else float(threshold)
    active: List[int] = []
    stack = list(reversed(hierarchy["roots"]))
    while stack:
        component_id = stack.pop()
        component = nodes[component_id]
        if component["kind"] == "leaf":
            active.append(component_id)
            continue
        if component["threshold"] < cut:
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
