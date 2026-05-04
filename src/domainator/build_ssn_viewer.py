"""Build a compact MST-based bundle for the standalone SSN viewer."""

import gzip
import json
import os
import sys
from os import PathLike
from pathlib import Path
from typing import Dict, List, Union

import pandas as pd
from jsonargparse import ActionConfigFile, ArgumentParser

from domainator import __version__, RawAndDefaultsFormatter
from domainator.build_ssn import subset_matrix_by_labels
from domainator.data_matrix import DataMatrix, MaxTree
from domainator.output_guardrails import (
    OutputSizeLimitExceeded,
    add_max_output_gb_argument,
    enforce_output_limit,
    make_temporary_output_path,
    max_output_gb_to_bytes,
)
from domainator.ssn_hierarchy import (
    DEFAULT_MAX_MERGE_EVENTS,
    MERGE_IMPACT_CHOICES,
    MERGE_IMPACT_MIN_CHILD,
    build_mst_component_hierarchy,
    component_size_summary_by_threshold,
    filter_merge_event_rows,
    threshold_merge_event_rows,
)
from domainator.ssn_viewer_html import write_ssn_viewer_html
from domainator.utils import list_and_file_to_dict_keys


SSN_VIEWER_BUNDLE_FORMAT = "domainator_ssn_viewer_bundle"
SSN_VIEWER_BUNDLE_VERSION = 1


def _infer_metadata_type(series: pd.Series) -> str:
    values = [value for value in series.tolist() if pd.notna(value)]
    if len(values) == 0:
        return "str"

    numeric_values = []
    for value in values:
        if isinstance(value, bool):
            return "str"
        if isinstance(value, (int, float)):
            numeric_values.append(float(value))
            continue
        return "str"

    if all(value.is_integer() for value in numeric_values):
        return "int"
    return "float"


def _normalize_metadata_value(value, column_type: str):
    if pd.isna(value):
        return None
    if column_type == "int":
        return int(value)
    if column_type == "float":
        return float(value)
    return str(value)


def _metadata_payload(node_data: pd.DataFrame) -> dict:
    columns = []
    rows = []
    for column_name in node_data.columns:
        columns.append({
            "name": column_name,
            "type": _infer_metadata_type(node_data[column_name]),
        })

    for _, row in node_data.iterrows():
        rows.append([
            _normalize_metadata_value(row[column["name"]], column["type"])
            for column in columns
        ])

    return {
        "columns": columns,
        "rows": rows,
    }


def _json_ready(value):
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    if hasattr(value, "tolist"):
        return _json_ready(value.tolist())
    if isinstance(value, float):
        if value == float("inf") or value == float("-inf"):
            return None
        return value
    return value


def _slider_stops(merge_event_rows):
    stops = [{
        "edge_index": -1,
        "threshold_label": "∞",
        "threshold_value": None,
    }]
    for merge_row in merge_event_rows:
        stops.append({
            "edge_index": int(merge_row["edge_index"]),
            "threshold_label": merge_row["threshold_to"],
            "threshold_value": float(merge_row["threshold_value"]),
        })
    return stops


def build_ssn_viewer_bundle(
    matrix: DataMatrix,
    metadata_files: List[Union[str, PathLike]] = None,
    subset_labels=None,
    merge_impact_metric: str = MERGE_IMPACT_MIN_CHILD,
    max_merge_events: int = DEFAULT_MAX_MERGE_EVENTS,
    color_by: str = None,
    label_by: str = None,
    name: str = None,
):
    if merge_impact_metric not in MERGE_IMPACT_CHOICES:
        raise ValueError(f"merge_impact_metric must be one of {sorted(MERGE_IMPACT_CHOICES)}")
    if max_merge_events < 0:
        raise ValueError("max_merge_events must be >= 0")

    matrix = subset_matrix_by_labels(matrix, subset_labels)
    if not matrix.symmetric_labels:
        raise ValueError("Input does not have symmetric axis labels. Can only build an SSN viewer bundle from a symmetric matrix.")

    node_data = pd.DataFrame(index=matrix.rows)
    if metadata_files is not None:
        for file in metadata_files:
            metadata = pd.read_csv(file, sep="\t", index_col=0)
            node_data = node_data.merge(metadata, how="left", left_index=True, right_index=True)

    if color_by is not None and color_by not in node_data.columns:
        raise ValueError(f"Requested color_by column '{color_by}' was not found in the merged metadata.")
    if label_by is not None and label_by not in node_data.columns:
        raise ValueError(f"Requested label_by column '{label_by}' was not found in the merged metadata.")

    tree = MaxTree(matrix)
    component_summary = component_size_summary_by_threshold(tree, merge_impact_metric=merge_impact_metric)
    merge_event_series = filter_merge_event_rows(
        threshold_merge_event_rows(component_summary),
        max_merge_events=max_merge_events,
    )
    hierarchy = build_mst_component_hierarchy(tree)

    bundle = {
        "format": SSN_VIEWER_BUNDLE_FORMAT,
        "version": SSN_VIEWER_BUNDLE_VERSION,
        "name": name,
        "domainator_version": __version__,
        "graph": {
            "nodes": list(matrix.rows),
            "mst_edges": tree.export_for_interactive_viz()["mst_edges"],
            "cluster_count_by_threshold": tree.cluster_count_by_threshold,
            "edges_by_threshold": tree.edges_by_threshold,
            "merge_impact_metric": merge_impact_metric,
            "merge_event_series": merge_event_series,
            "slider_stops": _slider_stops(merge_event_series),
            "hierarchy": hierarchy,
        },
        "metadata": _metadata_payload(node_data),
        "defaults": {
            "color_by": color_by,
            "label_by": label_by,
        },
    }
    return _json_ready(bundle)


def write_ssn_viewer_bundle(
    out_path: Union[str, PathLike],
    bundle: dict,
    max_output_bytes: int | None = None,
):
    json_bytes = json.dumps(bundle, separators=(",", ":"), sort_keys=True).encode("utf-8")
    compressed_bytes = gzip.compress(json_bytes, mtime=0)
    enforce_output_limit(
        projected_bytes=len(compressed_bytes),
        max_output_bytes=max_output_bytes,
        output_description=f"SSN viewer bundle output '{out_path}'",
        mitigation_options=["--subset", "--subset_file"],
        extra_guidance="This bundle stores only MST-derived hierarchy data; subset the network before export if it is still too large.",
    )

    temp_path = make_temporary_output_path(out_path)
    try:
        with open(temp_path, "wb") as out_handle:
            out_handle.write(compressed_bytes)
        os.replace(temp_path, out_path)
        temp_path = None
    finally:
        if temp_path is not None and Path(temp_path).exists():
            os.unlink(temp_path)


def main(argv):
    parser = ArgumentParser(description=f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Symmetric similarity matrix to convert into an SSN viewer bundle. Format can be tab-separated text or Domainator hdf5.")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Path to write the gzip-compressed SSN viewer bundle.")
    parser.add_argument("--viewer_html", type=str, default=None,
                        help="Optional path to write a standalone static HTML viewer shell that loads local bundle files.")
    parser.add_argument("--name", type=str, default=None,
                        help="Optional display name recorded in the bundle. Defaults to the output file stem.")
    parser.add_argument("--metadata", type=str, nargs="+", required=False, default=None,
                        help="Tab-separated metadata tables keyed by node id.")
    parser.add_argument("--subset", type=str, default=None, nargs='+',
                        help="Only consider matrix labels in this list. Additive with --subset_file.")
    parser.add_argument("--subset_file", type=str, default=None,
                        help="Text file containing matrix labels to retain.")
    parser.add_argument("--color_by", type=str, default=None,
                        help="Record a default metadata field for node coloring in the viewer bundle.")
    parser.add_argument("--label_by", type=str, default=None,
                        help="Record a default metadata field for node labels in the viewer bundle.")
    parser.add_argument("--merge_impact_metric", choices=list(MERGE_IMPACT_CHOICES), default=MERGE_IMPACT_MIN_CHILD,
                        help="Metric recorded for split events in the bundle.")
    parser.add_argument("--max_merge_events", type=int, default=DEFAULT_MAX_MERGE_EVENTS,
                        help="Maximum number of strongest merge events to embed in the viewer bundle threshold slider and split plot. Use 0 to include all merge events.")
    add_max_output_gb_argument(parser)
    parser.add_argument("--config", action=ActionConfigFile)

    params = parser.parse_args(argv)
    max_output_bytes = max_output_gb_to_bytes(params.max_output_gb)
    subset_labels = list_and_file_to_dict_keys(params.subset, params.subset_file, as_set=True)
    bundle_name = params.name if params.name is not None else Path(params.output).stem

    try:
        bundle = build_ssn_viewer_bundle(
            DataMatrix.from_file(params.input),
            metadata_files=params.metadata,
            subset_labels=subset_labels,
            merge_impact_metric=params.merge_impact_metric,
            max_merge_events=params.max_merge_events,
            color_by=params.color_by,
            label_by=params.label_by,
            name=bundle_name,
        )
        write_ssn_viewer_bundle(params.output, bundle, max_output_bytes=max_output_bytes)
        if params.viewer_html is not None:
            write_ssn_viewer_html(params.viewer_html, title=bundle_name)
    except OutputSizeLimitExceeded as exc:
        raise SystemExit(str(exc)) from None


def _entrypoint():
    main(sys.argv[1:])


if __name__ == '__main__':
    main(sys.argv[1:])