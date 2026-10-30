"""The Domainator Viewer bundle (``.dvb``) format, and the parts of it every view shares.

A bundle is gzip-compressed JSON that the Domainator Viewer (viewer_html.py) opens:
``{format, version, name, domainator_version, nodes, metadata, defaults, views}``, where
``format`` is :data:`VIEWER_BUNDLE_FORMAT` and ``views`` maps a view name ("ssn", written by
build_ssn; "tree", by build_tree; "projection", by build_projection) to that view's data about
``nodes``. See the .dvb section of docs/file_formats.md. ``metadata`` is the node metadata as
:func:`metadata_payload` writes it, and ``defaults`` the viewer's starting settings as
:func:`viewer_defaults` writes them. A bundle may also carry an ``app_state`` section written
by the viewer's "Save session".

Networks were written to a format of their own before, ``.dsnv``; see ssn_bundle.py.
"""

import gzip
import json
import os
from concurrent.futures import ThreadPoolExecutor
from os import PathLike
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Union

import pandas as pd

from domainator import __version__
from domainator.colors import normalize_color_hex
from domainator.output_guardrails import enforce_output_limit, make_temporary_output_path


VIEWER_BUNDLE_FORMAT = "domainator_viewer_bundle"
VIEWER_BUNDLE_VERSION = 2
SUPPORTED_VIEWER_BUNDLE_VERSIONS = (1, 2)


def read_node_metadata(paths: Optional[Iterable[Union[str, PathLike]]], index: Sequence[str]) -> pd.DataFrame:
    """The rows of ``index`` in tab-separated metadata tables keyed by their first column.

    Columns from each table are joined in turn; a node a table does not list gets missing
    values there.
    """
    node_data = pd.DataFrame(index=list(index))
    if paths:
        with ThreadPoolExecutor() as exe:
            tables = list(exe.map(lambda path: pd.read_csv(path, sep="\t", index_col=0), paths))
        for table in tables:
            node_data = node_data.merge(table, how="left", left_index=True, right_index=True)
    return node_data


def infer_metadata_type(series: pd.Series) -> str:
    """"int", "float" or "str": how the viewer reads a metadata column."""
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


def metadata_payload(node_data: pd.DataFrame) -> dict:
    """A bundle's ``metadata``: ``{columns: [{name, type}], rows}``, one row per node in
    ``node_data``'s order."""
    columns = []
    rows = []
    for column_name in node_data.columns:
        columns.append({
            "name": column_name,
            "type": infer_metadata_type(node_data[column_name]),
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


def viewer_defaults(
    node_data: pd.DataFrame,
    *,
    color_by: Optional[str] = None,
    color_table: Optional[Dict[str, str]] = None,
    label_by: Optional[str] = None,
    categorical: Optional[List[str]] = None,
) -> dict:
    """A bundle's ``defaults``: the settings the viewer opens with.

    Args:
        node_data: the node metadata, to check the named columns against.
        color_by: the metadata column the viewer colors nodes by at first.
        color_table: maps ``color_by`` values to hex colors, the column's starting palette in
            the viewer. Values not in it are colored gray. Requires color_by, and a numeric
            color_by column is colored as categories.
        label_by: the metadata column the viewer labels nodes by at first.
        categorical: numeric metadata columns to color as categories rather than as a gradient.
    """
    if color_by is not None and color_by not in node_data.columns:
        raise ValueError(f"Requested color_by column '{color_by}' was not found in the merged metadata.")
    if label_by is not None and label_by not in node_data.columns:
        raise ValueError(f"Requested label_by column '{label_by}' was not found in the merged metadata.")
    categorical_columns = list(categorical) if categorical else []
    for column_name in categorical_columns:
        if column_name not in node_data.columns:
            raise ValueError(f"Requested categorical column '{column_name}' was not found in the merged metadata.")
    if color_table is not None:
        if color_by is None:
            raise ValueError("color_table requires color_by, the metadata column the table colors.")
        # Missing values are colored gray whatever the table says.
        color_table = {str(value): normalize_color_hex(color) for value, color in color_table.items() if value is not None}
        # A color table is a discrete mapping, so a numeric column it colors is a categorical one.
        if infer_metadata_type(node_data[color_by]) != "str" and color_by not in categorical_columns:
            categorical_columns.append(color_by)

    defaults = {
        "color_by": color_by,
        "label_by": label_by,
        # Numeric columns the viewer should color with discrete category colors
        # instead of a gradient (e.g. integer cluster numbers). Toggleable in the
        # viewer's color picker; this only sets the initial state.
        "categorical_columns": categorical_columns,
    }
    if color_table is not None:
        # value -> "#RRGGBB" for the color_by column; the viewer starts from it in place of
        # its default palette. Optional, so bundles built without one are unchanged.
        defaults["color_table"] = color_table
    return defaults


def json_ready(value):
    """``value`` with numpy arrays and scalars as Python values, and infinities as None."""
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if hasattr(value, "tolist"):
        return json_ready(value.tolist())
    if isinstance(value, float):
        if value == float("inf") or value == float("-inf"):
            return None
        return value
    return value


def viewer_bundle(name: Optional[str], nodes: Sequence[str], node_data: pd.DataFrame, defaults: dict, views: dict) -> dict:
    """A ``domainator_viewer_bundle``: ``views`` maps view names to their data about ``nodes``,
    whose metadata rows are ``node_data``'s, in the same order."""
    return json_ready({
        "format": VIEWER_BUNDLE_FORMAT,
        "version": VIEWER_BUNDLE_VERSION,
        "name": name,
        "domainator_version": __version__,
        "nodes": list(nodes),
        "metadata": metadata_payload(node_data),
        "defaults": defaults,
        "views": views,
    })


def serialize_bundle(bundle: dict) -> bytes:
    """The bundle as a bundle file holds it: compact JSON, gzip-compressed without a timestamp."""
    json_bytes = json.dumps(bundle, separators=(",", ":"), sort_keys=True).encode("utf-8")
    return gzip.compress(json_bytes, mtime=0)


def write_bundle_bytes(
    out_path: Union[str, PathLike],
    compressed_bytes: bytes,
    *,
    description: str,
    max_output_bytes: Optional[int] = None,
    mitigation_options: Sequence[str] = (),
    extra_guidance: Optional[str] = None,
) -> None:
    """Write serialize_bundle's bytes through a temporary file, failing first if they exceed
    ``max_output_bytes``."""
    enforce_output_limit(
        projected_bytes=len(compressed_bytes),
        max_output_bytes=max_output_bytes,
        output_description=f"{description} '{out_path}'",
        mitigation_options=list(mitigation_options),
        extra_guidance=extra_guidance,
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
