import warnings
warnings.filterwarnings("ignore", module='numpy')
import numpy as np
import html
import seaborn as sns
import pandas as pd
from pandas.api.types import is_integer_dtype, is_float_dtype
import xml.etree.ElementTree as ET
from typing import Dict

from domainator.output_guardrails import build_output_limit_message, OutputSizeLimitExceeded

MAX_DOUBLE = np.finfo(np.float64).max
MIN_DOUBLE = np.finfo(np.float64).min
MAX_INT = np.iinfo(np.int64).max
MIN_INT = np.iinfo(np.int64).min

def indent(elem, level=0):
    """
        indent xml to make it more human-readable
    """
    # from: https://stackoverflow.com/questions/3095434/inserting-newlines-in-xml-file-generated-via-xml-etree-elementtree-in-python
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

def get_column_types(df):
    types = df.dtypes
    out = dict()
    for column, datatype in types.items():
        # Handle pandas extension dtypes like StringDtype that numpy cannot interpret.
        if is_integer_dtype(datatype):
            out[column] = 'integer'
        elif is_float_dtype(datatype):
            out[column] = 'real'
        else:
            out[column] = 'string'

    return out

def check_null(value, xgmml_type):
    """
        value: a value from a pandas dataframe
        xgmml_type: "string", "integer", "real"

        if value is not null, then return the value, otherwise return a default based on the type.
        string: ''
        integer: 0
        real: 0.0
    """
    defaults = {"string":'', "integer":0, "real":0.0}
    if pd.isna(value):
        return defaults[xgmml_type]
    if xgmml_type == "real":
        if value > MAX_DOUBLE:
            return MAX_DOUBLE
        if value < MIN_DOUBLE:
            return MIN_DOUBLE
    if xgmml_type == "integer":
        if value > MAX_INT:
            return MAX_INT
        if value < MIN_INT:
            return MIN_INT
    return value
        


#TODO: option to add edges for non-zero mutual-blast (or some threshold)
def _serialize_xgmml_element(element, tail):
    indent(element, level=1)
    element.tail = tail
    return ET.tostring(element, encoding="utf-8")


def _iter_cytoscape_xgmml_chunks(nodes, edges=None, name="", x_col=None, y_col=None, z_col=None, color_by=None, color_table:Dict[str,str]=None, shape="ELLIPSE", shape_x=25, shape_y=25):
    graph_name = name
    possible_shapes = {"DIAMOND", "ELLIPSE", "HEXAGON", "OCTAGON", "PARALLELOGRAM", "RECTANGLE", "ROUND_RECTANGLE", "TRIANGLE", "V"}
    if shape not in possible_shapes:
        raise ValueError(f"Unrecognized node shape specification {shape}, please choose one of: {possible_shapes}")

    palette = None
    if color_by is not None:
        if color_by not in nodes.columns:
            raise ValueError(f"color_by column {color_by} not found in nodes dataframe")
        if color_table is not None:
            palette = color_table
        else:
            raise ValueError("color_table must be provided if color_by is provided")

    metadata_types = get_column_types(nodes)
    node_names_to_ids = dict() if edges is not None else None
    col_to_idx = {v:i+1 for i,v in enumerate(nodes.columns)}
    node_count = len(nodes.index)
    edge_count = 0 if edges is None else len(edges.index)

    has_children = node_count + edge_count > 0

    yield b"<?xml version='1.0' encoding='utf-8'?>\n"
    yield f'<graph label="{html.escape(str(graph_name), quote=True)}" xmlns="http://www.cs.rpi.edu/XGMML">\n'.encode("utf-8")
    if has_children:
        yield b"  "

    for i, row in enumerate(nodes.itertuples(name=None)):
        name = str(row[0])
        node_id = f"n{i}"
        node = ET.Element("node", {"id":node_id, "label":node_id})
        if node_names_to_ids is not None:
            node_names_to_ids[name] = node_id

        ET.SubElement(node, "att", {"type":"string", "name":"Description", "value":html.escape(name)})
        for key, idx in col_to_idx.items():
            value = row[idx]
            if key != x_col and key != y_col and key != z_col:
                ET.SubElement(node, "att", {"type":metadata_types[key], "name":str(key), "value":html.escape(str(check_null(value, metadata_types[key])))})

        color = None
        if color_by is not None:
            if row[col_to_idx[color_by]] in palette:
                color = palette[row[col_to_idx[color_by]]]
            elif str(row[col_to_idx[color_by]]) in palette:
                color = palette[str(row[col_to_idx[color_by]])]
            else:
                color = palette[None]

        elem_dict = {"type":shape, "x":"0.0", "y":"0.0", "w":str(shape_x), "h":str(shape_y), "z":"0.0"}
        if color is not None:
            elem_dict["fill"] = color
        if x_col is not None and x_col in col_to_idx:
            elem_dict["x"] = str(row[col_to_idx[x_col]])
        if y_col is not None and y_col in col_to_idx:
            elem_dict["y"] = str(-1 * row[col_to_idx[y_col]])
        if z_col is not None and z_col in col_to_idx:
            elem_dict["z"] = str(row[col_to_idx[z_col]])
        ET.SubElement(node, "graphics", elem_dict)

        tail = "\n" if i == node_count - 1 and edge_count == 0 else "\n  "
        yield _serialize_xgmml_element(node, tail)

    if edges is not None:
        edge_metadata_types = get_column_types(edges)
        for i, row in enumerate(edges.itertuples()):
            source = node_names_to_ids[row[0][0]]
            target = node_names_to_ids[row[0][1]]
            edge = ET.Element("edge", {"id":f"e{i}", "label":f"e{i}", "source": source, "target": target})
            for key in row._fields[1:]:
                value = getattr(row, key)
                ET.SubElement(edge, "att", {"type":edge_metadata_types[key], "name":str(key), "value":str(check_null(value, edge_metadata_types[key]))})

            tail = "\n" if i == edge_count - 1 else "\n  "
            yield _serialize_xgmml_element(edge, tail)

    yield b"</graph>\n"


def estimate_cytoscape_xgmml_size(nodes, edges=None, name="", x_col=None, y_col=None, z_col=None, color_by=None, color_table:Dict[str,str]=None, shape="ELLIPSE", shape_x=25, shape_y=25):
    """Return the exact serialized XGMML size in bytes for the given node and edge tables."""
    return sum(
        len(chunk)
        for chunk in _iter_cytoscape_xgmml_chunks(
            nodes,
            edges=edges,
            name=name,
            x_col=x_col,
            y_col=y_col,
            z_col=z_col,
            color_by=color_by,
            color_table=color_table,
            shape=shape,
            shape_x=shape_x,
            shape_y=shape_y,
        )
    )


def write_cytoscape_xgmml(out_path, nodes, edges=None, name="", x_col=None, y_col=None, z_col=None, color_by=None,  color_table:Dict[str,str]=None, shape="ELLIPSE", shape_x=25, shape_y=25, max_output_bytes=None, output_description=None, mitigation_options=None, extra_guidance=None):
    """
        out_path: network file to write, should end in .xgmml
        nodes: pandas dataframe for nodes and their annotations. Should have a single index column which will be the node name.
        edges: pandas dataframe for edges and their annotations. Should have a double index, which will be the names of the start and end nodes.
        name: name of the graph
        x_col: which column in the nodes data frame to use as the x-coordinate for the node
        y_col: which column in the nodes data frame to use as the y-coordinate for the node
        color_by: give a different color to each value in this column
        color_table: dictionary mapping values in the color_by column to colors.
        shape: what shape the nodes should be
        shape_x: size of the nodes in the x-dimension
        shape_y: size of the nodes in the y-dimension
    """
    should_close = False
    if hasattr(out_path, "write"):
        out_handle = out_path
    else:
        out_handle = open(out_path, "wb")
        should_close = True

    bytes_written = 0
    try:
        for chunk in _iter_cytoscape_xgmml_chunks(
            nodes,
            edges=edges,
            name=name,
            x_col=x_col,
            y_col=y_col,
            z_col=z_col,
            color_by=color_by,
            color_table=color_table,
            shape=shape,
            shape_x=shape_x,
            shape_y=shape_y,
        ):
            out_handle.write(chunk)
            bytes_written += len(chunk)
            if max_output_bytes is not None and bytes_written > max_output_bytes:
                raise OutputSizeLimitExceeded(
                    build_output_limit_message(
                        output_description=output_description or "XGMML output",
                        projected_bytes=bytes_written,
                        max_output_bytes=max_output_bytes,
                        mitigation_options=mitigation_options,
                        extra_guidance=extra_guidance,
                    )
                )
    finally:
        if should_close:
            out_handle.close()
    