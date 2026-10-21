"""Tests for the XGMML writer."""
import xml.etree.ElementTree as ET

import pandas as pd

from domainator.cytoscape import _iter_cytoscape_xgmml_chunks

_NS = {"x": "http://www.cs.rpi.edu/XGMML"}


def _parse(nodes, **kwargs):
    return ET.fromstring(b"".join(_iter_cytoscape_xgmml_chunks(nodes, **kwargs)))


def test_values_with_markup_characters_round_trip():
    nodes = pd.DataFrame({"note": ['x<y & "z"']}, index=["A&B"])
    graph = _parse(nodes, name="g&h")
    assert graph.get("label") == "g&h"
    atts = {a.get("name"): a.get("value") for a in graph.find("x:node", _NS).findall("x:att", _NS)}
    assert atts == {"Description": "A&B", "note": 'x<y & "z"'}


def test_attribute_types_and_null_defaults():
    nodes = pd.DataFrame({"count": pd.array([3, None], dtype="Int64"), "score": [0.5, None],
                          "flag": [True, False], "label": ["a", None]}, index=["n1", "n2"])
    graph = _parse(nodes)
    second = graph.findall("x:node", _NS)[1]
    atts = {a.get("name"): (a.get("type"), a.get("value")) for a in second.findall("x:att", _NS)}
    assert atts == {"Description": ("string", "n2"), "count": ("integer", "0"), "score": ("real", "0.0"),
                    "flag": ("boolean", "false"), "label": ("string", "")}


def test_edges_reference_node_ids():
    nodes = pd.DataFrame(index=["A", "B"])
    edges = pd.DataFrame({"score": [1.5]}, index=pd.MultiIndex.from_tuples([("A", "B")]))
    graph = _parse(nodes, edges=edges)
    ids = {n.find("x:att", _NS).get("value"): n.get("id") for n in graph.findall("x:node", _NS)}
    edge = graph.find("x:edge", _NS)
    assert (edge.get("source"), edge.get("target")) == (ids["A"], ids["B"])
    assert edge.find("x:att", _NS).attrib == {"type": "real", "name": "score", "value": "1.5"}
