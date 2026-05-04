import base64
import gzip
import json
import os
import tempfile

import numpy as np
import pandas as pd
import pytest

from domainator import build_ssn_viewer
from domainator.data_matrix import DenseDataMatrix, MaxTree


def _read_bundle(path):
    with open(path, "rb") as handle:
        return json.loads(gzip.decompress(handle.read()).decode("utf-8"))


def _write_metadata(path, row_names):
    metadata = pd.DataFrame(
        {
            "category": ["alpha", "alpha", "beta", "gamma"],
            "count": [1, 2, 3, 4],
            "score": [1.5, 2.5, 3.5, 4.5],
        },
        index=row_names,
    )
    metadata.to_csv(path, sep="\t")


def test_build_ssn_viewer_writes_bundle_with_metadata_defaults():
    data = np.array([
        [0, 10, 6, 0],
        [10, 0, 7, 0],
        [6, 7, 0, 4],
        [0, 0, 4, 0],
    ], dtype=float)
    row_names = ["A", "B", "C", "D"]
    matrix = DenseDataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        metadata_file = os.path.join(output_dir, "metadata.tsv")
        bundle_file = os.path.join(output_dir, "test_bundle.ssnv")

        matrix.write(input_file, output_type="dense")
        _write_metadata(metadata_file, row_names)

        build_ssn_viewer.main([
            "-i", input_file,
            "-o", bundle_file,
            "--metadata", metadata_file,
            "--color_by", "category",
            "--label_by", "category",
        ])

        bundle = _read_bundle(bundle_file)

        assert bundle["format"] == build_ssn_viewer.SSN_VIEWER_BUNDLE_FORMAT
        assert bundle["version"] == build_ssn_viewer.SSN_VIEWER_BUNDLE_VERSION
        assert bundle["graph"]["nodes"] == row_names
        assert len(bundle["graph"]["mst_edges"]) == 3
        assert len(bundle["graph"]["merge_event_series"]) == 3
        assert [stop["edge_index"] for stop in bundle["graph"]["slider_stops"]] == [-1, 0, 1, 2]
        assert bundle["graph"]["hierarchy"]["roots"] == [6]
        assert bundle["graph"]["hierarchy"]["leaf_order"] == [0, 1, 2, 3]
        assert bundle["graph"]["hierarchy"]["nodes"][6]["leaf_count"] == 4
        assert bundle["defaults"] == {"color_by": "category", "label_by": "category"}
        assert bundle["metadata"]["columns"] == [
            {"name": "category", "type": "str"},
            {"name": "count", "type": "int"},
            {"name": "score", "type": "float"},
        ]
        assert bundle["metadata"]["rows"][0] == ["alpha", 1, 1.5]


def test_build_ssn_viewer_cluster_counts_match_maxtree():
    data = np.array([
        [0, 10, 0, 0, 0],
        [10, 0, 5, 0, 0],
        [0, 5, 0, 4, 0],
        [0, 0, 4, 0, 1],
        [0, 0, 0, 1, 0],
    ], dtype=float)
    row_names = ["A", "B", "C", "D", "E"]
    matrix = DenseDataMatrix(data, row_names, row_names)
    tree = MaxTree(matrix)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        bundle_file = os.path.join(output_dir, "test_bundle.ssnv")
        matrix.write(input_file, output_type="dense")

        build_ssn_viewer.main([
            "-i", input_file,
            "-o", bundle_file,
        ])

        bundle = _read_bundle(bundle_file)

        expected_counts = [[None, int(tree.n_nodes)]]
        expected_counts.extend([
            [float(row[0]), int(row[1])]
            for row in tree.cluster_count_by_threshold[1:]
        ])

        assert bundle["graph"]["cluster_count_by_threshold"] == expected_counts
        assert bundle["graph"]["edges_by_threshold"] == [
            [int(row[0]), float(row[1])]
            for row in tree.edges_by_threshold
        ]
        assert bundle["graph"]["hierarchy"]["nodes"][8]["size"] == 5
        assert bundle["graph"]["hierarchy"]["nodes"][8]["leaf_count"] == 5


def test_build_ssn_viewer_limits_merge_events_and_slider_stops():
    data = np.array([
        [0, 10, 0, 0, 0],
        [10, 0, 5, 0, 0],
        [0, 5, 0, 4, 0],
        [0, 0, 4, 0, 1],
        [0, 0, 0, 1, 0],
    ], dtype=float)
    row_names = ["A", "B", "C", "D", "E"]
    matrix = DenseDataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        bundle_file = os.path.join(output_dir, "test_bundle.ssnv")
        matrix.write(input_file, output_type="dense")

        build_ssn_viewer.main([
            "-i", input_file,
            "-o", bundle_file,
            "--max_merge_events", "2",
        ])

        bundle = _read_bundle(bundle_file)

        assert len(bundle["graph"]["merge_event_series"]) == 2
        assert len(bundle["graph"]["slider_stops"]) == 3
        assert bundle["graph"]["slider_stops"][0]["threshold_value"] is None


def test_build_ssn_viewer_subset_filters_nodes_and_metadata():
    data = np.array([
        [0, 10, 6, 1],
        [10, 0, 7, 1],
        [6, 7, 0, 4],
        [1, 1, 4, 0],
    ], dtype=float)
    row_names = ["A", "B", "C", "D"]
    matrix = DenseDataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        metadata_file = os.path.join(output_dir, "metadata.tsv")
        bundle_file = os.path.join(output_dir, "test_bundle.ssnv")

        matrix.write(input_file, output_type="dense")
        _write_metadata(metadata_file, row_names)

        build_ssn_viewer.main([
            "-i", input_file,
            "-o", bundle_file,
            "--metadata", metadata_file,
            "--subset", "A", "C", "D",
        ])

        bundle = _read_bundle(bundle_file)

        assert bundle["graph"]["nodes"] == ["A", "C", "D"]
        assert len(bundle["metadata"]["rows"]) == 3
        assert bundle["metadata"]["rows"][0][0] == "alpha"
        assert bundle["metadata"]["rows"][1][0] == "beta"
        assert bundle["metadata"]["rows"][2][0] == "gamma"


def test_build_ssn_viewer_writes_static_html_shell():
    data = np.array([
        [0, 10, 6, 0],
        [10, 0, 7, 0],
        [6, 7, 0, 4],
        [0, 0, 4, 0],
    ], dtype=float)
    row_names = ["A", "B", "C", "D"]
    matrix = DenseDataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        bundle_file = os.path.join(output_dir, "test_bundle.ssnv")
        html_file = os.path.join(output_dir, "viewer.html")
        matrix.write(input_file, output_type="dense")

        build_ssn_viewer.main([
            "-i", input_file,
            "-o", bundle_file,
            "--viewer_html", html_file,
            "--name", "Viewer Test",
        ])

        html_content = open(html_file, "r", encoding="utf-8").read()

        assert '<input id="bundle-file" type="file"' in html_content
        assert 'DecompressionStream' in html_content
        assert 'Cluster Splits vs Threshold' in html_content
        assert 'Hierarchy View' in html_content
        assert 'View Settings' in html_content
        assert 'Node Metadata' in html_content
        assert 'Export table TSV' in html_content
        assert 'threshold-min-label' in html_content
        assert 'threshold-max-label' in html_content
        assert 'threshold-input' in html_content
        assert 'Jump to threshold' in html_content
        assert 'Split impact' in html_content
        assert 'Threshold' in html_content
        assert 'componentMembers' in html_content
        assert 'mstLinksForActiveClusters' in html_content
        assert 'tidyForestLayout' in html_content
        assert 'forceDirectedForestLayout' in html_content
        assert 'organicForestLayout' in html_content
        assert 'gridClusterLayout' in html_content
        assert 'computeVisibleLayout' in html_content
        assert 'renderClusterView' in html_content
        assert 'refineLayoutGeometry' in html_content
        assert 'pointSegmentDistance' in html_content
        assert 'trimmedLinkEndpoints' in html_content
        assert 'renderedLinkSegments' in html_content
        assert 'sort-components-by-size' in html_content
        assert 'Sort components by size' in html_content
        assert 'sortComponentsBySizeEnabled' in html_content
        assert 'componentSelectionState' in html_content
        assert 'toggleSelectionForNode' in html_content
        assert 'hitTestNodeAt' in html_content
        assert 'metadata-sort-button' in html_content
        assert 'toggleMetadataSort' in html_content
        assert 'sortedMetadataNodeIndices' in html_content
        assert 'metadata-filter' in html_content
        assert 'metadata-select-nodes' in html_content
        assert 'metadata-reset-sort' in html_content
        assert 'metadata-null-order' in html_content
        assert 'filteredMetadataNodeIndices' in html_content
        assert 'formatMetadataDisplayValue' in html_content
        assert 'toggleMetadataRowSelection' in html_content
        assert 'metadataRowSelectionAnchor' in html_content
        assert 'selectNodesFromMetadataRows' in html_content
        assert 'layout-algorithm' in html_content
        assert 'Layout algorithm' in html_content
        assert '<option value="tree">Tree</option>' in html_content
        assert '<option value="force">Force-directed</option>' in html_content
        assert '<option value="organic">Organic</option>' in html_content
        assert '<option value="grid" selected>Grid (no edges)</option>' in html_content
        assert 'leaf-pruning-only' in html_content
        assert 'Minimum cluster size trims leaf clusters only' in html_content
        assert 'show-node-counts' in html_content
        assert 'Show node count labels' in html_content
        assert 'show-edge-scores' in html_content
        assert 'Show edge score labels' in html_content
        assert '<input id="exact-node-rendering" type="checkbox" checked />' in html_content
        assert 'Render every node and scale bubble area exactly' in html_content
        assert '<input id="leaf-pruning-only" type="checkbox" checked />' in html_content
        assert '<button id="sort-components-by-size" type="button" aria-pressed="true" disabled>' in html_content
        assert 'initialPosition: 0' in html_content
        assert 'reset-view' in html_content
        assert 'node-arrangement' in html_content
        assert 'Grouped subclusters' in html_content
        assert 'Radial split order' in html_content
        assert 'groupedDotLayoutForNode' in html_content
        assert 'radialDotPositions' in html_content
        assert 'orderedRadialPositions' in html_content
        assert 'positionClusterScore' in html_content
        assert 'splitRadialPositionsForSubclusters' in html_content
        assert 'assignMembersToRadialPositions' in html_content
        assert 'normalizedComponentDotLayout' in html_content