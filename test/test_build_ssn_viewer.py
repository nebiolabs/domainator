import base64
import gzip
import json
import os
import re
import tempfile

import numpy as np
import pandas as pd
import pytest

from domainator import build_ssn_viewer
from domainator.data_matrix import DataMatrix, DenseDataMatrix, MaxTree


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
        assert bundle["version"] == build_ssn_viewer.SSN_VIEWER_BUNDLE_VERSION == 4
        # app_state is written only by the HTML viewer's "Save session" button.
        assert "app_state" not in bundle
        assert bundle["graph"]["nodes"] == row_names
        assert len(bundle["graph"]["mst_edges"]) == 3
        assert len(bundle["graph"]["merge_event_series"]) == 3
        for event in bundle["graph"]["merge_event_series"]:
            counts = event["merge_size_counts"]
            assert sum(counts.values()) == event["merge_count"]
            assert max(int(size) for size in counts) == event["largest_merge"]
            assert sum(int(size) * n for size, n in counts.items()) == event["merge_impact"]
        moving_sum = bundle["graph"]["merge_moving_sum"]
        assert len(moving_sum["x"]) == len(moving_sum["y"]) > 0
        assert moving_sum["window"] > 0
        # A stop labelled T shows the `--lb T` cut, so edge_index is the last MST edge
        # scoring strictly above T (T's own tie group is excluded). The final stop is the
        # floor that keeps every MST edge, so the fully merged network is reachable --
        # without it the lowest stop still splits the weakest merge apart. It sits 1% of
        # the weight range below the weakest edge (4.0 - 0.01 * (10.0 - 4.0)).
        assert [stop["threshold_label"] for stop in bundle["graph"]["slider_stops"]] == [
            "∞", "10.00", "7.00", "4.00", "3.94"
        ]
        assert [stop["edge_index"] for stop in bundle["graph"]["slider_stops"]] == [-1, -1, 0, 1, 2]
        assert [stop["threshold_index"] for stop in bundle["graph"]["slider_stops"]] == [-1, 0, 1, 2, 3]
        assert bundle["graph"]["hierarchy"]["roots"] == [6]
        assert bundle["graph"]["hierarchy"]["leaf_order"] == [0, 1, 2, 3]
        assert bundle["graph"]["hierarchy"]["nodes"][6]["leaf_count"] == 4
        assert bundle["defaults"] == {
            "color_by": "category",
            "label_by": "category",
            "categorical_columns": [],
        }
        assert bundle["metadata"]["columns"] == [
            {"name": "category", "type": "str"},
            {"name": "count", "type": "int"},
            {"name": "score", "type": "float"},
        ]
        assert bundle["metadata"]["rows"][0] == ["alpha", 1, 1.5]


def test_build_ssn_viewer_records_categorical_columns():
    """--categorical marks numeric columns for discrete coloring in the viewer."""
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
            "--color_by", "count",
            "--categorical", "count",
        ])

        bundle = _read_bundle(bundle_file)
        assert bundle["defaults"]["categorical_columns"] == ["count"]
        # The column keeps its numeric type; only the viewer's coloring changes.
        assert {"name": "count", "type": "int"} in bundle["metadata"]["columns"]


def test_build_ssn_viewer_rejects_unknown_categorical_column():
    data = np.array([
        [0, 10],
        [10, 0],
    ], dtype=float)
    row_names = ["A", "B"]
    matrix = DenseDataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        bundle_file = os.path.join(output_dir, "test_bundle.ssnv")
        matrix.write(input_file, output_type="dense")

        with pytest.raises(ValueError, match="categorical column 'missing'"):
            build_ssn_viewer.main([
                "-i", input_file,
                "-o", bundle_file,
                "--categorical", "missing",
            ])


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

        # One row per distinct threshold now; no separate "infinite threshold" row.
        expected_counts = [
            [float(row[0]), int(row[1])]
            for row in tree.cluster_count_by_threshold
        ]

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
        # 2 capped merge events + the ∞ stop + the floor stop.
        assert len(bundle["graph"]["slider_stops"]) == 4
        assert bundle["graph"]["slider_stops"][-1]["threshold_value"] < min(
            edge[2] for edge in bundle["graph"]["mst_edges"]
        )
        assert bundle["graph"]["slider_stops"][0]["threshold_value"] is None

        # The moving sum must come from the UNFILTERED rows. The cap keeps only t=10.0 and
        # t=5.0, so a sum taken over the filtered series would span [5.0, 10.0] and never
        # see the dropped events at all.
        kept = {event["threshold_value"] for event in bundle["graph"]["merge_event_series"]}
        assert kept == {10.0, 5.0}
        moving_sum = bundle["graph"]["merge_moving_sum"]
        assert min(moving_sum["x"]) == pytest.approx(1.0)
        assert max(moving_sum["x"]) == pytest.approx(10.0)
        # The window at the bottom of the range still counts the dropped t=1.0 merge.
        assert moving_sum["y"][0] >= 1


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
            "--html", html_file,
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
        assert 'Largest single split' in html_content
        assert 'Split impact' not in html_content
        assert 'Threshold' in html_content
        assert 'componentMembers' in html_content
        assert 'mstLinksForActiveClusters' in html_content
        assert 'tidyForestLayout' in html_content
        assert 'tidyComponentLayout' in html_content
        assert 'forceDirectedForestLayout' in html_content
        assert 'gridClusterLayout' in html_content
        assert 'packedClusterLayout' in html_content
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
        assert 'metadata-rows-per-page' in html_content
        assert 'Rows per page' in html_content
        assert '<option value="all">All rows</option>' in html_content
        assert 'filteredMetadataNodeIndices' in html_content
        assert 'formatMetadataDisplayValue' in html_content
        assert 'toggleMetadataRowSelection' in html_content
        assert 'metadataRowSelectionAnchor' in html_content
        assert 'selectNodesFromMetadataRows' in html_content
        assert 'layout-algorithm' in html_content
        assert 'Layout algorithm' in html_content
        assert '<option value="tree">Tree</option>' in html_content
        assert '<option value="force">Force-directed</option>' in html_content
        assert '<option value="organic">' not in html_content
        assert '<option value="grid">Grid (no edges)</option>' in html_content
        assert '<option value="packed" selected>Packed (no edges)</option>' in html_content
        assert '<option value="treemap">Treemap (no edges)</option>' in html_content
        assert '<option value="treemap-recursive"' not in html_content
        assert 'treemapClusterLayout' in html_content
        assert 'componentSquareLayout' in html_content
        assert 'componentMemberLayout' in html_content
        assert 'gilbertCurve' in html_content
        assert 'ensureLatticeGlobal' in html_content
        assert 'latticeNodeAtWorld' in html_content
        assert 'rectIntersectsRect' in html_content
        assert '<input id="min-cluster-size" type="number" min="1" value="1" step="1" />' in html_content
        assert 'leaf-pruning-only' in html_content
        assert 'Minimum cluster size trims leaf clusters only' in html_content
        assert 'show-node-counts' in html_content
        assert 'Show node count labels' in html_content
        assert 'show-edge-scores' in html_content
        assert 'Show edge score labels' in html_content
        assert '<input id="leaf-pruning-only" type="checkbox" />' in html_content
        assert '<button id="sort-components-by-size" type="button" aria-pressed="true" disabled>' in html_content
        assert 'initialPosition: 0' in html_content
        assert 'reset-view' in html_content
        # The node-arrangement dropdown and exact-node-rendering toggle were removed:
        # dots are always grouped and every node is always rendered.
        assert 'node-arrangement' not in html_content
        assert 'exact-node-rendering' not in html_content
        assert 'groupedDotLayout' in html_content
        assert 'selectByCoord' in html_content
        assert 'collectMajorChildren' in html_content
        assert 'partitionAmongChildren' in html_content
        assert 'radialDotPositions' in html_content
        assert 'normalizedComponentDotLayout' in html_content
        assert 'reduce-elongation' in html_content
        assert 'reduceElongationEnabled' in html_content
        assert 'regionPrincipalAxis' in html_content
        assert 'selectByProjection' in html_content

        # Session save/load (bundle v4 app_state).
        assert '<button id="save-session" type="button" disabled' in html_content
        assert 'const VIEW_STATE_FIELDS = [' in html_content
        assert 'function collectSessionState()' in html_content
        assert 'function applySessionState(appState)' in html_content
        assert 'function saveSessionFile()' in html_content
        assert '<button id="save-extraction"' in html_content
        assert 'function saveExtractionFile()' in html_content
        assert 'function buildExtractionHierarchy(nodeCount, edges)' in html_content
        assert 'function extractionMergeEventRows(nodeCount, edges, metric)' in html_content
        assert 'function originalComponentByNode()' in html_content
        assert 'const SUPPORTED_BUNDLE_VERSIONS = [3, 4];' in html_content
        assert 'SUPPORTED_BUNDLE_VERSIONS.includes(bundle.version)' in html_content

        # Selection presets: ten slots, keyboard-addressable, hover preview.
        assert '<div id="preset-slots"' in html_content
        for slot in range(10):
            assert f'data-preset-slot="{slot}"' in html_content
        assert 'function storeSelectionPreset(slot)' in html_content
        assert "function recallSelectionPreset(slot, mode = 'replace')" in html_content
        assert 'function presetPreviewNodeSet()' in html_content
        assert 'function presetClickMode(event)' in html_content
        assert '/^Digit([0-9])$/' in html_content

        # Metadata editing.
        assert 'function setMetadataValue(nodeIndex, columnName, rawText)' in html_content
        assert 'function beginMetadataCellEdit(cell)' in html_content
        assert 'function addMetadataColumn(name, columnType)' in html_content
        assert 'function deleteMetadataColumn(columnName)' in html_content
        assert 'function applyBulkFill()' in html_content
        assert 'function applyPastedColumn()' in html_content
        assert '<input id="metadata-new-column-name"' in html_content
        assert '<select id="metadata-fill-target"' in html_content
        assert '<option value="all">All nodes</option>' in html_content
        assert '<div id="metadata-paste-overlay"' in html_content
        # Editing controls live behind three collapsed disclosure panels.
        for panel in ("add", "set", "rename", "delete"):
            assert f'<button id="metadata-panel-{panel}"' in html_content
            assert f'<div id="metadata-{panel}-panel" class="metadata-edit-panel" hidden>' in html_content
        assert 'function toggleMetadataEditPanel(name)' in html_content
        assert 'function renameMetadataColumn(oldName, newName)' in html_content
        assert 'function renameColumnInMenus(oldName, newName)' in html_content
        assert 'function addClusterColumn(name)' in html_content
        assert 'function clusterNumbersAtCurrentThreshold()' in html_content
        assert '<button id="metadata-add-cluster-column"' in html_content


def test_build_ssn_viewer_writes_static_html_without_input():
    with tempfile.TemporaryDirectory() as output_dir:
        html_file = os.path.join(output_dir, "viewer.html")

        build_ssn_viewer.main([
            "--html", html_file,
            "--name", "Viewer Only",
        ])

        html_content = open(html_file, "r", encoding="utf-8").read()

        assert '<title>Viewer Only</title>' in html_content
        # The tab keeps the bare name; the heading spells out the app.
        assert ('<h1 id="viewer-title">Domainator Similarity Network Viewer: Viewer Only</h1>'
                in html_content)
        assert 'const EMBEDDED_BUNDLE_BASE64 = null;' in html_content
        assert 'No bundle loaded.' in html_content


def test_build_ssn_viewer_embeds_data_in_viewer_html():
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
        html_file = os.path.join(output_dir, "viewer.html")
        matrix.write(input_file, output_type="dense")

        build_ssn_viewer.main([
            "-i", input_file,
            "--html", html_file,
            "--embed_data",
            "--name", "Embedded Viewer",
        ])

        html_content = open(html_file, "r", encoding="utf-8").read()

        assert 'const EMBEDDED_BUNDLE_BASE64 = ' in html_content
        assert 'Loading bundled data...' in html_content
        assert 'autoloadEmbeddedBundle()' in html_content

        match = re.search(r"const EMBEDDED_BUNDLE_BASE64 = \"([^\"]+)\";", html_content)
        assert match is not None

        embedded_bundle = json.loads(gzip.decompress(base64.b64decode(match.group(1))).decode("utf-8"))
        assert embedded_bundle["name"] == "Embedded Viewer"
        assert embedded_bundle["graph"]["nodes"] == row_names


def test_build_ssn_viewer_requires_output_or_embed_when_input_supplied():
    data = np.array([
        [0, 10],
        [10, 0],
    ], dtype=float)
    row_names = ["A", "B"]
    matrix = DenseDataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        matrix.write(input_file, output_type="dense")

        with pytest.raises(SystemExit, match="provide -o/--output or --embed_data"):
            build_ssn_viewer.main([
                "-i", input_file,
            ])


def test_build_ssn_viewer_rejects_output_without_input():
    with tempfile.TemporaryDirectory() as output_dir:
        bundle_file = os.path.join(output_dir, "test_bundle.ssnv")

        with pytest.raises(SystemExit, match="-o/--output requires -i/--input"):
            build_ssn_viewer.main([
                "-o", bundle_file,
            ])


def test_build_ssn_viewer_rejects_embed_without_viewer_html():
    data = np.array([
        [0, 10],
        [10, 0],
    ], dtype=float)
    row_names = ["A", "B"]
    matrix = DenseDataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        matrix.write(input_file, output_type="dense")

        with pytest.raises(SystemExit, match="--embed_data requires --html"):
            build_ssn_viewer.main([
                "-i", input_file,
                "--embed_data",
            ])


def test_build_ssn_viewer_requires_viewer_html_without_input():
    with pytest.raises(SystemExit, match="--html is required"):
        build_ssn_viewer.main([])

def test_load_bundle_accepts_every_supported_version():
    """v4 adds only the ignorable app_state section, so v3 files still load."""
    from domainator import ssn_bundle

    data = np.array([
        [0, 10, 6, 0],
        [10, 0, 7, 0],
        [6, 7, 0, 4],
        [0, 0, 4, 0],
    ], dtype=float)
    row_names = ["A", "B", "C", "D"]
    matrix = DenseDataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "matrix.hdf5")
        matrix.write(input_file, output_type="dense")
        bundle = build_ssn_viewer.build_ssn_viewer_bundle(
            DataMatrix.from_file(input_file), name="versions"
        )

        assert ssn_bundle.SUPPORTED_SSN_VIEWER_BUNDLE_VERSIONS == (3, 4)

        for version in ssn_bundle.SUPPORTED_SSN_VIEWER_BUNDLE_VERSIONS:
            path = os.path.join(output_dir, f"v{version}.ssnv")
            payload = dict(bundle, version=version)
            if version >= 4:
                # A saved session; the reader must ignore the extra section.
                payload["app_state"] = {"state_version": 1, "view": {"color_by": None}}
            build_ssn_viewer.write_ssn_viewer_bundle(path, payload)
            assert ssn_bundle.load_bundle(path)["version"] == version

        unsupported = os.path.join(output_dir, "v99.ssnv")
        build_ssn_viewer.write_ssn_viewer_bundle(unsupported, dict(bundle, version=99))
        with pytest.raises(ValueError, match="Unsupported SSN viewer bundle version 99"):
            ssn_bundle.load_bundle(unsupported)


def test_viewer_heading_names_the_network():
    """The <h1> is "<app name>: <network name>", with no name for a bare shell."""
    from domainator.ssn_viewer_html import VIEWER_APP_NAME, viewer_heading

    assert viewer_heading("GH17") == f"{VIEWER_APP_NAME}: GH17"
    assert viewer_heading("  GH17  ") == f"{VIEWER_APP_NAME}: GH17"
    # Titles that name no particular network collapse to the app name alone,
    # rather than "<app name>: Domainator SSN Viewer".
    for generic in (None, "", "   ", "Domainator SSN Viewer", VIEWER_APP_NAME):
        assert viewer_heading(generic) == VIEWER_APP_NAME


def test_lowest_slider_stop_reaches_the_fully_merged_network():
    """The lowest stop must keep every MST edge.

    Each stop excludes its own tie group (the strictly-above `--lb` convention),
    so without a floor stop below the weakest edge the slider could never show
    the true connected components -- it stopped one merge short.
    """
    from domainator import ssn_bundle

    # Two components: A-B-C joined at 10/6, and D-E joined at 8.
    data = np.array([
        [0, 10, 6, 0, 0],
        [10, 0, 7, 0, 0],
        [6, 7, 0, 0, 0],
        [0, 0, 0, 0, 8],
        [0, 0, 0, 8, 0],
    ], dtype=float)
    row_names = ["A", "B", "C", "D", "E"]
    matrix = DenseDataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "matrix.hdf5")
        matrix.write(input_file, output_type="dense")
        bundle = build_ssn_viewer.build_ssn_viewer_bundle(
            DataMatrix.from_file(input_file), name="floor"
        )

    hierarchy = bundle["graph"]["hierarchy"]
    stops = bundle["graph"]["slider_stops"]
    lowest = min(stop["threshold_value"] for stop in stops if stop["threshold_value"] is not None)
    weakest_mst_edge = min(edge[2] for edge in bundle["graph"]["mst_edges"])

    # Strictly below the weakest edge, so no merge is excluded at this cut, and
    # only 1% of the weight range below it so the slider track is not mostly empty.
    mst_weights = [edge[2] for edge in bundle["graph"]["mst_edges"]]
    span = max(mst_weights) - min(mst_weights)
    assert lowest < weakest_mst_edge
    assert lowest == pytest.approx(weakest_mst_edge - 0.01 * span)
    # At the floor stop every merge is applied, so the clusters are the components.
    assert ssn_bundle.clusters_at_threshold(hierarchy, lowest) == hierarchy["roots"]
    assert len(hierarchy["roots"]) == 2
    # Before the fix the lowest stop was the weakest edge itself, which the
    # strictly-above rule excludes; the viewer's cut there splits it back apart.
    assert stops[-1]["threshold_value"] == lowest
    assert stops[-2]["threshold_value"] == weakest_mst_edge
