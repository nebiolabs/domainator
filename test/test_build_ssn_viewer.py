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
from domainator import ssn_bundle


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
        bundle_file = os.path.join(output_dir, "test_bundle.dsnv")

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
        assert bundle["version"] == build_ssn_viewer.SSN_VIEWER_BUNDLE_VERSION == 6
        # app_state is written only by the HTML viewer's "Save session" button.
        assert "app_state" not in bundle
        assert bundle["graph"]["nodes"] == row_names
        assert len(bundle["graph"]["mst_edges"]) == 3

        # The split-event series, its moving sum and the slider's stops are derived from
        # the merge order rather than stored, so the bundle carries none of them.
        for derived_key in ("merge_event_series", "merge_event_total", "max_merge_events",
                            "merge_moving_sum", "slider_stops",
                            "cluster_count_by_threshold", "edges_by_threshold"):
            assert derived_key not in bundle["graph"]

        events = ssn_bundle.merge_event_rows(bundle)
        assert len(events) == 3
        for event in events:
            counts = event["merge_size_counts"]
            assert sum(counts.values()) == event["merge_count"]
            assert max(int(size) for size in counts) == event["largest_merge"]
            assert sum(int(size) * n for size, n in counts.items()) == event["merge_impact"]
        moving_sum = ssn_bundle.moving_sum(bundle, event_rows=events)
        assert len(moving_sum["x"]) == len(moving_sum["y"]) > 0
        assert moving_sum["window"] > 0
        # A stop labeled T shows the `--lb T` cut, so edge_index is the last MST edge
        # scoring strictly above T (T's own tie group is excluded). The final stop is the
        # floor that keeps every MST edge, so the fully merged network is reachable --
        # without it the lowest stop still splits the weakest merge apart. It sits 1% of
        # the weight range below the weakest edge (4.0 - 0.01 * (10.0 - 4.0)).
        stops = ssn_bundle.slider_stops(bundle, event_rows=events)
        assert [stop["threshold_label"] for stop in stops] == [
            "∞", "10.00", "7.00", "4.00", "3.94"
        ]
        assert [stop["edge_index"] for stop in stops] == [-1, -1, 0, 1, 2]
        # threshold_index indexed threshold tables describing the full input graph.
        # The bundle never carried those, so the stops no longer pretend to point at them.
        assert all("threshold_index" not in stop for stop in stops)
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
        bundle_file = os.path.join(output_dir, "test_bundle.dsnv")

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
        bundle_file = os.path.join(output_dir, "test_bundle.dsnv")
        matrix.write(input_file, output_type="dense")

        with pytest.raises(ValueError, match="categorical column 'missing'"):
            build_ssn_viewer.main([
                "-i", input_file,
                "-o", bundle_file,
                "--categorical", "missing",
            ])


def test_build_ssn_viewer_omits_threshold_tables_but_still_cuts_like_maxtree():
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
        bundle_file = os.path.join(output_dir, "test_bundle.dsnv")
        matrix.write(input_file, output_type="dense")

        build_ssn_viewer.main([
            "-i", input_file,
            "-o", bundle_file,
        ])

        bundle = _read_bundle(bundle_file)

        # The bundle no longer carries the threshold tables: they described the *full*
        # input graph, which a bundle never holds, and nothing in the viewer read them.
        assert "cluster_count_by_threshold" not in bundle["graph"]
        assert "edges_by_threshold" not in bundle["graph"]

        # What matters is that cutting the bundle still agrees with the matrix it came
        # from. MaxTree's table has one row per distinct MST weight, and cutting the
        # bundle's hierarchy at that weight must give the same component count.
        hierarchy = bundle["graph"]["hierarchy"]
        for threshold, expected_clusters in tree.cluster_count_by_threshold:
            active = ssn_bundle.clusters_at_threshold(hierarchy, float(threshold))
            assert len(active) == int(expected_clusters), f"clusters at lb={threshold}"

        assert bundle["graph"]["hierarchy"]["nodes"][8]["size"] == 5
        assert bundle["graph"]["hierarchy"]["nodes"][8]["leaf_count"] == 5


def test_bundle_merge_events_and_slider_stops_honour_a_caller_supplied_cap():
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
        bundle_file = os.path.join(output_dir, "test_bundle.dsnv")
        matrix.write(input_file, output_type="dense")

        build_ssn_viewer.main(["-i", input_file, "-o", bundle_file])

        bundle = _read_bundle(bundle_file)

        # The cap is applied where the series is consumed, not baked into the file, so
        # the same bundle answers for any cap.
        events = ssn_bundle.merge_event_rows(bundle)
        assert len(events) == 4

        # This network has four merge events spread over four different bands of the
        # axis, so the density back-fill restores every one the cap dropped -- the cap
        # and the back-fill are exercised against each other in test_ssn_hierarchy.py,
        # where the band layout can be controlled directly.
        selected = ssn_bundle.merge_event_series(bundle, max_merge_events=2, event_rows=events)
        assert len(selected) == 4

        stops = ssn_bundle.slider_stops(bundle, max_merge_events=2, event_rows=events)
        # 4 merge events + the ∞ stop + the floor stop.
        assert len(stops) == 6
        assert stops[-1]["threshold_value"] < min(edge[2] for edge in bundle["graph"]["mst_edges"])
        assert stops[0]["threshold_value"] is None

        # The moving sum must come from the UNFILTERED rows, whatever the cap kept.
        kept = {event["threshold_value"] for event in selected}
        assert kept == {10.0, 5.0, 4.0, 1.0}
        moving_sum = ssn_bundle.moving_sum(bundle, event_rows=events)
        assert min(moving_sum["x"]) == pytest.approx(1.0)
        assert max(moving_sum["x"]) == pytest.approx(10.0)
        # The window at the bottom of the range still counts the dropped t=1.0 merge.
        assert moving_sum["y"][0] >= 1


def test_bundle_cap_bites_on_a_network_large_enough_to_show_it():
    """The cap only shows on a network with more events than bands to spread them over.

    Three blocks of tightly connected nodes give ~60 merge events with a wide range of
    impacts, so a cap of 5 genuinely drops most of them and the back-fill adds at most
    one per 5% band.
    """
    from domainator.ssn_hierarchy import MERGE_EVENT_DENSITY_BINS

    rng = np.random.default_rng(0)
    node_count = 60
    data = np.zeros((node_count, node_count), dtype=float)
    for start, end in [(0, 22), (22, 40), (40, node_count)]:
        for i in range(start, end):
            for j in range(i + 1, end):
                data[i, j] = data[j, i] = rng.uniform(4, 12)
    row_names = [f"n{i:02d}" for i in range(node_count)]
    matrix = DenseDataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "matrix.hdf5")
        bundle_file = os.path.join(output_dir, "bundle.dsnv")
        matrix.write(input_file, output_type="dense")
        build_ssn_viewer.main(["-i", input_file, "-o", bundle_file])
        bundle = _read_bundle(bundle_file)
        graph = bundle["graph"]

        events = ssn_bundle.merge_event_rows(bundle)
        selected = ssn_bundle.merge_event_series(bundle, max_merge_events=5, event_rows=events)
        total = len(events)
        plotted = len(selected)
        assert total > 5 + MERGE_EVENT_DENSITY_BINS  # otherwise the cap proves nothing
        assert 5 <= plotted <= 5 + MERGE_EVENT_DENSITY_BINS
        assert plotted < total

        # The back-fill's purpose: the weakest plotted event reaches the bottom of the
        # axis rather than stopping wherever the strongest events happen to end.
        all_thresholds = [edge[2] for edge in graph["mst_edges"]]
        plotted_thresholds = [row["threshold_value"] for row in selected]
        axis_span = max(all_thresholds) - min(all_thresholds)
        assert (min(plotted_thresholds) - min(all_thresholds)) < 0.1 * axis_span

        # Every stop but ∞ and the floor is a plotted event.
        stops = ssn_bundle.slider_stops(bundle, max_merge_events=5, event_rows=events)
        assert len(stops) == plotted + 2

        # And the cap really is the caller's: raising it offers strictly more cut-points
        # off the very same file.
        assert len(ssn_bundle.slider_stops(bundle, max_merge_events=0, event_rows=events)) == total + 2


def test_uncapped_series_keeps_every_event():
    """max_merge_events=0 means no cap at all."""
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
        input_file = os.path.join(output_dir, "matrix.hdf5")
        bundle_file = os.path.join(output_dir, "bundle.dsnv")
        matrix.write(input_file, output_type="dense")
        build_ssn_viewer.main(["-i", input_file, "-o", bundle_file])
        bundle = _read_bundle(bundle_file)

        events = ssn_bundle.merge_event_rows(bundle)
        selected = ssn_bundle.merge_event_series(bundle, max_merge_events=0, event_rows=events)
        assert len(events) == len(selected) == 4


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
        bundle_file = os.path.join(output_dir, "test_bundle.dsnv")

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
        bundle_file = os.path.join(output_dir, "test_bundle.dsnv")
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
        assert 'View Settings' in html_content
        assert 'Node Metadata' in html_content
        assert 'Export table TSV' in html_content
        assert 'threshold-min-label' in html_content
        assert 'threshold-max-label' in html_content
        assert 'threshold-input' in html_content
        assert 'Jump to threshold' in html_content
        # Both metrics' axis titles ship, keyed by metric (see splitAxisLabels).
        from domainator.ssn_hierarchy import MERGE_IMPACT_CHOICES, merge_impact_axis_labels
        for metric in MERGE_IMPACT_CHOICES:
            assert merge_impact_axis_labels(metric)['largest'] in html_content
        # The y-axis was once titled bare "Split impact" and must not be again. The
        # product metric's hover label "Split impact within" is a different string,
        # and is meant to be there, so it is carved out rather than banned.
        assert 'Split impact' not in html_content.replace('Split impact within', '')
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

        # "Collapse long paths": the sub-option of leaf pruning that contracts a chain of
        # pass-through clusters into one dashed edge carrying the chain's weakest link.
        # It starts disabled because leaf pruning starts off.
        assert '<input id="collapse-long-paths" type="checkbox" disabled />' in html_content
        assert 'Collapse long paths' in html_content
        assert '.checkbox.sub-option {' in html_content
        assert 'function collapseLongPathsEnabled()' in html_content
        assert 'function updateCollapseLongPathsControl()' in html_content
        assert 'function collapseLongPathLinks(' in html_content
        assert (
            "domCheckboxField('view', 'collapse_long_paths', 'collapse-long-paths')"
            in html_content
        )
        assert "' path' + (collapsedPaths === 1 ? '' : 's') + ' collapsed'" in html_content
        # The dash pattern is shared by the canvas and the SVG export, so a collapsed
        # path looks the same in a figure as it did on screen.
        assert 'const LINK_DASH_ON = 6;' in html_content
        assert "clusterContext.setLineDash(link.collapsed" in html_content
        assert '''stroke-dasharray="' + LINK_DASH_ON''' in html_content
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
        assert '<canvas id="color-histogram"' in html_content
        assert '<div class="cp-range" id="color-range-slider">' in html_content
        assert '<div id="color-stop-list" class="cp-stop-list"></div>' in html_content
        assert '<button id="color-add-stop"' in html_content
        assert 'data-stop-field="hex"' in html_content
        assert 'function setGradientStopHexField(index, color)' in html_content
        assert 'function numericPaletteStops(palette, minValue, maxValue)' in html_content
        assert 'function addGradientStop()' in html_content
        assert 'function removeGradientStop(index)' in html_content
        assert 'function updateGradientSlider()' in html_content
        assert 'function setupGradientRangeSlider()' in html_content
        assert 'function columnHistogram(columnName)' in html_content
        assert 'function drawColorHistogram(columnName)' in html_content
        assert '<button id="focus-selection"' in html_content
        assert 'function selectedNodeBounds()' in html_content
        assert 'function focusSelection()' in html_content
        assert '<button id="save-extraction"' in html_content
        # The rename control lives beside the heading it edits, as a pencil glyph.
        assert '<button id="rename-network" type="button" class="title-edit"' in html_content
        assert '<div class="hero-title">' in html_content
        assert '<div id="name-overlay"' in html_content
        assert 'function openNameDialog(options)' in html_content
        assert 'function renameNetwork(name)' in html_content
        assert 'function saveExtractionFile()' in html_content
        assert 'function buildExtractionHierarchy(nodeCount, edges)' in html_content
        assert 'function mergeEventRows(nodeCount, edges, metric)' in html_content
        # The series is derived in the browser, so the shell has to carry the whole
        # pipeline and the control that sets its cap.
        assert 'function deriveMergeSeries(nodeCount, edges, metric)' in html_content
        assert 'function rebuildMergeSeries(maxMergeEvents)' in html_content
        assert '<input id="split-event-cap"' in html_content
        # The selection is taken against the chart's current window, so a zoom can
        # re-spend the cap on what is on screen.
        assert 'function selectMergeEvents(series, maxMergeEvents, window, pinnedThreshold)' in html_content
        assert 'function currentSplitWindow()' in html_content
        assert 'const DEFAULT_WINDOWED_MAX_MERGE_EVENTS = 50;' in html_content
        assert 'function originalComponentByNode()' in html_content
        assert 'const SUPPORTED_BUNDLE_VERSIONS = [3, 4, 5, 6];' in html_content

        # The split chart draws a capped selection of the split events, so it has to
        # say how much of the series that is.
        assert '<div class="note" id="split-event-count">' in html_content
        assert 'function updateSplitEventCount()' in html_content
        assert 'merge events plotted' in html_content
        assert 'const MERGE_EVENT_DENSITY_BINS = 20;' in html_content
        assert 'function mergeEventDensityBin(thresholdValue, lo, hi, densityBins)' in html_content
        assert 'function compareMergeEventRank(left, right)' in html_content

        # One geometry model behind the canvas painter, the SVG export and the
        # hit-test, so the three cannot place a mark or a tick differently.
        assert 'function splitChartLayout(viewWidth, viewHeight)' in html_content
        assert 'const SPLIT_CHART_COLORS = {' in html_content
        # Tick values on round numbers, printed to the precision their step needs:
        # a tick placed at 0.8736 and labeled "0.87" reads as misaligned beside a
        # lollipop whose own readout says 0.87.
        assert 'function splitAxisTicks(min, max, targetCount, options = {})' in html_content
        assert 'function decimalsForTickStep(step)' in html_content

        # Split-event chart export (PNG re-rasterized at the shared resolution
        # selector; SVG straight from the shared layout).
        assert '<button id="export-split-png" type="button" disabled' in html_content
        assert '<button id="export-split-svg" type="button" disabled' in html_content
        assert 'function buildSplitChartSVG()' in html_content
        # The threshold cursor is UI state, so the exports leave it out: the canvas
        # painter draws it only on screen, and the SVG builder never does.
        assert 'if (onScreen && layout.markerX !== null)' in html_content
        assert 'No threshold marker: this builder is only ever an export' in html_content
        assert 'function exportSplitChartSVG()' in html_content
        assert 'function exportSplitChartPNG()' in html_content
        assert "_split_events.svg'" in html_content
        assert "'_split_events@' + scaleFactor + 'x.png'" in html_content

        # Split-chart hover readout and click-to-jump. The geometry is recorded by
        # the draw itself, so a hit-test can never disagree with what was painted.
        assert '<div id="split-chart-tip" class="split-tip"' in html_content
        assert 'function recordSplitChartGeometry(geometry)' in html_content
        assert 'function splitChartHitAt(x, y)' in html_content
        assert 'function splitChartMovingSumAt(threshold)' in html_content
        assert 'function splitChartThresholdAt(x)' in html_content
        assert 'function splitChartTipLines(hit)' in html_content
        assert 'function splitImpactAmount(value)' in html_content
        assert 'function handleSplitChartClick(event)' in html_content
        assert 'function setupSplitChartHover()' in html_content
        assert 'Click to jump here' in html_content
        assert 'Click to jump to the nearest split' in html_content

        # Label level-of-detail is gated on each mark's on-screen size, not on a
        # single zoom threshold, and the canvas and SVG renderers share the rules.
        assert 'function itemScreenExtent(item)' in html_content
        assert 'function clusterCountLabelFits(text, item)' in html_content
        assert 'function edgeScoreLabelFits(text, linkScreenLength)' in html_content
        assert 'const MIN_LABELED_DOT_SCREEN_RADIUS' in html_content
        assert 'state.viewTransform.scale >= 0.11' not in html_content
        assert 'state.viewTransform.scale >= 0.16' not in html_content
        # A product merge impact is not a node count, in either chart's wording.
        assert 'Nodes displaced within' in html_content
        assert 'Split impact within' in html_content
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
        # Editing controls live behind collapsed disclosure panels.
        for panel in ("add", "set", "rename", "delete", "select"):
            assert f'<button id="metadata-panel-{panel}"' in html_content
            assert f'<div id="metadata-{panel}-panel" class="metadata-edit-panel" hidden>' in html_content
        assert 'function toggleMetadataEditPanel(name)' in html_content
        assert 'function renameMetadataColumn(oldName, newName)' in html_content
        assert 'function renameColumnInMenus(oldName, newName)' in html_content
        assert 'function addClusterColumn(name)' in html_content
        assert 'function clusterNumbersAtCurrentThreshold()' in html_content
        assert '<button id="metadata-add-cluster-column"' in html_content

        # Per-column charts and frequency tables: a glyph in each column
        # header opens a menu of chart kinds, and picking one opens a dialog.
        assert 'metadata-chart-button' in html_content
        assert 'data-chart-column=' in html_content
        assert 'metadata-chart-glyph' in html_content
        assert '<div id="column-chart-menu"' in html_content
        assert '<div id="column-chart-overlay"' in html_content
        assert '<div id="column-chart-preview"' in html_content
        assert '<select id="column-chart-kind"' in html_content
        assert '<button id="column-chart-copy-tsv"' in html_content
        assert '<button id="column-chart-download-tsv"' in html_content
        assert '<button id="column-chart-export-svg"' in html_content
        assert '<button id="column-chart-export-png"' in html_content
        assert 'function columnChartKindsFor(columnName)' in html_content
        assert 'function columnChartNodeIndices()' in html_content
        assert 'function columnValueDistribution(columnName, nodeIndices, options = {})' in html_content
        assert 'function columnNumericValues(columnName, nodeIndices)' in html_content
        assert 'function scopedColumnHistogram(columnName, nodeIndices)' in html_content
        assert 'function numericSummary(values)' in html_content
        assert 'function columnSummaryRows(columnName, nodeIndices)' in html_content
        assert 'function frequencyTableRows(model)' in html_content
        assert 'function buildBarChartSVG(model, meta)' in html_content
        assert 'function buildPieChartSVG(model, meta)' in html_content
        assert 'function buildHistogramSVG(histogram, meta)' in html_content
        assert 'function buildBoxPlotSVG(summary, meta)' in html_content
        assert 'function buildEcdfSVG(values, meta)' in html_content
        assert 'function buildTableSVG(rows, columns, meta)' in html_content
        assert 'function buildColumnChartArtifact(columnName, kind, nodeIndices)' in html_content
        assert 'function openColumnChartMenu(columnName, anchorButton)' in html_content
        assert 'function renderColumnChart()' in html_content
        assert 'function exportColumnChartSVG()' in html_content
        assert 'function columnChartTSV()' in html_content
        assert 'function setupColumnCharts()' in html_content
        # Alongside each value's share of the charted rows, its share of that
        # value's own network-wide population -- the question a selection raises.
        assert 'function columnGlobalCounts(columnName)' in html_content
        assert ('<th class="cc-num">All nodes</th><th class="cc-num">% of all</th>'
                in html_content)
        assert "{key: 'globalPercent', label: '% of all', align: 'right'}" in html_content
        assert ("['value', 'count', 'percent', 'count_all_nodes', 'percent_of_all', 'color']"
                in html_content)

        # Select by value: a metadata query that edits the node selection.
        assert '<select id="metadata-select-column"' in html_content
        assert '<select id="metadata-select-op"' in html_content
        assert '<input id="metadata-select-value"' in html_content
        assert '<input id="metadata-select-value2"' in html_content
        for action in ("add", "remove", "subset"):
            assert f'<button id="metadata-select-{action}"' in html_content
        assert 'function metadataSelectMatcher(opId, firstText, secondText)' in html_content
        assert 'function metadataSelectMatches(nodeIndices, field, matcher)' in html_content
        assert 'function applyMetadataSelectByValue(mode)' in html_content
        assert 'const METADATA_SELECT_OPS' in html_content

        # "Jump to threshold" steps between split-plot stops, so the field is a
        # plain text input (a number input's spinner would step by a constant).
        assert '<button id="threshold-step-down"' in html_content
        assert '<button id="threshold-step-up"' in html_content
        assert '<input id="threshold-input" type="text"' in html_content
        assert 'function stepThreshold(delta)' in html_content
        assert 'function updateThresholdStepButtons()' in html_content
        # An arrow press has to land on a threshold the view has not just been at.
        # Stops crowd onto shared slider positions, so the stop that was chosen is
        # remembered rather than re-derived from the position it sits at.
        assert 'function nextDistinctStopIndex(stops, index, delta)' in html_content
        # And a drag has to be able to land on every stop, which the plain
        # value-linear scale could not do where stops crowd: positions reserve one
        # slot per stop and spend the rest through a warp of the threshold axis,
        # weighted toward whatever stretch the split chart is zoomed into.
        assert 'function sliderValueWarp(lowValue, highValue)' in html_content
        assert 'function positionSliderStops(finiteStops)' in html_content
        assert 'function repositionSliderStops()' in html_content
        assert 'const SLIDER_MAX_FINITE_POSITION = 920;' in html_content
        assert 'const SLIDER_INFINITY_POSITION = 1000;' in html_content
        assert 'if (state.selectedStop && state.selectedStop.sliderPosition === sliderPosition)' in html_content

        # The split chart zooms and pans over its threshold axis.
        assert '<button id="split-chart-reset-zoom"' in html_content
        assert '<span id="split-chart-zoom-hint"' in html_content
        assert 'function splitChartVisibleWindow(dataMin, dataMax)' in html_content
        assert 'function setSplitChartWindow(min, max, dataMin, dataMax)' in html_content
        assert 'function zoomSplitChartAt(canvasX, factor)' in html_content
        assert 'function panSplitChartByPixels(dx)' in html_content
        assert 'function handleSplitChartWheel(event)' in html_content
        assert 'function resetSplitChartZoom()' in html_content
        assert 'function updateSplitChartZoomControls()' in html_content
        assert 'function splitChartWindowLabel(value, span)' in html_content
        assert 'function handleSplitChartDoubleClick()' in html_content
        # The page must not scroll out from under the zoom gesture.
        assert "splitCanvas.addEventListener('wheel', handleSplitChartWheel, {passive: false})" in html_content
        # A zoomed window means marks outside it, so both painters clip the plot box.
        assert '<clipPath id="split-plot-clip">' in html_content
        assert 'clip-path="url(#split-plot-clip)"' in html_content
        # Saved with the session, like the canvas's own pan and zoom.
        assert "key: 'split_chart_zoom'," in html_content

        # Categorical columns default to get_palette's own colors, so the
        # palette menu carries only real named palettes.
        assert "const DEFAULT_CATEGORICAL_PALETTE = 'domainator'" in html_content
        assert 'function defaultCategoricalPalette(columnName)' in html_content
        assert '__default__' not in html_content
        assert 'hashed hues' not in html_content

        # The page body is one big f-string, so a JS/CSS brace that was not
        # doubled -- or a doubled brace that leaked out of one of the plain
        # (non-f) string JS modules -- shows up as a literal '{{' here.
        assert '{{' not in html_content


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
        assert ('>Domainator Similarity Network Viewer: Viewer Only</h1>'
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
        bundle_file = os.path.join(output_dir, "test_bundle.dsnv")

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
    """Older bundles still load: everything this build reads has been there since v3.

    v4 added the ignorable app_state section and v5 two event-count keys, both
    additive. v6 *removed* the derived threshold data, so a v3/v4/v5 file carries keys
    a v6 reader simply ignores -- it derives those values from the merge order instead.
    """
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

        assert ssn_bundle.SUPPORTED_SSN_VIEWER_BUNDLE_VERSIONS == (3, 4, 5, 6)

        stale_derived_keys = {
            "cluster_count_by_threshold": [[10.0, 1.0]],
            "edges_by_threshold": [[0.0, 10.0]],
            "merge_event_series": [],
            "merge_event_total": 0,
            "max_merge_events": 500,
            "merge_moving_sum": {"window": 0.0, "x": [], "y": []},
            "slider_stops": [],
        }

        for version in ssn_bundle.SUPPORTED_SSN_VIEWER_BUNDLE_VERSIONS:
            path = os.path.join(output_dir, f"v{version}.dsnv")
            payload = dict(bundle, version=version)
            if version < 6:
                # A pre-v6 writer stored the derived data in the file. Loading must not
                # depend on it -- and, crucially, must not *use* it either: the stored
                # series here is empty and the real one has three events.
                payload["graph"] = dict(bundle["graph"], **stale_derived_keys)
            if version >= 4:
                # A saved session; the reader must ignore the extra section.
                payload["app_state"] = {"state_version": 1, "view": {"color_by": None}}
            build_ssn_viewer.write_ssn_viewer_bundle(path, payload)
            loaded = ssn_bundle.load_bundle(path)
            assert loaded["version"] == version
            # Derived from the merge order, so an old file's stale copy is ignored.
            assert len(ssn_bundle.merge_event_rows(loaded)) == 3
            assert len(ssn_bundle.slider_stops(loaded)) == 5

        unsupported = os.path.join(output_dir, "v99.dsnv")
        build_ssn_viewer.write_ssn_viewer_bundle(unsupported, dict(bundle, version=99))
        with pytest.raises(ValueError, match="Unsupported Domainator Similarity Network Viewer bundle version 99"):
            ssn_bundle.load_bundle(unsupported)


def test_viewer_heading_names_the_network():
    """The <h1> is "<app name>: <network name>", with no name for a bare shell."""
    from domainator.ssn_viewer_html import VIEWER_APP_NAME, viewer_heading

    assert viewer_heading("GH17") == f"{VIEWER_APP_NAME}: GH17"
    assert viewer_heading("  GH17  ") == f"{VIEWER_APP_NAME}: GH17"
    # Titles that name no particular network collapse to the app name alone,
    # rather than "<app name>: <app name>", including the pre-rename app name.
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
    stops = ssn_bundle.slider_stops(bundle)
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
