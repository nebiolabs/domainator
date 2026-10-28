import warnings
warnings.filterwarnings("ignore", module='numpy')
import io
from domainator import matrix_report
from domainator.data_matrix import DataMatrix
from domainator.ssn_edges import mst_knn_edge_counts_by_threshold
from domainator.matrix_report import MaxTree
import tempfile
import pytest
import numpy as np
import scipy.sparse
import os
from helpers import report_data


def _page_data(html_path):
    """The data the --html_out page was written with (see web/reports/matrix_report.js)."""
    envelope = report_data(open(html_path).read())
    assert envelope["report"] == "matrix_report"
    assert envelope["encoding"] == "gzip+base64"
    return envelope["data"]


@pytest.mark.parametrize("input_file",
[
    "scorefull.tsv",
    "scorefull.dense.hdf5"
])
def test_matrix_report_1(shared_datadir, input_file):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out_html = output_dir + "/matrix_report_test.html"
        out_txt = output_dir + "/matrix_report_test.txt"
        matrix_report.main(["-i", str(shared_datadir / input_file), "-o", out_txt, "--html_out", out_html])
        f_txt = open(out_txt).read()
        assert "Matrix Report" in f_txt
        assert "Min" in f_txt
        assert "Max" in f_txt
        # 20 nodes in the test data
        assert "Nodes: 20" in f_txt
        assert 'Projected MST + kNN edge counts' not in f_txt
        assert "<title>Matrix Report</title>" in open(out_html).read()
        data = _page_data(out_html)
        assert data["summary"]["nodes"] == 20
        assert data["knn_projection"] is None


def test_matrix_report_rejects_invalid_merge_impact_metric():
    matrix = DataMatrix(np.array([[0.0, 1.0], [1.0, 0.0]]), ['A', 'B'], ['A', 'B'])

    with pytest.raises(ValueError, match="merge_impact_metric must be one of"):
        matrix_report.matrix_report(matrix, text_out_handle=io.StringIO(), html_out_handle=None, merge_impact_metric="bad")



class TestInteractiveHTML:
    """Test suite for interactive HTML generation"""
    
    def test_export_for_interactive_viz(self):
        """Test that export_for_interactive_viz returns proper structure"""
        data = np.array([
            [0, 10, 5],
            [10, 0, 8],
            [5, 8, 0]
        ])
        row_names = ['A', 'B', 'C']
        matrix = DataMatrix(data, row_names, row_names)
        
        tree = MaxTree(matrix)
        viz_data = tree.export_for_interactive_viz()
        
        # Check structure
        assert 'n_nodes' in viz_data
        assert 'mst_edges' in viz_data
        assert viz_data['n_nodes'] == 3
        assert len(viz_data['mst_edges']) == 2
        
        # Check edge format
        for edge in viz_data['mst_edges']:
            assert len(edge) == 3
            assert isinstance(edge[0], int)
            assert isinstance(edge[1], int)
            assert isinstance(edge[2], float)
    
  
    def test_interactive_html_with_main(self):
        """Test generating interactive HTML through main function"""
        # Create a temporary symmetric matrix for testing
        data = np.array([
            [0, 10, 5, 3],
            [10, 0, 4, 6],
            [5, 4, 0, 9],
            [3, 6, 9, 0]
        ])
        row_names = ['A', 'B', 'C', 'D']
        matrix = DataMatrix(data, row_names, row_names)
        
        with tempfile.TemporaryDirectory() as output_dir:
            input_file = os.path.join(output_dir, "test_matrix.hdf5")
            output_file = os.path.join(output_dir, "interactive_test.html")
            
            # Save the matrix
            matrix.write(input_file, output_type="hdf")
            
            matrix_report.main([
                "-i", input_file,
                "--html_out", output_file,
                "--knn_projection"
            ])
            
            assert os.path.exists(output_file)
            
            with open(output_file, 'r') as f:
                content = f.read()
            
            assert 'Matrix Report' in content
            # Everything the page redraws the network from: the MST, the threshold
            # tables, and the MST + kNN edge counts for each threshold and k.
            data = _page_data(output_file)
            assert data["summary"]["nodes"] == 4
            assert len(data["mst_edges"]) == 3
            table = data["threshold_table"]
            assert len(table["threshold"]) == len(table["edges"]) == len(table["clusters"])
            knn = data["knn_projection"]
            assert (knn["min_k"], knn["max_k"], knn["default_k"]) == (1, 3, 3)
            assert len(knn["counts"]) == len(table["threshold"])
            assert all(len(row) == knn["max_k"] for row in knn["counts"])


def test_mst_knn_edge_counts_by_threshold_monotonic():
    data = np.array([
        [0, 10, 8, 7],
        [10, 0, 6, 5],
        [8, 6, 0, 9],
        [7, 5, 9, 0],
    ])
    row_names = ['A', 'B', 'C', 'D']
    matrix = DataMatrix(data, row_names, row_names)
    tree = MaxTree(matrix)

    counts = mst_knn_edge_counts_by_threshold(matrix, tree, 3)

    # One row per distinct threshold, each describing the `--lb <threshold>` cut
    assert counts.shape == (len(tree.thresholds), 3)
    assert np.all(counts[:, 0] >= tree.mst_edges_above_threshold)
    assert np.all(counts[:, 1] >= counts[:, 0])
    assert np.all(counts[:, 2] >= counts[:, 1])


def test_mst_knn_edge_counts_use_the_first_max_k_of_longer_rankings():
    from domainator.ssn_edges import build_symmetric_neighbor_rankings

    rng = np.random.default_rng(3)
    upper = np.triu(rng.integers(1, 9, size=(12, 12)).astype(float), k=1)
    data = upper + upper.T
    names = [str(i) for i in range(12)]
    matrix = DataMatrix(data, names, names)
    tree = MaxTree(matrix)
    expected = mst_knn_edge_counts_by_threshold(matrix, tree, 3)
    full_rankings = build_symmetric_neighbor_rankings(matrix)  # 11 neighbors per row
    assert np.array_equal(mst_knn_edge_counts_by_threshold(None, tree, 3, neighbor_rankings=full_rankings), expected)


def test_mst_knn_edge_counts_by_threshold_sparse_matches_dense_bruteforce():
    data = np.array([
        [0, 10, 8, 0, 4],
        [10, 0, 6, 6, 0],
        [8, 6, 0, 9, 9],
        [0, 6, 9, 0, 5],
        [4, 0, 9, 5, 0],
    ], dtype=float)
    labels = ['A', 'B', 'C', 'D', 'E']
    dense_matrix = DataMatrix(data, labels, labels)
    sparse_matrix = DataMatrix(scipy.sparse.csr_array(data), labels, labels)
    dense_tree = MaxTree(dense_matrix)
    sparse_tree = MaxTree(sparse_matrix)

    def brute_force_counts(array, tree, max_k):
        counts = np.zeros((len(tree.thresholds), max_k), dtype=int)

        # `--lb threshold` semantics throughout: only scores strictly above the
        # threshold survive, for both the MST prefix and the kNN candidates.
        for threshold_idx, threshold in enumerate(tree.thresholds):
            mst_prefix_edges = {
                (source_idx, target_idx) if source_idx < target_idx else (target_idx, source_idx)
                for source_idx, target_idx, weight in tree.mst_edges if weight > threshold
            }

            edge_min_rank = {}
            for row_idx in range(array.shape[0]):
                row_scores = np.maximum(array[row_idx, :], array[:, row_idx])
                candidates = [
                    (target_idx, row_scores[target_idx])
                    for target_idx in range(array.shape[0])
                    if target_idx != row_idx and row_scores[target_idx] > threshold
                ]
                candidates.sort(key=lambda item: (-item[1], item[0]))

                for rank, (target_idx, _) in enumerate(candidates[:max_k], start=1):
                    edge = (row_idx, target_idx) if row_idx < target_idx else (target_idx, row_idx)
                    previous_rank = edge_min_rank.get(edge)
                    if previous_rank is None or rank < previous_rank:
                        edge_min_rank[edge] = rank

            non_mst_rank_counts = np.zeros(max_k + 1, dtype=int)
            for edge, rank in edge_min_rank.items():
                if edge not in mst_prefix_edges:
                    non_mst_rank_counts[rank] += 1

            counts[threshold_idx, :] = len(mst_prefix_edges) + np.cumsum(non_mst_rank_counts)[1:]

        return counts

    expected = brute_force_counts(data, dense_tree, 4)
    dense_counts = mst_knn_edge_counts_by_threshold(dense_matrix, dense_tree, 4)
    sparse_counts = mst_knn_edge_counts_by_threshold(sparse_matrix, sparse_tree, 4)

    assert np.array_equal(dense_counts, expected)
    assert np.array_equal(sparse_counts, expected)


def test_matrix_report_text_includes_mst_knn_projection(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out_html = output_dir + "/matrix_report_test.html"
        out_txt = output_dir + "/matrix_report_test.txt"
        matrix_report.main(["-i", str(shared_datadir / "scorefull.tsv"), "-o", out_txt, "--html_out", out_html, "--knn_projection"])

        text_content = open(out_txt).read()

        assert 'Projected MST + kNN edge counts' in text_content
        assert 'MST + kNN edges' in text_content
        knn = _page_data(out_html)["knn_projection"]
        assert (knn["min_k"], knn["max_k"], knn["default_k"]) == (1, 19, 5)


def test_matrix_report_includes_merge_event_outputs():
    data = np.array([
        [0, 10, 6, 0],
        [10, 0, 7, 0],
        [6, 7, 0, 4],
        [0, 0, 4, 0],
    ], dtype=float)
    row_names = ['A', 'B', 'C', 'D']
    matrix = DataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        out_html = os.path.join(output_dir, "matrix_report_test.html")
        out_txt = os.path.join(output_dir, "matrix_report_test.txt")

        matrix.write(input_file, output_type="hdf")
        matrix_report.main([
            "-i", input_file,
            "-o", out_txt,
            "--html_out", out_html,
        ])

        text_content = open(out_txt).read()

        assert 'Strongest MST split events' in text_content
        assert '      From        To    Impact' in text_content
        assert 'Top-2 ratio' not in text_content
        assert 'Largest cluster changes (MST-derived)' not in text_content
        data = _page_data(out_html)
        events = data['summary']['split_events']
        assert len(events) == 3
        for event in events:
            counts = event['merge_size_counts']
            assert sum(counts.values()) == event['merge_count']
            assert max(int(size) for size in counts) == event['largest_merge']
            assert sum(int(size) * n for size, n in counts.items()) == event['merge_impact']
        moving_sum = data['merge_moving_sum']
        assert len(moving_sum['x']) == len(moving_sum['y']) > 0
        assert moving_sum['window'] > 0
        assert data['slider_stops'][0]['edge_index'] == -1
        # A stop labeled T shows the `--lb T` cut, so it excludes T's own tie group:
        # edge_index is the last MST edge scoring strictly above T. The list ends with
        # the `--lb 0` floor cut, which keeps every MST edge -- without it the slider
        # stopped one merge short of the fully merged network.
        assert [stop['threshold_label'] for stop in data['slider_stops']] == [
            '∞', '10.00', '7.00', '4.00', '3.94'
        ]
        assert [stop['edge_index'] for stop in data['slider_stops'][1:]] == [-1, 0, 1, 2]
        # threshold_index points at the matching row of threshold_table and the kNN counts
        assert [stop['threshold_index'] for stop in data['slider_stops']] == [-1, 0, 1, 2, 3]
        assert data['threshold_table'] == {
            'threshold': [10.0, 7.0, 4.0, 0.0], 'edges': [0, 1, 3, 4], 'clusters': [4, 3, 2, 1],
        }
        assert data['mst_edges'] == [[1, 0, 10.0], [2, 1, 7.0], [3, 2, 4.0]]


def test_matrix_report_supports_product_merge_impact_metric():
    data = np.array([
        [0, 10, 6, 0],
        [10, 0, 7, 0],
        [6, 7, 0, 4],
        [0, 0, 4, 0],
    ], dtype=float)
    row_names = ['A', 'B', 'C', 'D']
    matrix = DataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        out_html = os.path.join(output_dir, "matrix_report_test.html")
        out_txt = os.path.join(output_dir, "matrix_report_test.txt")

        matrix.write(input_file, output_type="hdf")
        matrix_report.main([
            "-i", input_file,
            "-o", out_txt,
            "--html_out", out_html,
            "--merge_impact_metric", "product",
        ])

        text_content = open(out_txt).read()

        assert 'Strongest MST split events (metric=product)' in text_content
        # The page labels the split chart by this metric (see test_report_html.py).
        data = _page_data(out_html)
        assert data['summary']['merge_impact_metric'] == 'product'
        # A product impact is the product of the two sizes, so B-C (2 x 1) merges at 7.
        assert {event['threshold_to']: event['merge_impact'] for event in data['summary']['split_events']}['7.00'] == 2


def test_matrix_report_default_excludes_mst_knn(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out_html = output_dir + "/matrix_report_test.html"
        out_txt = output_dir + "/matrix_report_test.txt"
        matrix_report.main(["-i", str(shared_datadir / "scorefull.tsv"), "-o", out_txt, "--html_out", out_html])

        text_content = open(out_txt).read()

        assert 'Projected MST + kNN edge counts' not in text_content
        assert 'MST + kNN edges' not in text_content
        assert _page_data(out_html)['knn_projection'] is None


def test_matrix_report_max_merge_events_filters_html_payload():
    data = np.array([
        [0, 10, 0, 0, 0, 0],
        [10, 0, 7, 0, 0, 0],
        [0, 7, 0, 9, 0, 0],
        [0, 0, 9, 0, 6, 0],
        [0, 0, 0, 6, 0, 8],
        [0, 0, 0, 0, 8, 0],
    ], dtype=float)
    row_names = ['A', 'B', 'C', 'D', 'E', 'F']
    matrix = DataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        out_html = os.path.join(output_dir, "matrix_report_test.html")

        matrix.write(input_file, output_type="hdf")
        matrix_report.main([
            "-i", input_file,
            "--html_out", out_html,
            "--max_merge_events", "2",
        ])

        summary = _page_data(out_html)['summary']
        stops = _page_data(out_html)['slider_stops']

        # The data says what the series was selected from, so the chart can caption
        # itself rather than leaving a filtered stretch of axis looking like a quiet one.
        assert summary['split_event_total'] == 5
        assert summary['max_merge_events'] == 2
        # This network's five events sit in five different bands of the axis, so the
        # density back-fill restores every one the cap dropped. The cap and the
        # back-fill are pitted against each other in test_ssn_hierarchy.py, where the
        # band layout can be set directly, and below on a network big enough to crowd
        # several events into one band.
        assert [row['edge_index'] for row in summary['split_events']] == [-1, 0, 1, 2, 3]
        # --max_merge_events caps the plotted events, but the floor stop is a cut
        # rather than an event, so it survives and the fully merged view stays reachable.
        assert [stop['edge_index'] for stop in stops] == [-1, -1, 0, 1, 2, 3, 4]
        assert [stop['threshold_index'] for stop in stops] == [-1, 0, 1, 2, 3, 4, 5]
        # The floor sits just below the weakest edge, so it does not strand the other
        # stops at one end of the track.
        assert [stop['threshold_label'] for stop in stops] == [
            '∞', '10.00', '9.00', '8.00', '7.00', '6.00', '5.96'
        ]




def test_matrix_report_caps_merge_events_on_a_network_large_enough_to_bite():
    """A cap only shows when a band has more than one event competing for it."""
    from domainator.ssn_hierarchy import MERGE_EVENT_DENSITY_BINS

    rng = np.random.default_rng(0)
    node_count = 60
    data = np.zeros((node_count, node_count), dtype=float)
    for start, end in [(0, 22), (22, 40), (40, node_count)]:
        for i in range(start, end):
            for j in range(i + 1, end):
                data[i, j] = data[j, i] = rng.uniform(4, 12)
    row_names = [f"n{i:02d}" for i in range(node_count)]
    matrix = DataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        out_html = os.path.join(output_dir, "matrix_report_test.html")
        matrix.write(input_file, output_type="hdf")
        matrix_report.main([
            "-i", input_file, "--html_out", out_html, "--max_merge_events", "5",
        ])
        summary = _page_data(out_html)['summary']

        total = summary['split_event_total']
        plotted = len(summary['split_events'])
        assert summary['max_merge_events'] == 5
        assert total > 5 + MERGE_EVENT_DENSITY_BINS  # otherwise the cap proves nothing
        assert 5 <= plotted <= 5 + MERGE_EVENT_DENSITY_BINS
        assert plotted < total
        # The page captions the chart with those numbers (test_report_browser.py).


def test_matrix_report_max_merge_events_zero_includes_all_html_events():
    data = np.array([
        [0, 10, 0, 0, 0, 0],
        [10, 0, 7, 0, 0, 0],
        [0, 7, 0, 9, 0, 0],
        [0, 0, 9, 0, 6, 0],
        [0, 0, 0, 6, 0, 8],
        [0, 0, 0, 0, 8, 0],
    ], dtype=float)
    row_names = ['A', 'B', 'C', 'D', 'E', 'F']
    matrix = DataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        out_html = os.path.join(output_dir, "matrix_report_test.html")

        matrix.write(input_file, output_type="hdf")
        matrix_report.main([
            "-i", input_file,
            "--html_out", out_html,
            "--max_merge_events", "0",
        ])

        data = _page_data(out_html)

        assert len(data['summary']['split_events']) == 5
        assert [stop['edge_index'] for stop in data['slider_stops']] == [-1, -1, 0, 1, 2, 3, 4]
        assert [stop['threshold_index'] for stop in data['slider_stops']] == [-1, 0, 1, 2, 3, 4, 5]


def test_matrix_report_text_merge_events_limited_to_top_25():
    n_nodes = 30
    data = np.zeros((n_nodes, n_nodes), dtype=float)
    for idx in range(n_nodes - 1):
        weight = float(n_nodes - idx)
        data[idx, idx + 1] = weight
        data[idx + 1, idx] = weight

    row_names = [f'N{idx}' for idx in range(n_nodes)]
    matrix = DataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        out_txt = os.path.join(output_dir, "matrix_report_test.txt")

        matrix.write(input_file, output_type="hdf")
        matrix_report.main([
            "-i", input_file,
            "-o", out_txt,
        ])

        text_lines = open(out_txt).read().splitlines()

        header_index = text_lines.index('Strongest MST split events (metric=min_child):')
        table_lines = []
        for line in text_lines[header_index + 2:]:
            if not line.strip():
                break
            table_lines.append(line)

        assert len(table_lines) == 25


def test_matrix_report_text_merge_events_are_sorted_by_from_threshold():
    data = np.array([
        [0, 10, 0, 0, 0, 0],
        [10, 0, 7, 0, 0, 0],
        [0, 7, 0, 9, 0, 0],
        [0, 0, 9, 0, 6, 0],
        [0, 0, 0, 6, 0, 8],
        [0, 0, 0, 0, 8, 0],
    ], dtype=float)
    row_names = ['A', 'B', 'C', 'D', 'E', 'F']
    matrix = DataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        out_txt = os.path.join(output_dir, "matrix_report_test.txt")

        matrix.write(input_file, output_type="hdf")
        matrix_report.main([
            "-i", input_file,
            "-o", out_txt,
        ])

        text_lines = open(out_txt).read().splitlines()

        header_index = text_lines.index('Strongest MST split events (metric=min_child):')
        from_thresholds = []
        for line in text_lines[header_index + 2:]:
            if not line.strip():
                break
            from_thresholds.append(line.split()[0])

        assert from_thresholds == ['∞', '10.00', '9.00', '8.00', '7.00']

def test_matrix_report_profile_stages_emits_timings(capsys):
    data = np.array([
        [0, 10, 6, 0],
        [10, 0, 7, 0],
        [6, 7, 0, 4],
        [0, 0, 4, 0],
    ], dtype=float)
    row_names = ['A', 'B', 'C', 'D']
    matrix = DataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        out_html = os.path.join(output_dir, "matrix_report_test.html")

        matrix.write(input_file, output_type="hdf")
        matrix_report.main([
            "-i", input_file,
            "--html_out", out_html,
            "--knn_projection",
            "--profile_stages",
        ])

        stderr_output = capsys.readouterr().err

        assert 'matrix_report stage timings:' in stderr_output
        assert 'load_matrix:' in stderr_output
        assert 'build_edge_table:' in stderr_output
        assert 'build_max_tree:' in stderr_output
        assert 'build_mst_knn_neighbor_rankings:' in stderr_output
        assert 'kept_directed_edges=' in stderr_output
        assert 'estimated_counts_bytes=' in stderr_output
        assert 'build_mst_knn_counts:' in stderr_output
        assert 'output_shape=' in stderr_output
        assert 'output_bytes=' in stderr_output
        assert 'render_outputs:' in stderr_output
        assert 'matrix_report_total:' in stderr_output


def test_matrix_report_progress_emits_live_updates(capsys):
    data = np.array([
        [0, 10, 6, 0],
        [10, 0, 7, 0],
        [6, 7, 0, 4],
        [0, 0, 4, 0],
    ], dtype=float)
    row_names = ['A', 'B', 'C', 'D']
    matrix = DataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        out_html = os.path.join(output_dir, "matrix_report_test.html")

        matrix.write(input_file, output_type="hdf")
        matrix_report.main([
            "-i", input_file,
            "--html_out", out_html,
            "--progress",
        ])

        stderr_output = capsys.readouterr().err

        assert 'loading matrix from' in stderr_output
        assert 'building edge table' in stderr_output
        assert 'building maximum spanning tree' in stderr_output
        assert 'rendering outputs' in stderr_output


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

def _tie_group_matrix():
    """Symmetric matrix whose MST contains a deliberate tie group (three edges at 5)."""
    data = np.array([
        [0.0, 9.0, 5.0, 5.0, 0.0, 0.0],
        [9.0, 0.0, 5.0, 0.0, 3.0, 0.0],
        [5.0, 5.0, 0.0, 5.0, 0.0, 2.0],
        [5.0, 0.0, 5.0, 0.0, 7.0, 0.0],
        [0.0, 3.0, 0.0, 7.0, 0.0, 2.0],
        [0.0, 0.0, 2.0, 0.0, 2.0, 0.0],
    ])
    labels = ['a', 'b', 'c', 'd', 'e', 'f']
    return data, labels


def test_threshold_tables_match_build_ssn_lb():
    """Every threshold row must equal what `build_ssn --lb <threshold>` actually emits.

    This is the invariant that keeps matrix_report's numbers honest: report rows are keyed
    by distinct threshold and count only scores strictly above it, exactly as --lb does.
    """
    from domainator.build_ssn import iter_default_ssn_edges, cluster_labels_from_graph

    data, labels = _tie_group_matrix()
    for matrix in (DataMatrix(data, labels, labels, data_type="score"),
                   DataMatrix(scipy.sparse.csr_array(data), labels, labels, data_type="score")):
        tree = MaxTree(matrix)

        # a tie group is present: two MST edges share a weight, so per-MST-edge rows
        # would have disagreed with each other at that threshold
        mst_weights = [weight for _, _, weight in tree.mst_edges]
        assert len(set(mst_weights)) < len(mst_weights)
        assert len(set(tree.thresholds.tolist())) == len(tree.thresholds)

        for row_idx, threshold in enumerate(tree.thresholds):
            expected_edges = len(list(iter_default_ssn_edges(matrix, threshold)))
            expected_clusters = len(set(cluster_labels_from_graph(matrix, threshold)))

            assert tree.edges_by_threshold[row_idx, 0] == expected_edges, f"edges at lb={threshold}"
            assert tree.edges_by_threshold[row_idx, 1] == threshold
            assert tree.cluster_count_by_threshold[row_idx, 0] == threshold
            assert tree.cluster_count_by_threshold[row_idx, 1] == expected_clusters, f"clusters at lb={threshold}"
            assert tree.cluster_count_by_edge_count[row_idx].tolist() == [expected_edges, expected_clusters]

        # the `--lb 0` row closes the table out with the complete graph
        assert tree.thresholds[-1] == 0
        assert tree.edges_by_threshold[-1, 0] == np.count_nonzero(np.tril(data, k=-1))
        assert tree.cluster_count_by_threshold[-1, 1] == tree.n_nodes - len(tree.mst_edges)


def test_mst_knn_counts_match_build_ssn_lb_mst_knn():
    """Projected counts must equal what `build_ssn --lb T --mst --knn k` emits."""
    from domainator.ssn_edges import mst_knn_edge_index_dict

    data, labels = _tie_group_matrix()
    matrix = DataMatrix(data, labels, labels, data_type="score")
    tree = MaxTree(matrix)
    max_k = 3

    counts = mst_knn_edge_counts_by_threshold(matrix, tree, max_k)
    assert counts.shape == (len(tree.thresholds), max_k)

    for row_idx, threshold in enumerate(tree.thresholds):
        for k in range(1, max_k + 1):
            expected = len(mst_knn_edge_index_dict(matrix, k, lower_bound=threshold))
            assert counts[row_idx, k - 1] == expected, f"lb={threshold} k={k}"


def test_slider_reaches_the_fully_merged_network():
    """The lowest cut must keep every MST edge.

    Regression: each stop excludes its own tie group under the strictly-above
    `--lb` convention, so the lowest event stop still split the weakest merge
    apart. The slider stopped one merge short and the true connected components
    were unreachable. A floor stop below the weakest edge closes the gap.
    """
    from domainator.ssn_bundle import clusters_at_threshold
    from domainator.ssn_hierarchy import build_mst_component_hierarchy
    from domainator.data_matrix import DataMatrix
    from domainator.ssn_edges import MaxTree

    # Two components: A-B-C and D-E.
    data = np.array([
        [0, 10, 6, 0, 0],
        [10, 0, 7, 0, 0],
        [6, 7, 0, 0, 0],
        [0, 0, 0, 0, 8],
        [0, 0, 0, 8, 0],
    ], dtype=float)
    row_names = ['A', 'B', 'C', 'D', 'E']

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        out_html = os.path.join(output_dir, "matrix_report_test.html")
        DataMatrix(data, row_names, row_names).write(input_file, output_type="hdf")
        matrix_report.main(["-i", input_file, "--html_out", out_html])
        stops = _page_data(out_html)['slider_stops']
        tree = MaxTree(DataMatrix.from_file(input_file))

    hierarchy = build_mst_component_hierarchy(tree)
    finite = [stop['threshold_value'] for stop in stops if stop['threshold_value'] is not None]
    weakest_mst_edge = min(float(edge[2]) for edge in tree.mst_edges)

    # The floor stop sits strictly below the weakest edge, so no merge is cut there,
    # and only 1% of the weight range below it so the track is not mostly empty.
    mst_weights = [float(edge[2]) for edge in tree.mst_edges]
    span = max(mst_weights) - min(mst_weights)
    assert stops[-1]['threshold_value'] == min(finite)
    assert min(finite) < weakest_mst_edge
    assert min(finite) == pytest.approx(weakest_mst_edge - 0.01 * span)
    assert clusters_at_threshold(hierarchy, min(finite)) == hierarchy['roots']
    assert len(hierarchy['roots']) == 2

    # The stop above it is the weakest edge itself, which that cut splits back apart.
    assert sorted(finite)[1] == weakest_mst_edge
    assert len(clusters_at_threshold(hierarchy, weakest_mst_edge)) > len(hierarchy['roots'])


def test_slider_stops_match_the_ssn_viewer_bundle():
    """matrix_report and the viewer bundle share one stop builder, so the two tools
    must offer the same cuts for the same matrix.

    The report embeds its stops; a bundle derives them on the way out (ssn_bundle), so
    this compares what each one actually hands its slider.
    """
    from domainator import build_ssn, ssn_bundle
    from domainator.ssn_bundle import load_bundle

    data = np.array([
        [0, 10, 6, 0, 0],
        [10, 0, 7, 0, 0],
        [6, 7, 0, 4, 0],
        [0, 0, 4, 0, 8],
        [0, 0, 0, 8, 0],
    ], dtype=float)
    row_names = ['A', 'B', 'C', 'D', 'E']

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        out_html = os.path.join(output_dir, "matrix_report_test.html")
        bundle_path = os.path.join(output_dir, "net.dvb")
        DataMatrix(data, row_names, row_names).write(input_file, output_type="hdf")
        matrix_report.main(["-i", input_file, "--html_out", out_html])
        build_ssn.main(["-i", input_file, "--dvb_out", bundle_path])
        report_stops = _page_data(out_html)['slider_stops']
        bundle_stops = ssn_bundle.slider_stops(load_bundle(bundle_path))

    keys = ('edge_index', 'threshold_label', 'threshold_value')
    assert [{k: stop[k] for k in keys} for stop in report_stops] == [
        {k: stop[k] for k in keys} for stop in bundle_stops
    ]
    # Only the report carries the threshold tables, so only the report's stops index
    # into them. A bundle's stops have no threshold_index at all.
    assert all(stop['threshold_index'] >= -1 for stop in report_stops)
    assert all('threshold_index' not in stop for stop in bundle_stops)
