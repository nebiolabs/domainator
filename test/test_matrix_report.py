import warnings
warnings.filterwarnings("ignore", module='numpy')
import base64
import gzip
import io
import json
import re
from domainator import matrix_report
from domainator.data_matrix import DataMatrix, DenseDataMatrix, SparseDataMatrix, mst_knn_edge_counts_by_threshold
from domainator.matrix_report import MaxTree
import tempfile
import pytest
import numpy as np
import scipy.sparse
import os


def _embedded_report_payload(html_content):
    match = re.search(r'const EMBEDDED_REPORT_DATA_GZIP_BASE64 = "([^"]+)";', html_content)
    assert match is not None
    return json.loads(gzip.decompress(base64.b64decode(match.group(1))).decode("utf-8"))


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
        matrix_report.main(["-i", str(shared_datadir / input_file), "-o", out_txt, "--html", out_html])
        for fh in (out_html, out_txt):
            f_txt = open(fh).read()
            assert "Matrix Report" in f_txt
            assert "Min" in f_txt or "min" in f_txt.lower()
            assert "Max" in f_txt or "max" in f_txt.lower()
            # Check that we have nodes count (20 nodes in the test data)
            assert "20" in f_txt
        assert 'Projected MST_KNN edge counts' not in open(out_txt).read()
        assert 'mst-knn-k-slider' not in open(out_html).read()


def test_matrix_report_rejects_invalid_merge_impact_metric():
    matrix = DenseDataMatrix(np.array([[0.0, 1.0], [1.0, 0.0]]), ['A', 'B'], ['A', 'B'])

    with pytest.raises(ValueError, match="merge_impact_metric must be one of"):
        matrix_report.matrix_report(matrix, io.StringIO(), None, merge_impact_metric="bad")


# def test_matrix_report_empty_input(shared_datadir):
#     pass
#     #TODO



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
        matrix = DenseDataMatrix(data, row_names, row_names)
        
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
        matrix = DenseDataMatrix(data, row_names, row_names)
        
        with tempfile.TemporaryDirectory() as output_dir:
            input_file = os.path.join(output_dir, "test_matrix.hdf5")
            output_file = os.path.join(output_dir, "interactive_test.html")
            
            # Save the matrix
            matrix.write(input_file, output_type="dense")
            
            matrix_report.main([
                "-i", input_file,
                "--html", output_file,
                "--include_mst_knn"
            ])
            
            assert os.path.exists(output_file)
            
            with open(output_file, 'r') as f:
                content = f.read()
            
            # Check that the HTML output contains interactive dashboard elements
            assert 'Matrix Report' in content
            assert 'VIZ_DATA' in content
            assert 'Plotly' in content or 'plotly' in content
            assert 'mst-knn-k-slider' in content
            assert 'MST_KNN_COUNTS' in content
            assert 'Projected MST_KNN Edges' in content


def test_mst_knn_edge_counts_by_threshold_monotonic():
    data = np.array([
        [0, 10, 8, 7],
        [10, 0, 6, 5],
        [8, 6, 0, 9],
        [7, 5, 9, 0],
    ])
    row_names = ['A', 'B', 'C', 'D']
    matrix = DenseDataMatrix(data, row_names, row_names)
    tree = MaxTree(matrix)

    counts = mst_knn_edge_counts_by_threshold(matrix, tree, 3)

    assert counts.shape == (len(tree.mst_edges), 2)
    assert np.all(counts[:, 0] >= np.arange(1, len(tree.mst_edges) + 1))
    assert np.all(counts[:, 1] >= counts[:, 0])


def test_mst_knn_edge_counts_by_threshold_sparse_matches_dense_bruteforce():
    data = np.array([
        [0, 10, 8, 0, 4],
        [10, 0, 6, 6, 0],
        [8, 6, 0, 9, 9],
        [0, 6, 9, 0, 5],
        [4, 0, 9, 5, 0],
    ], dtype=float)
    labels = ['A', 'B', 'C', 'D', 'E']
    dense_matrix = DenseDataMatrix(data, labels, labels)
    sparse_matrix = SparseDataMatrix(scipy.sparse.csr_array(data), labels, labels)
    dense_tree = MaxTree(dense_matrix)
    sparse_tree = MaxTree(sparse_matrix)

    def brute_force_counts(array, tree, max_k):
        counts = np.zeros((len(tree.mst_edges), max_k - 1), dtype=int)
        mst_prefix_edges = set()
        thresholds = tree.edges_by_threshold[:, 1]

        for threshold_idx, mst_edge in enumerate(tree.mst_edges):
            source_idx, target_idx, _ = mst_edge
            mst_prefix_edges.add((source_idx, target_idx) if source_idx < target_idx else (target_idx, source_idx))
            threshold = thresholds[threshold_idx]

            edge_min_rank = {}
            for row_idx in range(array.shape[0]):
                row_scores = np.maximum(array[row_idx, :], array[:, row_idx])
                candidates = [
                    (target_idx, row_scores[target_idx])
                    for target_idx in range(array.shape[0])
                    if target_idx != row_idx and row_scores[target_idx] >= threshold and row_scores[target_idx] > 0
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

            counts[threshold_idx, :] = len(mst_prefix_edges) + np.cumsum(non_mst_rank_counts)[2:]

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
        matrix_report.main(["-i", str(shared_datadir / "scorefull.tsv"), "-o", out_txt, "--html", out_html, "--include_mst_knn"])

        text_content = open(out_txt).read()
        html_content = open(out_html).read()

        assert 'Projected MST_KNN edge counts' in text_content
        assert 'MST_KNN edges' in text_content
        assert 'mst-knn-k-slider' in html_content


def test_matrix_report_includes_merge_event_outputs():
    data = np.array([
        [0, 10, 6, 0],
        [10, 0, 7, 0],
        [6, 7, 0, 4],
        [0, 0, 4, 0],
    ], dtype=float)
    row_names = ['A', 'B', 'C', 'D']
    matrix = DenseDataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        out_html = os.path.join(output_dir, "matrix_report_test.html")
        out_txt = os.path.join(output_dir, "matrix_report_test.txt")

        matrix.write(input_file, output_type="dense")
        matrix_report.main([
            "-i", input_file,
            "-o", out_txt,
            "--html", out_html,
        ])

        text_content = open(out_txt).read()
        html_content = open(out_html).read()

        assert 'Strongest MST split events' in text_content
        assert '      From        To    Impact' in text_content
        assert 'Top-2 ratio' not in text_content
        assert 'Largest cluster changes (MST-derived)' not in text_content
        assert 'Cluster Splits vs Threshold' in html_content
        assert 'Size of smallest new cluster' in html_content
        assert "orientation: 'h'" in html_content
        assert "y: 1.08" in html_content
        assert 'Cumulative split fraction' not in html_content
        assert 'range: [-0.2, 1.05]' not in html_content
        assert "zerolinecolor: '#9ca3af'" not in html_content
        assert 'Largest Merge Events (min_child)' not in html_content
        assert 'merge-events-table' not in html_content
        assert 'MERGE_EVENT_SERIES' in html_content
        payload = _embedded_report_payload(html_content)
        assert len(payload['merge_event_series']) == 3
        assert payload['slider_stops'][0]['edge_index'] == -1
        assert [stop['edge_index'] for stop in payload['slider_stops'][1:]] == [0, 1, 2]
        assert 'Clusters and Edge Count vs Threshold' in html_content
        assert 'id="edges-by-threshold"' not in html_content
        assert '<div class="chart chart-wide">\n        <div id="cluster-discontinuity-by-threshold"></div>' in html_content
        assert 'cdn.jsdelivr.net/npm/d3@7' in html_content
        assert 'cluster-size-bubble' in html_content
        assert 'updateClusterBubbleChart' in html_content
        assert 'num-singletons' in html_content
        assert 'Number of Singletons' in html_content
        assert 'Cluster Size Summaries vs Threshold' not in html_content
        assert 'cluster-size-hist' not in html_content
        assert "''.join(" not in html_content
        assert "row['threshold_from']" not in html_content
        assert '{{' not in html_content


def test_matrix_report_supports_product_merge_impact_metric():
    data = np.array([
        [0, 10, 6, 0],
        [10, 0, 7, 0],
        [6, 7, 0, 4],
        [0, 0, 4, 0],
    ], dtype=float)
    row_names = ['A', 'B', 'C', 'D']
    matrix = DenseDataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        out_html = os.path.join(output_dir, "matrix_report_test.html")
        out_txt = os.path.join(output_dir, "matrix_report_test.txt")

        matrix.write(input_file, output_type="dense")
        matrix_report.main([
            "-i", input_file,
            "-o", out_txt,
            "--html", out_html,
            "--merge_impact_metric", "product",
        ])

        text_content = open(out_txt).read()
        html_content = open(out_html).read()

        assert 'Strongest MST split events (metric=product)' in text_content
        assert 'Clusters and Edge Count vs Threshold' in html_content
        assert 'Cluster Splits vs Threshold' in html_content
        assert 'Largest Merge Events (product)' not in html_content


def test_matrix_report_default_excludes_mst_knn(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out_html = output_dir + "/matrix_report_test.html"
        out_txt = output_dir + "/matrix_report_test.txt"
        matrix_report.main(["-i", str(shared_datadir / "scorefull.tsv"), "-o", out_txt, "--html", out_html])

        text_content = open(out_txt).read()
        html_content = open(out_html).read()

        assert 'Projected MST_KNN edge counts' not in text_content
        assert 'MST_KNN edges' not in text_content
        assert 'mst-knn-k-slider' not in html_content
        payload = _embedded_report_payload(html_content)
        assert payload['has_mst_knn'] is False
        assert payload['mst_knn_counts'] == []
        assert 'Cluster Size Summaries vs Threshold' not in html_content


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
    matrix = DenseDataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        out_html = os.path.join(output_dir, "matrix_report_test.html")

        matrix.write(input_file, output_type="dense")
        matrix_report.main([
            "-i", input_file,
            "--html", out_html,
            "--max_merge_events", "2",
        ])

        html_content = open(out_html).read()
        payload = _embedded_report_payload(html_content)

        assert len(payload['merge_event_series']) == 2
        assert [row['edge_index'] for row in payload['merge_event_series']] == [3, 4]
        assert [stop['edge_index'] for stop in payload['slider_stops']] == [-1, 3, 4]
        assert 'id="threshold-slider" min="0" max="2" value="0" step="1"' in html_content


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
    matrix = DenseDataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        out_html = os.path.join(output_dir, "matrix_report_test.html")

        matrix.write(input_file, output_type="dense")
        matrix_report.main([
            "-i", input_file,
            "--html", out_html,
            "--max_merge_events", "0",
        ])

        html_content = open(out_html).read()
        payload = _embedded_report_payload(html_content)

        assert len(payload['merge_event_series']) == 5
        assert [stop['edge_index'] for stop in payload['slider_stops']] == [-1, 0, 1, 2, 3, 4]
        assert 'id="threshold-slider" min="0" max="5" value="0" step="1"' in html_content


def test_matrix_report_text_merge_events_limited_to_top_25():
    n_nodes = 30
    data = np.zeros((n_nodes, n_nodes), dtype=float)
    for idx in range(n_nodes - 1):
        weight = float(n_nodes - idx)
        data[idx, idx + 1] = weight
        data[idx + 1, idx] = weight

    row_names = [f'N{idx}' for idx in range(n_nodes)]
    matrix = DenseDataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        out_txt = os.path.join(output_dir, "matrix_report_test.txt")

        matrix.write(input_file, output_type="dense")
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
    matrix = DenseDataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        out_txt = os.path.join(output_dir, "matrix_report_test.txt")

        matrix.write(input_file, output_type="dense")
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
    matrix = DenseDataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        out_html = os.path.join(output_dir, "matrix_report_test.html")

        matrix.write(input_file, output_type="dense")
        matrix_report.main([
            "-i", input_file,
            "--html", out_html,
            "--include_mst_knn",
            "--profile_stages",
        ])

        stderr_output = capsys.readouterr().err

        assert 'matrix_report stage timings:' in stderr_output
        assert 'load_matrix:' in stderr_output
        assert 'build_edge_table:' in stderr_output
        assert 'build_max_tree:' in stderr_output
        assert 'build_mst_knn_neighbor_rankings:' in stderr_output
        assert 'build_mst_knn_counts:' in stderr_output
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
    matrix = DenseDataMatrix(data, row_names, row_names)

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = os.path.join(output_dir, "test_matrix.hdf5")
        out_html = os.path.join(output_dir, "matrix_report_test.html")

        matrix.write(input_file, output_type="dense")
        matrix_report.main([
            "-i", input_file,
            "--html", out_html,
            "--progress",
        ])

        stderr_output = capsys.readouterr().err

        assert 'loading matrix from' in stderr_output
        assert 'building edge table' in stderr_output
        assert 'building maximum spanning tree' in stderr_output
        assert 'rendering outputs' in stderr_output


if __name__ == '__main__':
    pytest.main([__file__, '-v'])