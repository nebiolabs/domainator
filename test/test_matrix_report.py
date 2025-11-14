import warnings
warnings.filterwarnings("ignore", module='numpy')
from domainator import matrix_report
from domainator.data_matrix import DataMatrix
from domainator.matrix_report import MaxTree
import tempfile
import pytest
import numpy as np
import scipy.sparse
import os


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


# def test_matrix_report_empty_input(shared_datadir):
#     pass
#     #TODO


class TestMaxTree:
    """Test suite for the MaxTree class"""
    
    def test_init_simple_triangle(self):
        """Test MST construction on a simple 3-node complete graph"""
        # Create a simple symmetric matrix
        # 3 nodes with edges: (0,1)=10, (0,2)=5, (1,2)=8
        data = np.array([
            [0, 10, 5],
            [10, 0, 8],
            [5, 8, 0]
        ])
        row_names = ['A', 'B', 'C']
        matrix = DataMatrix(data, row_names, row_names)
        
        tree = MaxTree(matrix)
        
        # Should have n-1 = 2 edges in MST
        assert len(tree.mst) == 2
        assert tree.n_nodes == 3
        
        # MST should contain edges with values 10 and 8 (highest two)
        mst_values = [tree.edges[tree.mst[i], 2] for i in range(len(tree.mst))]
        assert sorted(mst_values, reverse=True) == [10, 8]
    
    def test_init_four_nodes(self):
        """Test MST construction on a 4-node complete graph"""
        # 4 nodes in a square with diagonals
        data = np.array([
            [0, 10, 5, 3],
            [10, 0, 4, 6],
            [5, 4, 0, 9],
            [3, 6, 9, 0]
        ])
        row_names = ['A', 'B', 'C', 'D']
        matrix = DataMatrix(data, row_names, row_names)
        
        tree = MaxTree(matrix)
        
        # Should have n-1 = 3 edges in MST
        assert len(tree.mst) == 3
        assert tree.n_nodes == 4
        
        # MST should contain edges with values 10, 9, 6 (highest three without cycles)
        mst_values = sorted([tree.edges[tree.mst[i], 2] for i in range(len(tree.mst))], reverse=True)
        assert mst_values == [10, 9, 6]
    
    def test_init_disconnected_components(self):
        """Test MST on a matrix with zeros (disconnected if we filter)"""
        # Two separate triangles
        data = np.array([
            [0, 10, 5, 0, 0],
            [10, 0, 8, 0, 0],
            [5, 8, 0, 0, 0],
            [0, 0, 0, 0, 7],
            [0, 0, 0, 7, 0]
        ])
        row_names = ['A', 'B', 'C', 'D', 'E']
        matrix = DataMatrix(data, row_names, row_names)
        
        tree = MaxTree(matrix)
        
        # It's not a connected graph, so the "mst" is a multigraph, and not actually an mst
        assert len(tree.mst) == 3
        assert tree.n_nodes == 5
    
    def test_init_no_cycles(self):
        """Verify that MST doesn't create cycles"""
        # Create a graph where naive algorithm would create cycle
        data = np.array([
            [0, 10, 10, 0],
            [10, 0, 10, 0],
            [10, 10, 0, 10],
            [0, 0, 10, 0]
        ])
        row_names = ['A', 'B', 'C', 'D']
        matrix = DataMatrix(data, row_names, row_names)
        
        tree = MaxTree(matrix)
        
        # Should have exactly n-1 = 3 edges
        assert len(tree.mst) == 3
        
        # Verify it's a valid tree by checking connectivity using union-find
        parent = np.arange(4)
        
        def find(x):
            if parent[x] != x:
                parent[x] = find(parent[x])
            return parent[x]
        
        for i in range(len(tree.mst)):
            edge_idx = tree.mst[i]
            node1 = int(tree.edges[edge_idx, 0])
            node2 = int(tree.edges[edge_idx, 1])
            
            # These two nodes should not already be in the same component
            # (otherwise we'd have a cycle)
            root1 = find(node1)
            root2 = find(node2)
            
            # Union them
            parent[root1] = root2
        
        # All nodes should be in same component
        roots = [find(i) for i in range(4)]
        assert len(set(roots)) == 1
    
    def test_edges_by_threshold_simple(self):
        """Test edges_by_threshold on simple graph"""
        data = np.array([
            [0, 10, 5],
            [10, 0, 8],
            [5, 8, 0]
        ])
        row_names = ['A', 'B', 'C']
        matrix = DataMatrix(data, row_names, row_names)
        
        tree = MaxTree(matrix)
        result = tree.edges_by_threshold
        
        # Should have 2 rows (one per MST edge)
        assert result.shape[0] == 2
        assert result.shape[1] == 2
        
        # Check thresholds are in descending order
        thresholds = result[:, 1]
        assert all(thresholds[i] >= thresholds[i+1] for i in range(len(thresholds)-1))
    
    def test_edges_by_threshold_counts(self):
        """Test that edge counts in edges_by_threshold are cumulative"""
        data = np.array([
            [0, 10, 5, 3],
            [10, 0, 4, 6],
            [5, 4, 0, 9],
            [3, 6, 9, 0]
        ])
        row_names = ['A', 'B', 'C', 'D']
        matrix = DataMatrix(data, row_names, row_names)
        
        tree = MaxTree(matrix)
        result = tree.edges_by_threshold
        
        # Edge counts should be increasing
        edge_counts = result[:, 0]
        assert all(edge_counts[i] <= edge_counts[i+1] for i in range(len(edge_counts)-1))
        
        # Last edge count should be total number of edges
        assert edge_counts[-1] > 0
    
    def test_cluster_count_by_threshold_simple(self):
        """Test cluster_count_by_threshold on simple graph"""
        data = np.array([
            [0, 10, 5],
            [10, 0, 8],
            [5, 8, 0]
        ])
        row_names = ['A', 'B', 'C']
        matrix = DataMatrix(data, row_names, row_names)
        
        tree = MaxTree(matrix)
        result = tree.cluster_count_by_threshold
        
        # Should have n rows (one per MST edge plus initial state)
        assert result.shape[0] == 3  # n_nodes = 3, so 2 MST edges + 1 initial = 3
        assert result.shape[1] == 2
        
        # First row: infinite threshold, n clusters
        assert result[0, 0] == float('inf')
        assert result[0, 1] == 3
        
        # Last row: should have 1 cluster (fully connected)
        assert result[-1, 1] == 1
        
        # Cluster counts should be decreasing
        cluster_counts = result[:, 1]
        assert all(cluster_counts[i] >= cluster_counts[i+1] for i in range(len(cluster_counts)-1))
        
        # Thresholds should be decreasing (after the first inf)
        thresholds = result[1:, 0]
        assert all(thresholds[i] >= thresholds[i+1] for i in range(len(thresholds)-1))
    
    def test_cluster_count_by_threshold_four_nodes(self):
        """Test cluster_count_by_threshold on 4-node graph"""
        data = np.array([
            [0, 10, 5, 3],
            [10, 0, 4, 6],
            [5, 4, 0, 9],
            [3, 6, 9, 0]
        ])
        row_names = ['A', 'B', 'C', 'D']
        matrix = DataMatrix(data, row_names, row_names)
        
        tree = MaxTree(matrix)
        result = tree.cluster_count_by_threshold
        
        # Should have 5 rows (4 nodes -> 3 MST edges + 1 initial)
        assert result.shape[0] == 4
        
        # Check cluster count progression: 4 -> 3 -> 2 -> 1
        expected_clusters = [4, 3, 2, 1]
        actual_clusters = result[:, 1].astype(int).tolist()
        assert actual_clusters == expected_clusters
    
    def test_cluster_count_by_edge_count_simple(self):
        """Test cluster_count_by_edge_count on simple graph"""
        data = np.array([
            [0, 10, 5],
            [10, 0, 8],
            [5, 8, 0]
        ])
        row_names = ['A', 'B', 'C']
        matrix = DataMatrix(data, row_names, row_names)
        
        tree = MaxTree(matrix)
        result = tree.cluster_count_by_edge_count
        
        # Should have at least 3 rows (0 edges + 2 MST edges)
        assert result.shape[0] >= 3
        assert result.shape[1] == 2
        
        # First row: 0 edges, n clusters
        assert result[0, 0] == 0
        assert result[0, 1] == 3
        
        # Last row: should have 1 cluster
        assert result[-1, 1] == 1
        
        # Edge counts should be increasing
        edge_counts = result[:, 0]
        assert all(edge_counts[i] <= edge_counts[i+1] for i in range(len(edge_counts)-1))
        
        # Cluster counts should be decreasing when MST edges are added
        cluster_counts = result[:, 1]
        assert cluster_counts[0] == 3  # Start with 3 clusters
        assert cluster_counts[-1] == 1  # End with 1 cluster
    
    def test_cluster_count_by_edge_count_relationship(self):
        """Test that cluster_count_by_edge_count maintains n_nodes - n_components relationship"""
        data = np.array([
            [0, 10, 5, 3],
            [10, 0, 4, 6],
            [5, 4, 0, 9],
            [3, 6, 9, 0]
        ])
        row_names = ['A', 'B', 'C', 'D']
        matrix = DataMatrix(data, row_names, row_names)
        
        tree = MaxTree(matrix)
        result = tree.cluster_count_by_edge_count
        
        # For each MST edge added, cluster count should decrease by 1
        # Starting from n_nodes clusters at 0 edges
        assert result[0, 1] == 4  # 4 nodes, 0 edges -> 4 clusters
        
        # Find rows corresponding to MST edges
        edges_by_thresh = tree.edges_by_threshold
        mst_edge_counts = edges_by_thresh[:, 0]
        
        for i, mst_edge_count in enumerate(mst_edge_counts):
            # Find this edge count in result
            matching_rows = np.where(result[:, 0] == mst_edge_count)[0]
            if len(matching_rows) > 0:
                cluster_count = result[matching_rows[-1], 1]
                # After i+1 MST edges, should have n_nodes - (i+1) clusters
                expected_clusters = 4 - (i + 1)
                assert cluster_count == expected_clusters
    
    def test_single_node(self):
        """Test MaxTree with single node (edge case)"""
        data = np.array([[0]])
        row_names = ['A']
        matrix = DataMatrix(data, row_names, row_names)
        
        tree = MaxTree(matrix)
        
        # Single node: MST should be empty
        assert len(tree.mst) == 0
        assert tree.n_nodes == 1
    
    def test_two_nodes(self):
        """Test MaxTree with two nodes"""
        data = np.array([
            [0, 5],
            [5, 0]
        ])
        row_names = ['A', 'B']
        matrix = DataMatrix(data, row_names, row_names)
        
        tree = MaxTree(matrix)
        
        # Two nodes: MST should have 1 edge
        assert len(tree.mst) == 1
        assert tree.n_nodes == 2
        
        # Check that edge value is 5
        assert tree.edges[tree.mst[0], 2] == 5
        
        # Check cluster counts
        result = tree.cluster_count_by_threshold
        assert result.shape[0] == 2  # Initial state + 1 edge
        assert result[0, 1] == 2  # Start with 2 clusters
        assert result[1, 1] == 1  # End with 1 cluster
    
    def test_sparse_matrix(self):
        """Test MaxTree with sparse matrix"""
        # Create sparse matrix
        data = np.array([
            [0, 10, 0, 0],
            [10, 0, 8, 0],
            [0, 8, 0, 5],
            [0, 0, 5, 0]
        ])
        row_names = ['A', 'B', 'C', 'D']
        sparse_data = scipy.sparse.csr_array(data)
        matrix = DataMatrix(sparse_data, row_names, row_names)
        
        tree = MaxTree(matrix)
        
        # Should still construct valid MST
        assert len(tree.mst) == 3
        assert tree.n_nodes == 4
        
        # MST should only use non-zero edges
        mst_values = [tree.edges[tree.mst[i], 2] for i in range(len(tree.mst))]
        assert all(v > 0 for v in mst_values)
    
    def test_all_equal_weights(self):
        """Test MaxTree when all edge weights are equal"""
        data = np.array([
            [0, 5, 5],
            [5, 0, 5],
            [5, 5, 0]
        ])
        row_names = ['A', 'B', 'C']
        matrix = DataMatrix(data, row_names, row_names)
        
        tree = MaxTree(matrix)
        
        # Should still construct valid MST with 2 edges
        assert len(tree.mst) == 2
        assert tree.n_nodes == 3
        
        # All MST edges should have value 5
        mst_values = [tree.edges[tree.mst[i], 2] for i in range(len(tree.mst))]
        assert all(v == 5 for v in mst_values)
    
    def test_mst_values_sorted(self):
        """Test that MST edges are processed in descending value order"""
        data = np.array([
            [0, 10, 5, 3, 8],
            [10, 0, 4, 6, 2],
            [5, 4, 0, 9, 7],
            [3, 6, 9, 0, 1],
            [8, 2, 7, 1, 0]
        ])
        row_names = ['A', 'B', 'C', 'D', 'E']
        matrix = DataMatrix(data, row_names, row_names)
        
        tree = MaxTree(matrix)
        
        # Get MST edge values
        mst_values = [tree.edges[tree.mst[i], 2] for i in range(len(tree.mst))]
        
        # MST values should be in descending order (processed highest first)
        assert all(mst_values[i] >= mst_values[i+1] for i in range(len(mst_values)-1))


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
            matrix.write(input_file, output_type="dense")
            
            matrix_report.main([
                "-i", input_file,
                "--html", output_file
            ])
            
            assert os.path.exists(output_file)
            
            with open(output_file, 'r') as f:
                content = f.read()
            
            # Check that the HTML output contains interactive dashboard elements
            assert 'Matrix Report' in content
            assert 'VIZ_DATA' in content
            assert 'Plotly' in content or 'plotly' in content


if __name__ == '__main__':
    pytest.main([__file__, '-v'])