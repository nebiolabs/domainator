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