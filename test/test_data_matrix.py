import warnings
warnings.filterwarnings("ignore", module='numpy')
import pytest
from domainator.data_matrix import DataMatrix
import scipy.sparse
import numpy as np
import pytest_datadir

# Test initialization of DataMatrix
def test_init():
    # Test case 1: Initialize with data
    data = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    row_names = ['A', 'B', 'C']
    col_names = ['X', 'Y', 'Z']
    matrix = DataMatrix(data, row_names, col_names)
    assert matrix.shape == (3, 3)
    assert matrix.size == 9
    assert not matrix.sparse
    assert matrix.rows == row_names
    assert matrix.columns == col_names
    assert matrix.row_lengths is None
    assert matrix.column_lengths is None
    assert matrix.data_type == ""

    # Test case 2: Initialize without data
    matrix = DataMatrix()
    assert matrix.shape == (0, 0)
    assert matrix.size == 0
    assert not matrix.sparse
    assert matrix.rows is None
    assert matrix.columns is None
    assert matrix.row_lengths is None
    assert matrix.column_lengths is None
    assert matrix.data_type == ""

# Test from_file method of DataMatrix

@pytest.mark.parametrize("filename,sparse",
[
    ("test_matrix.dense.hdf5",False),
    ("test_matrix.dense.tsv",False),
    ("test_matrix.sparse.hdf5",True)
])
def test_from_file(shared_datadir, filename, sparse):
    # Test case 1: Read dense matrix from file
    matrix_file = shared_datadir / filename
    matrix = DataMatrix.from_file(matrix_file)
    assert matrix.shape == (3, 3)
    assert matrix.size == 9
    assert matrix.sparse is sparse
    assert matrix.rows == ['A', 'B', 'C']
    assert matrix.columns == ['X', 'Y', 'Z']
    assert matrix.row_lengths is None
    assert matrix.column_lengths is None
    assert matrix.data_type == ""


# Test convert_to_sparse method of DataMatrix
def test_convert_to_sparse():
    # Test case 1: Convert dense matrix to sparse
    data = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    row_names = ['A', 'B', 'C']
    col_names = ['X', 'Y', 'Z']
    matrix = DataMatrix(data, row_names, col_names)
    matrix.convert_to_sparse()
    assert matrix.sparse
    assert matrix.data.shape == (3, 3)

    # Test case 2: Convert already sparse matrix to sparse
    matrix = DataMatrix()
    matrix.sparse = True
    matrix.data = scipy.sparse.csr_matrix([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
    matrix.convert_to_sparse()
    assert matrix.sparse
    assert matrix.data.shape == (3, 3)

# Test iter_data order and zeros for synthetic data
def test_iter_data_order_and_zeros_fake():
    """
    Ensure iter_data for sparse and dense matrices never iterates over zeros and always matches order (synthetic data).
    """
    data = np.array([[1, 0, 2], [0, 3, 0], [4, 0, 5]])
    row_names = ['A', 'B', 'C']
    col_names = ['X', 'Y', 'Z']

    dense_matrix = DataMatrix(data, row_names, col_names)
    sparse_matrix = DataMatrix(scipy.sparse.csr_array(data), row_names, col_names)

    dense_iter = [(r, c, v) for r, c, v in dense_matrix.iter_data() if v != 0]
    sparse_iter = [(r, c, v) for r, c, v in sparse_matrix.iter_data()]

    # Check that sparse never iterates over zeros
    assert all(v != 0 for r, c, v in sparse_iter), "Sparse matrix should not iterate over zeros"
    # Check that dense (filtered) and sparse have same order and values
    assert dense_iter == sparse_iter, f"Dense and sparse iter_data results differ: {dense_iter} vs {sparse_iter}"

# Test iter_data order and zeros for real FeSOD files
def test_iter_data_order_and_zeros_real(shared_datadir):
    """
    Ensure iter_data for sparse and dense matrices never iterates over zeros and always matches order (real FeSOD files).
    """
    dense_file = shared_datadir / "FeSOD_dist.dense.hdf5"
    sparse_file = shared_datadir / "FeSOD_dist.sparse.hdf5"

    dense_matrix = DataMatrix.from_file(dense_file)
    sparse_matrix = DataMatrix.from_file(sparse_file)

    dense_iter = [(r, c, v) for r, c, v in dense_matrix.iter_data() if v != 0]
    sparse_iter = [(r, c, v) for r, c, v in sparse_matrix.iter_data()]

    # Check that sparse never iterates over zeros
    assert all(v != 0 for r, c, v in sparse_iter), "Sparse matrix should not iterate over zeros"
    # Check that dense (filtered) and sparse have same order and values
    assert dense_iter == sparse_iter, f"Dense and sparse iter_data results differ: {dense_iter[:10]} vs {sparse_iter[:10]}"


# Test triangular method with dense matrix
def test_triangular_dense():
    """Test triangular method with a dense matrix."""
    data = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    row_names = ['A', 'B', 'C']
    col_names = ['A', 'B', 'C']
    matrix = DataMatrix(data, row_names, col_names)
    
    # Test lower triangular without diagonal
    result = matrix.triangular(side='lower', include_diagonal=False, skip_zeros=True)
    expected = [('B', 'A', 4), ('C', 'A', 7), ('C', 'B', 8)]
    assert result == expected, f"Lower triangular (no diagonal) failed: {result}"
    
    # Test lower triangular with diagonal
    result = matrix.triangular(side='lower', include_diagonal=True, skip_zeros=True)
    expected = [('A', 'A', 1), ('B', 'A', 4), ('B', 'B', 5), ('C', 'A', 7), ('C', 'B', 8), ('C', 'C', 9)]
    assert result == expected, f"Lower triangular (with diagonal) failed: {result}"
    
    # Test upper triangular without diagonal
    result = matrix.triangular(side='upper', include_diagonal=False, skip_zeros=True)
    expected = [('A', 'B', 2), ('A', 'C', 3), ('B', 'C', 6)]
    assert result == expected, f"Upper triangular (no diagonal) failed: {result}"
    
    # Test upper triangular with diagonal
    result = matrix.triangular(side='upper', include_diagonal=True, skip_zeros=True)
    expected = [('A', 'A', 1), ('A', 'B', 2), ('A', 'C', 3), ('B', 'B', 5), ('B', 'C', 6), ('C', 'C', 9)]
    assert result == expected, f"Upper triangular (with diagonal) failed: {result}"


# Test triangular method with sparse matrix
def test_triangular_sparse():
    """Test triangular method with a sparse matrix."""
    data = np.array([[1, 0, 3], [0, 5, 0], [7, 0, 9]])
    row_names = ['A', 'B', 'C']
    col_names = ['A', 'B', 'C']
    matrix = DataMatrix(scipy.sparse.csr_array(data), row_names, col_names)
    
    # Test lower triangular without diagonal, skip zeros
    result = matrix.triangular(side='lower', include_diagonal=False, skip_zeros=True)
    expected = [('C', 'A', 7)]
    assert result == expected, f"Sparse lower triangular (no diagonal, skip zeros) failed: {result}"
    
    # Test lower triangular with diagonal, skip zeros
    result = matrix.triangular(side='lower', include_diagonal=True, skip_zeros=True)
    expected = [('A', 'A', 1), ('B', 'B', 5), ('C', 'A', 7), ('C', 'C', 9)]
    assert result == expected, f"Sparse lower triangular (with diagonal, skip zeros) failed: {result}"
    
    # Test upper triangular with diagonal, skip zeros
    result = matrix.triangular(side='upper', include_diagonal=True, skip_zeros=True)
    expected = [('A', 'A', 1), ('A', 'C', 3), ('B', 'B', 5), ('C', 'C', 9)]
    assert result == expected, f"Sparse upper triangular (with diagonal, skip zeros) failed: {result}"


# Test triangular method with skip_zeros=False
def test_triangular_include_zeros():
    """Test that triangular correctly includes zeros when skip_zeros=False."""
    data = np.array([[1, 0, 3], [0, 5, 0], [7, 0, 9]])
    row_names = ['A', 'B', 'C']
    col_names = ['A', 'B', 'C']
    
    dense_matrix = DataMatrix(data, row_names, col_names)
    sparse_matrix = DataMatrix(scipy.sparse.csr_array(data), row_names, col_names)
    
    # Test lower triangular with diagonal, including zeros
    dense_result = dense_matrix.triangular(side='lower', include_diagonal=True, skip_zeros=False)
    sparse_result = sparse_matrix.triangular(side='lower', include_diagonal=True, skip_zeros=False)
    
    expected = [
        ('A', 'A', 1),
        ('B', 'A', 0), ('B', 'B', 5),
        ('C', 'A', 7), ('C', 'B', 0), ('C', 'C', 9)
    ]
    
    assert len(dense_result) == 6, f"Dense result should have 6 items, got {len(dense_result)}"
    assert len(sparse_result) == 6, f"Sparse result should have 6 items, got {len(sparse_result)}"
    assert dense_result == expected, f"Dense result mismatch: {dense_result}"
    assert sparse_result == expected, f"Sparse result mismatch: {sparse_result}"


# Test triangular method with index_style parameter
def test_triangular_index_style():
    """Test that triangular correctly returns indices when index_style='index'."""
    data = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    row_names = ['A', 'B', 'C']
    col_names = ['A', 'B', 'C']
    matrix = DataMatrix(data, row_names, col_names)
    
    # Test with index_style='name' (default)
    result_name = matrix.triangular(side='lower', include_diagonal=True, skip_zeros=True, index_style='name')
    assert result_name[0] == ('A', 'A', 1)
    assert result_name[1] == ('B', 'A', 4)
    
    # Test with index_style='index'
    result_index = matrix.triangular(side='lower', include_diagonal=True, skip_zeros=True, index_style='index')
    assert result_index[0] == (0, 0, 1)
    assert result_index[1] == (1, 0, 4)
    
    # Verify that the values match
    assert len(result_name) == len(result_index)
    for (name_r, name_c, name_v), (idx_r, idx_c, idx_v) in zip(result_name, result_index):
        assert row_names[idx_r] == name_r
        assert col_names[idx_c] == name_c
        assert idx_v == name_v


# Test triangular consistency between dense and sparse
def test_triangular_dense_sparse_consistency():
    """Test that dense and sparse matrices return the same triangular results."""
    data = np.array([[1, 2, 3, 4], 
                     [5, 6, 7, 8],
                     [9, 10, 11, 12],
                     [13, 14, 15, 16]])
    row_names = ['A', 'B', 'C', 'D']
    col_names = ['A', 'B', 'C', 'D']
    
    dense = DataMatrix(data, row_names, col_names)
    sparse = DataMatrix(scipy.sparse.csr_array(data), row_names, col_names)
    
    # Test all combinations
    for side in ['lower', 'upper']:
        for include_diagonal in [True, False]:
            for skip_zeros in [True, False]:
                dense_result = dense.triangular(side=side, include_diagonal=include_diagonal, skip_zeros=skip_zeros)
                sparse_result = sparse.triangular(side=side, include_diagonal=include_diagonal, skip_zeros=skip_zeros)
                
                assert dense_result == sparse_result, \
                    f"Dense and sparse differ for side={side}, diagonal={include_diagonal}, skip_zeros={skip_zeros}"


# Test triangular with non-square matrix
def test_triangular_non_square_error():
    """Test that triangular raises an error for non-square matrices."""
    data = np.array([[1, 2, 3], [4, 5, 6]])
    row_names = ['A', 'B']
    col_names = ['X', 'Y', 'Z']
    matrix = DataMatrix(data, row_names, col_names)
    
    with pytest.raises(NotImplementedError, match="Triangular only supported for square matrices"):
        matrix.triangular()


# Test triangular with invalid parameters
def test_triangular_invalid_parameters():
    """Test that triangular raises errors for invalid parameters."""
    data = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    row_names = ['A', 'B', 'C']
    col_names = ['A', 'B', 'C']
    matrix = DataMatrix(data, row_names, col_names)
    
    # Test invalid side
    with pytest.raises(ValueError, match="side must be 'lower' or 'upper'"):
        matrix.triangular(side='middle')
    
    # Test invalid index_style
    with pytest.raises(ValueError, match="index_style must be 'name' or 'index'"):
        matrix.triangular(index_style='invalid')


# Test triangular sorting order
def test_triangular_sorting():
    """Test that triangular returns results sorted by row then column."""
    data = np.array([[1, 2, 3, 4], 
                     [5, 6, 7, 8],
                     [9, 10, 11, 12],
                     [13, 14, 15, 16]])
    row_names = ['A', 'B', 'C', 'D']
    col_names = ['A', 'B', 'C', 'D']
    matrix = DataMatrix(data, row_names, col_names)
    
    # Test lower triangular with diagonal
    result = matrix.triangular(side='lower', include_diagonal=True, skip_zeros=True, index_style='index')
    
    # Verify sorted order (row first, then column)
    for i in range(len(result) - 1):
        curr_row, curr_col, _ = result[i]
        next_row, next_col, _ = result[i + 1]
        
        # Next item should have row >= current row
        assert next_row >= curr_row, f"Row order violated at index {i}: {result[i]} -> {result[i+1]}"
        
        # If same row, next column should be > current column
        if next_row == curr_row:
            assert next_col > curr_col, f"Column order violated at index {i}: {result[i]} -> {result[i+1]}"


# Test symmetric_values with dense symmetric matrix
def test_symmetric_values_dense_symmetric():
    """Test symmetric_values returns True for a symmetric dense matrix."""
    data = np.array([[1, 2, 3], [2, 5, 6], [3, 6, 9]])
    row_names = ['A', 'B', 'C']
    col_names = ['A', 'B', 'C']
    matrix = DataMatrix(data, row_names, col_names)
    
    assert matrix.symmetric_values is True, "Symmetric matrix should return True"


# Test symmetric_values with dense non-symmetric matrix
def test_symmetric_values_dense_non_symmetric():
    """Test symmetric_values returns False for a non-symmetric dense matrix."""
    data = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    row_names = ['A', 'B', 'C']
    col_names = ['A', 'B', 'C']
    matrix = DataMatrix(data, row_names, col_names)
    
    assert matrix.symmetric_values is False, "Non-symmetric matrix should return False"


# Test symmetric_values with sparse symmetric matrix
def test_symmetric_values_sparse_symmetric():
    """Test symmetric_values returns True for a symmetric sparse matrix."""
    data = np.array([[1, 0, 3], [0, 5, 0], [3, 0, 9]])
    row_names = ['A', 'B', 'C']
    col_names = ['A', 'B', 'C']
    matrix = DataMatrix(scipy.sparse.csr_array(data), row_names, col_names)
    
    assert matrix.symmetric_values is True, "Symmetric sparse matrix should return True"


# Test symmetric_values with sparse non-symmetric matrix (different values)
def test_symmetric_values_sparse_non_symmetric_values():
    """Test symmetric_values returns False for sparse matrix with asymmetric values."""
    data = np.array([[1, 0, 3], [0, 5, 0], [7, 0, 9]])
    row_names = ['A', 'B', 'C']
    col_names = ['A', 'B', 'C']
    matrix = DataMatrix(scipy.sparse.csr_array(data), row_names, col_names)
    
    assert matrix.symmetric_values is False, "Sparse matrix with asymmetric values should return False"


# Test symmetric_values with sparse non-symmetric matrix (different structure)
def test_symmetric_values_sparse_non_symmetric_structure():
    """Test symmetric_values returns False for sparse matrix with asymmetric structure."""
    data = np.array([[1, 2, 3], [0, 5, 0], [0, 0, 9]])
    row_names = ['A', 'B', 'C']
    col_names = ['A', 'B', 'C']
    matrix = DataMatrix(scipy.sparse.csr_array(data), row_names, col_names)
    
    assert matrix.symmetric_values is False, "Sparse matrix with asymmetric structure should return False"


# Test symmetric_values with floating point imprecision
def test_symmetric_values_floating_point_tolerance():
    """Test that symmetric_values handles floating point imprecision correctly."""
    # Create a matrix that is symmetric within floating point tolerance
    data = np.array([[1.0, 2.5, 3.7], 
                     [2.5, 5.0, 6.3], 
                     [3.7, 6.3 + 1e-10, 9.0]])
    row_names = ['A', 'B', 'C']
    col_names = ['A', 'B', 'C']
    matrix = DataMatrix(data, row_names, col_names)
    
    assert matrix.symmetric_values is True, "Matrix with tiny floating point errors should be considered symmetric"
    
    # Create a matrix that is NOT symmetric beyond tolerance
    data = np.array([[1.0, 2.5, 3.7], 
                     [2.5, 5.0, 6.3], 
                     [3.7, 6.3 + 0.01, 9.0]])
    matrix = DataMatrix(data, row_names, col_names)
    
    assert matrix.symmetric_values is False, "Matrix with significant differences should not be considered symmetric"


# Test symmetric_values with non-symmetric labels
def test_symmetric_values_non_symmetric_labels():
    """Test that symmetric_values returns False immediately for non-symmetric labels."""
    data = np.array([[1, 2, 3], [2, 5, 6], [3, 6, 9]])
    row_names = ['A', 'B', 'C']
    col_names = ['X', 'Y', 'Z']
    matrix = DataMatrix(data, row_names, col_names)
    
    assert matrix.symmetric_values is False, "Matrix with non-symmetric labels should return False"


# Test symmetric_values consistency between dense and sparse
def test_symmetric_values_dense_sparse_consistency():
    """Test that dense and sparse matrices return the same symmetric_values result."""
    # Test with symmetric data
    data_symmetric = np.array([[1, 2, 3], [2, 5, 6], [3, 6, 9]])
    row_names = ['A', 'B', 'C']
    col_names = ['A', 'B', 'C']
    
    dense = DataMatrix(data_symmetric, row_names, col_names)
    sparse = DataMatrix(scipy.sparse.csr_array(data_symmetric), row_names, col_names)
    
    assert dense.symmetric_values == sparse.symmetric_values == True, \
        "Dense and sparse symmetric matrices should both return True"
    
    # Test with non-symmetric data
    data_non_symmetric = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    
    dense = DataMatrix(data_non_symmetric, row_names, col_names)
    sparse = DataMatrix(scipy.sparse.csr_array(data_non_symmetric), row_names, col_names)
    
    assert dense.symmetric_values == sparse.symmetric_values == False, \
        "Dense and sparse non-symmetric matrices should both return False"


# Test triangular with agg parameter - dense matrix
def test_triangular_agg_dense():
    """Test that triangular correctly applies aggregation function with dense matrices."""
    data = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    row_names = ['A', 'B', 'C']
    col_names = ['A', 'B', 'C']
    matrix = DataMatrix(data, row_names, col_names)
    
    # Test with sum aggregation
    result = matrix.triangular(side='lower', include_diagonal=False, skip_zeros=True, agg=lambda a, b: a + b)
    # Lower [1,0] = 4, transpose [0,1] = 2, sum = 6
    # Lower [2,0] = 7, transpose [0,2] = 3, sum = 10
    # Lower [2,1] = 8, transpose [1,2] = 6, sum = 14
    expected = [('B', 'A', 6), ('C', 'A', 10), ('C', 'B', 14)]
    assert result == expected, f"Sum aggregation failed: {result}"
    
    # Test with max aggregation
    result = matrix.triangular(side='upper', include_diagonal=False, skip_zeros=True, agg=max)
    # Upper [0,1] = 2, transpose [1,0] = 4, max = 4
    # Upper [0,2] = 3, transpose [2,0] = 7, max = 7
    # Upper [1,2] = 6, transpose [2,1] = 8, max = 8
    expected = [('A', 'B', 4), ('A', 'C', 7), ('B', 'C', 8)]
    assert result == expected, f"Max aggregation failed: {result}"
    
    # Test with min aggregation
    result = matrix.triangular(side='lower', include_diagonal=False, skip_zeros=True, agg=min)
    # Lower [1,0] = 4, transpose [0,1] = 2, min = 2
    # Lower [2,0] = 7, transpose [0,2] = 3, min = 3
    # Lower [2,1] = 8, transpose [1,2] = 6, min = 6
    expected = [('B', 'A', 2), ('C', 'A', 3), ('C', 'B', 6)]
    assert result == expected, f"Min aggregation failed: {result}"


# Test triangular with agg parameter - sparse matrix
def test_triangular_agg_sparse():
    """Test that triangular correctly applies aggregation function with sparse matrices."""
    data = np.array([[1, 0, 3], [0, 5, 0], [7, 0, 9]])
    row_names = ['A', 'B', 'C']
    col_names = ['A', 'B', 'C']
    matrix = DataMatrix(scipy.sparse.csr_array(data), row_names, col_names)
    
    # Test with sum aggregation
    result = matrix.triangular(side='lower', include_diagonal=False, skip_zeros=True, agg=lambda a, b: a + b)
    # Lower [2,0] = 7, transpose [0,2] = 3, sum = 10
    expected = [('C', 'A', 10)]
    assert result == expected, f"Sparse sum aggregation failed: {result}"
    
    # Test with max aggregation
    result = matrix.triangular(side='upper', include_diagonal=False, skip_zeros=True, agg=max)
    # Upper [0,2] = 3, transpose [2,0] = 7, max = 7
    expected = [('A', 'C', 7)]
    assert result == expected, f"Sparse max aggregation failed: {result}"


# Test triangular with agg parameter - averaging for symmetric matrix
def test_triangular_agg_average():
    """Test that triangular with average aggregation works correctly."""
    data = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    row_names = ['A', 'B', 'C']
    col_names = ['A', 'B', 'C']
    matrix = DataMatrix(data, row_names, col_names)
    
    # Test with average aggregation
    result = matrix.triangular(side='upper', include_diagonal=False, skip_zeros=True, 
                               agg=lambda a, b: (a + b) / 2)
    # Upper [0,1] = 2, transpose [1,0] = 4, avg = 3.0
    # Upper [0,2] = 3, transpose [2,0] = 7, avg = 5.0
    # Upper [1,2] = 6, transpose [2,1] = 8, avg = 7.0
    expected = [('A', 'B', 3.0), ('A', 'C', 5.0), ('B', 'C', 7.0)]
    assert result == expected, f"Average aggregation failed: {result}"


# Test triangular agg consistency between dense and sparse
def test_triangular_agg_dense_sparse_consistency():
    """Test that agg parameter works consistently for dense and sparse matrices."""
    data = np.array([[1, 2, 3, 4], 
                     [5, 6, 7, 8],
                     [9, 10, 11, 12],
                     [13, 14, 15, 16]])
    row_names = ['A', 'B', 'C', 'D']
    col_names = ['A', 'B', 'C', 'D']
    
    dense = DataMatrix(data, row_names, col_names)
    sparse = DataMatrix(scipy.sparse.csr_array(data), row_names, col_names)
    
    # Test with sum aggregation
    dense_result = dense.triangular(side='lower', include_diagonal=False, skip_zeros=True, 
                                    agg=lambda a, b: a + b)
    sparse_result = sparse.triangular(side='lower', include_diagonal=False, skip_zeros=True, 
                                      agg=lambda a, b: a + b)
    
    assert dense_result == sparse_result, \
        "Dense and sparse matrices with agg should produce same results"


# Test triangular agg with skip_zeros edge case
def test_triangular_agg_skip_zeros_edge_case():
    """Test that skip_zeros with agg only skips when BOTH values are zero."""
    data = np.array([[1, 2, 0], 
                     [0, 5, 6], 
                     [0, 0, 9]])
    row_names = ['A', 'B', 'C']
    col_names = ['A', 'B', 'C']
    matrix = DataMatrix(data, row_names, col_names)
    
    # Without agg: skips when value is zero
    result_no_agg = matrix.triangular(side='lower', include_diagonal=False, skip_zeros=True)
    # [1,0]=0 skipped, [2,0]=0 skipped, [2,1]=0 skipped
    # No items
    assert len(result_no_agg) == 0, f"Without agg should skip all zeros: {result_no_agg}"
    
    # With agg: only skips when BOTH value and symmetric complement are zero
    result_with_agg = matrix.triangular(side='lower', include_diagonal=False, skip_zeros=True, 
                                       agg=lambda a, b: a + b)
    # [1,0]=0, transpose [0,1]=2, both zero? No -> include, sum=2
    # [2,0]=0, transpose [0,2]=0, both zero? Yes -> SKIP
    # [2,1]=0, transpose [1,2]=6, both zero? No -> include, sum=6
    expected = [('B', 'A', 2), ('C', 'B', 6)]
    assert result_with_agg == expected, f"With agg should include when one is non-zero: {result_with_agg}"
    
    # Test upper triangular as well
    result_upper = matrix.triangular(side='upper', include_diagonal=False, skip_zeros=True, 
                                     agg=lambda a, b: a + b)
    # [0,1]=2, transpose [1,0]=0, both zero? No -> include, sum=2
    # [0,2]=0, transpose [2,0]=0, both zero? Yes -> SKIP
    # [1,2]=6, transpose [2,1]=0, both zero? No -> include, sum=6
    expected = [('A', 'B', 2), ('B', 'C', 6)]
    assert result_upper == expected, f"Upper triangular with agg edge case failed: {result_upper}"


# Test triangular agg with skip_zeros edge case sparse
def test_triangular_agg_skip_zeros_edge_case():
    """Test that skip_zeros with agg only skips when BOTH values are zero."""
    data = np.array([[1, 2, 0], 
                     [0, 5, 6], 
                     [0, 0, 9]])
    row_names = ['A', 'B', 'C']
    col_names = ['A', 'B', 'C']
    matrix = DataMatrix(data, row_names, col_names)

    matrix.convert_to_sparse()
    
    # Without agg: skips when value is zero
    result_no_agg = matrix.triangular(side='lower', include_diagonal=False, skip_zeros=True)
    # [1,0]=0 skipped, [2,0]=0 skipped, [2,1]=0 skipped
    # No items
    assert len(result_no_agg) == 0, f"Without agg should skip all zeros: {result_no_agg}"
    
    # With agg: only skips when BOTH value and symmetric complement are zero
    result_with_agg = matrix.triangular(side='lower', include_diagonal=False, skip_zeros=True, 
                                       agg=lambda a, b: a + b)
    # [1,0]=0, transpose [0,1]=2, both zero? No -> include, sum=2
    # [2,0]=0, transpose [0,2]=0, both zero? Yes -> SKIP
    # [2,1]=0, transpose [1,2]=6, both zero? No -> include, sum=6
    expected = [('B', 'A', 2), ('C', 'B', 6)]
    assert result_with_agg == expected, f"With agg should include when one is non-zero: {result_with_agg}"
    
    # Test upper triangular as well
    result_upper = matrix.triangular(side='upper', include_diagonal=False, skip_zeros=True, 
                                     agg=lambda a, b: a + b)
    # [0,1]=2, transpose [1,0]=0, both zero? No -> include, sum=2
    # [0,2]=0, transpose [2,0]=0, both zero? Yes -> SKIP
    # [1,2]=6, transpose [2,1]=0, both zero? No -> include, sum=6
    expected = [('A', 'B', 2), ('B', 'C', 6)]
    assert result_upper == expected, f"Upper triangular with agg edge case failed: {result_upper}"
