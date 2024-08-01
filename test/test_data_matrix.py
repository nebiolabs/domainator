import warnings
warnings.filterwarnings("ignore", message="The value of the smallest subnormal for <class 'numpy.float64'> type is zero.")
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
