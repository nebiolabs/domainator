import warnings
warnings.filterwarnings("ignore", message="The value of the smallest subnormal for <class 'numpy.float64'> type is zero.")
from domainator import transform_matrix
import tempfile
import pandas as pd
import numpy as np
import scipy
import scipy.sparse
from domainator.data_matrix import DataMatrix
from helpers import compare_iterables
import pytest

def compare_data_matrix(matrix_1:DataMatrix, matrix_2:DataMatrix):
    assert matrix_1.rows == matrix_2.rows
    assert matrix_1.columns == matrix_2.columns
    #assert np.testing.assert_array_equal(np.array(matrix_1.row_lengths), np.array(matrix_2.row_lengths))
    np.testing.assert_array_equal(matrix_1.row_lengths, matrix_2.row_lengths)
    np.testing.assert_array_equal(matrix_1.column_lengths, matrix_2.column_lengths)
    assert matrix_1.data_type == matrix_2.data_type


@pytest.mark.parametrize("output_mode,expected",
[("row_norm_score", np.array([
        [-1/3.,2/3.,1/3.,1,0],   
        [1/5.,2/5.,3/5.,4/5.,1],
        [1,4/5.,3/5.,2/5.,1/5.],
        [1/3.,1,1,1,2/3.],
        [1/3.,1/3.,1/3.,1,2/3.],
        ])
),
("norm_score", np.array([
     [1 + 1/3., 1/3.,   2/3.,  0.,    1.],   
     [4/5.,     1/2.,   0.,    0.,    0.],
     [0.,       0.,     0.,    1/2., 4/5.],
     [2/3.,     0.,     0.,    0.,   1/3.],
     [2/3.,     2/3.,   2/3.,  0.,   1/3.],
    ])
),
("efi_score", np.array([[-0., 29.19990958, 14.09725727, 44.15449935, -0., ],
        [13.90537175, 28.89887958, 43.89922684, 58.90496914, 73.91507624],
        [73.93527962, 58.82578789, 43.72313559, 28.62587831, 13.53298584],
        [13.60434175, 43.64934937, 43.59819685, 43.55243936, 28.45954689],
        [13.50743174, 13.44943979, 13.39828727, 43.45552935, 28.36263688]
    ])
)
])
def test_transform_matrix_1(output_mode, expected, shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        input = np.array([
        [-1,2,1,3,0],   
        [1,2,3,4,5],
        [5,4,3,2,1],
        [1,3,3,3,2],
        [1,1,1,3,2],
        ])

        input_file = output_dir + "/input.hdf5"
        dense_out = output_dir + "/dense_out.hdf5"
        sparse_out = output_dir + "/sparse_out.hdf5"

        rows =  ["a", "b", "c", "d", "e"]
        columns = ["A", "B", "C", "D", "E"]
        row_lengths = [1, 2, 3, 4, 5]
        col_lengths = [7, 8, 9, 10, 11]

        matrix = DataMatrix(input, rows, columns, row_lengths, col_lengths, "score")
        expected_matrix = DataMatrix(expected, rows, columns, row_lengths, col_lengths, output_mode)
        matrix.write(input_file, "dense")
        transform_matrix.main(["-i", input_file, "--dense", dense_out, "--sparse", sparse_out, "--mode", output_mode])
        dense_output_matrix = DataMatrix.from_file(dense_out)
        sparse_output_matrix = DataMatrix.from_file(sparse_out)

        compare_data_matrix(dense_output_matrix, expected_matrix)
        expected_matrix.convert_to_sparse()
        compare_data_matrix(sparse_output_matrix, expected_matrix)


def test_transform_matrix_2(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        input = np.array([
        [-1,2,1,3,0],   
        [1,2,3,4,5],
        [5,4,3,2,1],
        [1,3,3,3,2],
        [1,1,1,3,2],
        ])

        (output_mode, expected) = ("score_dist", np.array([
            [1 + 1/3., 1/3.,   2/3.,  0.,    1.],   
            [4/5.,     1/2.,   0.,    0.,    0.],
            [0.,       0.,     0.,    1/2., 4/5.],
            [2/3.,     0.,     0.,    0.,   1/3.],
            [2/3.,     2/3.,   2/3.,  0.,   1/3.],
            ]))

        input_file = output_dir + "/input.hdf5"
        dense_out = output_dir + "/dense_out.hdf5"

        rows =  ["a", "b", "c", "d", "e"]
        columns = ["A", "B", "C", "D", "E"]
        row_lengths = [1, 2, 3, 4, 5]
        col_lengths = [7, 8, 9, 10, 11]

        matrix = DataMatrix(input, rows, columns, row_lengths, col_lengths, "score")
        expected_matrix = DataMatrix(expected, rows, columns, row_lengths, col_lengths, output_mode)
        matrix.write(input_file, "dense")
        transform_matrix.main(["-i", input_file, "--dense", dense_out, "--mode", output_mode])
        dense_output_matrix = DataMatrix.from_file(dense_out)

        compare_data_matrix(dense_output_matrix, expected_matrix)


#TODO: test more conversion modes