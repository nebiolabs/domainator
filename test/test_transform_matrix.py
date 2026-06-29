import warnings
warnings.filterwarnings("ignore", module='numpy')
from domainator import transform_matrix
import tempfile
import pandas as pd
import numpy as np
import scipy
import scipy.sparse
from domainator.data_matrix import DataMatrix, DenseDataMatrix
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

        matrix = DenseDataMatrix(input, rows, columns, row_lengths, col_lengths, "score")
        expected_matrix = DenseDataMatrix(expected, rows, columns, row_lengths, col_lengths, output_mode)
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

        matrix = DenseDataMatrix(input, rows, columns, row_lengths, col_lengths, "score")
        expected_matrix = DenseDataMatrix(expected, rows, columns, row_lengths, col_lengths, output_mode)
        matrix.write(input_file, "dense")
        transform_matrix.main(["-i", input_file, "--dense", dense_out, "--mode", output_mode])
        dense_output_matrix = DataMatrix.from_file(dense_out)

        compare_data_matrix(dense_output_matrix, expected_matrix)


def test_transform_matrix_default_mode_preserves_values():
    with tempfile.TemporaryDirectory() as output_dir:
        input_data = np.array([
            [5.0, 1.0, 0.0],
            [2.0, 4.0, 3.0],
            [0.0, 6.0, 7.0],
        ])

        input_file = output_dir + "/input.hdf5"
        dense_out = output_dir + "/dense_out.hdf5"

        rows = ["a", "b", "c"]
        columns = ["A", "B", "C"]

        matrix = DenseDataMatrix(input_data, rows, columns, data_type="bool")
        matrix.write(input_file, "dense")

        transform_matrix.main(["-i", input_file, "--dense", dense_out])
        dense_output_matrix = DataMatrix.from_file(dense_out)

        np.testing.assert_array_equal(dense_output_matrix.data, input_data)
        assert dense_output_matrix.data_type == "bool"


def test_transform_matrix_lb_applied_after_mode_transform():
    with tempfile.TemporaryDirectory() as output_dir:
        input_data = np.array([
            [0.0, 0.2, 0.0],
            [0.3, 0.0, 0.4],
            [0.0, 0.0, 0.0],
        ])

        input_file = output_dir + "/input.hdf5"
        dense_out = output_dir + "/dense_out.hdf5"

        rows = ["a", "b", "c"]
        columns = ["A", "B", "C"]

        matrix = DenseDataMatrix(input_data, rows, columns, data_type="score")
        matrix.write(input_file, "dense")

        transform_matrix.main(["-i", input_file, "--dense", dense_out, "--mode", "bool", "--lb", "0.5"])
        dense_output_matrix = DataMatrix.from_file(dense_out)

        expected = np.array([
            [0, 1, 0],
            [1, 0, 1],
            [0, 0, 0],
        ], dtype=np.int8)

        np.testing.assert_array_equal(dense_output_matrix.data, expected)
        assert dense_output_matrix.data_type == "bool"


def test_transform_matrix_mst_knn_applied_after_mode_transform():
    with tempfile.TemporaryDirectory() as output_dir:
        input_data = np.array([
            [9.0, 2.0, 1.0],
            [2.0, 9.0, 8.0],
            [1.0, 8.0, 9.0],
        ])

        input_file = output_dir + "/input.hdf5"
        dense_out = output_dir + "/dense_out.hdf5"

        rows = ["a", "b", "c"]
        columns = ["a", "b", "c"]

        matrix = DenseDataMatrix(input_data, rows, columns, data_type="score")
        matrix.write(input_file, "dense")

        transform_matrix.main(["-i", input_file, "--dense", dense_out, "--mode", "bool", "--mst_knn", "2"])
        dense_output_matrix = DataMatrix.from_file(dense_out)

        expected = np.ones((3, 3), dtype=np.int8)
        np.testing.assert_array_equal(dense_output_matrix.data, expected)
        assert dense_output_matrix.data_type == "bool"


def _feed_matrix_entries(acc, data, order="forward"):
    """Feed every nonzero entry of a dense matrix into a StreamingMstKnnAccumulator."""
    entries = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            value = data[i, j]
            if value != 0:
                entries.append((i, j, float(value)))
    if order == "reversed":
        entries = entries[::-1]
    elif order == "chunked":
        # interleave to exercise the buffer flushing in a non-monotone order
        entries = entries[0::2] + entries[1::2]
    for i, j, value in entries:
        acc.add_edge(i, j, value)
    return acc


@pytest.mark.parametrize("order", ["forward", "reversed", "chunked"])
@pytest.mark.parametrize("buffer_cap,soft_cap", [(None, None), (1, 1), (2, 2)])
def test_streaming_mst_knn_matches_batch(order, buffer_cap, soft_cap):
    from domainator.data_matrix import SparseDataMatrix, StreamingMstKnnAccumulator
    from domainator.transform_matrix import apply_mst_knn_sparsification

    # symmetric matrix with all-distinct off-diagonal weights => unique MSF, no tie ambiguity.
    data = np.array([
        [5.0, 9.0, 1.0, 0.0, 2.0, 0.0],
        [9.0, 5.0, 8.0, 3.0, 0.0, 0.0],
        [1.0, 8.0, 5.0, 7.0, 0.0, 4.0],
        [0.0, 3.0, 7.0, 5.0, 6.0, 0.0],
        [2.0, 0.0, 0.0, 6.0, 5.0, 10.0],
        [0.0, 0.0, 4.0, 0.0, 10.0, 5.0],
    ])
    labels = ["a", "b", "c", "d", "e", "f"]
    k = 2

    matrix = SparseDataMatrix(scipy.sparse.csr_array(data), labels, labels, data_type="score")
    expected = apply_mst_knn_sparsification(matrix, k, lower_bound=0)

    acc = StreamingMstKnnAccumulator(len(labels), k, lower_bound=0,
                                     mst_buffer_cap=buffer_cap, knn_soft_cap=soft_cap)
    _feed_matrix_entries(acc, data, order=order)
    result = acc.to_csr()

    np.testing.assert_array_equal(result.toarray(), expected.toarray())


def test_streaming_mst_knn_respects_lower_bound():
    from domainator.data_matrix import SparseDataMatrix, StreamingMstKnnAccumulator
    from domainator.transform_matrix import apply_mst_knn_sparsification

    data = np.array([
        [0.0, 9.0, 1.0, 2.0],
        [9.0, 0.0, 8.0, 3.0],
        [1.0, 8.0, 0.0, 7.0],
        [2.0, 3.0, 7.0, 0.0],
    ])
    labels = ["a", "b", "c", "d"]
    k = 2
    lb = 2.0

    matrix = SparseDataMatrix(scipy.sparse.csr_array(data), labels, labels, data_type="score")
    expected = apply_mst_knn_sparsification(matrix, k, lower_bound=lb)

    acc = StreamingMstKnnAccumulator(len(labels), k, lower_bound=lb)
    _feed_matrix_entries(acc, data)
    result = acc.to_csr()

    np.testing.assert_array_equal(result.toarray(), expected.toarray())


def test_streaming_mst_knn_tie_tolerant_is_valid_forest():
    # Repeated weights => MSF edge selection may differ; assert it is still a valid
    # maximum spanning forest (same total weight + connectivity) and identical kNN edges.
    from domainator.data_matrix import SparseDataMatrix, StreamingMstKnnAccumulator
    from domainator.transform_matrix import apply_mst_knn_sparsification
    from scipy.sparse.csgraph import connected_components

    data = np.array([
        [0.0, 5.0, 5.0, 0.0],
        [5.0, 0.0, 5.0, 0.0],
        [5.0, 5.0, 0.0, 5.0],
        [0.0, 0.0, 5.0, 0.0],
    ])
    labels = ["a", "b", "c", "d"]
    k = 2

    matrix = SparseDataMatrix(scipy.sparse.csr_array(data), labels, labels, data_type="score")
    expected = apply_mst_knn_sparsification(matrix, k, lower_bound=0).toarray()

    acc = StreamingMstKnnAccumulator(len(labels), k, lower_bound=0)
    _feed_matrix_entries(acc, data)
    result = acc.to_csr().toarray()

    # same total retained weight and same connectivity
    assert np.isclose(result.sum(), expected.sum())
    n_exp = connected_components(expected > 0, directed=False)[0]
    n_res = connected_components(result > 0, directed=False)[0]
    assert n_exp == n_res


def test_transform_matrix_max_output_gb_blocks_dense_output():
    with tempfile.TemporaryDirectory() as output_dir:
        input_data = np.array([
            [5.0, 1.0, 0.0],
            [2.0, 4.0, 3.0],
            [0.0, 6.0, 7.0],
        ])

        input_file = output_dir + "/input.hdf5"
        dense_out = output_dir + "/dense_out.hdf5"

        rows = ["a", "b", "c"]
        columns = ["A", "B", "C"]

        matrix = DenseDataMatrix(input_data, rows, columns, data_type="score")
        matrix.write(input_file, "dense")

        with pytest.raises(SystemExit, match="--max_output_gb"):
            transform_matrix.main(["-i", input_file, "--dense", dense_out, "--max_output_gb", "0.000001"])


def test_transform_matrix_output_format_conversions():
    with tempfile.TemporaryDirectory() as output_dir:
        input_data = np.array([
            [5.0, 1.0, 0.0],
            [2.0, 4.0, 3.0],
            [0.0, 6.0, 7.0],
        ])

        input_file = output_dir + "/input.hdf5"
        dense_text_out = output_dir + "/dense_out.tsv"
        sparse_out = output_dir + "/sparse_out.hdf5"

        rows = ["a", "b", "c"]
        columns = ["A", "B", "C"]

        matrix = DenseDataMatrix(input_data, rows, columns, data_type="score")
        matrix.write(input_file, "dense")

        transform_matrix.main([
            "-i", input_file,
            "--dense_text", dense_text_out,
            "--sparse", sparse_out,
            "--mode", "bool",
        ])

        dense_text_matrix = DataMatrix.from_file(dense_text_out)
        sparse_matrix = DataMatrix.from_file(sparse_out)
        expected = np.array([
            [1, 1, 0],
            [1, 1, 1],
            [0, 1, 1],
        ], dtype=np.int8)

        np.testing.assert_array_equal(dense_text_matrix.data, expected)
        np.testing.assert_array_equal(sparse_matrix.data.toarray(), expected)
        assert sparse_matrix.data_type == "bool"