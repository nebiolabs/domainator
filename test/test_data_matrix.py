import warnings
warnings.filterwarnings("ignore", module='numpy')
import pytest
from domainator.edge_types import CompactNeighborRankings, SortedUndirectedEdges
from domainator.data_matrix import DataMatrix
from domainator.ssn_edges import MaxTree, build_symmetric_neighbor_rankings
import scipy.sparse
import numpy as np
import pytest_datadir

# Test initialization of DataMatrix
def test_init():
    # Test case 1: Initialize DataMatrix with data
    data = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    row_names = ['A', 'B', 'C']
    col_names = ['X', 'Y', 'Z']
    matrix = DataMatrix(data, row_names, col_names)
    assert matrix.shape == (3, 3)
    assert matrix.size == 9
    assert isinstance(matrix, DataMatrix)
    assert matrix.rows == row_names
    assert matrix.columns == col_names
    assert matrix.row_lengths is None
    assert matrix.column_lengths is None
    assert matrix.data_type == ""

    # Test case 2: Initialize DataMatrix with data
    sparse_data = scipy.sparse.csr_array(data)
    matrix = DataMatrix(sparse_data, row_names, col_names)
    assert matrix.shape == (3, 3)
    assert matrix.size == 9
    assert isinstance(matrix, DataMatrix)
    assert matrix.rows == row_names
    assert matrix.columns == col_names
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
    assert isinstance(matrix.data, scipy.sparse.csr_array)
    assert matrix.rows == ['A', 'B', 'C']
    assert matrix.columns == ['X', 'Y', 'Z']
    assert matrix.row_lengths is None
    assert matrix.column_lengths is None
    assert matrix.data_type == ""


def test_dense_text_reader_rejects_duplicate_row_labels(tmp_path):
    matrix_file = tmp_path / "duplicate_rows.tsv"
    matrix_file.write_text("\tA\tB\nrow1\t1\t2\nrow1\t3\t4\n")

    with pytest.raises(ValueError, match="row_names must be unique"):
        DataMatrix.from_file(matrix_file)


def test_dense_text_reader_rejects_duplicate_column_labels(tmp_path):
    matrix_file = tmp_path / "duplicate_columns.tsv"
    matrix_file.write_text("\tA\tA\nrow1\t1\t2\nrow2\t3\t4\n")

    with pytest.raises(ValueError, match="col_names must be unique"):
        DataMatrix.from_file(matrix_file)


def test_dense_text_reader_rejects_non_numeric_values(tmp_path):
    matrix_file = tmp_path / "non_numeric.tsv"
    matrix_file.write_text("\tA\tB\nrow1\t1\tnope\nrow2\t3\t4\n")

    with pytest.raises(ValueError, match="Dense text matrix values must be numeric"):
        DataMatrix.from_file(matrix_file)


def test_validate_matrix_data_keeps_prefix_labels_as_asymmetric(tmp_path):
    data = np.array([
        [1.0, 2.0, 3.0],
        [4.0, 5.0, 6.0],
    ])
    row_names = ['A', 'B']
    col_names = ['A', 'B', 'C']

    validated_rows, validated_cols = DataMatrix.validate_matrix_data(data, row_names, col_names)
    assert validated_rows == row_names
    assert validated_cols == col_names
    assert validated_cols is not validated_rows

    out_file = tmp_path / "prefix_labels.dense.hdf5"
    DataMatrix.write_hdf(data, out_file, row_names, col_names)
    matrix = DataMatrix.from_file(out_file)

    assert matrix.rows == row_names
    assert matrix.columns == col_names
    assert matrix.symmetric_labels is False


@pytest.mark.parametrize(
    "row_names,col_names,error_message",
    [
        (['A', 'A'], ['X', 'Y'], "row_names must be unique"),
        (['A', 'B'], ['X', 'X'], "col_names must be unique"),
    ],
)
def test_validate_matrix_data_rejects_duplicate_labels(row_names, col_names, error_message):
    data = np.ones((2, 2))

    with pytest.raises(ValueError, match=error_message):
        DataMatrix.validate_matrix_data(data, row_names, col_names)


def test_write_dense_validates_column_lengths():
    data = np.ones((2, 2))

    with pytest.raises(ValueError, match="col_lens size must match"):
        DataMatrix.write_hdf(
            data,
            "/tmp/unused_domainator_matrix.hdf5",
            ['A', 'B'],
            ['A', 'B'],
            row_lengths=np.array([100, 200]),
            col_lengths=np.array([100]),
        )


@pytest.mark.parametrize("input_kind", ["numpy", "csr"])
def test_hdf5_roundtrip_preserves_asymmetric_metadata(tmp_path, input_kind):
    data = np.array([
        [1.0, 0.0, 3.0],
        [0.0, 5.0, 6.0],
    ])
    row_names = ['row1', 'row2']
    col_names = ['col1', 'col2', 'col3']
    row_lengths = np.array([100, 200])
    col_lengths = np.array([10, 20, 30])
    out_file = tmp_path / f"metadata.{input_kind}.hdf5"

    matrix_data = scipy.sparse.csr_array(data) if input_kind == "csr" else data
    DataMatrix.write_hdf(matrix_data, out_file, row_names, col_names, row_lengths, col_lengths, data_type="score")

    matrix = DataMatrix.from_file(out_file)

    assert matrix.rows == row_names
    assert matrix.columns == col_names
    assert np.array_equal(matrix.row_lengths, row_lengths)
    assert np.array_equal(matrix.column_lengths, col_lengths)
    assert matrix.data_type == "score"
    assert np.array_equal(matrix.toarray(), data)


@pytest.mark.parametrize("input_kind", ["numpy", "csr"])
def test_hdf5_roundtrip_preserves_symmetric_lengths(tmp_path, input_kind):
    data = np.array([
        [1.0, 2.0],
        [2.0, 3.0],
    ])
    names = ['A', 'B']
    lengths = np.array([100, 200])
    out_file = tmp_path / f"symmetric_lengths.{input_kind}.hdf5"

    matrix_data = scipy.sparse.csr_array(data) if input_kind == "csr" else data
    # Equal copies of the lengths, not the same object: the labels are still symmetric.
    DataMatrix.write_hdf(matrix_data, out_file, names, list(names), lengths, lengths.copy(), data_type="efi_score")

    matrix = DataMatrix.from_file(out_file)

    assert matrix.rows == names
    assert matrix.columns == names
    assert np.array_equal(matrix.row_lengths, lengths)
    assert np.array_equal(matrix.column_lengths, lengths)
    assert matrix.data_type == "efi_score"
    assert matrix.symmetric_labels is True


def test_sparse_hdf5_roundtrip_preserves_zero_only_rows_and_columns(tmp_path):
    data = scipy.sparse.csr_array(np.array([
        [0.0, 0.0, 0.0],
        [0.0, 5.0, 0.0],
        [0.0, 0.0, 0.0],
    ]))
    row_names = ['row1', 'row2', 'row3']
    col_names = ['col1', 'col2', 'col3']
    out_file = tmp_path / "zero_only_axes.sparse.hdf5"

    DataMatrix.write_hdf(data, out_file, row_names, col_names)
    matrix = DataMatrix.from_file(out_file)

    assert isinstance(matrix, DataMatrix)
    assert matrix.shape == (3, 3)
    assert matrix.rows == row_names
    assert matrix.columns == col_names
    assert np.array_equal(matrix.toarray(), data.toarray())


def test_copy_data_false_reuses_csr_storage():
    data = scipy.sparse.csr_array(np.array([[1.0, 2.0], [3.0, 4.0]]))
    copied = DataMatrix(data, ['A', 'B'], ['A', 'B'])
    shared = DataMatrix(data, ['A', 'B'], ['A', 'B'], copy_data=False)

    data.data[0] = 99.0

    assert copied.data[0, 0] == 1.0
    assert shared.data[0, 0] == 99.0


def test_inputs_become_canonical_csr():
    # Dense input, other sparse formats, and explicit zeros all end up as canonical CSR.
    data = np.array([[1.0, 0.0, 3.0], [0.0, 0.0, 6.0]])
    explicit_zero = scipy.sparse.csr_array((np.array([0.0, 1.0]), np.array([1, 0]), np.array([0, 2, 2])), shape=(2, 3))
    for source in (data, scipy.sparse.coo_array(data), scipy.sparse.dok_array(data), explicit_zero):
        matrix = DataMatrix(source, ['A', 'B'], ['X', 'Y', 'Z'])
        assert isinstance(matrix.data, scipy.sparse.csr_array)
        assert matrix.data.has_sorted_indices
        assert not np.any(matrix.data.data == 0)


def test_assigning_data_canonicalizes_and_checks_shape():
    matrix = DataMatrix(np.eye(2), ['A', 'B'], ['A', 'B'])
    assert matrix.exactly_symmetric
    matrix.data = scipy.sparse.dok_array(np.array([[0.0, 2.0], [0.0, 1.0]]))
    assert isinstance(matrix.data, scipy.sparse.csr_array)
    assert not matrix.exactly_symmetric  # the cached answer is reset
    with pytest.raises(ValueError, match="does not match the labels"):
        matrix.data = np.zeros((3, 3))


def test_sorted_undirected_edges_sparse_matches_dense():
    data = np.array([
        [0, 10, 0, 7],
        [8, 0, 6, 0],
        [4, 9, 0, 5],
        [7, 0, 5, 0],
    ], dtype=float)
    row_names = ['A', 'B', 'C', 'D']

    dense_matrix = DataMatrix(data, row_names, row_names)
    sparse_matrix = DataMatrix(scipy.sparse.csr_array(data), row_names, row_names)

    dense_edges = dense_matrix.sorted_undirected_edges(skip_zeros=True, agg=max)
    sparse_edges = sparse_matrix.sorted_undirected_edges(skip_zeros=True, agg=max)

    assert dense_edges.n_nodes == sparse_edges.n_nodes == len(row_names)
    assert np.array_equal(dense_edges.source, sparse_edges.source)
    assert np.array_equal(dense_edges.target, sparse_edges.target)
    assert np.array_equal(dense_edges.score, sparse_edges.score)


def test_build_symmetric_neighbor_rankings_sparse_matches_dense():
    data = np.array([
        [0, 10, 0, 7],
        [8, 0, 6, 0],
        [4, 9, 0, 5],
        [7, 0, 5, 0],
    ], dtype=float)
    row_names = ['A', 'B', 'C', 'D']

    dense_matrix = DataMatrix(data, row_names, row_names)
    sparse_matrix = DataMatrix(scipy.sparse.csr_array(data), row_names, row_names)

    dense_rankings = build_symmetric_neighbor_rankings(dense_matrix, max_k=3)
    sparse_rankings = build_symmetric_neighbor_rankings(sparse_matrix, max_k=3)

    assert np.array_equal(dense_rankings.offsets, sparse_rankings.offsets)
    assert np.array_equal(dense_rankings.target, sparse_rankings.target)
    assert np.array_equal(dense_rankings.score, sparse_rankings.score)


def test_build_symmetric_neighbor_rankings_top_k_tie_order():
    data = np.array([
        [0, 10, 10, 8, 7],
        [10, 0, 5, 0, 0],
        [10, 5, 0, 0, 0],
        [8, 0, 0, 0, 0],
        [7, 0, 0, 0, 0],
    ], dtype=float)
    row_names = ['A', 'B', 'C', 'D', 'E']
    matrix = DataMatrix(data, row_names, row_names)

    rankings = build_symmetric_neighbor_rankings(matrix, max_k=3)
    start, end = rankings.row_bounds(0)

    assert np.array_equal(rankings.target[start:end], np.array([1, 2, 3], dtype=np.int32))
    assert np.allclose(rankings.score[start:end], np.array([10.0, 10.0, 8.0]))


@pytest.mark.parametrize("max_k", [1, 2, 5, 40])
def test_top_k_neighbor_rankings_are_the_full_rankings_cut_at_k(max_k):
    """The top-k kernel keeps each row's first k neighbors of the full ranking, ties included."""
    rng = np.random.default_rng(7)
    n = 30
    upper = np.triu(rng.integers(0, 6, size=(n, n)).astype(float), k=1)  # few distinct scores: many ties
    upper[rng.random((n, n)) < 0.4] = 0
    data = upper + upper.T
    names = [str(i) for i in range(n)]
    matrix = DataMatrix(data, names, names)

    full = build_symmetric_neighbor_rankings(matrix)
    top_k = build_symmetric_neighbor_rankings(matrix, max_k=max_k)
    for row in range(n):
        full_start, full_end = full.row_bounds(row)
        start, end = top_k.row_bounds(row)
        keep = min(max_k, full_end - full_start)
        assert np.array_equal(top_k.target[start:end], full.target[full_start:full_start + keep])
        assert np.array_equal(top_k.score[start:end], full.score[full_start:full_start + keep])


@pytest.mark.parametrize(
    "kwargs,error_message",
    [
        ({"n_nodes": -1, "source": np.array([], dtype=np.int32), "target": np.array([], dtype=np.int32), "score": np.array([], dtype=float)}, "n_nodes must be >= 0"),
        ({"n_nodes": 2, "source": np.array([0], dtype=np.int32), "target": np.array([2], dtype=np.int32), "score": np.array([1.0])}, "less than n_nodes"),
        ({"n_nodes": 2, "source": np.array([1], dtype=np.int32), "target": np.array([1], dtype=np.int32), "score": np.array([1.0])}, "self edges"),
        ({"n_nodes": 2, "source": np.array([[0]], dtype=np.int32), "target": np.array([1], dtype=np.int32), "score": np.array([1.0])}, "1-dimensional"),
        ({"n_nodes": 2, "source": np.array([0, 0], dtype=np.int32), "target": np.array([1], dtype=np.int32), "score": np.array([1.0])}, "same length"),
        ({"n_nodes": 2, "source": np.array([-1], dtype=np.int32), "target": np.array([1], dtype=np.int32), "score": np.array([1.0])}, ">= 0"),
        ({"n_nodes": 3, "source": np.array([0, 1], dtype=np.int32), "target": np.array([1, 2], dtype=np.int32), "score": np.array([1.0, 2.0])}, "descending"),
    ],
)
def test_sorted_undirected_edges_validates_invariants(kwargs, error_message):
    with pytest.raises(ValueError, match=error_message):
        SortedUndirectedEdges(**kwargs)


@pytest.mark.parametrize(
    "kwargs,error_message",
    [
        ({"offsets": np.array([], dtype=np.int64), "target": np.array([], dtype=np.int32), "score": np.array([], dtype=float)}, "at least one"),
        ({"offsets": np.array([1], dtype=np.int64), "target": np.array([], dtype=np.int32), "score": np.array([], dtype=float)}, "start at 0"),
        ({"offsets": np.array([0, 2, 1], dtype=np.int64), "target": np.array([1], dtype=np.int32), "score": np.array([1.0])}, "monotonically"),
        ({"offsets": np.array([0, 1], dtype=np.int64), "target": np.array([1], dtype=np.int32), "score": np.array([1.0])}, "less than the number of rows"),
        ({"offsets": np.array([[0]], dtype=np.int64), "target": np.array([], dtype=np.int32), "score": np.array([], dtype=float)}, "1-dimensional"),
        ({"offsets": np.array([0, -1, 0], dtype=np.int64), "target": np.array([], dtype=np.int32), "score": np.array([], dtype=float)}, "non-negative"),
        ({"offsets": np.array([0, 2], dtype=np.int64), "target": np.array([0], dtype=np.int32), "score": np.array([1.0])}, "end at the length"),
        ({"offsets": np.array([0, 1], dtype=np.int64), "target": np.array([-1], dtype=np.int32), "score": np.array([1.0])}, "target indices must be >= 0"),
    ],
)
def test_compact_neighbor_rankings_validates_invariants(kwargs, error_message):
    with pytest.raises(ValueError, match=error_message):
        CompactNeighborRankings(**kwargs)


def test_max_tree_filters_zero_edges_from_sorted_edge_table():
    edges = SortedUndirectedEdges(
        n_nodes=3,
        source=np.array([0, 1], dtype=np.int32),
        target=np.array([1, 2], dtype=np.int32),
        score=np.array([5.0, 0.0]),
    )

    with pytest.warns(RuntimeWarning, match="removed zero-score edges"):
        tree = MaxTree(edges, skip_zeros=True)

    assert tree.mst_edges == [(0, 1, 5.0)]


def test_max_tree_warns_when_zero_edges_are_included():
    edges = SortedUndirectedEdges(
        n_nodes=3,
        source=np.array([0, 1], dtype=np.int32),
        target=np.array([1, 2], dtype=np.int32),
        score=np.array([5.0, 0.0]),
    )

    with pytest.warns(RuntimeWarning, match="includes zero-score edges"):
        tree = MaxTree(edges, skip_zeros=False)

    assert tree.mst_edges == [(0, 1, 5.0), (1, 2, 0.0)]


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
    result = matrix.triangular(include_diagonal=False, skip_zeros=True)
    expected = [('B', 'A', 4), ('C', 'A', 7), ('C', 'B', 8)]
    assert result == expected, f"Lower triangular (no diagonal) failed: {result}"
    
    # Test lower triangular with diagonal
    result = matrix.triangular(include_diagonal=True, skip_zeros=True)
    expected = [('A', 'A', 1), ('B', 'A', 4), ('B', 'B', 5), ('C', 'A', 7), ('C', 'B', 8), ('C', 'C', 9)]
    assert result == expected, f"Lower triangular (with diagonal) failed: {result}"


# Test triangular method with sparse matrix
def test_triangular_sparse():
    """Test triangular method with a sparse matrix."""
    data = np.array([[1, 0, 3], [0, 5, 0], [7, 0, 9]])
    row_names = ['A', 'B', 'C']
    col_names = ['A', 'B', 'C']
    matrix = DataMatrix(scipy.sparse.csr_array(data), row_names, col_names)
    
    # Test lower triangular without diagonal, skip zeros
    result = matrix.triangular(include_diagonal=False, skip_zeros=True)
    expected = [('C', 'A', 7)]
    assert result == expected, f"Sparse lower triangular (no diagonal, skip zeros) failed: {result}"
    
    # Test lower triangular with diagonal, skip zeros
    result = matrix.triangular(include_diagonal=True, skip_zeros=True)
    expected = [('A', 'A', 1), ('B', 'B', 5), ('C', 'A', 7), ('C', 'C', 9)]
    assert result == expected, f"Sparse lower triangular (with diagonal, skip zeros) failed: {result}"


# Test triangular method with skip_zeros=False
def test_triangular_include_zeros():
    """Test that triangular correctly includes zeros when skip_zeros=False."""
    data = np.array([[1, 0, 3], [0, 5, 0], [7, 0, 9]])
    row_names = ['A', 'B', 'C']
    col_names = ['A', 'B', 'C']

    dense_matrix = DataMatrix(data, row_names, col_names)
    sparse_matrix = DataMatrix(scipy.sparse.csr_array(data), row_names, col_names)
    
    # Test lower triangular with diagonal, including zeros
    dense_result = dense_matrix.triangular(include_diagonal=True, skip_zeros=False)
    sparse_result = sparse_matrix.triangular(include_diagonal=True, skip_zeros=False)
    
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
    result_name = matrix.triangular(include_diagonal=True, skip_zeros=True, index_style='name')
    assert result_name[0] == ('A', 'A', 1)
    assert result_name[1] == ('B', 'A', 4)
    
    # Test with index_style='index'
    result_index = matrix.triangular(include_diagonal=True, skip_zeros=True, index_style='index')
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
    for include_diagonal in [True, False]:
        for skip_zeros in [True, False]:
            dense_result = dense.triangular(include_diagonal=include_diagonal, skip_zeros=skip_zeros)
            sparse_result = sparse.triangular(include_diagonal=include_diagonal, skip_zeros=skip_zeros)

            assert dense_result == sparse_result, \
                f"Dense and sparse differ for diagonal={include_diagonal}, skip_zeros={skip_zeros}"


def test_triangular_iter_matches_triangular():
    """Test that triangular_iter preserves triangular ordering and values."""
    data = np.array([[1, 0, 3], [4, 5, 0], [7, 8, 9]])
    row_names = ['A', 'B', 'C']
    col_names = ['A', 'B', 'C']

    dense = DataMatrix(data, row_names, col_names)
    sparse = DataMatrix(scipy.sparse.csr_array(data), row_names, col_names)

    for matrix in (dense, sparse):
        expected = matrix.triangular(include_diagonal=False, skip_zeros=True, agg=max)
        observed = list(matrix.triangular_iter(include_diagonal=False, skip_zeros=True, agg=max))
        assert observed == expected


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
    result = matrix.triangular(include_diagonal=True, skip_zeros=True, index_style='index')
    
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


def test_estimate_dense_text_size_matches_written_file(tmp_path):
    data = np.array([[1.5, 0.0], [2.0, 3.25]], dtype=np.float64)
    row_names = ['row1', 'row2']
    col_names = ['col1', 'col2']
    out_file = tmp_path / "matrix.tsv"

    estimated_size = DataMatrix.estimate_dense_text_size(data, row_names, col_names)
    DataMatrix.write_dense_text(data, out_file, row_names, col_names)

    assert estimated_size == out_file.stat().st_size


def test_estimate_hdf_size_matches_written_file_and_halves_for_symmetric(tmp_path):
    rng = np.random.default_rng(0)
    upper = np.triu(rng.random((60, 60)), k=1)
    symmetric = upper + upper.T + np.eye(60)
    asymmetric = symmetric.copy()
    asymmetric[0, 1] += 1.0
    names = [f'N{i}' for i in range(60)]

    symmetric_size = DataMatrix.estimate_hdf_size(symmetric, names, names)
    asymmetric_size = DataMatrix.estimate_hdf_size(asymmetric, names, names)
    forced_full_size = DataMatrix.estimate_hdf_size(symmetric, names, names, symmetric=False)

    assert forced_full_size == asymmetric_size
    # values + int32 column indices are stored for about half of the cells
    data_bytes = 60 * 60 * (8 + 4)
    assert asymmetric_size - symmetric_size == pytest.approx(data_bytes / 2, rel=0.05)

    for data, estimate in ((symmetric, symmetric_size), (asymmetric, asymmetric_size)):
        out_file = tmp_path / "estimate.hdf5"
        DataMatrix.write_hdf(data, out_file, names, names)
        assert estimate >= out_file.stat().st_size


def test_estimate_hdf_size_does_not_modify_input():
    data = scipy.sparse.csr_array((np.array([1.0, 0.0, 2.0]), np.array([1, 0, 0]), np.array([0, 2, 3])), shape=(2, 2))
    DataMatrix.estimate_hdf_size(data, ['A', 'B'], ['A', 'B'])
    assert data.nnz == 3  # the explicit zero and unsorted indices are left alone
    assert not data.has_sorted_indices


def test_estimate_write_size_rejects_unknown_output_type():
    data = np.array([[1.0]], dtype=np.float64)

    with pytest.raises(ValueError, match="Unrecognized output type"):
        DataMatrix.estimate_write_size("bogus", data, ['A'], ['A'])


def test_from_file_sparse_hdf5_applies_lower_bound(tmp_path):
    data = scipy.sparse.csr_array(np.array([
        [0.0, 0.2, 0.8],
        [0.2, 0.0, 0.5],
        [0.8, 0.5, 0.0],
    ], dtype=float))
    names = ['A', 'B', 'C']
    matrix_path = tmp_path / 'matrix.hdf5'

    DataMatrix.write_hdf(data, matrix_path, names, names)

    matrix = DataMatrix.from_file(matrix_path, lower_bound=0.5)

    assert isinstance(matrix, DataMatrix)
    assert matrix.toarray().tolist() == [
        [0.0, 0.0, 0.8],
        [0.0, 0.0, 0.0],
        [0.8, 0.0, 0.0],
    ]


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
    result = matrix.triangular(include_diagonal=False, skip_zeros=True, agg=lambda a, b: a + b)
    # Lower [1,0] = 4, transpose [0,1] = 2, sum = 6
    # Lower [2,0] = 7, transpose [0,2] = 3, sum = 10
    # Lower [2,1] = 8, transpose [1,2] = 6, sum = 14
    expected = [('B', 'A', 6), ('C', 'A', 10), ('C', 'B', 14)]
    assert result == expected, f"Sum aggregation failed: {result}"
    
    # Test with max aggregation
    result = matrix.triangular(include_diagonal=False, skip_zeros=True, agg=max)
    expected = [('B', 'A', 4), ('C', 'A', 7), ('C', 'B', 8)]
    assert result == expected, f"Max aggregation failed: {result}"
    
    # Test with min aggregation
    result = matrix.triangular(include_diagonal=False, skip_zeros=True, agg=min)
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
    result = matrix.triangular(include_diagonal=False, skip_zeros=True, agg=lambda a, b: a + b)
    # Lower [2,0] = 7, transpose [0,2] = 3, sum = 10
    expected = [('C', 'A', 10)]
    assert result == expected, f"Sparse sum aggregation failed: {result}"
    
    # Test with max aggregation
    result = matrix.triangular(include_diagonal=False, skip_zeros=True, agg=max)
    expected = [('C', 'A', 7)]
    assert result == expected, f"Sparse max aggregation failed: {result}"


# Test triangular with agg parameter - averaging for symmetric matrix
def test_triangular_agg_average():
    """Test that triangular with average aggregation works correctly."""
    data = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    row_names = ['A', 'B', 'C']
    col_names = ['A', 'B', 'C']
    matrix = DataMatrix(data, row_names, col_names)
    
    # Test with average aggregation
    result = matrix.triangular(include_diagonal=False, skip_zeros=True,
                               agg=lambda a, b: (a + b) / 2)
    expected = [('B', 'A', 3.0), ('C', 'A', 5.0), ('C', 'B', 7.0)]
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
    dense_result = dense.triangular(include_diagonal=False, skip_zeros=True, 
                                    agg=lambda a, b: a + b)
    sparse_result = sparse.triangular(include_diagonal=False, skip_zeros=True, 
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
    result_no_agg = matrix.triangular(include_diagonal=False, skip_zeros=True)
    # [1,0]=0 skipped, [2,0]=0 skipped, [2,1]=0 skipped
    # No items
    assert len(result_no_agg) == 0, f"Without agg should skip all zeros: {result_no_agg}"
    
    # With agg: only skips when BOTH value and symmetric complement are zero
    result_with_agg = matrix.triangular(include_diagonal=False, skip_zeros=True, 
                                       agg=lambda a, b: a + b)
    # [1,0]=0, transpose [0,1]=2, both zero? No -> include, sum=2
    # [2,0]=0, transpose [0,2]=0, both zero? Yes -> SKIP
    # [2,1]=0, transpose [1,2]=6, both zero? No -> include, sum=6
    expected = [('B', 'A', 2), ('C', 'B', 6)]
    assert result_with_agg == expected, f"With agg should include when one is non-zero: {result_with_agg}"


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

    def test_iter_mst_edges_matches_property(self):
        """Test that iter_mst_edges yields the same MST edge sequence as mst_edges."""
        data = np.array([
            [0, 10, 5],
            [10, 0, 8],
            [5, 8, 0]
        ])
        row_names = ['A', 'B', 'C']
        matrix = DataMatrix(data, row_names, row_names)

        tree = MaxTree(matrix)

        assert list(tree.iter_mst_edges()) == tree.mst_edges
    
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
        
        # One row per distinct MST weight (10, 8) plus the `--lb 0` full-graph row.
        # Each row holds the edges a `--lb <threshold>` cut keeps (scores strictly above).
        assert result.shape == (3, 2)
        assert result.tolist() == [[0.0, 10.0], [1.0, 8.0], [3.0, 0.0]]
        
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
        
        # One row per distinct MST weight (10, 8) plus the `--lb 0` full-graph row; each
        # row is the cluster count of the `--lb <threshold>` cut.
        assert result.shape == (3, 2)
        assert result.tolist() == [[10.0, 3.0], [8.0, 2.0], [0.0, 1.0]]
        
        # Highest threshold keeps nothing, lowest keeps everything
        assert result[0, 1] == 3
        assert result[-1, 1] == 1
        
        # Cluster counts should be decreasing
        cluster_counts = result[:, 1]
        assert all(cluster_counts[i] >= cluster_counts[i+1] for i in range(len(cluster_counts)-1))
        
        # Thresholds should be decreasing
        thresholds = result[:, 0]
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
        
        # Rows are the `--lb <threshold>` cuts, keyed by the edge count each one keeps.
        assert result[0].tolist() == [0.0, 4.0]   # nothing above the top weight
        assert result[-1].tolist() == [6.0, 1.0]  # `--lb 0`: every edge, one cluster
        
        # Cluster count is always n_nodes minus the MST edges strictly above the threshold
        for row_idx, threshold in enumerate(tree.thresholds):
            mst_above = sum(1 for _, _, weight in tree.mst_edges if weight > threshold)
            assert result[row_idx, 1] == tree.n_nodes - mst_above
            assert result[row_idx, 0] == tree.edges_by_threshold[row_idx, 0]
    
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

def test_sparse_write_narrows_index_dtype(tmp_path):
    """int64 index arrays are written as int32 when the matrix is small enough to address."""
    import h5py
    from domainator.data_matrix import _sparse_index_dtype

    data = np.array([
        [0.0, 5.0, 1.0],
        [5.0, 0.0, 2.0],
        [1.0, 2.0, 0.0],
    ])
    labels = ['a', 'b', 'c']

    wide = scipy.sparse.csr_array(data)
    wide.indices = wide.indices.astype(np.int64)
    wide.indptr = wide.indptr.astype(np.int64)

    out_file = tmp_path / "narrowed.sparse.hdf5"
    DataMatrix(wide, labels, labels, data_type="score").write(out_file, "hdf")

    with h5py.File(out_file, "r") as handle:
        assert handle["SPARSE_CSR_INDICES"].dtype == np.int32
        assert handle["SPARSE_CSR_INDPTR"].dtype == np.int32

    # values are untouched and the matrix round trips
    read_back = DataMatrix.from_file(out_file)
    np.testing.assert_array_equal(read_back.toarray(), data)
    assert read_back.rows == labels

    # the guardrail estimate must assume the same (narrowed) width as the write
    estimate = DataMatrix.estimate_hdf_size(wide, labels, labels, data_type="score")
    assert estimate >= out_file.stat().st_size

    class _FakeMatrix:
        def __init__(self, shape, nnz):
            self.shape = shape
            self.nnz = nnz

    int32_max = np.iinfo(np.int32).max
    assert _sparse_index_dtype(_FakeMatrix((10, 10), 50)) == np.int32
    # too many non-zeros for int32 indptr, or too many rows/columns for int32 indices
    assert _sparse_index_dtype(_FakeMatrix((10, 10), int32_max)) == np.int32
    assert _sparse_index_dtype(_FakeMatrix((10, 10), int32_max + 1)) == np.int64
    assert _sparse_index_dtype(_FakeMatrix((int32_max + 1, 10), 50)) == np.int64


def _symmetric_test_matrix():
    return np.array([
        [1.0, 2.0, 0.0, 4.0],
        [2.0, 0.0, 3.0, 0.0],
        [0.0, 3.0, 5.0, 0.5],
        [4.0, 0.0, 0.5, 0.0],
    ])


def test_symmetric_matrix_is_stored_as_upper_triangle(tmp_path):
    import h5py

    data = _symmetric_test_matrix()
    names = ['A', 'B', 'C', 'D']
    out_file = tmp_path / "symmetric.hdf5"
    DataMatrix.write_hdf(data, out_file, names, names, data_type="score")

    with h5py.File(out_file, "r") as handle:
        assert handle.attrs["ARRAY_TYPE"] == "SPARSE_CSR_UPPER_TRIANGLE"
        assert handle.attrs["MATRIX_FILE_VERSION"] == "2.0"
        assert len(handle["SPARSE_VALUES"]) == np.count_nonzero(np.triu(data))

    matrix = DataMatrix.from_file(out_file)
    np.testing.assert_array_equal(matrix.toarray(), data)
    assert matrix.exactly_symmetric
    assert matrix.symmetric


@pytest.mark.parametrize("case", ["asymmetric_values", "asymmetric_labels", "forced_full"])
def test_matrix_is_stored_full_unless_exactly_symmetric(tmp_path, case):
    import h5py

    data = _symmetric_test_matrix()
    rows = ['A', 'B', 'C', 'D']
    cols = rows
    kwargs = {}
    if case == "asymmetric_values":
        data[0, 1] = 2.0000001  # within np.allclose, but not exact
    elif case == "asymmetric_labels":
        cols = ['W', 'X', 'Y', 'Z']
    else:
        kwargs["symmetric"] = False
    out_file = tmp_path / f"{case}.hdf5"
    DataMatrix.write_hdf(data, out_file, rows, cols, **kwargs)

    with h5py.File(out_file, "r") as handle:
        assert handle.attrs["ARRAY_TYPE"] == "SPARSE_CSR"
    np.testing.assert_array_equal(DataMatrix.from_file(out_file).toarray(), data)


def test_matrix_write_reuses_known_symmetry(tmp_path):
    import h5py

    names = ['A', 'B', 'C', 'D']
    matrix = DataMatrix(_symmetric_test_matrix(), names, names)
    out_file = tmp_path / "written.hdf5"
    matrix.write(out_file, "hdf")
    with h5py.File(out_file, "r") as handle:
        assert handle.attrs["ARRAY_TYPE"] == "SPARSE_CSR_UPPER_TRIANGLE"
    np.testing.assert_array_equal(DataMatrix.from_file(out_file).toarray(), matrix.toarray())


def test_lower_bound_applies_to_triangle_files_and_text(tmp_path):
    data = _symmetric_test_matrix()
    names = ['A', 'B', 'C', 'D']
    hdf_file = tmp_path / "symmetric.hdf5"
    text_file = tmp_path / "symmetric.tsv"
    DataMatrix.write_hdf(data, hdf_file, names, names)
    DataMatrix.write_dense_text(data, text_file, names, names)

    expected = np.where(data <= 1.0, 0.0, data)
    for path in (hdf_file, text_file):
        matrix = DataMatrix.from_file(path, lower_bound=1.0)
        np.testing.assert_array_equal(matrix.toarray(), expected)
        assert matrix.exactly_symmetric


@pytest.mark.parametrize("filename", ["test_matrix.dense.hdf5", "test_matrix.sparse.hdf5", "FeSOD_dist.dense.hdf5", "FeSOD_dist.sparse.hdf5"])
def test_version_1_0_files_still_load(shared_datadir, filename):
    import h5py

    path = shared_datadir / filename
    with h5py.File(path, "r") as handle:
        assert handle.attrs.get("MATRIX_FILE_VERSION", handle.attrs.get("DOMAINATOR_MATRIX_FILE_VERSION")) == "1.0"
    matrix = DataMatrix.from_file(path)
    assert isinstance(matrix.data, scipy.sparse.csr_array)


def test_legacy_dense_and_sparse_files_hold_the_same_matrix(shared_datadir):
    dense = DataMatrix.from_file(shared_datadir / "FeSOD_dist.dense.hdf5")
    sparse = DataMatrix.from_file(shared_datadir / "FeSOD_dist.sparse.hdf5")
    text = DataMatrix.from_file(shared_datadir / "FeSOD_dist.tsv")
    assert dense.rows == sparse.rows == text.rows
    np.testing.assert_array_equal(dense.toarray(), sparse.toarray())
    np.testing.assert_allclose(dense.toarray(), text.toarray())


@pytest.mark.parametrize("version", ["1.0", "1.1", "2.0", "2.3"])
def test_major_versions_1_and_2_are_readable(tmp_path, version):
    """1.1 is the pre-release version of SPARSE_CSR_UPPER_TRIANGLE; 2.3 stands in for a future additive change."""
    import h5py

    data = _symmetric_test_matrix()
    names = ['A', 'B', 'C', 'D']
    out_file = tmp_path / "versioned.hdf5"
    DataMatrix.write_hdf(data, out_file, names, names)
    with h5py.File(out_file, "a") as handle:
        handle.attrs["MATRIX_FILE_VERSION"] = version
    np.testing.assert_array_equal(DataMatrix.from_file(out_file).toarray(), data)


@pytest.mark.parametrize("attribute,value,message", [
    ("MATRIX_FILE_VERSION", "3.0", "Unsupported matrix file version"),
    ("ARRAY_TYPE", "SPARSE_COO", "Unrecognized matrix ARRAY_TYPE"),
])
def test_unreadable_files_are_rejected(tmp_path, attribute, value, message):
    import h5py

    out_file = tmp_path / "future.hdf5"
    DataMatrix.write_hdf(np.eye(2), out_file, ['A', 'B'], ['A', 'B'])
    with h5py.File(out_file, "a") as handle:
        handle.attrs[attribute] = value
    with pytest.raises(ValueError, match=message):
        DataMatrix.from_file(out_file)


def test_symmetrize_matrix_methods():
    from domainator.data_matrix import symmetrize_matrix

    data = np.array([[1.0, 8.0, 0.0], [2.0, 0.0, 3.0], [0.0, 0.0, 4.0]])
    for array in (data, scipy.sparse.csr_array(data)):
        maxed = symmetrize_matrix(array, "max")
        maxed = maxed.toarray() if scipy.sparse.issparse(maxed) else maxed
        np.testing.assert_array_equal(maxed, [[1.0, 8.0, 0.0], [8.0, 0.0, 3.0], [0.0, 3.0, 4.0]])
        mined = symmetrize_matrix(array, "min")
        mined = mined.toarray() if scipy.sparse.issparse(mined) else mined
        np.testing.assert_array_equal(mined, [[1.0, 2.0, 0.0], [2.0, 0.0, 0.0], [0.0, 0.0, 4.0]])
        meaned = symmetrize_matrix(array, "mean")
        meaned = meaned.toarray() if scipy.sparse.issparse(meaned) else meaned
        np.testing.assert_array_equal(meaned, [[1.0, 5.0, 0.0], [5.0, 0.0, 1.5], [0.0, 1.5, 4.0]])
        assert symmetrize_matrix(array, "none") is array
    with pytest.raises(ValueError, match="square"):
        symmetrize_matrix(np.zeros((2, 3)), "max")


def test_negative_scores_rank_the_same_from_numpy_and_csr_input():
    # Before the dense class was removed, dense kNN rankings kept only scores > 0 while
    # sparse rankings kept any non-zero score.
    data = np.array([[0.0, -1.0, 0.0], [-1.0, 0.0, -2.0], [0.0, -2.0, 0.0]])
    names = ['A', 'B', 'C']
    from_numpy = build_symmetric_neighbor_rankings(DataMatrix(data, names, names))
    from_csr = build_symmetric_neighbor_rankings(DataMatrix(scipy.sparse.csr_array(data), names, names))
    np.testing.assert_array_equal(from_numpy.offsets, from_csr.offsets)
    np.testing.assert_array_equal(from_numpy.target, from_csr.target)
    np.testing.assert_array_equal(from_numpy.score, from_csr.score)


def test_require_similarity_rejects_distance_matrices():
    matrix = DataMatrix(np.eye(2), ['A', 'B'], ['A', 'B'], data_type="score_dist")
    with pytest.raises(ValueError, match="'score_dist' distance matrix"):
        matrix.require_similarity("build_ssn")
    DataMatrix(np.eye(2), ['A', 'B'], ['A', 'B'], data_type="score").require_similarity("build_ssn")


def test_duplicate_csr_entries_are_summed_before_zeros_are_dropped():
    # (0, 1) appears twice and sums to 0, so it must not be stored.
    data = scipy.sparse.csr_array((np.array([1.0, -1.0, 2.0]), np.array([1, 1, 0]), np.array([0, 2, 3])), shape=(2, 2))
    matrix = DataMatrix(data, ['A', 'B'], ['A', 'B'])
    assert matrix.nnz == 1
    np.testing.assert_array_equal(matrix.toarray(), [[0.0, 0.0], [2.0, 0.0]])
