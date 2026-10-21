import numpy as np
import pytest
import scipy.sparse

from domainator.sparse_scores import keep_top_k_per_row
from domainator.data_matrix import DataMatrix
from domainator.sparsify import (
    PerQueryTopK,
    check_sparsify_data_type,
    knn_arg,
    sparsify_matrix,
    sparsify_option_label,
)


def test_knn_arg_requires_at_least_one():
    assert knn_arg("3") == 3
    for bad in ("0", "-1"):
        with pytest.raises(ValueError):
            knn_arg(bad)


def test_sparsify_option_label():
    assert sparsify_option_label(True, None) == "--mst"
    assert sparsify_option_label(False, 4) == "--knn 4"
    assert sparsify_option_label(True, 4) == "--mst --knn 4"


@pytest.mark.parametrize("data_type", ["score_dist", "efi_score_dist", "tmscore_dist"])
def test_check_sparsify_data_type_rejects_distances(data_type):
    with pytest.raises(ValueError, match="distance"):
        check_sparsify_data_type(data_type, True, None)
    # nothing to check when no pruning was asked for
    check_sparsify_data_type(data_type, False, None)


@pytest.mark.parametrize("data_type", [None, "", "score", "bool", "efi_score", "tmscore"])
def test_check_sparsify_data_type_accepts_scores(data_type):
    check_sparsify_data_type(data_type, True, 2)


def test_sparsify_matrix_requires_square_and_a_selection():
    labels = ["a", "b", "c"]
    with pytest.raises(ValueError, match="square"):
        sparsify_matrix(DataMatrix(np.ones((3, 2)), labels, labels[:2], data_type="score"), mst=True, knn=None)
    square = DataMatrix(np.ones((3, 3)), labels, labels, data_type="score")
    with pytest.raises(ValueError):
        sparsify_matrix(square, mst=False, knn=None)
    with pytest.raises(ValueError):
        sparsify_matrix(square, mst=False, knn=0)


def test_sparsify_matrix_mst_union_knn():
    # Two triangles joined by a weak bridge that is nobody's nearest neighbor.
    labels = list("abcdef")
    weights = {(0, 1): 10.0, (0, 2): 9.0, (1, 2): 8.0, (3, 4): 7.0, (3, 5): 6.0, (4, 5): 5.0, (2, 3): 1.0}
    data = np.zeros((6, 6))
    for (i, j), value in weights.items():
        data[i, j] = data[j, i] = value
    matrix = DataMatrix(scipy.sparse.csr_array(data), labels, labels, data_type="score")

    def pairs(out):
        return {(i, j) for i, j in zip(*np.nonzero(out.toarray())) if i < j}

    mst = pairs(sparsify_matrix(matrix, mst=True, knn=None))
    knn = pairs(sparsify_matrix(matrix, mst=False, knn=1))
    both = pairs(sparsify_matrix(matrix, mst=True, knn=1))
    assert (2, 3) in mst and (2, 3) not in knn
    assert len(mst) == len(labels) - 1
    assert both == mst | knn


@pytest.mark.parametrize("cap", [0, 1, 2, 5])
def test_per_query_top_k_matches_batch(cap):
    rng = np.random.default_rng(7)
    n = 12
    # Rounded values give ties, which must break toward the lower column like the batch path.
    data = np.round(rng.random((n, n)), 1) * (rng.random((n, n)) < 0.6)

    top = PerQueryTopK(n, cap)
    hits = [(i, j, data[i, j]) for i in range(n) for j in range(n) if data[i, j] != 0]
    rng.shuffle(hits)
    for i, j, value in hits:
        # duplicates arrive too; the pair keeps its maximum
        top.add_edge(int(i), int(j), float(value) / 2)
        top.add_edge(int(i), int(j), float(value))
    streamed = np.zeros((n, n))
    for i, j, value in top.iter_edges():
        streamed[i, j] = value

    expected = keep_top_k_per_row(scipy.sparse.csr_array(data), cap).toarray()
    np.testing.assert_array_equal(streamed, expected)
