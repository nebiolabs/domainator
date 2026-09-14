import numpy as np
import pytest
import scipy.sparse

from domainator.data_matrix import SparseDataMatrix
from domainator.ssn_edges import StreamingMstKnnAccumulator
from domainator.transform_matrix import apply_mst_knn_sparsification


def _feed_dense(accumulator, data):
    """Feed every off-diagonal entry of a dense matrix into the accumulator."""
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if i != j and data[i, j] != 0:
                accumulator.add_edge(i, j, float(data[i, j]))
    return accumulator


def _as_matrix(data):
    labels = [f"s{i}" for i in range(data.shape[0])]
    return SparseDataMatrix(scipy.sparse.csr_array(data), labels, labels, data_type="score")


def test_streaming_knn_retains_negative_scores():
    accumulator = StreamingMstKnnAccumulator(
        n_nodes=3,
        k=1,
        lower_bound=-1.0,
        include_mst=False,
    )
    accumulator.add_edge(0, 1, -0.2)
    accumulator.add_edge(0, 2, -0.8)
    accumulator.add_edge(1, 2, -0.4)

    graph = accumulator.to_csr()

    np.testing.assert_array_equal(
        graph.toarray(),
        np.array(
            [
                [0.0, -0.2, 0.0],
                [-0.2, 0.0, -0.4],
                [0.0, -0.4, 0.0],
            ]
        ),
    )


def test_streaming_mst_knn_retains_negative_scores():
    # include_mst=True, so the forest seeds edge_dict before the kNN pass runs. The two
    # kNN-only edges here, (0,2) and (1,3), are the ones a "0.0 means absent" default in
    # finalize() would silently discard.
    data = np.array([
        [0.0, -0.2, -0.8, -1.5],
        [-0.2, 0.0, -0.4, -1.1],
        [-0.8, -0.4, 0.0, -0.6],
        [-1.5, -1.1, -0.6, 0.0],
    ])
    k = 2
    lb = -2.0

    accumulator = _feed_dense(
        StreamingMstKnnAccumulator(len(data), k, lower_bound=lb, include_mst=True), data
    )
    result = accumulator.to_csr().toarray()

    # MST keeps (0,1), (1,2), (2,3); k=2 adds (0,2) and (1,3); (0,3) stays out.
    np.testing.assert_array_equal(
        result,
        np.array([
            [0.0, -0.2, -0.8, 0.0],
            [-0.2, 0.0, -0.4, -1.1],
            [-0.8, -0.4, 0.0, -0.6],
            [0.0, -1.1, -0.6, 0.0],
        ]),
    )
    np.testing.assert_array_equal(
        result, apply_mst_knn_sparsification(_as_matrix(data), k, lower_bound=lb).toarray()
    )


def test_streaming_negative_scores_below_lower_bound_are_dropped():
    # A negative lower bound must still prune: only (0,1) clears -0.5.
    data = np.array([
        [0.0, -0.2, -0.8],
        [-0.2, 0.0, -0.9],
        [-0.8, -0.9, 0.0],
    ])

    accumulator = _feed_dense(
        StreamingMstKnnAccumulator(3, 2, lower_bound=-0.5, include_mst=True), data
    )

    np.testing.assert_array_equal(
        accumulator.to_csr().toarray(),
        np.array([
            [0.0, -0.2, 0.0],
            [-0.2, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ]),
    )


@pytest.mark.parametrize("include_mst", [True, False])
@pytest.mark.parametrize("k", [1, 2, 3])
@pytest.mark.parametrize("lb", [-10.0, -0.5])
def test_streaming_matches_batch_on_negative_scores(lb, k, include_mst):
    # Randomized parity against the batch path. Distinct weights everywhere, so the
    # maximum spanning forest is unique and there is no tie ambiguity to tolerate.
    rng = np.random.default_rng(20240917)
    for _ in range(10):
        n = int(rng.integers(4, 9))
        upper = np.round(rng.normal(0.0, 1.0, size=(n, n)), 3)
        data = np.triu(upper, 1)
        data = data + data.T

        expected = apply_mst_knn_sparsification(
            _as_matrix(data), k, lower_bound=lb, include_mst=include_mst
        ).toarray()
        result = _feed_dense(
            StreamingMstKnnAccumulator(n, k, lower_bound=lb, include_mst=include_mst), data
        ).to_csr().toarray()

        np.testing.assert_allclose(result, expected)
