"""Shared --mst / --knn graph sparsification for the matrix-producing and SSN tools.

Every tool that can prune a similarity graph exposes the same two independent options:

* ``--mst`` keeps the maximum spanning tree (forest), which preserves the connected
  components of the input.
* ``--knn K`` keeps the OR-symmetric K-nearest-neighbor edges, which need not preserve them.

Passing both keeps the union of the two.
"""

from typing import Dict, Iterator, Optional, Tuple

import numpy as np
import scipy.sparse

from domainator.ssn_edges import MaxTree, mst_edge_index_dict, mst_knn_edge_index_dict, symmetric_knn_edge_index_dict


MAX_HITS_PER_QUERY_SPARSIFY_NOTE = (
    "Applied before --mst and --knn, so a cap can split clusters that --mst would otherwise keep connected."
)


def knn_arg(value) -> int:
    """argparse type for --knn: an integer >= 1."""
    value = int(value)
    if value < 1:
        # k=0 selects no neighbors; --mst on its own is the way to ask for just the tree.
        raise ValueError("--knn must be an integer >= 1")
    return value


def add_sparsify_arguments(parser, *, requirements: str = "") -> None:
    """Add the shared --mst and --knn options to a CLI parser.

    ``requirements`` is appended to both help strings, for tool-specific constraints such
    as which --mode values can be pruned while streaming.
    """
    suffix = f" {requirements}" if requirements else ""
    parser.add_argument("--mst", action="store_true", default=False,
                        help="Keep only the edges of the maximum spanning tree (forest). The clusters "
                             "(connected components) are the same as the full graph's, with the "
                             "connections inside each cluster pruned to the minimum needed to keep it "
                             "connected. Combine with --knn to also keep each node's nearest neighbors." + suffix)
    parser.add_argument("--knn", type=knn_arg, default=None,
                        help="Keep only OR-symmetric K-nearest-neighbor edges (K >= 1): an edge is kept if "
                             "it is among the K best for either of its nodes. On its own this can split a "
                             "cluster of the full graph into several; add --mst to keep the clusters "
                             "connected." + suffix)


def sparsify_requested(mst: bool, knn: Optional[int]) -> bool:
    return bool(mst) or knn is not None


def sparsify_option_label(mst: bool, knn: Optional[int]) -> str:
    """The options as the user typed them, e.g. ``--mst --knn 5``, for error messages."""
    parts = []
    if mst:
        parts.append("--mst")
    if knn is not None:
        parts.append(f"--knn {knn}")
    return " ".join(parts)


def check_sparsify_data_type(data_type: Optional[str], mst: bool, knn: Optional[int]) -> None:
    """Reject --mst/--knn on distance-valued data.

    Both select the *highest* values, which are the weakest relationships in a distance
    matrix, so the result would be the opposite of what was asked for.
    """
    if not sparsify_requested(mst, knn) or not data_type:
        return
    if data_type.endswith("_dist"):
        score_type = data_type[:-len("_dist")]
        raise ValueError(
            f"{sparsify_option_label(mst, knn)} keeps the highest-valued edges, so it cannot be applied "
            f"to a '{data_type}' (distance) matrix. Sparsify the '{score_type}' matrix instead."
        )


def sparsify_edge_index_dict(matrix, *, mst: bool, knn: Optional[int], lower_bound: float = 0,
                             tree: Optional[MaxTree] = None) -> Dict[Tuple[int, int], float]:
    """Return the edges --mst and/or --knn keep, keyed by undirected index pair."""
    if knn is not None and knn < 1:
        raise ValueError("knn must be >= 1")
    if mst and knn is not None:
        return mst_knn_edge_index_dict(matrix, knn, lower_bound=lower_bound, tree=tree)
    if mst:
        return mst_edge_index_dict(tree if tree is not None else MaxTree(matrix), lower_bound=lower_bound)
    if knn is not None:
        return symmetric_knn_edge_index_dict(matrix, knn, lower_bound=lower_bound)
    raise ValueError("Select mst, knn, or both; otherwise no edges would be kept.")


def sparsify_matrix(matrix, *, mst: bool, knn: Optional[int], lower_bound: float = 0):
    """Keep only the maximum spanning tree and/or OR-symmetric kNN edges of a square matrix.

    Kept edges are written as a strict max over the diagonal: both cells of a pair get
    ``max(M[i, j], M[j, i])``, the same score the OR-symmetric kNN ranking and the maximum
    spanning tree select on. The output is therefore symmetric, matching
    :class:`~domainator.ssn_edges.StreamingMstKnnAccumulator`. The diagonal is copied
    unchanged. Returns the new data as a CSR array.
    """
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"{sparsify_option_label(mst, knn)} requires a square matrix.")

    edge_dict = sparsify_edge_index_dict(matrix, mst=mst, knn=knn, lower_bound=lower_bound)
    n = matrix.shape[0]
    data = matrix.data

    # One (low, high) row per undirected pair, so no cell is written twice (COO input sums duplicates).
    edges = np.sort(np.asarray(list(edge_dict.keys()), dtype=np.int64).reshape(-1, 2), axis=1)
    edges = np.unique(edges, axis=0)
    sources, targets = edges[:, 0], edges[:, 1]
    edge_values = np.asarray(data.maximum(data.transpose()).tocsr()[sources, targets]).ravel()

    diagonal = data.diagonal()
    diagonal_indices = np.flatnonzero(diagonal)

    rows = np.concatenate([diagonal_indices, sources, targets])
    cols = np.concatenate([diagonal_indices, targets, sources])
    values = np.concatenate([diagonal[diagonal_indices], edge_values, edge_values])
    out = scipy.sparse.csr_array((values, (rows, cols)), shape=(n, n))
    out.eliminate_zeros()
    out.sort_indices()
    return out


class PerQueryTopK:
    """Keep each row's ``cap`` highest-valued columns from a stream of (row, col, value) hits.

    The streaming equivalent of ``sparse_scores.keep_top_k_per_row`` for hits that arrive
    in no particular order: duplicate (row, col) pairs keep their maximum value, and ties
    are broken toward the lower column index. Memory is bounded at O(rows * cap).
    """

    def __init__(self, n_rows: int, cap: int):
        if cap < 0:
            raise ValueError("cap must be >= 0")
        self.cap = cap
        self._rows = [dict() for _ in range(n_rows)]

    def add_edge(self, row: int, col: int, value: float) -> None:
        if self.cap == 0:
            return
        hits = self._rows[row]
        if value > hits.get(col, -np.inf):
            hits[col] = value
        # Trimming is exact: anything dropped is outranked by `cap` values that can only grow.
        if len(hits) > 2 * self.cap:
            self._rows[row] = self._top(hits)

    def _top(self, hits: Dict[int, float]) -> Dict[int, float]:
        return dict(sorted(hits.items(), key=lambda item: (-item[1], item[0]))[:self.cap])

    def iter_edges(self) -> Iterator[Tuple[int, int, float]]:
        for row, hits in enumerate(self._rows):
            for col, value in self._top(hits).items():
                yield row, col, value
