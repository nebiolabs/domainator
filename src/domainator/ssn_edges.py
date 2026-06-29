"""Graph machinery built on top of DataMatrix: neighbor rankings, MST/kNN edges.

This module contains the maximum-spanning-tree (:class:`MaxTree`) and the neighbor
ranking / edge-selection functions used to build sequence similarity networks. It is
split out of ``data_matrix`` to keep that module focused on the matrix abstraction.

To avoid an import cycle, the concrete :class:`~domainator.data_matrix.DataMatrix` class
is imported lazily inside :meth:`MaxTree.__init__`; everything else references it only via
string type annotations (enabled by ``from __future__ import annotations``).
"""
from __future__ import annotations

import heapq
import warnings
from typing import Dict, Iterator, List, Optional, Tuple, Union

import numpy as np
import scipy.sparse

from .edge_types import (
    CompactNeighborRankings,
    NeighborRankings,
    SortedUndirectedEdges,
    _empty_compact_neighbor_rankings,
    _iter_neighbor_ranking_rows,
    _score_passes_lower_bound,
)

try:
    import numba as nb
except ImportError:
    nb = None


def _build_dense_neighbor_rankings(data: np.ndarray, max_k: Optional[int] = None) -> CompactNeighborRankings:
    offsets = [0]
    target_parts = []
    score_parts = []
    for row_idx in range(data.shape[0]):
        row_scores = np.maximum(data[row_idx, :], data[:, row_idx])
        candidate_indices = np.flatnonzero(row_scores > 0)
        candidate_indices = candidate_indices[candidate_indices != row_idx]

        if candidate_indices.size == 0:
            offsets.append(offsets[-1])
            continue

        candidate_scores = row_scores[candidate_indices]
        order = np.lexsort((candidate_indices, -candidate_scores))
        if max_k is not None:
            order = order[:max_k]
        row_target = candidate_indices[order].astype(np.int32, copy=False)
        row_score = candidate_scores[order].astype(float, copy=False)
        target_parts.append(row_target)
        score_parts.append(row_score)
        offsets.append(offsets[-1] + row_target.size)

    if offsets[-1] == 0:
        return _empty_compact_neighbor_rankings(data.shape[0])

    return CompactNeighborRankings(
        offsets=np.asarray(offsets, dtype=np.int64),
        target=np.concatenate(target_parts),
        score=np.concatenate(score_parts),
    )


def _build_top_k_neighbor_rankings_heap(edges: SortedUndirectedEdges, max_k: int) -> CompactNeighborRankings:
    top_k_by_row = [[] for _ in range(edges.n_nodes)]

    def add_candidate(row_idx: int, target_idx: int, score: float):
        heap = top_k_by_row[row_idx]
        candidate = (score, -target_idx, target_idx)
        if len(heap) < max_k:
            heapq.heappush(heap, candidate)
            return
        if candidate > heap[0]:
            heapq.heapreplace(heap, candidate)

    for source_idx, target_idx, score in zip(edges.source, edges.target, edges.score):
        add_candidate(int(source_idx), int(target_idx), float(score))
        add_candidate(int(target_idx), int(source_idx), float(score))

    offsets = np.zeros(edges.n_nodes + 1, dtype=np.int64)
    total_kept = 0
    for row_idx, heap in enumerate(top_k_by_row):
        total_kept += len(heap)
        offsets[row_idx + 1] = total_kept

    if total_kept == 0:
        return _empty_compact_neighbor_rankings(edges.n_nodes)

    targets = np.empty(total_kept, dtype=np.int32)
    scores = np.empty(total_kept, dtype=float)

    out_idx = 0
    for heap in top_k_by_row:
        if len(heap) == 0:
            continue
        heap.sort(key=lambda item: (-item[0], item[2]))
        for score, _, target_idx in heap:
            targets[out_idx] = target_idx
            scores[out_idx] = score
            out_idx += 1

    return CompactNeighborRankings(
        offsets=offsets,
        target=targets,
        score=scores,
    )


if nb is not None:
    @nb.njit(cache=True)
    def _is_better_top_k_candidate(candidate_score: float, candidate_target: int,
                                   current_score: float, current_target: int) -> bool:
        if candidate_score > current_score:
            return True
        if candidate_score < current_score:
            return False
        return candidate_target < current_target


    @nb.njit(cache=True)
    def _insert_top_k_candidate(top_targets: np.ndarray, top_scores: np.ndarray, row_counts: np.ndarray,
                                row_idx: int, candidate_target: int, candidate_score: float, max_k: int) -> None:
        count = int(row_counts[row_idx])
        if count == max_k and not _is_better_top_k_candidate(
            candidate_score,
            candidate_target,
            top_scores[row_idx, max_k - 1],
            top_targets[row_idx, max_k - 1],
        ):
            return

        insert_pos = count
        if insert_pos >= max_k:
            insert_pos = max_k - 1

        while insert_pos > 0 and _is_better_top_k_candidate(
            candidate_score,
            candidate_target,
            top_scores[row_idx, insert_pos - 1],
            top_targets[row_idx, insert_pos - 1],
        ):
            if insert_pos < max_k:
                top_scores[row_idx, insert_pos] = top_scores[row_idx, insert_pos - 1]
                top_targets[row_idx, insert_pos] = top_targets[row_idx, insert_pos - 1]
            insert_pos -= 1

        top_scores[row_idx, insert_pos] = candidate_score
        top_targets[row_idx, insert_pos] = candidate_target

        if count < max_k:
            row_counts[row_idx] = count + 1


    @nb.njit(cache=True)
    def _build_top_k_neighbor_arrays(source: np.ndarray, target: np.ndarray, score: np.ndarray,
                                     n_nodes: int, max_k: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        top_targets = np.full((n_nodes, max_k), -1, dtype=np.int32)
        top_scores = np.empty((n_nodes, max_k), dtype=np.float64)
        row_counts = np.zeros(n_nodes, dtype=np.int32)

        for edge_idx in range(source.shape[0]):
            source_idx = int(source[edge_idx])
            target_idx = int(target[edge_idx])
            edge_score = float(score[edge_idx])
            _insert_top_k_candidate(top_targets, top_scores, row_counts, source_idx, target_idx, edge_score, max_k)
            _insert_top_k_candidate(top_targets, top_scores, row_counts, target_idx, source_idx, edge_score, max_k)

        return top_targets, top_scores, row_counts


    def _build_top_k_neighbor_rankings_numba(edges: SortedUndirectedEdges, max_k: int) -> CompactNeighborRankings:
        top_targets, top_scores, row_counts = _build_top_k_neighbor_arrays(
            edges.source,
            edges.target,
            edges.score,
            edges.n_nodes,
            max_k,
        )

        offsets = np.zeros(edges.n_nodes + 1, dtype=np.int64)
        offsets[1:] = np.cumsum(row_counts, dtype=np.int64)
        total_kept = int(offsets[-1])
        if total_kept == 0:
            return _empty_compact_neighbor_rankings(edges.n_nodes)

        targets = np.empty(total_kept, dtype=np.int32)
        scores = np.empty(total_kept, dtype=float)

        out_idx = 0
        for row_idx in range(edges.n_nodes):
            keep_count = int(row_counts[row_idx])
            if keep_count == 0:
                continue
            next_idx = out_idx + keep_count
            targets[out_idx:next_idx] = top_targets[row_idx, :keep_count]
            scores[out_idx:next_idx] = top_scores[row_idx, :keep_count]
            out_idx = next_idx

        return CompactNeighborRankings(
            offsets=offsets,
            target=targets,
            score=scores,
        )
else:
    def _build_top_k_neighbor_rankings_numba(edges: SortedUndirectedEdges, max_k: int) -> CompactNeighborRankings:
        return _build_top_k_neighbor_rankings_heap(edges, max_k)


def _build_neighbor_rankings_from_sorted_edges(edges: SortedUndirectedEdges, max_k: Optional[int] = None) -> CompactNeighborRankings:
    if len(edges) == 0:
        return _empty_compact_neighbor_rankings(edges.n_nodes)

    if max_k is not None:
        return _build_top_k_neighbor_rankings_numba(edges, max_k)

    row_counts = np.bincount(edges.source, minlength=edges.n_nodes) + np.bincount(edges.target, minlength=edges.n_nodes)
    total_directed = int(row_counts.sum())

    row_indices = np.empty(total_directed, dtype=np.int32)
    targets = np.empty(total_directed, dtype=np.int32)
    scores = np.empty(total_directed, dtype=float)
    edge_count = len(edges)

    row_indices[:edge_count] = edges.source
    targets[:edge_count] = edges.target
    scores[:edge_count] = edges.score

    row_indices[edge_count:] = edges.target
    targets[edge_count:] = edges.source
    scores[edge_count:] = edges.score

    order = np.lexsort((targets, -scores, row_indices))
    row_indices = row_indices[order]
    targets = targets[order]
    scores = scores[order]

    row_offsets = np.zeros(edges.n_nodes + 1, dtype=np.int64)
    row_offsets[1:] = np.cumsum(np.bincount(row_indices, minlength=edges.n_nodes))

    if max_k is None:
        return CompactNeighborRankings(
            offsets=row_offsets,
            target=targets,
            score=scores,
        )

    keep_counts = np.minimum(np.diff(row_offsets), max_k)
    kept_offsets = np.zeros(edges.n_nodes + 1, dtype=np.int64)
    kept_offsets[1:] = np.cumsum(keep_counts)
    kept_total = int(kept_offsets[-1])
    kept_targets = np.empty(kept_total, dtype=np.int32)
    kept_scores = np.empty(kept_total, dtype=float)

    dest_start = 0
    for row_idx in range(edges.n_nodes):
        keep_count = int(keep_counts[row_idx])
        if keep_count == 0:
            continue
        src_start = int(row_offsets[row_idx])
        src_end = src_start + keep_count
        dest_end = dest_start + keep_count
        kept_targets[dest_start:dest_end] = targets[src_start:src_end]
        kept_scores[dest_start:dest_end] = scores[src_start:src_end]
        dest_start = dest_end

    return CompactNeighborRankings(
        offsets=kept_offsets,
        target=kept_targets,
        score=kept_scores,
    )


def _build_sparse_neighbor_rankings(matrix: 'DataMatrix', max_k: Optional[int] = None) -> CompactNeighborRankings:
    edges = matrix.sorted_undirected_edges(skip_zeros=True, agg=max)
    return _build_neighbor_rankings_from_sorted_edges(edges, max_k=max_k)


def build_symmetric_neighbor_rankings(matrix: Union['DataMatrix', SortedUndirectedEdges], max_k: Optional[int] = None) -> CompactNeighborRankings:
    """Return deterministic per-row neighbor rankings using max(row,col) scores."""
    if isinstance(matrix, SortedUndirectedEdges):
        return _build_neighbor_rankings_from_sorted_edges(matrix, max_k=max_k)
    if scipy.sparse.issparse(matrix.data):
        return _build_sparse_neighbor_rankings(matrix, max_k=max_k)
    return _build_dense_neighbor_rankings(matrix.data, max_k=max_k)


def symmetric_knn_edge_index_dict(matrix: 'DataMatrix', k: int, lower_bound: float = 0, include_equal: bool = False,
                                  neighbor_rankings: Optional[NeighborRankings] = None) -> Dict[Tuple[int, int], float]:
    """Return OR-symmetric kNN edges keyed by undirected index pairs."""
    if k < 1:
        raise ValueError("k must be >= 1")

    if neighbor_rankings is None:
        neighbor_rankings = build_symmetric_neighbor_rankings(matrix, max_k=k)

    edge_dict = dict()
    for source_idx, row_target, row_score in _iter_neighbor_ranking_rows(neighbor_rankings):
        selected = 0
        for target_idx, score in zip(row_target, row_score):
            if not _score_passes_lower_bound(score, lower_bound, include_equal=include_equal):
                break
            edge = (source_idx, target_idx) if source_idx < target_idx else (target_idx, source_idx)
            edge_dict[edge] = score
            selected += 1
            if selected >= k:
                break
    return edge_dict


def mst_edge_index_dict(tree: 'MaxTree', lower_bound: float = 0, include_equal: bool = False) -> Dict[Tuple[int, int], float]:
    """Return MST edges passing the threshold keyed by undirected index pair."""
    edge_dict = dict()
    for source_idx, target_idx, score in tree.mst_edges:
        if _score_passes_lower_bound(score, lower_bound, include_equal=include_equal):
            edge = (source_idx, target_idx) if source_idx < target_idx else (target_idx, source_idx)
            edge_dict[edge] = score
    return edge_dict


def mst_knn_edge_index_dict(matrix: 'DataMatrix', k: int, lower_bound: float = 0, include_equal: bool = False,
                            tree: Optional['MaxTree'] = None,
                            neighbor_rankings: Optional[NeighborRankings] = None) -> Dict[Tuple[int, int], float]:
    """Return the union of MST and OR-symmetric kNN edges keyed by undirected index pair."""
    if tree is None:
        tree = MaxTree(matrix)

    edge_dict = mst_edge_index_dict(tree, lower_bound=lower_bound, include_equal=include_equal)
    edge_dict.update(symmetric_knn_edge_index_dict(
        matrix,
        k,
        lower_bound=lower_bound,
        include_equal=include_equal,
        neighbor_rankings=neighbor_rankings,
    ))
    return edge_dict


def sorted_edges_from_edge_index_dict(n_nodes: int, edge_dict: Dict[Tuple[int, int], float]) -> SortedUndirectedEdges:
    """Return a descending score-sorted edge table from an undirected edge dict."""
    if len(edge_dict) == 0:
        return SortedUndirectedEdges(
            n_nodes=n_nodes,
            source=np.empty(0, dtype=np.int32),
            target=np.empty(0, dtype=np.int32),
            score=np.empty(0, dtype=float),
        )

    items = sorted(
        ((float(score), int(source_idx), int(target_idx)) for (source_idx, target_idx), score in edge_dict.items()),
        key=lambda item: (-item[0], item[1], item[2]),
    )

    return SortedUndirectedEdges(
        n_nodes=n_nodes,
        source=np.fromiter((source_idx for _, source_idx, _ in items), dtype=np.int32, count=len(items)),
        target=np.fromiter((target_idx for _, _, target_idx in items), dtype=np.int32, count=len(items)),
        score=np.fromiter((score for score, _, _ in items), dtype=float, count=len(items)),
    )


class StreamingMstKnnAccumulator:
    """Compute MST ∪ OR-symmetric kNN edges from a stream of edges with bounded memory.

    This is the streaming counterpart of
    :func:`~domainator.transform_matrix.apply_mst_knn_sparsification`. Instead of
    materializing the full all-vs-all matrix and pruning afterwards, edges are fed in
    one at a time via :meth:`add_edge` and only bounded state is retained:

    * **MST** — exploiting the matroid property ``MSF(E1 ∪ E2) = MSF(MSF(E1) ∪ E2)``,
      only the running maximum spanning forest (≤ ``n_nodes - 1`` edges) plus a bounded
      buffer of unprocessed edges is kept. When the buffer fills, the forest of
      ``running MSF ∪ buffer`` is recomputed with the existing batch :class:`MaxTree`
      (Kruskal/union-find) and the buffer is cleared. This is **exact**: the final
      forest equals a single batch ``MaxTree`` over all edges (ties among equal-weight
      edges may resolve to a different but equally-valid forest).

    * **kNN** — a per-node adjacency dict ``node -> {neighbor: [out_score, in_score]}``
      where ``out_score`` is ``M[node, neighbor]`` and ``in_score`` is
      ``M[neighbor, node]``; the symmetric ranking score is ``max(out_score, in_score)``,
      matching :func:`build_symmetric_neighbor_rankings`. Each node's dict is trimmed to
      ``knn_soft_cap`` neighbors (by descending symmetric score) when it grows past
      ``2 * knn_soft_cap``. For symmetric input (e.g. ``compare_contigs``) this is exact;
      for slightly-asymmetric input (e.g. ``seq_dist`` diamond bit scores) it is exact as
      long as no node retains more than ``knn_soft_cap`` neighbors above its true top-k —
      effectively always, except for pathological high-degree hubs.

    The stored ``[out_score, in_score]`` also lets :meth:`to_csr` reproduce the directional
    asymmetry of the batch path for the kept edges. If trimming has dropped both directions
    of a kept (MST-only) edge, the symmetric MST score is written to both cells as a fallback.
    """

    def __init__(self, n_nodes: int, k: int, lower_bound: float = 0.0, include_equal: bool = False,
                 mst_buffer_cap: Optional[int] = None, knn_soft_cap: Optional[int] = None):
        if k < 1:
            raise ValueError("k must be >= 1")
        self.n_nodes = n_nodes
        self.k = k
        self.lower_bound = lower_bound
        self.include_equal = include_equal
        self.mst_buffer_cap = max(1, mst_buffer_cap) if mst_buffer_cap is not None else max(4 * n_nodes, 100_000)
        # A soft cap below k could evict a true top-k neighbor, so never trim below k.
        self.knn_soft_cap = max(k, knn_soft_cap) if knn_soft_cap is not None else max(4 * k, k + 16)

        self._mst_edges: List[Tuple[int, int, float]] = []  # running maximum spanning forest
        self._mst_buffer: List[Tuple[int, int, float]] = []
        self._adj: List[Dict[int, List[float]]] = [dict() for _ in range(n_nodes)]  # node -> {nbr: [out, in]}
        self._diagonal: Dict[int, float] = dict()

    def add_edge(self, i: int, j: int, score: float) -> None:
        """Feed one directed edge ``i -> j`` (value ``M[i, j]``) into the accumulator."""
        score = float(score)
        if i == j:
            if score != 0:
                previous = self._diagonal.get(i)
                if previous is None or score > previous:
                    self._diagonal[i] = score
            return
        if not _score_passes_lower_bound(score, self.lower_bound, include_equal=self.include_equal):
            return

        self._mst_buffer.append((i, j, score))
        if len(self._mst_buffer) >= self.mst_buffer_cap:
            self._flush_mst()

        self._update_adj(i, j, score, direction=0)
        self._update_adj(j, i, score, direction=1)

    def _update_adj(self, node: int, neighbor: int, score: float, direction: int) -> None:
        row = self._adj[node]
        entry = row.get(neighbor)
        if entry is None:
            entry = [0.0, 0.0]
            row[neighbor] = entry
        if score > entry[direction]:
            entry[direction] = score
        if len(row) > 2 * self.knn_soft_cap:
            self._trim_adj(node)

    def _trim_adj(self, node: int) -> None:
        row = self._adj[node]
        # Keep the knn_soft_cap neighbors with the highest symmetric score.
        ranked = sorted(row.items(), key=lambda item: (-max(item[1]), item[0]))
        self._adj[node] = dict(ranked[:self.knn_soft_cap])

    def _flush_mst(self) -> None:
        if not self._mst_buffer and not self._mst_edges:
            return
        combined = self._mst_edges + self._mst_buffer
        # Drop any zero/below-threshold edges ourselves so MaxTree(skip_zeros=True) emits no warning.
        items = sorted(
            ((float(s), int(a), int(b)) for (a, b, s) in combined if s != 0),
            key=lambda item: (-item[0], item[1], item[2]),
        )
        if not items:
            self._mst_edges = []
            self._mst_buffer = []
            return
        edges = SortedUndirectedEdges(
            n_nodes=self.n_nodes,
            source=np.fromiter((a for _, a, _ in items), dtype=np.int32, count=len(items)),
            target=np.fromiter((b for _, _, b in items), dtype=np.int32, count=len(items)),
            score=np.fromiter((s for s, _, _ in items), dtype=float, count=len(items)),
        )
        tree = MaxTree(edges, skip_zeros=True)
        self._mst_edges = list(tree.iter_mst_edges())
        self._mst_buffer = []

    def finalize(self) -> Dict[Tuple[int, int], float]:
        """Return the union of MST and OR-symmetric kNN edges keyed by canonical ``(min, max)`` pair."""
        self._flush_mst()

        edge_dict: Dict[Tuple[int, int], float] = dict()
        for source_idx, target_idx, score in self._mst_edges:
            if source_idx == target_idx:
                continue
            edge = (source_idx, target_idx) if source_idx < target_idx else (target_idx, source_idx)
            # The MST tuple carries the surviving (max-direction) weight; trust it even if
            # adjacency trimming has since dropped this pair from both neighbor dicts.
            edge_dict[edge] = score

        for node, row in enumerate(self._adj):
            # Rank this node's neighbors exactly like symmetric_knn_edge_index_dict:
            # descending symmetric score, ties broken by target index.
            ranked = sorted(((max(entry), neighbor) for neighbor, entry in row.items()),
                            key=lambda item: (-item[0], item[1]))
            selected = 0
            for score, neighbor in ranked:
                if not _score_passes_lower_bound(score, self.lower_bound, include_equal=self.include_equal):
                    break
                edge = (node, neighbor) if node < neighbor else (neighbor, node)
                edge_dict[edge] = score
                selected += 1
                if selected >= self.k:
                    break

        return edge_dict

    def _directional_values(self, i: int, j: int, fallback: float) -> Tuple[float, float]:
        """Return ``(M[i, j], M[j, i])`` for a kept pair, falling back to ``fallback``."""
        entry = self._adj[i].get(j)
        if entry is not None:
            forward, reverse = entry[0], entry[1]
        else:
            mirror = self._adj[j].get(i)
            if mirror is not None:
                forward, reverse = mirror[1], mirror[0]
            else:
                forward = reverse = 0.0
        if forward == 0.0 and reverse == 0.0:
            # Both directions trimmed (only possible for a pure MST edge): use its weight.
            forward = reverse = fallback
        return forward, reverse

    def to_csr(self, dtype=np.float64) -> scipy.sparse.csr_array:
        """Assemble the pruned graph as a sparse CSR matrix.

        Off-diagonal kept edges are written in both directions (using the retained
        directional values ``M[i, j]`` / ``M[j, i]`` where available); the diagonal is
        populated from self-edges seen during streaming. Mirrors the sparse branch of
        :func:`~domainator.transform_matrix.apply_mst_knn_sparsification`.
        """
        edge_dict = self.finalize()
        out = scipy.sparse.dok_array((self.n_nodes, self.n_nodes), dtype=dtype)
        for index, value in self._diagonal.items():
            if value != 0:
                out[index, index] = value
        for (i, j), score in edge_dict.items():
            forward, reverse = self._directional_values(i, j, score)
            if forward != 0:
                out[i, j] = forward
            if reverse != 0:
                out[j, i] = reverse
        return scipy.sparse.csr_array(out)


def mst_knn_edge_counts_by_threshold(matrix: 'DataMatrix', tree: 'MaxTree', max_k: int,
                                     include_equal: bool = True,
                                     neighbor_rankings: Optional[NeighborRankings] = None) -> np.ndarray:
    """Return exact MST-kNN edge counts for each MST threshold and k in [2, max_k]."""
    if max_k < 2:
        raise ValueError("max_k must be >= 2")

    if neighbor_rankings is None:
        neighbor_rankings = build_symmetric_neighbor_rankings(matrix, max_k=max_k)

    thresholds = tree.edges_by_threshold[:, 1]
    counts = np.zeros((len(thresholds), max_k - 1), dtype=int)
    if len(thresholds) == 0:
        return counts

    if isinstance(neighbor_rankings, CompactNeighborRankings):
        total_events = len(neighbor_rankings.target)
    else:
        total_events = sum(min(len(ranked_neighbors), max_k) for ranked_neighbors in neighbor_rankings)
    # Pack all activation fields into one structured array so the sort is in-place
    # and avoids a separate O(total_events) argsort index array.
    act_dtype = np.dtype([('score', np.float64), ('rank', np.int16),
                          ('source', np.int32), ('target', np.int32)])
    activations = np.empty(total_events, dtype=act_dtype)

    activation_idx = 0
    for source_idx, row_target, row_score in _iter_neighbor_ranking_rows(neighbor_rankings):
        for rank, (target_idx, score) in enumerate(zip(row_target[:max_k], row_score[:max_k]), start=1):
            if source_idx < target_idx:
                edge_source, edge_target = source_idx, target_idx
            else:
                edge_source, edge_target = target_idx, source_idx
            activations['score'][activation_idx] = float(score)
            activations['rank'][activation_idx] = rank
            activations['source'][activation_idx] = edge_source
            activations['target'][activation_idx] = edge_target
            activation_idx += 1

    if total_events > 0:
        activations.sort(order='score', kind='stable')  # ascending in-place
        activations = activations[::-1]                  # descending view

    active_edge_min_rank = dict()
    non_mst_rank_counts = np.zeros(max_k + 1, dtype=int)
    mst_prefix_edges = set()
    activation_idx = 0

    def edge_key(source_idx: int, target_idx: int) -> int:
        return int(source_idx) * tree.n_nodes + int(target_idx)

    for threshold_idx, mst_edge in enumerate(tree.mst_edges):
        threshold = thresholds[threshold_idx]

        while activation_idx < total_events and _score_passes_lower_bound(activations['score'][activation_idx], threshold, include_equal=include_equal):
            rank = int(activations['rank'][activation_idx])
            edge = edge_key(activations['source'][activation_idx], activations['target'][activation_idx])
            current_rank = active_edge_min_rank.get(edge)

            if current_rank is None:
                active_edge_min_rank[edge] = rank
                if edge not in mst_prefix_edges:
                    non_mst_rank_counts[rank] += 1
            elif rank < current_rank:
                active_edge_min_rank[edge] = rank
                if edge not in mst_prefix_edges:
                    non_mst_rank_counts[current_rank] -= 1
                    non_mst_rank_counts[rank] += 1

            activation_idx += 1

        source_idx, target_idx, _ = mst_edge
        edge = edge_key(source_idx, target_idx) if source_idx < target_idx else edge_key(target_idx, source_idx)
        mst_prefix_edges.add(edge)
        current_rank = active_edge_min_rank.get(edge)
        if current_rank is not None:
            non_mst_rank_counts[current_rank] -= 1

        cumulative_non_mst = np.cumsum(non_mst_rank_counts)
        counts[threshold_idx, :] = len(mst_prefix_edges) + cumulative_non_mst[2:]

    return counts


class MaxTree():
    """Maximum Spanning Tree (MST) from a DataMatrix.

    A Maximum Spanning Tree is a subgraph that connects all nodes with the maximum
    total edge weight while avoiding cycles. This is useful for clustering and
    visualization of similarity/distance matrices.

    The MST is computed using Kruskal's algorithm with union-find for cycle detection.
    Edges are processed from highest to lowest weight, and only edges connecting
    different components are added to the tree.

    Use Cases:
        - Hierarchical clustering: Cut tree at different thresholds to create clusters
        - Network visualization: Display most significant relationships
        - Sequence similarity networks: Connect most similar sequences

    Args:
        matrix: DataMatrix or sorted edge table containing edge weights
        skip_zeros: If True, zero-weight edges are excluded. When there is no connected
                   MST without zero edges, the result will be a forest (multiple trees).

    Attributes:
        n_nodes: Number of nodes in the graph
        edges: Array of shape (n_edges, 4) containing [node_i, node_j, weight, gap_count]
               where gap_count is the number of edges between consecutive MST edges
        mst: Array of indices indicating which edges are in the MST
        edges_by_threshold: Array of (edge_count, threshold) pairs for MST edges
        cluster_count_by_threshold: Array of (threshold, cluster_count) showing how
                                   clusters form at different thresholds
        cluster_count_by_edge_count: Array of (edge_count, cluster_count) including
                                    non-MST edges between MST edges

    Properties:
        mst_edges: List of (node_i, node_j, weight) tuples for MST edges only

    Methods:
        export_for_interactive_viz: Export minimal structure for client-side clustering

    Example:
        >>> matrix = DataMatrix.from_file("similarities.hdf5")
        >>> tree = MaxTree(matrix, skip_zeros=True)
        >>> # Get clusters at a specific threshold
        >>> mst_edges = tree.mst_edges
        >>> # Find how many clusters at different thresholds
        >>> threshold_clusters = tree.cluster_count_by_threshold

    Note:
        When skip_zeros=True and the matrix is disconnected (contains at least one pair of nodes with no path),
        the resulting MST will be a forest of multiple disconnected trees, one per
        connected component.
    """

    def __init__(self, matrix: Union['DataMatrix', SortedUndirectedEdges], skip_zeros=True):
        from .data_matrix import DataMatrix  # local import to avoid an import cycle

        if isinstance(matrix, DataMatrix):
            sorted_edges = matrix.sorted_undirected_edges(skip_zeros=skip_zeros, agg=max)
        else:
            if skip_zeros:
                keep_mask = matrix.score != 0
                sorted_edges = SortedUndirectedEdges(
                    n_nodes=matrix.n_nodes,
                    source=matrix.source[keep_mask],
                    target=matrix.target[keep_mask],
                    score=matrix.score[keep_mask],
                )
            else:
                sorted_edges = matrix

        if not skip_zeros and len(sorted_edges) > 0 and np.any(sorted_edges.score == 0):
            warnings.warn(
                "MaxTree(skip_zeros=False) includes zero-score edges and may connect otherwise disconnected components.",
                RuntimeWarning,
                stacklevel=2,
            )
        elif not isinstance(matrix, DataMatrix) and skip_zeros and len(matrix) != len(sorted_edges):
            warnings.warn(
                "MaxTree(skip_zeros=True) removed zero-score edges from the provided SortedUndirectedEdges.",
                RuntimeWarning,
                stacklevel=2,
            )

        self.n_nodes = sorted_edges.n_nodes

        # calculate_max_spanning_tree using union-find
        parent = np.arange(self.n_nodes, dtype=int)  # Union-find parent array

        def find(x):
            """Find root with path compression"""
            root = x

            while parent[root] != root:
                root = parent[root]

            while parent[x] != x:
                next_x = parent[x]
                parent[x] = root
                x = next_x

            return root

        def union(x, y):
            """Union two sets"""
            root_x = find(x)
            root_y = find(y)
            if root_x != root_y:
                parent[root_x] = root_y
                return True
            return False

        mst_edges = []
        edge_count = 0
        for edge_i in range(len(sorted_edges)):
            node_1 = int(sorted_edges.source[edge_i])
            node_2 = int(sorted_edges.target[edge_i])
            score = float(sorted_edges.score[edge_i])
            edge_count += 1

            # Only add edge if it connects different components
            if union(node_1, node_2):
                mst_edges.append((node_1, node_2, score, edge_count))
                edge_count = 0

        if len(mst_edges) > 0:
            self.edges = np.asarray(mst_edges, dtype=float)
        else:
            self.edges = np.zeros((0, 4), dtype=float)
        self.mst = np.arange(len(self.edges), dtype=int)

        self.edges_by_threshold = self._edges_by_threshold()
        self.cluster_count_by_threshold = self._cluster_count_by_threshold()
        self.cluster_count_by_edge_count = self._cluster_count_by_edge_count()

    @property
    def mst_edges(self) -> List[Tuple[int,int,float]]: # (node_i, node_j, edge_value)
        """
        Returns the edges in the maximum spanning tree.

        Returns:
            List[Tuple[int,int,float]]: List of MST edges as (node_i, node_j, edge_value) tuples
        """
        return list(self.iter_mst_edges())

    def iter_mst_edges(self) -> Iterator[Tuple[int, int, float]]:
        """Iterate over the edges in the maximum spanning tree."""
        for mst_idx in self.mst:
            edge = self.edges[mst_idx]
            yield int(edge[0]), int(edge[1]), float(edge[2])


    def _edges_by_threshold(self):
        """
        Returns array of (edge_count, threshold) pairs.
        Each row represents the cumulative edge count and threshold at each MST edge.
        """
        edges_by_threshold = np.zeros((len(self.mst), 2)) # edge_count, threshold
        edges_count = 0
        for i in range(len(self.mst)):
            mst_edge_idx = self.mst[i]
            edges_count += int(self.edges[mst_edge_idx, 3])
            threshold = self.edges[mst_edge_idx, 2]
            edges_by_threshold[i, 0] = edges_count
            edges_by_threshold[i, 1] = threshold
        return edges_by_threshold

    def _cluster_count_by_threshold(self):
        """
        Returns array of (threshold, cluster_count) pairs.
        For each MST edge threshold, calculates how many clusters would result
        if we only include edges with values >= that threshold.
        """
        cluster_counts = np.zeros((len(self.mst) + 1, 2))  # threshold, cluster_count

        # Start with all nodes as separate clusters
        cluster_counts[0, 0] = float('inf')  # threshold = infinity
        cluster_counts[0, 1] = self.n_nodes  # all nodes separate

        # Process MST edges from highest to lowest threshold
        for i in range(len(self.mst)):
            mst_edge_idx = self.mst[i]
            threshold = self.edges[mst_edge_idx, 2]
            # Each MST edge reduces cluster count by 1
            cluster_counts[i + 1, 0] = threshold
            cluster_counts[i + 1, 1] = self.n_nodes - (i + 1)

        return cluster_counts

    def _cluster_count_by_edge_count(self):
        """
        Returns array of (edge_count, cluster_count) pairs.
        For each cumulative edge count, calculates the number of clusters.
        This includes non-MST edges between MST edges.
        """
        edges_by_thresh = self.edges_by_threshold

        if len(edges_by_thresh) == 0:
            return np.array([[0, self.n_nodes]])

        # Calculate cluster counts
        result = []
        result.append([0, self.n_nodes])  # Start with 0 edges, all nodes separate

        for i in range(len(edges_by_thresh)):
            edge_count = int(edges_by_thresh[i, 0])
            # Each MST edge reduces cluster count by 1
            cluster_count = self.n_nodes - (i + 1)
            result.append([edge_count, cluster_count])

        return np.array(result)

    def export_for_interactive_viz(self):
        """
        Export minimal data structure for client-side cluster computation.
        Returns dict with MST edges for JavaScript to rebuild clusters on-demand.

        Returns:
            dict with:
            - n_nodes: int
            - mst_edges: [[node1, node2, value], ...] sorted by value descending
        """
        mst_edges = []
        for i in range(len(self.mst)):
            edge_idx = self.mst[i]
            mst_edges.append([
                int(self.edges[edge_idx, 0]),
                int(self.edges[edge_idx, 1]),
                float(self.edges[edge_idx, 2])
            ])

        return {
            'n_nodes': int(self.n_nodes),
            'mst_edges': mst_edges
        }
