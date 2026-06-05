"""Pure data structures and helpers for undirected-edge / neighbor-ranking tables.

These types have no dependency on :class:`~domainator.data_matrix.DataMatrix` and are
imported by both ``data_matrix`` and ``ssn_edges`` to avoid an import cycle.
"""
from dataclasses import dataclass
from typing import List, Tuple, Union

import numpy as np


@dataclass(frozen=True)
class SortedUndirectedEdges:
    """Compact representation of sorted undirected edges."""

    n_nodes: int
    source: np.ndarray
    target: np.ndarray
    score: np.ndarray

    def __post_init__(self):
        if self.n_nodes < 0:
            raise ValueError("n_nodes must be >= 0")
        if self.source.ndim != 1 or self.target.ndim != 1 or self.score.ndim != 1:
            raise ValueError("source, target, and score arrays must be 1-dimensional")
        if not (len(self.source) == len(self.target) == len(self.score)):
            raise ValueError("source, target, and score arrays must have the same length")
        if len(self.source) == 0:
            return
        if np.any(self.source < 0) or np.any(self.target < 0):
            raise ValueError("source and target indices must be >= 0")
        if np.any(self.source >= self.n_nodes) or np.any(self.target >= self.n_nodes):
            raise ValueError("source and target indices must be less than n_nodes")
        if np.any(self.source == self.target):
            raise ValueError("self edges are not supported")
        if np.any(self.score[:-1] < self.score[1:]):
            raise ValueError("score array must be sorted in descending order")

    def __len__(self) -> int:
        return len(self.score)


@dataclass(frozen=True)
class CompactNeighborRankings:
    """Array-backed per-row neighbor rankings."""

    offsets: np.ndarray
    target: np.ndarray
    score: np.ndarray

    def __post_init__(self):
        if self.offsets.ndim != 1 or self.target.ndim != 1 or self.score.ndim != 1:
            raise ValueError("offsets, target, and score arrays must be 1-dimensional")
        if len(self.offsets) == 0:
            raise ValueError("offsets must contain at least one element")
        if self.offsets[0] != 0:
            raise ValueError("offsets must start at 0")
        if np.any(self.offsets < 0):
            raise ValueError("offsets must be non-negative")
        if np.any(self.offsets[:-1] > self.offsets[1:]):
            raise ValueError("offsets must be monotonically non-decreasing")
        if self.offsets[-1] != len(self.target) or self.offsets[-1] != len(self.score):
            raise ValueError("offsets must end at the length of target and score arrays")
        if len(self.target) > 0:
            if np.any(self.target < 0):
                raise ValueError("target indices must be >= 0")
            if np.any(self.target >= len(self)):
                raise ValueError("target indices must be less than the number of rows")

    def __len__(self) -> int:
        return len(self.offsets) - 1

    def row_bounds(self, row_idx: int) -> Tuple[int, int]:
        return int(self.offsets[row_idx]), int(self.offsets[row_idx + 1])


NeighborRankings = Union[CompactNeighborRankings, List[List[Tuple[int, float]]]]


def _empty_compact_neighbor_rankings(n_rows: int) -> CompactNeighborRankings:
    return CompactNeighborRankings(
        offsets=np.zeros(n_rows + 1, dtype=np.int64),
        target=np.empty(0, dtype=np.int32),
        score=np.empty(0, dtype=float),
    )


def _iter_neighbor_ranking_rows(neighbor_rankings: NeighborRankings):
    if isinstance(neighbor_rankings, CompactNeighborRankings):
        for row_idx in range(len(neighbor_rankings)):
            start, end = neighbor_rankings.row_bounds(row_idx)
            yield row_idx, neighbor_rankings.target[start:end], neighbor_rankings.score[start:end]
        return

    for row_idx, ranked_neighbors in enumerate(neighbor_rankings):
        if len(ranked_neighbors) == 0:
            yield row_idx, np.empty(0, dtype=np.int32), np.empty(0, dtype=float)
            continue
        row_target = np.fromiter((target_idx for target_idx, _ in ranked_neighbors), dtype=np.int32, count=len(ranked_neighbors))
        row_score = np.fromiter((score for _, score in ranked_neighbors), dtype=float, count=len(ranked_neighbors))
        yield row_idx, row_target, row_score


def _score_passes_lower_bound(score: float, lower_bound: float, include_equal: bool = False) -> bool:
    if include_equal:
        return score >= lower_bound
    return score > lower_bound
