"""Sparse score matrices built in chunks of query rows, shared by the matrix-producing tools."""
import math

import numpy as np
import scipy.sparse

_DEFAULT_CHUNKS_PER_WORKER = 16
_MAX_QUERY_ROWS_PER_CHUNK = 1024


def iter_query_chunks(n_rows, cpu, chunks_per_worker=_DEFAULT_CHUNKS_PER_WORKER, max_rows_per_chunk=_MAX_QUERY_ROWS_PER_CHUNK):
    if n_rows == 0:
        return []

    requested_workers = max(cpu, 1)
    chunk_size = max(1, math.ceil(n_rows / (requested_workers * chunks_per_worker)))
    if max_rows_per_chunk is not None:
        chunk_size = min(chunk_size, max_rows_per_chunk)
    return [(start, min(start + chunk_size, n_rows)) for start in range(0, n_rows, chunk_size)]


class StreamingCSRChunkBuilder:
    def __init__(self, shape, dtype=np.float64):
        self.shape = shape
        self.dtype = np.dtype(dtype)
        self._row_count = 0
        self._nnz = 0
        self._indices_chunks = list()
        self._data_chunks = list()
        self._indptr = [0]

    def append_chunk(self, start, end, chunk_matrix):
        if start != self._row_count:
            raise RuntimeError(
                f"Similarity chunks arrived out of order: expected row {self._row_count}, got {start}."
            )

        chunk_matrix = scipy.sparse.csr_array(chunk_matrix, dtype=self.dtype)
        expected_shape = (end - start, self.shape[1])
        if chunk_matrix.shape != expected_shape:
            raise RuntimeError(
                f"Similarity chunk shape mismatch for rows {start}:{end}: expected {expected_shape}, got {chunk_matrix.shape}."
            )

        if chunk_matrix.nnz > 0:
            self._indices_chunks.append(chunk_matrix.indices.astype(np.int32, copy=False))
            self._data_chunks.append(chunk_matrix.data.astype(self.dtype, copy=False))
        self._indptr.extend((chunk_matrix.indptr[1:] + self._nnz).tolist())
        self._nnz += chunk_matrix.nnz
        self._row_count = end

    def build(self):
        if self._row_count != self.shape[0]:
            raise RuntimeError(
                f"Similarity chunk assembly incomplete: expected {self.shape[0]} rows, assembled {self._row_count}."
            )

        if self._nnz == 0:
            return scipy.sparse.csr_array(self.shape, dtype=self.dtype)

        return scipy.sparse.csr_array(
            (
                np.concatenate(self._data_chunks),
                np.concatenate(self._indices_chunks),
                np.asarray(self._indptr, dtype=np.int64),
            ),
            shape=self.shape,
        )


def prune_scores_inplace(scores, threshold):
    if threshold <= 0 or scores.nnz == 0:
        return scores
    scores = scipy.sparse.csr_array(scores)
    scores.data[scores.data <= threshold] = 0
    scores.eliminate_zeros()
    return scores


def feed_chunk_to_accumulator(accumulator, start, chunk_matrix):
    """Feed the nonzero entries of a per-row chunk into a StreamingMstKnnAccumulator."""
    coo = scipy.sparse.coo_array(chunk_matrix)
    for local_row, col, value in zip(coo.row, coo.col, coo.data):
        accumulator.add_edge(int(start) + int(local_row), int(col), float(value))


def keep_top_k_per_row(scores, k):
    if k is None:
        return scores

    matrix = scipy.sparse.csr_array(scores)
    if k == 0 or matrix.nnz == 0:
        return scipy.sparse.csr_array(matrix.shape, dtype=matrix.dtype)

    row_idx = list()
    col_idx = list()
    data = list()

    for query_idx in range(matrix.shape[0]):
        row_start = matrix.indptr[query_idx]
        row_end = matrix.indptr[query_idx + 1]
        target_indices = matrix.indices[row_start:row_end]
        row_scores = matrix.data[row_start:row_end]

        if row_scores.size == 0:
            continue

        if row_scores.size > k:
            top_order = np.lexsort((target_indices, -row_scores))[:k]
            target_indices = target_indices[top_order]
            row_scores = row_scores[top_order]

        column_order = np.argsort(target_indices)
        target_indices = target_indices[column_order]
        row_scores = row_scores[column_order]

        row_idx.extend([query_idx] * target_indices.size)
        col_idx.extend(target_indices.tolist())
        data.extend(row_scores.tolist())

    if len(data) == 0:
        return scipy.sparse.csr_array(matrix.shape, dtype=matrix.dtype)

    return scipy.sparse.csr_array((np.asarray(data), (np.asarray(row_idx), np.asarray(col_idx))), shape=matrix.shape)
