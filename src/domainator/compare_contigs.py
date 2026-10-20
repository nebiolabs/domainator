"""Calculates similarity metrics for gene neighborhoods

Takes a genbank file annotated with domainator and calculates similarities between contigs and writes reports.

ji, ai, and dss metrics are similar to and inspired by the metrics used in the BigScape method:
Navarro-Muñoz, Jorge C., Nelly Selem-Mojica, Michael W. Mullowney, Satria A. Kautsar, James H. Tryon, Elizabeth I. Parkinson, Emmanuel L. C. De Los Santos, et al. "A Computational Framework to Explore Large-Scale Biosynthetic Diversity from Large-Scale Genomic Data." Nature Chemical Biology 16, no. 1 (January 2020): 60–68. https://doi.org/10.1038/s41589-019-0400-9.
"""
import warnings
warnings.filterwarnings("ignore", module='numpy')
import sys
from jsonargparse import ArgumentParser, ActionConfigFile
from domainator.utils import list_and_file_to_dict_keys, parse_seqfiles, write_genbank, copy_SeqRecord, make_pool
from abc import ABC, abstractmethod
from collections import OrderedDict
import scipy.cluster.hierarchy
import math
from domainator import __version__, DOMAIN_FEATURE_NAME, RawAndDefaultsFormatter
from domainator.data_matrix import DataMatrix, StreamingMstKnnAccumulator
from domainator.sparsify import MAX_HITS_PER_QUERY_SPARSIFY_NOTE, add_sparsify_arguments, sparsify_matrix, sparsify_requested
from domainator.filter_domains import filter_domains
from domainator.matrix_output import add_matrix_output_arguments, resolve_matrix_outputs, write_matrix_outputs
from domainator.output_guardrails import add_max_output_gb_argument, format_bytes, max_output_gb_to_bytes, OutputSizeLimitExceeded
from typing import Optional, Set
import numpy as np
import scipy.sparse
import psutil
import tqdm

_SIMILARITY_FEATURE_MATRIX = None
_SIMILARITY_ROW_COUNTS = None
_PREPARED_SIMILARITY_METRICS = None
_DEFAULT_CHUNKS_PER_WORKER = 16
_MAX_QUERY_ROWS_PER_CHUNK = 1024


def _progress(iterable, enabled, **kwargs):
    if not enabled:
        return iterable
    return tqdm.tqdm(iterable, leave=True, dynamic_ncols=True, **kwargs)


def _build_sparse_feature_matrix(contig_data):
    contig_ids = list(contig_data.keys())
    feature_to_idx = dict()
    row_idx = list()
    col_idx = list()

    for contig_idx, features in enumerate(contig_data.values()):
        for feature in features:
            if feature not in feature_to_idx:
                feature_to_idx[feature] = len(feature_to_idx)
            row_idx.append(contig_idx)
            col_idx.append(feature_to_idx[feature])

    if len(feature_to_idx) == 0:
        feature_matrix = scipy.sparse.csr_array((len(contig_ids), 0), dtype=np.float64)
    else:
        data = np.ones(len(row_idx), dtype=np.float64)
        feature_matrix = scipy.sparse.csr_array(
            (data, (row_idx, col_idx)),
            shape=(len(contig_ids), len(feature_to_idx)),
        )

    return contig_ids, feature_matrix


class _PreparedSimilarityMetric:
    def __init__(self, weight, feature_matrix, label):
        self.weight = weight
        self.feature_matrix = scipy.sparse.csr_array(feature_matrix, dtype=np.float64)
        self.row_counts = np.diff(self.feature_matrix.indptr)
        self.label = label


def _iter_query_chunks(n_rows, cpu, chunks_per_worker=_DEFAULT_CHUNKS_PER_WORKER, max_rows_per_chunk=_MAX_QUERY_ROWS_PER_CHUNK):
    if n_rows == 0:
        return []

    requested_workers = max(cpu, 1)
    chunk_size = max(1, math.ceil(n_rows / (requested_workers * chunks_per_worker)))
    if max_rows_per_chunk is not None:
        chunk_size = min(chunk_size, max_rows_per_chunk)
    return [(start, min(start + chunk_size, n_rows)) for start in range(0, n_rows, chunk_size)]


def _init_similarity_chunk_worker(feature_matrix, row_counts):
    global _SIMILARITY_FEATURE_MATRIX, _SIMILARITY_ROW_COUNTS
    _SIMILARITY_FEATURE_MATRIX = feature_matrix
    _SIMILARITY_ROW_COUNTS = row_counts


def _init_combined_similarity_chunk_worker(prepared_metrics):
    global _PREPARED_SIMILARITY_METRICS
    _PREPARED_SIMILARITY_METRICS = prepared_metrics


class _StreamingCSRChunkBuilder:
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


def _compute_similarity_chunk(feature_matrix, row_counts, start, end):
    overlaps = (feature_matrix[start:end] @ feature_matrix.transpose()).tocsr()
    overlaps.sort_indices()

    row_idx = list()
    col_idx = list()
    data = list()

    for local_query_idx in range(end - start):
        query_idx = start + local_query_idx
        row_idx.append(local_query_idx)
        col_idx.append(query_idx)
        data.append(1.0)

        row_start = overlaps.indptr[local_query_idx]
        row_end = overlaps.indptr[local_query_idx + 1]
        target_indices = overlaps.indices[row_start:row_end]
        shared_counts = overlaps.data[row_start:row_end].astype(np.float64, copy=False)

        off_diagonal_mask = target_indices != query_idx
        target_indices = target_indices[off_diagonal_mask]
        shared_counts = shared_counts[off_diagonal_mask]

        if target_indices.size > 0:
            union_counts = row_counts[query_idx] + row_counts[target_indices] - shared_counts
            scores = shared_counts / union_counts
            positive_mask = scores > 0
            target_indices = target_indices[positive_mask]
            scores = scores[positive_mask]
        else:
            scores = np.empty(0, dtype=np.float64)

        row_idx.extend([local_query_idx] * len(target_indices))
        col_idx.extend(target_indices.tolist())
        data.extend(scores.tolist())

    return scipy.sparse.csr_array(
        (
            np.asarray(data, dtype=np.float64),
            (np.asarray(row_idx, dtype=np.int32), np.asarray(col_idx, dtype=np.int32)),
        ),
        shape=(end - start, feature_matrix.shape[0]),
    )


def _prune_scores_inplace(scores, threshold):
    if threshold <= 0 or scores.nnz == 0:
        return scores
    scores = scipy.sparse.csr_array(scores)
    scores.data[scores.data <= threshold] = 0
    scores.eliminate_zeros()
    return scores


def _compute_combined_similarity_chunk(prepared_metrics, start, end, k, lb):
    chunk_shape = (end - start, prepared_metrics[0].feature_matrix.shape[0])
    combined_chunk = scipy.sparse.csr_array(chunk_shape, dtype=np.float64)
    remaining_weight = sum(metric.weight for metric in prepared_metrics)

    for metric in prepared_metrics:
        metric_chunk = _compute_similarity_chunk(metric.feature_matrix, metric.row_counts, start, end)
        combined_chunk = combined_chunk + (metric.weight * metric_chunk)
        remaining_weight -= metric.weight
        if lb > 0:
            combined_chunk = _prune_scores_inplace(combined_chunk, lb - remaining_weight)

    combined_chunk = _keep_top_k_per_row(combined_chunk, k)
    if lb > 0:
        combined_chunk = _prune_scores_inplace(combined_chunk, lb)
    return combined_chunk


def _compute_similarity_chunk_worker(args):
    start, end = args
    return start, end, _compute_similarity_chunk(_SIMILARITY_FEATURE_MATRIX, _SIMILARITY_ROW_COUNTS, start, end)


def _compute_combined_similarity_chunk_worker(args):
    start, end, k, lb = args
    return start, end, _compute_combined_similarity_chunk(_PREPARED_SIMILARITY_METRICS, start, end, k, lb)


def _compute_sparse_similarity_matrix(contig_data, cpu=1, progress=False, desc=None):
    contig_ids, feature_matrix = _build_sparse_feature_matrix(contig_data)
    matrix_shape = (len(contig_ids), len(contig_ids))

    if feature_matrix.shape[1] == 0:
        return scipy.sparse.identity(matrix_shape[0], dtype=np.float64, format="csr")

    row_counts = np.diff(feature_matrix.indptr)
    query_chunks = list(_iter_query_chunks(matrix_shape[0], cpu))
    chunk_builder = _StreamingCSRChunkBuilder(matrix_shape)

    if cpu > 1 and len(query_chunks) > 1:
        with make_pool(
            processes=min(cpu, len(query_chunks)),
            initializer=_init_similarity_chunk_worker,
            initargs=(feature_matrix, row_counts),
        ) as pool:
            for start, end, chunk_matrix in _progress(pool.imap(_compute_similarity_chunk_worker, query_chunks), progress, total=len(query_chunks), desc=desc):
                chunk_builder.append_chunk(start, end, chunk_matrix)
    else:
        for start, end in _progress(query_chunks, progress, total=len(query_chunks), desc=desc):
            chunk_builder.append_chunk(start, end, _compute_similarity_chunk(feature_matrix, row_counts, start, end))

    return chunk_builder.build()


def _compute_combined_sparse_similarity_matrix(prepared_metrics, k, lb, cpu=1, progress=False, desc=None):
    if len(prepared_metrics) == 0:
        return scipy.sparse.csr_array((0, 0), dtype=np.float64)

    matrix_shape = (
        prepared_metrics[0].feature_matrix.shape[0],
        prepared_metrics[0].feature_matrix.shape[0],
    )
    if k == 0:
        return scipy.sparse.csr_array(matrix_shape, dtype=np.float64)

    query_chunks = [(start, end, k, lb) for start, end in _iter_query_chunks(matrix_shape[0], cpu)]
    chunk_builder = _StreamingCSRChunkBuilder(matrix_shape)

    if cpu > 1 and len(query_chunks) > 1:
        with make_pool(
            processes=min(cpu, len(query_chunks)),
            initializer=_init_combined_similarity_chunk_worker,
            initargs=(prepared_metrics,),
        ) as pool:
            for start, end, chunk_matrix in _progress(pool.imap(_compute_combined_similarity_chunk_worker, query_chunks), progress, total=len(query_chunks), desc=desc):
                chunk_builder.append_chunk(start, end, chunk_matrix)
    else:
        for start, end, chunk_k, chunk_lb in _progress(query_chunks, progress, total=len(query_chunks), desc=desc):
            chunk_builder.append_chunk(start, end, _compute_combined_similarity_chunk(prepared_metrics, start, end, chunk_k, chunk_lb))

    return chunk_builder.build()


def _feed_chunk_to_accumulator(accumulator, start, chunk_matrix):
    """Feed the nonzero entries of a per-row chunk into a StreamingMstKnnAccumulator."""
    coo = scipy.sparse.coo_array(chunk_matrix)
    for local_row, col, value in zip(coo.row, coo.col, coo.data):
        accumulator.add_edge(int(start) + int(local_row), int(col), float(value))


def _compute_sparsified_similarity_matrix(prepared_metrics, *, mst, knn, k, lb, cpu=1, progress=False, desc=None):
    """Stream per-row similarity chunks into an MST and/or kNN accumulator.

    Equivalent to building the full combined similarity matrix (with the per-row top-k
    cap and ``--lb`` applied) and then calling :func:`~domainator.sparsify.sparsify_matrix`,
    but the full matrix is never materialized: each chunk's edges are fed into the
    accumulator and discarded.
    """
    n = prepared_metrics[0].feature_matrix.shape[0]
    accumulator = StreamingMstKnnAccumulator(n, mst=mst, knn=knn, lower_bound=lb)

    query_chunks = [(start, end, k, lb) for start, end in _iter_query_chunks(n, cpu)]

    if cpu > 1 and len(query_chunks) > 1:
        with make_pool(
            processes=min(cpu, len(query_chunks)),
            initializer=_init_combined_similarity_chunk_worker,
            initargs=(prepared_metrics,),
        ) as pool:
            for start, end, chunk_matrix in _progress(pool.imap(_compute_combined_similarity_chunk_worker, query_chunks), progress, total=len(query_chunks), desc=desc):
                _feed_chunk_to_accumulator(accumulator, start, chunk_matrix)
    else:
        for start, end, chunk_k, chunk_lb in _progress(query_chunks, progress, total=len(query_chunks), desc=desc):
            _feed_chunk_to_accumulator(accumulator, start, _compute_combined_similarity_chunk(prepared_metrics, start, end, chunk_k, chunk_lb))

    return accumulator.to_csr()


def _keep_top_k_per_row(scores, k):
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

class ContigMetric(ABC):
    weight: float = 0.0

    def prepare_sparse_metric(self, recs):
        return None


    @abstractmethod
    def compute(self, recs, cpu=1, progress=False):
        """
            computes the metric matrix for the contigs and returns a sparse matrix.
        """
        pass

class JaccardIndex(ContigMetric):
    def __init__(self, weight: float = 0.0, databases=None):
        self.weight = weight
        self.databases = databases
    
    def get_contig_data(self, recs, databases=None):
        """
            recs: a list of SeqRecords
            databases: a set of database names to keep, if None, keep all databases

            store the domain content for the contig
        """
        out = OrderedDict()
        for contig in recs:
            if contig.id not in out:
                out[contig.id] = set()
            
            for feature in contig.features:
                if feature.type == DOMAIN_FEATURE_NAME:
                    if databases is not None and feature.qualifiers["database"][0] not in databases:
                        continue
                    domain_name = feature.qualifiers["name"][0]
                    out[contig.id].add(domain_name)

        return out


    def compute(self, recs, cpu=1, progress=False):
        contig_data = self.get_contig_data(recs, self.databases)
        return _compute_sparse_similarity_matrix(contig_data, cpu=cpu, progress=progress, desc=f"Scoring {self.__class__.__name__}")

    def prepare_sparse_metric(self, recs):
        contig_data = self.get_contig_data(recs, self.databases)
        _contig_ids, feature_matrix = _build_sparse_feature_matrix(contig_data)
        return _PreparedSimilarityMetric(self.weight, feature_matrix, self.__class__.__name__)


class AdjacencyIndex(ContigMetric):
    def __init__(self, weight: float = 0.0, databases=None):
        self.weight = weight
        self.databases = databases
        #self.contig_data = OrderedDict()
    
    def get_contig_data(self, recs, databases=None):
        """
            store the domain content for the contig
        """
        out = OrderedDict()
        for contig in recs:
            if contig.id not in out:
                out[contig.id] = set()
            
            contig.features.sort(key=lambda x: x.location.start) #sort features by start location
            last_feature = None
            for feature in contig.features:
                if feature.type == DOMAIN_FEATURE_NAME:
                    if databases is not None and feature.qualifiers["database"][0] not in databases:
                        continue
                    domain_name = feature.qualifiers["name"][0]
                    if last_feature is not None:
                        #out[contig.id].add(tuple(sorted([domain_name,last_feature])))  
                        out[contig.id].add(tuple([domain_name, last_feature])) # If we don't sort alphabetically, we are assuming that the contigs are oriented in the same way, which might be a bad assumption?
                    last_feature = domain_name
        return out
    
    def compute(self, recs, cpu=1, progress=False):
        contig_data = self.get_contig_data(recs, self.databases)
        return _compute_sparse_similarity_matrix(contig_data, cpu=cpu, progress=progress, desc=f"Scoring {self.__class__.__name__}")

    def prepare_sparse_metric(self, recs):
        contig_data = self.get_contig_data(recs, self.databases)
        _contig_ids, feature_matrix = _build_sparse_feature_matrix(contig_data)
        return _PreparedSimilarityMetric(self.weight, feature_matrix, self.__class__.__name__)
        
class DomainSequenceSimilarity(ContigMetric):
    pass


LINKAGE_METHODS = ("average", "ward", "complete", "single", "weighted")
LEAF_ORDERINGS = ("greedy", "optimal", "none")
_OPTIMAL_LEAF_ORDERING_WARN_CONTIGS = 5000
_CONDENSED_BLOCK_NNZ = 1 << 20


class SortMemoryLimitExceeded(RuntimeError):
    """Raised when sorting the contigs for genbank output would need more memory than is available."""


def estimate_sort_memory_bytes(n):
    """Peak memory for hierarchical sorting: the condensed distance vector plus the copy made by scipy's linkage."""
    return (n * (n - 1) // 2) * np.dtype(np.float64).itemsize * 2


def _condensed_index(n, i, j):
    """Index of the pair (i, j), with i < j, in a scipy condensed distance vector. Works elementwise on int64 arrays."""
    return n * i - i * (i + 1) // 2 + (j - i - 1)


def _sparse_scores_to_condensed_distances(scores, max_score):
    """
        Converts a sparse score matrix to a scipy condensed distance vector, distance = 1 - score/max_score,
        without building a dense n x n matrix. Pairs with no score have distance 1.
        If the matrix is not symmetric (e.g. after --max_hits_per_query), each pair uses the higher of its two scores.
    """
    n = scores.shape[0]
    scores = scipy.sparse.csr_array(scores, dtype=np.float64)
    scores.sum_duplicates()
    condensed = np.ones(n * (n - 1) // 2, dtype=np.float64)
    indptr = scores.indptr
    # Rows are processed in order, so a pair's upper-triangle entry (from row min(i, j)) is always written
    # before its lower-triangle entry (from row max(i, j)) is merged in.
    start = 0
    while start < n:
        end = int(np.searchsorted(indptr, indptr[start] + _CONDENSED_BLOCK_NNZ, side="right")) - 1
        end = min(max(end, start + 1), n)
        lo, hi = indptr[start], indptr[end]
        rows = np.repeat(np.arange(start, end, dtype=np.int64), np.diff(indptr[start:end + 1]))
        cols = scores.indices[lo:hi].astype(np.int64)
        distances = 1.0 - scores.data[lo:hi] / max_score

        upper = cols > rows
        condensed[_condensed_index(n, rows[upper], cols[upper])] = distances[upper]
        lower = cols < rows
        lower_idx = _condensed_index(n, cols[lower], rows[lower])
        condensed[lower_idx] = np.minimum(condensed[lower_idx], distances[lower])
        start = end
    np.clip(condensed, 0.0, 1.0, out=condensed)
    return condensed


def _greedy_leaf_order(linkage_matrix, condensed, n):
    """
        Orders the leaves of a dendrogram by flipping subtrees: at each merge, the two child leaf orders are oriented
        so that the leaves meeting at the junction are as close as possible (Gruvaeus & Wainer, 1972).
        Linear time, and on contig data about as good as scipy's optimal_leaf_ordering.
    """
    def dist(i, j):
        if i > j:
            i, j = j, i
        return condensed[_condensed_index(n, i, j)]

    clusters = {i: [i] for i in range(n)}
    for merge_idx, (left, right, _height, _size) in enumerate(linkage_matrix):
        a = clusters.pop(int(left))
        b = clusters.pop(int(right))
        # orientations: a+b, a+rev(b), rev(a)+b, rev(a)+rev(b). Ties keep the earliest, so results are deterministic.
        options = (dist(a[-1], b[0]), dist(a[-1], b[-1]), dist(a[0], b[0]), dist(a[0], b[-1]))
        best = options.index(min(options))
        if best in (2, 3):
            a.reverse()
        if best in (1, 3):
            b.reverse()
        a.extend(b)
        clusters[n + merge_idx] = a
    return clusters[2 * n - 2]


def _hierarchical_order(scores, linkage_method="average", leaf_ordering="greedy"):
    """
        Runs hierarchical clustering on a square score matrix (higher is more similar) and
        returns a list of record indices ordered so that similar records are adjacent.
    """
    if linkage_method not in LINKAGE_METHODS:
        raise ValueError(f"Unknown linkage method {linkage_method!r}, expected one of {LINKAGE_METHODS}")
    if leaf_ordering not in LEAF_ORDERINGS:
        raise ValueError(f"Unknown leaf ordering {leaf_ordering!r}, expected one of {LEAF_ORDERINGS}")
    n = scores.shape[0]
    if n <= 2:
        return list(range(n))
    scores = scipy.sparse.csr_array(scores, dtype=np.float64)
    if scores.nnz == 0:
        return list(range(n))
    min_score = scores.data.min()
    if min_score < 0:
        raise ValueError(f"compare_contigs scores must be non-negative to sort contigs, but the minimum score is {min_score}.")
    max_score = scores.data.max()
    if max_score <= 0: # if all scores are zero, just return the original order
        return list(range(n))

    condensed = _sparse_scores_to_condensed_distances(scores, max_score)
    if not condensed.any(): # if all distances are zero, just return the original order
        return list(range(n))
    linkage_matrix = scipy.cluster.hierarchy.linkage(condensed, method=linkage_method)
    if leaf_ordering == "optimal":
        linkage_matrix = scipy.cluster.hierarchy.optimal_leaf_ordering(linkage_matrix, condensed)
        return scipy.cluster.hierarchy.leaves_list(linkage_matrix).tolist()
    if leaf_ordering == "greedy":
        return _greedy_leaf_order(linkage_matrix, condensed, n)
    return scipy.cluster.hierarchy.leaves_list(linkage_matrix).tolist()


def _check_sort_feasible(n, leaf_ordering, check_memory=True):
    """Fails early if sorting n contigs would not fit in available memory, and warns if it will be very slow."""
    needed = estimate_sort_memory_bytes(n)
    available = psutil.virtual_memory().available if check_memory else None
    if check_memory and needed > available:
        raise SortMemoryLimitExceeded(
            f"Sorting {n} contigs for genbank output (-o) needs about {format_bytes(needed)} of memory, "
            f"but only {format_bytes(available)} is available. "
            "To avoid this, drop -o and write a matrix with --hdf instead, sort a subset of contigs with --contigs or --contigs_file, "
            "or deduplicate the input first. Use --skip_sort_memory_check to try anyway."
        )
    if leaf_ordering == "optimal" and n > _OPTIMAL_LEAF_ORDERING_WARN_CONTIGS:
        warnings.warn(
            f"--leaf_ordering optimal on {n} contigs may take many minutes to hours. "
            "Consider --leaf_ordering greedy, which is much faster and usually about as good."
        )


class DistanceReport(ABC):
    outpath: str

    @abstractmethod
    def write(self, recs, scores, options):
        """
            recs: a list of SeqRecords
            scores: a dok_sparse matrix of scores representing how similar each pair of contigs are
            options: a dict of additional options
        """
        pass

class HdfTableOutput(DistanceReport):
    def __init__(self, outpath, max_output_bytes=None):
        self.outpath = outpath
        self.max_output_bytes = max_output_bytes
    
    def write(self, recs, scores, options):
        row_names = [x.id for x in recs]
        write_matrix_outputs(scores, row_names, row_names, hdf=self.outpath, dense_text=None,
                             max_output_bytes=self.max_output_bytes, data_type="score",
                             mitigation_options=["--max_hits_per_query", "--mst", "--knn"])


class ClinkerOutput(DistanceReport):

    def __init__(self, outpath):
        self.outpath = outpath
    
    def write(self, recs, scores, options):
        warnings.warn("Clinker output not implemented")

class SvgTreeOutput(DistanceReport):

    def __init__(self, outpath):
        self.outpath = outpath
    
    def write(self, recs, scores, options):
        warnings.warn("SVG Tree output not implemented")

class GenbankOutput(DistanceReport):

    def __init__(self, outpath, linkage_method="average", leaf_ordering="greedy"):
        self.outpath = outpath
        self.linkage_method = linkage_method
        self.leaf_ordering = leaf_ordering

    def write(self, recs, scores, options):
        if len(recs) > 0:
            order = _hierarchical_order(scores, linkage_method=self.linkage_method, leaf_ordering=self.leaf_ordering)
            digits = math.ceil(math.log10(len(order)+1))
            with open(self.outpath, "w") as outfile:
                print_i = 0
                for i in order:
                    rec = recs[i]
                    if "name_by_order" in options:
                        rec = copy_SeqRecord(rec)
                        rec.id = str(print_i).zfill(digits) + "_" + rec.id
                        print_i += 1
                    write_genbank([rec], outfile)
        else:
            with open(self.outpath, "w") as outfile:
                pass


class DenseTextTableOutput(DistanceReport):

    def __init__(self, outpath, max_output_bytes=None):
        self.outpath = outpath
        self.max_output_bytes = max_output_bytes
    
    def write(self, recs, scores, options):
        row_names = [x.id for x in recs]
        write_matrix_outputs(scores, row_names, row_names, hdf=None, dense_text=self.outpath,
                             max_output_bytes=self.max_output_bytes,
                             mitigation_options=["--max_hits_per_query", "--mst", "--knn"])

def compare_contigs(genbanks, metrics, reports, max_hits_per_query, contigs, name_by_order, databases:Optional[Set[str]]=None, lb: float = 0.0, progress: bool = False, cpu: int = 1, mst: bool = False, knn: Optional[int] = None, check_sort_memory: bool = True):
    sparsify = sparsify_requested(mst, knn)

    report_options = dict()
    if name_by_order:
        report_options["name_by_order"] = True


    name_to_idx = dict()
    recs = list()
    #collect data needed for the metrics
    for rec in _progress(parse_seqfiles(genbanks, contigs), progress, desc="Loading contigs"):
        if rec.id not in name_to_idx:
            name_to_idx[rec.id] = len(recs)
            if databases is not None:
                rec = tuple(filter_domains((rec,), evalue=float("inf"), max_overlap=1, databases_keep=databases))[0]
            recs.append(rec)
        else:
            warnings.warn(f"Duplicate contig names !!! :{rec.id}")

    # fail before scoring if the genbank output can't be sorted
    for report in reports:
        if isinstance(report, GenbankOutput):
            _check_sort_feasible(len(recs), report.leaf_ordering, check_memory=check_sort_memory)

    matrix_shape = (len(recs), len(recs))
    if max_hits_per_query == 0:
        results_matrix = scipy.sparse.csr_array(matrix_shape, dtype=np.float64)
        for report in reports:
            report.write(recs, results_matrix, report_options)
        return

    prepared_metrics = list()
    for metric in metrics:
        prepared_metric = metric.prepare_sparse_metric(recs)
        if prepared_metric is None:
            prepared_metrics = None
            break
        prepared_metrics.append(prepared_metric)

    if prepared_metrics is not None:
        if sparsify:
            results_matrix = _compute_sparsified_similarity_matrix(
                prepared_metrics,
                mst=mst,
                knn=knn,
                k=max_hits_per_query,
                lb=lb,
                cpu=cpu,
                progress=progress,
                desc="Scoring combined metrics (sparsified)",
            )
        else:
            results_matrix = _compute_combined_sparse_similarity_matrix(
                prepared_metrics,
                k=max_hits_per_query,
                lb=lb,
                cpu=cpu,
                progress=progress,
                desc="Scoring combined metrics",
            )
    else:
        results_matrix = scipy.sparse.csr_array(matrix_shape, dtype=np.float64)
        remaining_weight = sum(metric.weight for metric in metrics)
        # fallback compute path for metrics that do not support combined chunk scoring
        for metric in metrics:
            metric_matrix = metric.compute(recs, cpu=cpu, progress=progress)
            results_matrix = results_matrix + (metric.weight * metric_matrix)
            remaining_weight -= metric.weight

            if lb > 0:
                results_matrix = _prune_scores_inplace(results_matrix, lb - remaining_weight)

        results_matrix = _keep_top_k_per_row(results_matrix, max_hits_per_query)
        if lb > 0:
            results_matrix = _prune_scores_inplace(results_matrix, lb)

        if sparsify:
            row_names = [rec.id for rec in recs]
            wrapped = DataMatrix(scipy.sparse.csr_array(results_matrix), row_names, row_names, data_type="score")
            results_matrix = sparsify_matrix(wrapped, mst=mst, knn=knn, lower_bound=lb)

    # write reports
    for report in reports:
        report.write(recs, results_matrix, report_options)


def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument('-i', '--input', nargs='+', default=None, required=True, help="Genbank filenames.")

    parser.add_argument('-o', '--output', default=None, required=False,
                        help="Name of genbank output, contigs will be sorted hierarchically within the output genbank. If not supplied then no genbank output will be generated. To write to stdout, use '-'. default: None")
    
    parser.add_argument("--name_by_order", action="store_true",
                        help="Will rename contigs in the genbank file such that they have a prefix that can be sorted to restore the sorted order.")

    parser.add_argument("--linkage", default="average", choices=LINKAGE_METHODS,
                        help="Hierarchical clustering linkage method used to sort contigs for the genbank output (-o). Distances are 1 - score/max(score); pairs with no score have distance 1.")

    parser.add_argument("--leaf_ordering", default="greedy", choices=LEAF_ORDERINGS,
                        help="How to order the leaves of the clustering tree for the genbank output (-o). "
                        "greedy: flip subtrees so that the most similar contigs meet at each junction; fast and usually about as good as optimal. "
                        "optimal: scipy's optimal leaf ordering; very slow (minutes to hours) beyond a few thousand contigs. "
                        "none: plain dendrogram order.")

    parser.add_argument("--skip_sort_memory_check", action="store_true",
                        help="Don't fail early when sorting the contigs for the genbank output (-o) is estimated to need more memory than is available.")

    parser.add_argument('--contigs', default=[], nargs='+',
                        help="only annotate contigs with ids in this list. Additive with --contigs_file. default: all contigs.")

    parser.add_argument('--contigs_file', default=None, 
                        help="only annotate contigs with ids listed in this file (one per line). Additive with --contigs. default: all contigs.")

    parser.add_argument('--max_hits_per_query', default=None, type=int,
                        help="After combining all metric scores, keep this many of the top hits for each query and consider all other hits to be zero. " + MAX_HITS_PER_QUERY_SPARSIFY_NOTE + " Default: keep all hits.")

    parser.add_argument('--lb', default=0, type=float, required=False,
                        help="Round any final scores less than or equal to this down to zero. This is primarily useful for making sparse outputs smaller.")

    add_sparsify_arguments(parser, requirements="Computed while streaming, so the full matrix is never held in memory. Best paired with --hdf.")

    # Metrics
    parser.add_argument("--ji", default=0.5, required=False, type=float,
                        help="weighting for jaccard index metric, a comparison of all distinct types of domains in the contigs. When two contigs have no domains, their ji is 0.0.")
    parser.add_argument("--ai", default=0.5, required=False, type=float,
                        help="weighting for adjacency index metric, a comparison of shared domain pairs in the contigs. When two contigs have 0 or 1 domains, their ai is 0.0.")
#   parser.add_argument("--dss", default=0.0, required=False,
#                        help="weighting for domain sequence similarity, a comparison of shared domain pairs in the contigs.")
    parser.add_argument("--databases", default=None, required=False, type=str, nargs="+",
                        help="Domain databases to use for domain content comparison. default: all databases.")

    # Output
    # parser.add_argument("--clinker", type=str, default=None, help="Write an interactive html in clinker format to this path.")

    # parser.add_argument("--svg_tree", type=str, default=None, help="Write a tree graphic in svg format to this path.")
    
    add_matrix_output_arguments(parser, noun="contig similarity score matrix")

    parser.add_argument('--cpu', type=int, default=0, required=False,
                        help="How many cpu workers to use. Default: use all available physical cores.")

    parser.add_argument('--progress', action='store_true',
                        help="Show a progress bar for long-running steps.")

    add_max_output_gb_argument(parser)

    parser.add_argument('--config', action=ActionConfigFile)


    params = parser.parse_args(argv)
   
    ### Figure out what input and output files ####

    # if pipe_in:
    #     genbanks = [sys.stdin]
    # elif params.genbank is not None:
    genbanks = params.input
    max_output_bytes = max_output_gb_to_bytes(params.max_output_gb)
    if params.cpu <= 0:
        cpus = psutil.cpu_count(logical=False)
    else:
        cpus = params.cpu


    reports = list()
    # if params.clinker is not None:
    #     reports.append(ClinkerOutput(params.clinker))
    # if params.svg_tree is not None:
    #     reports.append(SvgTreeOutput(params.svg_tree))
    hdf, dense_text = resolve_matrix_outputs(params, require_output=False)
    if hdf is not None:
        reports.append(HdfTableOutput(hdf, max_output_bytes=max_output_bytes))
    if dense_text is not None:
        reports.append(DenseTextTableOutput(dense_text, max_output_bytes=max_output_bytes))
    if params.output is not None:
        output_path = '/dev/stdout' if params.output == "-" else params.output # TODO: not cross platform, will not work on windows.
        reports.append(GenbankOutput(output_path, linkage_method=params.linkage, leaf_ordering=params.leaf_ordering))

    databases=None
    if params.databases is not None:
        databases = set(params.databases)

    if params.ji < 0 or params.ai < 0:
        raise ValueError("Metric weights (--ji, --ai) must be non-negative.")
    metrics = list()
    if params.ji != 0.0:
        metrics.append(JaccardIndex(params.ji, databases=databases))
    if params.ai != 0.0:
        metrics.append(AdjacencyIndex(params.ai, databases=databases))
    # if params.dss != 0.0:
    #     reports["dss"] = params.dss
    if len(metrics) == 0:
        raise ValueError("Please specify a non-zero weight for at least one metric.")


        
    ## figure out if we are using every contig or just certain, named contigs ##
    contigs_needed = list_and_file_to_dict_keys(params.contigs, params.contigs_file)
    
    # Run
    try:
        compare_contigs(genbanks, metrics, reports, params.max_hits_per_query, contigs=contigs_needed, name_by_order=params.name_by_order, lb=params.lb, progress=params.progress, cpu=cpus, mst=params.mst, knn=params.knn, check_sort_memory=not params.skip_sort_memory_check)
    except (OutputSizeLimitExceeded, SortMemoryLimitExceeded) as exc:
        raise SystemExit(str(exc)) from None

def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])
