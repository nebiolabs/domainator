"""Calculates similarity scores between sequence records using exact k-mer set comparisons.

Typically the input and reference will be the same sequence file, creating a pairwise similarity matrix.
Protein inputs use protein k-mers. Nucleotide inputs use translated CDS protein k-mers by default,
or strand-invariant nucleotide k-mers when --nucleotide_kmers is supplied.
"""

import warnings
warnings.filterwarnings("ignore", module='numpy')

import sys
from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List, Optional, Set, Tuple

import numpy as np
import psutil
import scipy.sparse
import tqdm
from jsonargparse import ActionConfigFile, ArgumentParser

from domainator import RawAndDefaultsFormatter, __version__
from domainator.compare_contigs import (
    _StreamingCSRChunkBuilder,
    _iter_query_chunks,
    _keep_top_k_per_row,
    _prune_scores_inplace,
)
from domainator.data_matrix import DataMatrix
from domainator.domainate import (
    clean_rec,
    get_prot_list,
    prodigal_CDS_annotate,
)
from domainator.output_guardrails import (
    OutputSizeLimitExceeded,
    add_max_output_gb_argument,
    enforce_matrix_output_limit,
    max_output_gb_to_bytes,
)
from domainator.transform_matrix import transform_matrix
from domainator.utils import get_file_type, get_multiprocessing_context, parse_seqfiles


DEFAULT_PROTEIN_KSIZE = 5
DEFAULT_NUCLEOTIDE_KSIZE = 15
VALID_GENE_CALL_VALUES = {None, "all", "unannotated"}
VALID_METRICS = {"max_containment", "min_containment", "jaccard"}
VALID_MODES = {"score", "bool", "score_dist"}
_PROTEIN_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
_PROTEIN_BASE = len(_PROTEIN_ALPHABET)
_PROTEIN_CHAR_TO_INT = {char: idx for idx, char in enumerate(_PROTEIN_ALPHABET)}
_NUCLEOTIDE_ALPHABET = "ACGTRYMKBDHVNSW"
_NUCLEOTIDE_COMPLEMENTS = "TGCAYRKMVHDBNSW"
_NUCLEOTIDE_BASE = len(_NUCLEOTIDE_ALPHABET)
_NUCLEOTIDE_CHAR_TO_INT = {char: idx for idx, char in enumerate(_NUCLEOTIDE_ALPHABET)}
_NUCLEOTIDE_COMPLEMENT_INT = {
    _NUCLEOTIDE_CHAR_TO_INT[char]: _NUCLEOTIDE_CHAR_TO_INT[complement]
    for char, complement in zip(_NUCLEOTIDE_ALPHABET, _NUCLEOTIDE_COMPLEMENTS)
}

_QUERY_FEATURE_MATRIX = None
_REFERENCE_FEATURE_MATRIX = None
_QUERY_ROW_COUNTS = None
_REFERENCE_ROW_COUNTS = None
_METRIC_NAME = None
_SELF_COMPARISON = None


@dataclass(frozen=True)
class SequenceRepresentation:
    names: List[str]
    lengths: List[int]
    row_indices: List[int]
    col_indices: List[int]


def _progress(iterable, enabled, **kwargs):
    if not enabled:
        return iterable
    return tqdm.tqdm(iterable, leave=True, dynamic_ncols=True, **kwargs)


def _progress_checkpoint(enabled: bool, message: str) -> None:
    if enabled:
        tqdm.tqdm.write(message)


def _clean_sequence_letters(sequence: str, *, nucleotide: bool) -> str:
    letters = ''.join(filter(str.isalnum, str(sequence))).upper()
    if nucleotide:
        letters = letters.replace("U", "T")
    return letters


def _iter_encoded_kmers(sequence: str, ksize: int, char_to_int: Dict[str, int], base: int) -> Iterator[int]:
    if ksize < 1:
        return iter(())

    clean_sequence = _clean_sequence_letters(sequence, nucleotide=False)
    if len(clean_sequence) < ksize:
        return iter(())

    base_power = base ** (ksize - 1)

    def iterator() -> Iterator[int]:
        current_code = 0
        valid_run = 0
        for char in clean_sequence:
            value = char_to_int.get(char)
            if value is None:
                current_code = 0
                valid_run = 0
                continue

            if valid_run < ksize:
                current_code = (current_code * base) + value
                valid_run += 1
                if valid_run == ksize:
                    yield current_code
                continue

            current_code = ((current_code % base_power) * base) + value
            yield current_code

    return iterator()


def _iter_protein_kmers(sequence: str, ksize: int) -> Iterator[int]:
    return _iter_encoded_kmers(sequence, ksize, _PROTEIN_CHAR_TO_INT, _PROTEIN_BASE)


def _iter_nucleotide_kmers(sequence: str, ksize: int) -> Iterator[int]:
    clean_sequence = _clean_sequence_letters(sequence, nucleotide=True)
    if len(clean_sequence) < ksize:
        return iter(())

    base_power = _NUCLEOTIDE_BASE ** (ksize - 1)

    def iterator() -> Iterator[int]:
        forward_code = 0
        reverse_code = 0
        valid_run = 0
        for char in clean_sequence:
            value = _NUCLEOTIDE_CHAR_TO_INT.get(char)
            if value is None:
                forward_code = 0
                reverse_code = 0
                valid_run = 0
                continue

            complement_value = _NUCLEOTIDE_COMPLEMENT_INT[value]
            if valid_run < ksize:
                forward_code = (forward_code * _NUCLEOTIDE_BASE) + value
                reverse_code = reverse_code + (complement_value * (_NUCLEOTIDE_BASE ** valid_run))
                valid_run += 1
                if valid_run == ksize:
                    yield forward_code if forward_code <= reverse_code else reverse_code
                continue

            forward_code = ((forward_code % base_power) * _NUCLEOTIDE_BASE) + value
            reverse_code = (reverse_code // _NUCLEOTIDE_BASE) + (complement_value * base_power)
            yield forward_code if forward_code <= reverse_code else reverse_code

    return iterator()


def _normalize_record_for_translated_kmers(record, gene_call: Optional[str]) -> None:
    if gene_call not in VALID_GENE_CALL_VALUES:
        raise ValueError("gene_call must be one of None, 'all', or 'unannotated'")

    if record.annotations["molecule_type"] == "protein":
        return

    clear_cds_annotations = gene_call == "all"
    cds_count = clean_rec(record, clear_CDS_annotations=clear_cds_annotations)
    if cds_count == 0 and gene_call is not None:
        prodigal_CDS_annotate(record)
        clean_rec(record)


def _record_feature_set(record, *, gene_call: Optional[str], nucleotide_kmers: bool, ksize: int) -> Set[int]:
    molecule_type = record.annotations["molecule_type"]
    if molecule_type == "protein":
        return set(_iter_protein_kmers(str(record.seq), ksize))

    if nucleotide_kmers:
        return set(_iter_nucleotide_kmers(str(record.seq), ksize))

    _normalize_record_for_translated_kmers(record, gene_call)
    feature_set = set()
    for _name, protein_sequence in get_prot_list(record, 0):
        feature_set.update(_iter_protein_kmers(protein_sequence, ksize))
    return feature_set


def _load_sequence_representation(
    path: str,
    file_type: Optional[str],
    *,
    feature_to_idx: Dict[int, int],
    fasta_type: str,
    gene_call: Optional[str],
    nucleotide_kmers: bool,
    ksize: int,
    progress: bool,
    label: str,
) -> SequenceRepresentation:
    if fasta_type == "nucleotide" and file_type == "fasta" and not nucleotide_kmers and gene_call is None:
        raise ValueError(
            "Nucleotide FASTA input requires either --gene_call to create CDS annotations "
            "or --nucleotide_kmers to compare whole-contig nucleotide k-mers."
        )

    _progress_checkpoint(progress, f"Reading {label} records from {path}")
    names = list()
    lengths = list()
    row_indices = list()
    col_indices = list()
    seen_names = set()
    found_record = False

    for record_idx, record in enumerate(
        _progress(
            parse_seqfiles([path], filetype_override=file_type, default_molecule_type=fasta_type),
            progress,
            desc=f"Parsing {label} k-mers",
            unit="records",
        )
    ):
        found_record = True
        if record.id in seen_names:
            raise ValueError(f"Duplicate record ids found in {path}")
        seen_names.add(record.id)
        names.append(record.id)
        lengths.append(len(record))

        feature_set = _record_feature_set(
            record,
            gene_call=gene_call,
            nucleotide_kmers=nucleotide_kmers,
            ksize=ksize,
        )
        for feature in sorted(feature_set):
            if feature not in feature_to_idx:
                feature_to_idx[feature] = len(feature_to_idx)
            row_indices.append(record_idx)
            col_indices.append(feature_to_idx[feature])

    if not found_record:
        raise ValueError(f"No sequence records found in {path}")

    return SequenceRepresentation(
        names=names,
        lengths=lengths,
        row_indices=row_indices,
        col_indices=col_indices,
    )


def _build_sparse_feature_matrix(
    representation: SequenceRepresentation,
    feature_count: int,
) -> scipy.sparse.csr_array:
    row_count = len(representation.names)
    if feature_count == 0:
        return scipy.sparse.csr_array((row_count, 0), dtype=np.float64)
    if len(representation.row_indices) == 0:
        return scipy.sparse.csr_array((row_count, feature_count), dtype=np.float64)

    data = np.ones(len(representation.row_indices), dtype=np.float64)
    return scipy.sparse.csr_array(
        (
            data,
            (
                np.asarray(representation.row_indices, dtype=np.int32),
                np.asarray(representation.col_indices, dtype=np.int32),
            ),
        ),
        shape=(row_count, feature_count),
    )


def _build_aligned_feature_matrices(
    query_representation: SequenceRepresentation,
    reference_representation: SequenceRepresentation,
    feature_count: int,
    progress: bool,
    self_comparison: bool,
) -> Tuple[scipy.sparse.csr_array, scipy.sparse.csr_array]:
    _progress_checkpoint(progress, f"Materializing query sparse feature matrix ({len(query_representation.names)} records, {feature_count} unique k-mers)")
    query_matrix = _build_sparse_feature_matrix(query_representation, feature_count)
    if self_comparison:
        _progress_checkpoint(progress, "Reusing query sparse feature matrix for reference")
        return query_matrix, query_matrix
    _progress_checkpoint(progress, f"Materializing reference sparse feature matrix ({len(reference_representation.names)} records, {feature_count} unique k-mers)")
    reference_matrix = _build_sparse_feature_matrix(reference_representation, feature_count)
    return query_matrix, reference_matrix


def _score_shared_counts(shared_counts, query_count: float, reference_counts, metric: str):
    if shared_counts.size == 0:
        return np.empty(0, dtype=np.float64)

    if metric == "jaccard":
        denominator = query_count + reference_counts - shared_counts
    elif metric == "max_containment":
        denominator = np.minimum(query_count, reference_counts)
    elif metric == "min_containment":
        denominator = np.maximum(query_count, reference_counts)
    else:
        raise ValueError(f"Unsupported metric: {metric}")

    scores = np.divide(
        shared_counts,
        denominator,
        out=np.zeros(shared_counts.shape, dtype=np.float64),
        where=denominator > 0,
    )
    return scores


def _compute_metric_chunk(
    query_feature_matrix: scipy.sparse.csr_array,
    reference_feature_matrix: scipy.sparse.csr_array,
    query_row_counts: np.ndarray,
    reference_row_counts: np.ndarray,
    start: int,
    end: int,
    metric: str,
    self_comparison: bool,
    k: Optional[int],
    lb: float,
) -> scipy.sparse.csr_array:
    chunk_shape = (end - start, reference_feature_matrix.shape[0])
    if chunk_shape[0] == 0:
        return scipy.sparse.csr_array(chunk_shape, dtype=np.float64)

    overlaps = (query_feature_matrix[start:end] @ reference_feature_matrix.transpose()).tocsr()
    overlaps.sort_indices()

    row_idx = list()
    col_idx = list()
    data = list()

    for local_query_idx in range(end - start):
        query_idx = start + local_query_idx
        query_count = query_row_counts[query_idx]

        if self_comparison:
            row_idx.append(local_query_idx)
            col_idx.append(query_idx)
            data.append(1.0)

        row_start = overlaps.indptr[local_query_idx]
        row_end = overlaps.indptr[local_query_idx + 1]
        target_indices = overlaps.indices[row_start:row_end]
        shared_counts = overlaps.data[row_start:row_end].astype(np.float64, copy=False)

        if self_comparison:
            off_diagonal_mask = target_indices != query_idx
            target_indices = target_indices[off_diagonal_mask]
            shared_counts = shared_counts[off_diagonal_mask]

        if target_indices.size == 0:
            continue

        scores = _score_shared_counts(
            shared_counts,
            float(query_count),
            reference_row_counts[target_indices].astype(np.float64, copy=False),
            metric,
        )
        positive_mask = scores > 0
        if not np.any(positive_mask):
            continue

        target_indices = target_indices[positive_mask]
        scores = scores[positive_mask]
        row_idx.extend([local_query_idx] * target_indices.size)
        col_idx.extend(target_indices.tolist())
        data.extend(scores.tolist())

    if len(data) == 0:
        chunk_matrix = scipy.sparse.csr_array(chunk_shape, dtype=np.float64)
    else:
        chunk_matrix = scipy.sparse.csr_array(
            (
                np.asarray(data, dtype=np.float64),
                (np.asarray(row_idx, dtype=np.int32), np.asarray(col_idx, dtype=np.int32)),
            ),
            shape=chunk_shape,
        )

    chunk_matrix = _keep_top_k_per_row(chunk_matrix, k)
    if lb > 0:
        chunk_matrix = _prune_scores_inplace(chunk_matrix, lb)
    return chunk_matrix


def _init_metric_chunk_worker(
    query_feature_matrix,
    reference_feature_matrix,
    query_row_counts,
    reference_row_counts,
    metric,
    self_comparison,
):
    global _QUERY_FEATURE_MATRIX, _REFERENCE_FEATURE_MATRIX, _QUERY_ROW_COUNTS, _REFERENCE_ROW_COUNTS, _METRIC_NAME, _SELF_COMPARISON
    _QUERY_FEATURE_MATRIX = query_feature_matrix
    _REFERENCE_FEATURE_MATRIX = reference_feature_matrix
    _QUERY_ROW_COUNTS = query_row_counts
    _REFERENCE_ROW_COUNTS = reference_row_counts
    _METRIC_NAME = metric
    _SELF_COMPARISON = self_comparison


def _compute_metric_chunk_worker(args):
    start, end, k, lb = args
    return (
        start,
        end,
        _compute_metric_chunk(
            _QUERY_FEATURE_MATRIX,
            _REFERENCE_FEATURE_MATRIX,
            _QUERY_ROW_COUNTS,
            _REFERENCE_ROW_COUNTS,
            start,
            end,
            _METRIC_NAME,
            _SELF_COMPARISON,
            k,
            lb,
        ),
    )


def _compute_sparse_metric_matrix(
    query_feature_matrix: scipy.sparse.csr_array,
    reference_feature_matrix: scipy.sparse.csr_array,
    metric: str,
    k: Optional[int],
    lb: float,
    cpu: int,
    progress: bool,
    self_comparison: bool,
) -> scipy.sparse.csr_array:
    matrix_shape = (query_feature_matrix.shape[0], reference_feature_matrix.shape[0])
    if k == 0:
        return scipy.sparse.csr_array(matrix_shape, dtype=np.float64)

    if query_feature_matrix.shape[1] == 0 or reference_feature_matrix.shape[1] == 0:
        if self_comparison:
            return scipy.sparse.identity(matrix_shape[0], dtype=np.float64, format="csr")
        return scipy.sparse.csr_array(matrix_shape, dtype=np.float64)

    query_row_counts = np.diff(query_feature_matrix.indptr)
    reference_row_counts = np.diff(reference_feature_matrix.indptr)
    query_chunks = [(start, end, k, lb) for start, end in _iter_query_chunks(matrix_shape[0], cpu)]
    chunk_builder = _StreamingCSRChunkBuilder(matrix_shape)
    _progress_checkpoint(progress, f"Scoring {matrix_shape[0]} query records against {matrix_shape[1]} reference records")

    if cpu > 1 and len(query_chunks) > 1:
        with get_multiprocessing_context().Pool(
            processes=min(cpu, len(query_chunks)),
            initializer=_init_metric_chunk_worker,
            initargs=(
                query_feature_matrix,
                reference_feature_matrix,
                query_row_counts,
                reference_row_counts,
                metric,
                self_comparison,
            ),
        ) as pool:
            for start, end, chunk_matrix in _progress(
                pool.imap(_compute_metric_chunk_worker, query_chunks),
                progress,
                total=len(query_chunks),
                desc="Scoring",
            ):
                chunk_builder.append_chunk(start, end, chunk_matrix)
    else:
        for start, end, chunk_k, chunk_lb in _progress(query_chunks, progress, total=len(query_chunks), desc="Scoring"):
            chunk_builder.append_chunk(
                start,
                end,
                _compute_metric_chunk(
                    query_feature_matrix,
                    reference_feature_matrix,
                    query_row_counts,
                    reference_row_counts,
                    start,
                    end,
                    metric,
                    self_comparison,
                    chunk_k,
                    chunk_lb,
                ),
            )

    return chunk_builder.build()


def _resolve_sequence_file_type(path: str, file_type: Optional[str], label: str) -> str:
    resolved_type = file_type if file_type is not None else get_file_type(path)
    if resolved_type not in {"fasta", "genbank"}:
        raise ValueError(f"{label} must be a fasta or genbank file.")
    return resolved_type


def _resolve_molecule_type(path: str, file_type: str, fasta_type: str) -> str:
    molecule_types = {
        record.annotations["molecule_type"]
        for record in parse_seqfiles([path], filetype_override=file_type, default_molecule_type=fasta_type, max_recs=1)
    }
    if len(molecule_types) != 1:
        raise ValueError(f"Mixed molecule types within a single input are not supported: {path}")
    return next(iter(molecule_types))


def _resolve_ksize(molecule_type: str, nucleotide_kmers: bool, ksize: Optional[int]) -> int:
    if ksize is not None:
        if ksize < 1:
            raise ValueError("--ksize must be at least 1")
        return ksize
    if molecule_type == "protein" or not nucleotide_kmers:
        return DEFAULT_PROTEIN_KSIZE
    return DEFAULT_NUCLEOTIDE_KSIZE


def kmer_dist(
    input_path: str,
    input_type: str,
    reference_path: str,
    reference_type: str,
    k: Optional[int],
    metric: str,
    mode: str,
    cpu: int,
    dense: Optional[str],
    dense_text: Optional[str],
    sparse: Optional[str],
    lb: float,
    *,
    fasta_type: str,
    ksize: Optional[int],
    nucleotide_kmers: bool,
    gene_call: Optional[str],
    max_output_bytes: Optional[int],
    progress: bool,
) -> None:
    input_molecule_type = _resolve_molecule_type(input_path, input_type, fasta_type)
    reference_molecule_type = _resolve_molecule_type(reference_path, reference_type, fasta_type)
    if input_molecule_type != reference_molecule_type:
        raise ValueError("Input and reference must have the same molecule type.")

    self_comparison = input_path == reference_path and input_type == reference_type

    resolved_ksize = _resolve_ksize(input_molecule_type, nucleotide_kmers, ksize)
    feature_to_idx: Dict[int, int] = dict()
    query_representation = _load_sequence_representation(
        input_path,
        input_type,
        feature_to_idx=feature_to_idx,
        fasta_type=fasta_type,
        gene_call=gene_call,
        nucleotide_kmers=nucleotide_kmers,
        ksize=resolved_ksize,
        progress=progress,
        label="query",
    )
    if self_comparison:
        reference_representation = query_representation
    else:
        reference_representation = _load_sequence_representation(
            reference_path,
            reference_type,
            feature_to_idx=feature_to_idx,
            fasta_type=fasta_type,
            gene_call=gene_call,
            nucleotide_kmers=nucleotide_kmers,
            ksize=resolved_ksize,
            progress=progress,
            label="reference",
        )

    query_feature_matrix, reference_feature_matrix = _build_aligned_feature_matrices(
        query_representation,
        reference_representation,
        len(feature_to_idx),
        progress,
        self_comparison,
    )
    matrix = _compute_sparse_metric_matrix(
        query_feature_matrix,
        reference_feature_matrix,
        metric,
        k,
        lb,
        cpu,
        progress,
        self_comparison=self_comparison,
    )
    matrix = transform_matrix(
        matrix,
        mode,
        row_lengths=query_representation.lengths,
        col_lengths=reference_representation.lengths,
    )

    try:
        if dense:
            enforce_matrix_output_limit(
                output_type="dense",
                matrix=matrix,
                row_names=query_representation.names,
                col_names=reference_representation.names,
                row_lengths=query_representation.lengths,
                col_lengths=reference_representation.lengths,
                data_type=mode,
                max_output_bytes=max_output_bytes,
                output_path=dense,
                mitigation_options=["--sparse", "-k", "--lb"],
            )
            DataMatrix.write_dense(
                matrix,
                dense,
                query_representation.names,
                reference_representation.names,
                query_representation.lengths,
                reference_representation.lengths,
                mode,
            )
        if dense_text:
            enforce_matrix_output_limit(
                output_type="dense_text",
                matrix=matrix,
                row_names=query_representation.names,
                col_names=reference_representation.names,
                row_lengths=query_representation.lengths,
                col_lengths=reference_representation.lengths,
                data_type=mode,
                max_output_bytes=max_output_bytes,
                output_path=dense_text,
                mitigation_options=["--sparse", "-k", "--lb"],
            )
            DataMatrix.write_dense_text(
                matrix,
                dense_text,
                query_representation.names,
                reference_representation.names,
                query_representation.lengths,
                reference_representation.lengths,
                mode,
            )
        if sparse:
            enforce_matrix_output_limit(
                output_type="sparse",
                matrix=matrix,
                row_names=query_representation.names,
                col_names=reference_representation.names,
                row_lengths=query_representation.lengths,
                col_lengths=reference_representation.lengths,
                data_type=mode,
                max_output_bytes=max_output_bytes,
                output_path=sparse,
                mitigation_options=["-k", "--lb"],
            )
            DataMatrix.write_sparse(
                matrix,
                sparse,
                query_representation.names,
                reference_representation.names,
                query_representation.lengths,
                reference_representation.lengths,
                mode,
            )
    except OutputSizeLimitExceeded as exc:
        raise SystemExit(str(exc)) from None


def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)
    parser.add_argument('-i', '--input', type=str, required=True,
                        help="Input sequence file. This can be a fasta or genbank file.")
    parser.add_argument('--input_type', type=str, required=False, default=None, choices={"fasta", "genbank"},
                        help="File type for the input sequence file.")
    parser.add_argument('-r', '--reference', type=str, required=False,
                        help="Reference file. If not provided, the input file will be used as the reference.")
    parser.add_argument('--reference_type', type=str, required=False, default=None, choices={"fasta", "genbank"},
                        help="File type for the reference sequence file.")
    parser.add_argument('--fasta_type', type=str, default="protein", choices={"protein", "nucleotide"},
                        help="Whether fasta inputs should be interpreted as protein or nucleotide sequences.")
    parser.add_argument('--gene_call', type=str, default=None, choices={"all", "unannotated"}, required=False,
                        help="When activated, new CDS annotations will be added with Prodigal in Metagenomic mode before translated-CDS k-mer extraction. "
                        "If 'all', then any existing CDS annotations will be deleted and all contigs will be re-annotated. "
                        "If 'unannotated', then only contigs without CDS annotations will be annotated.")
    parser.add_argument('--nucleotide_kmers', action='store_true',
                        help="Use strand-invariant nucleotide k-mers instead of translated CDS protein k-mers for nucleotide inputs.")
    parser.add_argument('--ksize', type=int, default=None,
                        help="K-mer size. Defaults to 5 for protein/translated-CDS comparisons and 15 for nucleotide-kmer comparisons.")
    parser.add_argument('--metric', type=str, default="max_containment", choices=VALID_METRICS,
                        help="Similarity metric to compute from shared k-mer counts.")
    parser.add_argument('--mode', type=str, required=False, default="score", choices=VALID_MODES,
                        help="Output matrix values. score: raw similarity, bool: 1 if a similarity is non-zero, score_dist: 1 - (score / min(row_max, col_max)).")
    parser.add_argument('--dense', type=str, default=None, help="Write a dense distance matrix hdf5 file to this path.")
    parser.add_argument('--dense_text', type=str, default=None, help="Write a dense distance matrix tsv file to this path.")
    parser.add_argument('--sparse', type=str, default=None, help="Write a sparse distance matrix hdf5 file to this path.")
    parser.add_argument('--lb', default=0, type=float, required=False,
                        help="Round any similarities lower than or equal to this down to zero before mode transforms.")
    parser.add_argument('-k', type=int, required=False, default=None,
                        help="Include at most this many non-zero entries in the matrix for each input sequence.")
    parser.add_argument('--cpu', type=int, default=0, required=False,
                        help="How many cpu threads to use. Default: use all physical cores.")
    parser.add_argument('--progress', action='store_true',
                        help="Show a progress bar for long-running steps.")
    add_max_output_gb_argument(parser)
    parser.add_argument('--config', action=ActionConfigFile)

    params = parser.parse_args(argv)

    input_type = _resolve_sequence_file_type(params.input, params.input_type, "Input")
    if params.reference is None:
        reference_path = params.input
        reference_type = input_type
    else:
        reference_path = params.reference
        reference_type = _resolve_sequence_file_type(reference_path, params.reference_type, "Reference")

    if params.dense is not None and get_file_type(params.dense) != "hdf5":
        raise ValueError("Please use an hdf5 related extension for the --dense output, such as .h5, .hdf5, or .hdf.")
    if params.sparse is not None and get_file_type(params.sparse) != "hdf5":
        raise ValueError("Please use an hdf5 related extension for the --sparse output, such as .h5, .hdf5, or .hdf.")
    if params.mode == "score_dist" and params.sparse is not None:
        raise ValueError("Sparse distance matrices not implemented.")
    if not 0 <= params.lb <= 1:
        raise ValueError("--lb must be between 0 and 1 for k-mer similarities.")
    if params.k is not None and params.k < 0:
        raise ValueError("-k must be >= 0")
    if params.gene_call is not None and params.nucleotide_kmers:
        warnings.warn("--gene_call has no effect when --nucleotide_kmers is supplied.")

    dense = params.dense
    dense_text = params.dense_text
    sparse = params.sparse
    if dense is None and dense_text is None and sparse is None:
        raise ValueError("No output specified! Please specify at least one of: dense, dense_text, sparse")

    if params.cpu <= 0:
        cpu = psutil.cpu_count(logical=False)
    else:
        cpu = params.cpu

    max_output_bytes = max_output_gb_to_bytes(params.max_output_gb)
    kmer_dist(
        params.input,
        input_type,
        reference_path,
        reference_type,
        params.k,
        params.metric,
        params.mode,
        cpu,
        dense,
        dense_text,
        sparse,
        params.lb,
        fasta_type=params.fasta_type,
        ksize=params.ksize,
        nucleotide_kmers=params.nucleotide_kmers,
        gene_call=params.gene_call,
        max_output_bytes=max_output_bytes,
        progress=params.progress,
    )


def _entrypoint():
    main(sys.argv[1:])


if __name__ == '__main__':
    main(sys.argv[1:])