"""Calculates similarity or distance matrices between protein structures.

The structure-input counterpart of seq_dist.py. Typically only -i is supplied, giving an
all-vs-all matrix of the input structures; a separate -r produces a rectangular matrix of
input structures against reference structures.

Every chain of every input is a row, whether or not it had a hit, so the row and column
labels are exactly the record ids that structure_to_genbank.py writes for the same
structures, and the matrix feeds build_ssn.py, build_tree.py or build_ssn_viewer.py
directly. That does mean the input database is read end to end, unlike structure_search.py.

Alignment direction is the opposite of structure_search.py and structure_domainate.py: the
input is the QUERY and the reference is the TARGET, so rows are inputs and columns are
references, and --evalue is computed against the reference set. This matches seq_dist.py,
where diamond is run with the input as -q and the reference as -d.

--mode selects what goes in the cells. The score-based modes (score, bool, norm_score,
row_norm_score, score_dist, efi_score, efi_score_dist) are computed from the aligner's
alignment score and behave exactly as they do in seq_dist.py. The structural modes put a
similarity that is already bounded in [0, 1] in the cell -- tmscore, lddt, fident -- and
their *_dist forms are simply 1 - that value, needing none of the min(row_max, col_max)
normalization that score_dist does. tmscore and lddt require --algorithm foldseek and
C-alpha coordinates in both databases; fident works with either backend.

Note that with --algorithm reseek there is no bit score: structure_lib records -log10(p)
in its place. That is monotonic in significance, so score, bool and the normalized modes
remain meaningful, but the values are not comparable to foldseek's and efi_score, which
assumes a bit score in its 2^(-score) term, is rejected.
"""

import sys
import warnings
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import scipy.sparse
import tqdm
from jsonargparse import ActionConfigFile, ArgumentParser

from domainator import RawAndDefaultsFormatter, __version__, structure_lib
from domainator.compare_contigs import _keep_top_k_per_row
from domainator.data_matrix import DataMatrix, StreamingMstKnnAccumulator
from domainator.output_guardrails import (
    OutputSizeLimitExceeded,
    add_max_output_gb_argument,
    max_output_gb_to_bytes,
    enforce_matrix_output_limit,
)
from domainator.seq_dist import _SparseMaxResultBuilder
from domainator.transform_matrix import (
    LOG10_2,
    MODES,
    _knn_arg,
    _mst_knn_arg,
    transform_matrix,
)
from domainator.utils import get_file_type


# The seven score-based modes, taken from transform_matrix so the two tools cannot drift.
BITSCORE_MODES = frozenset(MODES)

# Structural modes: the StructureHit field each one reads.
STRUCTURE_VALUE_FIELDS = {"tmscore": "tmscore", "lddt": "lddt", "fident": "fident"}

# Distance forms of the structural modes. These are 1 - value rather than a normalized
# transform, because the underlying similarity is already bounded in [0, 1].
STRUCTURE_DIST_MODES = {"tmscore_dist": "tmscore", "lddt_dist": "lddt",
                        "fident_dist": "fident"}

ALL_MODES = BITSCORE_MODES | set(STRUCTURE_VALUE_FIELDS) | set(STRUCTURE_DIST_MODES)

# Modes whose cells are a distance, where an absent pair means "maximally distant" rather
# than zero. Sparse output would store that as an implicit zero and invert its meaning.
DIST_MODES = frozenset({"score_dist", "efi_score_dist"} | set(STRUCTURE_DIST_MODES))

# Per-mode structural metrics to request from the backend. The score-based and fident
# modes need nothing extra: fident is in foldseek's base output columns, and reseek
# derives it from pctid.
MODE_METRICS = {"tmscore": ("tmscore",), "tmscore_dist": ("tmscore",),
                "lddt": ("lddt",), "lddt_dist": ("lddt",)}

# Modes whose cell value is computable from a single hit, so --mst_knn/--knn can stream
# them. The rest normalize against global matrix maxima, which streaming cannot see.
MST_KNN_STREAMABLE_MODES = frozenset({"score", "bool", "efi_score",
                                      "tmscore", "lddt", "fident"})

# efi_score reads the cell as a bit score in its 2^(-score) term, which -log10(p) is not.
BITSCORE_ONLY_MODES = frozenset({"efi_score", "efi_score_dist"})

SYMMETRIZE_CHOICES = ("none", "min", "max", "mean")

# Above this many reference entries, say so before starting an uncapped all-vs-all.
SATURATION_PREFLIGHT_THRESHOLD = 10_000


def _progress(iterable, enabled, **kwargs):
    if not enabled:
        return iterable
    return tqdm.tqdm(iterable, leave=True, dynamic_ncols=True, **kwargs)


def _mode_value_field(mode: str) -> str:
    """The StructureHit attribute whose value fills the matrix for this mode."""
    if mode in STRUCTURE_VALUE_FIELDS:
        return STRUCTURE_VALUE_FIELDS[mode]
    if mode in STRUCTURE_DIST_MODES:
        return STRUCTURE_DIST_MODES[mode]
    return "bits"


def _mode_metrics(mode: str) -> List[str]:
    return list(MODE_METRICS.get(mode, ()))


def reject_mode_on_backend(aligner: structure_lib.StructureAligner, mode: str) -> None:
    """Fail on a --mode the selected backend cannot produce, naming --mode not --metrics.

    check_capabilities says the same thing in terms of metrics, which is confusing when
    the user never passed --metrics; this runs first so the error names what they typed.
    """
    if mode in BITSCORE_ONLY_MODES and aligner.name == "reseek":
        raise RuntimeError(
            f"--mode {mode} is computed from an alignment bit score, which backend "
            "'reseek' does not report (structure_lib records -log10(p-value) in its "
            "place, so the 2^(-score) term would be meaningless). Use --algorithm "
            "foldseek, or a --mode that only needs the values to be ordered, such as "
            "score, bool, norm_score, row_norm_score or score_dist."
        )
    missing = sorted(set(_mode_metrics(mode)) - set(aligner.supports_metrics))
    if missing:
        supported = sorted(aligner.supports_metrics)
        raise RuntimeError(
            f"--mode {mode} is computed from the {'/'.join(missing)} metric, which "
            f"backend '{aligner.name}' cannot produce (supported metrics: "
            f"{supported if supported else 'none'}). Use --algorithm foldseek, or a "
            "--mode based on the alignment score (score, bool, norm_score, "
            "row_norm_score, score_dist) or on sequence identity (fident, fident_dist), "
            "all of which every backend supports."
        )


def _enumerate_entries(aligner: structure_lib.StructureAligner,
                       db: structure_lib.StructureDB,
                       label: str) -> Tuple[Dict[str, int], List[str], List[int]]:
    """Every entry in a database, as (name -> index, names, lengths).

    Read from the database rather than from the hits, because a structure with no hits at
    all still needs its row and column, and efi_score needs a length for every one of
    them. Names are used verbatim: foldseek is built with --chain-name-mode 1, so they are
    already the extension-free '<basename>_<chain>' ids that structure_to_genbank.py
    writes, and stripping them again would corrupt a name like 'some.long.name_A'.
    """
    name_to_idx: Dict[str, int] = dict()
    names: List[str] = []
    lengths: List[int] = []
    for name, sequence, _source in aligner.iter_sequences(db):
        if name in name_to_idx:
            warnings.warn(
                f"Duplicate entry name '{name}' in the {label} database; keeping the "
                "first and skipping the rest.",
                RuntimeWarning,
            )
            continue
        name_to_idx[name] = len(names)
        names.append(name)
        lengths.append(len(sequence))
    if not names:
        raise ValueError(f"The {label} database contains no structures.")
    return name_to_idx, names, lengths


# foldseek's alntmscore normalization overshoots slightly on a perfect self-alignment --
# 1.01 to 1.05 is routine -- so clamping into [0, 1] is expected housekeeping, not a
# problem worth a warning on every run. Only a value well outside the range suggests the
# column was misread, and those are worth reporting.
METRIC_CLAMP_TOLERANCE = 0.1


class _MetricValueReader:
    """Reads one StructureHit field, reporting bad values once rather than per hit."""

    def __init__(self, field: str):
        self.field = field
        self.bounded = field in set(STRUCTURE_VALUE_FIELDS.values())
        self.missing = 0
        self.far_out_of_range = 0

    def value(self, hit) -> Optional[float]:
        raw = getattr(hit, self.field)
        if raw is None:
            # foldseek writes a non-float when it cannot compute a metric for a pair.
            self.missing += 1
            return None
        value = float(raw)
        if self.bounded and not 0.0 <= value <= 1.0:
            if not -METRIC_CLAMP_TOLERANCE <= value <= 1.0 + METRIC_CLAMP_TOLERANCE:
                self.far_out_of_range += 1
            value = min(1.0, max(0.0, value))
        return value

    def warn(self) -> None:
        if self.missing:
            warnings.warn(
                f"The aligner reported no {self.field} for {self.missing} alignment(s); "
                "those pairs are left as zero in the matrix.",
                RuntimeWarning,
            )
        if self.far_out_of_range:
            warnings.warn(
                f"{self.far_out_of_range} {self.field} value(s) fell well outside [0, 1] "
                f"(more than {METRIC_CLAMP_TOLERANCE} beyond it) and were clamped.",
                RuntimeWarning,
            )


def symmetrize_matrix(matrix, how: str):
    """Reconcile the two directions of every pair.

    Applied to the raw values, before the --mode transform, because it reconciles two
    measurements of one pair rather than adjusting a normalization -- and because 1 - x
    flips the ordering, so 'min' of a distance would be 'max' of a similarity.

    A pair the aligner reported in only one direction has a stored zero in the other,
    which means "no alignment reported", not "measured as zero". min therefore deletes
    such a pair and mean halves it; only max treats a one-directional hit as evidence.
    """
    if how == "none":
        return matrix
    sparse = scipy.sparse.issparse(matrix)
    transposed = matrix.T
    if how == "max":
        out = matrix.maximum(transposed) if sparse else np.maximum(matrix, transposed)
    elif how == "min":
        out = matrix.minimum(transposed) if sparse else np.minimum(matrix, transposed)
    elif how == "mean":
        out = (matrix + transposed) / 2.0
    else:
        raise ValueError(f"Unknown --symmetrize value: {how}")
    if scipy.sparse.issparse(out):
        out = scipy.sparse.csr_array(out)
        out.eliminate_zeros()
    return out


def _apply_mode(matrix, mode: str, row_lengths, col_lengths, self_comparison: bool):
    """Turn the raw value matrix into the requested mode."""
    if mode in BITSCORE_MODES:
        return transform_matrix(matrix, mode, row_lengths=row_lengths,
                                col_lengths=col_lengths)
    if mode in STRUCTURE_VALUE_FIELDS:
        return matrix  # the cell already holds the metric
    # A structural distance: 1 - similarity, necessarily dense because an absent pair is
    # the maximum distance rather than zero.
    dense = matrix.toarray() if scipy.sparse.issparse(matrix) else np.asarray(matrix)
    out = 1.0 - dense
    if self_comparison:
        # A structure is at distance 0 from itself by definition, even if its self-hit
        # fell below --evalue and never reached the matrix.
        np.fill_diagonal(out, 0.0)
    return out


def _streaming_value(raw: float, mode: str, log_lengths,
                     row_idx: int, col_idx: int) -> Optional[float]:
    """The cell value for one hit, for the --mst_knn/--knn path, which sees no maxima.

    Takes the already-read raw value rather than the hit, so the reader counts each
    alignment exactly once.
    """
    if mode == "bool":
        return 1.0
    if mode == "efi_score":
        value = raw * LOG10_2 - log_lengths[row_idx] - log_lengths[col_idx]
        return value if value > 0 else None
    return raw


def _write_matrix_outputs(matrix, mode, row_names, col_names, row_lengths, col_lengths,
                          dense, dense_text, sparse, max_output_bytes):
    """The enforce-then-write pair for each requested output format."""
    targets = (
        ("dense", dense, DataMatrix.write_dense, ["--sparse", "-k", "--lb", "-e/--evalue"]),
        ("dense_text", dense_text, DataMatrix.write_dense_text,
         ["--sparse", "-k", "--lb", "-e/--evalue"]),
        ("sparse", sparse, DataMatrix.write_sparse, ["-k", "--lb", "-e/--evalue"]),
    )
    try:
        for output_type, path, writer, mitigations in targets:
            if path is None:
                continue
            enforce_matrix_output_limit(
                output_type=output_type,
                matrix=matrix,
                row_names=row_names,
                col_names=col_names,
                row_lengths=row_lengths,
                col_lengths=col_lengths,
                data_type=mode,
                max_output_bytes=max_output_bytes,
                output_path=path,
                mitigation_options=mitigations,
            )
            writer(matrix, path, row_names, col_names, row_lengths, col_lengths, mode)
    except OutputSizeLimitExceeded as exc:
        raise SystemExit(str(exc)) from None


def structure_dist(input_values: Sequence[str],
                   reference_values: Optional[Sequence[str]],
                   aligner: structure_lib.StructureAligner,
                   mode: str = "score",
                   *,
                   dense: Optional[str] = None,
                   dense_text: Optional[str] = None,
                   sparse: Optional[str] = None,
                   evalue: float = 0.001,
                   k: Optional[int] = None,
                   lb: float = 0.0,
                   symmetrize: str = "none",
                   mst_knn: Optional[int] = None,
                   knn: Optional[int] = None,
                   alignment_type: Optional[int] = None,
                   max_seqs: Optional[int] = None,
                   tmp_dir: Optional[str] = None,
                   keep_db: Optional[str] = None,
                   hits_tsv: Optional[str] = None,
                   max_output_bytes: Optional[int] = None,
                   progress: bool = False) -> None:
    """Align every input structure against every reference structure, and write a matrix.

    reference_values of None compares the input against itself.
    """
    if mst_knn is not None and knn is not None:
        raise ValueError("--mst_knn already includes the kNN edges; pass either "
                         "--mst_knn or --knn, not both.")
    sparsify_k = mst_knn if mst_knn is not None else knn
    sparsify_include_mst = mst_knn is not None

    self_comparison = reference_values is None
    metrics = _mode_metrics(mode)
    field = _mode_value_field(mode)

    with structure_lib.prepared_databases(aligner, input_values, reference_values,
                                          tmp_dir=tmp_dir, keep_db=keep_db) as prepared:
        aligner.check_capabilities(alignment_type, metrics,
                                   [prepared.input_db, prepared.reference_db])

        row_name_to_idx, row_names, row_lengths = _enumerate_entries(
            aligner, prepared.input_db, "input")
        if self_comparison:
            # The same list objects, not equal copies: DataMatrix tests identity to decide
            # whether to record SYMMETRIC_LABELS, which build_ssn and build_ssn_viewer
            # then require.
            col_name_to_idx, col_names, col_lengths = row_name_to_idx, row_names, row_lengths
        else:
            col_name_to_idx, col_names, col_lengths = _enumerate_entries(
                aligner, prepared.reference_db, "reference")

        # foldseek's own default caps hits at 1000 per query, which would silently
        # truncate any all-vs-all larger than that. The reference count is the smallest
        # cap that cannot truncate anything.
        effective_max_seqs = max_seqs if max_seqs is not None else len(col_names)
        if max_seqs is None and effective_max_seqs > SATURATION_PREFLIGHT_THRESHOLD:
            print(
                f"Comparing {len(row_names)} structures against {len(col_names)}, with no "
                "--max_seqs cap so that no row is truncated. If this is too slow, bound "
                "the search with --evalue, or cap it with --max_seqs and accept that "
                "rows reaching the cap are incomplete.",
                file=sys.stderr, flush=True,
            )

        hits = aligner.search(
            prepared.input_db, prepared.reference_db,
            evalue=evalue,
            max_seqs=effective_max_seqs,
            want_tseq=False,
            alignment_type=alignment_type,
            metrics=metrics,
            work_dir=prepared.work_dir,
        )

        reader = _MetricValueReader(field)
        hit_counts: Dict[str, int] = dict()

        if sparsify_k is not None:
            log_lengths = (np.log10(np.asarray(row_lengths, dtype=np.float64))
                           if mode == "efi_score" else None)
            accumulator = StreamingMstKnnAccumulator(
                len(row_names), sparsify_k, lower_bound=0,
                include_mst=sparsify_include_mst)
            for hit in _progress(hits, progress, desc="Scoring pairs"):
                hit_counts[hit.query] = hit_counts.get(hit.query, 0) + 1
                raw = reader.value(hit)
                if raw is None or raw <= lb:
                    continue
                row_idx = row_name_to_idx[hit.query]
                col_idx = col_name_to_idx[hit.target]
                value = _streaming_value(raw, mode, log_lengths, row_idx, col_idx)
                if value is None:
                    continue
                accumulator.add_edge(row_idx, col_idx, value)
            matrix = accumulator.to_csr()
        else:
            # The order-agnostic builder: it makes no assumption about how the backend
            # groups its output, and duplicate (query, target) pairs -- which are what the
            # grouped builder exists to collapse early -- are rare in a structural search.
            builder = _SparseMaxResultBuilder(row_name_to_idx, col_name_to_idx)
            for hit in _progress(hits, progress, desc="Scoring pairs"):
                hit_counts[hit.query] = hit_counts.get(hit.query, 0) + 1
                value = reader.value(hit)
                if value is None or value <= lb:
                    continue
                builder.add_result(hit.query, hit.target, value)
            matrix = builder.build()

        reader.warn()
        # Only a cap the user actually asked for can be "saturated"; the derived default
        # is exactly the column count, which a fully-connected row reaches legitimately.
        structure_lib.warn_on_saturation(hit_counts, max_seqs, subject="input record")
        if hits_tsv is not None:
            # The work directory goes away when this context exits.
            structure_lib.copy_hits_tsv(prepared.work_dir, hits_tsv, aligner)

    if sparsify_k is None:
        # The streaming accumulator has already reconciled both directions itself.
        matrix = symmetrize_matrix(matrix, symmetrize)
    matrix = _keep_top_k_per_row(matrix, k)  # a no-op when k is None

    row_lengths_array = np.array(row_lengths)
    col_lengths_array = (row_lengths_array if self_comparison else np.array(col_lengths))

    if sparsify_k is None:
        matrix = _apply_mode(matrix, mode, row_lengths_array, col_lengths_array,
                             self_comparison)

    _write_matrix_outputs(matrix, mode, row_names, col_names,
                          row_lengths_array, col_lengths_array,
                          dense, dense_text, sparse, max_output_bytes)


def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__,
                            formatter_class=RawAndDefaultsFormatter)

    parser.add_argument('-i', '--input', required=True, nargs='+', type=str,
                        help="structures to compare: pdb/cif files (optionally gzipped), a directory searched recursively, a glob, or the prefix of a prebuilt aligner database.")
    parser.add_argument('-r', '--reference', default=None, nargs='+', type=str,
                        help="structures to compare the input against. Same accepted forms as --input. If not supplied, the input is compared against itself, giving a square all-vs-all matrix. Rows are inputs and columns are references.")

    parser.add_argument("--dense", type=str, default=None,
                        help="Write a dense distance matrix hdf5 file to this path.")
    parser.add_argument("--dense_text", type=str, default=None,
                        help="Write a dense distance matrix tsv file to this path.")
    parser.add_argument("--sparse", type=str, default=None,
                        help="Write a sparse distance matrix hdf5 file to this path. Not available for the *_dist modes, whose absent cells are the maximum distance rather than zero.")

    parser.add_argument('-e', '--evalue', type=float, default=0.001,
                        help="only pairs with an E-value strictly below this become non-zero cells. E-values are computed against the reference set, as they are in seq_dist. Raise it to fill in weaker relationships; lower it to shrink the matrix.")

    parser.add_argument('--mode', type=str, default="score", choices=sorted(ALL_MODES),
                        help="what kind of values should be in the matrix. Score-based, from the aligner's alignment score, exactly as in seq_dist: score: raw score, bool: 1 if a hit otherwise 0, norm_score: score/min(row_max, col_max), row_norm_score: score/row_max, score_dist: 1 - norm_score, efi_score: -log10[2^(-score) * (row_len * col_len)], efi_score_dist: 1 - (efi_score / min(row_max, col_max)). Structural, already in [0,1] so the distance forms are just 1 - value: tmscore, lddt, fident, tmscore_dist, lddt_dist, fident_dist. tmscore and lddt require --algorithm foldseek and C-alpha coordinates in both databases; fident works with either backend. efi_score and efi_score_dist require a real bit score, so they are rejected with --algorithm reseek. Default: score")

    parser.add_argument('--symmetrize', type=str, default="none", choices=list(SYMMETRIZE_CHOICES),
                        help="reconcile the two directions of each pair before the --mode transform. Requires a square self-comparison (no -r). none: leave the matrix as the aligner reported it, which is usually slightly asymmetric. max: the better of the two directions in both cells, the OR-symmetric convention used by --mst_knn and by build_ssn. min: the worse of the two -- note a pair found in only one direction has a stored zero in the other, so min DELETES it, making this a mutual-hit filter. mean: the average of the two -- note a pair found in only one direction is HALVED, because the missing direction contributes zero. Cannot be combined with --mst_knn/--knn, whose output is always max-symmetric.")

    parser.add_argument('--lb', type=float, default=0,
                        help="Round any values lower than or equal to this down to zero. This applies to the raw per-hit value before the --mode transform: the aligner's alignment score for the score-based modes, and the metric itself (0-1) for the tmscore/lddt/fident modes.")

    parser.add_argument('-k', type=int, default=None,
                        help="Include at most this many non-zero entries in the matrix for each input structure, keeping the highest-valued ones. Applied to the assembled matrix, after --symmetrize. Use --max_seqs to bound the search itself.")

    sparsify_group = parser.add_mutually_exclusive_group(required=False)
    sparsify_group.add_argument('--mst_knn', type=_mst_knn_arg, default=None,
                                help="Prune the output to the maximum spanning tree plus OR-symmetric k-nearest-neighbor edges (integer >= 0, where 0 keeps only the maximum spanning tree), computed as a streaming operation to keep memory and output size small. Requires a square self-comparison and a --mode whose value is computable from a single alignment. Best paired with --sparse.")
    sparsify_group.add_argument('--knn', type=_knn_arg, default=None,
                                help="Prune the output to OR-symmetric k-nearest-neighbor edges only (integer >= 1), computed as a streaming operation. Unlike --mst_knn this does not preserve the connected components of the full comparison. Same requirements as --mst_knn.")

    parser.add_argument('--keep_db', default=None, type=str,
                        help="when the input is structure files rather than a prebuilt database, also write the database built from them, using this path as its prefix.")
    parser.add_argument('--hits_tsv', default=None, type=str,
                        help="write the aligner's own hit table here, for debugging or for cross-checking against a hand-run search.")
    parser.add_argument('--progress', action='store_true',
                        help="Show a progress bar for long-running steps.")

    structure_lib.add_backend_arguments(parser)
    add_max_output_gb_argument(parser)
    parser.add_argument('--config', action=ActionConfigFile)

    params = parser.parse_args(argv)

    if params.metrics:
        raise ValueError(
            "structure_dist derives the structural metric from --mode (tmscore, lddt, "
            "fident), so --metrics is rejected here rather than looking as though it "
            "changes the matrix."
        )

    if params.evalue <= 0:
        raise ValueError("--evalue must be greater than 0.")
    if params.k is not None and params.k < 0:
        raise ValueError("-k must be >= 0")
    if params.dense is not None and get_file_type(params.dense) != "hdf5":
        raise ValueError("Please use an hdf5 related extension for the --dense output, such as .h5, .hdf5, or .hdf.")
    if params.sparse is not None and get_file_type(params.sparse) != "hdf5":
        raise ValueError("Please use an hdf5 related extension for the --sparse output, such as .h5, .hdf5, or .hdf.")
    if params.dense is None and params.dense_text is None and params.sparse is None:
        raise ValueError("No output specified! Please specify at least one of: dense, dense_text, sparse")
    if params.sparse is not None and params.mode in DIST_MODES:
        raise ValueError(
            f"Sparse distance matrices not implemented. With --mode {params.mode} an "
            "absent pair is the maximum distance, not zero, so a sparse file would "
            "invert the meaning of every cell it omits. Use --dense or --dense_text."
        )
    if params.mode in STRUCTURE_VALUE_FIELDS or params.mode in STRUCTURE_DIST_MODES:
        if not 0 <= params.lb <= 1:
            raise ValueError(f"--lb must be between 0 and 1 for --mode {params.mode}.")

    # None means "compare the input against itself", which also lets prepared_databases
    # build one database instead of two.
    reference_values = params.reference
    if reference_values is not None and list(reference_values) == list(params.input):
        reference_values = None
    self_comparison = reference_values is None

    sparsify_k = params.mst_knn if params.mst_knn is not None else params.knn
    if sparsify_k is not None:
        option = "--mst_knn" if params.mst_knn is not None else "--knn"
        if not self_comparison:
            raise ValueError(f"{option} requires a square, symmetric comparison: omit -r, "
                             "or pass the same value for -i and -r.")
        if params.mode not in MST_KNN_STREAMABLE_MODES:
            raise ValueError(
                f"{option} is only supported with --mode {sorted(MST_KNN_STREAMABLE_MODES)} "
                f"(modes that normalize against global matrix maxima cannot be streamed). "
                f"For '{params.mode}', run 'structure_dist --mode score --sparse out.hdf5' "
                f"then 'transform_matrix -i out.hdf5 --mode {params.mode} {option} "
                f"{sparsify_k} --sparse pruned.hdf5'."
            )
        if params.symmetrize != "none":
            raise ValueError(
                f"{option} output is always max-symmetric -- the accumulator selects "
                "edges on max(M[i,j], M[j,i]) and writes that value to both cells -- so "
                f"--symmetrize {params.symmetrize} cannot be applied on top of it. Drop "
                f"--symmetrize, or run 'structure_dist --mode score --sparse out.hdf5' "
                f"and prune with transform_matrix afterwards."
            )
    elif params.symmetrize != "none" and not self_comparison:
        raise ValueError("--symmetrize requires a square, symmetric comparison: omit -r, "
                         "or pass the same value for -i and -r.")

    aligner = structure_lib.build_aligner(params)
    reject_mode_on_backend(aligner, params.mode)
    aligner.effective_max_seqs(params.max_seqs)
    aligner.check_capabilities(params.alignment_type, _mode_metrics(params.mode), [])

    max_output_bytes = max_output_gb_to_bytes(params.max_output_gb)

    structure_dist(
        params.input,
        reference_values,
        aligner,
        params.mode,
        dense=params.dense,
        dense_text=params.dense_text,
        sparse=params.sparse,
        evalue=params.evalue,
        k=params.k,
        lb=params.lb,
        symmetrize=params.symmetrize,
        mst_knn=params.mst_knn,
        knn=params.knn,
        alignment_type=params.alignment_type,
        max_seqs=params.max_seqs,
        tmp_dir=params.tmp_dir,
        keep_db=params.keep_db,
        hits_tsv=params.hits_tsv,
        max_output_bytes=max_output_bytes,
        progress=params.progress,
    )


def _entrypoint():
    main(sys.argv[1:])


if __name__ == '__main__':
    main(sys.argv[1:])
