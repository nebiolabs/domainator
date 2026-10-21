"""Calculates similarity or distance matrices between protein structures.

The structure-input counterpart of domainator seq_dist. Typically only -i is supplied, giving an
all-vs-all matrix of the input structures; a separate -r produces a rectangular matrix of
input structures against reference structures.

Every chain of every input is a row, whether or not it had a hit, so the row and column
labels are exactly the record ids that domainator structure_to_genbank writes for the same
structures, and the matrix feeds domainator build_ssn, domainator build_tree or domainator build_ssn_viewer
directly. That does mean the input database is read end to end, unlike domainator structure_search.

Alignment direction is the opposite of domainator structure_search and domainator structure_domainate: the
input is the QUERY and the reference is the TARGET, so rows are inputs and columns are
references, and --evalue is computed against the reference set. This matches domainator seq_dist,
where diamond is run with the input as -q and the reference as -d.

--mode selects what goes in the cells. The score-based modes (score, bool, norm_score,
row_norm_score, score_dist, efi_score, efi_score_dist) are computed from the aligner's
alignment score and behave exactly as they do in domainator seq_dist. The structural modes put a
similarity that is already bounded in [0, 1] in the cell -- tmscore, lddt, fident -- and
their *_dist forms are simply 1 - that value, needing none of the min(row_max, col_max)
normalization that score_dist does. tmscore and lddt require --algorithm foldseek and
C-alpha coordinates in both databases; fident works with either backend.

Note that with --algorithm reseek there is no bit score: structure_lib records -log10(p)
in its place. That is monotonic in significance, so score, bool and the normalized modes
remain meaningful, but the values are not comparable to foldseek's and efi_score, which
assumes a bit score in its 2^(-score) term, is rejected.
"""

import logging
import warnings
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import scipy.sparse
import tqdm
from jsonargparse import ActionConfigFile
from jsonargparse.typing import PositiveFloat

from domainator.cli.common import add_input, tool_parser
from domainator import structure_lib
from domainator.sparse_scores import keep_top_k_per_row
from domainator.ssn_edges import StreamingMstKnnAccumulator
from domainator.data_matrix import DataMatrix, symmetrize_matrix
from domainator.matrix_output import (
    add_matrix_output_arguments,
    add_symmetrize_argument,
    resolve_matrix_outputs,
    write_matrix_outputs,
)
from domainator.output_guardrails import (
    add_max_output_gb_argument,
    max_output_gb_to_bytes,
)
from domainator.seq_dist import _SparseMaxResultBuilder
from domainator.sparsify import (
    MAX_HITS_PER_QUERY_SPARSIFY_NOTE,
    PerQueryTopK,
    add_sparsify_arguments,
    sparsify_option_label,
    sparsify_requested,
)
from domainator.transform_matrix import (
    LOG10_2,
    MODES,
    transform_matrix,
)

logger = logging.getLogger(__name__)


# The seven score-based modes, taken from transform_matrix so the two tools cannot drift.
BITSCORE_MODES = frozenset(MODES)

# Structural modes: the StructureHit field each one reads.
STRUCTURE_VALUE_FIELDS = {"tmscore": "tmscore", "lddt": "lddt", "fident": "fident"}

# Distance forms of the structural modes. These are 1 - value rather than a normalized
# transform, because the underlying similarity is already bounded in [0, 1].
STRUCTURE_DIST_MODES = {"tmscore_dist": "tmscore", "lddt_dist": "lddt",
                        "fident_dist": "fident"}

ALL_MODES = BITSCORE_MODES | set(STRUCTURE_VALUE_FIELDS) | set(STRUCTURE_DIST_MODES)

# Per-mode structural metrics to request from the backend. The score-based and fident
# modes need nothing extra: fident is in foldseek's base output columns, and reseek
# derives it from pctid.
MODE_METRICS = {"tmscore": ("tmscore",), "tmscore_dist": ("tmscore",),
                "lddt": ("lddt",), "lddt_dist": ("lddt",)}

# Modes whose cell value is computable from a single hit, so --mst/--knn can stream
# them. The rest normalize against global matrix maxima, which streaming cannot see.
SPARSIFY_STREAMABLE_MODES = frozenset({"score", "bool", "efi_score",
                                       "tmscore", "lddt", "fident"})

# efi_score reads the cell as a bit score in its 2^(-score) term, which -log10(p) is not.
BITSCORE_ONLY_MODES = frozenset({"efi_score", "efi_score_dist"})

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
    already the extension-free '<basename>_<chain>' ids that domainator structure_to_genbank
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


def _apply_mode(matrix, mode: str, row_lengths, col_lengths, self_comparison: bool):
    """Turn the raw value matrix into the requested mode."""
    if mode in BITSCORE_MODES:
        return transform_matrix(matrix, mode=mode, row_lengths=row_lengths,
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
    """The value --mst/--knn select on for one hit, since that path sees no maxima.

    bool selects on the raw value, since all-1 values would leave nothing to rank; the
    caller turns the kept cells into a 0/1 mask afterward. Takes the already-read raw
    value rather than the hit, so the reader counts each alignment exactly once.
    """
    if mode == "efi_score":
        value = raw * LOG10_2 - (log_lengths[row_idx] + log_lengths[col_idx])
        return value if value > 0 else None
    return raw


def structure_dist(input_values: Sequence[str],
                   *,
                   references: Optional[Sequence[str]],
                   aligner: structure_lib.StructureAligner,
                   mode: str = "score",
                   evalue: float = 0.001,
                   max_hits_per_query: Optional[int] = None,
                   threshold: float = 0.0,
                   symmetrize: str = "none",
                   mst: bool = False,
                   knn: Optional[int] = None,
                   alignment_type: Optional[int] = None,
                   max_seqs: Optional[int] = None,
                   tmp_dir: Optional[str] = None,
                   db_out: Optional[str] = None,
                   hits_tsv_out: Optional[str] = None,
                   progress: bool = False) -> DataMatrix:
    """Align every input structure against every reference structure.

    Returns a DataMatrix of ``mode`` values with a row per input and a column per reference.
    references of None compares the input against itself. ``db_out`` and ``hits_tsv_out``
    keep the aligner database and the raw hits table.
    """
    sparsify = sparsify_requested(mst, knn)

    self_comparison = references is None
    metrics = _mode_metrics(mode)
    field = _mode_value_field(mode)

    with structure_lib.prepared_databases(aligner, input_values, references,
                                          tmp_dir=tmp_dir, keep_db=db_out) as prepared:
        aligner.check_capabilities(alignment_type, metrics,
                                   [prepared.input_db, prepared.reference_db])

        row_name_to_idx, row_names, row_lengths = _enumerate_entries(
            aligner, prepared.input_db, "input")
        if self_comparison:
            col_name_to_idx, col_names, col_lengths = row_name_to_idx, row_names, row_lengths
        else:
            col_name_to_idx, col_names, col_lengths = _enumerate_entries(
                aligner, prepared.reference_db, "reference")

        # foldseek's own default caps hits at 1000 per query, which would silently
        # truncate any all-vs-all larger than that. The reference count is the smallest
        # cap that cannot truncate anything.
        effective_max_seqs = max_seqs if max_seqs is not None else len(col_names)
        if max_seqs is None and effective_max_seqs > SATURATION_PREFLIGHT_THRESHOLD:
            logger.info(
                f"Comparing {len(row_names)} structures against {len(col_names)}, with no "
                "--max_seqs cap so that no row is truncated. If this is too slow, bound "
                "the search with --evalue, or cap it with --max_seqs and accept that "
                "rows reaching the cap are incomplete.")

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

        if sparsify:
            log_lengths = (np.log10(np.asarray(row_lengths, dtype=np.float64))
                           if mode == "efi_score" else None)
            accumulator = StreamingMstKnnAccumulator(len(row_names), mst=mst, knn=knn, lower_bound=0)

            def add_raw(row_idx, col_idx, raw):
                value = _streaming_value(raw, mode, log_lengths, row_idx, col_idx)
                if value is not None:
                    accumulator.add_edge(row_idx, col_idx, value)

            # The per-query cap ranks raw values, as keep_top_k_per_row does on the
            # non-streaming path, and has to see a query's every hit before it can choose.
            top_hits = (PerQueryTopK(len(row_names), max_hits_per_query)
                        if max_hits_per_query is not None else None)
            for hit in _progress(hits, progress, desc="Scoring pairs"):
                hit_counts[hit.query] = hit_counts.get(hit.query, 0) + 1
                raw = reader.value(hit)
                if raw is None or raw <= threshold:
                    continue
                row_idx = row_name_to_idx[hit.query]
                col_idx = col_name_to_idx[hit.target]
                if top_hits is not None:
                    top_hits.add_edge(row_idx, col_idx, raw)
                else:
                    add_raw(row_idx, col_idx, raw)
            if top_hits is not None:
                for row_idx, col_idx, raw in top_hits.iter_edges():
                    add_raw(row_idx, col_idx, raw)
            matrix = accumulator.to_csr()
            if mode == "bool":
                matrix.data[:] = 1
        else:
            # The order-agnostic builder: it makes no assumption about how the backend
            # groups its output, and duplicate (query, target) pairs -- which are what the
            # grouped builder exists to collapse early -- are rare in a structural search.
            builder = _SparseMaxResultBuilder(row_name_to_idx, col_name_to_idx)
            for hit in _progress(hits, progress, desc="Scoring pairs"):
                hit_counts[hit.query] = hit_counts.get(hit.query, 0) + 1
                value = reader.value(hit)
                if value is None or value <= threshold:
                    continue
                builder.add_result(hit.query, hit.target, value)
            matrix = builder.build()

        reader.warn()
        # Only a cap the user actually asked for can be "saturated"; the derived default
        # is exactly the column count, which a fully-connected row reaches legitimately.
        structure_lib.warn_on_saturation(hit_counts, max_seqs, subject="input record")
        if hits_tsv_out is not None:
            # The work directory goes away when this context exits.
            structure_lib.copy_hits_tsv(prepared.work_dir, hits_tsv_out, aligner)

    if not sparsify:
        # The streaming accumulator has already reconciled both directions itself, and
        # applied the per-query cap before pruning.
        matrix = symmetrize_matrix(matrix, symmetrize)
        matrix = keep_top_k_per_row(matrix, max_hits_per_query)  # a no-op when None

    row_lengths_array = np.array(row_lengths)
    col_lengths_array = (row_lengths_array if self_comparison else np.array(col_lengths))

    if not sparsify:
        matrix = _apply_mode(matrix, mode, row_lengths_array, col_lengths_array,
                             self_comparison)

    return DataMatrix(matrix, row_names, col_names, row_lengths_array, col_lengths_array, data_type=mode, copy_data=False)


def build_parser():
    parser = tool_parser(__doc__)

    add_input(parser, "structures to compare: pdb/cif files (optionally gzipped), a directory searched recursively, a glob, or the prefix of a prebuilt aligner database.", stdin=False)
    parser.add_argument('-r', '--reference', default=None, nargs='+', type=str,
                        help="structures to compare the input against. Same accepted forms as --input. If not supplied, the input is compared against itself, giving a square all-vs-all matrix. Rows are inputs and columns are references.")

    add_matrix_output_arguments(parser, noun="similarity (or, for *_dist modes, distance) matrix")

    parser.add_argument('-e', '--evalue', type=PositiveFloat, default=0.001,
                        help="only pairs with an E-value strictly below this become non-zero cells. E-values are computed against the reference set, as they are in seq_dist. Raise it to fill in weaker relationships; lower it to shrink the matrix.")

    parser.add_argument('--mode', type=str, default="score", choices=sorted(ALL_MODES),
                        help="what kind of values should be in the matrix. Score-based, from the aligner's alignment score, exactly as in seq_dist: score: raw score, bool: 1 if a hit otherwise 0, norm_score: score/min(row_max, col_max), row_norm_score: score/row_max, score_dist: 1 - norm_score, efi_score: -log10[2^(-score) * (row_len * col_len)], efi_score_dist: 1 - (efi_score / min(row_max, col_max)). Structural, already in [0,1] so the distance forms are just 1 - value: tmscore, lddt, fident, tmscore_dist, lddt_dist, fident_dist. tmscore and lddt require --algorithm foldseek and C-alpha coordinates in both databases; fident works with either backend. efi_score and efi_score_dist require a real bit score, so they are rejected with --algorithm reseek.")

    add_symmetrize_argument(parser, extra_help="Requires a square self-comparison (no -r). Cannot be combined with --mst/--knn, whose output is always max-symmetric.")

    parser.add_argument('--threshold', type=float, default=0,
                        help="Round any values lower than or equal to this down to zero. This applies to the raw per-hit value before the --mode transform: the aligner's alignment score for the score-based modes, and the metric itself (0-1) for the tmscore/lddt/fident modes.")

    parser.add_argument('--max_hits_per_query', type=int, default=None,
                        help="Include at most this many non-zero entries in the matrix for each input structure, keeping the highest-valued ones. Applied after --symmetrize. " + MAX_HITS_PER_QUERY_SPARSIFY_NOTE + " Use --max_seqs to bound the search itself.")

    add_sparsify_arguments(parser, requirements="Computed while streaming, so the full matrix is never held in memory. Requires a square self-comparison and --mode in {score, bool, efi_score, tmscore, lddt, fident}; bool selects on the alignment scores and then keeps the selected cells as 1. Best paired with --hdf_out.")

    parser.add_argument('--db_out', default=None, type=str,
                        help="when the input is structure files rather than a prebuilt database, also write the database built from them, using this path as its prefix.")
    parser.add_argument('--hits_tsv_out', default=None, type=str,
                        help="write the aligner's own hit table here, for debugging or for cross-checking against a hand-run search.")
    parser.add_argument('--progress', action='store_true',
                        help="Show a progress bar for long-running steps.")

    structure_lib.add_backend_arguments(parser)
    add_max_output_gb_argument(parser)
    parser.add_argument('--config', action=ActionConfigFile)
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)

    if params.metrics:
        raise ValueError(
            "structure_dist derives the structural metric from --mode (tmscore, lddt, "
            "fident), so --metrics is rejected here rather than looking as though it "
            "changes the matrix."
        )

    if params.max_hits_per_query is not None and params.max_hits_per_query < 0:
        raise ValueError("--max_hits_per_query must be >= 0")
    hdf, dense_text = resolve_matrix_outputs(params)
    if params.mode in STRUCTURE_VALUE_FIELDS or params.mode in STRUCTURE_DIST_MODES:
        if not 0 <= params.threshold <= 1:
            raise ValueError(f"--threshold must be between 0 and 1 for --mode {params.mode}.")

    # None means "compare the input against itself", which also lets prepared_databases
    # build one database instead of two.
    reference_values = params.reference
    if reference_values is not None and list(reference_values) == list(params.input):
        reference_values = None
    self_comparison = reference_values is None

    if sparsify_requested(params.mst, params.knn):
        option = sparsify_option_label(params.mst, params.knn)
        if not self_comparison:
            raise ValueError(f"{option} requires a square, symmetric comparison: omit -r, "
                             "or pass the same value for -i and -r.")
        if params.mode not in SPARSIFY_STREAMABLE_MODES:
            raise ValueError(
                f"{option} is only supported with --mode {sorted(SPARSIFY_STREAMABLE_MODES)} "
                f"(modes that normalize against global matrix maxima cannot be streamed). "
                f"For '{params.mode}', run 'structure_dist --mode score --hdf_out out.hdf5' "
                f"then 'transform_matrix -i out.hdf5 --mode {params.mode} {option} "
                f"--hdf_out pruned.hdf5'."
            )
        if params.symmetrize != "none":
            raise ValueError(
                f"{option} output is always max-symmetric -- the accumulator selects "
                "edges on max(M[i,j], M[j,i]) and writes that value to both cells -- so "
                f"--symmetrize {params.symmetrize} cannot be applied on top of it. Drop "
                f"--symmetrize, or run 'structure_dist --mode score --hdf_out out.hdf5' "
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

    matrix = structure_dist(
        params.input,
        references=reference_values,
        aligner=aligner,
        mode=params.mode,
        evalue=params.evalue,
        max_hits_per_query=params.max_hits_per_query,
        threshold=params.threshold,
        symmetrize=params.symmetrize,
        mst=params.mst,
        knn=params.knn,
        alignment_type=params.alignment_type,
        max_seqs=params.max_seqs,
        tmp_dir=params.tmp_dir,
        db_out=params.db_out,
        hits_tsv_out=params.hits_tsv_out,
        progress=params.progress,
    )
    write_matrix_outputs(matrix, hdf=hdf, dense_text=dense_text, max_output_bytes=max_output_bytes,
                         mitigation_options=["--max_hits_per_query", "--threshold", "-e/--evalue", "--mst", "--knn"])
