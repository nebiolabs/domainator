"""Transforms a score matrix into various kinds of normalized matrices

   For matrices with any data type, can convert between hdf5 and dense text, or change the datatype.
"""

import warnings
warnings.filterwarnings("ignore", module='numpy')
from jsonargparse import ActionConfigFile

from domainator.cli.common import add_input, tool_parser
from domainator.data_matrix import DataMatrix, is_distance_data_type, symmetrize_matrix
from domainator.matrix_output import add_matrix_output_arguments, add_symmetrize_argument, resolve_matrix_outputs, write_matrix_outputs
from domainator.output_guardrails import add_max_output_gb_argument, max_output_gb_to_bytes, OutputSizeLimitExceeded
import scipy.sparse
import numpy as np
from domainator.sparsify import add_sparsify_arguments, check_sparsify_data_type, sparsify_matrix, sparsify_option_label, sparsify_requested

PASS_THROUGH_MODE = "none"
LOG10_2 = np.log10(2.0)

def _axis_maxes(array):
    """Return (row_max, col_max) as 1D float arrays; absent sparse cells count as 0."""
    if scipy.sparse.issparse(array):
        array = scipy.sparse.csr_array(array)
        return (array.max(axis=1).toarray().ravel().astype(np.float64),
                array.max(axis=0).toarray().ravel().astype(np.float64))
    array = np.asarray(array)
    return array.max(axis=1).astype(np.float64), array.max(axis=0).astype(np.float64)


def _scale_stored_values(array, row_divisor, col_divisor=None):
    """Divide each stored value by ``row_divisor[row]`` (or the min of row and column divisors).

    Returns a CSR array; cells that are absent in the input stay absent.
    """
    array = scipy.sparse.coo_array(array)
    divisor = row_divisor[array.row]
    if col_divisor is not None:
        divisor = np.minimum(divisor, col_divisor[array.col])
    return scipy.sparse.csr_array((array.data / divisor, (array.row, array.col)), shape=array.shape)


def norm_score(array, row_lengths=None, col_lengths=None):
    """
        score / min(row_max, col_max)
    """
    row_max, col_max = _axis_maxes(array)
    return _scale_stored_values(array, row_max, col_max)

def row_norm_score(array, row_lengths=None, col_lengths=None):
    """
        score / row_max
    """
    row_max, _ = _axis_maxes(array)
    return _scale_stored_values(array, row_max)

def efi_score(array, row_lengths=None, col_lengths=None):
    """
        -log10[2^(-score) * (row_seq_length * col_seq_length)]

        scores less than 0 are set to 0.

        where row_length and col_length are the lenghts of the sequences or profiles.

        This is the score used by EFI-EST.

    """
    if row_lengths is None or col_lengths is None:
        raise ValueError("Row and column lengths must be provided for EFI score.")

    row_logs = np.log10(np.asarray(row_lengths, dtype=np.float64))
    col_logs = np.log10(np.asarray(col_lengths, dtype=np.float64))

    # A cell with no score would get a negative EFI score, which is clamped to 0, so only
    # the stored cells need computing.
    array = scipy.sparse.coo_array(array)
    # Sum the two length terms first: addition commutes exactly, so a symmetric score
    # matrix gives an exactly symmetric result (and is stored as a triangle).
    scores = array.data * LOG10_2 - (row_logs[array.row] + col_logs[array.col])
    keep_mask = scores > 0
    return scipy.sparse.csr_array(
        (scores[keep_mask], (array.row[keep_mask], array.col[keep_mask])),
        shape=array.shape,
    )

def efi_score_dist(array, row_lengths=None, col_lengths=None):
    """
        1 - (efi_score / min(row_max, col_max))
    """
    return score_dist(efi_score(array, row_lengths=row_lengths, col_lengths=col_lengths))


def score_dist(array, row_lengths=None, col_lengths=None):
    """
    1 - (score / min(row_max, col_max))

    Returns a dense numpy array: a distance matrix has a value in (nearly) every cell.
    """
    #TODO: in https://doi.org/10.1093/nar/gkaa635, they use -ln( (score / min(row_max, col_max)) ), instead of 1 - ...
    #       unsure whether we should use the linear or log version.
    row_max, col_max = _axis_maxes(array)
    array = array.toarray() if scipy.sparse.issparse(array) else np.asarray(array)
    return 1 - (array / np.minimum(row_max[:, np.newaxis], col_max[np.newaxis, :]))


MODES = { 
        "score": None, 
        "bool": None, 
        "norm_score": norm_score,
        "row_norm_score": row_norm_score,
        "score_dist": score_dist,
        "efi_score": efi_score,
        "efi_score_dist": efi_score_dist}


def apply_lower_bound(array, lower_bound):
    if lower_bound is None:
        return array

    if scipy.sparse.issparse(array):
        out = scipy.sparse.csr_array(array, copy=True)
        out.data[out.data <= lower_bound] = 0
        out.eliminate_zeros()
        return out

    out = np.array(array, copy=True)
    out[out <= lower_bound] = 0
    return out


def transform_matrix(array, mode, row_lengths=None, col_lengths=None):
    """
    Transform an array to a different mode.
    """
    if mode == "score":
        return array
    elif mode == "bool":
        return (array > 0).astype(np.int8)
    else:
        return MODES[mode](array, row_lengths=row_lengths, col_lengths=col_lengths)


def build_parser():
    parser = tool_parser(__doc__)

    add_input(parser, "name of input matrix file.", multiple=False, stdin=False)

    parser.add_argument("--value_type", choices=sorted(MODES.keys()), default=None, help="What the values of the input matrix are, if the file does not say. Overrides the type stored in the file.")
    
    add_matrix_output_arguments(parser, noun="transformed matrix")

    parser.add_argument("--mode", type=str, required=False, default=PASS_THROUGH_MODE, choices=set(MODES.keys()).union({PASS_THROUGH_MODE}), 
                        help="what kind of values should be in the output matrix. none: preserve the input values and data type, score: raw score, bool: 1 if a hit otherwise 0, score_dist: 1 - (score / min(row_max, col_max)), norm_score: score/min(row_max, col_max), efi_score: -log10[2^(-score) * (input_seq_length * reference_seq_length)], efi_score_dist: 1 - (efi_score / min(row_max, col_max)). Default: none")

    parser.add_argument('--threshold', type=float, default=None, required=False,
                        help="Zero out all values less than or equal to this threshold after any mode transformation.")
    add_sparsify_arguments(parser, requirements="Uses the post-transform values, except that --mode bool selects on the scores and then keeps the selected cells as 1; requires a square, non-distance matrix.")

    add_symmetrize_argument(parser, extra_help="Applied before --mode; requires a square matrix with the same row and column labels, holding similarity (not *_dist) values. Cannot be combined with --mst/--knn, whose output is always max-symmetric.")

    add_max_output_gb_argument(parser)
    
    parser.add_argument('--config', action=ActionConfigFile)

    # TODO: option to supply row and/or column lengths if not supplied in the input matrix.

    # TODO: option to supply row and/or column names to override those in the input matrix.
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)



    hdf, dense_text = resolve_matrix_outputs(params)

    max_output_bytes = max_output_gb_to_bytes(params.max_output_gb)

    matrix = DataMatrix.from_file(params.input)

    if params.value_type is not None:
        if matrix.data_type == "" or matrix.data_type == params.value_type:
            pass
        else:
            warnings.warn(f"Input matrix format is '{matrix.data_type}', but --value_type is '{params.value_type}'. The input matrix format is overridden.")
        matrix.data_type = params.value_type

    requested_mode = None if params.mode == PASS_THROUGH_MODE else params.mode
    sparsify = sparsify_requested(params.mst, params.knn)

    if params.symmetrize != "none":
        if sparsify:
            raise ValueError(f"{sparsify_option_label(params.mst, params.knn)} output is always max-symmetric, "
                             f"so --symmetrize {params.symmetrize} cannot be applied on top of it.")
        if not matrix.symmetric_labels:
            raise ValueError("--symmetrize requires a square matrix with the same row and column labels.")
        if is_distance_data_type(matrix.data_type):
            raise ValueError(f"--symmetrize reconciles similarity values; the input is a '{matrix.data_type}' distance matrix. "
                             "Symmetrize the score matrix, then convert it.")
        matrix.data = symmetrize_matrix(matrix.data, params.symmetrize)

    if requested_mode is not None:
        if matrix.data_type == "":
            warnings.warn("Input matrix format not specified, assuming 'score'")
            matrix.data_type = "score"

        if matrix.data_type != "score" and matrix.data_type != requested_mode:
            raise ValueError(f"Input matrix format is '{matrix.data_type}', but --mode is '{requested_mode}'. Transforming between data types is only supported from a 'score' matrix.")

        if sparsify and requested_mode == "bool" and matrix.data_type == "score":
            # Every bool cell is 1, which would leave --mst/--knn nothing to rank, so select
            # on the scores and keep the mask afterward. The kept cells are exactly the
            # nonzero scores, so --threshold still applies to the bool values below.
            matrix.data = sparsify_matrix(matrix, mst=params.mst, knn=params.knn, lower_bound=0)
            sparsify = False

        if matrix.data_type == "score" and matrix.data_type != requested_mode:
            matrix.data = transform_matrix(matrix.data, requested_mode, matrix.row_lengths, matrix.column_lengths)

        matrix.data_type = requested_mode

    if sparsify:
        check_sparsify_data_type(matrix.data_type, params.mst, params.knn)
        matrix.data = sparsify_matrix(matrix, mst=params.mst, knn=params.knn,
                                      lower_bound=0 if params.threshold is None else params.threshold)

    if params.threshold is not None:
        if is_distance_data_type(matrix.data_type):
            raise ValueError(f"--threshold zeroes out low values, which in a '{matrix.data_type}' (distance) matrix are the "
                             "closest pairs (0 means identical). Apply --threshold to the score matrix instead.")
        matrix.data = apply_lower_bound(matrix.data, params.threshold)

    ### Run
    try:
        write_matrix_outputs(
            matrix.data, matrix.rows, matrix.columns,
            hdf=hdf, dense_text=dense_text, max_output_bytes=max_output_bytes,
            mitigation_options=["--threshold", "--mst", "--knn"],
            row_lengths=matrix.row_lengths, col_lengths=matrix.column_lengths, data_type=matrix.data_type,
            symmetric=matrix.exactly_symmetric,
        )
    except OutputSizeLimitExceeded as exc:
        raise SystemExit(str(exc)) from None
