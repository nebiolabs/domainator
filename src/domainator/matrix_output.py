"""Shared matrix output options (-o/--hdf_out, --dense_text_out) for the matrix-producing tools.

``--hdf_out`` writes Domainator hdf5 (CSR; the upper triangle only when the values are
exactly symmetric). ``--dense_text_out`` writes tab-separated text. ``--symmetrize`` is the
shared option for reconciling the two directions of a self-comparison.
"""

from typing import Optional, Sequence, Tuple

from domainator.data_matrix import DataMatrix, SYMMETRIZE_CHOICES, get_matrix_file_type
from domainator.output_guardrails import enforce_matrix_output_limit


def add_matrix_output_arguments(parser, *, noun: str = "matrix") -> None:
    """Add -o/--hdf_out and --dense_text_out."""
    parser.add_argument("-o", "--hdf_out", type=str, default=None,
                        help=f"Write the {noun} to this Domainator hdf5 file. Values are stored sparse; "
                             "a matrix with exactly symmetric values is stored as its upper triangle.")
    parser.add_argument("--dense_text_out", type=str, default=None,
                        help=f"Write the {noun} to this tab-separated text file, every cell included. Use - for standard output.")


def resolve_matrix_outputs(params, *, require_output: bool = True) -> Tuple[Optional[str], Optional[str]]:
    """Return ``(hdf_path, dense_text_path)`` from parsed arguments, checking the hdf5 file extension."""
    hdf = params.hdf_out
    if hdf is not None and get_matrix_file_type(hdf) != "hdf5":
        raise ValueError("Please use an hdf5 related extension for the --hdf_out output, such as .h5, .hdf5, or .hdf.")
    if require_output and hdf is None and params.dense_text_out is None:
        raise ValueError("No output specified! Please specify at least one of: -o/--hdf_out, --dense_text_out")
    return hdf, params.dense_text_out


def write_matrix_outputs(matrix, row_names, col_names, *, hdf: Optional[str], dense_text: Optional[str],
                         max_output_bytes: Optional[int], mitigation_options: Sequence[str] = (),
                         row_lengths=None, col_lengths=None, data_type: str = "",
                         symmetric: Optional[bool] = None) -> None:
    """Check the size guardrail for, then write, each requested matrix output.

    ``matrix`` is a numpy or scipy array. ``symmetric`` is passed to
    :meth:`DataMatrix.write_hdf`. Raises ``OutputSizeLimitExceeded`` before writing an
    output that would exceed ``max_output_bytes``.
    """
    mitigation_options = list(mitigation_options)
    if hdf is not None:
        enforce_matrix_output_limit(
            output_type="hdf", matrix=matrix, row_names=row_names, col_names=col_names,
            row_lengths=row_lengths, col_lengths=col_lengths, data_type=data_type,
            max_output_bytes=max_output_bytes, output_path=hdf, mitigation_options=mitigation_options,
            symmetric=symmetric,
        )
        DataMatrix.write_hdf(matrix, hdf, row_names, col_names, row_lengths, col_lengths, data_type, symmetric=symmetric)
    if dense_text is not None:
        enforce_matrix_output_limit(
            output_type="dense_text", matrix=matrix, row_names=row_names, col_names=col_names,
            row_lengths=row_lengths, col_lengths=col_lengths, data_type=data_type,
            max_output_bytes=max_output_bytes, output_path=dense_text,
            mitigation_options=["--hdf_out"] + mitigation_options,
        )
        DataMatrix.write_dense_text(matrix, dense_text, row_names, col_names, row_lengths, col_lengths, data_type)


def add_symmetrize_argument(parser, *, extra_help: str = "") -> None:
    """Add the shared --symmetrize option (see :func:`domainator.data_matrix.symmetrize_matrix`)."""
    suffix = f" {extra_help}" if extra_help else ""
    parser.add_argument('--symmetrize', type=str, default="none", choices=list(SYMMETRIZE_CHOICES),
                        help="reconcile the two directions of each pair of a square self-comparison, on the raw "
                             "values before any --mode transform. none: leave the matrix as computed, which may be "
                             "slightly asymmetric. max: the better of the two directions in both cells, the "
                             "OR-symmetric convention used by --mst/--knn and by build_ssn. min: the worse of the "
                             "two -- note a pair found in only one direction has a stored zero in the other, so min "
                             "DELETES it, making this a mutual-hit filter. mean: the average of the two -- note a "
                             "pair found in only one direction is HALVED, because the missing direction contributes "
                             "zero. A symmetric result is stored as its upper triangle by --hdf_out, about half the size."
                             + suffix)
