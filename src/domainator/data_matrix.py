"""
DataMatrix: labeled 2D matrices backed by scipy CSR arrays.

Every matrix, whatever its source format, is held in memory as a canonical
``scipy.sparse.csr_array`` (explicit zeros removed, indices sorted). A zero cell means
"no value" for similarity scores and "identical" for distances; ``data_type`` tells the
two apart (distance types end in ``_dist``).

Supported inputs (``DataMatrix.from_file``):
    - Domainator hdf5, sparse CSR storage (full matrix or upper triangle).
    - Domainator hdf5, legacy dense storage (read only; converted to CSR on load).
    - Tab-separated dense text, e.g. matrices produced outside Domainator.

Supported outputs (``DataMatrix.write``):
    - ``"hdf"``: Domainator hdf5. Matrices whose values are exactly symmetric (with
      matching row and column labels) are stored as their upper triangle, including the
      diagonal, and mirrored on load.
    - ``"dense_text"``: tab-separated text.

Usage:
    matrix = DataMatrix.from_file("similarities.hdf5")
    matrix = DataMatrix(scipy_or_numpy_array, row_names, col_names)

    for row, col, value in matrix.iter_data():
        print(f"{row} vs {col}: {value}")

    matrix.write("out.hdf5", "hdf")
"""
import contextlib
import sys
import csv
from os import PathLike
from pathlib import Path
import scipy.sparse
import numpy as np
import h5py

from domainator.errors import FormatVersionError
from typing import Callable, Dict, List, Tuple, Union, Optional, Iterator

from .edge_types import SortedUndirectedEdges


_HDF5_FILE_OVERHEAD_ESTIMATE = 8 * 1024
_HDF5_DATASET_OVERHEAD_ESTIMATE = 4 * 1024
_HDF5_ATTRIBUTE_OVERHEAD_ESTIMATE = 2 * 1024

MATRIX_OUTPUT_TYPES = ("hdf", "dense_text")


def _sparse_index_dtype(matrix: scipy.sparse.csr_array) -> np.dtype:
    """Return the narrowest dtype that can index this CSR matrix.

    scipy keeps whatever index width an input carried, so a matrix derived from an
    int64-indexed file (e.g. a large all-vs-all matrix passed through
    ``transform_matrix --threshold``) stays int64 even when int32 would address it. Narrowing
    saves 4 bytes per non-zero on write. ``indices`` holds column indices (< n_cols) and
    ``indptr`` counts up to nnz, so int32 is safe whenever both fit.
    """
    max_index = max(matrix.shape[0], matrix.shape[1], matrix.nnz)
    if max_index <= np.iinfo(np.int32).max:
        return np.dtype(np.int32)
    return np.dtype(np.int64)


def _narrowed_sparse_index_arrays(matrix: scipy.sparse.csr_array) -> Tuple[np.ndarray, np.ndarray]:
    """Return ``(indices, indptr)`` cast to :func:`_sparse_index_dtype` (no copy if already narrow)."""
    dtype = _sparse_index_dtype(matrix)
    return matrix.indices.astype(dtype, copy=False), matrix.indptr.astype(dtype, copy=False)


def _to_canonical_csr(matrix, *, in_place: bool) -> scipy.sparse.csr_array:
    """Return ``matrix`` as a CSR array with no explicit zeros and sorted indices.

    With ``in_place=False`` a CSR input that is already canonical is returned as is, and
    one that is not is copied before being fixed, so the caller's array is never changed.
    With ``in_place=True`` a CSR input is canonicalized in place (no copy).
    """
    if not isinstance(matrix, scipy.sparse.csr_array):
        # Any other input (numpy, another sparse format, nested lists) is converted to a
        # new array, which can then be modified freely.
        if scipy.sparse.issparse(matrix):
            matrix = scipy.sparse.csr_array(matrix)
        else:
            matrix = scipy.sparse.csr_array(np.asarray(matrix))
        in_place = True

    if not in_place and (not matrix.has_canonical_format or np.any(matrix.data == 0)):
        matrix = matrix.copy()
    matrix.sum_duplicates()  # also sorts indices; before eliminate_zeros, since duplicates can sum to 0
    matrix.eliminate_zeros()
    return matrix


def _is_exactly_symmetric(matrix: scipy.sparse.csr_array) -> bool:
    """True if a canonical CSR matrix is square and equal to its transpose, value for value.

    Exact equality (not a tolerance) is what makes upper-triangle storage lossless.
    """
    if matrix.shape[0] != matrix.shape[1]:
        return False
    transpose = matrix.transpose().tocsr()
    transpose.sort_indices()
    return (
        np.array_equal(matrix.indptr, transpose.indptr)
        and np.array_equal(matrix.indices, transpose.indices)
        and np.array_equal(matrix.data, transpose.data)
    )


def mirror_upper_triangle(upper) -> scipy.sparse.csr_array:
    """Return the symmetric matrix whose upper triangle (diagonal included) is ``upper``.

    ``upper`` must hold no values below the diagonal. Returns canonical CSR.
    """
    upper = scipy.sparse.csr_array(upper)
    full = (upper + scipy.sparse.triu(upper, k=1, format="csr").transpose()).tocsr()
    full.sum_duplicates()
    full.eliminate_zeros()
    return full


def _labels_are_symmetric(row_names: List[str], col_names: List[str],
                          row_lengths: Optional[np.ndarray], col_lengths: Optional[np.ndarray]) -> bool:
    """True if the column axis can be stored as a reference to the row axis."""
    if list(row_names) != list(col_names):
        return False
    if row_lengths is None or col_lengths is None:
        return row_lengths is None and col_lengths is None
    return np.array_equal(np.asarray(row_lengths), np.asarray(col_lengths))


def _encoded_string_bytes(values: List[str]) -> int:
    return sum(len(value.encode("utf-8")) for value in values)


def _estimate_label_dataset_size(values: List[str]) -> int:
    if len(values) == 0:
        return _HDF5_DATASET_OVERHEAD_ESTIMATE
    return _HDF5_DATASET_OVERHEAD_ESTIMATE + _encoded_string_bytes(values) + (16 * len(values))


def _estimate_optional_array_size(array: Optional[np.ndarray]) -> int:
    if array is None:
        return 0
    return _HDF5_DATASET_OVERHEAD_ESTIMATE + int(np.asarray(array).nbytes)


def _empty_sorted_undirected_edges(n_nodes: int) -> SortedUndirectedEdges:
    return SortedUndirectedEdges(
        n_nodes=n_nodes,
        source=np.empty(0, dtype=np.int32),
        target=np.empty(0, dtype=np.int32),
        score=np.empty(0, dtype=float),
    )


def _parse_index_style(index_style: str) -> bool:
    """Validate ``index_style`` and return True if names (rather than indices) are wanted."""
    index_style = index_style.lower()
    if index_style not in {"name", "index"}:
        raise ValueError("index_style must be 'name' or 'index'.")
    return index_style == "name"


#TODO: use hdf5plugin to blosc compress datasets: http://www.silx.org/doc/hdf5plugin/latest/usage.html

UTF8_h5py_encoding = h5py.string_dtype(encoding='utf-8')

def get_matrix_file_type(filename: Union[PathLike, str]) -> Optional[str]:
    """Detect file type from extension.

    Args:
        filename: Path to file

    Returns:
        File type string ('hdf5') or None if not recognized
    """
    extension_map = {
        "hdf5": "hdf5",
        "hdf": "hdf5",
        "h5": "hdf5",
    }
    file_extension = Path(filename).suffix[1:].lower()
    return extension_map.get(file_extension)


SYMMETRIZE_CHOICES = ("none", "min", "max", "mean")


def symmetrize_matrix(matrix, how: str):
    """Reconcile the two directions of every pair of a square matrix.

    Apply it to raw similarity values, before any distance transform: it reconciles two
    measurements of one pair rather than adjusting a normalization, and 1 - x flips the
    ordering, so 'min' of a distance would be 'max' of a similarity.

    A pair reported in only one direction has a stored zero in the other, which means "no
    hit", not "measured as zero". min therefore deletes such a pair and mean halves it;
    only max treats a one-directional hit as evidence. max is the convention every SSN
    tool (and --mst/--knn) uses.

    Returns a CSR array for sparse input, a numpy array otherwise ("none" returns the input).
    """
    if how == "none":
        return matrix
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Only a square matrix can be symmetrized.")
    sparse = scipy.sparse.issparse(matrix)
    transposed = matrix.T
    if how == "max":
        out = matrix.maximum(transposed) if sparse else np.maximum(matrix, transposed)
    elif how == "min":
        out = matrix.minimum(transposed) if sparse else np.minimum(matrix, transposed)
    elif how == "mean":
        out = (matrix + transposed) / 2.0
    else:
        raise ValueError(f"Unknown symmetrize method: {how!r}, expected one of {SYMMETRIZE_CHOICES}")
    if scipy.sparse.issparse(out):
        out = scipy.sparse.csr_array(out)
        out.eliminate_zeros()
        out.sort_indices()
    return out


def is_distance_data_type(data_type: Optional[str]) -> bool:
    """True for distance-valued matrices, where 0 means identical rather than "no value"."""
    return bool(data_type) and data_type.endswith("_dist")


class DataMatrix:
    """
    A labeled 2D matrix stored as a canonical ``scipy.sparse.csr_array``.

    Attributes:
        data: scipy.sparse.csr_array holding the values. Assigning to it re-canonicalizes
            the new value. Changing ``data.data`` in place is allowed, but then call
            ``eliminate_zeros()`` on it if any value became 0.
        rows: List of row labels
        columns: List of column labels
        row_to_idx: Dictionary mapping row names to indices
        column_to_idx: Dictionary mapping column names to indices (same object as row_to_idx if symmetric)
        row_lengths: Sequence lengths for rows as numpy array, or None
        column_lengths: Sequence lengths for columns as numpy array, or None
        data_type: String describing the kind of values, empty string if not specified
    """
    _ROW_LABELS_DATASET="ROW_LABELS"
    """string array"""
    _COL_LABELS_DATASET="COL_LABELS"
    """string array (not used if SYMMETRIC LABELS is True)"""

    _ROW_LENGTHS_DATASET="ROW_LENGTHS"
    """(optional) int array representing the lengths of the sequences in the rows"""
    _COL_LENGTHS_DATASET="COL_LENGTHS"
    """(optional) int array (not used if SYMMETRIC LABELS is True)"""

    _DENSE_DATASET="DENSE_DATA"
    """2d array of floats, for legacy dense files (read only)"""
    _SPARSE_VALUES_DATASET="SPARSE_VALUES"
    """values for CSR sparse"""
    _SPARSE_CSR_INDICES_DATASET="SPARSE_CSR_INDICES"
    """column indices for CSR sparse"""
    _SPARSE_CSR_INDPTR_DATASET="SPARSE_CSR_INDPTR"
    """row pointer offsets for CSR sparse"""

    _ARRAY_TYPE_ATTR="ARRAY_TYPE"
    """{DENSE, SPARSE_CSR, SPARSE_CSR_UPPER_TRIANGLE}"""
    _SYMMETRIC_LABELS_ATTR="SYMMETRIC_LABELS"
    """if the x-axis and y-axis have the same labels"""
    _MATRIX_FILE_VERSION_ATTR="MATRIX_FILE_VERSION"
    _LEGACY_MATRIX_FILE_VERSION_ATTR="DOMAINATOR_MATRIX_FILE_VERSION"
    _DATA_TYPE_ATTR="DATA_TYPE"
    """(optional) single string value describing the type of data in the matrix (e.g. 'score', 'norm_score', 'row_norm_score', 'score_dist', 'bool', 'efi_score')"""


    _DENSE = "DENSE"
    """Legacy dense storage. Read only: Domainator no longer writes it."""
    _SPARSE_CSR = "SPARSE_CSR"
    _SPARSE_CSR_UPPER_TRIANGLE = "SPARSE_CSR_UPPER_TRIANGLE"
    """CSR holding only the upper triangle and diagonal of an exactly symmetric matrix.

    A distinct ARRAY_TYPE (rather than a new attribute) makes Domainator versions that
    predate it fail loudly on these files instead of silently reading half the matrix.
    """

    _MATRIX_FILE_VERSION = "2.0"
    _SUPPORTED_MATRIX_FILE_MAJOR_VERSIONS = ("1", "2")
    """Major versions this reader accepts. See "Matrix file versions" in docs/file_formats.md."""

    def __init__(self, data: Union[np.ndarray, scipy.sparse.sparray, scipy.sparse.spmatrix],
                 row_names: List[str], col_names: List[str],
                 row_lengths: Optional[np.ndarray] = None, col_lengths: Optional[np.ndarray] = None,
                 data_type: Optional[str] = None, *, copy_data: bool = True):
        """
        Args:
            data: 2D numpy array, scipy sparse array/matrix, or array-like. Converted to a
                canonical CSR array.
            row_names: List of row labels, length must match data.shape[0]
            col_names: List of column labels, length must match data.shape[1]
            row_lengths: Optional sequence lengths for rows (array-like, stored as numpy array)
            col_lengths: Optional sequence lengths for columns (array-like, stored as numpy array)
            data_type: Optional string describing the data type (e.g., 'score', 'norm_score', 'score_dist')
            copy_data: If False and ``data`` is already a CSR array, reuse it (canonicalizing
                it in place) instead of copying it.

        Raises:
            ValueError: If data dimensions don't match label lengths
        """
        DataMatrix.validate_matrix_data(data, row_names, col_names, row_lengths, col_lengths)

        self.rows = list(row_names)
        self.columns = list(col_names)
        self.row_to_idx = {v:i for i,v in enumerate(self.rows)}
        self.data_type = str(data_type) if data_type is not None else ""

        if self.rows == self.columns:
            self.column_to_idx = self.row_to_idx
        else:
            self.column_to_idx = {v:i for i,v in enumerate(self.columns)}

        if row_lengths is not None:
            self.row_lengths = np.array(row_lengths, copy=True)
            if col_lengths is row_lengths or (col_lengths is not None and np.array_equal(row_lengths, col_lengths)):
                self.column_lengths = self.row_lengths
            else:
                self.column_lengths = np.array(col_lengths, copy=True)
        else:
            self.row_lengths = None
            self.column_lengths = None

        if copy_data and isinstance(data, scipy.sparse.csr_array):
            data = data.copy()
        self._data = _to_canonical_csr(data, in_place=True)
        self._exactly_symmetric: Optional[bool] = None

    @property
    def data(self) -> scipy.sparse.csr_array:
        return self._data

    @data.setter
    def data(self, value):
        value = _to_canonical_csr(value, in_place=True)
        if value.shape != (len(self.rows), len(self.columns)):
            raise ValueError(f"New data shape {value.shape} does not match the labels "
                             f"({len(self.rows)}, {len(self.columns)}).")
        self._data = value
        self._exactly_symmetric = None

    @property
    def shape(self) -> Tuple[int, int]:
        """Return shape of the matrix as (rows, columns)."""
        return self._data.shape

    @property
    def nnz(self) -> int:
        """Return number of stored (non-zero) values."""
        return self._data.nnz

    @property
    def size(self) -> int:
        """
            Returns:
                int: number of cells in the matrix
        """
        # Use Python ints to avoid int32 overflow on platforms where np.prod
        # would return a 32-bit integer for large matrices.
        return int(self.shape[0]) * int(self.shape[1])

    @property
    def symmetric_labels(self) -> bool:
        """Returns True if row and column labels are identical."""
        return self.rows == self.columns

    @property
    def symmetric(self) -> bool:
        """Returns True if both labels and values are symmetric."""
        return self.symmetric_labels and self.symmetric_values

    @property
    def symmetric_values(self) -> bool:
        """
            Returns True if the matrix has symmetric values.
            Uses np.allclose to allow for floating point imprecision.
        """
        if not self.symmetric_labels:
            return False
        if self._exactly_symmetric:
            return True

        # Zeros are eliminated, so structural symmetry of the non-zero pattern is
        # required, then values are compared element-wise with a tolerance.
        A = self._data
        AT = A.transpose().tocsr()
        AT.sort_indices()

        if A.nnz != AT.nnz:
            return False
        if not np.array_equal(A.indptr, AT.indptr) or not np.array_equal(A.indices, AT.indices):
            return False
        return np.allclose(A.data, AT.data)

    @property
    def exactly_symmetric(self) -> bool:
        """True if labels are symmetric and every value equals its transpose exactly.

        This is the condition for upper-triangle storage. The result is cached until
        ``data`` is reassigned.
        """
        if not self.symmetric_labels:
            return False
        if self._exactly_symmetric is None:
            self._exactly_symmetric = _is_exactly_symmetric(self._data)
        return self._exactly_symmetric

    def __len__(self) -> int:
        """Return number of rows in the matrix."""
        return self.shape[0]

    def __repr__(self) -> str:
        """Return detailed string representation for debugging."""
        return (f"{self.__class__.__name__}(shape={self.shape}, nnz={self.nnz}, "
                f"symmetric_labels={self.symmetric_labels}, data_type='{self.data_type}')")

    def __str__(self) -> str:
        """Return readable string representation."""
        return f"{self.__class__.__name__}{self.shape}"

    @classmethod
    def _validate_hdf5_attributes(cls, f: h5py.File):
        """Validate required HDF5 attributes are present and the file version is readable.

        Supports both old (DOMAINATOR_MATRIX_FILE_VERSION) and new (MATRIX_FILE_VERSION)
        attribute names for backward compatibility.

        Raises:
            ValueError: If required attributes are missing
            FormatVersionError: If the major version is unknown
        """
        version = f.attrs.get(cls._MATRIX_FILE_VERSION_ATTR, f.attrs.get(cls._LEGACY_MATRIX_FILE_VERSION_ATTR))
        if version is None:
            raise ValueError(f"Required attribute '{cls._MATRIX_FILE_VERSION_ATTR}' not found in HDF5 file")
        if isinstance(version, bytes):
            version = version.decode("utf-8")
        if str(version).split(".")[0] not in cls._SUPPORTED_MATRIX_FILE_MAJOR_VERSIONS:
            supported = ", ".join(f"{major}.x" for major in cls._SUPPORTED_MATRIX_FILE_MAJOR_VERSIONS)
            raise FormatVersionError(f"Unsupported matrix file version '{version}'. This version of Domainator reads "
                             f"version {supported} matrix files.")

        for attribute in (cls._ARRAY_TYPE_ATTR, cls._SYMMETRIC_LABELS_ATTR):
            if attribute not in f.attrs:
                raise ValueError(f"Required attribute '{attribute}' not found in HDF5 file")

    @classmethod
    def from_file(cls, matrix_file: Union[PathLike, str], lower_bound: Optional[float] = None) -> 'DataMatrix':
        """Load a matrix from a Domainator hdf5 file or a tab-separated dense text file.

        Args:
            matrix_file: path to the file. hdf5 is recognized by its extension (.hdf5, .hdf, .h5);
                anything else is read as dense text.
            lower_bound: if set, values <= lower_bound are dropped (set to 0) after loading.
                NaNs are kept.
        """
        file_type = get_matrix_file_type(matrix_file)
        if file_type == "hdf5":
            matrix = cls._from_hdf5(matrix_file)
        elif file_type is None:
            matrix = cls._from_text(matrix_file)
        else:
            raise ValueError(f"Unrecognized file type: {file_type}")

        if lower_bound is not None:
            values = matrix._data.data
            # NaNs are passed through unchanged (NaN <= lower_bound is False).
            dropped = values <= lower_bound
            if np.any(dropped):
                # A symmetric matrix stays symmetric: equal values are dropped together.
                values[dropped] = 0
                matrix._data.eliminate_zeros()
        return matrix

    @classmethod
    def _from_hdf5(cls, matrix_file: Union[PathLike, str]) -> 'DataMatrix':
        """Load a matrix from a Domainator hdf5 file."""
        with h5py.File(matrix_file, 'r') as f:
            cls._validate_hdf5_attributes(f)
            array_type = f.attrs[cls._ARRAY_TYPE_ATTR]

            row_names, _ = cls.parse_axis_index(f[cls._ROW_LABELS_DATASET])
            if f.attrs[cls._SYMMETRIC_LABELS_ATTR]:
                col_names = row_names
            else:
                col_names, _ = cls.parse_axis_index(f[cls._COL_LABELS_DATASET])

            row_lengths = f[cls._ROW_LENGTHS_DATASET][:] if cls._ROW_LENGTHS_DATASET in f else None
            col_lengths = f[cls._COL_LENGTHS_DATASET][:] if cls._COL_LENGTHS_DATASET in f else None
            if f.attrs[cls._SYMMETRIC_LABELS_ATTR] and row_lengths is not None and col_lengths is None:
                col_lengths = row_lengths
            data_type = f.attrs.get(cls._DATA_TYPE_ATTR, "")

            exactly_symmetric = None
            if array_type in (cls._SPARSE_CSR, cls._SPARSE_CSR_UPPER_TRIANGLE):
                data = scipy.sparse.csr_array(
                    (f[cls._SPARSE_VALUES_DATASET][:],
                     f[cls._SPARSE_CSR_INDICES_DATASET][:],
                     f[cls._SPARSE_CSR_INDPTR_DATASET][:]),
                    (len(row_names), len(col_names))
                )
                if array_type == cls._SPARSE_CSR_UPPER_TRIANGLE:
                    data = mirror_upper_triangle(data)
                    exactly_symmetric = True
            elif array_type == cls._DENSE:
                data = f[cls._DENSE_DATASET][:]
            else:
                raise ValueError(f"Unrecognized matrix ARRAY_TYPE '{array_type}'.")

        matrix = cls(data, row_names, col_names, row_lengths, col_lengths, data_type, copy_data=False)
        matrix._exactly_symmetric = exactly_symmetric
        return matrix

    @classmethod
    def _from_text(cls, matrix_file: Union[PathLike, str]) -> 'DataMatrix':
        """Load a matrix from a tab-separated dense text file."""
        try:
            with open(matrix_file, newline="", encoding="utf-8") as handle:
                reader = csv.reader(handle, delimiter="\t")
                try:
                    header = next(reader)
                except StopIteration:
                    raise ValueError("Dense text matrix file is empty.")

                col_names = header[1:]
                if len(set(col_names)) != len(col_names):
                    raise ValueError("col_names must be unique")

                row_names = []
                for line_number, row in enumerate(reader, start=2):
                    if len(row) == 0:
                        continue
                    row_names.append(row[0])
                    if len(row) != len(header):
                        raise ValueError(
                            f"Dense text matrix row {line_number} has {len(row) - 1} values; expected {len(col_names)}."
                        )

                if len(set(row_names)) != len(row_names):
                    raise ValueError("row_names must be unique")
        except UnicodeDecodeError:
            raise ValueError(f"Expecting utf-8 encoded matrix file, but found invalid unicode byte.")

        import pandas as pd  # only dense text matrices need it
        try:
            mat_file = pd.read_csv(matrix_file, sep="\t", index_col=0)
        except UnicodeDecodeError:
            raise ValueError(f"Expecting utf-8 encoded matrix file, but found invalid unicode byte.")
        except pd.errors.ParserError as err:
            raise ValueError(f"Dense text matrix could not be parsed as tab-separated data: {err}") from err

        row_names = list(mat_file.index)
        col_names = list(mat_file.columns)
        try:
            data = mat_file.to_numpy(dtype=float)
        except ValueError as err:
            raise ValueError(f"Dense text matrix values must be numeric: {err}") from err

        return cls(data, row_names, col_names)

    def require_similarity(self, tool: str):
        """Raise if this matrix holds distances, for tools where 0 must mean "no edge"."""
        if is_distance_data_type(self.data_type):
            score_type = self.data_type[:-len("_dist")]
            raise ValueError(f"{tool} needs a similarity matrix (higher values mean more similar, 0 means no edge), "
                             f"but the input is a '{self.data_type}' distance matrix. Use the '{score_type}' matrix instead.")

    def _require_square(self, operation: str):
        if self.shape[0] != self.shape[1]:
            raise NotImplementedError(f"{operation} only supported for square matrices")

    def _max_symmetrized(self) -> scipy.sparse.csr_array:
        """Return ``max(M, M.T)`` element-wise, as canonical CSR (absent cells count as 0)."""
        if self._exactly_symmetric:
            return self._data
        symmetrized = self._data.maximum(self._data.transpose()).tocsr()
        symmetrized.eliminate_zeros()
        symmetrized.sort_indices()
        return symmetrized

    def iter_data(self, index_style: str = "name", skip_zeros: bool = True
                  ) -> Iterator[Tuple[Union[str, int], Union[str, int], float]]:
        """
        Iterate over matrix values in row-major order, yielding (row_identifier, col_identifier, value).

        Args:
            index_style: Either "name" (return row/column names) or "index" (return integer indices)
            skip_zeros: If True, only yield stored (non-zero) values
        """
        use_names = _parse_index_style(index_style)
        data = self._data

        if skip_zeros:
            for r in range(data.shape[0]):
                for ind in range(data.indptr[r], data.indptr[r+1]):
                    c = int(data.indices[ind])
                    if use_names:
                        yield self.rows[r], self.columns[c], data.data[ind]
                    else:
                        yield r, c, data.data[ind]
        else:
            # Densify one row at a time to avoid per-element O(log nnz) CSR scalar indexing.
            n_cols = data.shape[1]
            for r in range(data.shape[0]):
                row = data[[r]].toarray().ravel()
                for c in range(n_cols):
                    if use_names:
                        yield self.rows[r], self.columns[c], row[c]
                    else:
                        yield r, c, row[c]

    def triangular_iter(self, include_diagonal: bool = False, skip_zeros: bool = False, index_style: str = "name",
                        agg: Optional[Callable[[float, float], float]] = None
                        ) -> Iterator[Tuple[Union[str, int], Union[str, int], float]]:
        """
        Iterate over the lower triangle of a square matrix, in row-major order.

        Args:
            include_diagonal: If True, include diagonal elements
            skip_zeros: If True, exclude zero values from results
            index_style: Either "name" (return row/column names) or "index" (return integer indices)
            agg: Optional function combining each value with its transpose, e.g. ``max``,
                yielding ``agg(M[r, c], M[c, r])``. With ``skip_zeros``, a position is
                skipped when both cells are 0; for the builtin ``max`` it is also skipped
                when the result is 0 (which only happens with negative values).

        Yields:
            (row_identifier, col_identifier, value) tuples.

        Raises:
            NotImplementedError: If matrix is not square
        """
        self._require_square("Triangular")
        use_names = _parse_index_style(index_style)
        k = 0 if include_diagonal else -1

        def emit(r, c, value):
            if use_names:
                return self.rows[r], self.columns[c], value
            return r, c, value

        if skip_zeros and (agg is None or agg is max):
            working = self._data if agg is None else self._max_symmetrized()
            lower = scipy.sparse.tril(working, k=k, format="csr")
            lower.sort_indices()
            for r in range(lower.shape[0]):
                for ind in range(lower.indptr[r], lower.indptr[r + 1]):
                    yield emit(r, int(lower.indices[ind]), lower.data[ind])

        elif skip_zeros:
            # Generic aggregation: visit every lower-triangle position where (r,c) or (c,r)
            # is stored, then look both values up in one vectorized pass.
            pattern = (self._data != 0).astype(np.int8) + (self._data.transpose() != 0).astype(np.int8)
            pattern = scipy.sparse.tril(pattern, k=k, format="coo")
            order = np.lexsort((pattern.col, pattern.row))
            rows = pattern.row[order]
            cols = pattern.col[order]
            values = np.asarray(self._data[rows, cols]).ravel()
            transpose_values = np.asarray(self._data[cols, rows]).ravel()
            for r, c, value, transpose_value in zip(rows, cols, values, transpose_values):
                yield emit(int(r), int(c), agg(value, transpose_value))

        else:
            # Need all values including zeros - densify one row at a time.
            transpose = self._data.transpose().tocsr() if agg is not None else None
            for r in range(self._data.shape[0]):
                end_c = r + 1 if include_diagonal else r
                if end_c <= 0:
                    continue
                row = self._data[[r]].toarray().ravel()
                if agg is not None:
                    transpose_row = transpose[[r]].toarray().ravel()
                for c in range(end_c):
                    value = row[c]
                    if agg is not None:
                        value = agg(value, transpose_row[c])
                    yield emit(r, c, value)

    def triangular(self, include_diagonal: bool = False, skip_zeros: bool = False, index_style: str = "name",
                   agg: Optional[Callable[[float, float], float]] = None
                   ) -> List[Tuple[Union[str, int], Union[str, int], float]]:
        """Return the lower triangle of the matrix as a list of tuples (see ``triangular_iter``)."""
        return list(self.triangular_iter(
            include_diagonal=include_diagonal,
            skip_zeros=skip_zeros,
            index_style=index_style,
            agg=agg,
        ))

    def sorted_undirected_edges(self, skip_zeros: bool = False,
                                agg: Optional[Callable[[float, float], float]] = None) -> SortedUndirectedEdges:
        """Return lower-triangular undirected edges sorted by descending score.

        ``agg=None`` or the builtin ``max`` take a vectorized path; any other callable is
        applied element-wise via ``triangular``.
        """
        self._require_square("sorted_undirected_edges")

        if agg is not None and agg is not max:
            triangular = self.triangular(include_diagonal=False, skip_zeros=skip_zeros, index_style="index", agg=agg)
            if len(triangular) == 0:
                return _empty_sorted_undirected_edges(len(self))
            edge_array = np.asarray(triangular, dtype=float)
            order = np.argsort(-edge_array[:, 2], kind="stable")
            return SortedUndirectedEdges(
                n_nodes=len(self),
                source=edge_array[order, 0].astype(np.int32, copy=False),
                target=edge_array[order, 1].astype(np.int32, copy=False),
                score=edge_array[order, 2].astype(float, copy=False),
            )

        working = self._data if agg is None else self._max_symmetrized()
        if skip_zeros:
            lower = scipy.sparse.tril(working, k=-1, format="csr")
            lower.sort_indices()
            coo = lower.tocoo()  # row-major, column-sorted
            source, target, score = coo.row, coo.col, coo.data
        else:
            # Every pair, zeros included: materialize the lower triangle.
            n = self.shape[0]
            source, target = np.tril_indices(n, k=-1)
            score = working.toarray()[source, target]

        if len(score) == 0:
            return _empty_sorted_undirected_edges(len(self))

        # Pack source/target/score into a single structured array so sorting is
        # in-place, then sort stable ascending and reverse, so ties keep row-major order
        # reversed (the order every consumer and test expects).
        dtype = np.dtype([('score', np.float64), ('source', np.int32), ('target', np.int32)])
        edges = np.empty(len(score), dtype=dtype)
        edges['score'] = score
        edges['source'] = source
        edges['target'] = target
        del source, target, score  # free before sort

        edges.sort(order='score', kind='stable')
        edges = edges[::-1]  # descending view; base array kept alive via field views below
        return SortedUndirectedEdges(
            n_nodes=len(self),
            source=edges['source'],
            target=edges['target'],
            score=edges['score'],
        )

    def itervalues(self) -> Iterator[float]:
        """Yield all values row by row, zeros included."""
        # Densify one row at a time to avoid per-element CSR scalar indexing.
        for r in range(self.shape[0]):
            yield from self._data[[r]].toarray().ravel()

    def toarray(self) -> np.ndarray:
        """Return matrix as a dense numpy array (zeros included)."""
        return self._data.toarray()

    def write(self, file_name: Union[PathLike, str], output_type: str = "hdf"):
        """
        Write the matrix to a file.

        Args:
            file_name: Path to the output file
            output_type: "hdf" (Domainator hdf5) or "dense_text" (TSV)

        Raises:
            ValueError: If output_type is not recognized
        """
        if output_type == "hdf":
            self.write_hdf(self._data, file_name, self.rows, self.columns, self.row_lengths, self.column_lengths,
                           self.data_type, symmetric=self.exactly_symmetric)
        elif output_type == "dense_text":
            self.write_dense_text(self._data, file_name, self.rows, self.columns, self.row_lengths, self.column_lengths, self.data_type)
        else:
            raise ValueError(f"Unrecognized output type: {output_type}")

    @staticmethod
    def parse_axis_index(axis_index):
        items_list = [x.decode("utf-8") for x in axis_index]
        items_dict = {v:i for i,v in enumerate(items_list)}

        return items_list, items_dict

    @staticmethod
    def names_dict_to_list(names_dict:Dict[str, int]) -> List[str]:
        """Converts a dict of indexes into a list of keys

        Args:
            names_dict (Dict[str, int]): {row_name: row_index}

        Returns:
            List[str]: List of strings where list[row_index] = row_name
        """
        out = [None] * len(names_dict)
        for name, idx in names_dict.items():
            out[idx] = name
        if any(v is None for v in out):
            raise ValueError("names_dict is not a unique mapping.")
        return out

    @classmethod
    def validate_matrix_data(cls, matrix: Union[np.array, scipy.sparse.sparray], row_names:List[str], col_names:List[str], row_lens:Optional[np.ndarray]=None, col_lens:Optional[np.ndarray]=None) -> Tuple[List[str], List[str]]:
        """Check for consistency between a matrix and row/column labels

        Args:
            matrix (Union[np.array, scipy.sparse.sparray]): dense or sparse data matrix
            row_names:List[str]: row labels
            col_names:List[str]: column labels

        Returns:
            (row_names, col_names) as lists; col_names is the same object as row_names if they are equal.
        """

        row_names = list(row_names)
        col_names = list(col_names)

        shape = matrix.shape if hasattr(matrix, "shape") else np.shape(matrix)
        if len(shape) != 2:
            raise ValueError(f"matrix must be 2-dimensional, not: {shape}")

        if len(row_names) != shape[0]:
            raise ValueError(f"row_names size must match matrix dim 0, not: {len(row_names)} vs. {shape[0]}")

        if len(col_names) != shape[1]:
            raise ValueError(f"col_names size must match matrix dim 1, not: {len(col_names)} vs. {shape[1]}")

        if len(set(row_names)) != len(row_names):
            raise ValueError("row_names must be unique")
        if len(set(col_names)) != len(col_names):
            raise ValueError("col_names must be unique")

        if row_lens is not None:
            if len(row_lens) != shape[0]:
                raise ValueError(f"row_lens size must match matrix dim 0, not: {len(row_lens)} vs. {shape[0]}")
        if col_lens is not None:
            if len(col_lens) != shape[1]:
                raise ValueError(f"col_lens size must match matrix dim 1, not: {len(col_lens)} vs. {shape[1]}")

        if row_lens is not None and col_lens is None:
            raise ValueError("If row_lens is provided, col_lens must also be provided.")
        if row_lens is None and col_lens is not None:
            raise ValueError("If col_lens is provided, row_lens must also be provided.")

        if row_names == col_names:
            col_names = row_names

        return row_names, col_names

    @classmethod
    def _hdf_payload(cls, matrix, row_names, col_names, row_lengths, col_lengths,
                     symmetric: Optional[bool]) -> Tuple[scipy.sparse.csr_array, str, bool]:
        """Return ``(csr_to_store, array_type, symmetric_labels)`` for an hdf write.

        The caller's matrix is never modified.
        """
        csr = _to_canonical_csr(matrix, in_place=False)
        symmetric_labels = _labels_are_symmetric(row_names, col_names, row_lengths, col_lengths)
        if symmetric_labels and symmetric is not False and (symmetric or _is_exactly_symmetric(csr)):
            upper = scipy.sparse.triu(csr, k=0, format="csr")
            upper.sort_indices()
            return upper, cls._SPARSE_CSR_UPPER_TRIANGLE, symmetric_labels
        return csr, cls._SPARSE_CSR, symmetric_labels

    @classmethod
    def write_hdf(cls, matrix: Union[np.ndarray, scipy.sparse.sparray, scipy.sparse.spmatrix],
                  file_name: Union[PathLike, str], row_names: List[str], col_names: List[str],
                  row_lengths: Optional[np.ndarray] = None, col_lengths: Optional[np.ndarray] = None,
                  data_type: str = "", symmetric: Optional[bool] = None) -> None:
        """Writes a matrix to a Domainator hdf5 file in CSR format.

        Zeros are dropped and indices sorted. When the row and column labels (and lengths)
        are the same and the values are exactly symmetric, only the upper triangle
        (including the diagonal) is stored.

        Args:
            matrix: numpy array or scipy sparse array/matrix. Not modified.
            symmetric: None (default) detects exact symmetry. True asserts the values are
                exactly symmetric (skips the check; a wrong assertion loses the lower
                triangle). False always stores the full matrix.
        """
        row_names, col_names = cls.validate_matrix_data(matrix, row_names, col_names, row_lengths, col_lengths)
        stored, array_type, symmetric_labels = cls._hdf_payload(matrix, row_names, col_names, row_lengths, col_lengths, symmetric)

        with h5py.File(file_name, 'w') as f:
            f.attrs[cls._MATRIX_FILE_VERSION_ATTR] = cls._MATRIX_FILE_VERSION
            f.attrs[cls._ARRAY_TYPE_ATTR] = array_type

            # Narrow the index arrays if they are wider than this matrix needs; the read
            # path hands them straight to csr_array, which accepts either width.
            indices, indptr = _narrowed_sparse_index_arrays(stored)

            f.create_dataset(cls._SPARSE_VALUES_DATASET, data=stored.data)
            f.create_dataset(cls._SPARSE_CSR_INDICES_DATASET, data=indices)
            f.create_dataset(cls._SPARSE_CSR_INDPTR_DATASET, data=indptr)

            f.create_dataset(cls._ROW_LABELS_DATASET, data=row_names, dtype=UTF8_h5py_encoding)

            if row_lengths is not None:
                f.create_dataset(cls._ROW_LENGTHS_DATASET, data=row_lengths)

            if symmetric_labels:
                # Store as int8 (0/1), not a Python/NumPy bool: HDF5 has no native
                # boolean type, so h5py writes numpy bool_ as an enum that non-h5py
                # readers may silently drop. int8 is unambiguously portable and is
                # still read back truthily. (0 == not symmetric, 1 == symmetric.)
                f.attrs[cls._SYMMETRIC_LABELS_ATTR] = np.int8(1)
            else:
                f.attrs[cls._SYMMETRIC_LABELS_ATTR] = np.int8(0)
                f.create_dataset(cls._COL_LABELS_DATASET, data=col_names, dtype=UTF8_h5py_encoding)
                if col_lengths is not None:
                    f.create_dataset(cls._COL_LENGTHS_DATASET, data=col_lengths)

            f.attrs[cls._DATA_TYPE_ATTR] = data_type

    @classmethod
    def write_dense_text(cls, matrix:Union[np.array, scipy.sparse.sparray], file_name: Union[PathLike, str], row_names:List[str], col_names:List[str], row_lengths:Optional[np.ndarray]=None, col_lengths:Optional[np.ndarray]=None, data_type:Optional[str]=""):
        """Writes a matrix to a tab-separated text file, every cell included.

        Args:
            matrix (Union[np.array, scipy.sparse.sparray]): 2d array
            file_name (Union[PathLike, str]): path to the target file to write the data to. "-" or
                "/dev/stdout" writes standard output.
            row_names:List[str]: row labels
            col_names:List[str]: column labels
        """

        row_names, col_names = cls.validate_matrix_data(matrix, row_names, col_names, row_lengths, col_lengths)
        sparse = scipy.sparse.issparse(matrix)
        if sparse:
            matrix = scipy.sparse.csr_array(matrix)
        to_stdout = str(file_name) in ("-", "/dev/stdout")
        with (contextlib.nullcontext(sys.stdout) if to_stdout else open(file_name, "w")) as outfile:
            print("\t" + "\t".join(col_names), file=outfile)
            for idx in range(len(row_names)):
                row = matrix[[idx], :].toarray() if sparse else matrix[idx]
                print(row_names[idx] + "\t" + "\t".join([str(x) for x in np.nditer(row)]), file=outfile)

    @classmethod
    def estimate_hdf_size(
        cls,
        matrix: Union[np.ndarray, scipy.sparse.sparray],
        row_names: List[str],
        col_names: List[str],
        row_lengths: Optional[np.ndarray] = None,
        col_lengths: Optional[np.ndarray] = None,
        data_type: str = "",
        symmetric: Optional[bool] = None,
    ) -> int:
        """Estimate the size of ``write_hdf`` output in bytes (``symmetric`` as in ``write_hdf``)."""
        row_names, col_names = cls.validate_matrix_data(matrix, row_names, col_names, row_lengths, col_lengths)
        stored, _, symmetric_labels = cls._hdf_payload(matrix, row_names, col_names, row_lengths, col_lengths, symmetric)

        estimated_size = _HDF5_FILE_OVERHEAD_ESTIMATE + _HDF5_ATTRIBUTE_OVERHEAD_ESTIMATE
        estimated_size += _HDF5_DATASET_OVERHEAD_ESTIMATE + int(stored.data.nbytes)
        # Match write_hdf, which narrows the index arrays before writing them.
        index_itemsize = _sparse_index_dtype(stored).itemsize
        estimated_size += _HDF5_DATASET_OVERHEAD_ESTIMATE + len(stored.indices) * index_itemsize
        estimated_size += _HDF5_DATASET_OVERHEAD_ESTIMATE + len(stored.indptr) * index_itemsize
        estimated_size += _estimate_label_dataset_size(row_names)
        estimated_size += _estimate_optional_array_size(row_lengths)

        if not symmetric_labels:
            estimated_size += _estimate_label_dataset_size(col_names)
            estimated_size += _estimate_optional_array_size(col_lengths)

        estimated_size += len(data_type.encode("utf-8"))
        return estimated_size

    @classmethod
    def estimate_dense_text_size(
        cls,
        matrix: Union[np.ndarray, scipy.sparse.sparray],
        row_names: List[str],
        col_names: List[str],
        row_lengths: Optional[np.ndarray] = None,
        col_lengths: Optional[np.ndarray] = None,
        data_type: str = "",
    ) -> int:
        row_names, col_names = cls.validate_matrix_data(matrix, row_names, col_names, row_lengths, col_lengths)

        total_bytes = len(("\t" + "\t".join(col_names) + "\n").encode("utf-8"))
        if scipy.sparse.issparse(matrix):
            matrix = scipy.sparse.csr_array(matrix)
            for idx, row_name in enumerate(row_names):
                row_text = row_name + "\t" + "\t".join(str(x) for x in np.nditer(matrix[[idx], :].toarray())) + "\n"
                total_bytes += len(row_text.encode("utf-8"))
        else:
            for idx, row_name in enumerate(row_names):
                row_text = row_name + "\t" + "\t".join(str(x) for x in np.nditer(matrix[idx])) + "\n"
                total_bytes += len(row_text.encode("utf-8"))
        return total_bytes

    @classmethod
    def estimate_write_size(
        cls,
        output_type: str,
        matrix: Union[np.ndarray, scipy.sparse.sparray],
        row_names: List[str],
        col_names: List[str],
        row_lengths: Optional[np.ndarray] = None,
        col_lengths: Optional[np.ndarray] = None,
        data_type: str = "",
        symmetric: Optional[bool] = None,
    ) -> int:
        """Estimate the size of a matrix output in bytes; ``symmetric`` only affects "hdf"."""
        if output_type == "hdf":
            return cls.estimate_hdf_size(matrix, row_names, col_names, row_lengths, col_lengths, data_type, symmetric)
        if output_type == "dense_text":
            return cls.estimate_dense_text_size(matrix, row_names, col_names, row_lengths, col_lengths, data_type)
        raise ValueError(f"Unrecognized output type: {output_type}")

