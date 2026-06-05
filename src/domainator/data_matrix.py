"""
DataMatrix: Unified interface for dense and sparse matrices.

This module provides an abstract base class DataMatrix with two concrete implementations:
- DenseDataMatrix: Uses numpy.ndarray for storage (best for matrices with many non-zeros)
- SparseDataMatrix: Uses scipy.sparse.csr_array for storage (best for matrices with many zeros)

The abstraction allows code to work with both dense and sparse matrices through a common
interface, with each implementation optimized for its storage format.

Key Features:
    - Polymorphic interface: Methods work on both dense and sparse matrices
    - Type-based dispatch: Use isinstance() to check matrix type
    - Factory pattern: DataMatrix.from_file() returns the appropriate subclass
    - Conversion methods: convert_to_sparse() and convert_to_dense()
    - HDF5 persistence: Efficient binary storage with metadata
    - Row/column labels: String identifiers for matrix axes
    - Optional sequence lengths: Support for biological sequence data
    - Value label: an string describing the kind of data stored in the matrix, (e.g. e-values, score, etc.)

Usage:
    # Load from file (returns appropriate subclass)
    matrix = DataMatrix.from_file("distances.hdf5")
    
    # Create directly
    dense = DenseDataMatrix(data_array, row_names, col_names)
    sparse = SparseDataMatrix(sparse_array, row_names, col_names)
    
    # Check type
    if isinstance(matrix, SparseDataMatrix):
        matrix = matrix.convert_to_dense()
    
    # Iterate over values
    for row, col, value in matrix.iter_data():
        print(f"{row} vs {col}: {value}")

"""
import warnings
warnings.filterwarnings("ignore", module='numpy')
import csv
from os import PathLike
from pathlib import Path
import pandas as pd
import scipy.sparse
import numpy as np
import h5py
from abc import ABC, abstractmethod
from typing import Callable, Dict, List, Tuple, Union, Optional, Iterator

# SortedUndirectedEdges / CompactNeighborRankings live in edge_types to avoid an
# import cycle with ssn_edges (MaxTree + the edge/neighbor-ranking functions).
from .edge_types import SortedUndirectedEdges, CompactNeighborRankings


_HDF5_FILE_OVERHEAD_ESTIMATE = 8 * 1024
_HDF5_DATASET_OVERHEAD_ESTIMATE = 4 * 1024
_HDF5_ATTRIBUTE_OVERHEAD_ESTIMATE = 2 * 1024


def _triangle_col_bounds(row_idx: int, side: str, include_diagonal: bool, n_cols: int) -> Tuple[int, int]:
    """Return the [start, end) column range for the triangular part of one row."""
    if side == "lower":
        return 0, (row_idx + 1 if include_diagonal else row_idx)
    return (row_idx if include_diagonal else row_idx + 1), n_cols


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


#TODO: use hdf5plugin to blosc compress datasets: http://www.silx.org/doc/hdf5plugin/latest/usage.html

UTF8_h5py_encoding = h5py.string_dtype(encoding='utf-8')

def _get_file_type(filename: Union[PathLike, str]) -> Optional[str]:
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

class DataMatrix(ABC):
    """
    Class used to abstract access to sparse and dense matrices

    Dense matrices are stored as numpy arrays
    Sparse matrices are stored as scipy.sparse.csr_array
    Intended to be read-only
    
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
    """2d array of floats, for dense data"""
    _SPARSE_VALUES_DATASET="SPARSE_VALUES" 
    """values for CSR sparse"""
    _SPARSE_CSR_INDICES_DATASET="SPARSE_CSR_INDICES"
    """column indices for CSR sparse"""
    _SPARSE_CSR_INDPTR_DATASET="SPARSE_CSR_INDPTR"
    """row pointer offsets for CSR sparse"""
    
    _ARRAY_TYPE_ATTR="ARRAY_TYPE"
    """{DENSE, SPARSE_CSR}"""
    _SYMMETRIC_LABELS_ATTR="SYMMETRIC_LABELS"
    """if the x-axis and y-axis have the same labels"""
    _MATRIX_FILE_VERSION_ATTR="MATRIX_FILE_VERSION"
    _DATA_TYPE_ATTR="DATA_TYPE"
    """(optional) single string value describing the type of data in the matrix (e.g. 'score', 'norm_score', 'row_norm_score', 'score_dist', 'bool', 'efi_score')"""


    _DENSE = "DENSE"
    _SPARSE_CSR = "SPARSE_CSR"
    _POSSIBLE_ARRAY_TYPES={_DENSE, _SPARSE_CSR}


    
    _MATRIX_FILE_VERSION = "1.0"  #TODO: increment after any breaking changes

    def __init__(self, row_names:List[str], col_names:List[str], row_lengths:Optional[np.ndarray]=None, col_lengths:Optional[np.ndarray]=None, data_type:str=None):
        """
        Base class initialization - common to all matrix types.
        Subclasses should call this via super().__init__()
        
        Args:
            row_names: List of row labels (typically gene/sequence names)
            col_names: List of column labels (typically gene/sequence names)
            row_lengths: Optional sequence lengths for rows as numpy array (will be converted if not)
            col_lengths: Optional sequence lengths for columns as numpy array (will be converted if not)
            data_type: Optional string describing the data type (e.g., 'score', 'norm_score', 'score_dist')
        
        Attributes:
            rows: List of row labels
            columns: List of column labels
            row_to_idx: Dictionary mapping row names to indices
            column_to_idx: Dictionary mapping column names to indices (same as row_to_idx if symmetric)
            row_lengths: Sequence lengths for rows as numpy array, or None
            column_lengths: Sequence lengths for columns as numpy array, or None
            data_type: String describing data type, empty string if not specified
        """
        # Validate and set up row/column names
        self.rows = list(row_names)
        self.columns = list(col_names)
        self.row_to_idx = {v:i for i,v in enumerate(row_names)}
        self.data_type = str(data_type) if data_type is not None else ""

        if row_names == col_names:
            self.column_to_idx = self.row_to_idx
        else:
            self.column_to_idx = {v:i for i,v in enumerate(col_names)}

        if row_lengths is not None:
            # Always store as numpy array for consistency and mathematical operations
            if not isinstance(row_lengths, np.ndarray):
                self.row_lengths = np.array(row_lengths)
            else:
                self.row_lengths = row_lengths.copy()
            
            if col_lengths is row_lengths or (col_lengths is not None and np.array_equal(row_lengths, col_lengths)):
                self.column_lengths = self.row_lengths
            else:
                if not isinstance(col_lengths, np.ndarray):
                    self.column_lengths = np.array(col_lengths)
                else:
                    self.column_lengths = col_lengths.copy()
        else:
            self.row_lengths = None
            self.column_lengths = None

    @property
    @abstractmethod
    def shape(self):
        """Return the shape of the matrix"""
        pass
    
    @property
    def size(self) -> int:
        """
            Returns:
                int: number of values in the matrix
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
    @abstractmethod
    def symmetric_values(self) -> bool:
        """
            Returns True if the matrix has symmetric values (i.e., the matrix is symmetric in its data).
            Uses np.isclose to allow for floating point imprecision.
        """
        pass

    def __len__(self) -> int:
        """Return number of rows in the matrix."""
        return self.shape[0]
    
    def __repr__(self) -> str:
        """Return detailed string representation for debugging."""
        return (f"{self.__class__.__name__}(shape={self.shape}, "
                f"symmetric_labels={self.symmetric_labels}, data_type='{self.data_type}')")
    
    def __str__(self) -> str:
        """Return readable string representation."""
        return f"{self.__class__.__name__}{self.shape}"

    @classmethod
    def _validate_hdf5_attributes(cls, f: h5py.File):
        """Validate required HDF5 attributes are present.
        
        Supports both old (DOMAINATOR_MATRIX_FILE_VERSION) and new (MATRIX_FILE_VERSION)
        attribute names for backward compatibility.
        
        Args:
            f: Open HDF5 file handle
        
        Raises:
            ValueError: If required attributes are missing
        """
        # Check for version attribute (support both old and new names)
        if cls._MATRIX_FILE_VERSION_ATTR not in f.attrs and "DOMAINATOR_MATRIX_FILE_VERSION" not in f.attrs:
            raise ValueError(f"Required attribute '{cls._MATRIX_FILE_VERSION_ATTR}' not found in HDF5 file")
        
        # Check other required attributes
        for attribute in (cls._ARRAY_TYPE_ATTR, cls._SYMMETRIC_LABELS_ATTR):
            if attribute not in f.attrs:
                raise ValueError(f"Required attribute '{attribute}' not found in HDF5 file")

    @classmethod
    def from_file(cls, matrix_file, lower_bound=None):
        """
        Factory method that returns the appropriate subclass (DenseDataMatrix or SparseDataMatrix)
        based on the file content.
        """
        file_type = _get_file_type(matrix_file)
        if file_type == "hdf5":
            # Peek at the file to determine type. The chosen subclass's _from_hdf5
            # validates required attributes (including ARRAY_TYPE), so a missing
            # ARRAY_TYPE here surfaces as a clear error there rather than a silent default.
            with h5py.File(matrix_file, 'r') as f:
                array_type = f.attrs.get(cls._ARRAY_TYPE_ATTR)

                if array_type == cls._SPARSE_CSR:
                    return SparseDataMatrix._from_hdf5(matrix_file, lower_bound=lower_bound)
                else:
                    return DenseDataMatrix._from_hdf5(matrix_file, lower_bound=lower_bound)
        elif file_type is None:
            # Assume dense tsv text file if extension not recognized
            return DenseDataMatrix._from_text(matrix_file)
        else:
            raise ValueError(f"Unrecognized file type: {file_type}")

    @abstractmethod
    def iter_data(self, index_style="name", skip_zeros=True):
        """
        Iterate over matrix values, yielding (row_identifier, col_identifier, value).
        
        Args:
            index_style: Either "name" (return row/column names) or "index" (return integer indices)
            skip_zeros: If True, skip zero values (efficient for sparse matrices)
        
        Yields:
            tuple: (row_identifier, col_identifier, value) where identifiers are strings 
                   if index_style="name" or integers if index_style="index"
        """
        pass

    @abstractmethod
    def triangular_iter(self, side="lower", include_diagonal=False, skip_zeros=False, index_style="name", agg=None):
        """
        Iterate over the triangular part of the matrix.
        
        Only works for square matrices. Returns values from either upper or lower triangle.
        
        Args:
            side: "lower" or "upper" triangle
            include_diagonal: If True, include diagonal elements
            skip_zeros: If True, exclude zero values from results
            index_style: Either "name" (return row/column names) or "index" (return integer indices)
            agg: Optional aggregation function to combine symmetric positions (e.g., lambda a,b: max(a,b)).
                 When provided, applies agg(value, transpose_value) for each position.
        
        Yields:
            tuple: (row_identifier, col_identifier, value)
        
        Raises:
            NotImplementedError: If matrix is not square
        """
        pass

    def triangular(self, side="lower", include_diagonal=False, skip_zeros=False, index_style="name", agg=None):
        """Return triangular part of the matrix as a list of tuples."""
        return list(self.triangular_iter(
            side=side,
            include_diagonal=include_diagonal,
            skip_zeros=skip_zeros,
            index_style=index_style,
            agg=agg,
        ))

    @abstractmethod
    def itervalues(self) -> Iterator[float]:
        """
        Yield all values from the matrix one by one, row by row, column by column.
        
        Yields:
            Numeric values from the matrix in row-major order
        """
        pass

    @abstractmethod
    def convert_to_sparse(self) -> 'SparseDataMatrix':
        """
        Convert to sparse representation.
        
        Returns:
            SparseDataMatrix: A new sparse matrix (if dense) or self (if already sparse)
        """
        pass

    @abstractmethod
    def convert_to_dense(self) -> 'DenseDataMatrix':
        """
        Convert to dense representation.
        
        Returns:
            DenseDataMatrix: A new dense matrix (if sparse) or self (if already dense)
        """
        pass

    @abstractmethod
    def toarray(self) -> np.ndarray:
        """
        Return matrix as a dense numpy array.
        
        Returns:
            np.ndarray: 2D numpy array with all matrix values (zeros included)
        """
        pass

    def write(self, file_name: Union[PathLike, str], output_type: str = "dense"):
        """
        Write the matrix to a file in the specified format.
        
        Args:
            file_name: Path to the output file
            output_type: Format to write - "dense" (HDF5), "sparse" (HDF5), or "dense_text" (TSV)
        
        Raises:
            ValueError: If output_type is not recognized
        """
        if output_type == "dense":
            self.write_dense(self.data, file_name, self.rows, self.columns, self.row_lengths, self.column_lengths, self.data_type)
        elif output_type == "sparse":
            self.write_sparse(self.data, file_name, self.rows, self.columns, self.row_lengths, self.column_lengths, self.data_type)
        elif output_type == "dense_text":
            self.write_dense_text(self.data, file_name, self.rows, self.columns, self.row_lengths, self.column_lengths, self.data_type)
        else:
            raise ValueError(f"Unrecognized output type: {output_type}")

    @classmethod
    def estimate_dense_hdf5_size(
        cls,
        matrix: Union[np.ndarray, scipy.sparse.sparray],
        row_names: List[str],
        col_names: List[str],
        row_lengths: Optional[np.ndarray] = None,
        col_lengths: Optional[np.ndarray] = None,
        data_type: str = "",
    ) -> int:
        row_names, col_names = cls.validate_matrix_data(matrix, row_names, col_names, row_lengths, col_lengths)
        dtype = matrix.dtype if hasattr(matrix, "dtype") else np.float64
        data_bytes = int(matrix.shape[0] * matrix.shape[1] * np.dtype(dtype).itemsize)

        estimated_size = _HDF5_FILE_OVERHEAD_ESTIMATE + _HDF5_ATTRIBUTE_OVERHEAD_ESTIMATE
        estimated_size += _HDF5_DATASET_OVERHEAD_ESTIMATE + data_bytes
        estimated_size += _estimate_label_dataset_size(row_names)
        estimated_size += _estimate_optional_array_size(row_lengths)

        if row_names is not col_names or row_lengths is not col_lengths:
            estimated_size += _estimate_label_dataset_size(col_names)
            estimated_size += _estimate_optional_array_size(col_lengths)

        estimated_size += len(data_type.encode("utf-8"))
        return estimated_size

    @classmethod
    def estimate_sparse_hdf5_size(
        cls,
        matrix: Union[np.ndarray, scipy.sparse.sparray],
        row_names: List[str],
        col_names: List[str],
        row_lengths: Optional[np.ndarray] = None,
        col_lengths: Optional[np.ndarray] = None,
        data_type: str = "",
    ) -> int:
        row_names, col_names = cls.validate_matrix_data(matrix, row_names, col_names, row_lengths, col_lengths)
        if not isinstance(matrix, scipy.sparse.csr_array):
            matrix = scipy.sparse.csr_array(matrix)

        matrix.eliminate_zeros()
        matrix.sort_indices()

        estimated_size = _HDF5_FILE_OVERHEAD_ESTIMATE + _HDF5_ATTRIBUTE_OVERHEAD_ESTIMATE
        estimated_size += _HDF5_DATASET_OVERHEAD_ESTIMATE + int(matrix.data.nbytes)
        estimated_size += _HDF5_DATASET_OVERHEAD_ESTIMATE + int(matrix.indices.nbytes)
        estimated_size += _HDF5_DATASET_OVERHEAD_ESTIMATE + int(matrix.indptr.nbytes)
        estimated_size += _estimate_label_dataset_size(row_names)
        estimated_size += _estimate_optional_array_size(row_lengths)

        if row_names is not col_names or row_lengths is not col_lengths:
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
    ) -> int:
        if output_type == "dense":
            return cls.estimate_dense_hdf5_size(matrix, row_names, col_names, row_lengths, col_lengths, data_type)
        if output_type == "sparse":
            return cls.estimate_sparse_hdf5_size(matrix, row_names, col_names, row_lengths, col_lengths, data_type)
        if output_type == "dense_text":
            return cls.estimate_dense_text_size(matrix, row_names, col_names, row_lengths, col_lengths, data_type)
        raise ValueError(f"Unrecognized output type: {output_type}")

    @staticmethod
    def parse_axis_index(axis_index):
        items_list = [x.decode("utf-8") for x in axis_index]
        items_dict = {v:i for i,v in enumerate(items_list)}
        
        return items_list, items_dict

    def sorted_undirected_edges(self, skip_zeros: bool = False, agg=None) -> SortedUndirectedEdges:
        """Return lower-triangular undirected edges sorted by descending score."""
        triangular = self.triangular(
            side="lower",
            include_diagonal=False,
            skip_zeros=skip_zeros,
            index_style="index",
            agg=agg,
        )
        if len(triangular) == 0:
            return SortedUndirectedEdges(
                n_nodes=len(self),
                source=np.empty(0, dtype=np.int32),
                target=np.empty(0, dtype=np.int32),
                score=np.empty(0, dtype=float),
            )

        edge_array = np.asarray(triangular, dtype=float)
        order = np.argsort(-edge_array[:, 2], kind="stable")
        return SortedUndirectedEdges(
            n_nodes=len(self),
            source=edge_array[order, 0].astype(np.int32, copy=False),
            target=edge_array[order, 1].astype(np.int32, copy=False),
            score=edge_array[order, 2].astype(float, copy=False),
        )


    @classmethod
    def write_sparse(cls, matrix, file_name, row_names, col_names, row_lengths:Optional[np.ndarray]=None, col_lengths:Optional[np.ndarray]=None, data_type=""):
        """Writes a matrix to an hdf5 file in CSR sparse format.

        The input is converted to scipy.sparse.csr_array if needed; zeros are
        eliminated and indices sorted before writing.
        """

        # check validity of input data
        row_names, col_names = cls.validate_matrix_data(matrix, row_names, col_names, row_lengths, col_lengths)
        
        if not isinstance(matrix, scipy.sparse.csr_array):
            matrix = scipy.sparse.csr_array(matrix)

        matrix.eliminate_zeros() # remove zeros, so that the sparse file is as small as possible
        matrix.sort_indices() # sort indices so iteration happens by row and column in order

        with h5py.File(file_name, 'w') as f:
            f.attrs[cls._MATRIX_FILE_VERSION_ATTR] = cls._MATRIX_FILE_VERSION
            f.attrs[cls._ARRAY_TYPE_ATTR] = cls._SPARSE_CSR

            f.create_dataset(cls._SPARSE_VALUES_DATASET, data=matrix.data)
            f.create_dataset(cls._SPARSE_CSR_INDICES_DATASET, data=matrix.indices)
            f.create_dataset(cls._SPARSE_CSR_INDPTR_DATASET, data=matrix.indptr)

            f.create_dataset(cls._ROW_LABELS_DATASET, data=row_names, dtype=UTF8_h5py_encoding)
            
            if row_lengths is not None:
                f.create_dataset(cls._ROW_LENGTHS_DATASET, data=row_lengths)


            if row_names is col_names and row_lengths is col_lengths:
                f.attrs[cls._SYMMETRIC_LABELS_ATTR] = True
            else:
                f.attrs[cls._SYMMETRIC_LABELS_ATTR] = False
                f.create_dataset(cls._COL_LABELS_DATASET, data=col_names, dtype=UTF8_h5py_encoding)
                if col_lengths is not None:
                    f.create_dataset(cls._COL_LENGTHS_DATASET, data=col_lengths)

            f.attrs[cls._DATA_TYPE_ATTR] = data_type


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
            matrix (Union[np.array, scipy.sparse.sparray]): dense or sparse (dok) data matrix
            row_names:List[str]: row labels
            col_names:List[str]: column labels

        """

        row_names = list(row_names)
        col_names = list(col_names)

        if len(matrix.shape) != 2:
            raise ValueError(f"matrix must be 2-dimensional, not: {matrix.shape}")
        
        if len(row_names) != matrix.shape[0]:
            raise ValueError(f"row_names size must match dense_matrix dim 0, not: {len(row_names)} vs. {matrix.shape[0]}")

        if len(col_names) != matrix.shape[1]:
            raise ValueError(f"col_names size must match dense_matrix dim 1, not: {len(col_names)} vs. {matrix.shape[1]}")

        if len(set(row_names)) != len(row_names):
            raise ValueError("row_names must be unique")
        if len(set(col_names)) != len(col_names):
            raise ValueError("col_names must be unique")
        
        if row_lens is not None:
            if len(row_lens) != matrix.shape[0]:
                raise ValueError(f"row_lens size must match dense_matrix dim 0, not: {len(row_lens)} vs. {matrix.shape[0]}")
        if col_lens is not None:
            if len(col_lens) != matrix.shape[1]:
                raise ValueError(f"col_lens size must match dense_matrix dim 1, not: {len(col_lens)} vs. {matrix.shape[1]}")
        
        if row_lens is not None and col_lens is None:
            raise ValueError("If row_lens is provided, col_lens must also be provided.")
        if row_lens is None and col_lens is not None:
            raise ValueError("If col_lens is provided, row_lens must also be provided.")

        if row_names == col_names:
            col_names = row_names

        return row_names, col_names        


    @classmethod
    def write_dense(cls, matrix:Union[np.array, scipy.sparse.sparray], file_name: Union[PathLike, str], row_names:List[str], col_names:List[str], row_lengths:Optional[np.ndarray]=None, col_lengths:Optional[np.ndarray]=None, data_type:Optional[str]=""):
        """Writes a dense matrix to an hdf5 file.
        
        Args:
            matrix (Union[np.array, scipy.sparse.sparray]): 2d array
            file_name (Union[PathLike, str]): path to the target file to write the data to
            row_names:List[str]: row labels
            col_names:List[str]: column labels
        """
        # check validity of input data
        if scipy.sparse.issparse(matrix):
            matrix = matrix.toarray()

        row_names, col_names = cls.validate_matrix_data(matrix, row_names, col_names, row_lengths, col_lengths)
        

        with h5py.File(file_name, 'w') as f:
            f.attrs[cls._MATRIX_FILE_VERSION_ATTR] = cls._MATRIX_FILE_VERSION
            f.attrs[cls._ARRAY_TYPE_ATTR] = cls._DENSE
            f.create_dataset(cls._DENSE_DATASET, data=matrix)
            
            
            f.create_dataset(cls._ROW_LABELS_DATASET, data=row_names, dtype=UTF8_h5py_encoding)
            
            if row_lengths is not None:
                f.create_dataset(cls._ROW_LENGTHS_DATASET, data=row_lengths)


            if row_names is col_names and row_lengths is col_lengths:
                f.attrs[cls._SYMMETRIC_LABELS_ATTR] = True
            else:
                f.attrs[cls._SYMMETRIC_LABELS_ATTR] = False
                f.create_dataset(cls._COL_LABELS_DATASET, data=col_names, dtype=UTF8_h5py_encoding)
                if col_lengths is not None:
                    f.create_dataset(cls._COL_LENGTHS_DATASET, data=col_lengths)

            f.attrs[cls._DATA_TYPE_ATTR] = data_type

    @classmethod
    def write_dense_text(cls, matrix:Union[np.array, scipy.sparse.sparray], file_name: Union[PathLike, str], row_names:List[str], col_names:List[str], row_lengths:Optional[np.ndarray]=None, col_lengths:Optional[np.ndarray]=None, data_type:Optional[str]=""):
        """Writes a dense matrix to a text file.
        
        Args:
            matrix (Union[np.array, scipy.sparse.sparray]): 2d array
            file_name (Union[PathLike, str]): path to the target file to write the data to
            row_names:List[str]: row labels
            col_names:List[str]: column labels
        """

        row_names, col_names = cls.validate_matrix_data(matrix, row_names, col_names, row_lengths, col_lengths)
        outfile = open(file_name, "w")
        print("\t" + "\t".join(col_names), file=outfile)
        
        if scipy.sparse.issparse(matrix):
            for idx in range(len(row_names)):
                print(row_names[idx] + "\t" + "\t".join([str(x) for x in np.nditer(matrix[[idx],:].toarray())]), file=outfile)
        else:
            for idx in range(len(row_names)):
                print(row_names[idx] + "\t" + "\t".join([str(x) for x in np.nditer(matrix[idx])]), file=outfile)

        outfile.close()


class DenseDataMatrix(DataMatrix):
    """
    Dense matrix implementation using numpy arrays.
    
    This class stores matrix data as a numpy.ndarray, optimized for matrices
    with many non-zero values. For matrices with many zeros, consider using
    SparseDataMatrix instead.
    
    Attributes:
        data: numpy.ndarray containing the matrix values
        rows: List of row labels (inherited from DataMatrix)
        columns: List of column labels (inherited from DataMatrix)
        row_to_idx: Dictionary mapping row names to indices (inherited)
        column_to_idx: Dictionary mapping column names to indices (inherited)
        row_lengths: Optional numpy array of sequence lengths for rows (inherited)
        column_lengths: Optional numpy array of sequence lengths for columns (inherited)
        data_type: String describing data type (inherited)
    """
    
    def __init__(self, data: np.ndarray, row_names: List[str], col_names: List[str], 
                 row_lengths: Optional[np.ndarray] = None, col_lengths: Optional[np.ndarray] = None, 
                 data_type: str = None, *, copy_data: bool = True):
        """
        Initialize a dense data matrix.
        
        Args:
            data: 2D numpy array or array-like object. If sparse, will be converted
                to dense via toarray(). Dense numpy inputs are copied by default.
            row_names: List of row labels, length must match data.shape[0]
            col_names: List of column labels, length must match data.shape[1]
            row_lengths: Optional sequence lengths for rows (array-like, stored as numpy array)
            col_lengths: Optional sequence lengths for columns (array-like, stored as numpy array)
            data_type: Optional string describing the data type (e.g., 'score', 'norm_score')
            copy_data: If True, copy dense input arrays. If False, reuse dense numpy
                storage when possible. Sparse input conversion still creates a dense array.
        
        Raises:
            ValueError: If data dimensions don't match label lengths
        """
        # Validate before calling super().__init__
        DataMatrix.validate_matrix_data(data, row_names, col_names, row_lengths, col_lengths)
        super().__init__(row_names, col_names, row_lengths, col_lengths, data_type)
        
        if scipy.sparse.issparse(data):
            self.data = data.toarray()
        elif copy_data:
            self.data = np.array(data, copy=True)
        else:
            self.data = np.asarray(data)
    
    @property
    def shape(self):
        return self.data.shape
    
    @property
    def symmetric_values(self) -> bool:
        """Check if matrix values are symmetric (np.isclose tolerance)."""
        if not self.symmetric_labels:
            return False

        return np.allclose(self.data, self.data.T)
    
    @classmethod
    def _from_hdf5(cls, matrix_file: Union[PathLike, str], lower_bound: Optional[float] = None) -> 'DenseDataMatrix':
        """Load a dense matrix from HDF5 file."""
        with h5py.File(matrix_file, 'r') as f:
            # Validate file format
            cls._validate_hdf5_attributes(f)

            if f.attrs[cls._ARRAY_TYPE_ATTR] != cls._DENSE:
                raise ValueError(f"Expected DENSE matrix, got {f.attrs[cls._ARRAY_TYPE_ATTR]}")

            # Load data
            data = f[cls._DENSE_DATASET][:]
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

            if lower_bound is not None:
                # NaNs are passed through unchanged (NaN <= lower_bound is False).
                data[data <= lower_bound] = 0

            return cls(data, row_names, col_names, row_lengths, col_lengths, data_type, copy_data=False)
    
    @classmethod
    def _from_text(cls, matrix_file: Union[PathLike, str]) -> 'DenseDataMatrix':
        """Load a dense matrix from text file."""
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

    def iter_data(self, index_style="name", skip_zeros=True):
        """Iterate over all values in the dense matrix."""
        index_style = index_style.lower()
        if index_style not in {"name", "index"}:
            raise ValueError("index_style must be 'name' or 'index'.")

        use_names = index_style == "name"
        for r in range(self.data.shape[0]):
            row = self.data[r]
            cols = np.flatnonzero(row) if skip_zeros else range(self.data.shape[1])
            for c in cols:
                c = int(c)
                value = row[c]
                if use_names:
                    yield self.rows[r], self.columns[c], value
                else:
                    yield r, c, value

    def triangular_iter(self, side="lower", include_diagonal=False, skip_zeros=False, index_style="name", agg=None):
        """Iterate over the triangular part of the matrix."""
        if self.data.shape[0] != self.data.shape[1]:
            raise NotImplementedError(f"Triangular only supported for square matrices")
        
        side = side.lower()
        if side not in {"upper", "lower"}:
            raise ValueError("side must be 'lower' or 'upper'")
        
        index_style = index_style.lower()
        if index_style not in {"name", "index"}:
            raise ValueError("index_style must be 'name' or 'index'.")
        
        n_cols = self.data.shape[1]
        for r in range(self.data.shape[0]):
            start_c, end_c = _triangle_col_bounds(r, side, include_diagonal, n_cols)

            for c in range(start_c, end_c):
                value = self.data[r, c]

                if agg is not None:
                    symmetric_value = self.data[c, r]
                    if skip_zeros and value == 0 and symmetric_value == 0:
                        continue
                    value = agg(value, symmetric_value)
                elif skip_zeros and value == 0:
                    continue

                if index_style == "name":
                    yield self.rows[r], self.columns[c], value
                else:
                    yield r, c, value

    def sorted_undirected_edges(self, skip_zeros: bool = False, agg=None) -> SortedUndirectedEdges:
        """Return lower-triangular undirected edges sorted by descending score."""
        if self.data.shape[0] != self.data.shape[1]:
            raise NotImplementedError("sorted_undirected_edges only supported for square matrices")

        # Fast path is an intentional optimization keyed on the builtin `max` object.
        # Any other agg callable (including an equivalent `lambda a, b: max(a, b)`) is
        # still handled correctly, just via the generic element-wise base-class path.
        if agg is not None and agg is not max:
            return super().sorted_undirected_edges(skip_zeros=skip_zeros, agg=agg)

        n = self.data.shape[0]
        # Pack source/target/score into a single structured array so sorting is
        # in-place and avoids a separate O(n²/2) argsort index array.
        dtype = np.dtype([('score', np.float64), ('source', np.int32), ('target', np.int32)])
        edges = np.empty(n * (n - 1) // 2, dtype=dtype)

        row_idx, col_idx = np.tril_indices(n, k=-1)
        if agg is max:
            edges['score'] = np.maximum(self.data[row_idx, col_idx], self.data[col_idx, row_idx])
        else:
            edges['score'] = self.data[row_idx, col_idx]
        edges['source'] = row_idx
        edges['target'] = col_idx
        del row_idx, col_idx  # free before sort

        if skip_zeros:
            edges = edges[edges['score'] != 0]

        if edges.size == 0:
            return SortedUndirectedEdges(
                n_nodes=len(self),
                source=np.empty(0, dtype=np.int32),
                target=np.empty(0, dtype=np.int32),
                score=np.empty(0, dtype=float),
            )

        edges.sort(order='score', kind='stable')  # ascending in-place
        edges = edges[::-1]  # descending view; base array kept alive via field views below
        return SortedUndirectedEdges(
            n_nodes=len(self),
            source=edges['source'],
            target=edges['target'],
            score=edges['score'],
        )
    
    def itervalues(self) -> Iterator[float]:
        """Yield all values row by row."""
        yield from self.data.flat

    def convert_to_sparse(self) -> 'SparseDataMatrix':
        """Convert to sparse representation."""
        return SparseDataMatrix(
            scipy.sparse.csr_array(self.data),
            self.rows, self.columns,
            self.row_lengths, self.column_lengths,
            self.data_type
        )
    
    def convert_to_dense(self) -> 'DenseDataMatrix':
        """Return self (already dense)."""
        return self
    
    def toarray(self) -> np.ndarray:
        """Return as numpy array."""
        return self.data


class SparseDataMatrix(DataMatrix):
    """
    Sparse matrix implementation using scipy CSR (Compressed Sparse Row) arrays.
    
    This class stores matrix data as a scipy.sparse.csr_array, optimized for matrices
    with many zero values. The sparse storage format saves memory and accelerates
    operations that skip zeros. For dense matrices (many non-zeros), use DenseDataMatrix.
    
    The CSR format is particularly efficient for row-based operations and arithmetic.
    Zeros are automatically eliminated and indices are sorted during initialization.
    
    Attributes:
        data: scipy.sparse.csr_array containing the sparse matrix
        rows: List of row labels (inherited from DataMatrix)
        columns: List of column labels (inherited from DataMatrix)
        row_to_idx: Dictionary mapping row names to indices (inherited)
        column_to_idx: Dictionary mapping column names to indices (inherited)
        row_lengths: Optional numpy array of sequence lengths for rows (inherited)
        column_lengths: Optional numpy array of sequence lengths for columns (inherited)
        data_type: String describing data type (inherited)
    """
    
    def __init__(self, data: scipy.sparse.csr_array, row_names: List[str], col_names: List[str],
                 row_lengths: Optional[np.ndarray] = None, col_lengths: Optional[np.ndarray] = None,
                 data_type: str = None, *, copy_data: bool = True):
        """
        Initialize a sparse data matrix.
        
        The input data is converted to CSR format if needed, and zeros are eliminated
        with indices sorted for efficiency.
        
        Args:
            data: Scipy sparse matrix or array-like object. Will be converted to 
                  scipy.sparse.csr_array format and copied. Dense arrays will be converted to sparse.
            row_names: List of row labels, length must match data.shape[0]
            col_names: List of column labels, length must match data.shape[1]
            row_lengths: Optional sequence lengths for rows (array-like, stored as numpy array)
            col_lengths: Optional sequence lengths for columns (array-like, stored as numpy array)
            data_type: Optional string describing the data type (e.g., 'score', 'norm_score')
        
        Raises:
            ValueError: If data dimensions don't match label lengths
        """
        # Validate before calling super().__init__
        DataMatrix.validate_matrix_data(data, row_names, col_names, row_lengths, col_lengths)
        super().__init__(row_names, col_names, row_lengths, col_lengths, data_type)
        
        if not scipy.sparse.issparse(data):
            self.data = scipy.sparse.csr_array(data)
        elif not isinstance(data, scipy.sparse.csr_array):
            self.data = scipy.sparse.csr_array(data)
        else:
            self.data = data.copy() if copy_data else data
        
        self.data.eliminate_zeros()
        self.data.sort_indices()
    
    @property
    def shape(self) -> Tuple[int, int]:
        """Return shape of the matrix as (rows, columns)."""
        return self.data.shape
    
    @property
    def nnz(self) -> int:
        """Return number of non-zero elements (efficient for sparse matrices)."""
        return self.data.nnz
    
    @property
    def symmetric_values(self) -> bool:
        """Check if matrix values are symmetric (optimized for sparse)."""
        if not self.symmetric_labels:
            return False

        # Compare the CSR structure and values against the transpose. Zeros are
        # eliminated, so structural symmetry of the non-zero pattern is required,
        # then values are compared element-wise with np.isclose tolerance.
        A = self.data
        A.eliminate_zeros()
        A.sort_indices()
        AT = A.transpose().tocsr()
        AT.sort_indices()

        if A.nnz != AT.nnz:
            return False
        if not np.array_equal(A.indptr, AT.indptr) or not np.array_equal(A.indices, AT.indices):
            return False
        return np.allclose(A.data, AT.data)
    
    @classmethod
    def _from_hdf5(cls, matrix_file: Union[PathLike, str], lower_bound: Optional[float] = None) -> 'SparseDataMatrix':
        """Load a sparse matrix from HDF5 file."""
        with h5py.File(matrix_file, 'r') as f:
            # Validate file format
            cls._validate_hdf5_attributes(f)
            
            if f.attrs[cls._ARRAY_TYPE_ATTR] != cls._SPARSE_CSR:
                raise ValueError(f"Expected SPARSE_CSR matrix, got {f.attrs[cls._ARRAY_TYPE_ATTR]}")
            
            # Load data
            row_names, _ = cls.parse_axis_index(f[cls._ROW_LABELS_DATASET])
            
            if f.attrs[cls._SYMMETRIC_LABELS_ATTR]:
                col_names = row_names
            else:
                col_names, _ = cls.parse_axis_index(f[cls._COL_LABELS_DATASET])
            
            data = scipy.sparse.csr_array(
                (f[cls._SPARSE_VALUES_DATASET][:], 
                 f[cls._SPARSE_CSR_INDICES_DATASET][:], 
                 f[cls._SPARSE_CSR_INDPTR_DATASET][:]), 
                (len(row_names), len(col_names))
            )
            if lower_bound is not None:
                # NaNs are passed through unchanged (NaN <= lower_bound is False).
                data.data[data.data <= lower_bound] = 0
                data.eliminate_zeros()
                data.sort_indices()
            
            row_lengths = f[cls._ROW_LENGTHS_DATASET][:] if cls._ROW_LENGTHS_DATASET in f else None
            col_lengths = f[cls._COL_LENGTHS_DATASET][:] if cls._COL_LENGTHS_DATASET in f else None
            if f.attrs[cls._SYMMETRIC_LABELS_ATTR] and row_lengths is not None and col_lengths is None:
                col_lengths = row_lengths
            data_type = f.attrs.get(cls._DATA_TYPE_ATTR, "")
            
            return cls(data, row_names, col_names, row_lengths, col_lengths, data_type, copy_data=False)

    def iter_data(self, index_style="name", skip_zeros=True):
        """Iterate over values in the sparse matrix."""
        index_style = index_style.lower()
        if index_style not in {"name", "index"}:
            raise ValueError("index_style must be 'name' or 'index'.")
        
        if skip_zeros:
            # Iterate only over non-zero values (efficient for sparse)
            self.data.eliminate_zeros()
            self.data.sort_indices()
            
            for r in range(self.data.shape[0]):
                for ind in range(self.data.indptr[r], self.data.indptr[r+1]):
                    c = self.data.indices[ind]
                    if index_style == "name":
                        yield self.rows[r], self.columns[c], self.data.data[ind]
                    else:
                        yield r, c, self.data.data[ind]
        else:
            # Have to iterate over entire matrix including zeros. Densify one row at a
            # time to avoid per-element O(log nnz) CSR scalar indexing.
            use_names = index_style == "name"
            n_cols = self.data.shape[1]
            for r in range(self.data.shape[0]):
                row = self.data[[r]].toarray().ravel()
                for c in range(n_cols):
                    value = row[c]
                    if use_names:
                        yield self.rows[r], self.columns[c], value
                    else:
                        yield r, c, value
    
    def triangular_iter(self, side="lower", include_diagonal=False, skip_zeros=False, index_style="name", agg=None):
        """Iterate over the triangular part of the matrix (optimized for sparse)."""
        if self.data.shape[0] != self.data.shape[1]:
            raise NotImplementedError(f"Triangular only supported for square matrices")
        
        side = side.lower()
        if side not in {"upper", "lower"}:
            raise ValueError("side must be 'lower' or 'upper'")
        
        index_style = index_style.lower()
        if index_style not in {"name", "index"}:
            raise ValueError("index_style must be 'name' or 'index'.")
        
        use_names = index_style == "name"
        n_cols = self.data.shape[1]

        if agg is max and skip_zeros:
            # Fast path keyed on the builtin `max`: symmetrize via element-wise maximum
            # (absent entries are treated as 0, matching max(value, symmetric_value)), then
            # walk the requested triangle's stored nonzeros in row-major, column order.
            symmetrized = self.data.maximum(self.data.transpose()).tocsr()
            symmetrized.eliminate_zeros()
            symmetrized.sort_indices()

            for r in range(symmetrized.shape[0]):
                start_c, end_c = _triangle_col_bounds(r, side, include_diagonal, n_cols)
                for ind in range(symmetrized.indptr[r], symmetrized.indptr[r + 1]):
                    c = int(symmetrized.indices[ind])
                    if start_c <= c < end_c:
                        value = symmetrized.data[ind]
                        if use_names:
                            yield self.rows[r], self.columns[c], value
                        else:
                            yield r, c, value

        elif agg is not None and skip_zeros:
            # Generic aggregation: visit every position where (r,c) or (c,r) is non-zero.
            self.data.eliminate_zeros()
            self.data.sort_indices()

            positions_to_process = set()
            for r in range(self.data.shape[0]):
                start_c, end_c = _triangle_col_bounds(r, side, include_diagonal, n_cols)
                for ind in range(self.data.indptr[r], self.data.indptr[r + 1]):
                    c = int(self.data.indices[ind])
                    if start_c <= c < end_c:
                        positions_to_process.add((r, c))
                    # The transpose entry (c, r) may also fall in the triangle.
                    comp_start_c, comp_end_c = _triangle_col_bounds(c, side, include_diagonal, n_cols)
                    if comp_start_c <= r < comp_end_c:
                        positions_to_process.add((c, r))

            for r, c in sorted(positions_to_process):
                value = self.data[r, c]
                symmetric_value = self.data[c, r]
                if value == 0 and symmetric_value == 0:
                    continue
                agg_value = agg(value, symmetric_value)
                if use_names:
                    yield self.rows[r], self.columns[c], agg_value
                else:
                    yield r, c, agg_value

        elif skip_zeros:
            # Sparse without agg - only iterate stored non-zeros.
            self.data.eliminate_zeros()
            self.data.sort_indices()

            for r in range(self.data.shape[0]):
                start_c, end_c = _triangle_col_bounds(r, side, include_diagonal, n_cols)
                for ind in range(self.data.indptr[r], self.data.indptr[r + 1]):
                    c = int(self.data.indices[ind])
                    if start_c <= c < end_c:
                        value = self.data.data[ind]
                        if use_names:
                            yield self.rows[r], self.columns[c], value
                        else:
                            yield r, c, value

        else:
            # Need all values including zeros - densify one row at a time.
            for r in range(self.data.shape[0]):
                start_c, end_c = _triangle_col_bounds(r, side, include_diagonal, n_cols)
                if start_c >= end_c:
                    continue
                row = self.data[[r]].toarray().ravel()
                if agg is not None:
                    transpose_col = self.data[:, [r]].toarray().ravel()
                for c in range(start_c, end_c):
                    value = row[c]
                    if agg is not None:
                        value = agg(value, transpose_col[c])
                    if use_names:
                        yield self.rows[r], self.columns[c], value
                    else:
                        yield r, c, value

    def sorted_undirected_edges(self, skip_zeros: bool = False, agg=None) -> SortedUndirectedEdges:
        """Return lower-triangular undirected edges sorted by descending score."""
        if self.data.shape[0] != self.data.shape[1]:
            raise NotImplementedError("sorted_undirected_edges only supported for square matrices")

        # Fast path keyed on the builtin `max` (or no aggregation). For agg=max we
        # symmetrize with an element-wise maximum (absent entries are treated as 0,
        # matching max(value, symmetric_value)); any other agg falls back to the
        # generic element-wise base-class implementation.
        if skip_zeros and (agg is None or agg is max):
            self.data.eliminate_zeros()
            self.data.sort_indices()

            working = self.data if agg is None else self.data.maximum(self.data.transpose())
            lower = scipy.sparse.tril(working, k=-1).tocsr()
            lower.eliminate_zeros()
            lower.sort_indices()
            coo = lower.tocoo()  # row-major, column-sorted (matches dense tril_indices order)

            if coo.nnz == 0:
                return SortedUndirectedEdges(
                    n_nodes=len(self),
                    source=np.empty(0, dtype=np.int32),
                    target=np.empty(0, dtype=np.int32),
                    score=np.empty(0, dtype=float),
                )

            # Sort exactly as DenseDataMatrix.sorted_undirected_edges does (stable
            # ascending then reversed) so dense and sparse outputs are byte-identical.
            dtype = np.dtype([('score', np.float64), ('source', np.int32), ('target', np.int32)])
            edges = np.empty(coo.nnz, dtype=dtype)
            edges['score'] = coo.data
            edges['source'] = coo.row
            edges['target'] = coo.col
            edges.sort(order='score', kind='stable')
            edges = edges[::-1]
            return SortedUndirectedEdges(
                n_nodes=len(self),
                source=edges['source'],
                target=edges['target'],
                score=edges['score'],
            )

        return super().sorted_undirected_edges(skip_zeros=skip_zeros, agg=agg)
    
    def itervalues(self) -> Iterator[float]:
        """Yield all values row by row."""
        # Densify one row at a time to avoid per-element CSR scalar indexing.
        for r in range(self.shape[0]):
            yield from self.data[[r]].toarray().ravel()

    def convert_to_sparse(self) -> 'SparseDataMatrix':
        """Return self (already sparse)."""
        return self
    
    def convert_to_dense(self) -> 'DenseDataMatrix':
        """Convert to dense representation."""
        return DenseDataMatrix(
            self.data.toarray(),
            self.rows, self.columns,
            self.row_lengths, self.column_lengths,
            self.data_type,
            copy_data=False,
        )
    
    def toarray(self) -> np.ndarray:
        """Return as numpy array."""
        return self.data.toarray()


# Backward-compatible re-exports. MaxTree and the edge/neighbor-ranking functions live
# in ssn_edges; importing here (after DataMatrix is defined) keeps existing
# `from domainator.data_matrix import MaxTree, ...` call sites working and avoids an
# import cycle (ssn_edges imports DataMatrix lazily).
from .ssn_edges import (  # noqa: E402
    MaxTree,
    build_symmetric_neighbor_rankings,
    symmetric_knn_edge_index_dict,
    mst_edge_index_dict,
    mst_knn_edge_index_dict,
    sorted_edges_from_edge_index_dict,
    mst_knn_edge_counts_by_threshold,
)
