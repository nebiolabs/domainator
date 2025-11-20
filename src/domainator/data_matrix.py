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
from os import PathLike
from pathlib import Path
import pandas as pd
import scipy.sparse
import numpy as np
import h5py
from abc import ABC, abstractmethod
from typing import Dict, List, Tuple, Union, Optional, Iterator

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
    """row index for CSR sparse"""
    _SPARSE_CSR_INDPTR_DATASET="SPARSE_CSR_INDPTR"
    """col index for CSR sparse"""
    
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
        return np.prod(self.shape)

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
    def from_file(cls, matrix_file):
        """
        Factory method that returns the appropriate subclass (DenseDataMatrix or SparseDataMatrix)
        based on the file content.
        """
        file_type = _get_file_type(matrix_file)
        if file_type == "hdf5":
            # Peek at the file to determine type
            with h5py.File(matrix_file, 'r') as f:
                array_type = f.attrs.get(cls._ARRAY_TYPE_ATTR, cls._DENSE)
                
                if array_type == cls._SPARSE_CSR:
                    return SparseDataMatrix._from_hdf5(matrix_file)
                else:
                    return DenseDataMatrix._from_hdf5(matrix_file)
        elif file_type is None:
            # Assume dense tsv text file if extension not recognized
            return DenseDataMatrix._from_text(matrix_file)
        else:
            raise ValueError(f"Unrecognized file type: {file_type}")

    @abstractmethod
    def _read_hdf5(self, matrix_file: Union[PathLike, str]):
        """
        Read matrix from HDF5 file - implemented by subclasses.
        
        This method should not be called directly. Use the factory method
        DataMatrix.from_file() instead, which returns the appropriate subclass.
        
        Args:
            matrix_file: Path to HDF5 file
        """
        pass

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
    def triangular(self, side="lower", include_diagonal=False, skip_zeros=False, index_style="name", agg=None):
        """
        Return triangular part of the matrix as a list of tuples.
        
        Only works for square matrices. Returns values from either upper or lower triangle.
        
        Args:
            side: "lower" or "upper" triangle
            include_diagonal: If True, include diagonal elements
            skip_zeros: If True, exclude zero values from results
            index_style: Either "name" (return row/column names) or "index" (return integer indices)
            agg: Optional aggregation function to combine symmetric positions (e.g., lambda a,b: max(a,b)).
                 When provided, applies agg(value, transpose_value) for each position.
        
        Returns:
            list: List of tuples (row_identifier, col_identifier, value)
        
        Raises:
            NotImplementedError: If matrix is not square
        """
        pass

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

    @staticmethod
    def parse_axis_index(axis_index):
        items_list = [x.decode("utf-8") for x in axis_index]
        items_dict = {v:i for i,v in enumerate(items_list)}
        
        return items_list, items_dict


    @classmethod
    def write_sparse(cls, matrix, file_name, row_names, col_names, row_lengths:Optional[np.ndarray]=None, col_lengths:Optional[np.ndarray]=None, data_type=""):
        """Writes a dok sparse matrix to an hdf5 file
        
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

        if len(matrix.shape) != 2:
            raise ValueError(f"matrix must be 2-dimensional, not: {matrix.shape}")
        
        if len(row_names) != matrix.shape[0]:
            raise ValueError(f"row_names size must match dense_matrix dim 0, not: {len(row_names)} vs. {matrix.shape[0]}")

        if len(col_names) != matrix.shape[1]:
            raise ValueError(f"col_names size must match dense_matrix dim 1, not: {len(col_names)} vs. {matrix.shape[1]}")
        
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



        if (row_names == col_names) or (all(a == b for a,b in zip(row_names, col_names))):
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

        row_names, col_names = cls.validate_matrix_data(matrix, row_names, col_names)
        

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

        row_names, col_names = cls.validate_matrix_data(matrix, row_names, col_names)
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
                 data_type: str = None):
        """
        Initialize a dense data matrix.
        
        Args:
            data: 2D numpy array or array-like object (will be copied). If sparse, 
                  will be converted to dense via toarray()
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
        
        if scipy.sparse.issparse(data):
            self.data = data.toarray()
        else:
            self.data = np.array(data, copy=True)
    
    @property
    def shape(self):
        return self.data.shape
    
    @property
    def symmetric_values(self) -> bool:
        """Check if matrix values are symmetric."""
        if not self.symmetric_labels:
            return False
        
        # Check upper triangular against lower triangular
        upper_triangular = self.triangular(side='upper', include_diagonal=False, skip_zeros=False, index_style='index')
        
        for row_idx, col_idx, value in upper_triangular:
            transpose_value = self.data[col_idx, row_idx]
            if not np.isclose(value, transpose_value):
                return False
        
        return True
    
    @classmethod
    def _from_hdf5(cls, matrix_file: Union[PathLike, str]) -> 'DenseDataMatrix':
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
            data_type = f.attrs.get(cls._DATA_TYPE_ATTR, "")
            
            return cls(data, row_names, col_names, row_lengths, col_lengths, data_type)
    
    @classmethod
    def _from_text(cls, matrix_file: Union[PathLike, str]) -> 'DenseDataMatrix':
        """Load a dense matrix from text file."""
        try:
            mat_file = pd.read_csv(matrix_file, sep="\t", index_col=0)
        except UnicodeDecodeError:
            raise ValueError(f"Expecting utf-8 encoded matrix file, but found invalid unicode byte.")
        
        row_names = list(mat_file.index)
        col_names = list(mat_file.columns)
        data = mat_file.to_numpy()
        
        return cls(data, row_names, col_names)
    
    def _read_hdf5(self, matrix_file: Union[PathLike, str]):
        """Not used - handled by _from_hdf5 classmethod"""
        raise NotImplementedError("Use DenseDataMatrix._from_hdf5() instead")
    
    def iter_data(self, index_style="name", skip_zeros=True):
        """Iterate over all values in the dense matrix."""
        index_style = index_style.lower()
        if index_style not in {"name", "index"}:
            raise ValueError("index_style must be 'name' or 'index'.")
        
        for r in range(self.data.shape[0]):
            for c in range(self.data.shape[1]):
                value = self.data[r, c]
                if not skip_zeros or value != 0:
                    if index_style == "name":
                        yield self.rows[r], self.columns[c], value
                    else:
                        yield r, c, value
    
    def triangular(self, side="lower", include_diagonal=False, skip_zeros=False, index_style="name", agg=None):
        """Get triangular part of the matrix."""
        if self.data.shape[0] != self.data.shape[1]:
            raise NotImplementedError(f"Triangular only supported for square matrices")
        
        side = side.lower()
        if side not in {"upper", "lower"}:
            raise ValueError("side must be 'lower' or 'upper'")
        
        index_style = index_style.lower()
        if index_style not in {"name", "index"}:
            raise ValueError("index_style must be 'name' or 'index'.")
        
        out = []
        
        for r in range(self.data.shape[0]):
            if side == "lower":
                start_c = 0
                end_c = r + 1 if include_diagonal else r
            else:  # upper
                start_c = r if include_diagonal else r + 1
                end_c = self.data.shape[1]
            
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
                    out.append((self.rows[r], self.columns[c], value))
                else:
                    out.append((r, c, value))
        
        return out
    
    def itervalues(self) -> Iterator[float]:
        """Yield all values row by row."""
        for r in range(self.shape[0]):
            for c in range(self.shape[1]):
                yield self.data[r, c]
    
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
                 data_type: str = None):
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
            self.data = data.copy()
        
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
        
        # Get upper and lower triangular non-zero elements
        upper_triangular = self.triangular(side='upper', include_diagonal=False, skip_zeros=True, index_style='index')
        lower_triangular = self.triangular(side='lower', include_diagonal=False, skip_zeros=True, index_style='index')
        
        if len(upper_triangular) != len(lower_triangular):
            return False
        
        # Create transposed lower triangular for comparison
        lower_transposed = [(l_c, l_r, l_v) for l_r, l_c, l_v in lower_triangular]
        lower_transposed.sort()
        
        # Compare element by element
        for (u_r, u_c, u_v), (lt_r, lt_c, lt_v) in zip(upper_triangular, lower_transposed):
            if (u_r, u_c) != (lt_r, lt_c):
                return False
            if not np.isclose(u_v, lt_v):
                return False
        
        return True
    
    @classmethod
    def _from_hdf5(cls, matrix_file: Union[PathLike, str]) -> 'SparseDataMatrix':
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
            
            row_lengths = f[cls._ROW_LENGTHS_DATASET][:] if cls._ROW_LENGTHS_DATASET in f else None
            col_lengths = f[cls._COL_LENGTHS_DATASET][:] if cls._COL_LENGTHS_DATASET in f else None
            data_type = f.attrs.get(cls._DATA_TYPE_ATTR, "")
            
            return cls(data, row_names, col_names, row_lengths, col_lengths, data_type)
    
    def _read_hdf5(self, matrix_file: Union[PathLike, str]):
        """Not used - handled by _from_hdf5 classmethod"""
        raise NotImplementedError("Use SparseDataMatrix._from_hdf5() instead")
    
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
            # Have to iterate over entire matrix including zeros
            for r in range(self.data.shape[0]):
                for c in range(self.data.shape[1]):
                    value = self.data[r, c]
                    if index_style == "name":
                        yield self.rows[r], self.columns[c], value
                    else:
                        yield r, c, value
    
    def triangular(self, side="lower", include_diagonal=False, skip_zeros=False, index_style="name", agg=None):
        """Get triangular part of the matrix (optimized for sparse)."""
        if self.data.shape[0] != self.data.shape[1]:
            raise NotImplementedError(f"Triangular only supported for square matrices")
        
        side = side.lower()
        if side not in {"upper", "lower"}:
            raise ValueError("side must be 'lower' or 'upper'")
        
        index_style = index_style.lower()
        if index_style not in {"name", "index"}:
            raise ValueError("index_style must be 'name' or 'index'.")
        
        out = []
        
        if agg is not None and skip_zeros:
            # Optimized path for sparse with aggregation
            self.data.eliminate_zeros()
            self.data.sort_indices()
            
            positions_to_process = set()
            
            # Collect positions where at least one of (r,c) or (c,r) is non-zero
            for r in range(self.data.shape[0]):
                for ind in range(self.data.indptr[r], self.data.indptr[r+1]):
                    c = self.data.indices[ind]
                    
                    if side == "lower":
                        start_c, end_c = 0, r + 1 if include_diagonal else r
                    else:
                        start_c, end_c = r if include_diagonal else r + 1, self.data.shape[1]
                    
                    if start_c <= c < end_c:
                        positions_to_process.add((r, c))
                    
                    # Check complement
                    if side == "lower":
                        comp_start_c, comp_end_c = 0, c + 1 if include_diagonal else c
                    else:
                        comp_start_c, comp_end_c = c if include_diagonal else c + 1, self.data.shape[1]
                    
                    if comp_start_c <= r < comp_end_c:
                        positions_to_process.add((c, r))
            
            for r, c in sorted(positions_to_process):
                value = self.data[r, c]
                symmetric_value = self.data[c, r]
                
                if skip_zeros and value == 0 and symmetric_value == 0:
                    continue
                
                agg_value = agg(value, symmetric_value)
                
                if index_style == "name":
                    out.append((self.rows[r], self.columns[c], agg_value))
                else:
                    out.append((r, c, agg_value))
        
        elif skip_zeros:
            # Sparse without agg - only iterate non-zeros
            self.data.eliminate_zeros()
            self.data.sort_indices()
            
            for r in range(self.data.shape[0]):
                if side == "lower":
                    start_c, end_c = 0, r + 1 if include_diagonal else r
                else:
                    start_c, end_c = r if include_diagonal else r + 1, self.data.shape[1]
                
                for ind in range(self.data.indptr[r], self.data.indptr[r+1]):
                    c = self.data.indices[ind]
                    
                    if start_c <= c < end_c:
                        value = self.data.data[ind]
                        
                        if index_style == "name":
                            out.append((self.rows[r], self.columns[c], value))
                        else:
                            out.append((r, c, value))
        
        else:
            # Need all values including zeros - use dense logic
            for r in range(self.data.shape[0]):
                if side == "lower":
                    start_c, end_c = 0, r + 1 if include_diagonal else r
                else:
                    start_c, end_c = r if include_diagonal else r + 1, self.data.shape[1]
                
                for c in range(start_c, end_c):
                    value = self.data[r, c]
                    
                    if agg is not None:
                        symmetric_value = self.data[c, r]
                        value = agg(value, symmetric_value)
                    
                    if index_style == "name":
                        out.append((self.rows[r], self.columns[c], value))
                    else:
                        out.append((r, c, value))
        
        return out
    
    def itervalues(self) -> Iterator[float]:
        """Yield all values row by row."""
        for r in range(self.shape[0]):
            for c in range(self.shape[1]):
                yield self.data[r, c]
    
    def convert_to_sparse(self) -> 'SparseDataMatrix':
        """Return self (already sparse)."""
        return self
    
    def convert_to_dense(self) -> 'DenseDataMatrix':
        """Convert to dense representation."""
        return DenseDataMatrix(
            self.data.toarray(),
            self.rows, self.columns,
            self.row_lengths, self.column_lengths,
            self.data_type
        )
    
    def toarray(self) -> np.ndarray:
        """Return as numpy array."""
        return self.data.toarray()


class MaxTree():
    """Maximum Spanning Tree (MST) from a DataMatrix.
    
    A Maximum Spanning Tree is a subgraph that connects all nodes with the maximum
    total edge weight while avoiding cycles. This is useful for clustering and
    visualization of similarity/distance matrices.
    
    The MST is computed using Kruskal's algorithm with union-find for cycle detection.
    Edges are processed from highest to lowest weight, and only edges connecting
    different components are added to the tree.
    
    Use Cases:
        - Hierarchical clustering: Cut tree at different thresholds to create clusters
        - Network visualization: Display most significant relationships
        - Sequence similarity networks: Connect most similar sequences
    
    Args:
        matrix: DataMatrix containing edge weights (typically similarity or distance scores)
        skip_zeros: If True, zero-weight edges are excluded. When there is no connected
                   MST without zero edges, the result will be a forest (multiple trees).
    
    Attributes:
        n_nodes: Number of nodes in the graph
        edges: Array of shape (n_edges, 4) containing [node_i, node_j, weight, gap_count]
               where gap_count is the number of edges between consecutive MST edges
        mst: Array of indices indicating which edges are in the MST
        edges_by_threshold: Array of (edge_count, threshold) pairs for MST edges
        cluster_count_by_threshold: Array of (threshold, cluster_count) showing how
                                   clusters form at different thresholds
        cluster_count_by_edge_count: Array of (edge_count, cluster_count) including
                                    non-MST edges between MST edges
    
    Properties:
        mst_edges: List of (node_i, node_j, weight) tuples for MST edges only
    
    Methods:
        export_for_interactive_viz: Export minimal structure for client-side clustering
    
    Example:
        >>> matrix = DataMatrix.from_file("similarities.hdf5")
        >>> tree = MaxTree(matrix, skip_zeros=True)
        >>> # Get clusters at a specific threshold
        >>> mst_edges = tree.mst_edges
        >>> # Find how many clusters at different thresholds
        >>> threshold_clusters = tree.cluster_count_by_threshold
    
    Note:
        When skip_zeros=True and the matrix is disconnected (contains at least one pair of nodes with no path),
        the resulting MST will be a forest of multiple disconnected trees, one per
        connected component.
    """

    def __init__(self, matrix:DataMatrix, skip_zeros=True):
        triangular = matrix.triangular(agg=max, index_style="index", skip_zeros=skip_zeros) #list of tuples: i, j, value 

       
        
        self.n_nodes = len(matrix)
        self.edges = np.zeros((len(triangular), 4)) # node_i, node_j, edge_value, edges between this and next mst_edge (if 0, then this is not an mst_edge)

        # Handle empty case (single node or empty matrix)
        if len(triangular) > 0:
            self.edges[:,0:3] = triangular

        del triangular
        np.take(self.edges, np.argsort(self.edges[:, 2])[::-1], axis=0, out=self.edges) # sort highest to lowest by value

        # calculate_max_spanning_tree using union-find
        self.mst = list() #np.zeros((max(0, self.n_nodes - 1),), dtype=int) # indexes of edges in the mst
        parent = np.arange(self.n_nodes, dtype=int)  # Union-find parent array
        
        def find(x):
            """Find root with path compression"""
            if parent[x] != x:
                parent[x] = find(parent[x])
            return parent[x]
        
        def union(x, y):
            """Union two sets"""
            root_x = find(x)
            root_y = find(y)
            if root_x != root_y:
                parent[root_x] = root_y
                return True
            return False
        
        edge_count = 0
        #mst_i = 0
        for edge_i in range(len(self.edges)):
            node_1 = int(self.edges[edge_i, 0])
            node_2 = int(self.edges[edge_i, 1])
            edge_count += 1
            
            # Only add edge if it connects different components
            if union(node_1, node_2):
                self.edges[edge_i, 3] = edge_count
                self.mst.append(edge_i)
                #mst_i += 1
                edge_count = 0
                
                # if mst_i == self.n_nodes - 1:
                #     break  # MST is complete
        self.mst = np.array(self.mst)

        self.edges_by_threshold = self._edges_by_threshold()
        self.cluster_count_by_threshold = self._cluster_count_by_threshold()
        self.cluster_count_by_edge_count = self._cluster_count_by_edge_count()

    @property
    def mst_edges(self) -> List[Tuple[int,int,float]]: # (node_i, node_j, edge_value)
        """
        Returns the edges in the maximum spanning tree.
        
        Returns:
            List[Tuple[int,int,float]]: List of MST edges as (node_i, node_j, edge_value) tuples
        """
        mst_edge_list = []
        for mst_idx in self.mst:
            edge = self.edges[mst_idx]
            mst_edge_list.append((int(edge[0]), int(edge[1]), float(edge[2])))
        return mst_edge_list
        

    def _edges_by_threshold(self):
        """
        Returns array of (edge_count, threshold) pairs.
        Each row represents the cumulative edge count and threshold at each MST edge.
        """
        edges_by_threshold = np.zeros((len(self.mst), 2)) # edge_count, threshold
        edges_count = 0
        for i in range(len(self.mst)):
            mst_edge_idx = self.mst[i]
            edges_count += int(self.edges[mst_edge_idx, 3])
            threshold = self.edges[mst_edge_idx, 2]
            edges_by_threshold[i, 0] = edges_count
            edges_by_threshold[i, 1] = threshold
        return edges_by_threshold

    def _cluster_count_by_threshold(self):
        """
        Returns array of (threshold, cluster_count) pairs.
        For each MST edge threshold, calculates how many clusters would result
        if we only include edges with values >= that threshold.
        """
        cluster_counts = np.zeros((len(self.mst) + 1, 2))  # threshold, cluster_count
        
        # Start with all nodes as separate clusters
        cluster_counts[0, 0] = float('inf')  # threshold = infinity
        cluster_counts[0, 1] = self.n_nodes  # all nodes separate
        
        # Process MST edges from highest to lowest threshold
        for i in range(len(self.mst)):
            mst_edge_idx = self.mst[i]
            threshold = self.edges[mst_edge_idx, 2]
            # Each MST edge reduces cluster count by 1
            cluster_counts[i + 1, 0] = threshold
            cluster_counts[i + 1, 1] = self.n_nodes - (i + 1)
        
        return cluster_counts

    def _cluster_count_by_edge_count(self):
        """
        Returns array of (edge_count, cluster_count) pairs.
        For each cumulative edge count, calculates the number of clusters.
        This includes non-MST edges between MST edges.
        """
        edges_by_thresh = self.edges_by_threshold
        
        if len(edges_by_thresh) == 0:
            return np.array([[0, self.n_nodes]])
        
        # Calculate cluster counts
        result = []
        result.append([0, self.n_nodes])  # Start with 0 edges, all nodes separate
        
        for i in range(len(edges_by_thresh)):
            edge_count = int(edges_by_thresh[i, 0])
            # Each MST edge reduces cluster count by 1
            cluster_count = self.n_nodes - (i + 1)
            result.append([edge_count, cluster_count])
        
        return np.array(result)

    def export_for_interactive_viz(self):
        """
        Export minimal data structure for client-side cluster computation.
        Returns dict with MST edges for JavaScript to rebuild clusters on-demand.
        
        Returns:
            dict with:
            - n_nodes: int
            - mst_edges: [[node1, node2, value], ...] sorted by value descending
        """
        mst_edges = []
        for i in range(len(self.mst)):
            edge_idx = self.mst[i]
            mst_edges.append([
                int(self.edges[edge_idx, 0]),
                int(self.edges[edge_idx, 1]),
                float(self.edges[edge_idx, 2])
            ])
        
        return {
            'n_nodes': int(self.n_nodes),
            'mst_edges': mst_edges
        }
