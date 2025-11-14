"""

The DataMatrix class is an abstraction around dense and sparse matrices

"""
import warnings
warnings.filterwarnings("ignore", module='numpy')
from os import PathLike
import pandas as pd
import scipy.sparse
import numpy as np
import h5py
from typing import Dict, List, Tuple, Union, Optional
from domainator.utils import get_file_type

#TODO: use hdf5plugin to blosc compress datasets: http://www.silx.org/doc/hdf5plugin/latest/usage.html

UTF8_h5py_encoding = h5py.string_dtype(encoding='utf-8')

#TODO: maybe sub-classes for dense vs sparse?
#TODO: maybe add an attr for value label

class DataMatrix():
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
    _DOMAINATOR_MATRIX_FILE_VERSION_ATTR="DOMAINATOR_MATRIX_FILE_VERSION"
    _DATA_TYPE_ATTR="DATA_TYPE"
    """(optional) single string value describing the type of data in the matrix (e.g. 'score', 'norm_score', 'row_norm_score', 'score_dist', 'bool', 'efi_score')"""


    _DENSE = "DENSE"
    _SPARSE_CSR = "SPARSE_CSR"
    _POSSIBLE_ARRAY_TYPES={_DENSE, _SPARSE_CSR}


    
    _DOMAINATOR_MATRIX_FILE_VERSION = "1.0"  #TODO: increment after any breaking changes

    def __init__(self, data=None, row_names:List[str]=None, col_names:List[str]=None, row_lengths:Optional[List[int]]=None, col_lengths:Optional[List[int]]=None, data_type:str=None):
        """
            
        
        """

        # define instance attributes
        self.sparse = False
        self.rows = None # row names
        self.columns = None # column names
        self.row_lengths = None # (optional) int list representing the lengths of the sequences in the rows. Used for some statistics, like efi_score
        self.column_lengths = None # (optional) int array (not used if SYMMETRIC LABELS is True)
        self.row_to_idx = None # row name to rows_idx (inverse of self.rows)
        self.column_to_idx = None # col name to rows_idx (inverse of self.cols)
        self.data_type = "" # single string value describing the type of data in the matrix (e.g. 'score', 'norm_score', 'row_norm_score', 'score_dist', 'bool', 'efi_score')
        self.data = None # numpy array or scipy sparse_csr array

        if data is None:
            self.data = np.zeros((0,0))
        else:
            rows_list, cols_list = self.validate_matrix_data(data, row_names, col_names, row_lengths, col_lengths)
            self.rows = rows_list
            self.columns = cols_list
            # self.row_lengths = list(row_lengths) if row_lengths is not None else None
            # self.column_lengths = list(col_lengths) if col_lengths is not None else None
            self.row_to_idx = {v:i for i,v in enumerate(row_names)}
            self.data_type = str(data_type) if data_type is not None else ""

            if row_names == col_names:
                self.column_to_idx = self.row_to_idx
            else:
                self.column_to_idx = {v:i for i,v in enumerate(col_names)} 

            if row_lengths is not None:
                self.row_lengths = row_lengths.copy()
                if col_lengths == row_lengths:
                    self.column_lengths = self.row_lengths
                else:
                    self.column_lengths = col_lengths.copy()

            if scipy.sparse.issparse(data):
                self.sparse = True
                if not isinstance(type(data), scipy.sparse.csr_array):
                    self.data = scipy.sparse.csr_array(data)
                else:
                    self.data = data.copy() 
                self.data.eliminate_zeros()
                self.data.sort_indices()
            else: #dense
                self.data = data.copy()

    @property
    def shape(self):
        return self.data.shape
    
    @property
    def size(self):
        """
            Returns:
                int: number of values in the matrix
        """
        return np.prod(self.shape)

    @property
    def symmetric(self):
        return self.rows == self.columns

    @property
    def symmetric_values(self) -> bool:
        """
            Returns True if the matrix has symmetric values (i.e., the matrix is symmetric in its data).
            Uses np.isclose to allow for floating point imprecision.
            
            For sparse matrices, uses an optimized approach that first checks if the 
            non-zero pattern is symmetric before comparing values.
        """
        if not self.symmetric:
            return False # If it doesn't have symmetric labels, then it can't have symmetric values.
        
        # Check if matrix values are symmetric: M[i,j] == M[j,i] for all i,j
        
        if self.sparse:
            # For sparse matrices, optimize by checking non-zero structure first
            # Get upper and lower triangular non-zero elements (excluding diagonal)
            upper_triangular = self.triangular(side='upper', include_diagonal=False, skip_zeros=True, index_style='index')
            lower_triangular = self.triangular(side='lower', include_diagonal=False, skip_zeros=True, index_style='index')
            
            # First, check if the number of non-zeros is the same
            if len(upper_triangular) != len(lower_triangular):
                return False
            
            # For symmetric matrix: upper[i,j] should equal lower[j,i]
            # Since both lists are sorted by (row, col), we can create a sorted list
            # of transposed lower triangular positions for comparison
            
            # Create list of (col, row, value) from lower triangular - this gives us the transpose
            # Then sort it the same way as upper triangular (by row, then col)
            lower_transposed = [(l_c, l_r, l_v) for l_r, l_c, l_v in lower_triangular]
            lower_transposed.sort()  # Sort by (row, col) which is now (l_c, l_r)
            
            # Now compare element by element
            for (u_r, u_c, u_v), (lt_r, lt_c, lt_v) in zip(upper_triangular, lower_transposed):
                # Check positions match
                if (u_r, u_c) != (lt_r, lt_c):
                    return False
                # Check values match (with floating point tolerance)
                if not np.isclose(u_v, lt_v):
                    return False
            
            return True
        else:
            # For dense matrices, check all upper triangular elements
            # We only need to check the upper triangular part against the lower triangular part
            upper_triangular = self.triangular(side='upper', include_diagonal=False, skip_zeros=False, index_style='index')
            
            # For each upper triangular element (i,j), check if it equals the corresponding lower element (j,i)
            for row_idx, col_idx, value in upper_triangular:
                # Get the corresponding element in the lower triangular (transpose position)
                transpose_value = self.data[col_idx, row_idx]
                
                # Use np.isclose to handle floating point imprecision
                if not np.isclose(value, transpose_value):
                    return False
            
            return True

    def __len__(self):
        return self.data.shape[0]

    @classmethod
    def from_file(cls, matrix_file):
        out = cls()
        file_type = get_file_type(matrix_file)
        if file_type == "hdf5": 
            out._read_hdf5(matrix_file)
        else: #assume dense tsv text file as parse as such
            out._read_text_dense_matrix(matrix_file) 
        return out

    def _read_hdf5(self, matrix_file: Union[PathLike, str]):
        
        f = h5py.File(matrix_file)
        
        ### validate file format  ###
        for attribute in (self._DOMAINATOR_MATRIX_FILE_VERSION_ATTR, self._ARRAY_TYPE_ATTR, self._SYMMETRIC_LABELS_ATTR):
            if attribute not in f.attrs:
                raise ValueError(f"attribute {attribute} not found in hdf5 file. Are you sure it is a Domainator formatted file?")
        file_version = f.attrs[self._DOMAINATOR_MATRIX_FILE_VERSION_ATTR]
        if file_version != self._DOMAINATOR_MATRIX_FILE_VERSION:
            warnings.warn(f"Matrix file version ({file_version}) in file {matrix_file} does not match program-supported version ({self._DOMAINATOR_MATRIX_FILE_VERSION})")
        array_type = f.attrs[self._ARRAY_TYPE_ATTR]
        if array_type not in self._POSSIBLE_ARRAY_TYPES:
            raise ValueError(f"Array type {array_type} not recognized")

        if array_type == "DENSE":
            for dataset_name in (self._ROW_LABELS_DATASET, self._DENSE_DATASET):
                if dataset_name not in f:
                    raise ValueError(f"expected dataset {self._DENSE_DATASET} not found in file {matrix_file}")
        elif array_type == "SPARSE_CSR": # sparse
            for dataset_name in (self._ROW_LABELS_DATASET, self._SPARSE_VALUES_DATASET, self._SPARSE_CSR_INDICES_DATASET, self._SPARSE_CSR_INDPTR_DATASET):
                if dataset_name not in f:
                    raise ValueError(f"expected dataset {dataset_name} not found in file {matrix_file}")

        if not f.attrs[self._SYMMETRIC_LABELS_ATTR]:
            if self._COL_LABELS_DATASET not in f:
                raise ValueError(f"expected dataset {self._COL_LABELS_DATASET} not found in file {matrix_file}")
        

        ### store data ###
        self.sparse = False
        self.rows = None
        self.columns = None
        self.row_to_idx = None
        self.column_to_idx = None
        self.data = None
        self.row_lengths = None
        self.column_lengths = None
        self.data_type = ""

        if self._DATA_TYPE_ATTR in f.attrs:
            self.data_type = f.attrs[self._DATA_TYPE_ATTR]

        self.rows, self.row_to_idx = self.parse_axis_index(f[self._ROW_LABELS_DATASET])
        if self._ROW_LENGTHS_DATASET in f:
            self.row_lengths = f[self._ROW_LENGTHS_DATASET][:]
            if f.attrs[self._SYMMETRIC_LABELS_ATTR]:
                self.column_lengths = self.row_lengths
            else:
                self.column_lengths = f[self._COL_LENGTHS_DATASET][:]
        

        if f.attrs[self._SYMMETRIC_LABELS_ATTR]:
            self.columns = self.rows
            self.column_to_idx = self.row_to_idx
        else:
            self.columns, self.column_to_idx = self.parse_axis_index(f[self._COL_LABELS_DATASET])

        if f.attrs[self._ARRAY_TYPE_ATTR] == "DENSE":
            self.sparse = False
            self.data = f[self._DENSE_DATASET][:]
        elif f.attrs[self._ARRAY_TYPE_ATTR] == "SPARSE_CSR":
            self.sparse = True
            self.data = scipy.sparse.csr_array((f[self._SPARSE_VALUES_DATASET][:], f[self._SPARSE_CSR_INDICES_DATASET][:], f[self._SPARSE_CSR_INDPTR_DATASET][:]), (len(self.rows), len(self.columns)) )

    
    def iter_data(self, index_style="name", skip_zeros=True):
        """
            yields (row_name, column_name, value).
            for dense matrices iterates over every value in the matrix
            for sparse matrices, skips zeros.

            Iterates in a row-first manner. That is, all values in row 0, then all values in row 1, etc.
        """
        
        index_style = index_style.lower()
        if index_style not in {"name", "index"}:
            raise ValueError("index_style must be 'name' or 'index'.")


        if not self.sparse: # it's a dense matrix
            for r in range(self.data.shape[0]):
                for c in range(self.data.shape[1]):
                    value = self.data[r,c]
                    if not skip_zeros or value != 0:
                        if index_style == "name":
                            yield self.rows[r], self.columns[c], value
                        else:
                            yield r, c, value
        else: # sparse
            if skip_zeros: # we can iterate over just the non-zeros
                self.data.eliminate_zeros()
                self.data.sort_indices()
                #https://stackoverflow.com/questions/4319014/iterating-through-a-scipy-sparse-vector-or-matrix
                for r in range(self.data.shape[0]):
                    for ind in range(self.data.indptr[r], self.data.indptr[r+1]):
                        c = self.data.indices[ind]
                        if index_style == "name":
                            yield self.rows[r], self.columns[c], self.data.data[ind]
                        else:
                            yield r, c, self.data.data[ind]
            else: # we have to iterate over the entire matrix
                for r in range(self.data.shape[0]):
                    for c in range(self.data.shape[1]):
                        value = self.data[r,c]
                        if index_style == "name":
                            yield self.rows[r], self.columns[c], value
                        else:
                            yield r, c, value
                    

    def triangular(self, side="lower", include_diagonal=False, skip_zeros=False, index_style="name", agg=None):
        """
            Iterate over the triangular part of the matrix.

            side: can be "lower" or "upper".
            index_style: can be "name" or "index".
            skip_zeros: if True then don't return rows with zero as the value.
            agg: a function taking two values as arguments, if supplied, then replaces each returned value with the 
                result of agg(value, symmetric_complement).
                if a value is zero, but the complement is not zero, then the result of agg is returned even when skip_zeros is true (unless it is also zero)

            results sorted by row_index, then column_index
            returns [(row_{name|col}, column_{name|col}, value),]
        
        """

        if self.data.shape[0] != self.data.shape[1]:
            raise NotImplementedError(f"Triangular only supported for square matrices. Not ({self.data.shape[0]}, {self.data.shape[1]}")

        side = side.lower()
        if side not in {"upper", "lower"}:
            raise ValueError("side must be 'lower' or 'upper'")

        index_style = index_style.lower()
        if index_style not in {"name", "index"}:
            raise ValueError("index_style must be 'name' or 'index'.")


        out = list()

        if not self.sparse or not skip_zeros or agg is not None:  # dense matrix or we want every value, inlcuding zeros, or we are doing an agg.
            for r in range(self.data.shape[0]):
                if side == "lower":
                    # Lower triangular: column index <= row index
                    start_c = 0
                    end_c = r + 1 if include_diagonal else r
                else:  # upper
                    # Upper triangular: column index >= row index
                    start_c = r if include_diagonal else r + 1
                    end_c = self.data.shape[1]
                
                for c in range(start_c, end_c):
                    value = self.data[r, c]
                    
                    # Apply aggregation if agg function is provided
                    if agg is not None:
                        # Get the symmetric complement (transpose position)
                        symmetric_value = self.data[c, r]
                        # When skip_zeros is True with agg, only skip if BOTH values are zero
                        if skip_zeros and value == 0 and symmetric_value == 0:
                            continue
                        value = agg(value, symmetric_value)
                    elif skip_zeros and value == 0:
                        # Without agg, skip if value is zero
                        continue
                    
                    if index_style == "name":
                        out.append((self.rows[r], self.columns[c], value))
                    else:
                        out.append((r, c, value))
        else:  # sparse matrix
            self.data.eliminate_zeros()
            self.data.sort_indices()

            for r in range(self.data.shape[0]):
                if side == "lower":
                    # Lower triangular: column index <= row index
                    start_c = 0
                    end_c = r + 1 if include_diagonal else r
                else:  # upper
                    # Upper triangular: column index >= row index
                    start_c = r if include_diagonal else r + 1
                    end_c = self.data.shape[1]
                
                if skip_zeros:
                    # Only iterate over non-zero values in the triangular region
                    for ind in range(self.data.indptr[r], self.data.indptr[r+1]):
                        c = self.data.indices[ind]
                        
                        # Check if this column is in the triangular region
                        if start_c <= c < end_c:
                            value = self.data.data[ind]
    
                            if index_style == "name":
                                out.append((self.rows[r], self.columns[c], value))
                            else:
                                out.append((r, c, value))
                else:
                    assert False # This case should have been handled above, treating the sparse matrix like a dense matrix.
        return out
        

    def itervalues(self):
        """ yield values from the matrix one by one by row then column.
        
        """

        for r in range(self.shape[0]):
            for c in range(self.shape[1]):
                yield self.data[r,c]
    
    def save(self, file):
        """write the matrix to an hdf5 file

        Args:
            file (path or string): 
        """

        if self.sparse:
            self.write_sparse(self.data, file, self.row_to_idx, self.column_to_idx)
        else: # dense
            self.write_dense(self.data, file, self.row_to_idx, self.column_to_idx)


    def convert_to_sparse(self):
        """Inplace conversion to csr sparse
           

        """
        if self.sparse:
            pass
        else:
            self.sparse = True
            self.data = scipy.sparse.csr_array(self.data)
            self.data.eliminate_zeros()

    def convert_to_dense(self):
        """Inplace conversion to dense

        """
        if self.sparse:
            self.sparse = False
            self.data = self.data.toarray()
        else:
            pass

    def toarray(self):
        if self.sparse:
            return self.data.toarray()
        else:
            return self.data

    def _read_text_dense_matrix(self, matrix_file: Union[PathLike, str]):
        self.sparse = False
        try:
            mat_file = pd.read_csv(matrix_file, sep="\t", index_col=0)
        except UnicodeDecodeError:
            raise ValueError(f"Expecting utf-8 encoded matrix file, but found invalid unicode byte. You may be trying to pass in an hdf5 file, but using an unrecognized extension.")
        self.rows = list(mat_file.index)
        self.columns = list(mat_file.columns)

        self.row_to_idx = {v:i for i,v in enumerate(self.rows)}
        self.column_to_idx = {v:i for i,v in enumerate(self.columns)}

        if all(a == b for a,b in zip(self.rows, self.columns)):
            self.columns = self.rows

        self.data = mat_file.to_numpy()

    def write(self, file_name: Union[PathLike, str], output_type: str = "dense"):
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
    def write_sparse(cls, matrix, file_name, row_names, col_names, row_lengths=None, col_lengths=None, data_type=""):
        """Writes a dok sparse matrix to an hdf5 file
        
        """

        # check validity of input data
        row_names, col_names = cls.validate_matrix_data(matrix, row_names, col_names, row_lengths, col_lengths)
        
        if not isinstance(matrix, scipy.sparse.csr_array):
            matrix = scipy.sparse.csr_array(matrix)

        matrix.eliminate_zeros() # remove zeros, so that the sparse file is as small as possible
        matrix.sort_indices() # sort indices so iteration happens by row and column in order

        f = h5py.File(file_name, 'w')
        f.attrs[cls._DOMAINATOR_MATRIX_FILE_VERSION_ATTR] = cls._DOMAINATOR_MATRIX_FILE_VERSION
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
        
        f.close()


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
    def validate_matrix_data(cls, matrix: Union[np.array, scipy.sparse.sparray], row_names:List[str], col_names:List[str], row_lens:Optional[List[int]]=None, col_lens:Optional[List[int]]=None) -> Tuple[List[str], List[str]]:
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
    def write_dense(cls, matrix:Union[np.array, scipy.sparse.sparray], file_name: Union[PathLike, str], row_names:List[str], col_names:List[str], row_lengths:Optional[List[int]]=None, col_lengths:Optional[List[int]]=None, data_type:Optional[str]=""):
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
        

        f = h5py.File(file_name, 'w')
        f.attrs[cls._DOMAINATOR_MATRIX_FILE_VERSION_ATTR] = cls._DOMAINATOR_MATRIX_FILE_VERSION
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
            
        f.close()

    @classmethod
    def write_dense_text(cls, matrix:Union[np.array, scipy.sparse.sparray], file_name: Union[PathLike, str], row_names:List[str], col_names:List[str], row_lengths:Optional[List[int]]=None, col_lengths:Optional[List[int]]=None, data_type:Optional[str]=""):
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