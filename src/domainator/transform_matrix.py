"""Transforms a score matrix into various kinds of normalized martices
   
   For matrices with any data type, can convert from dense to sparse, or sparse to dense, or change the datatype.
"""

import warnings
warnings.filterwarnings("ignore", module='numpy')
from jsonargparse import ArgumentParser, ActionConfigFile
import sys

from domainator.data_matrix import DataMatrix, mst_knn_edge_index_dict
from domainator import __version__, RawAndDefaultsFormatter
import scipy.sparse
import numpy as np
from domainator.utils import get_file_type

PASS_THROUGH_MODE = "none"

def is_sparse(matrix):
    """
    Check if a matrix is sparse.
    """
    return type(matrix) is scipy.sparse.dok_array or type(matrix) is scipy.sparse.coo_array or type(matrix) is scipy.sparse.csr_array

def norm_score(array, row_lengths=None, col_lengths=None):
    """
        score / min(row_max, col_max)
    """
    if is_sparse(array): #sparse matrix
        array = scipy.sparse.coo_matrix(array) 
        rowmaxes = dict()
        max_coo = array.max(axis=1)
        for r,c,v in zip(max_coo.row, max_coo.col, max_coo.data):
            rowmaxes[r] = v

        colmaxes = dict()        
        colmax_coo = array.max(axis=0)
        for r,c,v in zip(colmax_coo.row, colmax_coo.col, colmax_coo.data):
            colmaxes[c] = v

        out = scipy.sparse.dok_array(array.shape,dtype=np.float64)
        for r,c,v in zip(array.row, array.col, array.data):
            out[r, c] = v / min(rowmaxes[r], colmaxes[c])
        return out
    else: # dense matrix
        return array/np.minimum( array.max(axis=1)[:,np.newaxis], array.max(axis=0)[np.newaxis,:] )

def row_norm_score(array, row_lengths=None, col_lengths=None):
    """
        score / row_max
    """
    if is_sparse(array): #sparse matrix
        array = scipy.sparse.coo_matrix(array) 
        rowmaxes = dict()
        max_coo = array.max(axis=1)
        for r,c,v in zip(max_coo.row, max_coo.col, max_coo.data):
            rowmaxes[r] = v
        out = scipy.sparse.dok_array(array.shape,dtype=np.float64)
        for r,c,v in zip(array.row, array.col, array.data):
            out[r, c] = v / rowmaxes[r]
        return out
    else: # dense matrix
        return array/(array.max(axis=1)[:,np.newaxis])

def efi_score(array, row_lengths=None, col_lengths=None):
    """
        -log10[2^(-score) * (row_seq_length * col_seq_length)]

        scores less than 0 are set to 0.

        where row_length and col_length are the lenghts of the sequences or profiles.

        This is the score used by EFI-EST.


        
    """
    if row_lengths is None or col_lengths is None:
        raise ValueError("Row and column lengths must be provided for EFI score.")

    if is_sparse(array): #sparse matrix
        array = scipy.sparse.dok_array(array) #TODO: maybe there is a more efficient way to do this without duplicating the data
        out = scipy.sparse.dok_array(array.shape,dtype=np.float64)

        for (r,c),v in array.items():
            #score = -np.log10(2**(-v) * (row_lengths[r] * col_lengths[c]))
            score = v * np.log10(2) - np.log10(row_lengths[r] * col_lengths[c]) # refactored to avoid overflow
            if score > 0:
                out[r, c] = score
        return out
    
    else: # dense matrix
        out = -np.log10(2.0**(-array) * (row_lengths[:,np.newaxis] * col_lengths[np.newaxis,:]))
        return out * (out > 0)
    
def efi_score_dist(array, row_lengths=None, col_lengths=None):
    """
        1 - (efi_score / min(row_max, col_max))
    """
    out = efi_score(array, row_lengths=row_lengths, col_lengths=col_lengths)
    if is_sparse(out): #sparse matrix
            out = out.toarray()
    return 1 - (out/(np.minimum( out.max(axis=1)[:,np.newaxis], out.max(axis=0)[np.newaxis,:] ) ) )
        
        

def score_dist(array, row_lengths=None, col_lengths=None): 
    """
    1 - (score / min(row_max, col_max))
    """
    if is_sparse(array): #sparse matrix
        array = array.toarray()
        # raise ValueError("Sparse distance matrices not implemented.")
        # array = scipy.sparse.coo_matrix(array) 
        # rowmaxes = dict()
        # rowmax_coo = array.max(axis=1)
        # for r,c,v in zip(rowmax_coo.row, rowmax_coo.col, rowmax_coo.data):
        #     rowmaxes[r] = v

        # print(rowmaxes)
        # colmaxes = dict()        
        # colmax_coo = array.max(axis=0)
        # for r,c,v in zip(colmax_coo.row, colmax_coo.col, colmax_coo.data):
        #     colmaxes[c] = v
        # print(colmaxes)

        # out = scipy.sparse.dok_array(array.shape,dtype=np.float64)
        # for r,c,v in zip(array.row, array.col, array.data):
        #     print(r,c,v)
        #     new_val = 1 - (v / min(rowmaxes[r], colmaxes[c]))
        #     if new_val != 0:
        #         out[r, c] = new_val
        # return out
    
    #dense matrix
    #TODO: in https://doi.org/10.1093/nar/gkaa635, they use -ln( (score / min(row_max, col_max)) ), instead of 1 - ...
    #       unsure whether we should use the linear or log version.
    return 1 - (array/(np.minimum( array.max(axis=1)[:,np.newaxis], array.max(axis=0)[np.newaxis,:] ) ) )


MODES = { 
        "score": None, 
        "bool": None, 
        "norm_score": norm_score,
        "row_norm_score": row_norm_score,
        "score_dist": score_dist,
        "efi_score": efi_score,
        "efi_score_dist": efi_score_dist}


def _mst_knn_arg(value):
    value = int(value)
    if value <= 1:
        raise ValueError("--mst_knn must be an integer greater than 1")
    return value


def apply_lower_bound(array, lower_bound):
    if lower_bound is None:
        return array

    if scipy.sparse.issparse(array):
        out = scipy.sparse.csr_array(array)
        out.data[out.data <= lower_bound] = 0
        out.eliminate_zeros()
        return out

    out = np.array(array, copy=True)
    out[out <= lower_bound] = 0
    return out


def apply_mst_knn_sparsification(matrix, k, lower_bound=0):
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("--mst_knn requires a square matrix.")

    edge_dict = mst_knn_edge_index_dict(matrix, k, lower_bound=lower_bound)
    diagonal_indices = np.arange(matrix.shape[0])

    if scipy.sparse.issparse(matrix.data):
        out = scipy.sparse.dok_array(matrix.shape, dtype=matrix.data.dtype)
        for index in diagonal_indices:
            value = matrix.data[index, index]
            if value != 0:
                out[index, index] = value
        for source_idx, target_idx in edge_dict:
            forward_value = matrix.data[source_idx, target_idx]
            reverse_value = matrix.data[target_idx, source_idx]
            if forward_value != 0:
                out[source_idx, target_idx] = forward_value
            if reverse_value != 0:
                out[target_idx, source_idx] = reverse_value
        return scipy.sparse.csr_array(out)

    out = np.zeros_like(matrix.data)
    out[diagonal_indices, diagonal_indices] = matrix.data[diagonal_indices, diagonal_indices]
    for source_idx, target_idx in edge_dict:
        out[source_idx, target_idx] = matrix.data[source_idx, target_idx]
        out[target_idx, source_idx] = matrix.data[target_idx, source_idx]
    return out

def transform_matrix(array, mode, row_lengths=None, col_lengths=None):
    """
    Transform an array to a different mode.
    """
    if mode == "rank":
        return scipy.stats.rankdata(array, method='min', axis=1)
    elif mode == "score":
        return array
    elif mode == "bool":
        return (array > 0).astype(np.int8)
    else:
        return MODES[mode](array, row_lengths=row_lengths, col_lengths=col_lengths)
    
def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument("-i", "--input", default=None, required=True, type=str, help="name of input matrix file.")

    parser.add_argument("--input_type", choices=set(MODES.keys()).union({None}), default=None, help="If the type of the input matrix is not specified, specify the type of input matrix. Default: None")
    
    parser.add_argument("--dense", type=str, default=None, help="Write a dense distance matrix hdf5 file to this path.")
    
    parser.add_argument("--dense_text", type=str, default=None, help="Write a dense distance matrix tsv file to this path.")
    
    parser.add_argument("--sparse", type=str, default=None, help="Write a sparse distance matrix hdf5 file to this path.")

    parser.add_argument("--mode", type=str, required=False, default=PASS_THROUGH_MODE, choices=set(MODES.keys()).union({PASS_THROUGH_MODE}), 
                        help="what kind of values should be in the output matrix. none: preserve the input values and data type, score: raw score, bool: 1 if a hit otherwise 0, score_dist: 1 - (score / min(row_max, col_max)), norm_score: score/min(row_max, col_max), efi_score: -log10[2^(-score) * (input_seq_length * reference_seq_length)], efi_score_dist: 1 - (efi_score / min(row_max, col_max)). Default: none")

    parser.add_argument('--lb', type=float, default=None, required=False,
                        help="Zero out all values less than or equal to this threshold after any mode transformation.")
    parser.add_argument('--mst_knn', type=_mst_knn_arg, required=False, default=None,
                        help="Keep only the maximum spanning tree plus OR-symmetric k-nearest-neighbor edges, using the post-transform values.")
    
    parser.add_argument('--config', action=ActionConfigFile)

    # TODO: option to supply row and/or column lengths if not supplied in the input matrix.

    # TODO: option to supply row and/or column names to override those in the input matrix.
    params = parser.parse_args(argv)



    if params.dense is not None and get_file_type(params.dense) != "hdf5":
        raise ValueError("Please use an hdf5 related extension for the --dense output, such as .h5, .hdf5, or .hdf.")
    if params.sparse is not None and get_file_type(params.sparse) != "hdf5":
        raise ValueError("Please use an hdf5 related extension for the --sparse output, such as .h5, .hdf5, or .hdf.")
    dense = params.dense
    dense_text = params.dense_text
    sparse = params.sparse

    if ((dense is None) and (dense_text is None) and (sparse is None)):
        raise ValueError("No output specified! Please specify at least one of: dense, dense_text, sparse")

    if sparse is not None and params.mode == "score_dist":
        raise ValueError("Sparse distance matrices not implemented.")

    matrix = DataMatrix.from_file(params.input)

    if params.input_type is not None:
        if matrix.data_type == "" or matrix.data_type == params.input_type:
            pass
        else:
            warnings.warn(f"Input matrix format is '{matrix.data_type}', but --input_type is '{params.input_type}'. The input matrix format is overridden.")
        matrix.data_type = params.input_type

    requested_mode = None if params.mode == PASS_THROUGH_MODE else params.mode

    if requested_mode is not None:
        if matrix.data_type == "":
            warnings.warn("Input matrix format not specified, assuming 'score'")
            matrix.data_type = "score"

        if matrix.data_type != "score" and matrix.data_type != requested_mode:
            raise ValueError(f"Input matrix format is '{matrix.data_type}', but --mode is '{requested_mode}'. Transforming between data types is only supported from a 'score' matrix.")

        if matrix.data_type == "score" and matrix.data_type != requested_mode:
            matrix.data = transform_matrix(matrix.data, requested_mode, matrix.row_lengths, matrix.column_lengths)

        matrix.data_type = requested_mode

    if params.mst_knn is not None:
        matrix.data = apply_mst_knn_sparsification(matrix, params.mst_knn, lower_bound=0 if params.lb is None else params.lb)

    if params.lb is not None:
        matrix.data = apply_lower_bound(matrix.data, params.lb)

    ### Run
    if dense is not None:
        matrix.write(dense, "dense")
    if dense_text is not None:
        matrix.write(dense_text,"dense_text")
    if sparse:
        matrix.write(sparse, "sparse")

def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])
  