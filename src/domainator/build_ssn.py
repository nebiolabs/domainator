"""Generates Sequence Similarity Networks
    
Build a sequence similarity network and do analysis related to that.
"""

import warnings
warnings.filterwarnings("ignore", module='numpy')
from jsonargparse import ArgumentParser, ActionConfigFile
from concurrent.futures import ThreadPoolExecutor
from os import PathLike
import os
import sys
from domainator.cytoscape import write_cytoscape_xgmml
from domainator.data_matrix import DataMatrix, MaxTree, mst_knn_edge_index_dict
from domainator.output_guardrails import add_max_output_gb_argument, max_output_gb_to_bytes, OutputSizeLimitExceeded, make_temporary_output_path
from scipy.sparse.csgraph import connected_components
import pandas as pd
from pathlib import Path
from domainator import __version__, RawAndDefaultsFormatter
from typing import List, Union, Dict
import numpy as np
from domainator.utils import get_palette, list_and_file_to_dict_keys
from domainator.color_genbank import read_color_table

SCORE_COLUMN="SSN_SCORE"
CLUSTER_COLUMN="SSN_cluster"


def subset_matrix_by_labels(matrix: DataMatrix, subset_labels) -> DataMatrix:
    """Filter rows and columns to the labels present in subset_labels."""
    if subset_labels is None:
        return matrix

    row_labels = [label for label in matrix.rows if label in subset_labels]
    if matrix.symmetric_labels:
        col_labels = row_labels
    else:
        col_labels = [label for label in matrix.columns if label in subset_labels]

    row_indices = [matrix.row_to_idx[label] for label in row_labels]
    col_indices = row_indices if matrix.symmetric_labels else [matrix.column_to_idx[label] for label in col_labels]

    row_lengths = None if matrix.row_lengths is None else matrix.row_lengths[row_indices]
    if matrix.symmetric_labels:
        col_lengths = row_lengths
    else:
        col_lengths = None if matrix.column_lengths is None else matrix.column_lengths[col_indices]

    data = matrix.data[row_indices][:, col_indices]
    return matrix.__class__(data, row_labels, col_labels, row_lengths, col_lengths, matrix.data_type)

def rename_labels_by_frequency(labels:np.ndarray) -> np.ndarray:
    """rename labels by frequency, such that the most frequent label is 0, the second most frequent is 1, etc.

    Args:
        labels (np.ndarray): labels

    Returns:
        np.ndarray: labels
    """
    unique_labels, counts = np.unique(labels, return_counts=True)
    label_counts = dict(zip(unique_labels, counts))
    sorted_labels = sorted(label_counts, key=label_counts.get, reverse=True)
    new_labels = np.zeros(len(labels), dtype=int)
    for i, label in enumerate(sorted_labels):
        new_labels[labels == label] = i
    return new_labels


def _mst_knn_arg(value: Union[str, int]) -> int:
    value = int(value)
    if value <= 1:
        raise ValueError("--mst_knn must be an integer greater than 1")
    return value


def cluster_labels_from_tree(tree: MaxTree, lb: float) -> np.ndarray:
    """Return size-ranked connected component labels from MST edges above the threshold."""

    parent = np.arange(tree.n_nodes, dtype=int)

    def find(node_idx: int) -> int:
        root = node_idx
        while parent[root] != root:
            root = parent[root]

        while parent[node_idx] != node_idx:
            next_idx = parent[node_idx]
            parent[node_idx] = root
            node_idx = next_idx

        return root

    def union(left_idx: int, right_idx: int) -> None:
        left_root = find(left_idx)
        right_root = find(right_idx)
        if left_root != right_root:
            parent[left_root] = right_root

    for source_idx, target_idx, score in tree.iter_mst_edges():
        if score <= lb:
            break
        union(source_idx, target_idx)

    labels = np.empty(tree.n_nodes, dtype=int)
    root_to_label = dict()
    next_label = 0
    for node_idx in range(tree.n_nodes):
        root = find(node_idx)
        if root not in root_to_label:
            root_to_label[root] = next_label
            next_label += 1
        labels[node_idx] = root_to_label[root]

    return rename_labels_by_frequency(labels) + 1


def cluster_labels_from_graph(matrix: DataMatrix, lb: float) -> np.ndarray:
    """Return size-ranked connected component labels from the thresholded input graph."""
    labels = connected_components(matrix.data > lb, directed=False, return_labels=True)[1]
    return rename_labels_by_frequency(labels) + 1


def iter_default_ssn_edges(matrix: DataMatrix, lb: float):
    """Yield default-path SSN edges as index/index/score tuples in legacy output order."""
    for source_idx, target_idx, score in matrix.triangular_iter(skip_zeros=True, agg=max, index_style="index"):
        if score <= lb or source_idx == target_idx:
            continue
        if matrix.rows[target_idx] < matrix.rows[source_idx]:
            source_idx, target_idx = target_idx, source_idx
        yield source_idx, target_idx, score


def iter_mst_ssn_edges(tree: "MaxTree", rows, lb: float):
    """Yield MST edges above lb as (source_name, target_name, score) in canonical order."""
    for source_idx, target_idx, score in tree.iter_mst_edges():
        if score <= lb:
            break
        if source_idx == target_idx:
            continue
        source, target = rows[source_idx], rows[target_idx]
        if target < source:
            source, target = target, source
        yield source, target, score


def iter_mst_knn_ssn_edges(edge_dict, rows):
    """Yield MST+kNN edges as (source_name, target_name, score) in canonical order."""
    for (source_idx, target_idx), score in edge_dict.items():
        source, target = rows[source_idx], rows[target_idx]
        if target < source:
            source, target = target, source
        yield source, target, score


def build_ssn(matrix: DataMatrix, lb:float=0, metadata_files:List[Union[str, PathLike]]=None, color_by:str=None, color_table:Dict[str,str]=None, 
              xgmml:Union[str, PathLike]=None, cluster:bool=False, cluster_tsv:Union[str, PathLike]=None, no_cluster_header:bool=False, color_table_out:str = None, mst:bool = False, mst_knn: int = None, subset_labels=None, max_output_bytes=None):
    """build a sequence similarity network from a matrix

    Args:
        matrix (DataMatrix): symmetric matrix of edge values
        lb (float, optional): ignore scores lower than or equal to this. Defaults to 0.
        metadata_files (List[Union[str, PathLike]], optional): tab separated text files containg node annotations. Defaults to None.
        color_by (str, optional): color the nodes based on their value in this metadata field. Set to SSN_cluster to color by cluster. Defaults to None.
        xgmml (Union[str, PathLike], optional): write a Cytoscape xgmml file to this path. Defaults to None.
        cluster (bool, optional): If True, then add new annotations with connected cluster IDs. Defaults to False.
        cluster_tsv (Union[str, PathLike], optional): write a metadata tab-separated file with the cluster ids to this path. Defaults to None.
        no_cluster_header (bool, optional): If True, then the cluster tsv file will not have a header. Defaults to False.
        color_table_out (str, optional): write a color table to this path. Defaults to None.
        mst (bool, optional): If True, then write 
        mst_knn (int, optional): If set, emit the union of the MST and an OR-symmetric kNN graph.
        subset_labels (set, optional): keep only rows and columns with labels in this set. Defaults to None.

    Raises:
        ValueError: if the input does not have symmetric row and column labels
    """
    matrix = subset_matrix_by_labels(matrix, subset_labels)

    if not matrix.symmetric_labels:
        raise ValueError("Input does not have symmetric axis labels. Can only build an SSN from a symmetric matrix.")

    # if not matrix.symmetric_values:
    #     raise ValueError("Input not symmetric data values. Can only build an SSN from a symmetric matrix.")

    if cluster_tsv is not None:
        cluster = True
    tree = MaxTree(matrix) if mst or mst_knn is not None else None

    node_data = pd.DataFrame(index=matrix.rows)
    if cluster:
        if tree is None:
            node_data[CLUSTER_COLUMN] = cluster_labels_from_graph(matrix, lb)
        else:
            node_data[CLUSTER_COLUMN] = cluster_labels_from_tree(tree, lb)

    # Determine the edge source for the XGMML write.  For mst/mst_knn we free
    # the large matrix.data as soon as it is no longer needed.
    if mst:
        rows = matrix.rows
        del matrix  # free O(n²) data array; rows is a plain list
        xgmml_edge_rows = iter_mst_ssn_edges(tree, rows, lb)
        xgmml_edge_kwargs = dict(
            edge_rows=xgmml_edge_rows,
            edge_column_names=[SCORE_COLUMN],
            edge_metadata_types={SCORE_COLUMN: "real"},
            edge_endpoint_style="name",
        )
    elif mst_knn is not None:
        edge_dict = mst_knn_edge_index_dict(matrix, mst_knn, lower_bound=lb, tree=tree)
        rows = matrix.rows
        del matrix  # free O(n²) data array
        xgmml_edge_rows = iter_mst_knn_ssn_edges(edge_dict, rows)
        xgmml_edge_kwargs = dict(
            edge_rows=xgmml_edge_rows,
            edge_column_names=[SCORE_COLUMN],
            edge_metadata_types={SCORE_COLUMN: "real"},
            edge_endpoint_style="name",
        )
    else:
        xgmml_edge_kwargs = dict(
            edge_rows=iter_default_ssn_edges(matrix, lb),
            edge_column_names=[SCORE_COLUMN],
            edge_metadata_types={SCORE_COLUMN: "real"},
            edge_endpoint_style="index",
        )

    if metadata_files is not None:
        with ThreadPoolExecutor() as exe:
            dfs = list(exe.map(lambda f: pd.read_csv(f, sep="\t", index_col=0), metadata_files))
        for df in dfs:
            node_data = node_data.merge(df, how="left", left_index=True, right_index=True)

    if color_by is not None:
        if color_table is None:
            color_table = get_palette(node_data[color_by].unique())

    if xgmml is not None:
        temp_xgmml_path = make_temporary_output_path(xgmml)
        try:
            write_cytoscape_xgmml(
                temp_xgmml_path,
                node_data,
                **xgmml_edge_kwargs,
                name=Path(xgmml).stem,
                color_by=color_by,
                color_table=color_table,
                max_output_bytes=max_output_bytes,
                output_description=f"XGMML network output '{xgmml}'",
                mitigation_options=["--lb", "--mst", "--mst_knn", "--subset"],
            )
            os.replace(temp_xgmml_path, xgmml)
            temp_xgmml_path = None
        finally:
            if temp_xgmml_path is not None and Path(temp_xgmml_path).exists():
                os.unlink(temp_xgmml_path)

    if cluster_tsv:
        index_label = None
        header = None
        if not no_cluster_header:
            header = ["cluster"]
            index_label = "contig"

        node_data[CLUSTER_COLUMN].to_csv(cluster_tsv, index_label=index_label, header=header, sep="\t")
    
    if color_table_out is not None:
        with open(color_table_out, "w") as out:
            color_table_items = list(color_table.items())
            color_table_items.sort(key=lambda x: (x[0] is None, x[0])) # sort lexically
            for domain, color in color_table_items:
                out.write(f"{domain}\t{color}\n")


def main(argv):
    parser = ArgumentParser(description=f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)
    
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-i', '--input', type=str, required=False,
                        help="pairwise similarity (or distance) scores. In the network, nodes will be derived from row and column names, and edges will be added between pairs of rows and columns meeting the threshold criteria. Format can be tab-separated text, or Domainator hdf5.")

    parser.add_argument('--lb', type=float, default=0, required=False,
                        help="exclude all edges with weights less than or equal to this value. This should be >= 0.")

    # parser.add_argument('--metric', type=str, default="score", required=False,
    #                     help="what to call the metric in the edge annotations") #TODO: reactivate this? or add "value" annotation in DataMatrix

    parser.add_argument('--metadata', type=str, nargs="+", required=False, default=None, 
                        help="tab separated files of sequence metadata.")

    parser.add_argument('--subset', type=str, default=None, nargs='+',
                        help="Only consider matrix labels in this list. Additive with --subset_file.")
    parser.add_argument('--subset_file', type=str, default=None,
                        help="text file containing matrix labels to retain.")
    
    parser.add_argument('--color_by', type=str, required=False, default=None,
                        help="Color the points in the output image based on this column of the metadata table.")
    parser.add_argument("--color_table", required=False, default=None, type=str, help="tab separated file with two columns and no header, columns are: annotation, hex color. For example: CCDB   cc0000")
    parser.add_argument("--color_table_out", required=False, default=None, type=str, help="tab separated file with two columns and no header, columns are: annotation, hex color. Written after the color table is updated with new colors, for example if using --color_by, but not supplying an external color table.")

    parser.add_argument('--xgmml', type=str, required=True, default=None, #TODO: what other kinds of output might be useful?
                        help="write a cytoscape xgmml file of the projection.")

    parser.add_argument('--cluster', required=False, default=False, action="store_true",
                        help="Creates a new metadata column called 'SSN_cluster', replacing that column from the metadata file if already present. Finds clusters by following connectivity, assigns each cluster a distinct integer, based on size rank, bigger clusters get smaller numbers.")

    parser.add_argument('--cluster_tsv', required=False, type=str, default=None,
                        help="Path to write a tsv file will mapping each sequence to its cluster. If set, then --cluster is automatically activated. If --no_cluster_header not set, then the file will have header contig cluster.")

    parser.add_argument('--no_cluster_header', action="store_true", required=False, default=False,
                        help="If set, then the tsv file will not have a header. Only relevant if --cluster_tsv is set.")

    add_max_output_gb_argument(parser)

    sparsify_group = parser.add_mutually_exclusive_group(required=False)
    sparsify_group.add_argument('--mst', action="store_true", required=False, default=False,
                                help="If set, then only include edges that are part of the maximum spanning tree of the graph. " \
                                "The clusters will be the same as the full graph, but the intra-cluster connections will be pruned to " \
                                "the minimum necessary to preserve the clusters.")
    sparsify_group.add_argument('--mst_knn', type=_mst_knn_arg, required=False, default=None,
                                help="Include the maximum spanning tree plus OR-symmetric k-nearest-neighbor edges, where K must be greater than 1.")

    parser.add_argument('--config', action=ActionConfigFile)


    # TODO: color by connected component (use a single color for all singletons)
    # TODO: generate MSAs for connected clusters (write to files and make a cytoscape column with the alinged sequences)
    # TODO: generate hmms for connected clusters (write to files)
    # TODO: 
    # parser.add_argument('--gradient', type=str, required=False, default=None,
    #                     help="Color nodes based on a gradient instead of discrete value mapping.")


    params = parser.parse_args(argv)
    max_output_bytes = max_output_gb_to_bytes(params.max_output_gb)


    input_file = params.input
    
    cluster = False
    if params.cluster:
        cluster=True
    if params.cluster_tsv:
        cluster=True

    color_table = None
    if params.color_table:
        color_table = read_color_table(params.color_table)

    subset_labels = list_and_file_to_dict_keys(params.subset, params.subset_file, as_set=True)
    
    # Run
    # params.metric,
    try:
        load_lower_bound = params.lb if (params.mst_knn is None and not params.mst) else None
        build_ssn(DataMatrix.from_file(input_file, lower_bound=load_lower_bound), params.lb, params.metadata, params.color_by, color_table,
                  params.xgmml, cluster, params.cluster_tsv, params.no_cluster_header, params.color_table_out,
                  params.mst, params.mst_knn, subset_labels, max_output_bytes=max_output_bytes)
    except OutputSizeLimitExceeded as exc:
        raise SystemExit(str(exc)) from None


def _entrypoint():
    main(sys.argv[1:])
    
if __name__ == '__main__':
    main(sys.argv[1:])
