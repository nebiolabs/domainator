"""Generates Sequence Similarity Networks
    
Build a sequence similarity network and do analysis related to that.
"""

from jsonargparse import ActionConfigFile
from concurrent.futures import ThreadPoolExecutor
from os import PathLike
import os
from domainator.cli.common import add_input, add_output, check_outputs, declare_outputs, tool_parser
from domainator.cytoscape import write_cytoscape_xgmml
from domainator.data_matrix import DataMatrix
from domainator.ssn_edges import MaxTree, mst_edge_index_dict
from domainator.sparsify import add_sparsify_arguments, check_sparsify_data_type, sparsify_edge_index_dict
from domainator.output_guardrails import add_max_output_gb_argument, max_output_gb_to_bytes, make_temporary_output_path
import scipy.sparse
from scipy.sparse.csgraph import connected_components
import pandas as pd
from pathlib import Path
from typing import Collection, List, Optional, Union, Dict
import numpy as np
from domainator.colors import get_palette, sort_palette_values
from domainator.cli.common import list_and_file_to_dict_keys
from domainator.colors import read_color_table

SCORE_COLUMN="SSN_SCORE"
CLUSTER_COLUMN="SSN_cluster"
MST_COLUMN="MST"


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


def cluster_labels_from_edge_dict(edge_dict, n_nodes: int) -> np.ndarray:
    """Return size-ranked connected component labels from a set of selected edges.

    Needed for a pure --knn graph, whose components are its own rather than the input
    graph's: dropping the spanning tree can leave a cluster split into several pieces.
    """
    if len(edge_dict) == 0:
        return rename_labels_by_frequency(np.arange(n_nodes, dtype=int)) + 1

    pairs = np.fromiter((index for pair in edge_dict for index in pair), dtype=int,
                        count=2 * len(edge_dict)).reshape(-1, 2)
    graph = scipy.sparse.coo_array(
        (np.ones(len(pairs), dtype=np.int8), (pairs[:, 0], pairs[:, 1])),
        shape=(n_nodes, n_nodes),
    )
    labels = connected_components(graph, directed=False, return_labels=True)[1]
    return rename_labels_by_frequency(labels) + 1


def _normalized_pair(source_idx, target_idx):
    """Key an undirected edge the way the edge index dicts in ssn_edges do."""
    return (source_idx, target_idx) if source_idx < target_idx else (target_idx, source_idx)


def iter_default_ssn_edges(matrix: DataMatrix, lb: float, mst_edges=None):
    """Yield default-path SSN edges as index/index/score tuples in legacy output order.

    mst_edges, when given, is a container of normalized (low_idx, high_idx) pairs;
    each edge then carries a trailing bool saying whether it is an MST edge.
    """
    for source_idx, target_idx, score in matrix.triangular_iter(skip_zeros=True, agg=max, index_style="index"):
        if score <= lb or source_idx == target_idx:
            continue
        if matrix.rows[target_idx] < matrix.rows[source_idx]:
            source_idx, target_idx = target_idx, source_idx
        if mst_edges is None:
            yield source_idx, target_idx, score
        else:
            yield source_idx, target_idx, score, _normalized_pair(source_idx, target_idx) in mst_edges


def iter_mst_ssn_edges(tree: "MaxTree", rows, lb: float, mark_mst: bool = False):
    """Yield MST edges above lb as (source_name, target_name, score) in canonical order.

    Every edge on this path is an MST edge, so mark_mst just appends a constant True.
    """
    for source_idx, target_idx, score in tree.iter_mst_edges():
        if score <= lb:
            break
        if source_idx == target_idx:
            continue
        source, target = rows[source_idx], rows[target_idx]
        if target < source:
            source, target = target, source
        if mark_mst:
            yield source, target, score, True
        else:
            yield source, target, score


def iter_edge_dict_ssn_edges(edge_dict, rows, mst_edges=None):
    """Yield --knn (optionally plus --mst) edges as (source_name, target_name, score) in canonical order.

    mst_edges, when given, is a container of normalized (low_idx, high_idx) pairs;
    each edge then carries a trailing bool separating the MST backbone from the
    kNN edges added around it.
    """
    for (source_idx, target_idx), score in edge_dict.items():
        source, target = rows[source_idx], rows[target_idx]
        if target < source:
            source, target = target, source
        if mst_edges is None:
            yield source, target, score
        else:
            yield source, target, score, _normalized_pair(source_idx, target_idx) in mst_edges


def build_ssn(matrix: DataMatrix, *, threshold:float=0, metadata:List[Union[str, PathLike]]=None, color_by:str=None, color_table:Dict[str,str]=None, 
              xgmml_out:Union[str, PathLike]=None, cluster:bool=False, cluster_tsv_out:Union[str, PathLike]=None, no_cluster_header:bool=False, color_table_out:str = None, mst:bool = False, knn: int = None, subset_labels: Optional[Collection[str]] = None, max_output_bytes: Optional[int] = None, max_color_groups:int=None, mark_mst:bool=False) -> None:
    """build a sequence similarity network from a matrix

    Args:
        matrix (DataMatrix): symmetric matrix of edge values
        threshold (float, optional): ignore scores lower than or equal to this. Defaults to 0.
        metadata (List[Union[str, PathLike]], optional): tab separated text files containg node annotations. Defaults to None.
        color_by (str, optional): color the nodes based on their value in this metadata field. Set to SSN_cluster to color by cluster. Defaults to None.
        xgmml_out (Union[str, PathLike], optional): write a Cytoscape XGMML file to this path. Defaults to None.
        cluster (bool, optional): If True, then add new annotations with connected cluster IDs. Defaults to False.
        cluster_tsv_out (Union[str, PathLike], optional): write a metadata tab-separated file with the cluster ids to this path. Defaults to None.
        no_cluster_header (bool, optional): If True, then the cluster tsv file will not have a header. Defaults to False.
        color_table_out (str, optional): write a color table to this path. Defaults to None.
        mst (bool, optional): If True, emit only maximum spanning tree edges, which preserves the input's connected components. Defaults to False.
        knn (int, optional): If set, emit OR-symmetric k-nearest-neighbor edges. Combined with mst=True the result is the union of the two; on its own the emitted graph need not preserve the input's connected components. Defaults to None.
        subset_labels (set, optional): keep only rows and columns with labels in this set. Defaults to None.
        max_color_groups (int, optional): if set, only this many of the largest color_by groups get distinct colors, the rest share a neutral gray. Defaults to None, meaning every group gets a color.
        mark_mst (bool, optional): If True, give every edge a boolean MST property saying whether it belongs to the maximum spanning tree. Defaults to False.

    Raises:
        ValueError: if the input does not have symmetric row and column labels
    """
    if knn is not None and knn < 1:
        raise ValueError("knn must be >= 1")
    check_sparsify_data_type(matrix.data_type, mst, knn)

    matrix = subset_matrix_by_labels(matrix, subset_labels)

    if not matrix.symmetric_labels:
        raise ValueError("Input does not have symmetric axis labels. Can only build an SSN from a symmetric matrix.")
    matrix.require_similarity("build_ssn")

    if cluster_tsv_out is not None:
        cluster = True
    # --mst (with or without --knn) changes which edges are emitted while preserving the
    # input's components, so its clusters are read off the tree. --mark_mst only annotates
    # edges, so it must keep the graph-derived clustering: the two agree on the partition,
    # but rename_labels_by_frequency can number equal-sized clusters differently depending
    # on which one produced the raw labels.
    tree = MaxTree(matrix) if (mst or mark_mst) else None

    # A pure --knn graph has components of its own, so they have to be read off the
    # selected edges -- computed here, while matrix.data is still around.
    knn_edge_dict = None
    if knn is not None and not mst:
        knn_edge_dict = sparsify_edge_index_dict(matrix, mst=False, knn=knn, lower_bound=threshold)

    node_data = pd.DataFrame(index=matrix.rows)
    if cluster:
        if mst:
            node_data[CLUSTER_COLUMN] = cluster_labels_from_tree(tree, threshold)
        elif knn_edge_dict is not None:
            node_data[CLUSTER_COLUMN] = cluster_labels_from_edge_dict(knn_edge_dict, len(matrix.rows))
        else:
            node_data[CLUSTER_COLUMN] = cluster_labels_from_graph(matrix, threshold)

    # Determine the edge source for the XGMML write.  For --mst/--knn we free
    # the large matrix.data as soon as it is no longer needed.
    edge_column_names = [SCORE_COLUMN]
    edge_metadata_types = {SCORE_COLUMN: "real"}
    if mark_mst:
        edge_column_names.append(MST_COLUMN)
        edge_metadata_types[MST_COLUMN] = "boolean"

    if mst and knn is None:
        rows = matrix.rows
        del matrix  # free O(n²) data array; rows is a plain list
        # every edge on this path is an MST edge
        xgmml_edge_rows = iter_mst_ssn_edges(tree, rows, threshold, mark_mst=mark_mst)
        xgmml_edge_kwargs = dict(
            edge_rows=xgmml_edge_rows,
            edge_column_names=edge_column_names,
            edge_metadata_types=edge_metadata_types,
            edge_endpoint_style="name",
        )
    elif knn is not None:
        # mst + knn is the union of the two; knn on its own was selected above.
        edge_dict = (sparsify_edge_index_dict(matrix, mst=True, knn=knn, lower_bound=threshold, tree=tree)
                     if mst else knn_edge_dict)
        mst_edges = set(mst_edge_index_dict(tree, lower_bound=threshold)) if mark_mst else None
        rows = matrix.rows
        del matrix  # free O(n²) data array
        xgmml_edge_rows = iter_edge_dict_ssn_edges(edge_dict, rows, mst_edges=mst_edges)
        xgmml_edge_kwargs = dict(
            edge_rows=xgmml_edge_rows,
            edge_column_names=edge_column_names,
            edge_metadata_types=edge_metadata_types,
            edge_endpoint_style="name",
        )
    else:
        mst_edges = set(mst_edge_index_dict(tree, lower_bound=threshold)) if mark_mst else None
        xgmml_edge_kwargs = dict(
            edge_rows=iter_default_ssn_edges(matrix, threshold, mst_edges=mst_edges),
            edge_column_names=edge_column_names,
            edge_metadata_types=edge_metadata_types,
            edge_endpoint_style="index",
        )

    if metadata is not None:
        with ThreadPoolExecutor() as exe:
            dfs = list(exe.map(lambda f: pd.read_csv(f, sep="\t", index_col=0), metadata))
        for df in dfs:
            node_data = node_data.merge(df, how="left", left_index=True, right_index=True)

    if color_by is not None:
        if color_table is None:
            color_table = get_palette(node_data[color_by], max_groups=max_color_groups)

    if xgmml_out is not None:
        temp_xgmml_path = make_temporary_output_path(xgmml_out)
        try:
            write_cytoscape_xgmml(
                temp_xgmml_path,
                node_data,
                **xgmml_edge_kwargs,
                name=Path(xgmml_out).stem,
                color_by=color_by,
                color_table=color_table,
                max_output_bytes=max_output_bytes,
                output_description=f"XGMML network output '{xgmml_out}'",
                mitigation_options=["--threshold", "--mst", "--knn", "--subset"],
            )
            os.replace(temp_xgmml_path, xgmml_out)
            temp_xgmml_path = None
        finally:
            if temp_xgmml_path is not None and Path(temp_xgmml_path).exists():
                os.unlink(temp_xgmml_path)

    if cluster_tsv_out:
        index_label = None
        header = None
        if not no_cluster_header:
            header = ["cluster"]
            index_label = "contig"

        node_data[CLUSTER_COLUMN].to_csv(cluster_tsv_out, index_label=index_label, header=header, sep="\t")
    
    if color_table_out is not None:
        if color_table is None:
            raise ValueError("--color_table_out requires either --color_by or --color_table, otherwise there are no colors to write.")
        with open(color_table_out, "w") as out:
            for value in sort_palette_values(color_table.keys()):
                out.write(f"{value}\t{color_table[value]}\n")


def build_parser():
    parser = tool_parser(__doc__)
    
    add_input(parser, "pairwise similarity (or distance) scores. In the network, nodes will be derived from row and column names, and edges will be added between pairs of rows and columns meeting the threshold criteria. Format can be tab-separated text, or Domainator hdf5.", multiple=False, stdin=False)

    parser.add_argument('--threshold', type=float, default=0, required=False,
                        help="keep only edges with weights strictly above this value. This should be >= 0.")

    parser.add_argument('--metadata', type=str, nargs="+", required=False, default=None, 
                        help="tab separated files of sequence metadata.")

    parser.add_argument('--subset', type=str, default=None, nargs='+',
                        help="Only consider matrix labels in this list. Additive with --subset_file.")
    parser.add_argument('--subset_file', type=str, default=None,
                        help="text file containing matrix labels to retain.")
    
    parser.add_argument('--color_by', type=str, required=False, default=None,
                        help="Color the points in the output image based on this column of the metadata table.")
    parser.add_argument("--color_table", required=False, default=None, type=str, help="tab separated file with two columns and no header, columns are: annotation, hex color. For example: CCDB   cc0000")
    parser.add_argument("--max_color_groups", required=False, default=None, type=int, help="Give distinct colors only to this many of the largest --color_by groups, all remaining groups share a neutral gray. By default every group gets a color, cycling over a palette of 64 distinguishable colors.")
    parser.add_argument("--mark_mst", action="store_true", default=False, help=f"Give every edge a boolean '{MST_COLUMN}' property saying whether it belongs to the maximum spanning tree. Useful for styling or filtering the MST backbone within a larger network in Cytoscape.")
    parser.add_argument("--color_table_out", required=False, default=None, type=str, help="tab separated file with two columns and no header, columns are: annotation, hex color. Written after the color table is updated with new colors, for example if using --color_by, but not supplying an external color table.")

    add_output(parser, "xgmml", "Write the network to this Cytoscape xgmml file.", primary=True, required=True, stdout=False)

    parser.add_argument('--cluster', required=False, default=False, action="store_true",
                        help="Creates a new metadata column called 'SSN_cluster', replacing that column from the metadata file if already present. Finds clusters by following connectivity, assigns each cluster a distinct integer, based on size rank, bigger clusters get smaller numbers.")

    parser.add_argument('--cluster_tsv_out', required=False, type=str, default=None,
                        help="Path to write a tsv file mapping each sequence to its cluster. If set, then --cluster is automatically activated. If --no_cluster_header is not set, then the file will have the header: contig cluster.")

    parser.add_argument('--no_cluster_header', action="store_true", required=False, default=False,
                        help="If set, then the tsv file will not have a header. Only relevant if --cluster_tsv_out is set.")

    add_max_output_gb_argument(parser)

    add_sparsify_arguments(parser)

    parser.add_argument('--config', action=ActionConfigFile)


    # TODO: color by connected component (use a single color for all singletons)
    # TODO: generate MSAs for connected clusters (write to files and make a cytoscape column with the alinged sequences)
    # TODO: generate hmms for connected clusters (write to files)
    # TODO: 
    declare_outputs(parser, ["xgmml_out"])
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)
    check_outputs(parser, params)
    max_output_bytes = max_output_gb_to_bytes(params.max_output_gb)


    input_file = params.input
    
    cluster = False
    if params.cluster:
        cluster=True
    if params.cluster_tsv_out:
        cluster=True

    color_table = None
    if params.color_table:
        color_table = read_color_table(params.color_table)

    subset_labels = list_and_file_to_dict_keys(params.subset, params.subset_file, as_set=True)
    
    # Run
    # params.metric,
    # The MST paths need the unthresholded matrix; --threshold is applied later by the tree
    # walk. Pure --knn and the default path can prune on load and save the memory.
    load_lower_bound = None if params.mst else params.threshold
    build_ssn(DataMatrix.from_file(input_file, lower_bound=load_lower_bound), threshold=params.threshold, metadata=params.metadata, color_by=params.color_by, color_table=color_table,
              xgmml_out=params.xgmml_out, cluster=cluster, cluster_tsv_out=params.cluster_tsv_out, no_cluster_header=params.no_cluster_header, color_table_out=params.color_table_out,
              mst=params.mst, knn=params.knn, subset_labels=subset_labels, max_output_bytes=max_output_bytes,
              max_color_groups=params.max_color_groups, mark_mst=params.mark_mst)
