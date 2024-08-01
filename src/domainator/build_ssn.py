"""Generates Sequence Similarity Networks
    
Build a sequence similarity network and do analysis related to that.
"""

# TODO: MSA column for cystoscape cluster
# TODO: nothing here is very memory efficient. How big of a problem is that?
import warnings
warnings.filterwarnings("ignore", module='numpy')
import argparse
from os import PathLike
import sys
from domainator.cytoscape import write_cytoscape_xgmml
from domainator.data_matrix import DataMatrix
from scipy.sparse.csgraph import connected_components
import pandas as pd
from pathlib import Path
from domainator import __version__, RawAndDefaultsFormatter
from typing import List, Union, Dict
import numpy as np
from domainator.utils import get_palette
from domainator.color_genbank import read_color_table


SCORE_COLUMN="SSN_SCORE"
CLUSTER_COLUMN="SSN_cluster"

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

def build_ssn(matrix: DataMatrix, lb:float=0, metadata_files:List[Union[str, PathLike]]=None, color_by:str=None, color_table:Dict[str,str]=None, 
              xgmml:Union[str, PathLike]=None, cluster:bool=False, cluster_tsv:Union[str, PathLike]=None, no_cluster_header:bool=False, color_table_out:str = None):
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


    Raises:
        ValueError: if the input does not have symmetric row and column labels
    """
    
    if not matrix.symmetric:
        raise ValueError("Input not symmetric. Can only build an SSN from a symmetric matrix.")
    
    if cluster_tsv is not None:
        cluster = True
    # find edges
    edge_data=dict() # (source, target): score
    for source,target,score in matrix.iter_data():
        if score > lb:
            if source != target: #don't write self-edges
                if target < source:
                    tmp = source
                    source = target
                    
                    target = tmp
                edge_data[(source,target)] = score

    #TODO: could save a lot of memory by streaming: writing the xgmml row by row instead of loading the entire matrix and metadata into memory.
    edge_data = pd.DataFrame.from_dict(edge_data, orient="index", columns=[SCORE_COLUMN]) # (source_node, dest_node): score
    node_data = pd.DataFrame()
    node_data.index = matrix.rows # node: []

    if cluster:
        n_components, labels =  connected_components( (matrix.data > lb), directed=False, return_labels=True)
        labels = rename_labels_by_frequency(labels)
        labels = labels + 1 # make the cluster labels start at 1 instead of 0
        node_data[CLUSTER_COLUMN] = labels

    if metadata_files is not None:
        for file in metadata_files:
            metadata = pd.read_csv(file, sep="\t", index_col=0)
            node_data = node_data.merge(metadata,how="left", left_index=True, right_index=True)
    
    if color_by is not None:
        if color_table is None:
            color_table = get_palette(node_data[color_by].unique())

    if xgmml is not None:
        write_cytoscape_xgmml(xgmml, node_data, edge_data, name=Path(xgmml).stem, color_by=color_by, color_table=color_table)

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
            color_table_items.sort(key=lambda x: x[0]) # sort lexically
            for domain, color in color_table_items:
                out.write(f"{domain}\t{color}\n")


def main(argv):
    parser = argparse.ArgumentParser(description=f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)
    
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-i', '--input', type=str, required=False,
                        help="pairwise similarity (or distance) scores. In the network, nodes will be derived from row and column names, and edges will be added between pairs of rows and columns meeting the threshold criteria. Format can be tab-separated text, or Domainator hdf5.")

    #TODO: add a separate program for matrix analysis, such as density, non-zeros, histograms like they do in EFI
    #           maybe just a 'hist' option that write a histogram and then quits?
    #TODO: tree output (or as a separate script, or as part of the ssn script?), in the ssn 
    #TODO: auto-generate MSAs and hmms and logos and stuff from clusters

    parser.add_argument('--lb', type=float, default=0, required=False,
                        help="exclude all edges with weights less than or equal to this value. This should be >= 0.")

    # parser.add_argument('--metric', type=str, default="score", required=False,
    #                     help="what to call the metric in the edge annotations") #TODO: reactivate this? or add "value" annotation in DataMatrix

    parser.add_argument('--metadata', type=str, nargs="+", required=False, default=None, 
                        help="tab separated files of sequence metadata.")
    
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




    # TODO: add annotation for connected clusters
    # TODO: color by connected cluster (use a single color for all singletons)
    # TODO: generate MSAs for connected clusters (write to files and make a cytoscape column with the alinged sequences)
    # TODO: generate hmms for connected clusters (write to files)
    # TODO: 
    # parser.add_argument('--gradient', type=str, required=False, default=None,
    #                     help="Color nodes based on a gradient instead of discrete value mapping.")


    params = parser.parse_args(argv)


    input_file = params.input
    
    cluster = False
    if params.cluster:
        cluster=True
    if params.cluster_tsv:
        cluster=True

    color_table = None
    if params.color_table:
        color_table = read_color_table(params.color_table)
    # Run
    # params.metric,
    build_ssn(DataMatrix.from_file(input_file), params.lb, params.metadata, params.color_by, color_table, params.xgmml, cluster, params.cluster_tsv, params.no_cluster_header, params.color_table_out)


def _entrypoint():
    main(sys.argv[1:])
    
if __name__ == '__main__':
    main(sys.argv[1:])
