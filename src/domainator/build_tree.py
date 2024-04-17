"""Generates Trees from distance matrices
    
"""

import argparse
from os import PathLike
import sys
from domainator.build_projection import write_cytoscape_xgmml
from domainator.data_matrix import DataMatrix
import pandas as pd
from pathlib import Path
from domainator import __version__, RawAndDefaultsFormatter
from typing import List, Union
import warnings
from scipy.cluster import hierarchy
from domainator.build_projection import write_cytoscape_xgmml
from scipy.spatial.distance import squareform


SCORE_COLUMN="TREE_DIST"

# TODO: instructions for how to get cytoscape formatted files to look nice (y-files tree layout and some others?).
#see: https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format
#see also: https://github.com/biopython/biopython/blob/master/Bio/Phylo/TreeConstruction.py
#see also: https://dendropy.org/library/phylogeneticdistance.html?highlight=upgma#dendropy.calculate.phylogeneticdistance.PhylogeneticDistanceMatrix.upgma_tree
# https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html


def clean_name_for_newick(name, quiet=False):
  """
    converts all: " " (space), ";" (semicolon), ":" (colon), "," (comma), "()" (parentheses), "'" (quote) characters to "_" in a string
  """

  bad_chars = " ;:,()'"
  chars = ["_" if x in bad_chars else x for x in name]
  out = "".join(chars)
  if out != name:
    if not quiet:
        warnings.warn(f"sequence name {name} contains characters that don't play well with newick, so for newick export, it has been changed to {out}")
  return out

def get_tree_edges(node, parent_dist, leaf_names, parent_name=None, edges=None, internal_nodes=None) -> dict: #adapted from: https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format
    """
    Convert scipy.cluster.hierarchy.to_tree()-output to an edge list

    node: output of scipy.cluster.hierarchy.to_tree()
    parent_dist: output of scipy.cluster.hierarchy.to_tree().dist
    leaf_names: list of leaf names
    parent_name: the name of the parent node
    edges: leave empty, this variable will be populated by an empty dict, and is used in recursion.
    
    returns:
        dictionary of edges: [(start,end)]:dist
        internal_nodes: names of internal nodes
    """
    #TODO: make this not recursive

    if edges is None:
        edges = dict()
        internal_nodes = list()

    if node.is_leaf():
        name = leaf_names[node.id]
        edges[(parent_name,name)] = parent_dist - node.dist
    else:
        name = f"__INTERNAL_NODE_{node.id}__"
        internal_nodes.append(name)
        if parent_name is not None:
            edges[(parent_name,name)] = parent_dist - node.dist
        get_tree_edges(node.get_left(), node.dist, leaf_names, parent_name=name, edges=edges, internal_nodes=internal_nodes)
        get_tree_edges(node.get_right(), node.dist, leaf_names, parent_name=name, edges=edges, internal_nodes=internal_nodes)
        return edges, internal_nodes


def get_newick(node, parent_dist, leaf_names, newick='') -> str: #from: https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format
    """ 
    Convert scipy.cluster.hierarchy.to_tree()-output to Newick format.

    node: output of scipy.cluster.hierarchy.to_tree()
    parent_dist: output of scipy.cluster.hierarchy.to_tree().dist
    leaf_names: list of leaf names
    newick: leave empty, this variable is used in recursion.
    
    returns:
        tree in Newick format
    """
    #TODO: make this not recursive

    if len(leaf_names) == 1:
        return f"({leaf_names[0]}:0.00);" #if there's only one leaf, the tree is just a leaf 

    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parent_dist - node.dist, newick)
        else:
            newick = ");"
        newick = get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
        newick = get_newick(node.get_right(), node.dist, leaf_names, newick=",%s" % (newick))
        newick = "(%s" % (newick)
        return newick


def build_tree(data_matrix, algorithm, metadata_files, xgmml_file, newick_file, quiet=False):
    """ 
        Build a tree from a distance matrix.

        data_matrix: a DataMatrix object
        algorithm: the algorithm to use to build the tree. Currently only supports "upgma"
        metadata_files: metadata files to add to the XGMML file
        xgmml_file: path to write the XGMML file to
        newick_file: path to write the newick file to
    """

    if data_matrix.sparse:
        warnings.warn(f"Input is a sparse matrix. Distance matrices are usually not sparse. You may want to rerun seq_dist with --mode score_dist and output to --dense or --dense_text.")
        data_matrix.convert_to_dense()
    if algorithm == "upgma":
        linkage_matrix = hierarchy.average(squareform(data_matrix.data, checks=False))
    
    tree = hierarchy.to_tree(linkage_matrix, False)

    if newick_file is not None:
        names = list(map(lambda x: clean_name_for_newick(x, quiet), data_matrix.rows))

        newick_string = get_newick(tree, tree.dist, names)

        with open(newick_file, "w") as outfile:
            print(newick_string, file=outfile)

    if xgmml_file is not None:
        edge_data, internal_nodes = get_tree_edges(tree, tree.dist, data_matrix.rows)
        edge_data = pd.DataFrame.from_dict(edge_data, orient="index", columns=[SCORE_COLUMN])

        node_data = pd.DataFrame()
        node_data.index = data_matrix.rows + internal_nodes
        if metadata_files is not None:
            for file in metadata_files:
                metadata = pd.read_csv(file, sep="\t", index_col=0)
                node_data = node_data.merge(metadata,how="left", left_index=True, right_index=True)
        write_cytoscape_xgmml(xgmml_file, node_data, edge_data, name=Path(xgmml_file).stem)


def main(argv):
    parser = argparse.ArgumentParser(description=f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)
    
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-i', '--input', type=str, required=False,
                        help="pairwise distances. Format can be tab-separated text, or Domainator dense hdf5.")

    parser.add_argument('--algorithm', type=str, required=False, default="upgma", choices={"upgma"},
                        help="Which clustering algorithm to use.")

    parser.add_argument('--metadata', type=str, nargs="+", required=False, default=None, 
                        help="Tab separated files of sequence metadata. These data will be added to the xgmml output.")
    
    parser.add_argument('--xgmml', type=str, required=False, default=None, 
                        help="write a cytoscape xgmml file of the tree.")
    
    parser.add_argument('--newick', type=str, required=False, default=None, #TODO: what other kinds of output might be useful?
                        help="write a newick file of the tree.")
    
    parser.add_argument('--quiet', action="store_true", required=False, default=False, 
                        help="Suppress warnings.")

    # parser.add_argument('--svg', type=str, required=True, default=None, #TODO: what other kinds of output might be useful?
    #                     help="write an image of the tree in svg format.")

    # parser.add_argument('--png', type=str, required=True, default=None, #TODO: what other kinds of output might be useful?
    #                     help="write an image of the tree in png format.")

    params = parser.parse_args(argv)


    input_file = params.input
    

    # Run
    # params.metric,
    build_tree(DataMatrix.from_file(input_file), params.algorithm, params.metadata, params.xgmml, params.newick, params.quiet)

def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])