"""Generates trees from distance or similarity matrices

    Distance matrices (data types ending in _dist, or dense text with no data type) are used
    as is. Similarity matrices (e.g. score, efi_score) are first converted to distances with
    1 - (score / min(row_max, col_max)), so a pair with no score has distance 1.
    The two directions of each pair are reconciled by taking the smaller distance.
"""
import warnings
warnings.filterwarnings("ignore", module='numpy')
from jsonargparse import ActionConfigFile
from domainator.cli.common import add_input, add_output, check_outputs, tool_parser
from domainator.build_projection import write_cytoscape_xgmml
from domainator.data_matrix import DataMatrix, is_distance_data_type
from domainator.transform_matrix import score_dist
import numpy as np
import pandas as pd
from pathlib import Path
import warnings
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform


SCORE_COLUMN="TREE_DIST"

# TODO: instructions for how to get cytoscape formatted files to look nice (y-files tree layout and some others?).
#see: https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format
#see also: https://github.com/biopython/biopython/blob/master/Bio/Phylo/TreeConstruction.py
#see also: https://dendropy.org/library/phylogeneticdistance.html?highlight=upgma#dendropy.calculate.phylogeneticdistance.PhylogeneticDistanceMatrix.upgma_tree
# https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html


# TODO: add options for coloring the xgmml file, like domainator build_ssn

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


DISTANCE_CONVERSIONS = ("auto", "precomputed", "score_dist")


def condensed_distances(data_matrix: DataMatrix, distance: str = "auto") -> np.ndarray:
    """
        Return a scipy condensed distance vector for a square matrix.

        distance:
            precomputed: the values are distances.
            score_dist: the values are similarity scores; distance = 1 - (score / min(row_max, col_max)).
            auto: precomputed if the data type is a distance (ends in _dist) or not recorded
                  (e.g. dense text), otherwise score_dist.

        Each pair uses the smaller of its two directional distances (the higher similarity),
        the same convention as the SSN tools.
    """
    if distance not in DISTANCE_CONVERSIONS:
        raise ValueError(f"distance must be one of {DISTANCE_CONVERSIONS}, not {distance!r}")
    if not data_matrix.symmetric_labels:
        raise ValueError("Input does not have symmetric axis labels. Can only build a tree from a square matrix comparing a set of sequences to itself.")
    if distance == "auto":
        distance = "precomputed" if (data_matrix.data_type == "" or is_distance_data_type(data_matrix.data_type)) else "score_dist"

    if distance == "precomputed":
        dense = data_matrix.toarray().astype(np.float64, copy=False)
    else:
        dense = score_dist(data_matrix.data)
        # A row or column with no scores at all gives 0/0; it is unrelated to everything.
        dense[np.isnan(dense)] = 1.0
    dense = np.minimum(dense, dense.T)
    return squareform(dense, checks=False)


def build_tree(data_matrix, algorithm, metadata_files, xgmml_file, newick_file, quiet=False, distance="auto"):
    """ 
        Build a tree from a distance or similarity matrix.

        data_matrix: a DataMatrix object
        algorithm: the algorithm to use to build the tree. Currently only supports "upgma"
        metadata_files: metadata files to add to the XGMML file
        xgmml_file: path to write the XGMML file to
        newick_file: path to write the newick file to
        distance: how to get distances from the matrix values, see condensed_distances
    """

    if algorithm == "upgma":
        linkage_matrix = hierarchy.average(condensed_distances(data_matrix, distance))
    
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


def build_parser():
    parser = tool_parser(__doc__)
    
    add_input(parser, "Square matrix of pairwise distances or similarity scores. Format can be tab-separated text or Domainator hdf5.", multiple=False, stdin=False)

    parser.add_argument('--distance', type=str, required=False, default="auto", choices=DISTANCE_CONVERSIONS,
                        help="How to turn matrix values into distances. precomputed: the values are distances. score_dist: the values are similarity scores, distance = 1 - (score / min(row_max, col_max)), so a pair with no score has distance 1. auto: precomputed if the matrix data type is a distance (ends in _dist) or not recorded (e.g. tab-separated text), otherwise score_dist.")

    parser.add_argument('--algorithm', type=str, required=False, default="upgma", choices={"upgma"},
                        help="Which clustering algorithm to use.")

    parser.add_argument('--metadata', type=str, nargs="+", required=False, default=None, 
                        help="Tab separated files of sequence metadata. These data will be added to the xgmml output.")
    
    add_output(parser, "newick", "Write the tree to this Newick file.", primary=True)
    add_output(parser, "xgmml", "Write the tree to this Cytoscape xgmml file.", stdout=False)
    
    parser.add_argument('--config', action=ActionConfigFile)
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)
    check_outputs(parser, params, ["newick_out", "xgmml_out"])
    build_tree(DataMatrix.from_file(params.input), params.algorithm, params.metadata, params.xgmml_out, params.newick_out, params.quiet, params.distance)
