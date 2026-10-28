"""Generates trees from distance or similarity matrices, or reads them from Newick files

    Distance matrices (data types ending in _dist, or dense text with no data type) are used
    as is. Similarity matrices (e.g. score, efi_score) are first converted to distances with
    1 - (score / min(row_max, col_max)), so a pair with no score has distance 1.
    The two directions of each pair are reconciled by taking the smaller distance.

    Instead of a matrix, --newick reads a tree built elsewhere (for example a maximum-likelihood
    tree), to write it with metadata, colors, and an alignment.

    Outputs are a Newick tree, a Cytoscape xgmml file (laid out as a rectangular tree), and an
    interactive Domainator Viewer page, which can show a multiple sequence alignment beside the
    tree.
"""
from dataclasses import dataclass
from os import PathLike
from typing import Dict, List, Optional, Set, Union
import warnings
from jsonargparse import ActionConfigFile
from domainator.cli.common import add_input, add_output, check_outputs, declare_outputs, list_and_file_to_dict_keys, tool_parser
from domainator.build_ssn import subset_matrix_by_labels
from domainator.colors import get_palette, read_color_table, sort_palette_values
from domainator.cytoscape import write_cytoscape_xgmml
from domainator.errors import InputError
from domainator.outputs import output_file
from domainator.output_guardrails import DEFAULT_MAX_OUTPUT_GB, add_max_output_gb_argument, max_output_gb_to_bytes
from domainator.data_matrix import DataMatrix, is_distance_data_type
from domainator.seqio import parse_seqfiles
from domainator.transform_matrix import score_dist
from domainator.viewer_bundle import read_node_metadata, viewer_bundle, viewer_defaults
from domainator.viewer_html import write_viewer_outputs
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform


SCORE_COLUMN="TREE_DIST"

# The xgmml layout: the deepest leaf is this far from the root, and leaves are this far apart.
XGMML_TREE_WIDTH = 1000
XGMML_LEAF_SPACING = 30
_X_COLUMN = "__tree_x__"
_Y_COLUMN = "__tree_y__"

# Characters other than these (and gaps) in at least this fraction of an alignment's residues
# make it protein rather than nucleotide.
_NUCLEOTIDE_LETTERS = set("ACGTUN")
_NUCLEOTIDE_MIN_FRACTION = 0.9
_GAP_CHARACTERS = set("-.")


# TODO: instructions for how to get cytoscape formatted files to look nice (y-files tree layout and some others?).
#see: https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format
#see also: https://github.com/biopython/biopython/blob/master/Bio/Phylo/TreeConstruction.py
#see also: https://dendropy.org/library/phylogeneticdistance.html?highlight=upgma#dendropy.calculate.phylogeneticdistance.PhylogeneticDistanceMatrix.upgma_tree
# https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html


@dataclass
class Tree:
    """A rooted tree as arrays over its nodes in preorder, so node 0 is the root and every
    node comes after its parent; children are in the order they were written.

    parent: each node's parent, -1 for the root.
    length: each node's branch length to its parent (0 for the root).
    name: each leaf's name; an internal node's label, or None.
    support: each internal node's support value, or None.
    """
    parent: List[int]
    length: List[float]
    name: List[Optional[str]]
    support: List[Optional[float]]

    def children(self) -> List[List[int]]:
        children = [[] for _ in self.parent]
        for node, parent in enumerate(self.parent):
            if parent >= 0:
                children[parent].append(node)
        return children

    def leaves(self) -> List[int]:
        """The leaf nodes, in preorder (the order the tree draws them)."""
        children = self.children()
        return [node for node in range(len(self.parent)) if not children[node]]

    def leaf_names(self) -> List[str]:
        return [self.name[node] for node in self.leaves()]


def _tree_from_linkage(root, leaf_names) -> Tree:
    """The tree of scipy.cluster.hierarchy.to_tree's ``root``, each node's right subtree before
    its left one."""
    tree = Tree([], [], [], [])
    stack = [(root, -1, root.dist)]
    while stack:
        node, parent, parent_dist = stack.pop()
        index = len(tree.parent)
        tree.parent.append(parent)
        tree.length.append(parent_dist - node.dist if parent >= 0 else 0.0)
        tree.name.append(leaf_names[node.id] if node.is_leaf() else None)
        tree.support.append(None)
        if not node.is_leaf():
            stack.append((node.get_left(), index, node.dist))
            stack.append((node.get_right(), index, node.dist))
    return tree


def _tree_from_newick(path: Union[str, PathLike]) -> Tree:
    """The first tree in a Newick file.

    Numeric internal labels are read as support values. A tree with no branch lengths at all
    gets length 1 everywhere, so it draws as a cladogram; otherwise a missing length is 0.
    """
    from Bio import Phylo
    try:
        phylo_tree = Phylo.read(path, "newick")
    except Exception as error:
        raise InputError(f"Could not read a Newick tree from {path}: {error}") from error
    clades = list(phylo_tree.find_clades(order="preorder"))
    has_lengths = any(clade.branch_length is not None for clade in clades[1:])
    index_of = {id(clade): i for i, clade in enumerate(clades)}
    tree = Tree([], [], [], [])
    for i, clade in enumerate(clades):
        tree.parent.append(-1)
        if i == 0:
            tree.length.append(0.0)
        elif has_lengths:
            tree.length.append(float(clade.branch_length or 0.0))
        else:
            tree.length.append(1.0)
        tree.name.append(clade.name)
        tree.support.append(None if clade.is_terminal() or clade.confidence is None else float(clade.confidence))
    for clade in clades:
        for child in clade.clades:
            tree.parent[index_of[id(child)]] = index_of[id(clade)]
    names = tree.leaf_names()
    if any(name is None or name == "" for name in names):
        raise InputError(f"Every leaf of the tree in {path} needs a name.")
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise InputError(f"Leaf names in {path} must be unique; repeated: {', '.join(duplicates[:10])}")
    return tree


def subset_tree(tree: Tree, keep: Set[str]) -> Tree:
    """The tree of the leaves named in ``keep``: other leaves are removed, then internal nodes
    left with no children, and an internal node left with one child is merged into it (adding
    their branch lengths)."""
    children = tree.children()
    kept = [False] * len(tree.parent)
    for node in reversed(range(len(tree.parent))):
        kept[node] = (tree.name[node] in keep) if not children[node] else any(kept[child] for child in children[node])
    if not kept[0]:
        raise InputError("None of the tree's leaves are in the subset.")

    out = Tree([], [], [], [])
    # (node, new parent, length carried down from merged single-child ancestors)
    stack = [(0, -1, 0.0)]
    while stack:
        node, parent, carried = stack.pop()
        kept_children = [child for child in children[node] if kept[child]]
        length = carried + (tree.length[node] if node else 0.0)
        if len(kept_children) == 1 and children[node]:
            stack.append((kept_children[0], parent, length))
            continue
        index = len(out.parent)
        out.parent.append(parent)
        out.length.append(length if parent >= 0 else 0.0)
        out.name.append(tree.name[node])
        out.support.append(tree.support[node])
        for child in reversed(kept_children):
            stack.append((child, index, 0.0))
    return out


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


def get_newick(tree: Tree, quiet: bool = False) -> str:
    """The tree in Newick format, with branch lengths to six significant digits.

    Internal nodes are written with their labels, or failing that their support values.
    """
    children = tree.children()
    if len(children[0]) == 0:
        return f"({clean_name_for_newick(tree.name[0], quiet)}:0);"  # if there's only one leaf, the tree is just a leaf

    out = list()
    stack = [0]  # a node to write, or a string to write as is
    while stack:
        item = stack.pop()
        if isinstance(item, str):
            out.append(item)
            continue
        node = item
        branch = "" if node == 0 else ":%.6g" % tree.length[node]
        if not children[node]:
            out.append(clean_name_for_newick(tree.name[node], quiet) + branch)
            continue
        if tree.name[node] is not None:
            label = clean_name_for_newick(tree.name[node], quiet)
        elif tree.support[node] is not None:
            label = "%.6g" % tree.support[node]
        else:
            label = ""
        out.append("(")
        stack.append(")" + label + branch)
        for i, child in enumerate(reversed(children[node])):
            if i:
                stack.append(",")
            stack.append(child)
    return "".join(out) + ";"


def tree_layout(tree: Tree):
    """Rectangular layout coordinates: (x, y) for every node, x the distance from the root and y
    the leaf's place in preorder, an internal node halfway between its first and last child."""
    children = tree.children()
    x = [0.0] * len(tree.parent)
    for node in range(1, len(tree.parent)):
        x[node] = x[tree.parent[node]] + tree.length[node]
    y = [0.0] * len(tree.parent)
    leaf_rank = 0
    for node in range(len(tree.parent)):
        if not children[node]:
            y[node] = float(leaf_rank)
            leaf_rank += 1
    for node in reversed(range(len(tree.parent))):
        if children[node]:
            y[node] = (y[children[node][0]] + y[children[node][-1]]) / 2
    return x, y


def internal_node_name(node: int) -> str:
    return f"__INTERNAL_NODE_{node}__"


def _xgmml_frames(tree: Tree, node_data: pd.DataFrame):
    """The xgmml writer's node and edge tables: every node, laid out as a rectangular tree, with
    its metadata, and every branch with its length."""
    children = tree.children()
    names = [tree.name[node] if not children[node] else internal_node_name(node) for node in range(len(tree.parent))]
    edges = {(names[tree.parent[node]], names[node]): tree.length[node] for node in range(1, len(tree.parent))}
    edge_data = pd.DataFrame.from_dict(edges, orient="index", columns=[SCORE_COLUMN])
    x, y = tree_layout(tree)
    depth = max(x) or 1.0
    # Leaves first, then internal nodes.
    order = tree.leaves() + [node for node in range(len(tree.parent)) if children[node]]
    nodes = pd.DataFrame(index=[names[node] for node in order])
    nodes = nodes.merge(node_data, how="left", left_index=True, right_index=True)
    nodes[_X_COLUMN] = [x[node] * XGMML_TREE_WIDTH / depth for node in order]
    nodes[_Y_COLUMN] = [-y[node] * XGMML_LEAF_SPACING for node in order]  # the writer negates y
    return nodes, edge_data


def _guess_alphabet(sequences) -> str:
    """"nucleotide" when almost every residue is A, C, G, T, U or N, otherwise "protein"."""
    residues = 0
    nucleotides = 0
    for sequence in sequences:
        for letter in sequence.upper():
            if letter in _GAP_CHARACTERS:
                continue
            residues += 1
            nucleotides += letter in _NUCLEOTIDE_LETTERS
    return "nucleotide" if residues and nucleotides >= _NUCLEOTIDE_MIN_FRACTION * residues else "protein"


def read_msa(path: Union[str, PathLike], leaf_names: List[str]) -> dict:
    """The tree view's ``msa``: ``{alphabet, length, sequences}``, ``sequences`` holding each
    leaf's aligned sequence (None for a leaf the alignment lacks)."""
    sequences = {}
    for record in parse_seqfiles([path], filetype_override="fasta"):
        if record.id in sequences:
            raise InputError(f"The alignment {path} has more than one sequence named {record.id}.")
        sequences[record.id] = str(record.seq)
    if not sequences:
        raise InputError(f"The alignment {path} has no sequences.")
    lengths = {len(sequence) for sequence in sequences.values()}
    if len(lengths) != 1:
        raise InputError(f"The sequences in {path} are not all the same length ({min(lengths)} to {max(lengths)}), so they are not an alignment.")
    leaves = set(leaf_names)
    missing = [name for name in leaf_names if name not in sequences]
    extra = [name for name in sequences if name not in leaves]
    if missing:
        warnings.warn(f"{len(missing)} of the tree's {len(leaf_names)} leaves have no sequence in {path}, for example {missing[0]}.")
    if extra:
        warnings.warn(f"{len(extra)} sequences in {path} are not leaves of the tree and are left out, for example {extra[0]}.")
    return {
        "alphabet": _guess_alphabet(sequences.values()),
        "length": lengths.pop(),
        "sequences": [sequences.get(name) for name in leaf_names],
    }


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


def tree_from_matrix(data_matrix: DataMatrix, algorithm: str = "upgma", distance: str = "auto") -> Tree:
    """The tree of a distance or similarity matrix's rows (see condensed_distances)."""
    if algorithm != "upgma":
        raise ValueError(f"algorithm must be upgma, not {algorithm!r}")
    if len(data_matrix.rows) == 1:
        return Tree([-1], [0.0], [data_matrix.rows[0]], [None])
    linkage_matrix = hierarchy.average(condensed_distances(data_matrix, distance))
    return _tree_from_linkage(hierarchy.to_tree(linkage_matrix, False), data_matrix.rows)


def tree_view(tree: Tree, msa: Optional[dict] = None) -> dict:
    """The bundle's ``views.tree``: the tree's arrays, which node each of the bundle's nodes (the
    leaves, in preorder) is, and the alignment."""
    children = tree.children()
    return {
        "parent": tree.parent,
        "length": tree.length,
        # Leaves are named by the bundle's nodes; only internal labels are kept here.
        "label": [tree.name[node] if children[node] else None for node in range(len(tree.parent))],
        "support": tree.support,
        "leaf_nodes": tree.leaves(),
        "msa": msa,
    }


def build_tree(data_matrix: Optional[DataMatrix] = None, *, newick: Optional[Union[str, PathLike]] = None,
               algorithm: str = "upgma", distance: str = "auto", metadata: Optional[List[str]] = None,
               subset_labels: Optional[Set[str]] = None, msa: Optional[Union[str, PathLike]] = None,
               color_by: Optional[str] = None, color_table: Optional[Dict[str, str]] = None,
               color_table_out: Optional[str] = None, palette_offset: int = 0, max_color_groups: Optional[int] = None,
               label_by: Optional[str] = None, categorical: Optional[List[str]] = None, name: Optional[str] = None,
               newick_out: Optional[str] = None, xgmml_out: Optional[str] = None, html_out: Optional[str] = None,
               dvb_out: Optional[str] = None, max_output_gb: float = DEFAULT_MAX_OUTPUT_GB, quiet: bool = False) -> Tree:
    """Build a tree from a distance or similarity matrix, or read one from a Newick file, and write it.

    Args:
        data_matrix: a square distance or similarity matrix. Exactly one of data_matrix and newick.
        newick: read the tree from this Newick file instead of building one.
        algorithm: the algorithm to build the tree with. Only "upgma".
        distance: how to get distances from the matrix values, see condensed_distances.
        metadata: tab-separated tables of leaf metadata, keyed by their first column.
        subset_labels: only include these leaves (all of them if None).
        msa: an aligned FASTA file of the leaves' sequences, shown beside the tree in the viewer.
        color_by: color the leaves in the xgmml and, at first, the viewer by this metadata column.
        color_table: maps color_by values to colors. None: assign colors automatically.
        color_table_out: write the color table, with any colors assigned, to this path: two
            tab-separated columns (annotation, hex color) and no header.
        palette_offset: skip this many colors at the start of the palette when assigning colors
            to color_by values. Ignored when color_table is given.
        max_color_groups: give distinct colors only to this many of the largest color_by groups.
        label_by: the metadata column the viewer labels leaves by at first.
        categorical: numeric metadata columns the viewer colors as categories, not a gradient.
        name: the tree's name, shown in the viewer. Defaults to the html or dvb output's stem.
        newick_out: write the tree to this Newick file.
        xgmml_out: write the tree to this Cytoscape xgmml file, laid out as a rectangular tree.
        html_out: write the Domainator Viewer, opening this tree, to this HTML file.
        dvb_out: write the viewer bundle (gzip-compressed JSON, .dvb) to this file.
        max_output_gb: refuse to write an html or dvb output larger than this (0: no limit).
        quiet: do not warn about leaf names changed for Newick.

    Returns:
        The tree.
    """
    if (data_matrix is None) == (newick is None):
        raise ValueError("Give exactly one of a matrix and a Newick tree.")
    if msa is not None and html_out is None and dvb_out is None:
        raise ValueError("An alignment is only shown in the viewer, so it requires html_out or dvb_out.")
    if color_table_out is not None and color_by is None and color_table is None:
        raise ValueError("color_table_out requires either color_by or color_table, otherwise there are no colors to write.")

    if data_matrix is not None:
        tree = tree_from_matrix(subset_matrix_by_labels(data_matrix, subset_labels), algorithm, distance)
    else:
        tree = _tree_from_newick(newick)
        if subset_labels is not None:
            tree = subset_tree(tree, set(subset_labels))
    leaf_names = tree.leaf_names()

    if newick_out is not None:
        with output_file(newick_out) as outfile:
            print(get_newick(tree, quiet), file=outfile)

    node_data = read_node_metadata(metadata, leaf_names)
    if color_by is not None and color_by not in node_data.columns:
        raise ValueError(f"Requested color_by column '{color_by}' was not found in the merged metadata.")
    xgmml_colors = color_table
    if color_by is not None and xgmml_colors is None:
        xgmml_colors = get_palette(node_data[color_by], max_groups=max_color_groups, offset=palette_offset)

    if xgmml_out is not None:
        nodes, edges = _xgmml_frames(tree, node_data)
        write_cytoscape_xgmml(xgmml_out, nodes, edges, name=Path(xgmml_out).stem, x_col=_X_COLUMN, y_col=_Y_COLUMN,
                              color_by=color_by, color_table=xgmml_colors)

    if color_table_out is not None:
        with open(color_table_out, "w") as out:
            for value in sort_palette_values(xgmml_colors.keys()):
                out.write(f"{'' if value is None else value}\t{xgmml_colors[value]}\n")

    if html_out is not None or dvb_out is not None:
        if name is None:
            name = Path(html_out if html_out is not None else dvb_out).stem
        # The viewer assigns default colors as get_palette does; an offset or capped palette goes in as a table.
        viewer_colors = color_table if color_table is not None or (palette_offset == 0 and max_color_groups is None) else xgmml_colors
        defaults = viewer_defaults(node_data, color_by=color_by, color_table=viewer_colors, label_by=label_by, categorical=categorical)
        msa_view = read_msa(msa, leaf_names) if msa is not None else None
        bundle = viewer_bundle(name, leaf_names, node_data, defaults, {"tree": tree_view(tree, msa_view)})
        write_viewer_outputs(bundle, html_out=html_out, dvb_out=dvb_out, max_output_bytes=max_output_gb_to_bytes(max_output_gb),
                             mitigation_options=["--subset", "--subset_file"])
    return tree


def build_parser():
    parser = tool_parser(__doc__)

    add_input(parser, "Square matrix of pairwise distances or similarity scores. Format can be tab-separated text or Domainator hdf5. Give this or --newick.", multiple=False, required=False, stdin=False)
    parser.add_argument("--newick", type=str, required=False, default=None,
                        help="Read the tree from this Newick file instead of building one from a matrix. Numeric internal node labels are read as support values.")

    parser.add_argument('--distance', type=str, required=False, default="auto", choices=DISTANCE_CONVERSIONS,
                        help="How to turn matrix values into distances. precomputed: the values are distances. score_dist: the values are similarity scores, distance = 1 - (score / min(row_max, col_max)), so a pair with no score has distance 1. auto: precomputed if the matrix data type is a distance (ends in _dist) or not recorded (e.g. tab-separated text), otherwise score_dist.")

    parser.add_argument('--algorithm', type=str, required=False, default="upgma", choices=["upgma"],
                        help="Which clustering algorithm to use.")

    parser.add_argument('--metadata', type=str, nargs="+", required=False, default=None,
                        help="Tab separated files of sequence metadata, keyed by their first column. These data are added to the xgmml and html outputs.")
    parser.add_argument("--subset", type=str, default=None, nargs='+',
                        help="Only include these leaves. Additive with --subset_file.")
    parser.add_argument("--subset_file", type=str, default=None,
                        help="Text file of leaf names to include.")
    parser.add_argument("--msa", type=str, required=False, default=None,
                        help="An aligned FASTA file of the leaves' sequences, keyed by name, shown beside the tree in the html output.")

    parser.add_argument('--color_by', type=str, required=False, default=None,
                        help="Color the leaves in the xgmml, and at first in the html viewer, by this column of the metadata table.")
    parser.add_argument("--color_table", required=False, default=None, type=str, help="tab separated file with two columns and no header, columns are: annotation, hex color. For example: CCDB   cc0000. Requires --color_by.")
    parser.add_argument("--color_table_out", required=False, default=None, type=str, help="tab separated file with two columns and no header, columns are: annotation, hex color. Written after the color table is updated with new colors, for example if using --color_by, but not supplying an external color table.")
    parser.add_argument("--palette_offset", required=False, default=0, type=int, help="Skip this many colors at the start of the palette before assigning colors to --color_by values. Ignored with --color_table.")
    parser.add_argument("--max_color_groups", required=False, default=None, type=int, help="Give distinct colors only to this many of the largest --color_by groups, all remaining groups share a neutral gray.")
    parser.add_argument("--label_by", type=str, default=None,
                        help="Label the leaves in the html viewer by this metadata column at first.")
    parser.add_argument("--categorical", type=str, nargs="+", required=False, default=None,
                        help="Numeric metadata columns (for example integer cluster numbers) that the html viewer should color as discrete categories instead of a gradient.")
    parser.add_argument("--name", type=str, default=None,
                        help="The tree's name, shown in the html viewer. Defaults to the html (or dvb) output's file name.")

    add_output(parser, "newick", "Write the tree to this Newick file.", primary=True)
    add_output(parser, "xgmml", "Write the tree to this Cytoscape xgmml file, laid out as a rectangular tree.", stdout=False)
    add_output(parser, "html", "Write the Domainator Viewer, showing this tree, to this HTML file.", stdout=False)
    add_output(parser, "dvb", "Write the gzip-compressed JSON Domainator Viewer bundle (.dvb) to this file. The viewer opens it, and saves sessions in the same format.", stdout=False)
    add_max_output_gb_argument(parser)

    parser.add_argument('--config', action=ActionConfigFile)
    declare_outputs(parser, ["newick_out", "xgmml_out", "html_out", "dvb_out"])
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)
    check_outputs(parser, params)
    if (params.input is None) == (params.newick is None):
        parser.error("Give exactly one of -i/--input and --newick.")
    if params.newick is not None and (params.distance != "auto" or params.algorithm != "upgma"):
        parser.error("--distance and --algorithm apply to a matrix, not to a --newick tree.")
    if params.msa is not None and params.html_out is None and params.dvb_out is None:
        parser.error("--msa is shown only in the viewer, so it requires --html_out or --dvb_out.")
    if params.color_table is not None and params.color_by is None:
        parser.error("--color_table requires --color_by.")

    build_tree(DataMatrix.from_file(params.input) if params.input is not None else None,
               newick=params.newick, algorithm=params.algorithm, distance=params.distance, metadata=params.metadata,
               subset_labels=list_and_file_to_dict_keys(params.subset, params.subset_file, as_set=True), msa=params.msa,
               color_by=params.color_by, color_table=read_color_table(params.color_table) if params.color_table else None,
               color_table_out=params.color_table_out, palette_offset=params.palette_offset,
               max_color_groups=params.max_color_groups, label_by=params.label_by, categorical=params.categorical,
               name=params.name, newick_out=params.newick_out, xgmml_out=params.xgmml_out, html_out=params.html_out,
               dvb_out=params.dvb_out, max_output_gb=params.max_output_gb, quiet=params.quiet)
