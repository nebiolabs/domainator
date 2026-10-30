import gzip
import json
import tempfile
import warnings

import pytest
from domainator.build_tree import Tree, _tree_from_linkage, get_newick
from domainator import build_tree
from helpers import compare_files, compare_xgmml, viewer_bundle

# Simple tree structure for testing purposes
class TreeNode:
    def __init__(self, id, dist=0):
        self.id = id
        self.dist = dist
        self.left = None
        self.right = None

    def is_leaf(self):
        return self.left is None and self.right is None

    def get_left(self):
        return self.left

    def get_right(self):
        return self.right

@pytest.fixture
def leaf_names():
    return ['A', 'B', 'C', 'D', 'E']

@pytest.fixture
def tree():
    root = TreeNode(-1, 4)
    root.left = TreeNode(-1, 2)
    root.right = TreeNode(-1, 1)
    root.left.left = TreeNode(0, 0)
    root.left.right = TreeNode(1, 0)
    root.right.left = TreeNode(2, 0)
    root.right.right = TreeNode(3, 0)
    return root

def test_get_newick(tree, leaf_names):
    result = get_newick(_tree_from_linkage(tree, leaf_names))
    assert result == '((D:1,C:1):3,(B:2,A:2):2);'

def test_get_newick_single_node():
    assert get_newick(Tree([-1], [0.0], ['A'], [None])) == '(A:0);'

def test_get_newick_single_level(leaf_names):
    root = TreeNode(-1, 3)
    root.left = TreeNode(0, 0)
    root.right = TreeNode(1, 0)
    result = get_newick(_tree_from_linkage(root, leaf_names))
    assert result == '(B:3,A:3);'

def test_get_newick_with_empty_leaf_names(tree, leaf_names):
    result = get_newick(_tree_from_linkage(tree, [''] * len(leaf_names)))
    assert result == '((:1,:1):3,(:2,:2):2);'


def test_deep_trees_do_not_recurse():
    # A caterpillar tree is as deep as it has leaves, deeper than Python's recursion limit.
    n = 5000
    node = TreeNode(0, 0)
    for i in range(1, n):
        parent = TreeNode(n + i, i)
        parent.left, parent.right = node, TreeNode(i, 0)
        node = parent
    names = [f"s{i}" for i in range(n)]
    tree = _tree_from_linkage(node, names)
    assert len(tree.parent) == 2 * n - 1 and len(tree.leaves()) == n
    newick = get_newick(tree)
    assert newick.count("(") == n - 1 and newick.endswith(");")
    import pandas as pd
    nodes, edges = build_tree._xgmml_frames(tree, pd.DataFrame(index=names))
    assert len(edges) == 2 * n - 2 and len(nodes) == 2 * n - 1


def test_newick_output(shared_datadir):

    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        newick_out = str(output_dir + "/test.newick")
        xgmm_out = str(output_dir + "/test.xgmml")
        metadata = str(shared_datadir / 'FeSOD_metadata.tsv')

        build_tree.main(['--input', str(shared_datadir / 'FeSOD_score_dist.tsv'), '--newick_out', newick_out, '--xgmml_out', xgmm_out, "--metadata", metadata])
        compare_files(newick_out, str(shared_datadir / 'FeSOD_score_dist.newick'))
        compare_xgmml(xgmm_out, str(shared_datadir / 'FeSOD_score_dist.xgmml'))


def _condensed_to_square(condensed):
    from scipy.spatial.distance import squareform
    return squareform(condensed)


def test_condensed_distances_from_similarity_treats_missing_pairs_as_unrelated():
    import numpy as np
    from domainator.data_matrix import DataMatrix

    scores = np.array([
        [10.0, 8.0, 0.0],
        [6.0, 10.0, 0.0],
        [0.0, 0.0, 10.0],
    ])
    names = ["a", "b", "c"]
    matrix = DataMatrix(scores, names, names, data_type="score")
    distances = _condensed_to_square(build_tree.condensed_distances(matrix))
    # a-b: the better direction (8/10) is used; a-c and b-c have no score, so distance 1.
    np.testing.assert_allclose(distances, [[0.0, 0.2, 1.0], [0.2, 0.0, 1.0], [1.0, 1.0, 0.0]])


def test_condensed_distances_precomputed_uses_smaller_distance_of_each_pair():
    import numpy as np
    from domainator.data_matrix import DataMatrix

    distances = np.array([
        [0.0, 0.3, 0.9],
        [0.1, 0.0, 0.5],
        [0.7, 0.6, 0.0],
    ])
    names = ["a", "b", "c"]
    expected = [[0.0, 0.1, 0.7], [0.1, 0.0, 0.5], [0.7, 0.5, 0.0]]
    for data_type in ("", "score_dist"):  # auto picks precomputed for both
        matrix = DataMatrix(distances, names, names, data_type=data_type)
        np.testing.assert_allclose(_condensed_to_square(build_tree.condensed_distances(matrix)), expected)
    # explicit override on a matrix whose data type says similarity
    matrix = DataMatrix(distances, names, names, data_type="score")
    np.testing.assert_allclose(_condensed_to_square(build_tree.condensed_distances(matrix, "precomputed")), expected)


def test_condensed_distances_requires_matching_labels():
    import numpy as np
    from domainator.data_matrix import DataMatrix

    matrix = DataMatrix(np.eye(2), ["a", "b"], ["x", "y"], data_type="score")
    with pytest.raises(ValueError, match="symmetric axis labels"):
        build_tree.condensed_distances(matrix)


def test_build_tree_from_hdf_distance_matrix_matches_text(shared_datadir, tmp_path):
    # A distance matrix stored in hdf5 drops its zeros (identical pairs); they must come back.
    from domainator.data_matrix import DataMatrix

    text_matrix = DataMatrix.from_file(shared_datadir / 'FeSOD_score_dist.tsv')
    hdf_path = tmp_path / "score_dist.hdf5"
    DataMatrix.write_hdf(text_matrix.data, hdf_path, text_matrix.rows, text_matrix.columns, data_type="score_dist")
    newick_out = tmp_path / "from_hdf.newick"
    build_tree.main(['--input', str(hdf_path), '--newick_out', str(newick_out), '--quiet'])
    compare_files(str(newick_out), str(shared_datadir / 'FeSOD_score_dist.newick'))


def _read_bundle(path):
    with gzip.open(path, "rt") as handle:
        return json.load(handle)


def test_newick_input_round_trips(shared_datadir, tmp_path):
    out = tmp_path / "out.newick"
    build_tree.main(["--newick", str(shared_datadir / "tree_support.newick"), "--newick_out", str(out)])
    assert out.read_text().strip() == (shared_datadir / "tree_support.newick").read_text().strip().replace("G:1.0", "G:1")


def test_newick_subset_prunes_leaves_and_merges_single_child_nodes(shared_datadir, tmp_path):
    out = tmp_path / "out.newick"
    build_tree.main(["--newick", str(shared_datadir / "tree_support.newick"), "--newick_out", str(out),
                     "--subset", "A", "C", "D", "G"])
    # (A,B)90 loses B, so A takes its branch; clade_x keeps only (C,D)75, which takes its branch.
    assert out.read_text().strip() == "(A:0.4,(C:0.1,D:0.15)75:0.25,G:1);"


def test_newick_without_branch_lengths_draws_as_a_cladogram(tmp_path):
    newick = tmp_path / "plain.newick"
    newick.write_text("((A,B),C);\n")
    tree = build_tree._tree_from_newick(newick)
    assert tree.length == [0.0, 1.0, 1.0, 1.0, 1.0]


@pytest.mark.parametrize("text, message", [
    ("((A:1,A:1):1,B:1);", "must be unique"),
    ("((A:1,:1):1,B:1);", "needs a name"),
])
def test_newick_leaf_names_must_be_unique_and_present(tmp_path, text, message):
    from domainator.errors import InputError
    newick = tmp_path / "bad.newick"
    newick.write_text(text)
    with pytest.raises(InputError, match=message):
        build_tree._tree_from_newick(newick)


def test_bundle_holds_the_tree_metadata_and_alignment(shared_datadir, tmp_path):
    dvb = tmp_path / "tree.dvb"
    build_tree.main(["-i", str(shared_datadir / "FeSOD_score_dist.tsv"), "--dvb_out", str(dvb),
                     "--metadata", str(shared_datadir / "FeSOD_metadata.tsv"), "--msa", str(shared_datadir / "FeSOD_20.afa"),
                     "--color_by", "Status", "--label_by", "Accession", "--categorical", "active", "--name", "FeSOD"])
    bundle = _read_bundle(dvb)
    assert bundle["format"] == "domainator_viewer_bundle" and bundle["version"] == 2
    assert bundle["name"] == "FeSOD"
    assert bundle["defaults"] == {"color_by": "Status", "label_by": "Accession", "categorical_columns": ["active"]}
    view = bundle["views"]["tree"]
    tree = build_tree.tree_from_matrix(build_tree.DataMatrix.from_file(shared_datadir / "FeSOD_score_dist.tsv"))
    # The bundle's nodes are the leaves in the order the tree draws them.
    assert bundle["nodes"] == tree.leaf_names()
    assert view["parent"] == tree.parent
    assert [view["label"][node] for node in view["leaf_nodes"]] == [None] * 20
    assert all(view["parent"][node] < node for node in range(1, len(view["parent"])))
    # Metadata rows follow the nodes.
    columns = [column["name"] for column in bundle["metadata"]["columns"]]
    accession = columns.index("Accession")
    assert [row[accession] for row in bundle["metadata"]["rows"]] == [name.split("|")[0].removeprefix("FeSOD_") for name in bundle["nodes"]]
    msa = view["msa"]
    assert msa["alphabet"] == "protein" and msa["length"] == 234
    fasta = {record.id: str(record.seq) for record in build_tree.parse_seqfiles([shared_datadir / "FeSOD_20.afa"], filetype_override="fasta")}
    assert msa["sequences"] == [fasta[name] for name in bundle["nodes"]]


def test_newick_bundle_keeps_internal_labels_and_support(shared_datadir, tmp_path):
    dvb = tmp_path / "tree.dvb"
    build_tree.build_tree(newick=shared_datadir / "tree_support.newick", dvb_out=str(dvb))
    view = _read_bundle(dvb)["views"]["tree"]
    assert view["label"] == [None, None, None, None, "clade_x", None, None, None, None, None, None]
    assert view["support"][1] == 90 and view["support"][5] == 75
    assert view["msa"] is None


def test_html_embeds_the_bundle(shared_datadir, tmp_path):
    html = tmp_path / "tree.html"
    build_tree.main(["--newick", str(shared_datadir / "tree_support.newick"), "--html_out", str(html),
                     "--metadata", str(shared_datadir / "tree_support_metadata.tsv")])
    page = html.read_text()
    assert "<title>tree</title>" in page and ">Domainator Viewer: tree</h1>" in page
    bundle = viewer_bundle(page)
    assert bundle["nodes"] == list("ABCDEFG") and "tree" in bundle["views"]


def test_alignment_rows_follow_the_leaves(tmp_path):
    msa = tmp_path / "aln.fasta"
    msa.write_text(">B\nAC-GT\n>A\nACGGT\n>X\nTTTTT\n")
    with pytest.warns(UserWarning) as record:
        result = build_tree.read_msa(msa, ["A", "B", "C"])
    messages = " ".join(str(warning.message) for warning in record)
    assert "1 of the tree's 3 leaves have no sequence" in messages
    assert "1 sequences in" in messages and "are left out" in messages
    assert result == {"alphabet": "nucleotide", "length": 5, "sequences": ["ACGGT", "AC-GT", None]}


@pytest.mark.parametrize("text, message", [
    (">A\nACGT\n>B\nACG\n", "not all the same length"),
    (">A\nACGT\n>A\nACGA\n", "more than one sequence named A"),
])
def test_alignment_errors(tmp_path, text, message):
    from domainator.errors import InputError
    msa = tmp_path / "aln.fasta"
    msa.write_text(text)
    with pytest.raises(InputError, match=message):
        build_tree.read_msa(msa, ["A", "B"])


def test_xgmml_colors_leaves_and_lays_out_the_tree(shared_datadir, tmp_path):
    import xml.etree.ElementTree as ET
    from domainator.colors import OTHER_COLOR, get_palette
    xgmml = tmp_path / "tree.xgmml"
    table = tmp_path / "colors.tsv"
    build_tree.main(["--newick", str(shared_datadir / "tree_support.newick"), "--xgmml_out", str(xgmml),
                     "--metadata", str(shared_datadir / "tree_support_metadata.tsv"), "--color_by", "group",
                     "--color_table_out", str(table)])
    ns = {"x": "http://www.cs.rpi.edu/XGMML"}
    nodes = {}
    for node in ET.parse(xgmml).getroot().findall("x:node", ns):
        name = next(att.get("value") for att in node.findall("x:att", ns) if att.get("name") == "Description")
        nodes[name] = node.find("x:graphics", ns)
    import pandas as pd
    palette = get_palette(pd.read_csv(shared_datadir / "tree_support_metadata.tsv", sep="\t", index_col=0)["group"])
    assert nodes["A"].get("fill") == palette["alpha"]
    assert nodes["G"].get("fill") == palette["beta"]
    assert nodes["F"].get("fill") == OTHER_COLOR
    assert nodes["__INTERNAL_NODE_0__"].get("fill") == OTHER_COLOR
    # Leaves are spaced evenly down the page in preorder; the root is at x 0, the deepest leaf at the far right.
    ys = [float(nodes[leaf].get("y")) for leaf in "ABCDEFG"]
    assert ys == sorted(ys) and len(set(ys)) == 7
    assert float(nodes["__INTERNAL_NODE_0__"].get("x")) == 0
    assert max(float(graphics.get("x")) for graphics in nodes.values()) == pytest.approx(build_tree.XGMML_TREE_WIDTH)
    rows = [line.split("\t") for line in table.read_text().splitlines()]
    assert rows == [["alpha", palette["alpha"]], ["beta", palette["beta"]], ["gamma", palette["gamma"]], ["", OTHER_COLOR]]


@pytest.mark.parametrize("argv, message", [
    ([], "exactly one of"),
    (["-i", "MATRIX", "--newick", "NEWICK"], "exactly one of"),
    (["--newick", "NEWICK", "--distance", "score_dist"], "apply to a matrix"),
    (["--newick", "NEWICK", "--msa", "MSA", "--newick_out", "OUT"], "requires --html_out or --dvb_out"),
    (["--newick", "NEWICK", "--color_table", "TABLE"], "--color_table requires --color_by"),
])
def test_argument_errors(shared_datadir, tmp_path, capsys, argv, message):
    replacements = {"MATRIX": str(shared_datadir / "FeSOD_score_dist.tsv"), "NEWICK": str(shared_datadir / "tree_support.newick"),
                    "MSA": str(shared_datadir / "FeSOD_20.afa"), "OUT": str(tmp_path / "out.nwk"),
                    "TABLE": str(shared_datadir / "color_table_123.tsv")}
    argv = [replacements.get(arg, arg) for arg in argv]
    if "--newick_out" not in argv:
        argv += ["--html_out", str(tmp_path / "out.html")]
    with pytest.raises(SystemExit):
        build_tree.main(argv)
    assert message in capsys.readouterr().err
