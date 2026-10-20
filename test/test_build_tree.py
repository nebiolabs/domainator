import pytest
import tempfile
from domainator.build_tree import get_newick
from domainator import build_tree
from helpers import compare_files

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
    result = get_newick(tree, tree.dist, leaf_names)
    assert result == '((D:1.00,C:1.00):3.00,(B:2.00,A:2.00):2.00);'

def test_get_newick_single_node():
    single_node_tree = TreeNode(0, 0)
    leaf_names = ['A']
    result = get_newick(single_node_tree, 0, leaf_names)
    assert result == '(A:0.00);'

def test_get_newick_single_level(leaf_names):
    root = TreeNode(-1, 3)
    root.left = TreeNode(0, 0)
    root.right = TreeNode(1, 0)
    result = get_newick(root, root.dist, leaf_names)
    assert result == '(B:3.00,A:3.00);'

def test_get_newick_with_empty_leaf_names(tree, leaf_names):
    result = get_newick(tree, tree.dist, [''] * len(leaf_names))
    assert result == '((:1.00,:1.00):3.00,(:2.00,:2.00):2.00);'


def test_newick_output(shared_datadir):

    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        newick_out = str(output_dir + "/test.newick")
        xgmm_out = str(output_dir + "/test.xgmml")
        metadata = str(shared_datadir / 'FeSOD_metadata.tsv')

        build_tree.main(['--input', str(shared_datadir / 'FeSOD_score_dist.tsv'), '--newick', newick_out, '--xgmml', xgmm_out, "--metadata", metadata])
        compare_files(newick_out, str(shared_datadir / 'FeSOD_score_dist.newick'))
        compare_files(xgmm_out, str(shared_datadir / 'FeSOD_score_dist.xgmml'))


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
    build_tree.main(['--input', str(hdf_path), '--newick', str(newick_out), '--quiet'])
    compare_files(str(newick_out), str(shared_datadir / 'FeSOD_score_dist.newick'))
