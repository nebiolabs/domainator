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
