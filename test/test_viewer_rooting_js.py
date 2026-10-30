"""Topology invariants for the viewer's iterative JavaScript rooting algorithms."""

import json
from io import StringIO
import shutil
import subprocess

import pytest
from Bio import Phylo

from domainator.viewer_html import _script

pytestmark = pytest.mark.skipif(shutil.which("node") is None, reason="node is not installed")

_RUN = """
const input = JSON.parse(process.argv[2]);
const state = {bundle: {nodes: input.view.leaf_nodes.map((_, i) => String(i)), views: {tree: input.view}}};
const nodeId = i => state.bundle.nodes[i];
state.tree = createTreeState(state.bundle);
const results = [];
for (const action of input.actions) {
    try {
        const old = state.tree;
        const selected = new Set(action.selection || []);
        const clade = treeSelectedClade(old, selected);
        if (action.mode === 'selection' && clade === null) { throw new Error('not monophyletic'); }
        const location = action.mode === 'original' ? null : action.mode === 'selection'
            ? treeSelectionRootLocation(old, clade) : treeMidpointLocation(old.rootOriginal);
        const result = treeRootedTopology(old.rootOriginal, location);
        state.bundle.views.tree = result.view;
        state.tree = createTreeState(state.bundle);
        state.tree.rootOriginal = old.rootOriginal;
        state.tree.rootLocation = location;
        state.tree.rootOrigins = result.origins;
        results.push({view: result.view, location, newick: treeNewick(null),
            distances: Array.from(state.tree.rootDistance)});
    } catch (error) {
        results.push({error: error.message});
    }
}
process.stdout.write(JSON.stringify(results));
"""


def _run(tmp_path, view, actions):
    script = tmp_path / "rooting.js"
    script.write_text(_script("views/tree/tree.js") + _script("views/tree/rooting.js") + _RUN)
    result = subprocess.run(
        ["node", str(script), json.dumps({"view": view, "actions": actions})],
        capture_output=True, text=True, check=True,
    )
    return json.loads(result.stdout)


def _view(parent, length, support=None, labels=None):
    parents = set(parent)
    return {
        "parent": parent,
        "length": length,
        "support": support or [None] * len(parent),
        "label": labels or [None] * len(parent),
        "leaf_nodes": [node for node in range(len(parent)) if node not in parents],
    }


def _distances(view):
    graph = [[] for _ in view["parent"]]
    for node, parent in enumerate(view["parent"][1:], 1):
        length = view["length"][node]
        graph[node].append((parent, length))
        graph[parent].append((node, length))
    out = {}
    for index, leaf in enumerate(view["leaf_nodes"]):
        distances = {}
        stack = [(leaf, -1, 0)]
        while stack:
            node, parent, distance = stack.pop()
            distances[node] = distance
            stack.extend((other, node, distance + length) for other, length in graph[node] if other != parent)
        for other, target in enumerate(view["leaf_nodes"]):
            out[index, other] = distances[target]
    return out


def _supported_splits(view):
    descendants = [set() for _ in view["parent"]]
    for index, node in enumerate(view["leaf_nodes"]):
        descendants[node].add(index)
    for node in range(len(view["parent"]) - 1, 0, -1):
        descendants[view["parent"][node]].update(descendants[node])
    all_leaves = descendants[0]
    return sorted(
        (support, frozenset((frozenset(descendants[node]), frozenset(all_leaves - descendants[node]))))
        for node, support in enumerate(view["support"]) if node and support is not None
    )


@pytest.mark.parametrize("view", [
    _view([-1, 0, 1, 1, 0], [0, 2, 1, 3, 8], [None, 90, None, None, None], [None, "clade", None, None, None]),
    _view([-1, 0, 0], [0, 2, 2]),
    _view([-1, 0, 0, 0], [0, 0, 3, 7]),
    _view([-1, 0, 1, 1], [0, 9, 2, 4], labels=["stem", None, None, None]),
    _view([-1, 0, 1, 2, 2, 1, 0], [0, 1, 2, 3, 0, 4, 6]),
])
def test_midpoint_preserves_distances_splits_and_original(tmp_path, view):
    results = _run(tmp_path, view, [{"mode": "midpoint"}, {"mode": "midpoint"}, {"mode": "original"}])
    first, second, restored = results
    assert "error" not in first, first
    assert first == second
    assert restored["view"] == view
    assert _distances(first["view"]) == pytest.approx(_distances(view))
    assert _supported_splits(first["view"]) == _supported_splits(view)
    parsed = Phylo.read(StringIO(first["newick"]), "newick")
    assert sorted(leaf.name for leaf in parsed.get_terminals()) == [str(i) for i in range(len(view["leaf_nodes"]))]
    for (left, right), distance in _distances(view).items():
        assert parsed.distance(str(left), str(right)) == pytest.approx(distance, rel=1e-5, abs=1e-6)
    diameter = max(_distances(view).values())
    assert max(first["distances"]) == pytest.approx(diameter / 2)


@pytest.mark.parametrize("selection,valid", [
    ([], False), ([0], True), ([0, 1], True), ([0, 2], False), ([0, 1, 2], False),
])
def test_selection_requires_exact_current_root_clade(tmp_path, selection, valid):
    view = _view([-1, 0, 1, 1, 0], [0, 2, 1, 3, 8])
    result, = _run(tmp_path, view, [{"mode": "selection", "selection": selection}])
    if not valid:
        assert result == {"error": "not monophyletic"}
    else:
        assert "error" not in result, result
        assert _distances(result["view"]) == pytest.approx(_distances(view))
        assert result["view"]["parent"].count(0) == 2


def test_repeated_selection_roots_do_not_accumulate_nodes(tmp_path):
    view = _view([-1, 0, 1, 1, 0], [0, 0, 1, 3, 8])
    actions = [{"mode": "selection", "selection": [0]}] * 25
    results = _run(tmp_path, view, actions)
    for result in results:
        assert "error" not in result, result
        assert len(result["view"]["parent"]) <= len(view["parent"]) + 1
        assert _distances(result["view"]) == pytest.approx(_distances(view))


def test_zero_length_selection_boundary_is_rootable(tmp_path):
    view = _view([-1, 0, 0], [0, 0, 0])
    rooted, restored = _run(tmp_path, view, [{"mode": "selection", "selection": [0]}, {"mode": "original"}])
    assert "error" not in rooted
    assert _distances(rooted["view"]) == _distances(view)
    assert restored["view"] == view


def test_large_finite_lengths_do_not_overflow_edge_midpoints(tmp_path):
    view = _view([-1, 0, 0], [0, 1e308, 0])
    actions = [{"mode": "midpoint"}] + [{"mode": "selection", "selection": [0]}] * 5
    for result in _run(tmp_path, view, actions):
        assert "error" not in result, result
        assert _distances(result["view"]) == pytest.approx(_distances(view))


def test_named_original_root_is_not_lost_from_full_newick(tmp_path):
    view = _view([-1, 0, 1, 1, 0], [0, 2, 1, 3, 8], labels=["named_root", "clade", None, None, None])
    result, = _run(tmp_path, view, [{"mode": "midpoint"}])
    assert "named_root" in result["newick"]
    assert "clade" in result["newick"]


@pytest.mark.parametrize("view,reason", [
    (_view([-1], [0]), "at least two leaves"),
    (_view([-1, 0, 0], [0, 0, 0]), "positive tree diameter"),
    (_view([-1, 0, 0], [0, -1, 2]), "nonnegative"),
    (_view([-1, 0, 1, 1, 0], [0, 1e308, 1e308, 1, 1]), "overflow"),
])
def test_midpoint_refuses_undefined_distances(tmp_path, view, reason):
    result, = _run(tmp_path, view, [{"mode": "midpoint"}])
    assert reason in result["error"]


def test_deep_tree_uses_iterative_traversals(tmp_path):
    depth = 5000
    parent = [-1] + list(range(depth)) + [depth]
    view = _view(parent, [0] + [1] * (len(parent) - 1))
    # Two leaves below a long unary stem.
    view["parent"].append(depth)
    view["length"].append(3)
    view["label"].append(None)
    view["support"].append(None)
    view["leaf_nodes"] = [depth + 1, depth + 2]
    result, = _run(tmp_path, view, [{"mode": "midpoint"}])
    assert "error" not in result, result
    assert len(result["view"]["leaf_nodes"]) == 2
    assert max(result["distances"]) == 2


def test_deep_branching_tree_does_not_need_recursion(tmp_path):
    depth = 2000
    parent = [-1] + list(range(depth)) + list(reversed(range(depth)))
    view = _view(parent, [0] + [1] * (len(parent) - 1))
    result, = _run(tmp_path, view, [{"mode": "midpoint"}])
    assert "error" not in result, result
    assert len(result["view"]["leaf_nodes"]) == depth + 1
    assert len(result["view"]["parent"]) <= len(view["parent"]) + 1
    assert max(result["distances"]) == (depth + 1) / 2
