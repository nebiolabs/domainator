"""Browser-driven tests of the Domainator Viewer's tree view (build_tree --html_out).

Like test_ssn_viewer_browser.py, these load the self-contained page in headless Chromium via
Playwright and skip when it is not installed (pip install -e ".[test,browser]" and playwright
install chromium).
"""

import gzip
import json
import shutil

import numpy as np
import pandas as pd
import pytest

from domainator import build_ssn, build_tree
from domainator.colors import get_palette
from domainator.data_matrix import DataMatrix
from domainator.seqio import parse_seqfiles

sync_api = pytest.importorskip("playwright.sync_api")


def _loaded_page(html_path, viewport=None):
    with sync_api.sync_playwright() as p:
        try:
            browser = p.chromium.launch()
        except sync_api.Error as exc:
            pytest.skip(f"Chromium not available for Playwright; run 'playwright install chromium'. ({exc})")
        page = browser.new_page(viewport=viewport or {"width": 1400, "height": 1100}, accept_downloads=True)
        page.pageerrors = []
        page.requests = []
        page.on("pageerror", lambda exc: page.pageerrors.append(str(exc)))
        page.on("request", lambda request: page.requests.append(request.url) if not request.url.startswith(("file:", "data:", "blob:")) else None)
        page.goto(html_path.as_uri())
        page.wait_for_function("() => document.getElementById('bundle-status').textContent.startsWith('Loaded')")
        try:
            yield page
        finally:
            browser.close()


@pytest.fixture(scope="module")
def small_html(tmp_path_factory, shared_datadir_module):
    out = tmp_path_factory.mktemp("tree_small") / "small.html"
    build_tree.main(["--newick", str(shared_datadir_module / "tree_support.newick"), "--html_out", str(out),
                     "--metadata", str(shared_datadir_module / "tree_support_metadata.tsv"), "--color_by", "group"])
    return out


@pytest.fixture(scope="module")
def msa_html(tmp_path_factory, shared_datadir_module):
    out = tmp_path_factory.mktemp("tree_msa") / "fesod.html"
    build_tree.main(["-i", str(shared_datadir_module / "FeSOD_score_dist.tsv"), "--html_out", str(out),
                     "--metadata", str(shared_datadir_module / "FeSOD_metadata.tsv"), "--msa", str(shared_datadir_module / "FeSOD_20.afa"),
                     "--quiet"])
    return out


@pytest.fixture(scope="module")
def shared_datadir_module(tmp_path_factory, request):
    """test/data, copied once for the module's page builders."""
    source = request.path.parent / "data"
    target = tmp_path_factory.mktemp("data")
    for name in ("tree_support.newick", "tree_support_metadata.tsv", "FeSOD_score_dist.tsv", "FeSOD_metadata.tsv", "FeSOD_20.afa"):
        shutil.copy(source / name, target / name)
    return target


@pytest.fixture
def page(small_html):
    yield from _loaded_page(small_html)


@pytest.fixture
def msa_page(msa_html):
    yield from _loaded_page(msa_html)


def _leaf_point(page, name):
    """Where leaf `name`'s tip is on the screen, in client coordinates."""
    return page.evaluate("""name => {
        const rect = document.getElementById('tree-canvas').getBoundingClientRect();
        const g = treeGeometry(rect.width, rect.height);
        const node = state.tree.view.leaf_nodes[state.bundle.nodes.indexOf(name)];
        return {x: rect.left + g.xPx(node), y: rect.top + g.yPx(state.tree.layout.nodeRow[node])};
    }""", name)


def _node_point(page, node):
    return page.evaluate("""node => {
        const rect = document.getElementById('tree-canvas').getBoundingClientRect();
        const g = treeGeometry(rect.width, rect.height);
        return {x: rect.left + g.xPx(node), y: rect.top + g.yPx(state.tree.layout.nodeRow[node])};
    }""", node)


def _selected(page):
    return sorted(page.evaluate("Array.from(state.selectedNodeIndices).map(i => state.bundle.nodes[i])"))


def _pixel(page, canvas_id, x, y):
    """The color of the canvas pixel at client point (x, y), as "#rrggbb"."""
    return page.evaluate("""([id, x, y]) => {
        const canvas = document.getElementById(id);
        const rect = canvas.getBoundingClientRect();
        const ratio = canvas.width / rect.width;
        const [r, g, b] = canvas.getContext('2d').getImageData(Math.round((x - rect.left) * ratio), Math.round((y - rect.top) * ratio), 1, 1).data;
        return '#' + [r, g, b].map(v => v.toString(16).padStart(2, '0')).join('').toUpperCase();
    }""", [canvas_id, x, y])


def test_tree_page_opens_the_tree_view(page):
    assert page.evaluate("document.body.dataset.activeView") == "tree"
    assert page.evaluate("state.tree.layout.rowNode.length") == 7
    assert page.is_visible("#tree-canvas")
    assert not page.is_visible("#cluster-view") and not page.is_visible("#split-chart")
    assert not page.is_visible("#msa-pane")
    assert page.locator("#metadata-table tbody tr").count() == 7
    assert page.pageerrors == [] and page.requests == []


def test_click_leaf_and_clade(page):
    point = _leaf_point(page, "C")
    page.mouse.click(point["x"] + 15, point["y"])
    assert _selected(page) == ["C"]
    assert page.locator("#metadata-table tbody tr").count() == 1
    clade = page.evaluate("state.tree.view.parent[state.tree.view.leaf_nodes[state.bundle.nodes.indexOf('C')]]")
    point = _node_point(page, clade)
    page.mouse.click(point["x"], point["y"])
    assert _selected(page) == ["C", "D"]
    page.mouse.click(point["x"], point["y"])
    assert _selected(page) == []
    assert page.pageerrors == []


def test_double_click_collapses_and_expand_all_restores(page):
    clade_x = page.evaluate("state.tree.view.label.indexOf('clade_x')")
    point = _node_point(page, clade_x)
    page.mouse.dblclick(point["x"], point["y"])
    assert page.evaluate("Array.from(state.tree.collapsed)") == [clade_x]
    assert page.evaluate("state.tree.layout.rowNode.length") == 4
    # The collapsed clade's row toggles all of its leaves.
    point = _node_point(page, clade_x)
    page.mouse.click(point["x"] + 30, point["y"])
    assert _selected(page) == ["C", "D", "E", "F"]
    page.click("#tree-expand-all")
    assert page.evaluate("state.tree.layout.rowNode.length") == 7


def test_shift_drag_selects_and_alt_deselects(page):
    a, b = _leaf_point(page, "A"), _leaf_point(page, "B")
    page.keyboard.down("Shift")
    page.mouse.move(a["x"] - 30, a["y"] - 6)
    page.mouse.down()
    page.mouse.move(b["x"] + 40, b["y"] + 6, steps=4)
    page.mouse.up()
    page.keyboard.up("Shift")
    assert _selected(page) == ["A", "B"]
    page.keyboard.down("Shift")
    page.keyboard.down("Alt")
    page.mouse.move(a["x"] - 30, a["y"] - 6)
    page.mouse.down()
    page.mouse.move(a["x"] + 40, a["y"] + 6, steps=4)
    page.mouse.up()
    page.keyboard.up("Alt")
    page.keyboard.up("Shift")
    assert _selected(page) == ["B"]


def test_leaves_take_the_color_by_palette(page, shared_datadir_module):
    palette = get_palette(pd.read_csv(shared_datadir_module / "tree_support_metadata.tsv", sep="\t", index_col=0)["group"])
    for leaf, group in (("A", "alpha"), ("G", "beta"), ("E", "gamma")):
        point = _leaf_point(page, leaf)
        assert _pixel(page, "tree-canvas", point["x"], point["y"]) == palette[group].upper()
    page.select_option("#color-by", "")
    point = _leaf_point(page, "A")
    assert _pixel(page, "tree-canvas", point["x"], point["y"]) != palette["alpha"].upper()


def test_ladderize_puts_smaller_clades_first(page):
    leaf_order = "state.tree.layout.rowNode.map(node => state.bundle.nodes[state.tree.treeToBundle[node]])"
    assert page.evaluate(leaf_order) == list("ABCDEFG")
    page.check("#tree-ladderize")
    assert page.evaluate(leaf_order) == ["G", "A", "B", "E", "F", "C", "D"]


def test_cut_fills_cluster_numbers_and_exports_them(page):
    page.evaluate("document.getElementById('tree-cut').value = 250; document.getElementById('tree-cut').dispatchEvent(new Event('input'))")
    assert page.text_content("#tree-cut-label") == "0.25"
    page.click("#metadata-panel-add")
    page.click("#metadata-add-cluster-column")
    column = page.evaluate("state.metadataColumnIndexByName.get('tree_cluster')")
    numbers = page.evaluate(f"state.metadataByNodeIndex.map(row => row[{column}])")
    assert numbers == [1, 1, 2, 2, 3, 4, 5]
    with page.expect_download() as download:
        page.click("#export-selected")
    lines = open(download.value.path()).read().splitlines()
    assert lines[0].split("\t")[:2] == ["node_id", "tree_cluster"]
    assert [line.split("\t")[1] for line in lines[1:]] == ["1", "1", "2", "2", "3", "4", "5"]


def test_newick_export_matches_build_tree(page, shared_datadir_module, tmp_path):
    with page.expect_download() as download:
        page.click("#export-newick")
    assert download.value.suggested_filename == "small.nwk"
    expected = build_tree.get_newick(build_tree._tree_from_newick(shared_datadir_module / "tree_support.newick"))
    assert open(download.value.path()).read().strip() == expected
    for leaf in ("A", "C", "D", "G"):
        point = _leaf_point(page, leaf)
        page.mouse.click(point["x"] + 15, point["y"])
    with page.expect_download() as download:
        page.click("#export-newick")
    assert open(download.value.path()).read().strip() == "(A:0.4,(C:0.1,D:0.15)75:0.25,G:1);"


def test_view_exports(page):
    with page.expect_download() as download:
        page.click("#export-png")
    assert download.value.suggested_filename == "small_view.png"
    assert open(download.value.path(), "rb").read(8) == b"\x89PNG\r\n\x1a\n"
    with page.expect_download() as download:
        page.click("#export-svg")
    svg = open(download.value.path()).read()
    assert svg.startswith("<svg") and ">G</text>" in svg and "90</text>" in svg


def test_session_round_trip(page, small_html):
    clade_x = page.evaluate("state.tree.view.label.indexOf('clade_x')")
    point = _node_point(page, clade_x)
    page.mouse.dblclick(point["x"], point["y"])
    point = _leaf_point(page, "G")
    page.mouse.click(point["x"] + 15, point["y"])
    page.evaluate("document.getElementById('tree-cut').value = 500; document.getElementById('tree-cut').dispatchEvent(new Event('input'))")
    page.check("#tree-ladderize")
    with page.expect_download() as download:
        page.click("#save-session")
    assert download.value.suggested_filename == "small_session.dvb"
    session = download.value.path()
    bundle = json.loads(gzip.open(session).read())
    assert bundle["format"] == "domainator_viewer_bundle" and bundle["app_state"]["view"]["tree_collapsed"] == [clade_x]

    page.goto(small_html.as_uri())
    page.wait_for_function("() => document.getElementById('bundle-status').textContent.startsWith('Loaded')")
    page.set_input_files("#bundle-file", session)
    page.wait_for_function("() => state.tree && state.tree.collapsed.size === 1")
    assert _selected(page) == ["G"]
    assert page.input_value("#tree-cut") == "500"
    assert page.is_checked("#tree-ladderize")
    assert page.evaluate("state.tree.layout.rowNode.length") == 4
    assert "skipped" not in page.text_content("#bundle-status")
    assert page.pageerrors == []


def test_alignment_rows_line_up_with_the_tree(msa_page, shared_datadir_module):
    page = msa_page
    assert page.is_visible("#msa-pane") and page.is_visible("#export-fasta")
    sequences = {record.id: str(record.seq) for record in parse_seqfiles([shared_datadir_module / "FeSOD_20.afa"], filetype_override="fasta")}
    rows = page.evaluate("state.tree.layout.rowNode.map(node => [state.bundle.nodes[state.tree.treeToBundle[node]], msaRowSequence(node)])")
    assert [sequence for _, sequence in rows] == [sequences[name] for name, _ in rows]
    # The cell of a row's first residue has that residue's color, at the row's height in the tree.
    name = rows[3][0]
    column = next(i for i, letter in enumerate(sequences[name]) if letter != "-")
    color = page.evaluate("([letter, column]) => msaResidueColor(letter, column, 'protein')", [sequences[name][column], column])
    leaf = _leaf_point(page, name)
    cell = page.evaluate("""column => {
        const rect = document.getElementById('msa-canvas').getBoundingClientRect();
        return rect.left + (column + 0.5) * msaColumnWidth() - state.tree.msaScrollX;
    }""", column)
    assert _pixel(page, "msa-canvas", cell, leaf["y"]) == color.upper()
    assert page.pageerrors == []


def test_alignment_scrolls_and_zooms(msa_page):
    page = msa_page
    box = page.locator("#msa-canvas").bounding_box()
    page.mouse.move(box["x"] + 100, box["y"] + 200)
    page.keyboard.down("Shift")
    page.mouse.wheel(0, 300)
    page.keyboard.up("Shift")
    page.wait_for_timeout(50)
    assert page.evaluate("state.tree.msaScrollX") == 300
    page.keyboard.down("Control")
    page.mouse.wheel(0, -200)
    page.keyboard.up("Control")
    page.wait_for_timeout(50)
    assert float(page.input_value("#msa-column-width")) > 10


def test_fasta_export(msa_page, shared_datadir_module):
    page = msa_page
    sequences = {record.id: str(record.seq) for record in parse_seqfiles([shared_datadir_module / "FeSOD_20.afa"], filetype_override="fasta")}
    name = page.evaluate("state.bundle.nodes[0]")
    point = _leaf_point(page, name)
    page.mouse.click(point["x"] + 15, point["y"])
    with page.expect_download() as download:
        page.click("#export-fasta")
    assert download.value.suggested_filename == "fesod_selected.afa"
    exported = {record.id: str(record.seq) for record in parse_seqfiles([download.value.path()], filetype_override="fasta")}
    assert exported == {name: sequences[name]}
    page.select_option("#export-fasta-gaps", "ungapped")
    with page.expect_download() as download:
        page.click("#export-fasta")
    exported = {record.id: str(record.seq) for record in parse_seqfiles([download.value.path()], filetype_override="fasta")}
    assert exported == {name: sequences[name].replace("-", "")}


def test_ctrl_wheel_zooms_rows(msa_page):
    page = msa_page
    before = page.evaluate("state.tree.rowHeight")
    box = page.locator("#tree-canvas").bounding_box()
    page.mouse.move(box["x"] + 100, box["y"] + 100)
    page.keyboard.down("Control")
    page.mouse.wheel(0, -300)
    page.keyboard.up("Control")
    page.wait_for_timeout(50)
    assert page.evaluate("state.tree.rowHeight") > before


def test_one_page_opens_network_and_tree_bundles(tmp_path, shared_datadir_module):
    names = ["a", "b", "c"]
    matrix = DataMatrix(np.array([[0, 5, 1], [5, 0, 0], [1, 0, 0]], dtype=float), names, names, data_type="score")
    matrix_path = tmp_path / "m.hdf5"
    DataMatrix.write_hdf(matrix.data, matrix_path, names, names, data_type="score")
    network = tmp_path / "net.dvb"
    html = tmp_path / "net.html"
    build_ssn.main(["-i", str(matrix_path), "--dvb_out", str(network), "--html_out", str(html)])
    dvb = tmp_path / "small.dvb"
    build_tree.main(["--newick", str(shared_datadir_module / "tree_support.newick"), "--dvb_out", str(dvb)])
    for page in _loaded_page(html):
        assert page.evaluate("document.body.dataset.activeView") == "ssn"
        page.set_input_files("#bundle-file", str(dvb))
        page.wait_for_function("() => document.body.dataset.activeView === 'tree'")
        assert page.is_visible("#tree-canvas") and not page.is_visible("#cluster-view")
        assert page.evaluate("state.tree.layout.rowNode.length") == 7
        page.set_input_files("#bundle-file", str(network))
        page.wait_for_function("() => document.body.dataset.activeView === 'ssn'")
        page.wait_for_function("() => document.getElementById('stat-clusters').textContent !== '0'")
        assert page.is_visible("#cluster-view") and not page.is_visible("#tree-canvas")
        assert page.locator("#metadata-table tbody tr").count() == 3
        assert page.pageerrors == []
