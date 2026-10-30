"""Property-mapped node shapes and independent tree annotation slots."""

import base64
import gzip
import json
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import pytest

from domainator import build_ssn
from domainator.data_matrix import DataMatrix
from domainator.viewer import _read_bundle
from domainator.viewer_bundle import serialize_bundle, viewer_bundle, viewer_defaults
from domainator.viewer_html import write_viewer_html
from test_viewer_image_export_browser import _pixels

sync_api = pytest.importorskip("playwright.sync_api")
SVG = "{http://www.w3.org/2000/svg}"


@pytest.fixture(scope="module")
def appearance_pages(tmp_path_factory):
    directory = tmp_path_factory.mktemp("appearance")
    names = list("ABCDEFGH")
    metadata = pd.DataFrame({
        "group": ["a", "b", "a", "b", "c", "c", None, ""],
        "score": list(range(8)),
        "description": ["description_" + name * 12 for name in names],
    }, index=names)
    pages = {}
    for mode in ("tree", "projection"):
        if mode == "tree":
            view = {
                "parent": [-1, 0, 1, 1, 1, 1, 0, 0, 0, 0],
                "length": [0, 0.3, 0.7, 0.5, 0.4, 0.6, 1, 1, 1, 1],
                "label": [None, "clade_abcd"] + [None] * 8,
                "support": [None, 90] + [None] * 8,
                "leaf_nodes": [2, 3, 4, 5, 6, 7, 8, 9], "msa": None,
            }
        else:
            view = {"algorithm": "pca", "axes": ["X", "Y"],
                    "coords": [list(range(8)), [0, 3, 1, 4, 2, 6, 5, 7]]}
        bundle = viewer_bundle(mode, names, metadata, viewer_defaults(metadata, color_by="group"), {mode: view})
        pages[mode] = directory / f"{mode}.html"
        write_viewer_html(pages[mode], embedded_bundle=serialize_bundle(bundle))
    metadata_path = directory / "metadata.tsv"
    metadata.to_csv(metadata_path, sep="\t")
    scores = np.zeros((8, 8))
    for i in range(7):
        scores[i, i + 1] = scores[i + 1, i] = 10 - i
    matrix_path = directory / "scores.hdf5"
    DataMatrix(scores, names, names).write(matrix_path, output_type="hdf")
    pages["ssn"] = directory / "ssn.html"
    build_ssn.main(["-i", str(matrix_path), "--html_out", str(pages["ssn"]), "--metadata",
                    str(metadata_path), "--color_by", "group"])
    return pages


@pytest.fixture
def open_appearance(appearance_pages):
    with sync_api.sync_playwright() as playwright:
        try:
            browser = playwright.chromium.launch()
        except sync_api.Error as exc:
            pytest.skip(f"Chromium unavailable; run 'playwright install chromium': {exc}")

        def load(mode):
            page = browser.new_page(viewport={"width": 1500, "height": 1200}, accept_downloads=True)
            page.pageerrors = []
            page.requests = []
            page.on("pageerror", lambda exc: page.pageerrors.append(str(exc)))
            page.on("request", lambda request: page.requests.append(request.url)
                    if not request.url.startswith(("file:", "data:", "blob:")) else None)
            page.goto(appearance_pages[mode].as_uri())
            page.wait_for_function("() => state.bundle && (state.activeViewName !== 'ssn' || state.visibleLayout.length > 0)")
            return page

        yield load
        browser.close()


def _download(page, button):
    with page.expect_download() as downloaded:
        page.click(button)
    return downloaded.value.path()


def _node(page, index):
    return page.locator("#node-styles .appearance-entry").nth(index)


def _label(page, index):
    return page.locator("#label-styles .appearance-entry").nth(index)


def _three_nodes_two_labels(page):
    page.select_option("#shape-by", "score")
    page.click("#add-node-style")
    _node(page, 1).locator('[data-appearance-key="color_by"]').select_option("score")
    _node(page, 1).locator('[data-appearance-key="shape_by"]').select_option("group")
    page.click("#add-node-style")
    _node(page, 2).locator('[data-appearance-key="color_by"]').select_option("group")
    _node(page, 2).locator('[data-appearance-key="shape_by"]').select_option("score")
    page.click("#add-label-style")
    _label(page, 1).locator("select").select_option("description")


@pytest.mark.parametrize("mode", ["tree", "projection", "ssn"])
def test_shape_mapping_is_categorical_without_changing_gradient(open_appearance, mode):
    page = open_appearance(mode)
    assert page.locator("legend", has_text="Nodes").count() == 1
    assert page.locator("legend", has_text="Labels").count() == 1
    assert page.is_visible("#add-node-style") is (mode == "tree")
    page.select_option("#color-by", "score")
    colors = page.evaluate("state.nodeColorCache.slice()")
    page.select_option("#shape-by", "score")
    assert page.evaluate("state.nodeColorCache") == colors
    assert page.evaluate("colorInfo('score').type") == "numeric"
    assert page.evaluate("state.allNodeIndices.map(i => shapeForNode(state.nodeStyles[0], i))") == [
        "circle", "square", "diamond", "triangle_up", "triangle_down", "hexagon", "circle", "square",
    ]
    _node(page, 0).get_by_text("Customize shapes...", exact=True).click()
    assert "Shapes repeat" in page.text_content("#shape-picker-note")
    page.get_by_label("Shape for 1", exact=True).select_option("hexagon")
    assert page.evaluate("shapeForNode(state.nodeStyles[0], 1)") == "hexagon"
    page.click("#shape-picker-reset")
    assert page.evaluate("shapeForNode(state.nodeStyles[0], 1)") == "square"
    page.click("#shape-picker-close")
    assert page.pageerrors == [] and page.requests == []


def test_tree_multiple_slots_render_and_export_without_overlaps(open_appearance):
    page = open_appearance("tree")
    _three_nodes_two_labels(page)
    page.select_option("#export-extent", "entire")
    svg = ET.fromstring(_download(page, "#export-svg").read_bytes())
    texts = [element.text for element in svg.iter(SVG + "text")]
    for name in "ABCDEFGH":
        assert name in texts and "description_" + name * 12 in texts
    filled = [element for element in svg.iter(SVG + "path") if element.get("fill") not in (None, "none")]
    assert len(filled) == 24
    assert sum("L" in element.get("d") for element in filled) >= 12
    assert page.evaluate("""() => {
        const size = treeExportSize();
        const g = treeGeometry(size.treeWidth, size.height, size.options);
        const widths = g.labelWidths;
        return state.tree.layout.rowNode.every(node =>
            g.xPx(node) + g.markerWidth + widths.reduce((a, b) => a + b, 0) <= size.treeWidth);
    }""")
    assert _download(page, "#export-png").read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    page.locator("#tree-canvas").scroll_into_view_if_needed()
    point = page.evaluate("""() => {
        const rect = document.getElementById('tree-canvas').getBoundingClientRect();
        const g = treeGeometry(rect.width, rect.height);
        const node = state.tree.view.leaf_nodes[0];
        return {x: rect.x + g.xPx(node) + 2 * g.markerStep, y: rect.y + g.yPx(state.tree.layout.nodeRow[node])};
    }""")
    page.mouse.click(**point)
    assert page.evaluate("Array.from(state.selectedNodeIndices)") == [0]
    assert page.pageerrors == []


def test_tree_reordering_primary_changes_branches_and_empty_lists_round_trip(open_appearance):
    page = open_appearance("tree")
    _three_nodes_two_labels(page)
    ids = page.evaluate("state.nodeStyles.map(style => style.id)")
    _node(page, 1).get_by_role("button", name="Move up node 2", exact=True).click()
    assert page.evaluate("state.nodeStyles.map(style => style.id)") == [ids[1], ids[0], ids[2]]
    assert page.input_value("#color-by") == "score"
    assert page.evaluate("treeCladeSummaries().color[state.tree.view.leaf_nodes[0]]") == page.evaluate("nodeColor(0, 'score')")
    for _ in range(3):
        _node(page, 0).get_by_role("button", name="Remove node 1", exact=True).click()
    for _ in range(2):
        _label(page, 0).get_by_role("button", name="Remove label 1", exact=True).click()
    assert page.evaluate("treeCladeSummaries().coloring") is False
    svg = ET.fromstring(_download(page, "#export-svg").read_bytes())
    assert not [element for element in svg.iter(SVG + "path") if element.get("fill") not in (None, "none")]
    assert not set("ABCDEFGH") & {element.text for element in svg.iter(SVG + "text")}
    saved = _download(page, "#save-session")
    _, payload = _read_bundle(saved)
    assert payload["version"] == 2
    assert payload["app_state"]["view"]["node_styles"] == []
    assert payload["app_state"]["view"]["label_styles"] == []
    page.set_input_files("#bundle-file", saved)
    page.wait_for_function("() => state.nodeStyles.length === 0 && state.labelStyles.length === 0")
    assert "skipped" not in page.text_content("#bundle-status")
    assert page.pageerrors == []


def test_secondary_color_editor_and_metadata_lifecycle(open_appearance):
    page = open_appearance("tree")
    _three_nodes_two_labels(page)
    primary = page.evaluate("state.nodeColorCache.slice()")
    _node(page, 1).get_by_text("Customize colors...", exact=True).click()
    assert "score" in page.text_content("#color-picker-title")
    page.check("#color-as-categorical")
    assert page.evaluate("state.nodeColorCache") == primary
    page.click("#color-picker-close")
    _node(page, 1).get_by_text("Customize shapes...", exact=True).click()
    page.get_by_label("Shape for a", exact=True).select_option("hexagon")
    page.click("#shape-picker-close")
    page.evaluate("renameMetadataColumn('group', 'family')")
    assert page.evaluate("state.nodeStyles.map(s => s.color_by)") == ["family", "score", "family"]
    assert page.evaluate("shapeForNode(state.nodeStyles[1], 0)") == "hexagon"
    page.evaluate("setMetadataValue(0, 'family', 'new category')")
    assert page.evaluate("metadataValue(0, 'family')") == "new category"
    page.on("dialog", lambda dialog: dialog.accept())
    page.evaluate("deleteMetadataColumn('family')")
    assert page.evaluate("state.nodeStyles.map(s => s.color_by)") == ["", "score", ""]
    assert page.evaluate("state.nodeStyles[1].shape_by") == ""
    assert "Cleared" in page.text_content("#bundle-status")
    assert page.pageerrors == []


@pytest.mark.parametrize("mode", ["tree", "projection", "ssn"])
def test_new_session_round_trip_and_old_session_migration(open_appearance, mode, tmp_path):
    page = open_appearance(mode)
    page.select_option("#shape-by", "group")
    _node(page, 0).get_by_text("Customize shapes...", exact=True).click()
    page.get_by_label("Shape for a", exact=True).select_option("hexagon")
    page.get_by_label("Shape for No value", exact=True).select_option("diamond")
    page.click("#shape-picker-close")
    if mode == "tree":
        page.click("#add-label-style")
    saved = _download(page, "#save-session")
    _, payload = _read_bundle(saved)
    assert payload["app_state"]["shapes"]["custom_mappings"]["group"]["values"]["a"] == "hexagon"
    view = payload["app_state"]["view"]
    view["color_by"] = "description"
    # New arrays must win even when legacy aliases occur later in the file.
    payload["app_state"]["view"] = dict(reversed(list(view.items())))
    reordered = tmp_path / "reordered.dvb"
    reordered.write_bytes(serialize_bundle(payload))
    page.set_input_files("#bundle-file", reordered)
    page.wait_for_function("() => currentColorField() === 'group' && shapeForNode(state.nodeStyles[0], 0) === 'hexagon'")
    assert page.evaluate("shapeForNode(state.nodeStyles[0], 6)") == "diamond"
    assert "skipped" not in page.text_content("#bundle-status")
    payload["version"] = 1
    payload["app_state"] = {"view": {"color_by": "score", "label_by": "", "tree_show_labels": False}} if mode == "tree" else {
        "view": {"color_by": "score", "label_by": "description"},
    }
    legacy = tmp_path / "legacy.dvb"
    legacy.write_bytes(gzip.compress(json.dumps(payload).encode()))
    assert _read_bundle(legacy)[1]["version"] == 1
    page.set_input_files("#bundle-file", legacy)
    page.wait_for_function("() => currentColorField() === 'score'")
    assert page.evaluate("state.nodeStyles.length") == 1
    assert page.evaluate("state.nodeStyles[0].shape_by") == ""
    assert page.evaluate("Object.keys(state.customShapes)") == []
    if mode == "tree":
        assert not page.is_checked("#tree-show-labels")
        assert page.input_value("#label-by") == "__node_id__"
    assert page.pageerrors == []


@pytest.mark.parametrize("layout", ["packed", "tree", "force", "grid", "treemap"])
def test_network_shape_svg_and_png_in_every_layout(open_appearance, layout):
    page = open_appearance("ssn")
    page.select_option("#layout-algorithm", layout)
    page.wait_for_function("() => !state.layoutComputing && state.visibleLayout.length > 0")
    page.select_option("#shape-by", "score")
    page.select_option("#export-extent", "entire")
    svg = ET.fromstring(_download(page, "#export-svg").read_bytes())
    assert len(list(svg.iter(SVG + "polygon"))) == 6
    assert _download(page, "#export-png").read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert page.pageerrors == []


@pytest.mark.parametrize("mode", ["tree", "projection", "ssn"])
@pytest.mark.parametrize("extent", ["viewport", "entire"])
def test_mapped_png_and_svg_match_actual_canvas(open_appearance, mode, extent):
    page = open_appearance(mode)
    if mode == "tree":
        _three_nodes_two_labels(page)
    else:
        page.select_option("#shape-by", "score")
    page.select_option("#export-extent", extent)
    page.select_option("#export-png-scale", "1")
    samples = page.evaluate("""() => {
        const canvas = viewerCanvas();
        const rect = canvas.getBoundingClientRect();
        const screen = [], exported = [];
        if (state.activeViewName === 'tree') {
            const size = treeExportSize();
            const a = treeGeometry(rect.width, rect.height);
            const b = treeGeometry(size.treeWidth, size.height, size.options);
            for (let i = 0; i < 6; i++) {
                const node = state.tree.view.leaf_nodes[i];
                for (let j = 0; j < state.nodeStyles.length; j++) {
                    screen.push([a.xPx(node) + j * a.markerStep, a.yPx(state.tree.layout.nodeRow[node])]);
                    exported.push([b.xPx(node) + j * b.markerStep, b.yPx(state.tree.layout.nodeRow[node])]);
                }
            }
        } else if (state.activeViewName === 'projection') {
            const size = projectionExportSize();
            const a = projectionScales(rect.width, rect.height);
            const b = projectionScales(size.width, size.height, size.options.view, size.options.margin);
            const [xs, ys] = projectionPair().map(axis => state.projection.view.coords[axis]);
            for (let i = 0; i < 6; i++) {
                screen.push([a.sx(xs[i]), a.sy(ys[i])]);
                exported.push([b.sx(xs[i]), b.sy(ys[i])]);
            }
        } else {
            const size = clusterExportSize();
            for (const item of state.visibleLayout) {
                const component = state.bundle.views.ssn.hierarchy.nodes[item.componentId];
                for (const dot of componentMemberLayout(component, item)) {
                    const a = worldToScreenPoint(dot.x, dot.y);
                    const b = worldToScreenPoint(dot.x, dot.y, size.transform);
                    screen.push([a.x, a.y]); exported.push([b.x, b.y]);
                }
            }
        }
        return {screen, exported, image: canvas.toDataURL().split(',')[1]};
    }""")
    before = page.evaluate("""() => [state.viewTransform, state.tree && [state.tree.rowHeight, state.tree.scrollY],
        state.projection && [state.projection.center, state.projection.unitsPerPixel]]""")
    expected = _pixels(page, base64.b64decode(samples["image"]), "image/png", samples["screen"])
    assert all(min(color) < 225 for color in expected)
    for button, mime in (("#export-png", "image/png"), ("#export-svg", "image/svg+xml")):
        data = _download(page, button).read_bytes()
        observed = _pixels(page, data, mime, samples["exported"])
        for actual, wanted in zip(observed, expected):
            assert actual == pytest.approx(wanted, abs=30)
    assert page.evaluate("""() => [state.viewTransform, state.tree && [state.tree.rowHeight, state.tree.scrollY],
        state.projection && [state.projection.center, state.projection.unitsPerPixel]]""") == before
    assert page.pageerrors == []


def test_dense_and_collapsed_tree_annotations(open_appearance):
    page = open_appearance("tree")
    _three_nodes_two_labels(page)
    dense_svg = page.evaluate("""() => {
        const old = state.tree.rowHeight;
        state.tree.rowHeight = 1;
        const ctx = createSvgContext(900, 400);
        drawTree(ctx, treeGeometry(900, 400), treeCladeSummaries());
        state.tree.rowHeight = old;
        return ctx.toSVG();
    }""")
    dense = ET.fromstring(dense_svg)
    assert len([element for element in dense.iter(SVG + "rect") if element.get("width") == "8"]) == 24
    page.evaluate("treeToggleCollapsed(1)")
    svg = ET.fromstring(_download(page, "#export-svg").read_bytes())
    filled = [element for element in svg.iter(SVG + "path") if element.get("fill") not in (None, "none")]
    assert len(filled) == 13  # Four visible leaves times three markers, plus the wedge.
    texts = [element.text for element in svg.iter(SVG + "text")]
    assert texts.count("clade_abcd (4 leaves)") == 1
    assert "description_" + "A" * 12 not in texts
    assert page.pageerrors == []


def test_invalid_appearance_fields_are_atomic_and_reported(open_appearance, tmp_path):
    page = open_appearance("tree")
    saved = _download(page, "#save-session")
    _, payload = _read_bundle(saved)
    view = payload["app_state"]["view"]
    view["color_by"] = "score"
    view["node_styles"] = [
        {"id": "node-1", "color_by": "group", "shape_by": ""},
        {"id": "node-2", "color_by": "missing", "shape_by": ""},
    ]
    payload["app_state"]["shapes"]["custom_mappings"] = {"group": {"values": {"a": "not-a-shape"}, "null_shape": "circle"}}
    path = tmp_path / "invalid.dvb"
    path.write_bytes(serialize_bundle(payload))
    page.set_input_files("#bundle-file", path)
    page.wait_for_function("() => document.getElementById('bundle-status').textContent.includes('skipped 2')")
    assert page.input_value("#color-by") == "score"
    assert page.evaluate("state.nodeStyles.length") == 1
    assert "view.node_styles" in page.text_content("#bundle-status")
    assert "shapes.custom_mappings" in page.text_content("#bundle-status")
    assert page.pageerrors == []


def test_network_extraction_retains_shape_assignments(open_appearance, tmp_path):
    page = open_appearance("ssn")
    page.select_option("#shape-by", "score")
    _node(page, 0).get_by_text("Customize shapes...", exact=True).click()
    page.get_by_label("Shape for 0", exact=True).select_option("hexagon")
    page.click("#shape-picker-close")
    page.evaluate("state.selectedNodeIndices = new Set(state.allNodeIndices)")
    result = page.evaluate("buildExtractionBundle()")
    assert "error" not in result
    bundle = result["bundle"]
    assert bundle["app_state"]["view"]["node_styles"][0]["shape_by"] == "score"
    assert bundle["app_state"]["shapes"]["custom_mappings"]["score"]["values"]["0"] == "hexagon"
    path = tmp_path / "extracted.dvb"
    path.write_bytes(serialize_bundle(bundle))
    page.set_input_files("#bundle-file", path)
    page.wait_for_function("() => shapeForNode(state.nodeStyles[0], 0) === 'hexagon'")
    assert "skipped" not in page.text_content("#bundle-status")
    assert page.pageerrors == []


def test_loading_different_view_resets_annotation_and_editor_state(open_appearance, tmp_path):
    page = open_appearance("tree")
    _three_nodes_two_labels(page)
    saved = _download(page, "#save-session")
    _, payload = _read_bundle(saved)
    _node(page, 1).get_by_text("Customize shapes...", exact=True).click()
    page.get_by_label("Shape for a", exact=True).select_option("hexagon")
    payload.pop("app_state")
    payload["views"] = {"projection": {"algorithm": "pca", "axes": ["X", "Y"],
                                      "coords": [list(range(8)), list(range(8))]}}
    path = tmp_path / "other_view.dvb"
    path.write_bytes(serialize_bundle(payload))
    page.set_input_files("#bundle-file", path)
    page.wait_for_function("() => state.activeViewName === 'projection'")
    assert page.evaluate("state.nodeStyles.length") == 1
    assert page.evaluate("state.labelStyles.length") == 1
    assert page.evaluate("state.nodeStyles[0].shape_by") == ""
    assert page.evaluate("Object.keys(state.customShapes)") == []
    assert not page.is_visible("#shape-picker-overlay")
    assert not page.is_visible("#add-node-style")
    assert page.pageerrors == []
