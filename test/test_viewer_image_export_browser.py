"""Verify complete-scene PNG/SVG exports, not just successful downloads."""

import base64
import gzip
import json
import struct
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import pytest

from domainator import build_ssn
from domainator.data_matrix import DataMatrix
from domainator.viewer_bundle import serialize_bundle, viewer_bundle, viewer_defaults
from domainator.viewer_html import write_viewer_html

sync_api = pytest.importorskip("playwright.sync_api")


@pytest.fixture(scope="module")
def image_pages(tmp_path_factory):
    directory = tmp_path_factory.mktemp("entire_content")
    pages = {}
    for mode, count in (("tree", 60), ("projection", 6)):
        names = [f"leaf_{i:03d}" for i in range(count)]
        names[-1] += "_" + "long_label_" * 8
        if mode == "tree":
            sequence = "ACDEFGHIKLMNPQRSTVWY" * 20
            view = {
                "parent": [-1] + [0] * count, "length": [0] + [1] * count,
                "label": [None] * (count + 1), "support": [None] * (count + 1),
                "leaf_nodes": list(range(1, count + 1)),
                "msa": {"alphabet": "protein", "length": len(sequence), "sequences": [sequence] * count},
            }
        else:
            view = {"algorithm": "pca", "axes": ["X", "Y", "Z"],
                    "coords": [[0, 1, 2, 3, 4, 5], [2, 0, 4, 1, 3, 5], [0, 2, 1, 3, 5, 4]]}
        metadata = pd.DataFrame(index=names)
        bundle = viewer_bundle(mode, names, metadata, viewer_defaults(metadata), {mode: view})
        path = directory / f"{mode}.html"
        write_viewer_html(path, embedded_bundle=serialize_bundle(bundle))
        pages[mode] = path
    names = list("ABCDEF")
    matrix = np.zeros((6, 6))
    for i in range(5):
        matrix[i, i + 1] = matrix[i + 1, i] = 10 - i
    matrix_path = directory / "matrix.hdf5"
    DataMatrix(matrix, names, names).write(matrix_path, output_type="hdf")
    pages["ssn"] = directory / "network.html"
    build_ssn.main(["-i", str(matrix_path), "--html_out", str(pages["ssn"]), "--name", "network"])
    return pages


@pytest.fixture
def open_image_page(image_pages):
    with sync_api.sync_playwright() as playwright:
        try:
            browser = playwright.chromium.launch()
        except sync_api.Error as exc:
            pytest.skip(f"Chromium unavailable; run 'playwright install chromium': {exc}")

        def load(mode):
            page = browser.new_page(viewport={"width": 1400, "height": 1100}, accept_downloads=True)
            page.pageerrors = []
            page.requests = []
            page.on("pageerror", lambda exc: page.pageerrors.append(str(exc)))
            page.on("request", lambda request: page.requests.append(request.url)
                    if not request.url.startswith(("file:", "data:", "blob:")) else None)
            page.goto(image_pages[mode].as_uri())
            page.wait_for_function("() => document.getElementById('bundle-status').textContent.startsWith('Loaded')")
            page.wait_for_function("() => state.activeViewName !== 'ssn' || state.visibleLayout.length > 0")
            return page

        try:
            yield load
        finally:
            browser.close()


def _download(page, button):
    with page.expect_download() as downloaded:
        page.click(button)
    return downloaded.value.path().read_bytes(), downloaded.value.suggested_filename


def _pixels(page, data, mime, points):
    return page.evaluate("""async ({data, mime, points}) => {
        const image = new Image();
        await new Promise((resolve, reject) => {
            image.onload = resolve; image.onerror = reject;
            image.src = 'data:' + mime + ';base64,' + data;
        });
        const canvas = document.createElement('canvas');
        canvas.width = image.width; canvas.height = image.height;
        const ctx = canvas.getContext('2d');
        ctx.drawImage(image, 0, 0);
        return points.map(([x, y]) => Array.from(ctx.getImageData(Math.round(x), Math.round(y), 1, 1).data).slice(0, 3));
    }""", {"data": base64.b64encode(data).decode("ascii"), "mime": mime, "points": points})


_SNAPSHOT = """() => ({
    transform: state.viewTransform, tree: state.tree && [state.tree.rowHeight, state.tree.scrollY,
        state.tree.msaScrollX, Array.from(state.tree.collapsed)],
    projection: state.projection && [state.projection.center, state.projection.unitsPerPixel],
    selected: Array.from(state.selectedNodeIndices),
    layout: state.visibleLayout,
    canvases: ['tree-canvas', 'msa-canvas', 'cluster-view', 'projection-canvas'].map(id => {
        const c = document.getElementById(id); return [c.width, c.height, c.toDataURL()];
    })
})"""


@pytest.mark.parametrize("mode,layout", [
    ("tree", None), ("projection", "0,1"), ("projection", "0,2"), ("projection", "1,2"),
    ("ssn", "packed"), ("ssn", "tree"), ("ssn", "force"), ("ssn", "grid"), ("ssn", "treemap"),
])
def test_entire_content_exports_offscreen_geometry_without_changing_view(open_image_page, mode, layout):
    page = open_image_page(mode)
    if mode == "projection":
        page.select_option("#projection-axes", layout)
    elif mode == "ssn":
        page.evaluate("""() => {
            const slider = document.getElementById('threshold-slider');
            slider.value = slider.max; slider.dispatchEvent(new Event('input'));
        }""")
        page.select_option("#layout-algorithm", layout)
        page.wait_for_function("() => !state.layoutComputing")
    page.select_option("#label-by", "__node_id__")
    page.evaluate("""() => {
        state.selectedNodeIndices = new Set([0]);
        updateMetadataTable();
        if (state.activeViewName === 'tree') {
            state.tree.rowHeight = 16; state.tree.scrollY = 100000; state.tree.msaScrollX = 100000;
        } else if (state.activeViewName === 'projection') {
            state.projection.unitsPerPixel /= 2;
            state.projection.center = {x: 1000, y: 1000};
        } else {
            state.viewTransform.scale *= 2;
            state.viewTransform.offsetX = state.viewTransform.offsetY = -100000;
        }
        repaintNow();
    }""")
    before = page.evaluate(_SNAPSHOT)
    viewport_svg, viewport_name = _download(page, "#export-svg")
    first = page.evaluate("state.bundle.nodes[0]")
    assert f">{first}</text>".encode() not in viewport_svg
    assert viewport_name.endswith("_view.svg")
    page.select_option("#export-extent", "entire")
    description = page.evaluate("""() => {
        const mode = state.activeViewName;
        const size = mode === 'tree' ? treeExportSize() : mode === 'ssn' ? clusterExportSize() : projectionExportSize();
        let points;
        if (mode === 'tree') {
            const g = treeGeometry(size.treeWidth, size.height, size.options);
            const rows = [0, state.tree.layout.rowNode.length - 1];
            points = rows.map(row => [g.xPx(state.tree.layout.rowNode[row]), g.yPx(row)]);
            for (const row of rows) {
                for (const col of [0, state.tree.msa.length - 1]) {
                    points.push([size.treeWidth + (col + 0.2) * msaColumnWidth(), g.yPx(row) - 4]);
                }
            }
        } else if (mode === 'projection') {
            const scales = projectionScales(size.width, size.height, size.options.view, size.options.margin);
            points = [];
            forEachProjectionPoint(scales, projectionPointRadius(), (_, x, y) => points.push([x, y]));
        } else {
            points = state.visibleLayout.flatMap(item =>
                Array.from(componentMemberLayout(state.bundle.views.ssn.hierarchy.nodes[item.componentId], item),
                    dot => { const p = worldToScreenPoint(dot.x + dot.radius * 0.45, dot.y + dot.radius * 0.45, size.transform);
                        return [p.x, p.y]; }));
        }
        return {width: Math.round(size.width), height: Math.round(size.height), points, nodes: state.bundle.nodes};
    }""")
    svg, name = _download(page, "#export-svg")
    assert name.endswith("_entire.svg")
    root = ET.fromstring(svg)
    assert int(root.attrib["width"]) == description["width"]
    assert int(root.attrib["height"]) == description["height"]
    assert root.attrib["viewBox"] == f"0 0 {description['width']} {description['height']}"
    texts = {element.text for element in root.iter("{http://www.w3.org/2000/svg}text")}
    if mode != "ssn":
        assert set(description["nodes"]) <= texts
    clipped_text = page.evaluate("""svg => {
        const host = document.createElement('div');
        host.style.cssText = 'position:fixed;left:-100000px;top:0';
        host.innerHTML = svg; document.body.appendChild(host);
        try {
            const image = host.querySelector('svg');
            const width = Number(image.getAttribute('width')), height = Number(image.getAttribute('height'));
            return Array.from(image.querySelectorAll('text')).filter(text => {
                const b = text.getBBox();
                return b.x < -2 || b.y < -2 || b.x + b.width > width + 2 || b.y + b.height > height + 2;
            }).map(text => text.textContent);
        } finally { host.remove(); }
    }""", svg.decode())
    assert clipped_text == []
    for x, y in description["points"]:
        assert 0 < x < description["width"] and 0 < y < description["height"]
    svg_pixels = _pixels(page, svg, "image/svg+xml", description["points"])
    assert all(min(pixel) < 245 for pixel in svg_pixels), svg_pixels
    for scale in (1, 2):
        page.select_option("#export-png-scale", str(scale))
        png, name = _download(page, "#export-png")
        assert "_entire" in name
        assert struct.unpack(">II", png[16:24]) == (description["width"] * scale, description["height"] * scale)
        points = [[x * scale, y * scale] for x, y in description["points"]]
        png_pixels = _pixels(page, png, "image/png", points)
        assert all(min(pixel) < 245 for pixel in png_pixels), png_pixels
        assert all(max(abs(a - b) for a, b in zip(left, right)) < 45
                   for left, right in zip(svg_pixels, png_pixels)), (scale, svg_pixels, png_pixels)
    assert page.evaluate(_SNAPSHOT) == before
    assert page.pageerrors == [] and page.requests == []


@pytest.mark.parametrize("mode", ["tree", "ssn", "projection"])
def test_export_extent_session_and_old_session_default(open_image_page, mode):
    page = open_image_page(mode)
    page.select_option("#export-extent", "entire")
    data, _ = _download(page, "#save-session")
    bundle = json.loads(gzip.decompress(data))
    assert bundle["app_state"]["view"]["export_extent"] == "entire"
    page.evaluate("bundle => installBundle(bundle)", bundle)
    assert page.input_value("#export-extent") == "entire"
    del bundle["app_state"]["view"]["export_extent"]
    page.evaluate("bundle => installBundle(bundle)", bundle)
    assert page.input_value("#export-extent") == "viewport"
    assert page.pageerrors == []


@pytest.mark.parametrize("failure", ["size", "context", "encoding"])
def test_png_failure_is_explicit_and_never_downloads(open_image_page, failure):
    page = open_image_page("projection")
    downloads = []
    page.on("download", lambda download: downloads.append(download))
    page.evaluate("""failure => {
        if (failure === 'size') {
            projectionExportSize = () => ({width: 40000, height: 40000, options: {}});
        } else if (failure === 'context') {
            HTMLCanvasElement.prototype.getContext = () => null;
        } else {
            HTMLCanvasElement.prototype.toBlob = callback => callback(null);
        }
    }""", failure)
    page.click("#export-png")
    page.wait_for_function("() => document.getElementById('bundle-status').textContent.startsWith('Image export failed:')")
    assert "SVG" in page.text_content("#bundle-status")
    assert not downloads
    assert page.pageerrors == []


def test_png_allocation_limits_are_checked_exactly(open_image_page):
    page = open_image_page("projection")
    result = page.evaluate("""() => [[8000, 8000], [8001, 8000], [32767, 1], [32768, 1], [0, 1], [1.5, 1]].map(([w, h]) => {
        try { validateImageExportSize(w, h, true); return true; } catch (_) { return false; }
    })""")
    assert result == [True, False, True, False, False, False]


@pytest.mark.parametrize("count", [0, 1, 3])
def test_entire_projection_handles_empty_singleton_and_coincident_points(open_image_page, count):
    page = open_image_page("projection")
    page.evaluate("""count => {
        const bundle = JSON.parse(JSON.stringify(state.bundle));
        bundle.nodes = Array.from({length: count}, (_, i) => 'p' + i);
        bundle.metadata.rows = Array.from({length: count}, () => []);
        bundle.views.projection.coords = bundle.views.projection.axes.map(() => Array(count).fill(5));
        installBundle(bundle);
    }""", count)
    page.select_option("#export-extent", "entire")
    svg, _ = _download(page, "#export-svg")
    root = ET.fromstring(svg)
    assert int(root.attrib["width"]) > 0 and int(root.attrib["height"]) > 0
    png, _ = _download(page, "#export-png")
    assert png.startswith(b"\x89PNG\r\n\x1a\n")
    assert page.pageerrors == []


def test_entire_network_does_not_reveal_hidden_nodes(open_image_page):
    page = open_image_page("ssn")
    page.select_option("#label-by", "__node_id__")
    page.fill("#min-cluster-size", "100")
    page.locator("#min-cluster-size").dispatch_event("change")
    page.wait_for_function("() => state.visibleLayout.length === 0")
    page.select_option("#export-extent", "entire")
    svg, _ = _download(page, "#export-svg")
    root = ET.fromstring(svg)
    assert not list(root.iter("{http://www.w3.org/2000/svg}circle"))
    assert not list(root.iter("{http://www.w3.org/2000/svg}text"))
    assert page.pageerrors == []


@pytest.mark.parametrize("extent", ["viewport", "entire"])
@pytest.mark.parametrize("position", [0, 1, 250, 999, 1000])
def test_cut_line_hidden_only_at_endpoints_in_canvas_png_and_svg(open_image_page, extent, position):
    page = open_image_page("tree")
    assert page.input_value("#tree-cut") == "0"
    assert page.locator("#tree-show-cut-line").count() == 0
    page.evaluate("""position => {
        const cut = document.getElementById('tree-cut');
        cut.value = position; cut.dispatchEvent(new Event('input'));
    }""", position)
    visible = 0 < position < 1000
    assignments = page.evaluate("state.bundle.nodes.map((_, i) => treeClusterNumbers().numberFor(i))")
    assert len(set(assignments)) == (1 if position == 0 else 60)
    page.select_option("#export-extent", extent)
    page.select_option("#export-png-scale", "1")
    point = page.evaluate("""() => {
        const size = treeExportSize();
        const g = treeGeometry(size.treeWidth, size.height, size.options);
        const x = Number.isFinite(treeCutX()) ? treeCutX() : state.tree.layout.maxX / 4;
        return [g.xOf(x), g.top + 2];
    }""")
    live_pixel = """() => {
        const canvas = document.getElementById('tree-canvas');
        const r = canvas.getBoundingClientRect(), g = treeGeometry(r.width, r.height);
        const ratio = canvas.width / r.width;
        const x = Number.isFinite(treeCutX()) ? treeCutX() : state.tree.layout.maxX / 4;
        return Array.from(canvas.getContext('2d').getImageData(
            Math.round(g.xOf(x) * ratio), Math.round((g.top + 2) * ratio), 1, 1).data).slice(0, 3);
    }"""
    assert (min(page.evaluate(live_pixel)) < 245) is visible
    svg, _ = _download(page, "#export-svg")
    assert (b'stroke-dasharray="5 4"' in svg) is visible
    assert (min(_pixels(page, svg, "image/svg+xml", [point])[0]) < 245) is visible
    png, _ = _download(page, "#export-png")
    assert (min(_pixels(page, png, "image/png", [point])[0]) < 245) is visible
    assert page.evaluate("state.bundle.nodes.map((_, i) => treeClusterNumbers().numberFor(i))") == assignments
    assert page.pageerrors == []
