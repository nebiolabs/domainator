"""Browser-driven tests of the Domainator Viewer's projection view (build_projection --html_out).

Like test_ssn_viewer_browser.py, these load the self-contained page in headless Chromium via
Playwright and skip when it is not installed (pip install -e ".[test,browser]" and playwright
install chromium).
"""

import gzip
import json

import pandas as pd
import pytest

from domainator import build_projection
from domainator.colors import get_palette

sync_api = pytest.importorskip("playwright.sync_api")


@pytest.fixture
def open_page():
    with sync_api.sync_playwright() as p:
        try:
            browser = p.chromium.launch()
        except sync_api.Error as exc:
            pytest.skip(f"Chromium not available for Playwright; run 'playwright install chromium'. ({exc})")

        def load(path):
            page = browser.new_page(viewport={"width": 1400, "height": 1100}, accept_downloads=True)
            page.pageerrors = []
            page.requests = []
            page.on("pageerror", lambda exc: page.pageerrors.append(str(exc)))
            page.on("request", lambda request: page.requests.append(request.url) if not request.url.startswith(("file:", "data:", "blob:")) else None)
            page.goto(path.as_uri())
            page.wait_for_function("() => document.getElementById('bundle-status').textContent.startsWith('Loaded')")
            return page

        try:
            yield load
        finally:
            browser.close()


def _projection_page(path, *extra):
    build_projection.main(["--html_out", str(path), *extra])
    return path


def _colored_page(tmp_path, shared_datadir, *extra):
    return _projection_page(tmp_path / "projection.html", "-i", str(shared_datadir / "scorefull.tsv"),
                            "--metadata", str(shared_datadir / "metadata_FeSOD_20.tsv"), "--color_by", "group", *extra)


def _point(page, name):
    """Where point `name` is on the screen, in client coordinates."""
    return page.evaluate("""name => {
        const rect = document.getElementById('projection-canvas').getBoundingClientRect();
        const scales = projectionScales(rect.width, rect.height);
        const i = state.bundle.nodes.indexOf(name);
        const [xs, ys] = projectionPair().map(axis => state.projection.view.coords[axis]);
        return {x: rect.left + scales.sx(xs[i]), y: rect.top + scales.sy(ys[i])};
    }""", name)


def _pixel(page, x, y):
    return page.evaluate("""([x, y]) => {
        const canvas = document.getElementById('projection-canvas');
        const rect = canvas.getBoundingClientRect();
        const ratio = canvas.width / rect.width;
        const [r, g, b] = canvas.getContext('2d').getImageData(Math.round((x - rect.left) * ratio), Math.round((y - rect.top) * ratio), 1, 1).data;
        return [r, g, b];
    }""", [x, y])


def _domain(page):
    return page.evaluate("""() => {
        const rect = document.getElementById('projection-canvas').getBoundingClientRect();
        const s = projectionScales(rect.width, rect.height);
        return {x0: s.x0, x1: s.x1, y0: s.y0, y1: s.y1};
    }""")


def _selected(page):
    return sorted(page.evaluate("Array.from(state.selectedNodeIndices).map(i => state.bundle.nodes[i])"))


def test_projection_draws_every_point_in_its_color(open_page, tmp_path, shared_datadir):
    page = open_page(_colored_page(tmp_path, shared_datadir))
    assert page.evaluate("document.body.dataset.activeView") == "projection"
    assert not page.is_visible("#tree-canvas") and not page.is_visible("#cluster-view")
    # Projections have no clusters to fill a column with.
    assert page.evaluate("getComputedStyle(document.getElementById('metadata-add-cluster-column')).display") == "none"
    metadata = pd.read_csv(shared_datadir / "metadata_FeSOD_20.tsv", sep="\t", index_col=0)["group"]
    palette = get_palette(metadata)
    for name in page.evaluate("state.bundle.nodes"):
        point = _point(page, name)
        red, green, blue = _pixel(page, point["x"], point["y"])
        expected = [int(palette[metadata[name]][i:i + 2], 16) for i in (1, 3, 5)]
        # Points are drawn at 75% opacity over white, and may overlap.
        assert max(abs(channel - (0.75 * want + 0.25 * 255)) for channel, want in zip((red, green, blue), expected)) < 40, name
    assert page.pageerrors == [] and page.requests == []


def test_projection_downloads(open_page, tmp_path, shared_datadir):
    page = open_page(_colored_page(tmp_path, shared_datadir))
    with page.expect_download() as download:
        page.click("#export-svg")
    assert download.value.suggested_filename == "projection_view.svg"
    svg = open(download.value.path()).read()
    assert svg.startswith("<svg") and ">projection_X</text>" in svg and ">projection_Y</text>" in svg
    with page.expect_download() as download:
        page.click("#export-png")
    assert open(download.value.path(), "rb").read(8) == b"\x89PNG\r\n\x1a\n"


def test_projection_hover_zoom_and_pan(open_page, tmp_path, shared_datadir):
    page = open_page(_colored_page(tmp_path, shared_datadir))
    name = page.evaluate("state.bundle.nodes[3]")
    point = _point(page, name)
    page.mouse.move(point["x"], point["y"])
    assert page.is_visible("#view-tip") and page.text_content("#view-tip").startswith(name)
    fitted = _domain(page)
    # Equal units on both axes.
    rect = page.locator("#projection-canvas").bounding_box()
    assert (fitted["x1"] - fitted["x0"]) / (fitted["y1"] - fitted["y0"]) == pytest.approx(
        page.evaluate("(() => { const r = document.getElementById('projection-canvas').getBoundingClientRect(); const a = projectionPlotArea(r.width, r.height); return a.width / a.height; })()"))
    page.mouse.wheel(0, -400)
    page.wait_for_timeout(50)
    zoomed = _domain(page)
    assert zoomed["x1"] - zoomed["x0"] < fitted["x1"] - fitted["x0"]
    # The point under the pointer stays there.
    after = _point(page, name)
    assert after["x"] == pytest.approx(point["x"], abs=1) and after["y"] == pytest.approx(point["y"], abs=1)
    page.click("#reset-view")
    assert _domain(page) == fitted
    page.mouse.move(rect["x"] + 300, rect["y"] + 300)
    page.mouse.down()
    page.mouse.move(rect["x"] + 360, rect["y"] + 300, steps=4)
    page.mouse.up()
    assert _domain(page)["x0"] < fitted["x0"]
    page.mouse.dblclick(rect["x"] + 300, rect["y"] + 300)
    assert _domain(page)["x1"] - _domain(page)["x0"] == pytest.approx((fitted["x1"] - fitted["x0"]) / 2)


def test_projection_selection_drives_the_table(open_page, tmp_path, shared_datadir):
    page = open_page(_colored_page(tmp_path, shared_datadir))
    names = page.evaluate("state.bundle.nodes")
    point = _point(page, names[0])
    page.mouse.click(point["x"], point["y"])
    assert _selected(page) == [names[0]]
    assert page.locator("#metadata-table tbody tr").count() == 1
    page.mouse.click(point["x"], point["y"])
    assert _selected(page) == []
    # A box around every point selects them all; Alt takes them away again.
    rect = page.locator("#projection-canvas").bounding_box()
    for keys, expected in ((["Shift"], 20), (["Shift", "Alt"], 0)):
        for key in keys:
            page.keyboard.down(key)
        page.mouse.move(rect["x"] + 66, rect["y"] + 18)
        page.mouse.down()
        page.mouse.move(rect["x"] + rect["width"] - 18, rect["y"] + rect["height"] - 48, steps=5)
        page.mouse.up()
        for key in reversed(keys):
            page.keyboard.up(key)
        assert len(_selected(page)) == expected


def test_projection_labels_and_focus(open_page, tmp_path, shared_datadir):
    page = open_page(_colored_page(tmp_path, shared_datadir, "--label_by", "group"))
    with page.expect_download() as download:
        page.click("#export-svg")
    assert ">g1</text>" in open(download.value.path()).read()
    name = page.evaluate("state.bundle.nodes[5]")
    point = _point(page, name)
    page.mouse.click(point["x"], point["y"])
    page.mouse.wheel(0, -1200)
    page.wait_for_timeout(50)
    rect = page.locator("#projection-canvas").bounding_box()
    page.mouse.move(rect["x"] + 200, rect["y"] + 200)
    page.mouse.down()
    page.mouse.move(rect["x"] + 600, rect["y"] + 500, steps=4)
    page.mouse.up()
    page.click("#focus-selection")
    centered = _point(page, name)
    area = page.evaluate("(() => { const r = document.getElementById('projection-canvas').getBoundingClientRect(); return projectionPlotArea(r.width, r.height); })()")
    assert centered["x"] == pytest.approx(rect["x"] + area["left"] + area["width"] / 2, abs=1)


def test_projection_3d_switches_axis_pairs(open_page, tmp_path, shared_datadir):
    page = open_page(_projection_page(tmp_path / "projection3d.html", "-i", str(shared_datadir / "scorefull.tsv"), "--three_d"))
    assert page.is_visible("#projection-axes")
    page.select_option("#projection-axes", "1,2")
    assert page.evaluate("projectionPair()") == [1, 2]
    with page.expect_download() as download:
        page.click("#export-svg")
    svg = open(download.value.path()).read()
    assert ">projection_Z</text>" in svg and ">projection_Y</text>" in svg and ">projection_X</text>" not in svg


def test_projection_session_round_trip(open_page, tmp_path, shared_datadir):
    path = _colored_page(tmp_path, shared_datadir)
    page = open_page(path)
    names = page.evaluate("state.bundle.nodes")
    point = _point(page, names[2])
    page.mouse.click(point["x"], point["y"])
    page.mouse.wheel(0, -300)
    page.fill("#projection-point-size", "7")
    page.dispatch_event("#projection-point-size", "input")
    view = page.evaluate("({center: state.projection.center, u: state.projection.unitsPerPixel})")
    with page.expect_download() as download:
        page.click("#save-session")
    assert download.value.suggested_filename == "projection_session.dvb"
    session = download.value.path()
    assert json.loads(gzip.open(session).read())["app_state"]["view"]["projection_point_size"] == "7"
    page = open_page(path)
    page.set_input_files("#bundle-file", session)
    page.wait_for_function("() => state.selectedNodeIndices.size === 1")
    assert _selected(page) == [names[2]]
    assert page.input_value("#projection-point-size") == "7"
    assert page.evaluate("({center: state.projection.center, u: state.projection.unitsPerPixel})") == view
    assert page.pageerrors == []
