"""Browser-driven regression tests for the generated SSN viewer HTML.

The other viewer tests (test_build_ssn_viewer.py) only assert that the
generated HTML *string* contains the expected markup and function names. They
cannot catch runtime regressions -- e.g. a broken event listener that stops the
viewport from repainting -- because nothing executes the JavaScript.

These tests load the self-contained ``--embed_data`` viewer in a headless
Chromium via Playwright and assert on real rendered behavior. They are gated
behind an optional dependency:

    pip install -e ".[test,browser]"
    playwright install chromium

When the ``playwright`` package is not importable the whole module is skipped;
when the package is present but the browser binary has not been downloaded the
individual tests skip with a hint.
"""

import numpy as np
import pytest

from domainator import build_ssn_viewer
from domainator.data_matrix import DenseDataMatrix

# Skip the entire module if Playwright isn't installed (it lives in the
# optional `browser` extra, not the default `test` extra).
sync_api = pytest.importorskip("playwright.sync_api")


def _build_embedded_viewer(out_dir):
    """Build a self-contained viewer HTML with data embedded, return its path."""
    data = np.array([
        [0, 10, 6, 0, 0, 0],
        [10, 0, 7, 0, 0, 0],
        [6, 7, 0, 4, 0, 0],
        [0, 0, 4, 0, 8, 5],
        [0, 0, 0, 8, 0, 9],
        [0, 0, 0, 5, 9, 0],
    ], dtype=float)
    row_names = ["A", "B", "C", "D", "E", "F"]
    matrix = DenseDataMatrix(data, row_names, row_names)

    input_file = out_dir / "matrix.hdf5"
    html_file = out_dir / "viewer.html"
    matrix.write(str(input_file), output_type="dense")

    build_ssn_viewer.main([
        "-i", str(input_file),
        "--html", str(html_file),
        "--embed_data",
        "--name", "Browser Test Viewer",
    ])
    return html_file


def _build_embedded_viewer_with_metadata(out_dir, color_by="family"):
    """Build a viewer carrying metadata so coloring can be exercised.

    The matrix is the same 6-node A-F graph as ``_build_embedded_viewer``; a
    sidecar metadata TSV adds a categorical column (``family``) and a numeric
    column (``score``). The first TSV column is the node id (consumed as the
    index, left-merged on the matrix row names A-F).
    """
    data = np.array([
        [0, 10, 6, 0, 0, 0],
        [10, 0, 7, 0, 0, 0],
        [6, 7, 0, 4, 0, 0],
        [0, 0, 4, 0, 8, 5],
        [0, 0, 0, 8, 0, 9],
        [0, 0, 0, 5, 9, 0],
    ], dtype=float)
    row_names = ["A", "B", "C", "D", "E", "F"]
    matrix = DenseDataMatrix(data, row_names, row_names)

    input_file = out_dir / "matrix.hdf5"
    html_file = out_dir / "viewer_meta.html"
    meta_file = out_dir / "meta.tsv"
    matrix.write(str(input_file), output_type="dense")
    meta_file.write_text(
        "node_id\tfamily\tscore\n"
        "A\talpha\t1\n"
        "B\talpha\t2\n"
        "C\tbeta\t5\n"
        "D\tbeta\t8\n"
        "E\tgamma\t3\n"
        "F\tgamma\t9\n"
    )

    build_ssn_viewer.main([
        "-i", str(input_file),
        "--html", str(html_file),
        "--embed_data",
        "--name", "Color Test Viewer",
        "--metadata", str(meta_file),
        "--color_by", color_by,
    ])
    return html_file


def _build_embedded_viewer_many_categories(out_dir, node_count=120):
    """Build a viewer whose categorical column has more values than one picker
    page (default page size is 100), so swatch pagination can be exercised."""
    data = np.zeros((node_count, node_count), dtype=float)
    # A simple path so every node is connected into one component.
    for i in range(node_count - 1):
        data[i, i + 1] = 8.0
        data[i + 1, i] = 8.0
    row_names = [f"n{i:03d}" for i in range(node_count)]
    matrix = DenseDataMatrix(data, row_names, row_names)

    input_file = out_dir / "matrix.hdf5"
    html_file = out_dir / "viewer_many.html"
    meta_file = out_dir / "meta.tsv"
    matrix.write(str(input_file), output_type="dense")
    # One unique family per node -> node_count distinct categorical values.
    lines = ["node_id\tfamily"]
    lines += [f"n{i:03d}\tfam_{i:03d}" for i in range(node_count)]
    meta_file.write_text("\n".join(lines) + "\n")

    build_ssn_viewer.main([
        "-i", str(input_file),
        "--html", str(html_file),
        "--embed_data",
        "--name", "Many Category Viewer",
        "--metadata", str(meta_file),
        "--color_by", "family",
    ])
    return html_file


@pytest.fixture(scope="module")
def viewer_html(tmp_path_factory):
    return _build_embedded_viewer(tmp_path_factory.mktemp("ssn_viewer_browser"))


@pytest.fixture(scope="module")
def meta_viewer_html(tmp_path_factory):
    return _build_embedded_viewer_with_metadata(
        tmp_path_factory.mktemp("ssn_viewer_meta")
    )


@pytest.fixture(scope="module")
def many_cat_viewer_html(tmp_path_factory):
    return _build_embedded_viewer_many_categories(
        tmp_path_factory.mktemp("ssn_viewer_many")
    )


def _yield_loaded_page(html_path):
    """Launch headless Chromium, load the viewer, and yield the ready page.

    Exposes ``page.pageerrors`` -- a list of uncaught JS errors captured during
    the test. The previous renderClusterView() event-listener regression threw
    ``ctx.clearRect is not a function`` on every toggle, so an empty list is
    itself a meaningful assertion.
    """
    with sync_api.sync_playwright() as p:
        try:
            browser = p.chromium.launch()
        except sync_api.Error as exc:  # binary not downloaded / sandbox issue
            pytest.skip(
                "Chromium not available for Playwright; run "
                f"'playwright install chromium'. ({exc})"
            )
        page = browser.new_page()
        pageerrors = []
        page.on("pageerror", lambda exc: pageerrors.append(str(exc)))
        page.pageerrors = pageerrors
        page.goto(html_path.as_uri())
        # Default layout is the synchronous "grid" algorithm, so the first
        # render lands as soon as the embedded bundle finishes decompressing.
        # stat-clusters flips from "0" once applyComputedLayout runs.
        page.wait_for_function(
            "() => document.getElementById('stat-clusters').textContent !== '0'"
        )
        try:
            yield page
        finally:
            browser.close()


@pytest.fixture
def page(viewer_html):
    yield from _yield_loaded_page(viewer_html)


@pytest.fixture
def meta_page(meta_viewer_html):
    """A loaded page for the metadata-bearing viewer (categorical default color)."""
    yield from _yield_loaded_page(meta_viewer_html)


@pytest.fixture
def many_cat_page(many_cat_viewer_html):
    """A loaded page whose categorical column spans multiple picker pages."""
    yield from _yield_loaded_page(many_cat_viewer_html)


def _canvas_snapshot(page):
    """Return the cluster canvas pixels as a PNG data URL."""
    return page.eval_on_selector("#cluster-view", "c => c.toDataURL()")


def _wait_for_canvas_change(page, before):
    """Wait until the cluster canvas differs from ``before``.

    Color edits repaint through ``scheduleClusterRender`` (requestAnimationFrame
    coalesced), so the canvas updates on the next frame rather than synchronously;
    polling avoids a race with an immediate snapshot.
    """
    page.wait_for_function(
        "prev => document.getElementById('cluster-view').toDataURL() !== prev",
        arg=before,
    )


def test_viewport_renders_on_load(page):
    snapshot = _canvas_snapshot(page)
    assert snapshot.startswith("data:image/png;base64,")
    # A rendered canvas with several clusters has substantial PNG payload; a
    # blank canvas data URL is tiny by comparison.
    assert len(snapshot) > 5000
    assert page.pageerrors == []


def test_render_nodes_toggle_updates_viewport_immediately(page):
    """Unchecking 'Render nodes' must repaint without needing a pan.

    This is the regression that motivated these tests: the change listener was
    passing the DOM Event as renderClusterView's ctx parameter, so toggles
    silently did nothing (and threw) until a pan called renderClusterView()
    with no args.
    """
    before = _canvas_snapshot(page)
    page.uncheck("#render-nodes")
    after = _canvas_snapshot(page)
    assert after != before, "viewport did not update when toggling 'Render nodes'"

    page.check("#render-nodes")
    restored = _canvas_snapshot(page)
    assert restored != after, "viewport did not update when re-checking 'Render nodes'"
    assert page.pageerrors == []


def test_render_cluster_bounds_toggle_updates_viewport_immediately(page):
    before = _canvas_snapshot(page)
    page.uncheck("#render-cluster-bounds")
    after = _canvas_snapshot(page)
    assert after != before, "viewport did not update when toggling 'Render cluster bounds'"
    assert page.pageerrors == []


def test_split_chart_renders_with_moving_sum(page):
    """The split chart (lollipops + moving-sum trace) renders without errors."""
    # The legend advertises the moving-sum trace.
    assert page.is_visible(".legend-sum")
    snapshot = page.eval_on_selector("#split-chart", "c => c.toDataURL()")
    assert snapshot.startswith("data:image/png;base64,")
    assert len(snapshot) > 5000, "split chart appears blank"
    assert page.pageerrors == []


def test_view_setting_toggles_do_not_throw(page):
    """Toggling every View Settings checkbox should never raise a JS error."""
    for checkbox_id in (
        "show-node-counts",
        "show-edge-scores",
        "render-cluster-bounds",
        "render-nodes",
    ):
        page.click(f"#{checkbox_id}")
        page.click(f"#{checkbox_id}")
    assert page.pageerrors == []


# --- Custom color palette: picker, color-table TSV, and legend export ---


def _open_color_picker(page):
    page.click("#customize-colors")
    page.wait_for_selector("#color-picker-overlay:not([hidden])")


def test_discrete_swatch_recolors(meta_page):
    """Editing a categorical swatch repaints the viewport."""
    page = meta_page
    _open_color_picker(page)
    page.wait_for_selector("#color-picker-discrete:not([hidden])")
    before = _canvas_snapshot(page)
    swatch = page.locator("#color-picker-swatch-list input[type=color]").first
    swatch.fill("#ff0000")
    swatch.dispatch_event("input")
    _wait_for_canvas_change(page, before)
    assert page.pageerrors == []


def test_numeric_stops_recolor(meta_page):
    """Editing low/high colors for a numeric column repaints the viewport."""
    page = meta_page
    page.select_option("#color-by", "score")
    _open_color_picker(page)
    page.wait_for_selector("#color-picker-continuous:not([hidden])")
    before = _canvas_snapshot(page)
    for input_id, value in (("#color-low", "#0000ff"), ("#color-high", "#ffff00")):
        node = page.locator(input_id)
        node.fill(value)
        node.dispatch_event("input")
    _wait_for_canvas_change(page, before)
    assert page.pageerrors == []


def test_discrete_picker_paginates(many_cat_page):
    """A column with >100 values shows a pager that navigates swatch pages."""
    page = many_cat_page
    _open_color_picker(page)
    page.wait_for_selector("#color-picker-discrete:not([hidden])")
    page.wait_for_selector("#color-picker-pager:not([hidden])")
    status = page.locator("#color-picker-page-status").inner_text()
    assert "page 1/2" in status, status
    assert page.is_disabled("#color-picker-prev")
    assert not page.is_disabled("#color-picker-next")
    first_page_top = page.locator("#color-picker-swatch-list .cp-swatch-label").first.inner_text()

    page.click("#color-picker-next")
    page.wait_for_function(
        "() => document.getElementById('color-picker-page-status').textContent.includes('page 2/2')"
    )
    assert page.is_disabled("#color-picker-next")
    assert not page.is_disabled("#color-picker-prev")
    second_page_top = page.locator("#color-picker-swatch-list .cp-swatch-label").first.inner_text()
    assert second_page_top != first_page_top, "Next did not change the visible swatches"
    assert page.pageerrors == []


def test_numeric_bounds_recolor(meta_page):
    """Editing the low/high domain-bound values remaps the gradient."""
    page = meta_page
    page.select_option("#color-by", "score")
    _open_color_picker(page)
    page.wait_for_selector("#color-picker-continuous:not([hidden])")
    before = _canvas_snapshot(page)
    # Tighten the upper bound so mid-range scores saturate differently.
    node = page.locator("#color-high-value")
    node.fill("4")
    node.dispatch_event("input")
    _wait_for_canvas_change(page, before)
    assert page.pageerrors == []


def test_midpoint_shift_without_mid_color(meta_page):
    """Moving the mid value recolors even when 'Use midpoint color' is off."""
    page = meta_page
    page.select_option("#color-by", "score")
    _open_color_picker(page)
    page.wait_for_selector("#color-picker-continuous:not([hidden])")
    # Mid color stays off; the mid-value input must remain editable.
    assert not page.is_checked("#color-mid-toggle")
    assert not page.is_disabled("#color-mid-value")
    before = _canvas_snapshot(page)
    node = page.locator("#color-mid-value")
    node.fill("8")  # shove the midpoint toward the high end
    node.dispatch_event("input")
    _wait_for_canvas_change(page, before)
    assert page.pageerrors == []


def test_reset_values_button(meta_page):
    """Reset values restores the low/high inputs to the data range."""
    page = meta_page
    page.select_option("#color-by", "score")
    _open_color_picker(page)
    page.wait_for_selector("#color-picker-continuous:not([hidden])")
    page.fill("#color-low-value", "2")
    page.locator("#color-low-value").dispatch_event("input")
    page.fill("#color-high-value", "4")
    page.locator("#color-high-value").dispatch_event("input")
    page.click("#color-reset-values")
    # score column ranges 1..9 across nodes A-F.
    assert page.input_value("#color-low-value") == "1"
    assert page.input_value("#color-high-value") == "9"
    assert page.pageerrors == []


def test_save_color_table_download(meta_page):
    """Save color table emits a value<TAB>#HEX TSV."""
    page = meta_page
    _open_color_picker(page)
    with page.expect_download() as download_info:
        page.click("#save-color-table")
    download = download_info.value
    text = open(download.path(), encoding="utf-8").read()
    data_lines = [ln for ln in text.splitlines() if ln.strip()]
    assert data_lines, "saved color table was empty"
    for line in data_lines:
        value, _, hexcode = line.rpartition("\t")
        assert value != "", f"missing value column in {line!r}"
        assert len(hexcode) == 7 and hexcode[0] == "#", f"bad color in {line!r}"
    assert page.pageerrors == []


def test_load_color_table_recolors(meta_page, tmp_path):
    """Loading a color-table TSV recolors the categorical view without errors."""
    page = meta_page
    tsv = tmp_path / "palette.tsv"
    tsv.write_text("alpha\t#ff0000\nbeta\t#00ff00\ngamma\t#0000ff\n")
    _open_color_picker(page)
    page.wait_for_selector("#color-picker-discrete:not([hidden])")
    before = _canvas_snapshot(page)
    page.set_input_files("#color-table-file", str(tsv))
    page.wait_for_function(
        "() => document.querySelector('#color-picker-swatch-list input[type=color]')"
        " && document.querySelector('#color-picker-swatch-list input[type=color]').value === '#ff0000'"
    )
    after = _canvas_snapshot(page)
    assert after != before, "viewport did not update after loading a color table"
    assert page.pageerrors == []


def test_export_legend_svg(meta_page):
    page = meta_page
    _open_color_picker(page)
    with page.expect_download() as download_info:
        page.click("#export-legend-svg")
    download = download_info.value
    assert download.suggested_filename.endswith(".svg")
    content = open(download.path(), encoding="utf-8").read()
    assert "<svg" in content
    assert page.pageerrors == []


def test_export_legend_png(meta_page):
    page = meta_page
    _open_color_picker(page)
    with page.expect_download() as download_info:
        page.click("#export-legend-png")
    download = download_info.value
    assert download.suggested_filename.endswith(".png")
    import os
    assert os.path.getsize(download.path()) > 0
    assert page.pageerrors == []


# --- Treemap layout ---


def _switch_to_treemap(page):
    """Select the treemap layout and wait for the canvas to repaint."""
    before = _canvas_snapshot(page)
    page.select_option("#layout-algorithm", "treemap")
    _wait_for_canvas_change(page, before)


def test_treemap_layout_renders(page):
    """Selecting the treemap layout repaints the viewport without JS errors."""
    _switch_to_treemap(page)
    snapshot = _canvas_snapshot(page)
    assert snapshot.startswith("data:image/png;base64,")
    assert len(snapshot) > 5000, "treemap canvas appears blank"
    assert page.pageerrors == []


def _fragment_threshold_to_multiple_clusters(page):
    """Push the threshold to its maximum so the network splits into several clusters."""
    page.evaluate(
        """() => {
            const slider = document.getElementById('threshold-slider');
            slider.value = slider.max;
            slider.dispatchEvent(new Event('input', {bubbles: true}));
        }"""
    )
    page.wait_for_function(
        "() => Number(document.getElementById('stat-clusters').textContent) > 1"
    )


def test_treemap_toggles_do_not_throw(page):
    """Treemap-mode toggles (render nodes/bounds) repaint cleanly with no JS errors."""
    _switch_to_treemap(page)
    before = _canvas_snapshot(page)
    page.uncheck("#render-nodes")
    _wait_for_canvas_change(page, before)
    page.check("#render-nodes")
    page.click("#render-cluster-bounds")
    page.click("#render-cluster-bounds")
    assert page.pageerrors == []


def test_treemap_nodes_are_fixed_size(page):
    """Every treemap node square is the same fixed world size across all clusters, and
    each tile renders squares for all its members (regression: nodes must always render)."""
    _switch_to_treemap(page)
    stats = page.evaluate(
        """() => {
            const radii = new Set();
            let tilesWithoutSquares = 0;
            state.visibleLayout.forEach(item => {
                const component = state.bundle.graph.hierarchy.nodes[item.componentId];
                const layout = componentMemberLayout(component, item);
                if (layout.length === 0) tilesWithoutSquares += 1;
                layout.forEach(dot => radii.add(Math.round(dot.radius * 1000)));
            });
            return {distinctRadii: radii.size, tilesWithoutSquares};
        }"""
    )
    assert stats["tilesWithoutSquares"] == 0
    assert stats["distinctRadii"] == 1, "treemap node squares must all be the same size"
    assert page.pageerrors == []


def test_treemap_nodes_fixed_across_threshold(page):
    """The core lattice invariant: node positions do not move when the threshold changes;
    only which clusters they belong to changes."""
    _switch_to_treemap(page)

    def lattice_positions():
        return page.evaluate(
            """() => {
                const g = state.latticeGlobal;
                return {w: g.width, h: g.height,
                        cols: Array.from(g.colOf), rows: Array.from(g.rowOf)};
            }"""
        )

    before = lattice_positions()
    cluster_count_before = int(page.locator("#stat-clusters").inner_text())
    _fragment_threshold_to_multiple_clusters(page)
    after = lattice_positions()
    cluster_count_after = int(page.locator("#stat-clusters").inner_text())

    assert cluster_count_after > cluster_count_before, "threshold change should split clusters"
    assert after == before, "node lattice positions must not move when the threshold changes"
    assert page.pageerrors == []


def test_treemap_min_size_never_hides_nodes(page):
    """Raising the minimum cluster size must not hide any nodes (option a): every node is
    still laid out; min size only suppresses small-cluster outlines."""
    _switch_to_treemap(page)
    _fragment_threshold_to_multiple_clusters(page)

    def rendered_node_count():
        return page.evaluate(
            """() => {
                let total = 0;
                state.visibleLayout.forEach(item => {
                    const component = state.bundle.graph.hierarchy.nodes[item.componentId];
                    total += componentMemberLayout(component, item).length;
                });
                return total;
            }"""
        )

    total_nodes = page.evaluate("() => state.bundle.graph.nodes.length")
    assert rendered_node_count() == total_nodes
    # Raising the minimum cluster size well past the largest cluster must not remove any nodes
    # from the layout (it only affects which outlines are drawn).
    page.fill("#min-cluster-size", "9999")
    page.wait_for_timeout(300)
    assert rendered_node_count() == total_nodes, "min cluster size must not hide nodes"
    assert page.pageerrors == []


def test_treemap_svg_export_has_rects(page):
    """SVG export under treemap emits <rect> tiles/members and no JS errors."""
    _switch_to_treemap(page)
    with page.expect_download() as download_info:
        page.click("#export-svg")
    download = download_info.value
    assert download.suggested_filename.endswith(".svg")
    content = open(download.path(), encoding="utf-8").read()
    assert "<svg" in content
    assert "<rect" in content
    assert page.pageerrors == []


def test_treemap_click_selection(page):
    """Ctrl-clicking a treemap node square selects it via the rect hit-test."""
    _switch_to_treemap(page)
    assert page.is_disabled("#clear-selection")
    # Locate the first member square in canvas-CSS coordinates (relative to the
    # canvas element) so the click lands on a node rather than in the inter-square
    # gap of the grid. locator.click(position=...) auto-scrolls the canvas into view.
    point = page.evaluate(
        """() => {
            const item = state.visibleLayout[0];
            const component = state.bundle.graph.hierarchy.nodes[item.componentId];
            const dot = componentMemberLayout(component, item)[0];
            const sp = worldToScreenPoint(dot.x, dot.y);
            const canvas = document.getElementById('cluster-view');
            const rect = canvas.getBoundingClientRect();
            return {
                x: sp.x * (rect.width / canvas.width),
                y: sp.y * (rect.height / canvas.height),
            };
        }"""
    )
    page.locator("#cluster-view").click(
        position={"x": point["x"], "y": point["y"]}, modifiers=["Control"]
    )
    # A successful rect hit-test selects a member, which enables Clear selection.
    page.wait_for_function(
        "() => !document.getElementById('clear-selection').disabled"
    )
    assert page.pageerrors == []
