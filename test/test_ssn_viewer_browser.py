"""Browser-driven regression tests for the generated Domainator Similarity Network Viewer HTML.

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

import re

import numpy as np
import pandas as pd
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


def _build_embedded_viewer_with_metadata(out_dir, color_by="family", categorical=None):
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

    args = [
        "-i", str(input_file),
        "--html", str(html_file),
        "--embed_data",
        "--name", "Color Test Viewer",
        "--metadata", str(meta_file),
        "--color_by", color_by,
    ]
    if categorical:
        args += ["--categorical", *categorical]
    build_ssn_viewer.main(args)
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
def numeric_categorical_viewer_html(tmp_path_factory):
    """A viewer whose numeric color column was marked categorical at build time."""
    return _build_embedded_viewer_with_metadata(
        tmp_path_factory.mktemp("ssn_viewer_numcat"),
        color_by="score",
        categorical=["score"],
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
        # Default layout is the synchronous "packed" algorithm, so the first
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
def numeric_categorical_page(numeric_categorical_viewer_html):
    """A loaded page whose numeric column defaults to discrete coloring."""
    yield from _yield_loaded_page(numeric_categorical_viewer_html)


@pytest.fixture
def many_cat_page(many_cat_viewer_html):
    """A loaded page whose categorical column spans multiple picker pages."""
    yield from _yield_loaded_page(many_cat_viewer_html)



def _build_embedded_viewer_dense_random(out_dir, node_count=60):
    """A 3-block network with many distinct edge weights.

    The equivalence tests need a network whose MST edges mostly have *different*
    weights, so the per-threshold grouping produces many rows and the moving sum
    slides a real window -- the numerically delicate part of the JS port.
    """
    rng = np.random.default_rng(0)
    data = np.zeros((node_count, node_count), dtype=float)
    for start, end in [(0, 22), (22, 40), (40, node_count)]:
        for i in range(start, end):
            for j in range(i + 1, end):
                data[i, j] = data[j, i] = rng.uniform(4, 12)
    for i in range(0, node_count - 1, 17):
        data[i, i + 1] = data[i + 1, i] = 3.0
    row_names = [f"seq_{i:03d}" for i in range(node_count)]
    matrix = DenseDataMatrix(data, row_names, row_names)

    input_file = out_dir / "matrix.hdf5"
    html_file = out_dir / "viewer_dense.html"
    matrix.write(str(input_file), output_type="dense")
    build_ssn_viewer.main([
        "-i", str(input_file),
        "--html", str(html_file),
        "--embed_data",
        "--name", "Dense Test Viewer",
    ])
    return html_file


@pytest.fixture(scope="module")
def dense_viewer_html(tmp_path_factory):
    return _build_embedded_viewer_dense_random(tmp_path_factory.mktemp("ssn_viewer_dense"))


@pytest.fixture
def dense_page(dense_viewer_html):
    """A loaded page whose MST edge weights are nearly all distinct."""
    yield from _yield_loaded_page(dense_viewer_html)


CAPPED_PAGE_MAX_MERGE_EVENTS = 5


def _build_embedded_viewer_capped(out_dir, node_count=60):
    """A dense network with enough split events that a small cap drops most of them.

    Every other fixture sits under the default cap of 500, so nothing in the suite
    would otherwise exercise a filtered split chart -- which is the state every large
    network is in. The cap itself is applied in the viewer (see capped_page): it is a
    control there, not something the bundle was built with.
    """
    rng = np.random.default_rng(0)
    data = np.zeros((node_count, node_count), dtype=float)
    for start, end in [(0, 22), (22, 40), (40, node_count)]:
        for i in range(start, end):
            for j in range(i + 1, end):
                data[i, j] = data[j, i] = rng.uniform(4, 12)
    row_names = [f"seq_{i:03d}" for i in range(node_count)]
    matrix = DenseDataMatrix(data, row_names, row_names)

    input_file = out_dir / "matrix.hdf5"
    html_file = out_dir / "viewer_capped.html"
    matrix.write(str(input_file), output_type="dense")
    build_ssn_viewer.main([
        "-i", str(input_file),
        "--html", str(html_file),
        "--embed_data",
        "--name", "Capped Test Viewer",
    ])
    return html_file


def _set_split_event_cap(page, cap):
    """Move the viewer's split-event cap -- the control that replaced the old
    build-time ``--max_merge_events`` flag."""
    page.fill("#split-event-cap", str(cap))
    page.dispatch_event("#split-event-cap", "change")
    page.wait_for_function("cap => state.maxMergeEvents === cap", arg=cap)


@pytest.fixture(scope="module")
def capped_viewer_html(tmp_path_factory):
    return _build_embedded_viewer_capped(tmp_path_factory.mktemp("ssn_viewer_capped"))


@pytest.fixture
def capped_page(capped_viewer_html):
    """A loaded page whose split chart shows only a capped slice of its events."""
    for page in _yield_loaded_page(capped_viewer_html):
        _set_split_event_cap(page, CAPPED_PAGE_MAX_MERGE_EVENTS)
        yield page


def _build_embedded_viewer_two_bead(out_dir):
    """A network where one threshold carries merges of two different sizes.

    At t=11 two pairs form (a 1-node merge each); at t=10 those two pairs join
    (a 2-node merge) while A-B forms (a 1-node merge). So the t=10 event draws two
    beads on one stem, which is what the hover's bead-picking has to resolve.
    """
    names = ["A", "B", "E", "F", "G", "H"]
    index = {name: position for position, name in enumerate(names)}
    data = np.zeros((6, 6), dtype=float)
    for left, right, weight in [("E", "F", 11.0), ("G", "H", 11.0),
                                ("F", "G", 10.0), ("A", "B", 10.0)]:
        data[index[left], index[right]] = weight
        data[index[right], index[left]] = weight
    matrix = DenseDataMatrix(data, names, names)

    input_file = out_dir / "matrix.hdf5"
    html_file = out_dir / "viewer_two_bead.html"
    matrix.write(str(input_file), output_type="dense")
    build_ssn_viewer.main([
        "-i", str(input_file), "--html", str(html_file), "--embed_data",
        "--name", "Two Bead Viewer",
    ])
    return html_file


@pytest.fixture(scope="module")
def two_bead_viewer_html(tmp_path_factory):
    return _build_embedded_viewer_two_bead(tmp_path_factory.mktemp("ssn_viewer_two_bead"))


@pytest.fixture
def two_bead_page(two_bead_viewer_html):
    """A loaded page with two beads stacked on one split event."""
    yield from _yield_loaded_page(two_bead_viewer_html)


def _build_embedded_viewer_product_metric(out_dir):
    """The dense network under --merge_impact_metric product.

    A product impact is not a count of nodes, so the hover readout must not call it
    one -- the distinction merge_impact_axis_labels exists to keep.
    """
    rng = np.random.default_rng(0)
    node_count = 40
    data = np.zeros((node_count, node_count), dtype=float)
    for start, end in [(0, 15), (15, 28), (28, node_count)]:
        for i in range(start, end):
            for j in range(i + 1, end):
                data[i, j] = data[j, i] = rng.uniform(4, 12)
    names = [f"seq_{i:03d}" for i in range(node_count)]
    matrix = DenseDataMatrix(data, names, names)

    input_file = out_dir / "matrix.hdf5"
    html_file = out_dir / "viewer_product.html"
    matrix.write(str(input_file), output_type="dense")
    build_ssn_viewer.main([
        "-i", str(input_file), "--html", str(html_file), "--embed_data",
        "--name", "Product Metric Viewer", "--merge_impact_metric", "product",
    ])
    return html_file


@pytest.fixture(scope="module")
def product_metric_viewer_html(tmp_path_factory):
    return _build_embedded_viewer_product_metric(
        tmp_path_factory.mktemp("ssn_viewer_product"))


@pytest.fixture
def product_metric_page(product_metric_viewer_html):
    """A loaded page whose merge impacts are size products, not node counts."""
    yield from _yield_loaded_page(product_metric_viewer_html)


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
    """The split chart (stems + beads + moving-sum trace) renders without errors."""
    # The legend advertises both series.
    assert page.is_visible(".legend-sum")
    assert page.is_visible(".legend-line")
    # The moving sum is derived in the browser from the bundle's merge order.
    assert page.evaluate("() => Array.isArray(state.series.movingSum.y)")
    # Stems are drawn at the largest single merge, so every event carries the new keys.
    assert page.evaluate(
        "() => state.series.selectedRows.every("
        "e => typeof e.largest_merge === 'number' && e.merge_size_counts !== undefined)"
    )
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
        "reduce-elongation",
        "leaf-pruning-only",
    ):
        page.click(f"#{checkbox_id}")
        page.click(f"#{checkbox_id}")
    # "Collapse long paths" is only clickable while leaf pruning is on.
    page.check("#leaf-pruning-only")
    page.click("#collapse-long-paths")
    page.click("#collapse-long-paths")
    page.uncheck("#leaf-pruning-only")
    assert page.pageerrors == []


# --- Custom color palette: picker, color-table TSV, and legend export ---


def _open_color_picker(page):
    page.click("#customize-colors")
    page.wait_for_selector("#color-picker-overlay:not([hidden])")


def _stop_field(page, index, field):
    """Locator for one gradient stop row's color or value input."""
    return page.locator(
        f'.cp-stop-row[data-stop-index="{index}"] [data-stop-field="{field}"]')


def _set_stop(page, index, field, value, commit=False):
    """Type into a stop row the way a user would, firing input (and change)."""
    node = _stop_field(page, index, field)
    node.fill(value)
    node.dispatch_event("input")
    if commit:
        node.dispatch_event("change")


def _stop_values(page):
    return page.evaluate("() => state.gradientStops.map(stop => stop.value)")


def _stop_colors(page):
    return page.evaluate("() => state.gradientStops.map(stop => stop.color)")


def test_default_color_by_paints_on_load(meta_page):
    """--color_by colors the nodes on load, before any color control is touched.

    Regression: the node color cache used to be built before the "Color by" menu was
    populated from bundle.defaults, so every node painted the fallback color until
    something else triggered a rebuild.
    """
    page = meta_page
    assert page.input_value("#color-by") == "family"
    colors = page.evaluate("state.nodeColorCache.slice()")
    assert len(set(colors)) == 3, colors
    assert colors == page.evaluate("state.bundle.graph.nodes.map((_, i) => nodeColor(i))")
    assert page.pageerrors == []


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
    _set_stop(page, 0, "color", "#0000ff")
    _set_stop(page, 1, "color", "#ffff00")
    _wait_for_canvas_change(page, before)
    assert _stop_colors(page) == ["#0000ff", "#ffff00"]
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
    _set_stop(page, 1, "value", "4")
    _wait_for_canvas_change(page, before)
    assert _stop_values(page) == [1, 4]
    assert page.pageerrors == []


def test_an_intermediate_stop_shifts_the_ramp(meta_page):
    """Moving an added stop recolors the network."""
    page = meta_page
    page.select_option("#color-by", "score")
    _open_color_picker(page)
    page.wait_for_selector("#color-picker-continuous:not([hidden])")
    page.click("#color-add-stop")
    before = _canvas_snapshot(page)

    _set_stop(page, 1, "value", "8", commit=True)  # shove it toward the high end

    _wait_for_canvas_change(page, before)
    assert _stop_values(page) == [1, 8, 9]
    assert page.pageerrors == []


def test_two_stops_leave_a_straight_ramp(meta_page):
    """With only the two ends, the ramp is a plain interpolation.

    Regression: the gradient used to bend whenever a mid VALUE was stored, and
    that value survived switching the midpoint off. Narrowing the high bound
    below the stale value pinned the bend against the top of the range, so the
    ramp only ever reached the low/high average before jumping to the high
    color. There is now no state to go stale -- a removed stop is gone.
    """
    page = meta_page
    page.select_option("#color-by", "score")
    _open_color_picker(page)
    page.wait_for_selector("#color-picker-continuous:not([hidden])")
    page.click("#color-add-stop")
    page.click('.cp-stop-row[data-stop-index="1"] [data-stop-remove]')
    assert len(_stop_values(page)) == 2

    # Pull the top of the ramp below where the removed stop used to sit.
    _set_stop(page, 1, "value", "3", commit=True)

    high_color = page.evaluate("() => state.gradientStops[1].color")
    assert page.evaluate(
        """() => numericColor(9, state.colorHistogram.min, state.colorHistogram.max,
                              customPalette('score'))""") == high_color
    background = page.eval_on_selector(
        "#color-gradient-preview", "e => getComputedStyle(e).backgroundImage")
    assert background.count("rgb") == 2, background
    assert page.pageerrors == []


def test_reset_values_button(meta_page):
    """Reset values refits the ramp onto the data range."""
    page = meta_page
    page.select_option("#color-by", "score")
    _open_color_picker(page)
    page.wait_for_selector("#color-picker-continuous:not([hidden])")
    _set_stop(page, 0, "value", "2", commit=True)
    _set_stop(page, 1, "value", "4", commit=True)

    page.click("#color-reset-values")

    # score column ranges 1..9 across nodes A-F.
    assert _stop_values(page) == [1, 9]
    assert page.input_value('.cp-stop-row[data-stop-index="0"] [data-stop-field="value"]') == "1"
    assert page.pageerrors == []


def test_reset_values_keeps_the_shape_of_a_tuned_ramp(meta_page):
    """Refitting rescales intermediate stops instead of discarding them."""
    page = meta_page
    page.select_option("#color-by", "score")
    _open_color_picker(page)
    page.wait_for_selector("#color-picker-continuous:not([hidden])")
    page.click("#color-add-stop")
    # Squeeze the whole ramp into 1..5 with the intermediate a quarter of the way in.
    _set_stop(page, 1, "value", "2", commit=True)
    _set_stop(page, 2, "value", "5", commit=True)
    assert _stop_values(page) == [1, 2, 5]

    page.click("#color-reset-values")

    # 1..9 now, with the intermediate still a quarter of the way along.
    assert _stop_values(page) == [1, 3, 9]
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


def test_numeric_column_switches_to_categorical(meta_page):
    """A numeric column can be recolored as discrete categories and switched back."""
    page = meta_page
    page.select_option("#color-by", "score")
    _open_color_picker(page)
    page.wait_for_selector("#color-picker-continuous:not([hidden])")
    page.wait_for_selector("#color-categorical-control:not([hidden])")
    assert not page.is_checked("#color-as-categorical")

    gradient_snapshot = _canvas_snapshot(page)
    page.check("#color-as-categorical")
    page.wait_for_selector("#color-picker-discrete:not([hidden])")
    assert page.is_hidden("#color-picker-continuous")
    # One swatch per distinct score (the fixture has 6 unique values).
    assert page.locator("#color-picker-swatch-list input[type=color]").count() == 6
    categorical_snapshot = _canvas_snapshot(page)
    assert categorical_snapshot != gradient_snapshot, "discrete coloring did not repaint"

    # Recolor one category, then flip back to the gradient: the gradient returns
    # unchanged and the discrete palette survives a second flip.
    swatch = page.locator("#color-picker-swatch-list input[type=color]").first
    swatch.fill("#ff0000")
    swatch.dispatch_event("input")
    _wait_for_canvas_change(page, categorical_snapshot)
    edited_snapshot = _canvas_snapshot(page)

    page.uncheck("#color-as-categorical")
    page.wait_for_selector("#color-picker-continuous:not([hidden])")
    assert _canvas_snapshot(page) == gradient_snapshot

    page.check("#color-as-categorical")
    page.wait_for_selector("#color-picker-discrete:not([hidden])")
    assert page.locator("#color-picker-swatch-list input[type=color]").first.input_value() == "#ff0000"
    assert _canvas_snapshot(page) == edited_snapshot
    assert page.pageerrors == []


def test_bundle_default_categorical_column(numeric_categorical_page):
    """--categorical in the bundle opens the numeric column in discrete mode."""
    page = numeric_categorical_page
    _open_color_picker(page)
    page.wait_for_selector("#color-picker-discrete:not([hidden])")
    assert page.is_checked("#color-as-categorical")
    assert page.locator("#color-picker-swatch-list input[type=color]").count() == 6
    assert page.pageerrors == []


def test_categorical_toggle_hidden_for_string_column(meta_page):
    """The gradient/discrete switch only applies to numeric columns."""
    page = meta_page
    _open_color_picker(page)
    page.wait_for_selector("#color-picker-discrete:not([hidden])")
    assert page.is_hidden("#color-categorical-control")
    assert page.pageerrors == []


def test_save_color_table_for_categorical_numeric(meta_page):
    """A numeric column colored as categories can save/load a value->color TSV."""
    page = meta_page
    page.select_option("#color-by", "score")
    _open_color_picker(page)
    page.check("#color-as-categorical")
    page.wait_for_selector("#color-picker-discrete:not([hidden])")
    with page.expect_download() as download_info:
        page.click("#save-color-table")
    text = open(download_info.value.path(), encoding="utf-8").read()
    keys = [line.split("\t")[0] for line in text.splitlines() if line.strip()]
    assert "1" in keys and "9" in keys, keys
    assert page.pageerrors == []


def _saved_color_table(page):
    """Click "Save color table" and return the downloaded value -> color mapping."""
    with page.expect_download() as download_info:
        page.click("#save-color-table")
    text = open(download_info.value.path(), encoding="utf-8").read()
    table = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        value, _, hexcode = line.rpartition("\t")
        table[value] = hexcode
    return table


def test_named_palette_matches_get_palette(meta_page):
    """The "Domainator distinct" palette assigns the colors build_ssn.py would.

    It is also the viewer's default, so an untouched column already carries those
    colors and picking it from the menu only writes them down -- which is why this
    checks the color table before and after, and not the canvas (it cannot change).
    """
    from domainator.utils import get_palette

    page = meta_page
    _open_color_picker(page)
    page.wait_for_selector("#color-picker-discrete:not([hidden])")

    expected = get_palette(pd.Series(["alpha", "alpha", "beta", "beta", "gamma", "gamma"]))
    assert page.evaluate("() => Object.keys(state.customPalettes)") == []
    default_table = _saved_color_table(page)
    assert {key: value.upper() for key, value in default_table.items()
            if key != "\u2014"} == expected
    assert page.input_value("#color-palette") == "domainator"

    page.select_option("#color-palette", "domainator")
    page.wait_for_function(
        "() => !!(state.customPalettes.family && state.customPalettes.family.scheme)"
    )
    saved = _saved_color_table(page)
    assert {key: value.upper() for key, value in saved.items() if key != "\u2014"} == expected
    assert page.pageerrors == []


def test_named_palette_orders_numeric_values_numerically(meta_page):
    """A numeric column colored as categories matches get_palette's sorted assignment."""
    from domainator.utils import get_palette

    page = meta_page
    page.select_option("#color-by", "score")
    _open_color_picker(page)
    page.check("#color-as-categorical")
    page.wait_for_selector("#color-picker-discrete:not([hidden])")
    page.select_option("#color-palette", "domainator")
    page.wait_for_function(
        "() => document.getElementById('color-palette-note').textContent.includes('64 colors')"
    )

    expected = {str(key): color for key, color in get_palette(pd.Series([1, 2, 5, 8, 3, 9])).items()}
    saved = _saved_color_table(page)
    assert {key: value.upper() for key, value in saved.items() if key != "\u2014"} == expected
    assert page.pageerrors == []


def test_named_palette_cycles_and_reverts_to_default(many_cat_page):
    """A short palette cycles across many values, and the default palette comes back.

    The menu carries no pseudo-entry for the default any more, so "Reset to
    defaults" is the only way back -- and what it comes back to is
    DEFAULT_CATEGORICAL_PALETTE, not a hash of each value.
    """
    page = many_cat_page
    _open_color_picker(page)
    page.wait_for_selector("#color-picker-discrete:not([hidden])")
    # The per-node color cache is the layout-independent record of what is painted.
    default_colors = page.evaluate("state.nodeColorCache.slice()")
    before = _canvas_snapshot(page)

    page.select_option("#color-palette", "okabe_ito")
    _wait_for_canvas_change(page, before)
    note = page.locator("#color-palette-note").inner_text()
    assert "colors repeat" in note, note
    saved = _saved_color_table(page)
    # 120 values over an 8-color palette: 8 distinct colors, each reused.
    assert len(saved) == 120
    assert len(set(saved.values())) == 8
    assert page.evaluate("state.nodeColorCache.slice()") != default_colors

    # A hand edit stops the menu from claiming the column is still that palette.
    swatch = page.locator("#color-picker-swatch-list input[type=color]").first
    swatch.fill("#123456")
    swatch.dispatch_event("input")
    page.click("#color-picker-close")
    _open_color_picker(page)
    assert page.input_value("#color-palette") == ""

    assert page.evaluate(
        "() => Array.from(document.getElementById('color-palette').options)"
        ".map(option => option.value)"
    ) == ["", "domainator", "tab10", "okabe_ito", "brewer_dark2", "brewer_set2", "brewer_paired"]

    page.click("#color-picker-reset")
    page.wait_for_selector("#color-picker-discrete:not([hidden])")
    assert page.evaluate("Object.keys(state.customPalettes)") == []
    assert page.evaluate("state.nodeColorCache.slice()") == default_colors
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


def test_treemap_click_in_node_gap_still_selects(page):
    """A click in the padding gap between nodes still selects (the whole lattice cell is
    clickable, not just the drawn node square) -- regression for gap clicks falling through."""
    _switch_to_treemap(page)
    assert page.is_disabled("#clear-selection")
    # Shift the click into the padding: a node square is TREEMAP_NODE (13) wide, so a point
    # 7.5px right of a member center lands in the inter-node gap, outside the drawn square.
    point = page.evaluate(
        """() => {
            const item = state.visibleLayout[0];
            const component = state.bundle.graph.hierarchy.nodes[item.componentId];
            const dot = componentMemberLayout(component, item)[0];
            const sp = worldToScreenPoint(dot.x + 7.5, dot.y);
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
    page.wait_for_function(
        "() => !document.getElementById('clear-selection').disabled"
    )
    assert page.pageerrors == []


# ---------------------------------------------------------------------------
# Session save/load (bundle v4 ``app_state``)
# ---------------------------------------------------------------------------


def _save_session(page, out_dir, filename="session.dsnv"):
    """Click "Save session…" and return the path the download was saved to."""
    with page.expect_download() as download_info:
        page.click("#save-session")
    target = out_dir / filename
    download_info.value.save_as(str(target))
    return target


def _load_bundle_file(page, path):
    """Feed a bundle file to the viewer's file picker and wait for it to install."""
    page.set_input_files("#bundle-file", str(path))
    page.wait_for_function(
        "() => document.getElementById('bundle-status').textContent.startsWith('Loaded ')"
    )
    page.wait_for_function(
        "() => document.getElementById('stat-clusters').textContent !== '0'"
    )


def _read_session_bundle(path):
    """Decompress a saved .dsnv session and return the parsed JSON."""
    import gzip
    import json

    raw = open(path, "rb").read()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return json.loads(raw.decode("utf-8"))


def test_saved_session_is_a_valid_current_bundle(meta_page, tmp_path):
    """Saving produces a bundle Python readers accept, carrying an app_state."""
    from domainator.ssn_bundle import (
        SSN_VIEWER_BUNDLE_VERSION,
        load_bundle,
    )

    page = meta_page
    saved = _save_session(page, tmp_path)

    bundle = load_bundle(saved)  # would raise on a bad format/version
    assert bundle["version"] == SSN_VIEWER_BUNDLE_VERSION == 6
    assert bundle["graph"]["nodes"] == ["A", "B", "C", "D", "E", "F"]

    app_state = bundle["app_state"]
    assert app_state["state_version"] == 1
    assert app_state["saved_by"]
    # The registry groups fields into sections; spot-check one from each.
    assert app_state["view"]["color_by"] == "family"
    assert app_state["table"]["rows_per_page"] == "250"
    assert app_state["colors"]["categorical_columns"] == []
    assert app_state["selection"]["node_ids"] == []
    assert page.pageerrors == []


def test_session_round_trip_restores_view_and_selection(meta_page, tmp_path):
    """View settings, the threshold, and the selection survive save -> load."""
    page = meta_page

    page.select_option("#color-by", "score")
    page.select_option("#label-by", "family")
    page.select_option("#layout-algorithm", "grid")
    page.uncheck("#show-node-counts")
    page.fill("#metadata-filter", "alpha")
    page.select_option("#metadata-rows-per-page", "500")
    # Select every node in the graph so the selection has something to restore.
    page.evaluate(
        "() => { state.selectedNodeIndices = new Set([0, 1, 2]);"
        " renderClusterView(); updateMetadataTable(); }"
    )
    page.wait_for_function("() => state.selectedNodeIndices.size === 3")

    saved = _save_session(page, tmp_path)
    _load_bundle_file(page, saved)

    assert page.eval_on_selector("#color-by", "e => e.value") == "score"
    assert page.eval_on_selector("#label-by", "e => e.value") == "family"
    assert page.eval_on_selector("#layout-algorithm", "e => e.value") == "grid"
    assert page.eval_on_selector("#show-node-counts", "e => e.checked") is False
    assert page.eval_on_selector("#metadata-filter", "e => e.value") == "alpha"
    assert page.eval_on_selector("#metadata-rows-per-page", "e => e.value") == "500"
    assert page.evaluate("() => Array.from(state.selectedNodeIndices).sort()") == [0, 1, 2]
    assert page.pageerrors == []


def test_session_restores_threshold_by_value_not_slider_position(meta_page, tmp_path):
    """The threshold persists as a value, so it survives independently of stops."""
    page = meta_page
    page.evaluate(
        "() => { const stops = state.sliderModel.stops.filter(s => s.threshold_value !== null);"
        " snapSliderToStop(stops[stops.length - 1]); updateThresholdUI(false); }"
    )
    before = page.evaluate("() => selectedThresholdValue()")
    assert before != float("inf")

    saved = _save_session(page, tmp_path)
    assert _read_session_bundle(saved)["app_state"]["view"]["threshold_value"] == before

    _load_bundle_file(page, saved)
    assert page.evaluate("() => selectedThresholdValue()") == before
    assert page.pageerrors == []


def test_session_restores_pan_and_zoom(meta_page, tmp_path):
    """A restored view transform is not clobbered by the layout's auto-fit."""
    page = meta_page
    page.evaluate(
        "() => { state.viewTransform.scale = 2.5;"
        " state.viewTransform.offsetX = 123; state.viewTransform.offsetY = -45;"
        " renderClusterView(); }"
    )
    saved = _save_session(page, tmp_path)
    _load_bundle_file(page, saved)

    transform = page.evaluate(
        "() => ({scale: state.viewTransform.scale,"
        " offsetX: state.viewTransform.offsetX, offsetY: state.viewTransform.offsetY})"
    )
    assert transform == {"scale": 2.5, "offsetX": 123, "offsetY": -45}
    assert page.evaluate("() => state.pendingRestoreViewTransform") is None
    assert page.pageerrors == []


def test_session_with_unknown_keys_loads_and_reports(meta_page, tmp_path):
    """Unknown or unusable saved settings are skipped and reported, not fatal.

    This is the forward/backward-compatibility contract: the viewer's options
    are expected to change between versions, so a session written by a
    different build must still open.
    """
    import gzip
    import json

    page = meta_page
    saved = _save_session(page, tmp_path)
    bundle = _read_session_bundle(saved)
    bundle["app_state"]["view"]["a_setting_from_the_future"] = 42
    bundle["app_state"]["view"]["color_by"] = "no_such_column"
    bundle["app_state"]["an_entire_unknown_section"] = {"x": 1}
    edited = tmp_path / "edited.dsnv"
    edited.write_bytes(gzip.compress(json.dumps(bundle).encode("utf-8")))

    _load_bundle_file(page, edited)

    status = page.eval_on_selector("#bundle-status", "e => e.textContent")
    assert "skipped 3 saved settings" in status
    assert "view.a_setting_from_the_future (unknown)" in status
    assert 'view.color_by (no option "no_such_column")' in status
    # A bad color_by must not blank the control -- it keeps the bundle default.
    assert page.eval_on_selector("#color-by", "e => e.value") == "family"
    assert page.pageerrors == []


def test_session_reports_node_ids_missing_from_bundle(meta_page, tmp_path):
    """Selections persist as node ids; ids absent from the bundle are counted."""
    import gzip
    import json

    page = meta_page
    saved = _save_session(page, tmp_path)
    bundle = _read_session_bundle(saved)
    bundle["app_state"]["selection"]["node_ids"] = ["A", "C", "GONE_1", "GONE_2"]
    edited = tmp_path / "missing.dsnv"
    edited.write_bytes(gzip.compress(json.dumps(bundle).encode("utf-8")))

    _load_bundle_file(page, edited)

    status = page.eval_on_selector("#bundle-status", "e => e.textContent")
    assert "2 saved node ids are not in this bundle" in status
    assert page.evaluate("() => Array.from(state.selectedNodeIndices).sort()") == [0, 2]
    assert page.pageerrors == []


def test_unknown_bundle_version_loads_with_a_warning(meta_page, tmp_path):
    """A wrong format is fatal, but an unrecognized version only warns."""
    import gzip
    import json

    page = meta_page
    saved = _save_session(page, tmp_path)
    bundle = _read_session_bundle(saved)
    bundle["version"] = 99
    edited = tmp_path / "v99.dsnv"
    edited.write_bytes(gzip.compress(json.dumps(bundle).encode("utf-8")))

    _load_bundle_file(page, edited)

    status = page.eval_on_selector("#bundle-status", "e => e.textContent")
    assert "bundle version 99 is not one this viewer knows about" in status
    assert page.evaluate("() => state.bundle.graph.nodes.length") == 6
    assert page.pageerrors == []


# ---------------------------------------------------------------------------
# Selection presets
# ---------------------------------------------------------------------------


def _select_nodes(page, node_indices):
    """Set the graph selection directly and wait for it to take effect."""
    page.evaluate(
        "indices => { state.selectedNodeIndices = new Set(indices);"
        " renderClusterView(); updateMetadataTable(); }",
        node_indices,
    )
    page.wait_for_function(
        "n => state.selectedNodeIndices.size === n", arg=len(node_indices)
    )


def test_preset_store_and_recall_by_keyboard(page):
    """Shift+digit stores the selection; the bare digit recalls it."""
    _select_nodes(page, [0, 1, 2])
    page.keyboard.press("Shift+Digit3")
    assert page.evaluate("() => state.selectionPresets.get(3).nodeIndices.size") == 3

    page.click("#clear-selection")
    page.wait_for_function("() => state.selectedNodeIndices.size === 0")

    page.keyboard.press("Digit3")
    assert page.evaluate("() => Array.from(state.selectedNodeIndices).sort()") == [0, 1, 2]
    assert page.pageerrors == []


def test_preset_slot_button_shows_filled_state_and_recalls_on_click(page):
    _select_nodes(page, [4, 5])
    page.keyboard.press("Shift+Digit7")

    slot = page.query_selector('[data-preset-slot="7"]')
    assert "preset-slot-filled" in slot.get_attribute("class")
    assert "2 nodes" in slot.get_attribute("title")
    # An untouched slot stays empty-looking.
    empty = page.query_selector('[data-preset-slot="1"]')
    assert "preset-slot-filled" not in empty.get_attribute("class")

    page.click("#clear-selection")
    page.wait_for_function("() => state.selectedNodeIndices.size === 0")
    page.click('[data-preset-slot="7"]')
    assert page.evaluate("() => Array.from(state.selectedNodeIndices).sort()") == [4, 5]
    assert page.pageerrors == []


def test_storing_an_empty_selection_clears_the_slot(page):
    _select_nodes(page, [0])
    page.keyboard.press("Shift+Digit2")
    assert page.evaluate("() => state.selectionPresets.has(2)") is True

    page.click("#clear-selection")
    page.wait_for_function("() => state.selectedNodeIndices.size === 0")
    page.keyboard.press("Shift+Digit2")

    assert page.evaluate("() => state.selectionPresets.has(2)") is False
    assert "Cleared preset 2" in page.eval_on_selector("#bundle-status", "e => e.textContent")
    assert page.pageerrors == []


def test_digit_shortcuts_ignore_typing_in_inputs(page):
    """Typing a digit into a text field must not recall or store a preset."""
    _select_nodes(page, [0, 1])
    page.keyboard.press("Shift+Digit5")
    page.click("#clear-selection")

    page.click("#metadata-filter")
    page.keyboard.press("Digit5")
    page.keyboard.press("Shift+Digit5")

    assert page.eval_on_selector("#metadata-filter", "e => e.value") == "5%"
    assert page.evaluate("() => state.selectedNodeIndices.size") == 0
    assert page.evaluate("() => state.selectionPresets.get(5).nodeIndices.size") == 2
    assert page.pageerrors == []


def test_hovering_a_preset_slot_outlines_its_nodes(page):
    """Hovering previews the preset in gray without changing the selection."""
    _select_nodes(page, [0, 1, 2])
    page.keyboard.press("Shift+Digit4")
    page.click("#clear-selection")
    page.wait_for_function("() => state.selectedNodeIndices.size === 0")

    before = _canvas_snapshot(page)
    page.hover('[data-preset-slot="4"]')
    page.wait_for_function("() => state.presetPreviewSlot === 4")
    _wait_for_canvas_change(page, before)
    # Previewing must not touch the live selection.
    assert page.evaluate("() => state.selectedNodeIndices.size") == 0

    hovered = _canvas_snapshot(page)
    page.hover("#hidden-summary")
    page.wait_for_function("() => state.presetPreviewSlot === null")
    _wait_for_canvas_change(page, hovered)
    assert _canvas_snapshot(page) == before
    assert page.pageerrors == []


def test_hovering_an_empty_preset_slot_does_not_repaint(page):
    before = _canvas_snapshot(page)
    page.hover('[data-preset-slot="8"]')
    page.wait_for_function("() => state.presetPreviewSlot === 8")
    assert _canvas_snapshot(page) == before
    assert page.pageerrors == []


def test_png_export_excludes_the_hover_preview(page):
    """The preview is a transient cue, so exports must not bake it in."""
    _select_nodes(page, [0, 1, 2])
    page.keyboard.press("Shift+Digit6")
    page.click("#clear-selection")
    page.wait_for_function("() => state.selectedNodeIndices.size === 0")

    plain = page.evaluate(
        "() => { const c = document.createElement('canvas');"
        " c.width = clusterCanvas.width; c.height = clusterCanvas.height;"
        " renderClusterView(c.getContext('2d'), c.width, c.height, {preview: false});"
        " return c.toDataURL(); }"
    )
    page.hover('[data-preset-slot="6"]')
    page.wait_for_function("() => state.presetPreviewSlot === 6")
    with_preview_suppressed = page.evaluate(
        "() => { const c = document.createElement('canvas');"
        " c.width = clusterCanvas.width; c.height = clusterCanvas.height;"
        " renderClusterView(c.getContext('2d'), c.width, c.height, {preview: false});"
        " return c.toDataURL(); }"
    )
    assert with_preview_suppressed == plain
    assert page.pageerrors == []


def test_presets_survive_a_session_round_trip(meta_page, tmp_path):
    """Presets persist as node ids inside app_state."""
    page = meta_page
    _select_nodes(page, [1, 3])
    page.keyboard.press("Shift+Digit9")

    saved = _save_session(page, tmp_path, "presets.dsnv")
    presets = _read_session_bundle(saved)["app_state"]["selection"]["presets"]
    assert presets["9"]["node_ids"] == ["B", "D"]

    _load_bundle_file(page, saved)
    assert page.evaluate("() => Array.from(state.selectionPresets.keys())") == [9]
    assert page.evaluate("() => Array.from(state.selectionPresets.get(9).nodeIndices).sort()") == [1, 3]
    slot = page.query_selector('[data-preset-slot="9"]')
    assert "preset-slot-filled" in slot.get_attribute("class")
    assert page.pageerrors == []


# ---------------------------------------------------------------------------
# Metadata editing
# ---------------------------------------------------------------------------


def _cell(page, node_index, column):
    """Return the <td> for one node/column of the rendered table."""
    return page.query_selector(
        f'#metadata-table tbody tr[data-node-index="{node_index}"] '
        f'td[data-column-key="{column}"]'
    )


def _edit_cell(page, node_index, column, text, commit_key="Enter"):
    _cell(page, node_index, column).dblclick()
    page.wait_for_selector(".metadata-cell-input")
    page.fill(".metadata-cell-input", text)
    page.keyboard.press(commit_key)


def _row_values(page, node_index):
    return page.evaluate(
        "i => state.metadataByNodeIndex[i].slice()", node_index
    )


def _open_edit_panel(page, name):
    """Expand one of the Add/Set/Delete column disclosure panels."""
    if page.eval_on_selector(f"#metadata-panel-{name}", "e => e.getAttribute('aria-expanded')") != "true":
        page.click(f"#metadata-panel-{name}")
    page.wait_for_selector(f"#metadata-{name}-panel:not([hidden])")


def _add_column(page, name, column_type="text"):
    _open_edit_panel(page, "add")
    page.fill("#metadata-new-column-name", name)
    page.select_option("#metadata-new-column-type", column_type)
    page.click("#metadata-add-column")


def test_double_click_edits_a_metadata_cell(meta_page):
    page = meta_page
    assert _row_values(page, 0) == ["alpha", 1]

    _edit_cell(page, 0, "family", "delta")

    assert _row_values(page, 0) == ["delta", 1]
    assert _cell(page, 0, "family").inner_text() == "delta"
    # node_id is the graph key and must stay read-only.
    assert _cell(page, 0, "node_id") is None
    assert page.pageerrors == []


def test_escape_cancels_a_cell_edit(meta_page):
    page = meta_page
    _edit_cell(page, 1, "family", "should_not_stick", commit_key="Escape")
    assert _row_values(page, 1) == ["alpha", 2]
    assert page.pageerrors == []


def test_editing_a_numeric_cell_rejects_non_numbers(meta_page):
    page = meta_page
    _edit_cell(page, 2, "score", "77")
    assert _row_values(page, 2) == ["beta", 77]

    _edit_cell(page, 2, "score", "not a number")
    assert _row_values(page, 2) == ["beta", 77]
    status = page.eval_on_selector("#bundle-status", "e => e.textContent")
    assert '"not a number" is not a number' in status
    assert page.pageerrors == []


def test_fractional_input_widens_an_int_column_to_float(meta_page):
    """The int/float split is inferred from the source values, so a fractional
    edit widens the column rather than silently truncating to an integer."""
    page = meta_page
    assert page.evaluate("() => metadataColumnType('score')") == "int"

    _edit_cell(page, 2, "score", "42.5")

    assert _row_values(page, 2) == ["beta", 42.5]
    assert page.evaluate("() => metadataColumnType('score')") == "float"
    assert "widened to float" in page.eval_on_selector("#bundle-status", "e => e.textContent")
    assert page.pageerrors == []


def test_blanking_a_cell_clears_it_to_null(meta_page):
    page = meta_page
    _edit_cell(page, 3, "family", "")
    assert _row_values(page, 3)[0] is None
    assert _cell(page, 3, "family").inner_text() == "—"
    assert page.pageerrors == []


def test_editing_a_cell_updates_node_colors(meta_page):
    """An edit must invalidate the color cache, not just the table.

    The default palette hands out colors in sorted value order, so introducing
    "a_brand_new_family" -- which sorts first -- shifts every other family along
    by one. Node 0 keeps the first color precisely because it is the renamed one;
    node 1, still "alpha", is the node whose color has to move.
    """
    page = meta_page
    before = page.evaluate("() => state.nodeColorCache.slice()")
    _edit_cell(page, 0, "family", "a_brand_new_family")
    after = page.evaluate("() => state.nodeColorCache.slice()")
    assert after != before
    assert after[1] != before[1]
    assert page.pageerrors == []


def test_editing_a_cell_updates_the_search_index(meta_page):
    """rebuildMetadataCaches also maintains the filter's search text."""
    page = meta_page
    _edit_cell(page, 0, "family", "zzz_unique_marker")
    page.fill("#metadata-filter", "zzz_unique_marker")
    page.wait_for_function(
        "() => document.querySelectorAll('#metadata-table tbody tr[data-node-index]').length === 1"
    )
    assert _cell(page, 0, "family").inner_text() == "zzz_unique_marker"
    assert page.pageerrors == []


def test_add_and_delete_a_metadata_column(meta_page):
    page = meta_page
    _add_column(page, "my_notes")

    assert page.evaluate(
        "() => state.metadataColumns.map(c => c.name)"
    ) == ["family", "score", "my_notes"]
    assert page.evaluate("() => state.metadataByNodeIndex[0].length") == 3
    # New columns join every consumer of the column list.
    assert "my_notes" in page.eval_on_selector_all(
        "#color-by option", "opts => opts.map(o => o.value)"
    )
    assert _cell(page, 0, "my_notes") is not None

    _edit_cell(page, 0, "my_notes", "a note")
    assert _row_values(page, 0)[2] == "a note"

    _open_edit_panel(page, "delete")
    page.select_option("#metadata-delete-column", "my_notes")
    page.click("#metadata-delete-column-apply")
    assert page.evaluate("() => state.metadataColumns.map(c => c.name)") == ["family", "score"]
    assert page.evaluate("() => state.metadataByNodeIndex[0].length") == 2
    assert page.pageerrors == []


def test_delete_menu_lists_every_column_labelled_by_origin(meta_page):
    """Bundle columns are deletable too; the menu says where each came from."""
    page = meta_page
    _add_column(page, "scratch")
    _open_edit_panel(page, "delete")

    assert page.eval_on_selector_all(
        "#metadata-delete-column option", "opts => opts.map(o => o.value)"
    ) == ["family", "score", "scratch"]
    assert page.eval_on_selector_all(
        "#metadata-delete-column option", "opts => opts.map(o => o.textContent)"
    ) == ["family (from bundle)", "score (from bundle)", "scratch (added here)"]

    page.select_option("#metadata-delete-column", "scratch")
    assert page.eval_on_selector("#metadata-delete-note", "e => e.textContent") == (
        "Created in the viewer."
    )
    page.select_option("#metadata-delete-column", "family")
    assert "reload it to restore" in page.eval_on_selector(
        "#metadata-delete-note", "e => e.textContent"
    )
    assert page.pageerrors == []


def test_deleting_a_bundle_column_asks_first(meta_page):
    """Dropping source data is the one un-undoable edit, so it prompts."""
    page = meta_page
    _open_edit_panel(page, "delete")
    page.select_option("#metadata-delete-column", "family")

    page.once("dialog", lambda dialog: dialog.dismiss())
    page.click("#metadata-delete-column-apply")
    assert page.evaluate("() => state.metadataColumns.map(c => c.name)") == ["family", "score"]

    page.once("dialog", lambda dialog: dialog.accept())
    page.click("#metadata-delete-column-apply")
    assert page.evaluate("() => state.metadataColumns.map(c => c.name)") == ["score"]
    assert page.evaluate("() => state.metadataByNodeIndex[0]") == [1]
    # The column is gone from every consumer, not just the table.
    assert "family" not in page.eval_on_selector_all(
        "#color-by option", "opts => opts.map(o => o.value)"
    )
    assert page.pageerrors == []


def test_deleting_a_viewer_added_column_does_not_prompt(meta_page):
    """Columns created here are cheap to recreate, so they go without a prompt."""
    page = meta_page
    _add_column(page, "scratch")
    _open_edit_panel(page, "delete")
    page.select_option("#metadata-delete-column", "scratch")

    dialogs = []
    page.on("dialog", lambda dialog: (dialogs.append(dialog), dialog.accept()))
    page.click("#metadata-delete-column-apply")

    assert dialogs == []
    assert page.evaluate("() => state.metadataColumns.map(c => c.name)") == ["family", "score"]
    assert page.pageerrors == []


def test_bulk_fill_applies_to_all_nodes(meta_page):
    """The "All nodes" target ignores the selection and the current page."""
    page = meta_page
    _add_column(page, "everywhere")
    _select_nodes(page, [0])

    _open_edit_panel(page, "set")
    page.select_option("#metadata-fill-column", "everywhere")
    page.select_option("#metadata-fill-target", "all")
    page.fill("#metadata-fill-value", "all_of_them")
    page.click("#metadata-fill-apply")

    assert [_row_values(page, i)[2] for i in range(6)] == ["all_of_them"] * 6
    assert "on 6 nodes in the network" in page.eval_on_selector(
        "#bundle-status", "e => e.textContent"
    )
    assert page.pageerrors == []


def test_edit_panels_open_one_at_a_time(meta_page):
    """Only one editing panel is expanded at a time, and clicking again closes it."""
    page = meta_page

    def expanded():
        return [
            name for name in ("add", "set", "delete")
            if page.eval_on_selector(f"#metadata-panel-{name}", "e => e.getAttribute('aria-expanded')") == "true"
        ]

    assert expanded() == []
    assert page.eval_on_selector("#metadata-add-panel", "e => e.hidden") is True

    page.click("#metadata-panel-add")
    assert expanded() == ["add"]
    assert page.eval_on_selector("#metadata-add-panel", "e => e.hidden") is False

    page.click("#metadata-panel-set")
    assert expanded() == ["set"]
    assert page.eval_on_selector("#metadata-add-panel", "e => e.hidden") is True

    page.click("#metadata-panel-set")
    assert expanded() == []
    assert page.pageerrors == []


def test_loading_a_bundle_collapses_the_edit_panels(meta_page, tmp_path):
    page = meta_page
    page.click("#metadata-panel-add")
    saved = _save_session(page, tmp_path, "panels.dsnv")
    _load_bundle_file(page, saved)

    assert page.eval_on_selector("#metadata-panel-add", "e => e.getAttribute('aria-expanded')") == "false"
    assert page.eval_on_selector("#metadata-add-panel", "e => e.hidden") is True
    assert page.pageerrors == []


def test_adding_a_column_keeps_the_current_color_by(meta_page):
    """populateMetadataControls resets to bundle defaults; adding must not."""
    page = meta_page
    page.select_option("#color-by", "score")
    _add_column(page, "extra")
    assert page.eval_on_selector("#color-by", "e => e.value") == "score"
    assert page.pageerrors == []


def test_duplicate_column_names_are_rejected(meta_page):
    page = meta_page
    _add_column(page, "family")
    assert page.evaluate("() => state.metadataColumns.length") == 2
    assert 'A column named "family" already exists' in page.eval_on_selector(
        "#bundle-status", "e => e.textContent"
    )
    assert page.pageerrors == []


def test_bulk_fill_applies_to_the_graph_selection(meta_page):
    """The main annotation workflow: select a cluster, then label it."""
    page = meta_page
    _add_column(page, "label")
    _select_nodes(page, [0, 1])

    _open_edit_panel(page, "set")
    page.select_option("#metadata-fill-column", "label")
    page.select_option("#metadata-fill-target", "selection")
    page.fill("#metadata-fill-value", "candidate")
    page.click("#metadata-fill-apply")

    assert [_row_values(page, i)[2] for i in range(6)] == [
        "candidate", "candidate", None, None, None, None
    ]
    assert page.pageerrors == []


def test_bulk_fill_applies_to_staged_table_rows(meta_page):
    page = meta_page
    _add_column(page, "label")
    page.evaluate(
        "() => { state.selectedMetadataNodeIndices = new Set([2, 4]);"
        " applyMetadataTableRowHighlights(); }"
    )

    _open_edit_panel(page, "set")
    page.select_option("#metadata-fill-column", "label")
    page.select_option("#metadata-fill-target", "rows")
    page.fill("#metadata-fill-value", "staged")
    page.click("#metadata-fill-apply")

    assert [_row_values(page, i)[2] for i in range(6)] == [
        None, None, "staged", None, "staged", None
    ]
    assert page.pageerrors == []


def test_bulk_fill_reports_an_empty_target(meta_page):
    page = meta_page
    _add_column(page, "label")
    _open_edit_panel(page, "set")
    page.select_option("#metadata-fill-column", "label")
    page.select_option("#metadata-fill-target", "rows")
    page.fill("#metadata-fill-value", "x")
    page.click("#metadata-fill-apply")
    assert "No staged table rows to fill" in page.eval_on_selector(
        "#bundle-status", "e => e.textContent"
    )
    assert page.pageerrors == []


def test_paste_column_fills_the_visible_page(meta_page):
    page = meta_page
    _add_column(page, "pasted")

    _open_edit_panel(page, "set")
    page.click("#metadata-paste-open")
    page.wait_for_selector("#metadata-paste-overlay:not([hidden])")
    assert page.eval_on_selector("#metadata-paste-count", "e => e.textContent") == "6"
    page.select_option("#metadata-paste-column", "pasted")
    page.fill("#metadata-paste-values", "p0\np1\np2\np3\np4\np5\n")
    page.click("#metadata-paste-apply")

    page.wait_for_function("() => document.getElementById('metadata-paste-overlay').hidden")
    assert [_row_values(page, i)[2] for i in range(6)] == ["p0", "p1", "p2", "p3", "p4", "p5"]
    assert page.pageerrors == []


def test_paste_column_refuses_a_count_mismatch(meta_page):
    """A short paste means the page moved since the copy; filling a prefix
    would silently mislabel rows, so it is refused."""
    page = meta_page
    _add_column(page, "pasted")

    _open_edit_panel(page, "set")
    page.click("#metadata-paste-open")
    page.select_option("#metadata-paste-column", "pasted")
    page.fill("#metadata-paste-values", "only\ntwo")
    page.click("#metadata-paste-apply")

    status = page.eval_on_selector("#bundle-status", "e => e.textContent")
    assert "Paste has 2 values but this page shows 6 rows" in status
    assert [_row_values(page, i)[2] for i in range(6)] == [None] * 6
    # The dialog stays open so the paste can be corrected.
    assert page.eval_on_selector("#metadata-paste-overlay", "e => e.hidden") is False
    assert page.pageerrors == []


def test_paste_round_trips_a_copied_column(meta_page):
    """Paste is the inverse of the per-column Copy button."""
    page = meta_page
    _add_column(page, "copy_target")

    copied = page.evaluate(
        "() => metadataColumnValues("
        " metadataDisplayNodeIndices(metadataBaseNodeIndices()), 'family').join('\\n')"
    )
    _open_edit_panel(page, "set")
    page.click("#metadata-paste-open")
    page.select_option("#metadata-paste-column", "copy_target")
    page.fill("#metadata-paste-values", copied)
    page.click("#metadata-paste-apply")

    assert [_row_values(page, i)[2] for i in range(6)] == [
        "alpha", "alpha", "beta", "beta", "gamma", "gamma"
    ]
    assert page.pageerrors == []


def test_edits_and_added_columns_survive_a_session_round_trip(meta_page, tmp_path):
    page = meta_page
    _add_column(page, "annotation")
    _edit_cell(page, 0, "annotation", "kept")
    _edit_cell(page, 1, "family", "edited")

    saved = _save_session(page, tmp_path, "edits.dsnv")
    bundle = _read_session_bundle(saved)
    # The saved file is an ordinary bundle: the edits live in `metadata`, so
    # Python readers see them without knowing anything about app_state.
    assert [c["name"] for c in bundle["metadata"]["columns"]] == ["family", "score", "annotation"]
    assert bundle["metadata"]["columns"][2]["origin"] == "viewer"
    assert bundle["metadata"]["rows"][0] == ["alpha", 1, "kept"]
    assert bundle["metadata"]["rows"][1] == ["edited", 2, None]

    _load_bundle_file(page, saved)
    assert _row_values(page, 0) == ["alpha", 1, "kept"]
    assert _row_values(page, 1) == ["edited", 2, None]
    # `origin` round-trips, so the reloaded column is still known to be viewer-made.
    _open_edit_panel(page, "delete")
    assert page.eval_on_selector_all(
        "#metadata-delete-column option", "opts => opts.map(o => o.textContent)"
    ) == ["family (from bundle)", "score (from bundle)", "annotation (added here)"]
    assert page.pageerrors == []


def test_edited_bundle_is_readable_by_the_python_reader(meta_page, tmp_path):
    """A saved session stays a valid bundle for ssn_bundle/ssn_navigator."""
    from domainator.ssn_bundle import load_bundle, summarize_cluster_metadata

    page = meta_page
    _add_column(page, "cluster_label")
    _select_nodes(page, [0, 1, 2])
    _open_edit_panel(page, "set")
    page.select_option("#metadata-fill-column", "cluster_label")
    page.fill("#metadata-fill-value", "group_one")
    page.click("#metadata-fill-apply")

    bundle = load_bundle(_save_session(page, tmp_path, "annotated.dsnv"))
    summaries = summarize_cluster_metadata([0, 1, 2], bundle["metadata"])
    label_summary = next(s for s in summaries if s["name"] == "cluster_label")
    assert label_summary["count"] == 3
    assert label_summary["top"] == [{"value": "group_one", "count": 3}]
    assert page.pageerrors == []


def test_columns_can_be_added_to_a_bundle_with_no_metadata(page):
    """The editing controls bootstrap a bundle that shipped without metadata."""
    assert page.evaluate("() => state.metadataColumns.length") == 0
    assert page.eval_on_selector("#metadata-fill-apply", "e => e.disabled") is True

    _add_column(page, "notes")

    assert page.eval_on_selector("#metadata-fill-apply", "e => e.disabled") is False
    _select_nodes(page, [0, 1])
    _open_edit_panel(page, "set")
    page.select_option("#metadata-fill-column", "notes")
    page.fill("#metadata-fill-value", "hello")
    page.click("#metadata-fill-apply")
    assert page.evaluate("() => state.metadataByNodeIndex.map(r => r[0])") == [
        "hello", "hello", None, None, None, None
    ]
    assert page.pageerrors == []


def test_added_columns_and_edits_reach_the_tsv_export(meta_page, tmp_path):
    page = meta_page
    _add_column(page, "note")
    _edit_cell(page, 0, "note", "exported")
    _edit_cell(page, 0, "family", "renamed")

    with page.expect_download() as download_info:
        page.click("#export-selected")
    rows = [
        line.split("\t")
        for line in open(download_info.value.path(), encoding="utf-8").read().splitlines()
    ]
    assert rows[0] == ["node_id", "SSN_cluster", "family", "score", "note"]
    assert rows[1][0] == "A" and rows[1][2] == "renamed" and rows[1][4] == "exported"
    assert page.pageerrors == []


def test_session_round_trips_custom_colors(meta_page, tmp_path):
    """Palette edits made in the color picker persist with the session."""
    page = meta_page
    page.evaluate(
        "() => { state.customPalettes[paletteKey('family')] = {colors: {alpha: '#112233'}};"
        " rebuildNodeColorCache(); renderClusterView(); }"
    )
    page.evaluate("() => setColumnCategorical('score', true)")
    assert page.evaluate("() => state.nodeColorCache[0]") == "#112233"

    saved = _save_session(page, tmp_path, "colors.dsnv")
    colors = _read_session_bundle(saved)["app_state"]["colors"]
    assert colors["custom_palettes"]["family"] == {"colors": {"alpha": "#112233"}}
    assert colors["categorical_columns"] == ["score"]

    _load_bundle_file(page, saved)
    assert page.evaluate("() => state.nodeColorCache[0]") == "#112233"
    assert page.evaluate("() => Array.from(state.categoricalColumns)") == ["score"]
    assert page.pageerrors == []


def test_deleting_a_column_drops_its_custom_palette(meta_page):
    """Palettes are keyed by paletteKey(), which suffixes categorical numerics."""
    page = meta_page
    _add_column(page, "grouping", "number")
    page.evaluate(
        "() => { setColumnCategorical('grouping', true);"
        " state.customPalettes[paletteKey('grouping')] = {'1': '#445566'};"
        " state.customPalettes['grouping'] = {'1': '#778899'}; }"
    )
    assert page.evaluate("() => Object.keys(state.customPalettes).length") == 2

    _open_edit_panel(page, "delete")
    page.select_option("#metadata-delete-column", "grouping")
    page.click("#metadata-delete-column-apply")

    assert page.evaluate("() => Object.keys(state.customPalettes)") == []
    assert page.evaluate("() => Array.from(state.categoricalColumns)") == []
    assert page.pageerrors == []


def test_shift_click_adds_a_preset_to_the_selection(page):
    """Shift-click unions the preset into the current selection."""
    _select_nodes(page, [0, 1])
    page.keyboard.press("Shift+Digit1")
    _select_nodes(page, [4, 5])
    page.keyboard.press("Shift+Digit2")

    _select_nodes(page, [0, 1])
    page.click('[data-preset-slot="2"]', modifiers=["Shift"])

    assert page.evaluate("() => Array.from(state.selectedNodeIndices).sort((a,b)=>a-b)") == [0, 1, 4, 5]
    assert "Added preset 2 to the selection (+2 nodes, 4 total)" in page.eval_on_selector(
        "#bundle-status", "e => e.textContent"
    )
    assert page.pageerrors == []


def test_alt_shift_click_subtracts_a_preset_from_the_selection(page):
    _select_nodes(page, [1, 2])
    page.keyboard.press("Shift+Digit3")
    _select_nodes(page, [0, 1, 2, 3])

    page.click('[data-preset-slot="3"]', modifiers=["Shift", "Alt"])

    assert page.evaluate("() => Array.from(state.selectedNodeIndices).sort((a,b)=>a-b)") == [0, 3]
    assert "Subtracted preset 3 from the selection (-2 nodes, 2 remaining)" in page.eval_on_selector(
        "#bundle-status", "e => e.textContent"
    )
    assert page.pageerrors == []


def test_plain_click_still_replaces_the_selection(page):
    """An unmodified click is unchanged: it replaces, it does not union."""
    _select_nodes(page, [0, 1])
    page.keyboard.press("Shift+Digit4")
    _select_nodes(page, [3, 4, 5])

    page.click('[data-preset-slot="4"]')

    assert page.evaluate("() => Array.from(state.selectedNodeIndices).sort((a,b)=>a-b)") == [0, 1]
    assert page.pageerrors == []


def test_adding_a_preset_that_overlaps_does_not_double_count(page):
    _select_nodes(page, [0, 1, 2])
    page.keyboard.press("Shift+Digit5")
    _select_nodes(page, [2, 3])

    page.click('[data-preset-slot="5"]', modifiers=["Shift"])

    assert page.evaluate("() => Array.from(state.selectedNodeIndices).sort((a,b)=>a-b)") == [0, 1, 2, 3]
    assert "+2 nodes, 4 total" in page.eval_on_selector("#bundle-status", "e => e.textContent")
    assert page.pageerrors == []


def test_modified_click_on_an_empty_slot_leaves_the_selection_alone(page):
    _select_nodes(page, [0, 1])
    page.click('[data-preset-slot="6"]', modifiers=["Shift"])
    assert page.evaluate("() => Array.from(state.selectedNodeIndices).sort((a,b)=>a-b)") == [0, 1]
    assert "Preset 6 is empty" in page.eval_on_selector("#bundle-status", "e => e.textContent")
    assert page.pageerrors == []


# ---------------------------------------------------------------------------
# Page title
# ---------------------------------------------------------------------------


def test_heading_names_the_loaded_network(meta_page):
    page = meta_page
    assert page.eval_on_selector("#viewer-title", "e => e.textContent") == (
        "Domainator Similarity Network Viewer: Color Test Viewer"
    )
    # The didactic blurb that used to sit under the heading is gone.
    assert page.query_selector(".hero p") is None
    assert "runs entirely in the browser" not in page.content()
    assert page.pageerrors == []


def test_heading_follows_a_newly_loaded_bundle(meta_page, tmp_path):
    """Loading a different bundle must retitle the page, not keep the build-time name."""
    import gzip
    import json

    page = meta_page
    saved = _save_session(page, tmp_path, "renamed.dsnv")
    bundle = _read_session_bundle(saved)
    bundle["name"] = "Some Other Network"
    renamed = tmp_path / "renamed_bundle.dsnv"
    renamed.write_bytes(gzip.compress(json.dumps(bundle).encode("utf-8")))

    _load_bundle_file(page, renamed)

    assert page.eval_on_selector("#viewer-title", "e => e.textContent") == (
        "Domainator Similarity Network Viewer: Some Other Network"
    )
    assert page.title() == "Some Other Network"
    assert page.pageerrors == []


def _rename_column(page, old_name, new_name):
    _open_edit_panel(page, "rename")
    page.select_option("#metadata-rename-column", old_name)
    page.fill("#metadata-rename-value", new_name)
    page.click("#metadata-rename-apply")


def test_rename_a_metadata_column(meta_page):
    page = meta_page
    _rename_column(page, "family", "clan")

    assert page.evaluate("() => state.metadataColumns.map(c => c.name)") == ["clan", "score"]
    # Values are positional, so the row contents are untouched by a rename.
    assert _row_values(page, 0) == ["alpha", 1]
    assert _cell(page, 0, "clan").inner_text() == "alpha"
    options = page.eval_on_selector_all("#color-by option", "opts => opts.map(o => o.value)")
    assert "clan" in options and "family" not in options
    assert 'Renamed column "family" to "clan"' in page.eval_on_selector(
        "#bundle-status", "e => e.textContent"
    )
    assert page.pageerrors == []


def test_renaming_the_color_by_column_keeps_it_selected(meta_page):
    """A renamed column keeps its role rather than falling back to the default."""
    page = meta_page
    assert page.eval_on_selector("#color-by", "e => e.value") == "family"
    colors_before = page.evaluate("() => state.nodeColorCache.slice()")

    _rename_column(page, "family", "clan")

    assert page.eval_on_selector("#color-by", "e => e.value") == "clan"
    assert page.evaluate("() => state.nodeColorCache.slice()") == colors_before
    assert page.pageerrors == []


def test_rename_carries_column_state_across(meta_page):
    """Width, palette, categorical flag and sort are keyed by name; all follow."""
    page = meta_page
    page.evaluate(
        "() => { setColumnCategorical('score', true);"
        " state.customPalettes[paletteKey('score')] = {colors: {'1': '#010203'}};"
        " state.metadataColumnWidths.set('score', 321);"
        " state.metadataSort = {columnKey: 'score', direction: 'desc'}; }"
    )

    _rename_column(page, "score", "rating")

    assert page.evaluate("() => Array.from(state.categoricalColumns)") == ["rating"]
    # Palettes are stored under paletteKey(), which for a categorical numeric
    # column suffixes the name -- that spelling has to move too.
    assert page.evaluate("() => state.customPalettes[paletteKey('rating')]") == {
        "colors": {"1": "#010203"}
    }
    assert page.evaluate("() => Object.keys(state.customPalettes).length") == 1
    assert page.evaluate(
        "() => Object.keys(state.customPalettes)[0].startsWith('rating')"
    ) is True
    assert page.evaluate("() => state.metadataColumnWidths.get('rating')") == 321
    assert page.evaluate("() => state.metadataSort") == {"columnKey": "rating", "direction": "desc"}
    assert page.pageerrors == []


def test_rename_rejects_duplicate_and_empty_names(meta_page):
    page = meta_page
    _open_edit_panel(page, "rename")
    page.select_option("#metadata-rename-column", "family")

    page.fill("#metadata-rename-value", "score")
    page.click("#metadata-rename-apply")
    assert page.evaluate("() => state.metadataColumns.map(c => c.name)") == ["family", "score"]
    assert 'A column named "score" already exists' in page.eval_on_selector(
        "#bundle-status", "e => e.textContent"
    )

    page.fill("#metadata-rename-value", "   ")
    page.click("#metadata-rename-apply")
    assert page.evaluate("() => state.metadataColumns.map(c => c.name)") == ["family", "score"]
    assert "Enter a column name" in page.eval_on_selector("#bundle-status", "e => e.textContent")

    # node_id is the table's synthetic key column and is reserved.
    page.fill("#metadata-rename-value", "node_id")
    page.click("#metadata-rename-apply")
    assert page.evaluate("() => state.metadataColumns.map(c => c.name)") == ["family", "score"]
    assert page.pageerrors == []


def test_renamed_columns_survive_a_session_round_trip(meta_page, tmp_path):
    page = meta_page
    _rename_column(page, "family", "clan")

    saved = _save_session(page, tmp_path, "renamed_col.dsnv")
    bundle = _read_session_bundle(saved)
    assert [c["name"] for c in bundle["metadata"]["columns"]] == ["clan", "score"]
    assert bundle["metadata"]["rows"][0] == ["alpha", 1]

    _load_bundle_file(page, saved)
    assert page.evaluate("() => state.metadataColumns.map(c => c.name)") == ["clan", "score"]
    assert page.eval_on_selector("#color-by", "e => e.value") == "clan"
    assert page.pageerrors == []


def test_rename_panel_joins_the_accordion(meta_page):
    page = meta_page
    page.click("#metadata-panel-rename")
    assert page.eval_on_selector("#metadata-rename-panel", "e => e.hidden") is False
    page.click("#metadata-panel-delete")
    assert page.eval_on_selector("#metadata-rename-panel", "e => e.hidden") is True
    assert page.eval_on_selector("#metadata-delete-panel", "e => e.hidden") is False
    assert page.pageerrors == []


def test_rename_retargets_every_column_menu(meta_page):
    """Which column each menu points at lives in the DOM, so it is one more
    thing keyed by column name that a rename has to carry across."""
    page = meta_page
    page.select_option("#label-by", "score")
    _open_edit_panel(page, "set")
    page.select_option("#metadata-fill-column", "score")
    _open_edit_panel(page, "delete")
    page.select_option("#metadata-delete-column", "score")

    _rename_column(page, "score", "rating")

    def value(selector):
        return page.eval_on_selector(selector, "e => e.value")

    assert value("#label-by") == "rating"
    assert value("#metadata-fill-column") == "rating"
    assert value("#metadata-delete-column") == "rating"
    assert value("#metadata-rename-column") == "rating"
    # Menus that were pointing elsewhere are left alone.
    assert value("#color-by") == "family"
    assert value("#metadata-paste-column") == "family"
    assert page.pageerrors == []


# ---------------------------------------------------------------------------
# Cluster-number column (the viewer's `build_ssn.py --lb N --cluster`)
# ---------------------------------------------------------------------------


def _set_threshold(page, threshold_value):
    """Snap the slider to the stop nearest `threshold_value` and redraw."""
    page.evaluate(
        "v => { snapSliderToStop(nearestStopForThreshold(v)); updateThresholdUI(false); }",
        threshold_value,
    )


def test_cluster_column_matches_build_ssn(meta_page):
    """The column the button writes is the one `build_ssn.py --cluster` writes."""
    from domainator.build_ssn import cluster_labels_from_tree
    from domainator.data_matrix import DenseDataMatrix, MaxTree

    page = meta_page
    _set_threshold(page, 6.0)
    threshold = page.evaluate("() => selectedThresholdValue()")

    _open_edit_panel(page, "add")
    page.click("#metadata-add-cluster-column")

    assert page.evaluate("() => state.metadataColumns.map(c => c.name)") == [
        "family", "score", "SSN_cluster"
    ]
    viewer_numbers = page.evaluate("() => state.metadataByNodeIndex.map(r => r[2])")

    # Same matrix the meta_viewer_html fixture was built from.
    data = np.array([
        [0, 10, 6, 0, 0, 0],
        [10, 0, 7, 0, 0, 0],
        [6, 7, 0, 4, 0, 0],
        [0, 0, 4, 0, 8, 5],
        [0, 0, 0, 8, 0, 9],
        [0, 0, 0, 5, 9, 0],
    ], dtype=float)
    names = ["A", "B", "C", "D", "E", "F"]
    expected = cluster_labels_from_tree(MaxTree(DenseDataMatrix(data, names, names)), threshold)

    assert viewer_numbers == [int(v) for v in expected]
    assert page.pageerrors == []


def test_cluster_column_tracks_the_threshold(meta_page):
    """Re-running at a finer threshold rewrites the column with more clusters."""
    page = meta_page
    _open_edit_panel(page, "add")

    # The lowest stop is the weakest MST edge, i.e. the coarsest partition.
    _set_threshold(page, 0.0)
    page.click("#metadata-add-cluster-column")
    coarse = page.evaluate("() => state.metadataByNodeIndex.map(r => r[2])")

    _set_threshold(page, 8.5)
    page.click("#metadata-add-cluster-column")
    fine = page.evaluate("() => state.metadataByNodeIndex.map(r => r[2])")

    # Raising the threshold cuts more edges, so the partition only gets finer.
    assert len(set(fine)) > len(set(coarse))
    assert set(coarse) == set(range(1, len(set(coarse)) + 1))
    assert set(fine) == set(range(1, len(set(fine)) + 1))
    # Overwritten in place rather than added a second time.
    assert page.evaluate("() => state.metadataColumns.map(c => c.name)") == [
        "family", "score", "SSN_cluster"
    ]
    assert page.pageerrors == []


def test_cluster_numbers_are_ranked_by_size(meta_page):
    """Largest cluster is 1, as in build_ssn.py's size-rank numbering."""
    page = meta_page
    _open_edit_panel(page, "add")
    _set_threshold(page, 6.0)
    page.click("#metadata-add-cluster-column")

    numbers = page.evaluate("() => state.metadataByNodeIndex.map(r => r[2])")
    sizes = {}
    for number in numbers:
        sizes[number] = sizes.get(number, 0) + 1
    ranked = [sizes[n] for n in sorted(sizes)]
    assert ranked == sorted(ranked, reverse=True)
    assert min(sizes) == 1
    assert page.pageerrors == []


def test_cluster_column_takes_a_custom_name_and_is_categorical(meta_page):
    page = meta_page
    _open_edit_panel(page, "add")
    page.fill("#metadata-new-column-name", "clusters_at_6")
    page.click("#metadata-add-cluster-column")

    assert page.evaluate("() => state.metadataColumns.map(c => c.name)") == [
        "family", "score", "clusters_at_6"
    ]
    assert page.evaluate("() => metadataColumnType('clusters_at_6')") == "int"
    # Cluster ids are labels, not magnitudes, so they color as discrete categories.
    assert "clusters_at_6" in page.evaluate("() => Array.from(state.categoricalColumns)")
    assert page.eval_on_selector("#metadata-new-column-name", "e => e.value") == ""
    # A second name gives a second snapshot rather than replacing the first.
    page.fill("#metadata-new-column-name", "clusters_at_9")
    page.click("#metadata-add-cluster-column")
    assert page.evaluate("() => state.metadataColumns.length") == 4
    assert page.pageerrors == []


def test_cluster_column_confirms_before_overwriting_bundle_data(meta_page):
    page = meta_page
    _open_edit_panel(page, "add")
    page.fill("#metadata-new-column-name", "family")

    page.once("dialog", lambda dialog: dialog.dismiss())
    page.click("#metadata-add-cluster-column")
    assert _row_values(page, 0) == ["alpha", 1]

    page.fill("#metadata-new-column-name", "family")
    page.once("dialog", lambda dialog: dialog.accept())
    page.click("#metadata-add-cluster-column")
    assert isinstance(_row_values(page, 0)[0], int)
    assert page.pageerrors == []


def test_cluster_column_agrees_with_the_tsv_export(meta_page, tmp_path):
    """Both read the same helper, so the column and the export cannot drift."""
    page = meta_page
    _set_threshold(page, 6.0)
    _open_edit_panel(page, "add")
    page.click("#metadata-add-cluster-column")

    with page.expect_download() as download_info:
        page.click("#export-selected")
    rows = [
        line.split("\t")
        for line in open(download_info.value.path(), encoding="utf-8").read().splitlines()
    ]
    header, body = rows[0], rows[1:]
    exported = {row[0]: row[header.index("SSN_cluster")] for row in body}

    in_column = page.evaluate(
        "() => Object.fromEntries(state.bundle.graph.nodes.map("
        "(id, i) => [id, String(metadataValue(i, 'SSN_cluster'))]))"
    )
    assert exported == in_column
    assert page.pageerrors == []


# ---------------------------------------------------------------------------
# Save extraction
# ---------------------------------------------------------------------------


def _save_extraction(page, out_dir, filename="extraction.dsnv", name=None):
    """Click "Save extraction…", name it in the dialog, and save the download."""
    page.click("#save-extraction")
    page.wait_for_selector("#name-overlay:not([hidden])")
    if name is not None:
        page.fill("#name-value", name)
    with page.expect_download() as download_info:
        page.click("#name-apply")
    target = out_dir / filename
    download_info.value.save_as(str(target))
    return target


def _assert_js_series_matches_python(page, graph):
    """The viewer derives its split-event series in JavaScript; ssn_hierarchy derives
    the same thing in Python. Same bundle, same numbers -- that is the port contract,
    and from bundle v6 it is the only thing keeping the two in step, since neither
    reads a stored copy any more.

    The cap and the chart's window are read off the page, so this holds zoomed as well
    as flat.
    """
    from domainator import ssn_bundle
    from domainator.ssn_hierarchy import filter_merge_event_rows, threshold_slider_stops

    bundle = {"graph": graph}
    max_merge_events = page.evaluate("() => state.maxMergeEvents")
    window = page.evaluate("() => currentSplitWindow()")
    window = None if window is None else (window["min"], window["max"])
    pinned = page.evaluate("() => selectedThresholdValue()")
    pinned = pinned if isinstance(pinned, (int, float)) else None
    py_events = ssn_bundle.merge_event_rows(bundle)
    js_events = page.evaluate("() => state.series.eventRows")

    assert len(js_events) == len(py_events)
    for got, want in zip(js_events, py_events):
        for field in ("edge_index", "threshold_to", "threshold_from", "merge_count",
                      "summary_row_from", "summary_row_to"):
            assert got[field] == want[field], field
        # JSON object keys are strings; Python keys the histogram by int.
        assert ({str(k): v for k, v in got["merge_size_counts"].items()}
                == {str(k): v for k, v in want["merge_size_counts"].items()})
        for field in ("threshold_value", "merge_impact", "largest_merge",
                      "delta_largest", "delta_avg_non_singleton"):
            assert got[field] == pytest.approx(want[field]), field

    py_selected = filter_merge_event_rows(
        py_events, max_merge_events=max_merge_events, window=window, pinned_threshold=pinned)
    js_selected = page.evaluate("() => state.series.selectedRows")
    assert ([row["edge_index"] for row in js_selected]
            == [row["edge_index"] for row in py_selected])

    py_stops = threshold_slider_stops(py_selected, tree=ssn_bundle.bundle_tree(bundle))
    js_stops = page.evaluate("() => state.series.sliderStops")
    assert len(js_stops) == len(py_stops)
    for got, want in zip(js_stops, py_stops):
        assert got["edge_index"] == want["edge_index"]
        assert got["threshold_label"] == want["threshold_label"]
        # threshold_index went with the threshold tables it indexed.
        assert "threshold_index" not in got
        if want["threshold_value"] is None:
            assert got["threshold_value"] is None
        else:
            assert got["threshold_value"] == pytest.approx(want["threshold_value"])

    py_sum = ssn_bundle.moving_sum(bundle, event_rows=py_events)
    js_sum = page.evaluate("() => state.series.movingSum")
    assert js_sum["window"] == pytest.approx(py_sum["window"])
    assert js_sum["x"] == pytest.approx(py_sum["x"], rel=1e-9)
    assert js_sum["y"] == pytest.approx(py_sum["y"], rel=1e-9)


def test_extraction_of_everything_reproduces_the_python_bundle(meta_page, tmp_path):
    """The JS hierarchy/merge-series port must agree with ssn_hierarchy.py.

    Selecting every node makes the extraction a rebuild of the whole network, so its
    structure can be compared field-for-field against what build_ssn_viewer.py wrote.
    The derived series is checked against Python separately, since from v6 neither the
    original nor the extraction stores one.
    """
    page = meta_page
    original = page.evaluate("() => state.bundle.graph")
    _assert_js_series_matches_python(page, original)
    _select_nodes(page, list(range(6)))

    extracted = _read_session_bundle(_save_extraction(page, tmp_path))["graph"]

    assert extracted["nodes"] == original["nodes"]
    assert extracted["mst_edges"] == original["mst_edges"]
    assert extracted["hierarchy"] == original["hierarchy"]
    assert extracted["merge_impact_metric"] == original["merge_impact_metric"]
    # A v6 graph carries nothing derived, so an extraction writes nothing derived.
    assert set(extracted) == {"nodes", "mst_edges", "hierarchy", "merge_impact_metric"}
    assert page.pageerrors == []


def test_extraction_refuses_a_selection_split_across_the_mst(meta_page, tmp_path):
    """Two pieces with no MST path between them cannot be extracted faithfully."""
    page = meta_page
    # A and F sit at opposite ends of the 6-node path, joined only through B-E.
    _select_nodes(page, [0, 5])

    messages = []
    page.on("dialog", lambda dialog: (messages.append(dialog.message), dialog.accept()))
    page.click("#save-extraction")
    page.wait_for_function(
        "() => document.getElementById('bundle-status').textContent.includes('pieces')"
    )

    assert len(messages) == 1
    assert "splits 1 network component into 2 pieces" in messages[0]
    assert "leaving out the nodes that join them in the MST" in messages[0]
    assert "--subset" in messages[0]
    assert page.pageerrors == []


def test_extraction_writes_a_loadable_subset(meta_page, tmp_path):
    """A connected selection extracts to a bundle the viewer can open."""
    from domainator.ssn_bundle import clusters_at_threshold, coarsest_threshold, load_bundle

    page = meta_page
    _select_nodes(page, [0, 1, 2])          # A-B-C, connected in the MST
    saved = _save_extraction(page, tmp_path, "abc.dsnv")

    bundle = load_bundle(saved)             # valid for the Python readers too
    assert bundle["graph"]["nodes"] == ["A", "B", "C"]
    assert bundle["metadata"]["rows"] == [["alpha", 1], ["alpha", 2], ["beta", 5]]
    assert bundle["name"] == "Color Test Viewer_extraction"
    # Connected, so the hierarchy has a single root spanning all three nodes.
    assert len(bundle["graph"]["hierarchy"]["roots"]) == 1
    assert sorted(bundle["graph"]["hierarchy"]["leaf_order"]) == [0, 1, 2]
    hierarchy = bundle["graph"]["hierarchy"]
    assert clusters_at_threshold(hierarchy, coarsest_threshold(hierarchy)) == hierarchy["roots"]

    _load_bundle_file(page, saved)
    assert page.evaluate("() => state.bundle.graph.nodes") == ["A", "B", "C"]
    assert page.eval_on_selector("#viewer-title", "e => e.textContent").endswith("_extraction")
    assert page.pageerrors == []


def test_extraction_reindexes_mst_edges(meta_page, tmp_path):
    """Edge endpoints are renumbered into the extraction's own node list."""
    page = meta_page
    _select_nodes(page, [3, 4, 5])          # D-E-F
    bundle = _read_session_bundle(_save_extraction(page, tmp_path, "def.dsnv"))

    assert bundle["graph"]["nodes"] == ["D", "E", "F"]
    for source, target, score in bundle["graph"]["mst_edges"]:
        assert 0 <= source < 3 and 0 <= target < 3
        assert score > 0
    assert page.pageerrors == []


def test_extraction_carries_session_state_for_surviving_nodes_only(meta_page, tmp_path):
    """Presets and the selection are filtered to what the extraction contains."""
    page = meta_page
    _select_nodes(page, [4, 5])
    page.keyboard.press("Shift+Digit1")     # a preset that survives
    _select_nodes(page, [0, 1])
    page.keyboard.press("Shift+Digit2")     # a preset that does not
    page.select_option("#color-by", "score")

    _select_nodes(page, [3, 4, 5])
    app_state = _read_session_bundle(_save_extraction(page, tmp_path, "state.dsnv"))["app_state"]

    assert app_state["view"]["color_by"] == "score"
    assert app_state["selection"]["node_ids"] == ["D", "E", "F"]
    assert app_state["selection"]["presets"]["1"]["node_ids"] == ["E", "F"]
    # Slot 2 held only A and B, so it is dropped rather than saved empty.
    assert "2" not in app_state["selection"]["presets"]
    assert page.pageerrors == []


def test_extraction_loads_without_reporting_missing_node_ids(meta_page, tmp_path):
    """Filtering the saved ids means the extraction opens cleanly."""
    page = meta_page
    _select_nodes(page, [0, 1, 2])
    page.keyboard.press("Shift+Digit4")
    saved = _save_extraction(page, tmp_path, "clean.dsnv")

    _load_bundle_file(page, saved)
    status = page.eval_on_selector("#bundle-status", "e => e.textContent")
    assert "not in this bundle" not in status
    assert "skipped" not in status
    assert page.evaluate("() => Array.from(state.selectedNodeIndices).sort()") == [0, 1, 2]
    assert page.evaluate("() => state.selectionPresets.get(4).nodeIndices.size") == 3
    assert page.pageerrors == []


def test_save_extraction_button_needs_a_selection(meta_page):
    page = meta_page
    assert page.eval_on_selector("#save-extraction", "e => e.disabled") is True
    _select_nodes(page, [0, 1])
    assert page.eval_on_selector("#save-extraction", "e => e.disabled") is False
    page.click("#clear-selection")
    assert page.eval_on_selector("#save-extraction", "e => e.disabled") is True
    assert page.pageerrors == []


def test_extraction_of_everything_matches_python_on_a_tied_network(many_cat_page, tmp_path):
    """The 120-node fixture is a path of equal-weight edges, so every MST edge
    lands in one tie group -- the case the per-threshold grouping exists for."""
    page = many_cat_page
    original = page.evaluate("() => state.bundle.graph")
    _select_nodes(page, list(range(120)))

    extracted = _read_session_bundle(_save_extraction(page, tmp_path, "all120.dsnv"))["graph"]

    assert extracted["hierarchy"] == original["hierarchy"]
    assert extracted["mst_edges"] == original["mst_edges"]
    _assert_js_series_matches_python(page, original)
    # One tie group covering all 119 merges.
    events = page.evaluate("() => state.series.eventRows")
    assert len(events) == 1
    assert events[0]["merge_count"] == 119
    assert page.pageerrors == []


def test_extraction_of_a_contiguous_run_is_connected(many_cat_page, tmp_path):
    """A contiguous stretch of the path is MST-connected; a gapped one is not."""
    page = many_cat_page
    _select_nodes(page, list(range(10, 30)))
    bundle = _read_session_bundle(_save_extraction(page, tmp_path, "run.dsnv"))
    assert len(bundle["graph"]["nodes"]) == 20
    assert len(bundle["graph"]["mst_edges"]) == 19
    assert len(bundle["graph"]["hierarchy"]["roots"]) == 1

    _select_nodes(page, list(range(10, 20)) + list(range(40, 50)))
    messages = []
    page.on("dialog", lambda dialog: (messages.append(dialog.message), dialog.accept()))
    page.click("#save-extraction")
    page.wait_for_function(
        "() => document.getElementById('bundle-status').textContent.includes('pieces')"
    )
    assert "splits 1 network component into 2 pieces" in messages[0]
    assert page.pageerrors == []


def test_extraction_of_everything_matches_python_on_distinct_weights(dense_page, tmp_path):
    """Same equivalence check where nearly every MST edge has its own weight, so
    the merge series has many rows and the moving sum window actually slides."""
    page = dense_page
    original = page.evaluate("() => state.bundle.graph")
    # Many distinct thresholds, so the moving-sum window actually slides.
    assert page.evaluate("() => state.series.eventRows.length") > 40
    _assert_js_series_matches_python(page, original)
    _select_nodes(page, list(range(60)))

    extracted = _read_session_bundle(_save_extraction(page, tmp_path, "all60.dsnv"))["graph"]

    assert extracted["hierarchy"] == original["hierarchy"]
    assert extracted["mst_edges"] == original["mst_edges"]
    # Reloading the extraction must derive the very same series the original did.
    _load_bundle_file(page, tmp_path / "all60.dsnv")
    _assert_js_series_matches_python(page, extracted)
    assert page.pageerrors == []


def test_extraction_allows_separate_network_components(dense_page, tmp_path):
    """Components that were already unrelated stay apart without complaint.

    The check is per original component, not global: keeping two components
    apart invents nothing, so a whole multi-component network extracts cleanly.
    """
    page = dense_page
    roots = page.evaluate("() => state.bundle.graph.hierarchy.roots.length")
    assert roots == 3          # the fixture is three unconnected blocks

    _select_nodes(page, list(range(60)))
    bundle = _read_session_bundle(_save_extraction(page, tmp_path, "multi.dsnv"))
    assert len(bundle["graph"]["hierarchy"]["roots"]) == 3
    assert len(bundle["graph"]["nodes"]) == 60
    assert page.pageerrors == []


def test_extraction_still_refuses_a_component_broken_in_two(dense_page, tmp_path):
    """Exempting separate components must not exempt a component split by the
    selection -- that is the case the MST cannot speak to."""
    page = dense_page
    # Two halves of the first block with the linking nodes left out.
    _select_nodes(page, list(range(0, 5)) + list(range(15, 22)))
    messages = []
    page.on("dialog", lambda dialog: (messages.append(dialog.message), dialog.accept()))
    page.click("#save-extraction")
    page.wait_for_function(
        "() => document.getElementById('bundle-status').textContent.includes('pieces')"
    )
    assert "network component" in messages[0]
    assert page.pageerrors == []


def test_leftmost_slider_stop_shows_the_true_components(dense_page):
    """The lowest stop must merge everything the MST joins.

    Regression: every stop excludes its own tie group, so the lowest one used to
    be the weakest MST edge -- which that rule then split back apart. The viewer
    showed one cluster too many and could never reach the real components.
    """
    page = dense_page
    roots = page.evaluate("() => state.bundle.graph.hierarchy.roots.length")
    assert roots == 3

    page.evaluate(
        "() => { document.getElementById('threshold-slider').value = '0';"
        " snapSliderToStop(currentSliderStop()); updateThresholdUI(false); }"
    )
    page.wait_for_function(
        "n => document.getElementById('stat-clusters').textContent === String(n)", arg=roots
    )

    stop = page.evaluate("() => currentSliderStop()")
    weights = page.evaluate("() => state.bundle.graph.mst_edges.map(e => e[2])")
    weakest, span = min(weights), max(weights) - min(weights)
    # Strictly below the weakest merge, but only by 1% of the weight range: a floor
    # at 0 would leave most of the slider track empty on a network scoring 350-650.
    assert stop["threshold_value"] < weakest
    assert stop["threshold_value"] == pytest.approx(weakest - 0.01 * span)
    assert page.evaluate("() => activeClustersAtThreshold(selectedThresholdValue()).length") == roots
    assert page.pageerrors == []


def test_extraction_slider_keeps_the_floor_stop(meta_page, tmp_path):
    """The extraction port mirrors the floor stop, so extractions are reachable
    at their own fully merged view too."""
    page = meta_page
    _select_nodes(page, [0, 1, 2])
    saved = _save_extraction(page, tmp_path, "floor.dsnv")

    graph = _read_session_bundle(saved)["graph"]
    from domainator import ssn_bundle
    stops = ssn_bundle.slider_stops({"graph": graph})
    weights = [edge[2] for edge in graph["mst_edges"]]
    weakest, span = min(weights), max(weights) - min(weights)
    # The JS port applies the same 1%-of-range offset as ssn_hierarchy does.
    assert stops[-1]["threshold_value"] < weakest
    assert stops[-1]["threshold_value"] == pytest.approx(weakest - 0.01 * span)
    assert stops[-1]["edge_index"] == 1          # both induced MST edges are above it

    _load_bundle_file(page, saved)
    page.evaluate(
        "() => { document.getElementById('threshold-slider').value = '0';"
        " snapSliderToStop(currentSliderStop()); updateThresholdUI(false); }"
    )
    page.wait_for_function("() => document.getElementById('stat-clusters').textContent === '1'")
    assert page.pageerrors == []


def test_session_restores_categorical_colouring_of_a_numeric_column(meta_page, tmp_path):
    """A numeric column marked categorical must come back *colored* that way.

    Regression: whether a numeric column colors as discrete categories is derived
    from state.categoricalColumns by rebuildMetadataCaches, which runs before a
    session is applied. Restoring the set alone left the derived flag stale, so the
    checkbox read categorical while the colors and the picker stayed on the
    gradient until the box was toggled by hand.
    """
    page = meta_page
    page.select_option("#color-by", "score")
    page.evaluate(
        "() => { setColumnCategorical('score', true);"
        " rebuildNodeColorCache(); renderClusterView(); }"
    )
    assert page.evaluate("() => colorInfo('score').type") == "categorical"
    categorical_colors = page.evaluate("() => state.nodeColorCache.slice()")

    saved = _save_session(page, tmp_path, "categorical.dsnv")
    _load_bundle_file(page, saved)

    assert page.evaluate("() => state.categoricalColumns.has('score')") is True
    assert page.evaluate("() => colorInfo('score').type") == "categorical"
    assert page.evaluate("() => state.nodeColorCache.slice()") == categorical_colors

    _open_color_picker(page)
    assert page.eval_on_selector("#color-as-categorical", "e => e.checked") is True
    assert page.eval_on_selector("#color-picker-discrete", "e => e.hidden") is False
    assert page.eval_on_selector("#color-picker-continuous", "e => e.hidden") is True
    assert page.pageerrors == []


def test_session_restores_gradient_colouring_when_the_flag_was_cleared(numeric_categorical_page, tmp_path):
    """The same derivation has to run the other way too.

    This bundle ships --categorical score, so a session that turned it *off* must
    come back on the gradient rather than inheriting the bundle's default.
    """
    page = numeric_categorical_page
    assert page.evaluate("() => colorInfo('score').type") == "categorical"
    page.evaluate(
        "() => { setColumnCategorical('score', false);"
        " rebuildNodeColorCache(); renderClusterView(); }"
    )
    gradient_colors = page.evaluate("() => state.nodeColorCache.slice()")

    saved = _save_session(page, tmp_path, "gradient.dsnv")
    _load_bundle_file(page, saved)

    assert page.evaluate("() => state.categoricalColumns.has('score')") is False
    assert page.evaluate("() => colorInfo('score').type") == "numeric"
    assert page.evaluate("() => state.nodeColorCache.slice()") == gradient_colors

    _open_color_picker(page)
    assert page.eval_on_selector("#color-as-categorical", "e => e.checked") is False
    assert page.eval_on_selector("#color-picker-continuous", "e => e.hidden") is False
    assert page.pageerrors == []


# ---------------------------------------------------------------------------
# Naming
# ---------------------------------------------------------------------------


def test_rename_network_updates_heading_tab_and_bundle(meta_page):
    page = meta_page
    page.click("#rename-network")
    page.wait_for_selector("#name-overlay:not([hidden])")
    assert page.eval_on_selector("#name-value", "e => e.value") == "Color Test Viewer"

    page.fill("#name-value", "GH17 family")
    page.click("#name-apply")

    assert page.eval_on_selector("#name-overlay", "e => e.hidden") is True
    assert page.evaluate("() => state.bundle.name") == "GH17 family"
    assert page.eval_on_selector("#viewer-title", "e => e.textContent") == (
        "Domainator Similarity Network Viewer: GH17 family"
    )
    assert page.title() == "GH17 family"
    assert page.pageerrors == []


def test_double_clicking_the_heading_opens_the_rename_dialog(meta_page):
    page = meta_page
    page.dblclick("#viewer-title")
    page.wait_for_selector("#name-overlay:not([hidden])")
    assert page.eval_on_selector("#name-dialog-title", "e => e.textContent") == "Rename network"
    assert page.eval_on_selector("#name-value", "e => e.value") == "Color Test Viewer"
    assert page.pageerrors == []


def test_rename_is_carried_into_saved_files(meta_page, tmp_path):
    """The name is stored in the bundle and names the file it is saved to."""
    page = meta_page
    page.click("#rename-network")
    page.fill("#name-value", "renamed net")
    page.click("#name-apply")

    with page.expect_download() as download_info:
        page.click("#save-session")
    download = download_info.value
    assert download.suggested_filename == "renamed_net_session.dsnv"
    saved = tmp_path / "renamed.dsnv"
    download.save_as(str(saved))
    assert _read_session_bundle(saved)["name"] == "renamed net"

    _load_bundle_file(page, saved)
    assert page.eval_on_selector("#viewer-title", "e => e.textContent").endswith("renamed net")
    assert page.pageerrors == []


def test_rename_dialog_rejects_an_empty_name(meta_page):
    page = meta_page
    page.click("#rename-network")
    page.fill("#name-value", "   ")
    page.click("#name-apply")

    # The dialog stays open so the name can be corrected.
    assert page.eval_on_selector("#name-overlay", "e => e.hidden") is False
    assert page.evaluate("() => state.bundle.name") == "Color Test Viewer"
    assert "Enter a name" in page.eval_on_selector("#bundle-status", "e => e.textContent")
    page.click("#name-cancel")
    assert page.pageerrors == []


def test_extraction_is_named_through_the_dialog(meta_page, tmp_path):
    page = meta_page
    _select_nodes(page, [0, 1, 2])
    page.click("#save-extraction")
    page.wait_for_selector("#name-overlay:not([hidden])")

    assert page.eval_on_selector("#name-dialog-title", "e => e.textContent") == "Save extraction"
    # Pre-filled with a sensible default, and says what is being extracted.
    assert page.eval_on_selector("#name-value", "e => e.value") == "Color Test Viewer_extraction"
    assert "3 of 6 nodes" in page.eval_on_selector("#name-dialog-note", "e => e.textContent")

    page.fill("#name-value", "clade A")
    with page.expect_download() as download_info:
        page.click("#name-apply")
    download = download_info.value
    assert download.suggested_filename == "clade_A.dsnv"
    saved = tmp_path / "cladeA.dsnv"
    download.save_as(str(saved))

    bundle = _read_session_bundle(saved)
    assert bundle["name"] == "clade A"
    assert bundle["graph"]["nodes"] == ["A", "B", "C"]
    # The source network keeps its own name.
    assert page.evaluate("() => state.bundle.name") == "Color Test Viewer"

    _load_bundle_file(page, saved)
    assert page.eval_on_selector("#viewer-title", "e => e.textContent").endswith("clade A")
    assert page.pageerrors == []


def test_cancelling_the_extraction_dialog_saves_nothing(meta_page):
    page = meta_page
    _select_nodes(page, [0, 1, 2])
    page.click("#save-extraction")
    page.wait_for_selector("#name-overlay:not([hidden])")

    downloads = []
    page.on("download", lambda download: downloads.append(download))
    page.click("#name-cancel")

    assert page.eval_on_selector("#name-overlay", "e => e.hidden") is True
    assert downloads == []
    assert page.pageerrors == []


def test_a_refused_extraction_never_asks_for_a_name(meta_page):
    """Validation runs first, so an impossible extraction is not named first."""
    page = meta_page
    _select_nodes(page, [0, 5])          # two pieces, no MST path between them
    page.once("dialog", lambda dialog: dialog.accept())
    page.click("#save-extraction")
    page.wait_for_function(
        "() => document.getElementById('bundle-status').textContent.includes('pieces')"
    )
    assert page.eval_on_selector("#name-overlay", "e => e.hidden") is True
    assert page.pageerrors == []


def test_digit_shortcuts_are_inert_while_naming(meta_page):
    """Typing a digit into the name field must not recall a preset."""
    page = meta_page
    _select_nodes(page, [0, 1])
    page.keyboard.press("Shift+Digit3")
    page.click("#clear-selection")

    page.click("#rename-network")
    page.fill("#name-value", "")
    page.keyboard.press("Digit3")

    assert page.eval_on_selector("#name-value", "e => e.value") == "3"
    assert page.evaluate("() => state.selectedNodeIndices.size") == 0
    page.click("#name-cancel")
    assert page.pageerrors == []


def _selection_screen_box(page):
    """Screen-space box of the selection, plus the scale and canvas size."""
    return page.evaluate(
        """() => {
            const bounds = selectedNodeBounds();
            if (!bounds) { return null; }
            const view = state.viewTransform;
            const canvas = document.getElementById('cluster-view');
            return {
                left: (bounds.minX * view.scale) + view.offsetX,
                top: (bounds.minY * view.scale) + view.offsetY,
                right: (bounds.maxX * view.scale) + view.offsetX,
                bottom: (bounds.maxY * view.scale) + view.offsetY,
                scale: view.scale,
                width: canvas.width,
                height: canvas.height,
            };
        }"""
    )


def test_focus_selection_disabled_without_a_selection(meta_page):
    """The button is selection-driven, so it stays off until something is picked."""
    page = meta_page
    assert page.eval_on_selector("#focus-selection", "e => e.disabled") is True
    _select_nodes(page, [0, 1, 2])
    assert page.eval_on_selector("#focus-selection", "e => e.disabled") is False
    page.click("#clear-selection")
    assert page.eval_on_selector("#focus-selection", "e => e.disabled") is True
    assert page.pageerrors == []


def test_focus_selection_centers_without_changing_zoom_when_it_fits(meta_page):
    """A selection that already fits is centered at the zoom the user set."""
    page = meta_page
    _select_nodes(page, [0, 1, 2])
    before_scale = page.evaluate("() => state.viewTransform.scale")

    page.click("#focus-selection")
    box = _selection_screen_box(page)

    assert box["scale"] == before_scale
    assert abs(((box["left"] + box["right"]) / 2) - (box["width"] / 2)) < 0.5
    assert abs(((box["top"] + box["bottom"]) / 2) - (box["height"] / 2)) < 0.5
    assert page.pageerrors == []


def test_focus_selection_zooms_out_until_everything_fits(meta_page):
    """Zoomed far in on the whole network, focusing must pull back until it fits."""
    page = meta_page
    # A small canvas makes the fit-to-width/height calculation the binding
    # constraint rather than viewTransform.maxScale, which is what we want to test.
    page.evaluate(
        """() => {
            const canvas = document.getElementById('cluster-view');
            canvas.width = 320;
            canvas.height = 240;
            state.viewTransform.scale = 30;
        }"""
    )
    _select_nodes(page, list(range(60)))

    page.click("#focus-selection")
    box = _selection_screen_box(page)

    assert box["scale"] < 30
    assert box["left"] >= 0 and box["top"] >= 0
    assert box["right"] <= box["width"] and box["bottom"] <= box["height"]
    assert abs(((box["left"] + box["right"]) / 2) - (box["width"] / 2)) < 0.5
    assert abs(((box["top"] + box["bottom"]) / 2) - (box["height"] / 2)) < 0.5
    assert page.pageerrors == []


def test_focus_selection_never_zooms_in(meta_page):
    """Zoomed out, a tiny selection is centered but the view is not magnified."""
    page = meta_page
    page.evaluate("() => { state.viewTransform.scale = 0.5; }")
    _select_nodes(page, [0])

    page.click("#focus-selection")

    assert page.evaluate("() => state.viewTransform.scale") == 0.5
    assert page.pageerrors == []


def test_focus_selection_reports_when_nothing_selected_is_visible(meta_page):
    """Selected nodes in clusters the current view hides have no position."""
    page = meta_page
    _select_nodes(page, [0, 1])
    page.evaluate("() => { state.visibleLayout = []; }")

    page.click("#focus-selection")

    status = page.eval_on_selector("#bundle-status", "e => e.textContent")
    assert "No selected nodes are visible" in status
    assert page.pageerrors == []


def test_focus_selection_accounts_for_node_radius(meta_page):
    """Bounds grow by each node's radius, so the drawn dot fits, not just its center."""
    page = meta_page
    _select_nodes(page, [0])
    measured = page.evaluate(
        """() => {
            const bounds = selectedNodeBounds();
            const item = state.visibleLayout.find(
                entry => componentMembers(entry.componentId).includes(0));
            const component = state.bundle.graph.hierarchy.nodes[item.componentId];
            const dot = componentMemberLayout(component, item)
                .find(entry => entry.memberIndex === 0);
            return {span: bounds.maxX - bounds.minX, radius: dot.radius};
        }"""
    )
    assert measured["radius"] > 0
    assert abs(measured["span"] - (measured["radius"] * 2)) < 1e-6
    assert page.pageerrors == []


def _open_gradient_picker(page, column="score"):
    """Open the color picker on a numeric column, in gradient (not discrete) mode."""
    page.select_option("#color-by", column)
    _open_color_picker(page)
    page.wait_for_selector("#color-picker-continuous:not([hidden])")


def _histogram(page):
    return page.evaluate("() => state.colorHistogram")


def test_gradient_histogram_bins_integers_one_per_value(meta_page):
    """An int column gets one bin per integer, not evenly divided bins.

    The fixture's scores are 1, 2, 5, 8, 3, 9. Splitting that range into e.g. 8
    even bins leaves empty gaps and doubled-up bars that read as structure in the
    data when they are only an artifact of the binning.
    """
    page = meta_page
    _open_gradient_picker(page)

    histogram = _histogram(page)
    assert histogram["binCount"] == 9
    # The axis is exactly the data range, so the knobs below can reach both ends.
    assert [histogram["lowEdge"], histogram["highEdge"]] == [1, 9]
    #                  1  2  3  4  5  6  7  8  9
    assert histogram["counts"] == [1, 1, 1, 0, 1, 0, 0, 1, 1]
    assert histogram["total"] == 6
    assert histogram["missing"] == 0
    assert page.is_visible("#color-histogram-wrap")
    assert page.eval_on_selector("#color-histogram-note", "e => e.textContent") == (
        "6 values in 9 bins · peak bin 1"
    )
    assert page.pageerrors == []


def test_gradient_histogram_paints_bars(meta_page):
    """The canvas is actually drawn on, not left blank."""
    page = meta_page
    _open_gradient_picker(page)
    snapshot = page.eval_on_selector("#color-histogram", "c => c.toDataURL()")
    assert snapshot.startswith("data:image/png;base64,")
    assert len(snapshot) > 2000
    assert page.pageerrors == []


def test_gradient_histogram_recolors_when_the_ramp_changes(meta_page):
    """Bars carry the color their values get, so narrowing the ramp repaints them."""
    page = meta_page
    _open_gradient_picker(page)
    before = page.eval_on_selector("#color-histogram", "c => c.toDataURL()")

    _set_stop(page, 0, "value", "6")
    page.wait_for_function(
        "prev => document.getElementById('color-histogram').toDataURL() !== prev",
        arg=before,
    )
    # Binning is over the data range, so clamped values stay on the chart -- that
    # they now share one flat color is the point.
    assert _histogram(page)["lowEdge"] == 1
    assert page.pageerrors == []


def test_gradient_histogram_counts_values_with_no_data(meta_page):
    """Nodes with no value are excluded from the bins and reported separately."""
    page = meta_page
    page.evaluate(
        """() => {
            const columnIndex = state.metadataColumnIndexByName.get('score');
            state.metadataByNodeIndex[0][columnIndex] = null;
            state.metadataByNodeIndex[1][columnIndex] = null;
            rebuildMetadataCaches();
        }"""
    )
    _open_gradient_picker(page)

    histogram = _histogram(page)
    assert histogram["total"] == 4
    assert histogram["missing"] == 2
    assert sum(histogram["counts"]) == 4
    assert "2 with no value" in page.eval_on_selector(
        "#color-histogram-note", "e => e.textContent")
    assert page.pageerrors == []


def test_gradient_histogram_handles_one_repeated_value(meta_page):
    """A zero-width range collapses to a single centered bin rather than dividing by zero."""
    page = meta_page
    page.evaluate(
        """() => {
            const columnIndex = state.metadataColumnIndexByName.get('score');
            state.metadataByNodeIndex.forEach(row => { row[columnIndex] = 4; });
            rebuildMetadataCaches();
        }"""
    )
    _open_gradient_picker(page)

    histogram = _histogram(page)
    assert histogram["binCount"] == 1
    assert histogram["counts"] == [6]
    assert [histogram["lowEdge"], histogram["highEdge"]] == [3.5, 4.5]
    assert page.pageerrors == []


def test_gradient_histogram_bins_floats_across_the_data_range(meta_page):
    """A float column falls back to even bins spanning exactly min..max."""
    page = meta_page
    page.evaluate(
        """() => {
            const columnIndex = state.metadataColumnIndexByName.get('score');
            state.metadataColumnByName.get('score').type = 'float';
            [0.5, 1.25, 2.75, 3.0, 4.5, 6.25].forEach((value, nodeIndex) => {
                state.metadataByNodeIndex[nodeIndex][columnIndex] = value;
            });
            rebuildMetadataCaches();
        }"""
    )
    _open_gradient_picker(page)

    histogram = _histogram(page)
    # sqrt(6) rounds up to 3, below the floor of 8 bins.
    assert histogram["binCount"] == 8
    assert [histogram["lowEdge"], histogram["highEdge"]] == [0.5, 6.25]
    assert sum(histogram["counts"]) == 6
    assert page.pageerrors == []


def test_gradient_histogram_absent_for_a_categorical_column(meta_page):
    """The chart belongs to the gradient panel, which discrete mode hides entirely."""
    page = meta_page
    page.select_option("#color-by", "family")
    _open_color_picker(page)
    page.wait_for_selector("#color-picker-discrete:not([hidden])")

    assert page.is_hidden("#color-picker-continuous")
    assert page.is_hidden("#color-histogram-wrap")
    assert page.pageerrors == []


def _knob_percent(page, index):
    """A stop knob's left offset as a float percentage of the shared axis."""
    left = page.eval_on_selector(
        f'#color-range-slider .cp-range-knob[data-stop-index="{index}"]', "e => e.style.left")
    return float(left.rstrip("%"))


def _ramp_geometry(page):
    """Content-box left edge and width of each layer of the ramp stack."""
    return page.evaluate(
        """() => ['color-histogram', 'color-gradient-preview', 'color-range-slider']
            .map(id => {
                const element = document.getElementById(id);
                const rect = element.getBoundingClientRect();
                const border = parseFloat(getComputedStyle(element).borderLeftWidth) || 0;
                return [rect.left + border, rect.width - (border * 2)];
            })"""
    )


def test_gradient_ramp_layers_share_one_axis(meta_page):
    """Histogram, gradient bar and slider must line up, or the knobs lie."""
    page = meta_page
    _open_gradient_picker(page)
    histogram, gradient, slider = _ramp_geometry(page)
    assert abs(histogram[0] - gradient[0]) < 0.5
    assert abs(histogram[0] - slider[0]) < 0.5
    assert abs(histogram[1] - gradient[1]) < 0.5
    assert abs(histogram[1] - slider[1]) < 0.5
    assert page.pageerrors == []


def test_gradient_bounds_start_at_the_ends_of_the_axis(meta_page):
    """Untouched, the ramp spans the data, so the knobs sit at 0% and 100%."""
    page = meta_page
    _open_gradient_picker(page)
    assert _knob_percent(page, 0) == 0
    assert _knob_percent(page, 1) == 100
    # The axis labels describe the data extent, not the ramp's own ends.
    assert page.eval_on_selector("#color-min-label", "e => e.textContent") == "1"
    assert page.eval_on_selector("#color-max-label", "e => e.textContent") == "9"
    assert page.pageerrors == []


def test_narrowing_the_ramp_moves_the_knobs_not_the_axis(meta_page):
    """Bounds inside the data range pull the knobs in and leave the axis alone."""
    page = meta_page
    _open_gradient_picker(page)
    _set_stop(page, 0, "value", "3")
    _set_stop(page, 1, "value", "7")

    domain = page.evaluate("() => gradientRangeDomain()")
    assert [domain["low"], domain["high"]] == [1, 9]
    # (3 - 1) / 8 and (7 - 1) / 8.
    assert _knob_percent(page, 0) == pytest.approx(25, abs=0.02)
    assert _knob_percent(page, 1) == pytest.approx(75, abs=0.02)
    assert page.eval_on_selector("#color-min-label", "e => e.textContent") == "1"
    assert page.eval_on_selector("#color-max-label", "e => e.textContent") == "9"
    assert page.pageerrors == []


def test_gradient_bar_holds_its_end_colors_outside_the_bounds(meta_page):
    """The bar spans the data axis, so clamped stretches show as flat color.

    The stops carry the bounds' real positions on that axis; CSS holds the first
    and last colors flat outside them, which is what numericColor does to values
    beyond the bounds.
    """
    page = meta_page
    _open_gradient_picker(page)
    _set_stop(page, 0, "value", "3")
    _set_stop(page, 1, "value", "7")

    background = page.eval_on_selector(
        "#color-gradient-preview", "e => getComputedStyle(e).backgroundImage")
    assert "25%" in background
    assert "75%" in background
    assert page.pageerrors == []


def test_adding_and_removing_stops_tracks_the_knobs(meta_page):
    """Rows and knobs are two views of one list, so they change together."""
    page = meta_page
    _open_gradient_picker(page)
    assert page.locator(".cp-stop-row").count() == 2
    assert page.locator("#color-range-slider .cp-range-knob").count() == 2
    # The two ends define the ramp, so neither carries a remove button.
    assert page.locator("[data-stop-remove]").count() == 0

    page.click("#color-add-stop")

    assert page.locator(".cp-stop-row").count() == 3
    assert page.locator("#color-range-slider .cp-range-knob").count() == 3
    assert page.locator("[data-stop-remove]").count() == 1
    # Split into the widest gap, so the new stop lands mid-range.
    assert _stop_values(page) == [1, 5, 9]
    assert _knob_percent(page, 1) == pytest.approx(50, abs=0.02)
    assert page.eval_on_selector_all(
        ".cp-stop-role", "els => els.map(e => e.textContent)") == ["Min", "Stop 2", "Max"]

    page.click('.cp-stop-row[data-stop-index="1"] [data-stop-remove]')

    assert page.locator(".cp-stop-row").count() == 2
    assert page.locator("#color-range-slider .cp-range-knob").count() == 2
    assert _stop_values(page) == [1, 9]
    assert page.pageerrors == []


def test_a_new_stop_takes_the_color_the_ramp_already_has(meta_page):
    """Adding a stop should shape the ramp you have, not restate it.

    Sampled across the whole range, no value shifts by more than one step per
    channel -- the residue of lerpHexColor rounding the new stop's color to 8
    bits, not a change in the ramp. (The *first* add is exempt: it replaces the
    built-in HSL sweep with an interpolation between its two ends.)
    """
    page = meta_page
    _open_gradient_picker(page)
    page.click("#color-add-stop")          # leaves the built-in ramp behind

    sample = """() => {
        const palette = customPalette('score');
        const colors = [];
        for (let value = 1; value <= 9; value += 0.1) {
            colors.push(numericColor(value, 1, 9, palette));
        }
        return colors;
    }"""
    before = page.evaluate(sample)

    page.click("#color-add-stop")

    after = page.evaluate(sample)
    assert len(after) == len(before)
    for old_hex, new_hex in zip(before, after):
        for channel in range(3):
            old_value = int(old_hex[1 + (channel * 2):3 + (channel * 2)], 16)
            new_value = int(new_hex[1 + (channel * 2):3 + (channel * 2)], 16)
            assert abs(old_value - new_value) <= 1, (old_hex, new_hex)
    # The added stop splits [1, 5] and takes that segment's midpoint color.
    assert _stop_values(page) == [1, 3, 5, 9]
    assert page.evaluate(
        """() => lerpHexColor(state.gradientStops[0].color,
                              state.gradientStops[2].color, 0.5)""") == _stop_colors(page)[1]
    assert page.pageerrors == []


def test_an_intermediate_stop_colors_its_own_value(meta_page):
    """A stop's color is what its value gets, whatever the ends are."""
    page = meta_page
    _open_gradient_picker(page)
    page.click("#color-add-stop")
    _set_stop(page, 1, "color", "#ff0000")

    assert page.evaluate(
        """() => numericColor(5, state.colorHistogram.min, state.colorHistogram.max,
                              customPalette('score'))""") == "#ff0000"
    # ...and the ends are untouched by it.
    assert page.evaluate(
        """() => numericColor(1, state.colorHistogram.min, state.colorHistogram.max,
                              customPalette('score'))""") == _stop_colors(page)[0]
    assert page.pageerrors == []


def test_the_two_end_stops_cannot_be_removed(meta_page):
    """The ramp needs two ends, so only intermediates are removable."""
    page = meta_page
    _open_gradient_picker(page)
    page.click("#color-add-stop")

    removable = page.eval_on_selector_all(
        ".cp-stop-row",
        "els => els.map(row => row.querySelector('[data-stop-remove]') !== null)")
    assert removable == [False, True, False]
    # Nothing in the API removes one either.
    page.evaluate("() => { removeGradientStop(0); removeGradientStop(2); }")
    assert len(_stop_values(page)) == 3
    assert page.pageerrors == []


def test_stop_count_is_capped(meta_page):
    """Add is refused past the cap rather than crowding the rows and track."""
    page = meta_page
    _open_gradient_picker(page)
    cap = page.evaluate("() => MAX_GRADIENT_STOPS")
    # Two stops already exist, so `cap - 2` adds reach the limit exactly.
    for _ in range(cap - 2):
        page.click("#color-add-stop")

    assert len(_stop_values(page)) == cap
    assert page.eval_on_selector("#color-add-stop", "e => e.disabled") is True
    page.evaluate("() => addGradientStop()")
    assert len(_stop_values(page)) == cap
    assert "at most" in page.eval_on_selector("#bundle-status", "e => e.textContent")
    assert page.pageerrors == []


def test_typing_a_stop_past_its_neighbour_reorders_on_commit(meta_page):
    """Re-sorting mid-keystroke would yank the row out from under the cursor."""
    page = meta_page
    _open_gradient_picker(page)
    page.click("#color-add-stop")

    _set_stop(page, 1, "value", "12")
    assert _stop_values(page) == [1, 12, 9], "reordered while still typing"

    _stop_field(page, 1, "value").dispatch_event("change")
    assert _stop_values(page) == [1, 9, 12]
    assert page.eval_on_selector_all(
        ".cp-stop-role", "els => els.map(e => e.textContent)") == ["Min", "Stop 2", "Max"]
    assert page.pageerrors == []


def _drag_knob_to_fraction(page, index, fraction):
    """Drag a stop's knob to a fraction of the way along the slider track."""
    knob = page.query_selector(
        f'#color-range-slider .cp-range-knob[data-stop-index="{index}"]')
    track = page.query_selector("#color-range-slider").bounding_box()
    box = knob.bounding_box()
    middle = box["y"] + (box["height"] / 2)
    page.mouse.move(box["x"] + (box["width"] / 2), middle)
    page.mouse.down()
    page.mouse.move(track["x"] + (track["width"] * fraction), middle, steps=8)
    page.mouse.up()


def test_dragging_a_knob_recolors_the_network(meta_page):
    """A knob drag runs the same apply() path as typing in the number input."""
    page = meta_page
    _open_gradient_picker(page)
    before = _canvas_snapshot(page)

    _drag_knob_to_fraction(page, 0, 0.5)

    _wait_for_canvas_change(page, before)
    # An int column snaps to whole numbers: halfway along 1..9 is 5.
    assert _stop_values(page)[0] == 5
    assert page.input_value(
        '.cp-stop-row[data-stop-index="0"] [data-stop-field="value"]') == "5"
    assert page.evaluate(
        "() => state.customPalettes[paletteKey('score')].stops[0].value") == 5
    assert page.pageerrors == []


def test_a_dragged_knob_cannot_cross_its_neighbour(meta_page):
    """Low is held at high (and vice versa) so the ramp cannot invert."""
    page = meta_page
    _open_gradient_picker(page)
    _set_stop(page, 1, "value", "6", commit=True)

    _drag_knob_to_fraction(page, 0, 1.0)
    assert _stop_values(page) == [6, 6]

    _set_stop(page, 0, "value", "3", commit=True)
    _drag_knob_to_fraction(page, 1, 0.0)
    assert _stop_values(page) == [3, 3]
    assert page.pageerrors == []


def test_knobs_are_keyboard_operable(meta_page):
    """Arrow keys step the bound; Home/End jump to the ends of the axis."""
    page = meta_page
    _open_gradient_picker(page)
    page.eval_on_selector(
        '#color-range-slider .cp-range-knob[data-stop-index="0"]', "e => e.focus()")

    page.keyboard.press("ArrowRight")
    assert _stop_values(page)[0] == 2
    page.keyboard.press("ArrowLeft")
    assert _stop_values(page)[0] == 1
    page.keyboard.press("Home")
    assert _stop_values(page)[0] == 1

    page.eval_on_selector(
        '#color-range-slider .cp-range-knob[data-stop-index="1"]', "e => e.focus()")
    page.keyboard.press("End")
    assert _stop_values(page)[1] == 9
    assert page.evaluate(
        "() => document.activeElement.dataset.stopIndex") == "1"
    assert page.pageerrors == []


def test_delete_key_removes_the_focused_intermediate_knob(meta_page):
    """Delete on a knob drops that stop; on an end knob it does nothing."""
    page = meta_page
    _open_gradient_picker(page)
    page.click("#color-add-stop")
    page.eval_on_selector(
        '#color-range-slider .cp-range-knob[data-stop-index="1"]', "e => e.focus()")

    page.keyboard.press("Delete")
    assert _stop_values(page) == [1, 9]

    page.eval_on_selector(
        '#color-range-slider .cp-range-knob[data-stop-index="0"]', "e => e.focus()")
    page.keyboard.press("Delete")
    assert _stop_values(page) == [1, 9]
    assert page.pageerrors == []


def test_reset_values_returns_the_knobs_to_the_axis_ends(meta_page):
    """"Reset values to data range" puts the ramp back over the whole axis."""
    page = meta_page
    _open_gradient_picker(page)
    _set_stop(page, 0, "value", "4", commit=True)
    _set_stop(page, 1, "value", "6", commit=True)
    assert _knob_percent(page, 0) > 0

    page.click("#color-reset-values")

    assert _knob_percent(page, 0) == 0
    assert _knob_percent(page, 1) == 100
    assert page.pageerrors == []


def test_float_column_knobs_do_not_snap_to_integers(meta_page):
    """Only int columns snap; a float column keeps fractional bounds."""
    page = meta_page
    page.evaluate(
        """() => {
            const columnIndex = state.metadataColumnIndexByName.get('score');
            state.metadataColumnByName.get('score').type = 'float';
            [0.0, 2.0, 4.0, 6.0, 8.0, 10.0].forEach((value, nodeIndex) => {
                state.metadataByNodeIndex[nodeIndex][columnIndex] = value;
            });
            rebuildMetadataCaches();
        }"""
    )
    _open_gradient_picker(page)
    assert page.evaluate("() => state.colorHistogram.integer") is False

    _drag_knob_to_fraction(page, 0, 0.25)

    value = float(_stop_values(page)[0])
    assert 2.0 < value < 3.0, value
    assert page.pageerrors == []


def test_a_pre_stops_session_keeps_its_colors(meta_page):
    """Sessions saved before gradients were stop lists must still load.

    They carry {low, mid, high} with separate lowValue/midValue/highValue
    bounds; converting on read is what lets the rest of the viewer know only
    about stops.
    """
    page = meta_page
    converted = page.evaluate(
        """() => {
            const legacy = {'score': {type: 'numeric', low: '#000080', mid: '#00ff00',
                                      high: '#ff0000', lowValue: 10, midValue: 40,
                                      highValue: 80, nullColor: '#cccccc'}};
            normalizeCustomPalettes(legacy);
            const palette = legacy['score'];
            return {
                stops: palette.stops,
                nullColor: palette.nullColor,
                at_low: numericColor(10, 1, 99, palette),
                at_mid: numericColor(40, 1, 99, palette),
                at_high: numericColor(80, 1, 99, palette),
                below: numericColor(5, 1, 99, palette),
                above: numericColor(90, 1, 99, palette),
            };
        }"""
    )
    assert converted["stops"] == [
        {"value": 10, "color": "#000080"},
        {"value": 40, "color": "#00ff00"},
        {"value": 80, "color": "#ff0000"},
    ]
    assert converted["nullColor"] == "#cccccc"
    # Each old color still lands on the value it was bound to...
    assert converted["at_low"] == "#000080"
    assert converted["at_mid"] == "#00ff00"
    assert converted["at_high"] == "#ff0000"
    # ...and the old bounds still clamp.
    assert converted["below"] == "#000080"
    assert converted["above"] == "#ff0000"
    assert page.pageerrors == []


def test_a_two_color_pre_stops_session_converts_to_two_stops(meta_page):
    """No midpoint color in the old shape means no third stop."""
    page = meta_page
    stops = page.evaluate(
        """() => {
            const legacy = {'score': {type: 'numeric', low: '#000000', mid: null,
                                      high: '#ffffff', lowValue: null, midValue: 4,
                                      highValue: null, nullColor: null}};
            normalizeCustomPalettes(legacy);
            return legacy['score'].stops;
        }"""
    )
    # Null bounds fall back to the data range, and the stale mid VALUE that used
    # to bend the ramp is dropped along with the absent mid color.
    assert [stop["color"] for stop in stops] == ["#000000", "#ffffff"]
    assert len(stops) == 2
    assert page.pageerrors == []


def test_gradient_stops_survive_a_session_round_trip(meta_page, tmp_path):
    """A tuned multi-stop ramp reloads exactly."""
    page = meta_page
    _open_gradient_picker(page)
    page.click("#color-add-stop")
    page.click("#color-add-stop")
    _set_stop(page, 1, "color", "#ff0000")
    _set_stop(page, 2, "color", "#00ff00")
    before = _stop_values(page)
    colors = _stop_colors(page)
    page.click("#color-picker-close")

    saved = _save_session(page, tmp_path)
    _load_bundle_file(page, saved)
    _open_gradient_picker(page)

    assert _stop_values(page) == before
    assert _stop_colors(page) == colors
    assert page.locator(".cp-stop-row").count() == 4
    assert page.pageerrors == []


def test_legend_gradient_carries_one_svg_stop_per_palette_stop(meta_page):
    """The exported legend is the ramp exactly, not a sampled approximation."""
    page = meta_page
    _open_gradient_picker(page)
    page.click("#color-add-stop")
    _set_stop(page, 1, "color", "#ff0000")
    _set_stop(page, 1, "value", "3", commit=True)

    svg = page.evaluate("() => buildLegendSVG().svg")
    offsets = re.findall(r'<stop offset="([\d.]+)%" stop-color="([^"]+)"', svg)
    assert len(offsets) == 3
    # Stop values 1, 3, 9 across a 1..9 bar put the middle stop a quarter along.
    assert [round(float(offset), 4) for offset, _ in offsets] == [0.0, 25.0, 100.0]
    assert offsets[1][1] == "#ff0000"
    assert page.pageerrors == []


def test_legend_ticks_drop_only_where_they_would_collide(meta_page):
    """A crowded ramp still labels its ends, whatever it does in between."""
    page = meta_page
    _open_gradient_picker(page)
    cap = page.evaluate("() => MAX_GRADIENT_STOPS")
    for _ in range(cap - 2):
        page.click("#color-add-stop")

    svg = page.evaluate("() => buildLegendSVG().svg")
    ticks = re.findall(r'text-anchor="(\w+)">([^<]*)</text>', svg)
    assert ticks[0] == ("start", "1")
    assert ticks[-1] == ("end", "9")
    assert 2 <= len(ticks) <= cap
    assert page.pageerrors == []


def _hex_field(page, index):
    return page.locator(
        f'.cp-stop-row[data-stop-index="{index}"] [data-stop-field="hex"]')


def _type_hex(page, index, text):
    """Type into a stop's hex field and leave it, as a user would."""
    node = _hex_field(page, index)
    node.fill(text)
    node.dispatch_event("input")
    node.dispatch_event("change")


def test_hex_field_shows_the_canonical_spelling(meta_page):
    """Each stop carries a copy-pastable hex code beside its swatch."""
    page = meta_page
    _open_gradient_picker(page)
    for index in range(2):
        shown = _hex_field(page, index).input_value()
        assert re.fullmatch(r"#[0-9A-F]{6}", shown), shown
        # Same color as the swatch, just spelled canonically.
        assert shown.lower() == _stop_colors(page)[index]
    assert page.pageerrors == []


@pytest.mark.parametrize("typed, expected", [
    ("ff0000", "#FF0000"),          # bare six digits
    ("#00ff00", "#00FF00"),         # lower case
    ("  #0000FF  ", "#0000FF"),     # surrounding whitespace
    ("abc", "#AABBCC"),             # three-digit shorthand, no hash
    ("#F0F", "#FF00FF"),            # three-digit shorthand with hash
])
def test_hex_field_normalizes_on_leaving_the_field(meta_page, typed, expected):
    """Input is taken loosely and rewritten canonically once the field is left."""
    page = meta_page
    _open_gradient_picker(page)

    _type_hex(page, 0, typed)

    assert _hex_field(page, 0).input_value() == expected
    assert _stop_colors(page)[0] == expected.lower()
    # The swatch is the same color by another name, so it follows.
    assert page.input_value(
        '.cp-stop-row[data-stop-index="0"] [data-stop-field="color"]') == expected.lower()
    assert page.pageerrors == []


def test_hex_field_recolors_the_network(meta_page):
    """Typing a color is a gradient edit like any other."""
    page = meta_page
    _open_gradient_picker(page)
    before = _canvas_snapshot(page)

    _type_hex(page, 0, "#ff0000")

    _wait_for_canvas_change(page, before)
    assert page.evaluate(
        "() => state.customPalettes[paletteKey('score')].stops[0].color") == "#ff0000"
    assert page.pageerrors == []


def test_hex_field_reverts_text_that_is_not_a_color(meta_page):
    """A rejected entry restores the stop's color rather than being left standing."""
    page = meta_page
    _open_gradient_picker(page)
    _type_hex(page, 0, "#123456")

    _type_hex(page, 0, "aabbccdd")

    assert _hex_field(page, 0).input_value() == "#123456"
    assert _stop_colors(page)[0] == "#123456"
    assert "not a hex color" in page.eval_on_selector("#bundle-status", "e => e.textContent")
    assert page.pageerrors == []


def test_clearing_the_hex_field_cancels_quietly(meta_page):
    """Emptying the field and tabbing away reads as canceling, not as an error."""
    page = meta_page
    _open_gradient_picker(page)
    _type_hex(page, 0, "#123456")
    page.evaluate("() => setStatus('unchanged')")

    _type_hex(page, 0, "")

    assert _hex_field(page, 0).input_value() == "#123456"
    assert page.eval_on_selector("#bundle-status", "e => e.textContent") == "unchanged"
    assert page.pageerrors == []


def test_editing_the_swatch_updates_the_hex_field(meta_page):
    """The swatch and the hex code are two views of one color."""
    page = meta_page
    _open_gradient_picker(page)

    _set_stop(page, 1, "color", "#00ff7f")

    assert _hex_field(page, 1).input_value() == "#00FF7F"
    assert page.pageerrors == []


def test_added_stops_get_a_hex_field_too(meta_page):
    """Rows are generated, so a new stop is not a special case."""
    page = meta_page
    _open_gradient_picker(page)
    page.click("#color-add-stop")

    assert page.locator("[data-stop-field='hex']").count() == 3
    assert _hex_field(page, 1).input_value().lower() == _stop_colors(page)[1]
    assert page.pageerrors == []


def test_color_table_load_accepts_shorthand_hex(meta_page, tmp_path):
    """The looser parser is shared, so color tables take shorthand as well."""
    page = meta_page
    page.select_option("#color-by", "family")
    table = tmp_path / "colors.tsv"
    table.write_text("alpha\t#f00\nbeta\t0f0\ngamma\t#0000ff\n")

    page.set_input_files("#color-table-file", str(table))
    page.wait_for_function(
        "() => document.getElementById('bundle-status').textContent.startsWith('Loaded 3')"
    )

    assert page.evaluate(
        "() => state.customPalettes[paletteKey('family')].colors") == {
            "alpha": "#FF0000", "beta": "#00FF00", "gamma": "#0000FF"}
    assert page.pageerrors == []


# ---------------------------------------------------------------------------
# Per-column charts and frequency tables
# ---------------------------------------------------------------------------


def _open_column_chart_menu(page, column):
    """Click a column header's chart glyph and wait for its menu."""
    page.click(f'.metadata-chart-button[data-chart-column="{column}"]')
    page.wait_for_selector("#column-chart-menu:not([hidden])")


def _column_chart_menu_kinds(page):
    return page.eval_on_selector_all(
        "#column-chart-menu button", "els => els.map(e => e.dataset.chartKind)"
    )


def _open_column_chart(page, column, kind):
    """Open the chart dialog for one column/kind through the header menu."""
    _open_column_chart_menu(page, column)
    page.click(f'#column-chart-menu button[data-chart-kind="{kind}"]')
    page.wait_for_selector("#column-chart-overlay:not([hidden])")
    page.wait_for_function(
        "kind => state.columnChartArtifact && state.columnChartArtifact.kind === kind",
        arg=kind,
    )


def _set_column_chart_kind(page, kind):
    page.select_option("#column-chart-kind", kind)
    page.wait_for_function(
        "kind => state.columnChartArtifact && state.columnChartArtifact.kind === kind",
        arg=kind,
    )


def _column_chart_svg(page):
    return page.inner_html("#column-chart-preview")


def _column_chart_rows(page):
    """The rendered frequency/summary table as a list of row cell lists."""
    return page.eval_on_selector_all(
        "#column-chart-preview tbody tr",
        "rows => rows.map(r => Array.from(r.cells).map(c => c.textContent))",
    )


def test_column_chart_menu_adapts_to_column_type(meta_page):
    """A text column offers no histogram; a numeric one offers every kind."""
    page = meta_page
    _open_column_chart_menu(page, "family")
    assert _column_chart_menu_kinds(page) == ["frequency", "bar", "pie", "summary"]

    page.keyboard.press("Escape")
    page.wait_for_selector("#column-chart-menu", state="hidden")
    assert page.get_attribute(
        '.metadata-chart-button[data-chart-column="family"]', "aria-expanded"
    ) == "false"

    _open_column_chart_menu(page, "score")
    assert _column_chart_menu_kinds(page) == [
        "frequency", "bar", "pie", "histogram", "box", "ecdf", "summary",
    ]
    assert page.pageerrors == []


def test_node_id_header_has_no_chart_button(meta_page):
    """node_id is unique by construction, so every chart of it is degenerate."""
    page = meta_page
    assert page.locator(".metadata-chart-button").count() == 2
    assert page.locator(
        "#metadata-table thead th:nth-child(1) .metadata-chart-button"
    ).count() == 0
    assert page.pageerrors == []


def test_chart_glyph_click_does_not_sort(meta_page):
    """The glyph sits inside the header cell; clicking it must not sort."""
    page = meta_page
    _open_column_chart_menu(page, "family")

    assert page.evaluate("() => state.metadataSort.columnKey") is None
    assert page.eval_on_selector(
        "#metadata-table thead th:nth-child(2)", "e => e.getAttribute('aria-sort')"
    ) == "none"
    assert page.pageerrors == []


def test_frequency_table_counts_and_inherits_colors(meta_page):
    """Counts come from the scoped rows; swatches from the column's palette."""
    page = meta_page
    _open_column_chart(page, "family", "frequency")

    assert _column_chart_rows(page) == [
        ["", "alpha", "2", "33.3%", "2", "100.0%"],
        ["", "beta", "2", "33.3%", "2", "100.0%"],
        ["", "gamma", "2", "33.3%", "2", "100.0%"],
    ]
    expected = page.evaluate(
        "() => ['alpha', 'beta', 'gamma'].map("
        "  value => cssToHex(categoricalColor(value, customPalette('family'))))"
    )
    swatches = page.eval_on_selector_all(
        ".cc-swatch", "els => els.map(e => cssToHex(e.style.background))"
    )
    assert swatches == expected
    assert page.pageerrors == []


def test_frequency_table_inherits_a_custom_palette(meta_page):
    """A palette edit reaches an open chart, not just the network canvas."""
    page = meta_page
    _open_column_chart(page, "family", "frequency")
    page.evaluate("""() => {
        ensureCategoricalPalette('family').colors = {alpha: '#112233'};
        rebuildNodeColorCache();
    }""")

    assert page.eval_on_selector(
        ".cc-swatch", "e => cssToHex(e.style.background)"
    ) == "#112233"
    assert "#112233" in page.evaluate("() => columnChartTSV()")
    assert page.pageerrors == []


def test_column_chart_follows_the_selection(meta_page):
    """Charts summarize the rows the table shows, i.e. the selection."""
    page = meta_page
    _open_column_chart(page, "family", "frequency")
    page.evaluate("() => { state.selectedNodeIndices = new Set([0, 1, 2]); updateMetadataTable(); }")

    assert _column_chart_rows(page) == [
        ["", "alpha", "2", "66.7%", "2", "100.0%"],
        ["", "beta", "1", "33.3%", "2", "50.0%"],
    ]
    assert page.text_content("#column-chart-note").startswith("3 of 6 nodes")
    assert "selection" in page.text_content("#column-chart-note")
    assert page.pageerrors == []


def test_frequency_table_reports_each_value_global_share(meta_page):
    """The second percentage answers "how much of this value did I catch?".

    A/B are alpha and C/D are beta, so selecting A, B and C makes the two
    denominators differ: beta is a third of the selection but half of the betas
    in the network, and only the global column can say the second.
    """
    page = meta_page
    _open_column_chart(page, "family", "frequency")
    page.evaluate("() => { state.selectedNodeIndices = new Set([0, 1, 2]); updateMetadataTable(); }")

    assert _column_chart_rows(page) == [
        ["", "alpha", "2", "66.7%", "2", "100.0%"],
        ["", "beta", "1", "33.3%", "2", "50.0%"],
    ]
    assert "% of all" in page.text_content("#column-chart-note")
    assert page.evaluate("() => columnChartTSV()").splitlines()[:3] == [
        "value\tcount\tpercent\tcount_all_nodes\tpercent_of_all\tcolor",
        "alpha\t2\t66.6667\t2\t100.0000\t#1f77b4",
        "beta\t1\t33.3333\t2\t50.0000\t#ff7f0e",
    ]
    # The charted kinds share the row model, so their TSVs carry it too, and the
    # gloss stays out of a chart that does not draw the column.
    _set_column_chart_kind(page, "bar")
    assert page.evaluate("() => columnChartTSV()").splitlines()[1].endswith("2\t100.0000\t#1f77b4")
    assert "% of all" not in page.text_content("#column-chart-note")
    assert page.pageerrors == []


def test_frequency_table_global_share_is_whole_when_nothing_is_selected(meta_page):
    """With the whole bundle in scope every value is entirely in scope."""
    page = meta_page
    _open_column_chart(page, "family", "frequency")

    assert [row[-1] for row in _column_chart_rows(page)] == ["100.0%"] * 3
    # Nothing to gloss when both denominators are the same number.
    assert "% of all" not in page.text_content("#column-chart-note")
    assert page.pageerrors == []


def test_column_chart_honours_the_table_filter(meta_page):
    """The same row set the header's Copy button exports, filter included."""
    page = meta_page
    _open_column_chart(page, "family", "frequency")
    page.fill("#metadata-filter", "alpha")
    # The filter is debounced, so wait for the re-render to reach the chart.
    page.wait_for_function("() => state.columnChartArtifact.notes[0].startsWith('2 of 6')")

    assert _column_chart_rows(page) == [["", "alpha", "2", "100.0%", "2", "100.0%"]]
    assert 'filter "alpha"' in page.text_content("#column-chart-note")
    assert page.evaluate("() => columnChartNodeIndices().length") == 2
    assert page.pageerrors == []


def test_pie_chart_single_value_draws_a_circle(meta_page):
    """An arc whose start and end angles are equal draws nothing at all."""
    page = meta_page
    page.fill("#metadata-filter", "alpha")
    # The filter is debounced, and the re-render it eventually runs rebuilds the
    # header -- closing an open chart menu with it. Waiting only for the scope to
    # narrow would open the menu inside that window and have it shut mid-click.
    page.wait_for_function(
        "() => state.pendingMetadataFilterTimer === null"
        " && document.querySelectorAll('#metadata-table tbody tr[data-node-index]').length === 2"
    )
    _open_column_chart(page, "family", "pie")

    svg = _column_chart_svg(page)
    assert "<circle" in svg
    assert "<path" not in svg
    assert "NaN" not in svg
    assert page.pageerrors == []


def test_bar_chart_draws_one_bar_per_value(meta_page):
    page = meta_page
    _open_column_chart(page, "family", "bar")

    # One bar rect per value, plus the SVG's own white background rect.
    assert page.eval_on_selector_all("#column-chart-preview rect", "els => els.length") == 4
    assert page.pageerrors == []


def test_histogram_bins_scoped_integers_one_per_value(meta_page):
    """score is 1,2,3,5,8,9: nine bins, one per integer, three of them empty."""
    page = meta_page
    _open_column_chart(page, "score", "histogram")

    histogram = page.evaluate(
        "() => scopedColumnHistogram('score', columnChartNodeIndices())"
    )
    assert histogram["binCount"] == 9
    assert histogram["counts"] == [1, 1, 1, 0, 1, 0, 0, 1, 1]
    assert histogram["total"] == 6
    # Empty bins draw nothing, so six bars plus the background rect.
    assert page.eval_on_selector_all("#column-chart-preview rect", "els => els.length") == 7
    # One bin per integer means the bins are values, not edges.
    assert page.evaluate("() => columnChartTSV()").splitlines()[:3] == [
        "value\tcount", "1\t1", "2\t1",
    ]
    assert page.pageerrors == []


def test_histogram_bars_take_the_column_gradient(meta_page):
    """Bars are colored the way the nodes are, so the chart doubles as ramp
    feedback -- which means the whole-column range, not the scoped one."""
    page = meta_page
    _open_column_chart(page, "score", "histogram")
    before = _column_chart_svg(page)

    page.evaluate("""() => {
        ensureNumericPalette('score').stops = [
            {value: 1, color: '#000000'}, {value: 9, color: '#FFFFFF'}];
        rebuildNodeColorCache();
    }""")

    after = _column_chart_svg(page)
    assert after != before
    # Every bar now sits on a black-to-white ramp, so every fill is a gray.
    # They are bin-*center* colors, so none is pure black or pure white.
    bar_fills = re.findall(r'<rect [^>]*fill="(#[0-9a-f]{6})"', after)
    assert len(bar_fills) == 6
    assert all(fill[1:3] == fill[3:5] == fill[5:7] for fill in bar_fills)
    assert bar_fills == sorted(bar_fills)
    assert page.pageerrors == []


def test_summary_statistics_match_the_five_number_summary(meta_page):
    """score is 1,2,3,5,8,9; quartiles interpolate like numpy's default."""
    page = meta_page
    _open_column_chart(page, "score", "summary")

    summary = page.evaluate(
        "() => numericSummary(columnNumericValues('score', columnChartNodeIndices()).values)"
    )
    assert summary["min"] == 1
    assert summary["q1"] == 2.25
    assert summary["median"] == 4
    assert summary["q3"] == 7.25
    assert summary["max"] == 9
    assert summary["outliers"] == []

    rows = dict(row[:2] for row in _column_chart_rows(page))
    assert rows["Rows"] == "6"
    assert rows["Distinct values"] == "6"
    assert rows["Median"] == "4"
    assert rows["1st quartile"] == "2.25"
    assert page.pageerrors == []


def test_box_plot_whiskers_stop_at_the_data(meta_page):
    """No outliers in this column, so the whiskers reach min and max."""
    page = meta_page
    _open_column_chart(page, "score", "box")

    assert page.evaluate("() => columnChartTSV()").splitlines() == [
        "statistic\tvalue",
        "count\t6",
        "minimum\t1",
        "low_whisker\t1",
        "q1\t2.25",
        "median\t4",
        "q3\t7.25",
        "high_whisker\t9",
        "maximum\t9",
        "iqr\t5",
        "outliers\t",
    ]
    # No outlier dots to draw when nothing falls outside the fences.
    assert page.eval_on_selector_all("#column-chart-preview circle", "els => els.length") == 0
    assert page.pageerrors == []


def test_ecdf_rises_monotonically_to_one(meta_page):
    page = meta_page
    _open_column_chart(page, "score", "ecdf")

    tsv = page.evaluate("() => columnChartTSV()").splitlines()[1:]
    fractions = [float(line.split("\t")[1]) for line in tsv]
    assert fractions == sorted(fractions)
    assert fractions[-1] == pytest.approx(1.0)

    # The curve is one path whose y coordinates never increase (y grows downward).
    path = re.search(r'<path d="([^"]+)"', _column_chart_svg(page)).group(1)
    ys = [float(value) for value in re.findall(r"[ML] [\d.]+ ([\d.]+)", path)]
    assert ys == sorted(ys, reverse=True)
    assert page.pageerrors == []


def test_many_categories_roll_into_other(many_cat_page):
    """120 families: the pie caps its slices, the table lists them all."""
    page = many_cat_page
    _open_column_chart(page, "family", "pie")

    # Twelve most common plus one "Other" slice.
    assert page.eval_on_selector_all("#column-chart-preview path", "els => els.length") == 13
    assert "108 more rolled into" in page.text_content("#column-chart-note")
    # The rollup takes get_palette's neutral gray, not a category hue.
    assert 'fill="#bfbfbf"' in _column_chart_svg(page)

    _set_column_chart_kind(page, "frequency")
    assert len(_column_chart_rows(page)) == 120
    assert "rolled into" not in page.text_content("#column-chart-note")
    assert page.pageerrors == []


def test_numeric_categorical_column_uses_its_discrete_palette(numeric_categorical_page):
    """A numeric column colored as categories still charts as numbers, but its
    category colors come from the paletteKey-suffixed discrete palette."""
    page = numeric_categorical_page
    _open_column_chart_menu(page, "score")
    assert "histogram" in _column_chart_menu_kinds(page)

    page.click('#column-chart-menu button[data-chart-kind="pie"]')
    page.wait_for_selector("#column-chart-overlay:not([hidden])")
    page.evaluate("""() => {
        ensureCategoricalPalette('score').colors = {'5': '#445566'};
        rebuildNodeColorCache();
    }""")

    assert page.evaluate("() => paletteKey('score')") == "score\u0000categorical"
    assert 'fill="#445566"' in _column_chart_svg(page)
    assert page.pageerrors == []


def test_column_chart_exports_svg_png_and_tsv(meta_page, tmp_path):
    page = meta_page
    _open_column_chart(page, "family", "frequency")

    with page.expect_download() as download_info:
        page.click("#column-chart-export-svg")
    svg_download = download_info.value
    assert svg_download.suggested_filename == "Color_Test_Viewer_family_frequency.svg"
    svg_path = tmp_path / "chart.svg"
    svg_download.save_as(str(svg_path))
    # The XML prolog is added only on export: the same string is assigned to
    # innerHTML for the preview, which refuses a processing instruction.
    assert svg_path.read_text().startswith("<?xml version=")
    assert "<?xml" not in _column_chart_svg(page)

    with page.expect_download() as download_info:
        page.click("#column-chart-download-tsv")
    tsv_download = download_info.value
    assert tsv_download.suggested_filename == "Color_Test_Viewer_family_frequency.tsv"
    tsv_path = tmp_path / "chart.tsv"
    tsv_download.save_as(str(tsv_path))
    assert tsv_path.read_text().splitlines()[0] == (
        "value\tcount\tpercent\tcount_all_nodes\tpercent_of_all\tcolor")

    with page.expect_download() as download_info:
        page.click("#column-chart-export-png")
    png_download = download_info.value
    assert png_download.suggested_filename.endswith(".png")
    png_path = tmp_path / "chart.png"
    png_download.save_as(str(png_path))
    assert png_path.read_bytes().startswith(b"\x89PNG")
    assert page.pageerrors == []


def test_column_chart_tracks_a_metadata_edit(meta_page):
    page = meta_page
    _open_column_chart(page, "family", "frequency")
    # The dialog is modal, so drive the edit through the same functions the
    # cell editor calls rather than clicking the table underneath it.
    page.evaluate("""() => {
        setMetadataValue(0, 'family', 'delta');
        refreshAfterMetadataEdit();
    }""")
    page.wait_for_function(
        "() => state.columnChartArtifact.notes.some(n => n.includes('4 distinct'))"
    )

    assert _column_chart_rows(page) == [
        ["", "beta", "2", "33.3%", "2", "100.0%"],
        ["", "gamma", "2", "33.3%", "2", "100.0%"],
        ["", "alpha", "1", "16.7%", "1", "100.0%"],
        ["", "delta", "1", "16.7%", "1", "100.0%"],
    ]
    assert page.pageerrors == []


def test_column_chart_closes_when_its_column_is_deleted(meta_page):
    """A column can be deleted from under an open dialog."""
    page = meta_page
    _add_column(page, "notes")
    _open_column_chart(page, "notes", "frequency")

    page.evaluate("() => deleteMetadataColumn('notes')")
    page.wait_for_selector("#column-chart-overlay", state="hidden")

    assert page.evaluate("() => state.columnChart") is None
    assert page.pageerrors == []


def test_chart_dialog_closes_on_escape_and_on_the_backdrop(meta_page):
    """The two dismissals every other overlay in the viewer supports."""
    page = meta_page
    _open_column_chart(page, "family", "frequency")
    page.keyboard.press("Escape")
    page.wait_for_selector("#column-chart-overlay", state="hidden")
    assert page.evaluate("() => state.columnChart") is None
    assert page.evaluate("() => state.columnChartArtifact") is None

    _open_column_chart(page, "family", "frequency")
    # Click the overlay itself, outside the dialog, at the very top of the page.
    page.mouse.click(4, 4)
    page.wait_for_selector("#column-chart-overlay", state="hidden")
    assert page.pageerrors == []


def test_a_second_click_on_the_glyph_closes_the_menu(meta_page):
    page = meta_page
    _open_column_chart_menu(page, "family")

    page.click('.metadata-chart-button[data-chart-column="family"]')
    page.wait_for_selector("#column-chart-menu", state="hidden")
    assert page.evaluate("() => state.columnChartMenu") is None
    assert page.pageerrors == []


def test_clicking_another_glyph_reopens_the_menu_for_that_column(meta_page):
    page = meta_page
    _open_column_chart_menu(page, "family")

    page.click('.metadata-chart-button[data-chart-column="score"]')
    page.wait_for_function("() => state.columnChartMenu?.columnName === 'score'")
    assert "histogram" in _column_chart_menu_kinds(page)
    assert page.get_attribute(
        '.metadata-chart-button[data-chart-column="family"]', "aria-expanded"
    ) == "false"
    assert page.pageerrors == []


def test_chart_says_why_when_there_is_nothing_to_draw(meta_page):
    """A filter that matches nothing must explain itself, not draw a blank box."""
    page = meta_page
    _open_column_chart(page, "score", "histogram")
    page.fill("#metadata-filter", "zzzzz")
    page.wait_for_function("() => state.columnChartArtifact.empty")

    assert page.eval_on_selector("#column-chart-preview", "e => e.textContent") == (
        "No numeric values in these rows."
    )
    assert page.eval_on_selector_all("#column-chart-preview svg", "els => els.length") == 0
    for button in ("copy-tsv", "download-tsv", "export-svg", "export-png"):
        assert page.is_disabled(f"#column-chart-{button}")
    assert page.pageerrors == []


def test_an_all_null_column_charts_as_one_no_value_entry(meta_page):
    """A column nobody has filled in yet is a legitimate thing to chart."""
    page = meta_page
    page.evaluate("() => addMetadataColumn('empty', 'number')")
    _open_column_chart(page, "empty", "frequency")

    assert _column_chart_rows(page) == [["", "—", "6", "100.0%", "6", "100.0%"]]
    assert "6 with no value" in page.text_content("#column-chart-note")
    # The no-value color is categoricalColor's own fallback, so an empty cell
    # is the same color in the chart as its node is on the canvas.
    assert page.eval_on_selector(
        ".cc-swatch", "e => cssToHex(e.style.background)"
    ) == page.evaluate("() => cssToHex(categoricalColor(null, customPalette('empty')))")

    # ...but the numeric kinds have nothing to bin.
    _set_column_chart_kind(page, "histogram")
    assert page.evaluate("() => state.columnChartArtifact.empty") is True
    assert page.pageerrors == []


def test_a_single_repeated_value_does_not_break_the_numeric_charts(meta_page):
    """A zero-width numeric domain would divide by zero in every axis."""
    page = meta_page
    page.evaluate("""() => {
        addMetadataColumn('flat', 'number');
        for (let i = 0; i < 6; i++) { setMetadataValue(i, 'flat', '7'); }
        refreshAfterMetadataEdit({columnsChanged: true});
    }""")

    for kind in ("histogram", "box", "ecdf", "summary"):
        _open_column_chart(page, "flat", kind)
        svg = _column_chart_svg(page)
        assert "NaN" not in svg
        assert "Infinity" not in svg
        assert page.evaluate("() => state.columnChartArtifact.empty") is False
        page.keyboard.press("Escape")
        page.wait_for_selector("#column-chart-overlay", state="hidden")

    assert page.evaluate(
        "() => scopedColumnHistogram('flat', columnChartNodeIndices()).binCount"
    ) == 1
    assert page.pageerrors == []


def test_chart_menu_is_keyboard_operable(meta_page):
    """Opening focuses the first item; the arrows and Home/End walk the rest."""
    page = meta_page
    _open_column_chart_menu(page, "score")

    focused = "() => document.activeElement.dataset.chartKind"
    assert page.evaluate(focused) == "frequency"
    page.keyboard.press("ArrowDown")
    assert page.evaluate(focused) == "bar"
    page.keyboard.press("ArrowUp")
    assert page.evaluate(focused) == "frequency"
    # The arrows wrap, so ArrowUp from the first item lands on the last.
    page.keyboard.press("ArrowUp")
    assert page.evaluate(focused) == "summary"
    page.keyboard.press("Home")
    assert page.evaluate(focused) == "frequency"
    page.keyboard.press("End")
    assert page.evaluate(focused) == "summary"

    page.keyboard.press("Enter")
    page.wait_for_selector("#column-chart-overlay:not([hidden])")
    assert page.evaluate("() => state.columnChart.kind") == "summary"
    assert page.pageerrors == []


def test_closing_the_chart_returns_focus_to_its_glyph(meta_page):
    page = meta_page
    _open_column_chart(page, "score", "histogram")
    page.keyboard.press("Escape")
    page.wait_for_selector("#column-chart-overlay", state="hidden")

    assert page.evaluate(
        "() => document.activeElement.dataset.chartColumn"
    ) == "score"
    assert page.pageerrors == []


def test_preset_digits_do_not_fire_while_a_chart_is_open(meta_page):
    """A digit shortcut would change the selection out from under the chart."""
    page = meta_page
    page.evaluate("() => { state.selectedNodeIndices = new Set([0, 1]); updateMetadataTable(); }")
    page.keyboard.press("Shift+Digit1")
    page.evaluate("() => { state.selectedNodeIndices = new Set(); updateMetadataTable(); }")
    _open_column_chart(page, "family", "frequency")

    page.keyboard.press("Digit1")
    assert page.evaluate("() => state.selectedNodeIndices.size") == 0

    # ...and it works again once the dialog is closed.
    page.keyboard.press("Escape")
    page.wait_for_selector("#column-chart-overlay", state="hidden")
    page.keyboard.press("Digit1")
    assert page.evaluate("() => state.selectedNodeIndices.size") == 2
    assert page.pageerrors == []


# ---------------------------------------------------------------------------
# "Jump to threshold" arrows
# ---------------------------------------------------------------------------


def _finite_stop_values(page):
    return page.evaluate(
        "() => state.sliderModel.stops.filter(stop => stop.threshold_value !== null)"
        ".map(stop => stop.threshold_value)"
    )


def _snap_to_stop(page, stop_index):
    page.evaluate(
        "i => { const stops = state.sliderModel.stops;"
        " snapSliderToStop(stops[i < 0 ? stops.length + i : i]); updateThresholdUI(false); }",
        stop_index,
    )


def _wait_for_threshold_field(page, value):
    page.wait_for_function(
        "v => document.getElementById('threshold-input').value === String(v)", arg=value
    )


def test_threshold_arrows_step_between_split_stops(meta_page):
    """The arrows walk the slider's stop list, not a fixed numeric step.

    A number input's spinner stepped by 1.0, which on a distance axis is either a
    no-op or a leap over several splits depending on the metric.
    """
    page = meta_page
    finite = _finite_stop_values(page)
    assert len(finite) >= 3
    # The premise: consecutive splits are not one unit apart, so a constant step
    # could not land on them.
    assert finite[1] - finite[0] != 1

    _snap_to_stop(page, 0)
    _wait_for_threshold_field(page, finite[0])
    assert page.is_disabled("#threshold-step-down")
    assert page.is_enabled("#threshold-step-up")

    page.click("#threshold-step-up")
    _wait_for_threshold_field(page, finite[1])
    assert page.evaluate("() => selectedThresholdValue()") == finite[1]

    page.click("#threshold-step-up")
    _wait_for_threshold_field(page, finite[2])

    page.click("#threshold-step-down")
    _wait_for_threshold_field(page, finite[1])
    assert page.evaluate("() => selectedThresholdValue()") == finite[1]
    assert page.pageerrors == []


def test_threshold_arrow_keys_step_stops_too(meta_page):
    """ArrowUp/ArrowDown in the jump field are the keyboard form of the buttons."""
    page = meta_page
    finite = _finite_stop_values(page)
    _snap_to_stop(page, 0)
    _wait_for_threshold_field(page, finite[0])

    page.click("#threshold-input")
    page.keyboard.press("ArrowUp")
    _wait_for_threshold_field(page, finite[1])
    page.keyboard.press("ArrowDown")
    _wait_for_threshold_field(page, finite[0])
    assert page.pageerrors == []


def test_threshold_arrows_stop_at_the_ends_of_the_stop_list(meta_page):
    """The last stop is the infinity stop; there is nothing past either end."""
    page = meta_page
    _snap_to_stop(page, -1)
    page.wait_for_selector("#threshold-step-up[disabled]")
    assert page.is_enabled("#threshold-step-down")
    assert page.evaluate("() => selectedThresholdValue() === Infinity") is True
    # Calling past the end is a no-op rather than an error or a wrap-around.
    assert page.evaluate(
        "() => { stepThreshold(1); return selectedThresholdValue() === Infinity; }") is True

    _snap_to_stop(page, 0)
    page.wait_for_selector("#threshold-step-down[disabled]")
    lowest = page.evaluate("() => selectedThresholdValue()")
    assert page.evaluate(
        "() => { stepThreshold(-1); return selectedThresholdValue(); }") == lowest
    assert page.pageerrors == []


def test_typing_a_threshold_still_snaps_to_the_nearest_stop(meta_page):
    """The field is type="text" now, so its own parsing has to keep working."""
    page = meta_page
    finite = _finite_stop_values(page)
    page.fill("#threshold-input", str(finite[-1] + 100))
    page.press("#threshold-input", "Enter")
    _wait_for_threshold_field(page, finite[-1])
    assert page.evaluate("() => selectedThresholdValue()") == finite[-1]
    assert page.pageerrors == []


# ---------------------------------------------------------------------------
# Select by value
# ---------------------------------------------------------------------------


def _selection(page):
    return page.evaluate("() => Array.from(state.selectedNodeIndices).sort((a, b) => a - b)")


def _select_by_value(page, column, op, value, second=None, action="add"):
    """Fill in the Select-by-value panel and click one of its three buttons."""
    _open_edit_panel(page, "select")
    page.select_option("#metadata-select-column", column)
    page.select_option("#metadata-select-op", op)
    page.fill("#metadata-select-value", value)
    if second is not None:
        page.fill("#metadata-select-value2", second)
    page.click(f"#metadata-select-{action}")


def test_select_by_value_contains_searches_the_whole_network(meta_page):
    page = meta_page
    _select_by_value(page, "family", "contains", "alph")
    assert _selection(page) == [0, 1]
    assert "Matched 2 of 6 nodes in the network" in page.text_content("#metadata-select-note")
    assert page.pageerrors == []


def test_select_by_value_exact_is_case_insensitive(meta_page):
    page = meta_page
    _select_by_value(page, "family", "exact", "BETA")
    assert _selection(page) == [2, 3]
    # ...and it is exact: a prefix of a real value matches nothing.
    page.evaluate("() => { state.selectedNodeIndices = new Set(); updateMetadataTable(); }")
    _select_by_value(page, "family", "exact", "bet")
    assert _selection(page) == []
    assert page.pageerrors == []


def test_select_by_value_exact_matches_a_number_however_it_is_typed(meta_page):
    """Node C's score is 5; "5.0" is the same number and has to match it."""
    page = meta_page
    for typed in ("5", "5.0", "05"):
        page.evaluate("() => { state.selectedNodeIndices = new Set(); updateMetadataTable(); }")
        _select_by_value(page, "score", "exact", typed)
        assert _selection(page) == [2], typed
    assert page.pageerrors == []


def test_select_by_value_compares_numeric_columns_as_numbers(meta_page):
    """The scores are 1, 2, 5, 8, 3, 9: "less than 10" is all of them.

    Text order would put "10" between "1" and "2" and match only score 1, so this
    is the case that distinguishes numeric from lexicographic comparison.
    """
    page = meta_page
    _select_by_value(page, "score", "lt", "10")
    assert _selection(page) == [0, 1, 2, 3, 4, 5]

    page.evaluate("() => { state.selectedNodeIndices = new Set(); updateMetadataTable(); }")
    _select_by_value(page, "score", "gt", "5")
    assert _selection(page) == [3, 5]
    assert page.pageerrors == []


def test_select_by_value_between_is_inclusive_and_order_tolerant(meta_page):
    page = meta_page
    _select_by_value(page, "score", "between", "2", second="5")
    assert _selection(page) == [1, 2, 4]

    # The same query with the bounds typed the wrong way round.
    page.evaluate("() => { state.selectedNodeIndices = new Set(); updateMetadataTable(); }")
    _select_by_value(page, "score", "between", "5", second="2")
    assert _selection(page) == [1, 2, 4]
    assert page.pageerrors == []


def test_select_by_value_between_asks_for_both_bounds(meta_page):
    page = meta_page
    _select_by_value(page, "score", "between", "2")
    assert _selection(page) == []
    assert page.text_content("#metadata-select-note") == (
        "Type both bounds, or pick another comparison."
    )
    # The second field only appears for a two-operand comparison.
    assert page.is_visible("#metadata-select-value2")
    page.select_option("#metadata-select-op", "contains")
    assert page.is_hidden("#metadata-select-value2")
    assert page.pageerrors == []


def test_select_by_value_regex_and_its_error_path(meta_page):
    page = meta_page
    _select_by_value(page, "family", "regex", "^g")
    assert _selection(page) == [4, 5]

    # A malformed pattern reports itself and leaves the selection alone.
    _select_by_value(page, "family", "regex", "[")
    assert _selection(page) == [4, 5]
    assert page.text_content("#metadata-select-note").startswith(
        "Not a valid regular expression"
    )
    assert page.pageerrors == []


def test_select_by_value_can_query_node_id(meta_page):
    """node_id is offered alongside the metadata columns."""
    page = meta_page
    _select_by_value(page, "__node_id__", "exact", "c")
    assert _selection(page) == [2]
    assert page.pageerrors == []


def test_select_by_value_remove_and_subset_never_grow_the_selection(meta_page):
    page = meta_page
    _select_nodes(page, [0, 1, 2, 3])

    # F (index 5) matches "score greater than 4" but is not selected, so neither
    # narrowing action may pull it in.
    _select_by_value(page, "score", "gt", "4", action="subset")
    assert _selection(page) == [2, 3]
    assert "in the selection" in page.text_content("#metadata-select-note")

    _select_nodes(page, [0, 1, 2, 3])
    _select_by_value(page, "family", "contains", "beta", action="remove")
    assert _selection(page) == [0, 1]
    assert page.pageerrors == []


def test_select_by_value_narrowing_needs_a_selection(meta_page):
    page = meta_page
    _open_edit_panel(page, "select")
    assert page.is_enabled("#metadata-select-add")
    assert page.is_disabled("#metadata-select-remove")
    assert page.is_disabled("#metadata-select-subset")

    _select_nodes(page, [0, 1])
    assert page.is_enabled("#metadata-select-remove")
    assert page.is_enabled("#metadata-select-subset")
    assert page.pageerrors == []


def test_select_by_value_ignores_blank_cells(meta_page):
    """An empty cell is an absence, so no comparison may match it."""
    page = meta_page
    _add_column(page, "empty")
    _select_by_value(page, "empty", "contains", "a")
    assert _selection(page) == []
    assert "Matched 0 of 6 nodes" in page.text_content("#metadata-select-note")

    # ...including an ordering comparison, which would otherwise treat null as 0.
    _select_by_value(page, "empty", "lt", "1000")
    assert _selection(page) == []
    assert page.pageerrors == []


def test_select_by_value_enter_adds_to_the_selection(meta_page):
    page = meta_page
    _open_edit_panel(page, "select")
    page.select_option("#metadata-select-column", "family")
    page.select_option("#metadata-select-op", "contains")
    page.fill("#metadata-select-value", "gamma")
    page.press("#metadata-select-value", "Enter")
    assert _selection(page) == [4, 5]
    assert page.pageerrors == []


def test_select_by_value_menu_follows_a_column_rename(meta_page):
    page = meta_page
    _open_edit_panel(page, "select")
    page.select_option("#metadata-select-column", "family")

    _open_edit_panel(page, "rename")
    page.select_option("#metadata-rename-column", "family")
    page.fill("#metadata-rename-value", "clan")
    page.click("#metadata-rename-apply")

    _open_edit_panel(page, "select")
    assert page.input_value("#metadata-select-column") == "clan"
    _select_by_value(page, "clan", "contains", "alpha")
    assert _selection(page) == [0, 1]
    assert page.pageerrors == []


# ---------------------------------------------------------------------------
# Default categorical palette
# ---------------------------------------------------------------------------


def test_an_untouched_column_is_already_on_the_default_palette(meta_page):
    """No stored palette, yet every value has a get_palette color."""
    from domainator.utils import NAMED_CATEGORICAL_PALETTES

    page = meta_page
    distinct_colors = next(
        palette["colors"] for palette in NAMED_CATEGORICAL_PALETTES
        if palette["name"] == "domainator"
    )
    assert page.evaluate("() => Object.keys(state.customPalettes)") == []
    # alpha, beta, gamma in sorted value order take the first three colors.
    assert page.evaluate("() => state.nodeColorCache.slice()") == [
        distinct_colors[0], distinct_colors[0],
        distinct_colors[1], distinct_colors[1],
        distinct_colors[2], distinct_colors[2],
    ]
    assert page.pageerrors == []


def test_a_hand_edited_swatch_keeps_the_defaults_of_the_other_values(meta_page):
    """Storing one swatch must not shadow the default palette for the rest.

    The stored palette is what customPalette() returns once it exists, so it has
    to be seeded with the default assignment rather than start empty.
    """
    page = meta_page
    before = page.evaluate("() => state.nodeColorCache.slice()")
    _open_color_picker(page)
    page.wait_for_selector("#color-picker-discrete:not([hidden])")

    swatch = page.locator("#color-picker-swatch-list input[type=color]").first
    swatch.fill("#123456")
    swatch.dispatch_event("input")
    page.wait_for_function("() => state.nodeColorCache[0] === '#123456'")

    after = page.evaluate("() => state.nodeColorCache.slice()")
    assert after[0] == after[1] == "#123456"
    assert after[2:] == before[2:]
    assert page.pageerrors == []


def test_reset_to_defaults_returns_a_column_to_the_default_palette(meta_page):
    page = meta_page
    before = page.evaluate("() => state.nodeColorCache.slice()")
    _open_color_picker(page)
    page.wait_for_selector("#color-picker-discrete:not([hidden])")
    page.select_option("#color-palette", "okabe_ito")
    page.wait_for_function("() => state.customPalettes.family !== undefined")
    assert page.evaluate("() => state.nodeColorCache.slice()") != before

    page.click("#color-picker-reset")
    page.wait_for_selector("#color-picker-discrete:not([hidden])")
    assert page.evaluate("() => Object.keys(state.customPalettes)") == []
    assert page.evaluate("() => state.nodeColorCache.slice()") == before
    assert page.input_value("#color-palette") == "domainator"
    assert page.pageerrors == []


def test_default_palette_cache_is_rebuilt_after_a_metadata_edit(meta_page):
    """The default assignment depends on the whole value set, so it must not stick."""
    page = meta_page
    # Already populated by the first paint: every node's color goes through it.
    assert page.evaluate("() => state.defaultPalettes.size") == 1
    alpha_before = page.evaluate("() => customPalette('family').colors.alpha")

    page.evaluate(
        "() => { setMetadataValue(0, 'family', 'aardvark'); refreshAfterMetadataEdit(); }")

    colors = page.evaluate("() => customPalette('family').colors")
    assert set(colors) == {"aardvark", "alpha", "beta", "gamma"}
    # Colors are handed out in sorted value order, and aardvark sorts first, so it
    # takes the color alpha had and pushes every other family along by one. A stale
    # cache would have left alpha where it was and given aardvark nothing at all.
    assert colors["aardvark"] == alpha_before
    assert colors["alpha"] != alpha_before
    assert page.evaluate("() => nodeColor(0)") == colors["aardvark"]
    assert page.evaluate("() => nodeColor(1)") == colors["alpha"]
    assert page.pageerrors == []


def test_no_hidden_element_is_actually_rendered(meta_page):
    """A `display` rule on a class overrides the `hidden` attribute silently.

    The viewer hides most of its panels with `hidden` and styles them with
    `display: flex`, which wins -- so every such rule needs a `[hidden]` guard,
    and forgetting one leaves a dead control on screen with nothing to say. This
    sweeps the whole page instead of trusting each rule to be remembered.
    """
    page = meta_page
    visible_but_hidden = "() => Array.from(document.querySelectorAll('[hidden]'))" \
        ".filter(el => getComputedStyle(el).display !== 'none')" \
        ".map(el => el.id || el.className)"
    assert page.evaluate(visible_but_hidden) == []

    _open_color_picker(page)
    page.wait_for_selector("#color-picker-discrete:not([hidden])")
    assert page.evaluate(visible_but_hidden) == []
    page.click("#color-picker-close")

    for panel in ("add", "set", "rename", "delete", "select"):
        _open_edit_panel(page, panel)
        assert page.evaluate(visible_but_hidden) == [], panel
    assert page.pageerrors == []


# ---------------------------------------------------------------------------
# "N merge events plotted": the split chart is a capped selection
# ---------------------------------------------------------------------------


def _split_event_caption(page):
    return page.text_content("#split-event-count")


def test_split_event_caption_says_the_whole_series_is_drawn(page):
    """Small networks are under the cap, and the caption has to say so plainly."""
    series = page.evaluate("() => ({plotted: state.series.selectedRows.length,"
                           " total: state.series.total})")
    assert series["plotted"] == series["total"]
    assert _split_event_caption(page) == (
        f"All {series['total']:,} merge events plotted."
    )
    assert "Every split event" in page.get_attribute("#split-event-count", "title")
    assert page.pageerrors == []


def test_split_event_caption_reports_a_capped_series(capped_page):
    """The number that makes a blank stretch of axis readable as "filtered"."""
    page = capped_page
    series = page.evaluate("() => ({plotted: state.series.selectedRows.length,"
                           " total: state.series.total, cap: state.series.cap})")
    plotted, total, cap = series["plotted"], series["total"], series["cap"]
    assert cap == CAPPED_PAGE_MAX_MERGE_EVENTS == 5
    assert plotted < total
    # The back-fill: between the cap and the cap plus one per 5% band.
    assert cap <= plotted <= cap + 20

    caption = _split_event_caption(page)
    assert caption.startswith(f"{plotted:,} of {total:,} merge events plotted")
    assert f"the strongest {cap:,} by impact" in caption
    assert f"plus {plotted - cap:,}" in caption
    assert "every 5% of the axis" in caption
    assert '"Split events" box' in page.get_attribute("#split-event-count", "title")
    assert page.pageerrors == []


def test_capped_series_still_reaches_the_bottom_of_the_axis(capped_page):
    """What the back-fill is for: ranking by impact alone strands the left end."""
    page = capped_page
    weakest_edge = page.evaluate(
        "() => Math.min(...state.bundle.graph.mst_edges.map(e => e[2]))")
    strongest_edge = page.evaluate(
        "() => Math.max(...state.bundle.graph.mst_edges.map(e => e[2]))")
    weakest_plotted = page.evaluate(
        "() => Math.min(...state.series.selectedRows.map(e => e.threshold_value))")
    span = strongest_edge - weakest_edge
    assert (weakest_plotted - weakest_edge) < 0.1 * span

    # And the slider gets stops across the whole track rather than only the right end.
    positions = page.evaluate(
        "() => state.sliderModel.stops.map(stop => stop.sliderPosition)")
    assert min(p for p in positions if p > 0) < 200
    assert page.pageerrors == []


def test_split_event_caption_ignores_a_stale_stored_series(page):
    """A pre-v6 bundle carries its own capped copy of the series. The viewer derives
    its own, so a stale stored one can neither shrink the chart nor skew the caption."""
    page.evaluate("""() => {
        state.bundle.graph.merge_event_series = [];
        state.bundle.graph.merge_event_total = 999999;
        state.bundle.graph.max_merge_events = 1;
        state.bundle.graph.merge_moving_sum = {window: 0, x: [], y: []};
        state.bundle.graph.slider_stops = [];
        updateSplitEventCount();
    }""")
    total = page.evaluate("() => state.series.total")
    assert total > 0
    assert _split_event_caption(page) == f"All {total:,} merge events plotted."
    assert page.pageerrors == []


def test_js_merge_event_filter_matches_python_including_the_backfill(page):
    """The JS port is what viewer-built extractions use; it must not drift.

    Both implementations are handed the same rows -- big events crowded into the
    top of the axis, tiny ones spread below -- which is the distribution that makes
    the cap and the back-fill disagree.
    """
    from domainator.ssn_hierarchy import MERGE_EVENT_DENSITY_BINS, filter_merge_event_rows

    rows = []
    for index in range(200):
        rows.append({"edge_index": index, "threshold_value": 0.95 + (0.05 * index / 200),
                     "merge_impact": 100.0 + index, "delta_largest": 0.0,
                     "delta_avg_non_singleton": 0.0})
    for index in range(200, 300):
        rows.append({"edge_index": index, "threshold_value": 0.95 * (index - 200) / 100.0,
                     "merge_impact": 1.0, "delta_largest": 0.0,
                     "delta_avg_non_singleton": 0.0})

    for cap in (1, 5, 25, 200):
        got = page.evaluate(
            "([rows, cap]) => filterMergeEventRows(rows, cap)"
            ".map(row => row.edge_index)",
            [rows, cap],
        )
        want = [row["edge_index"] for row in filter_merge_event_rows(rows, max_merge_events=cap)]
        assert got == want, cap
        assert cap <= len(got) <= cap + MERGE_EVENT_DENSITY_BINS

    # The knob itself ports too.
    assert page.evaluate(
        "rows => filterMergeEventRows(rows, 5, 0).length", rows) == 5
    assert page.pageerrors == []


def test_split_event_cap_control_rederives_the_series(dense_page):
    """The cap is a control, not something the bundle was built with: moving it
    re-derives the series and re-lays-out the slider, off the same file."""
    page = dense_page
    total = page.evaluate("() => state.series.total")
    assert total > 25, "fixture needs more events than the caps under test"

    _set_split_event_cap(page, 5)
    capped = page.evaluate("() => state.series.selectedRows.length")
    capped_stops = page.evaluate("() => state.sliderModel.stops.length")
    assert 5 <= capped <= 5 + 20
    assert capped < total

    _set_split_event_cap(page, 0)
    assert page.evaluate("() => state.series.selectedRows.length") == total
    # Every event gets a stop, plus the floor; the infinity stop is in the model too.
    assert page.evaluate("() => state.sliderModel.stops.length") > capped_stops
    assert _split_event_caption(page) == f"All {total:,} merge events plotted."
    assert page.pageerrors == []


def test_zooming_reveals_more_split_events_in_the_window(crowded_viewer_html):
    """The point of the windowed cap: the same small number buys more detail as you
    zoom, because it is spent on the events actually on screen."""
    for page in _yield_loaded_page(crowded_viewer_html):
        _set_split_event_cap(page, 10)

        def in_window(window):
            return page.evaluate(
                """window => state.series.selectedRows.filter(
                       row => row.threshold_value >= window.min
                           && row.threshold_value <= window.max).length""",
                window)

        x, y = _split_chart_fraction_point(page, fraction=0.85)
        page.mouse.move(x, y)
        page.mouse.wheel(0, -1500)
        page.wait_for_function("() => Boolean(state.splitChartZoom)")
        window = page.evaluate("() => state.splitChartZoom")
        zoomed = in_window(window)

        page.click("#split-chart-reset-zoom")
        page.wait_for_function("() => state.splitChartZoom === null")
        flat = in_window(window)

        # Zoomed, the cap is spent inside the window; flat, that same stretch competes
        # with the whole axis for the same ten slots.
        assert zoomed > flat
        assert zoomed <= 10
        assert page.pageerrors == []


def test_zooming_keeps_whole_range_coverage_and_the_current_cut(crowded_viewer_html):
    """Zooming must not strand the rest of the axis, nor move the cut in effect.

    The band back-fill runs over the whole range however narrow the window is, so the
    slider still spans everything; and the selected threshold is pinned through the
    re-selection, so a zoom can never re-cluster the network under the user.
    """
    for page in _yield_loaded_page(crowded_viewer_html):
        _set_split_event_cap(page, 10)
        # Sit on a stop well away from where we are about to zoom.
        page.evaluate("""() => {
            const finite = state.sliderModel.stops.filter(s => s.threshold_value !== null);
            snapSliderToStop(finite[Math.floor(finite.length / 2)]);
            updateThresholdUI();
        }""")
        before = page.evaluate("() => selectedThresholdValue()")
        clusters_before = page.evaluate("() => activeClustersAtThreshold(selectedThresholdValue()).length")

        x, y = _split_chart_fraction_point(page, fraction=0.05)
        page.mouse.move(x, y)
        page.mouse.wheel(0, -2000)
        page.wait_for_function("() => Boolean(state.splitChartZoom)")

        # The cut is untouched, so the clustering on the canvas is untouched.
        assert page.evaluate("() => selectedThresholdValue()") == pytest.approx(before)
        assert page.evaluate(
            "() => activeClustersAtThreshold(selectedThresholdValue()).length") == clusters_before
        # And the selected stop is still in the list, not merely nearest to it.
        assert page.evaluate(
            "value => state.sliderModel.stops.some(s => s.threshold_value === value)", before)

        # Stops still reach across the whole axis, not just the zoomed window.
        coverage = page.evaluate("""() => {
            const finite = state.sliderModel.stops.filter(s => s.threshold_value !== null);
            const all = state.series.eventRows.map(r => r.threshold_value);
            const lo = Math.min(...all), hi = Math.max(...all);
            // The production bin, clamp included: the topmost event lands in the last
            // band rather than one past the end.
            const bin = value => mergeEventDensityBin(value, lo, hi, 20);
            return {
                stopBins: [...new Set(finite.map(s => bin(s.threshold_value)))].sort((a, b) => a - b),
                eventBins: [...new Set(all.map(bin))].sort((a, b) => a - b),
                windowSpan: state.splitChartZoom.max - state.splitChartZoom.min,
                dataSpan: hi - lo,
            };
        }""")
        # Every 5% band that has an event to offer still has a stop, however narrow the
        # window is -- that is what the back-fill running over the whole range buys. The
        # band's representative is its strongest event, not its highest-scoring one, so
        # this is coverage of the axis rather than a claim about particular thresholds.
        assert set(coverage["eventBins"]) <= set(coverage["stopBins"])
        assert len(coverage["eventBins"]) > 1, "fixture needs events in more than one band"
        # ...and the window really was narrow, or the coverage claim proves nothing.
        assert coverage["windowSpan"] < 0.5 * coverage["dataSpan"]
        assert page.pageerrors == []


def test_split_event_cap_keeps_the_threshold_you_are_looking_at(dense_page):
    """Raising the cap re-lays-out the track, so the *value* has to be what survives."""
    page = dense_page
    _set_split_event_cap(page, 5)
    page.evaluate("""() => {
        const stop = state.sliderModel.stops.filter(s => s.threshold_value !== null)[1];
        snapSliderToStop(stop);
        updateThresholdUI();
    }""")
    before = page.evaluate("() => selectedThresholdValue()")

    _set_split_event_cap(page, 0)
    after = page.evaluate("() => selectedThresholdValue()")
    # The old stop is still a stop at a larger cap, so the threshold is unchanged.
    assert after == pytest.approx(before)
    assert page.pageerrors == []


def test_split_event_cap_is_saved_with_a_session(dense_page, tmp_path):
    """A session records the cap, so reopening it shows the chart that was saved."""
    from domainator.ssn_bundle import SSN_VIEWER_BUNDLE_VERSION, load_bundle

    page = dense_page
    _set_split_event_cap(page, 7)
    saved = _save_session(page, tmp_path, "capped_session.dsnv")
    bundle = load_bundle(saved)

    assert bundle["version"] == SSN_VIEWER_BUNDLE_VERSION == 6
    assert bundle["app_state"]["view"]["max_merge_events"] == 7
    # Saving re-serializes the loaded graph, which carries nothing derived.
    for derived_key in ("merge_event_series", "merge_event_total", "max_merge_events",
                        "merge_moving_sum", "slider_stops",
                        "cluster_count_by_threshold", "edges_by_threshold"):
        assert derived_key not in bundle["graph"]

    _load_bundle_file(page, saved)
    assert page.evaluate("() => state.maxMergeEvents") == 7
    assert page.evaluate("() => state.series.cap") == 7
    assert page.pageerrors == []


# ---------------------------------------------------------------------------
# Split chart hover readout and click-to-jump
# ---------------------------------------------------------------------------


def _split_chart_scale(page):
    rect = page.eval_on_selector(
        "#split-chart",
        "c => { const r = c.getBoundingClientRect();"
        " return {left: r.left, top: r.top, width: r.width}; }",
    )
    return rect, rect["width"] / 1100


def _split_chart_point(page, x, y):
    """Client coordinates for a point given in split-canvas coordinates."""
    rect, scale = _split_chart_scale(page)
    return rect["left"] + (x * scale), rect["top"] + (y * scale)


def _hover_split_chart(page, x, y):
    page.mouse.move(*_split_chart_point(page, x, y))


def _click_split_chart(page, x, y):
    page.mouse.click(*_split_chart_point(page, x, y))


def _split_tip_lines(page):
    return page.eval_on_selector_all(
        "#split-chart-tip > div", "els => els.map(el => el.textContent)")


def _split_hit_events(page):
    return page.evaluate("() => state.splitChartHit.events")


def test_split_chart_records_the_geometry_it_painted(page):
    """Hit-testing reads what the draw recorded, so the two cannot disagree."""
    hit = page.evaluate("() => state.splitChartHit")
    events = page.evaluate("() => state.series.selectedRows")

    assert len(hit["events"]) == len(events)
    for entry, event in zip(hit["events"], events):
        assert entry["event"]["edge_index"] == event["edge_index"]
        # One bead per distinct merge size, exactly as the chart draws.
        assert len(entry["beads"]) == len(event["merge_size_counts"])
        assert hit["plotLeft"] <= entry["x"] <= hit["plotLeft"] + hit["plotWidth"]
        for bead in entry["beads"]:
            assert hit["plotTop"] - 1 <= bead["y"] <= hit["plotTop"] + hit["plotHeight"] + 1
    assert page.pageerrors == []


def test_hovering_a_split_event_reports_the_same_numbers_as_the_report(page):
    """The readout carries matrix_report's hovertemplate fields."""
    entry = max(_split_hit_events(page), key=lambda e: e["event"]["largest_merge"])
    bead = entry["beads"][0]
    _hover_split_chart(page, entry["x"], bead["y"])
    page.wait_for_selector("#split-chart-tip:not([hidden])")

    lines = _split_tip_lines(page)
    event = entry["event"]
    assert lines[0].startswith(f"Threshold {event['threshold_to']}")
    assert lines[1] == (
        f"{bead['count']:,} split{'' if bead['count'] == 1 else 's'} of "
        f"{bead['size']:,} node{'' if bead['size'] == 1 else 's'}"
    )
    assert lines[2] == (
        f"{event['merge_count']:,} split{'' if event['merge_count'] == 1 else 's'} here, "
        f"{int(event['merge_impact']):,} node"
        f"{'' if event['merge_impact'] == 1 else 's'} total"
    )
    # The moving-sum label is metric-aware, so it comes from the same helper
    # matrix_report's hovertemplate reads.
    assert lines[3].startswith("Nodes displaced within \u00b1")
    assert lines[-1] == "Click to jump here"
    assert page.eval_on_selector("#split-chart", "c => c.style.cursor") == "pointer"
    assert page.pageerrors == []


def test_hovering_picks_the_bead_nearest_the_cursor(two_bead_page):
    """One stem, two merge sizes: the readout has to name the one under the pointer."""
    page = two_bead_page
    entry = next(e for e in _split_hit_events(page) if len(e["beads"]) == 2)
    beads = sorted(entry["beads"], key=lambda b: b["size"])
    assert [b["size"] for b in beads] == [1, 2]

    for bead in beads:
        _hover_split_chart(page, entry["x"], bead["y"])
        page.wait_for_selector("#split-chart-tip:not([hidden])")
        assert _split_tip_lines(page)[1] == (
            f"1 split of {bead['size']} node{'' if bead['size'] == 1 else 's'}")
    assert page.pageerrors == []


def test_hovering_between_beads_falls_back_to_the_largest_split(two_bead_page):
    """Away from every bead the stem still means something: its own height."""
    page = two_bead_page
    entry = next(e for e in _split_hit_events(page) if len(e["beads"]) == 2)
    first, second = (bead["y"] for bead in entry["beads"])
    # Halfway between them, which on this network is ~58px from each -- well past
    # the bead radius, so the fallback is genuinely what is being exercised.
    assert abs(first - second) > 40
    _hover_split_chart(page, entry["x"], (first + second) / 2)
    page.wait_for_selector("#split-chart-tip:not([hidden])")

    largest = int(entry["event"]["largest_merge"])
    assert _split_tip_lines(page)[1] == f"Largest single split: {largest} nodes"
    assert page.pageerrors == []


def test_hovering_the_moving_sum_line_reads_the_window(page):
    """Off every event, the readout is the moving sum at that threshold."""
    hit = page.evaluate("() => state.splitChartHit")
    xs = sorted(entry["x"] for entry in hit["events"])
    # The middle of the widest gap between events is as far from a stem as it gets.
    gap_x, widest = xs[0], 0.0
    for left, right in zip(xs, xs[1:]):
        if right - left > widest:
            widest, gap_x = right - left, (left + right) / 2
    assert widest > 40, "fixture has no gap wide enough to land between events"

    _hover_split_chart(page, gap_x, hit["plotTop"] + (hit["plotHeight"] / 2))
    page.wait_for_selector("#split-chart-tip:not([hidden])")
    lines = _split_tip_lines(page)
    assert lines[0].startswith("Threshold ")
    assert lines[1].startswith("Nodes displaced within \u00b1")
    assert lines[-1] == "Click to jump to the nearest split"

    # The value is the step function's, read at the cursor's threshold.
    threshold = page.evaluate("x => splitChartThresholdAt(x)", gap_x)
    expected = page.evaluate("t => splitChartMovingSumAt(t)", threshold)
    assert f"{int(expected):,}" in lines[1]
    assert page.pageerrors == []


def test_moving_sum_lookup_matches_a_linear_scan(page):
    """The binary search has to agree with the step function actually drawn.

    Checked against a scan rather than against itself, at sample points and at the
    exact sample thresholds where the step jumps.
    """
    xs = page.evaluate("() => state.series.movingSum.x")
    ys = page.evaluate("() => state.series.movingSum.y")
    assert len(xs) > 100

    probes = [xs[0], xs[1], xs[len(xs) // 2], xs[-1]]
    for left, right in zip(xs, xs[1:]):
        probes.append((left + right) / 2)
        if len(probes) > 40:
            break
    for threshold in probes:
        expected = ys[max(i for i, x in enumerate(xs) if x <= threshold)]
        assert page.evaluate("t => splitChartMovingSumAt(t)", threshold) == expected, threshold

    # Below the first sample there is no line to read.
    assert page.evaluate("t => splitChartMovingSumAt(t)", xs[0] - 1) is None
    assert page.pageerrors == []


def test_hover_readout_hides_when_the_pointer_leaves(page):
    entry = _split_hit_events(page)[0]
    _hover_split_chart(page, entry["x"], entry["beads"][0]["y"])
    page.wait_for_selector("#split-chart-tip:not([hidden])")

    page.mouse.move(5, 5)
    page.wait_for_selector("#split-chart-tip", state="hidden")
    assert page.eval_on_selector("#split-chart", "c => c.style.cursor") == ""
    assert page.pageerrors == []


def test_hover_is_inert_outside_the_plot_area(page):
    """The axes, the titles and the margins are not data."""
    hit = page.evaluate("() => state.splitChartHit")
    for x, y in [(4, 4), (hit["plotLeft"] - 30, hit["plotTop"] + 10),
                 (hit["plotLeft"] + 10, hit["plotTop"] + hit["plotHeight"] + 40)]:
        _hover_split_chart(page, x, y)
        assert page.eval_on_selector("#split-chart-tip", "e => e.hidden") is True, (x, y)
    assert page.pageerrors == []


def test_clicking_a_split_event_jumps_to_its_threshold(page):
    """Every plotted event is also a slider stop, so the landing is exact."""
    entry = max(_split_hit_events(page), key=lambda e: e["event"]["largest_merge"])
    target = entry["event"]["threshold_value"]
    page.evaluate("() => { snapSliderToStop(state.sliderModel.stops[0]); updateThresholdUI(false); }")
    assert page.evaluate("() => selectedThresholdValue()") != target

    _click_split_chart(page, entry["x"], entry["beads"][0]["y"])
    page.wait_for_function("v => selectedThresholdValue() === v", arg=target)
    # Waited for rather than asserted: the click snaps the slider immediately but
    # repaints the readouts on the next animation frame, so asserting here would be
    # asserting against whichever of the two won the race.
    page.wait_for_function(
        "v => document.getElementById('threshold-input').value === String(v)", arg=target)
    assert page.pageerrors == []


def test_clicking_the_moving_sum_line_jumps_to_the_nearest_stop(page):
    """A click out on the line has no event of its own, so it snaps."""
    hit = page.evaluate("() => state.splitChartHit")
    xs = sorted(entry["x"] for entry in hit["events"])
    gap_x, widest = xs[0], 0.0
    for left, right in zip(xs, xs[1:]):
        if right - left > widest:
            widest, gap_x = right - left, (left + right) / 2

    threshold = page.evaluate("x => splitChartThresholdAt(x)", gap_x)
    expected = page.evaluate("t => nearestStopForThreshold(t).threshold_value", threshold)
    _click_split_chart(page, gap_x, hit["plotTop"] + (hit["plotHeight"] / 2))
    page.wait_for_function("v => selectedThresholdValue() === v", arg=expected)
    assert page.pageerrors == []


def test_clicking_the_split_chart_keeps_the_pan_and_zoom(page):
    """Stepping along the chart to watch one region break up is the point."""
    page.evaluate("() => { state.viewTransform.scale = 2.5;"
                  " state.viewTransform.offsetX = 40; state.viewTransform.offsetY = -25; }")
    entry = max(_split_hit_events(page), key=lambda e: e["event"]["largest_merge"])
    _click_split_chart(page, entry["x"], entry["beads"][0]["y"])
    page.wait_for_function("v => selectedThresholdValue() === v",
                           arg=entry["event"]["threshold_value"])

    assert page.evaluate("() => state.viewTransform.scale") == 2.5
    assert page.evaluate("() => state.viewTransform.offsetX") == 40
    assert page.evaluate("() => state.viewTransform.offsetY") == -25
    assert page.pageerrors == []


def test_hover_readout_does_not_call_a_product_impact_a_node_count(product_metric_page):
    """min_child impacts are nodes; product impacts are not, and must not say so."""
    page = product_metric_page
    assert page.evaluate("() => state.bundle.graph.merge_impact_metric") == "product"
    entry = max(_split_hit_events(page), key=lambda e: e["event"]["largest_merge"])
    _hover_split_chart(page, entry["x"], entry["beads"][0]["y"])
    page.wait_for_selector("#split-chart-tip:not([hidden])")

    lines = _split_tip_lines(page)
    assert "node" not in " ".join(lines)
    assert lines[1].startswith("1 split of ") or " splits of " in lines[1]
    # The axis title makes the same distinction.
    assert page.evaluate("() => splitAxisLabels().largest") == (
        "Largest single split (size product)")
    assert page.pageerrors == []


def test_split_chart_hover_survives_a_threshold_change(page):
    """Every redraw re-records the geometry; a stale record would mis-hit."""
    before = page.evaluate("() => state.splitChartHit.events.map(e => e.x)")
    page.click("#threshold-step-up")
    page.wait_for_function("n => state.splitChartHit.events.length === n", arg=len(before))
    # The marks do not move with the threshold, only the dashed cursor line does.
    assert page.evaluate("() => state.splitChartHit.events.map(e => e.x)") == before

    entry = _split_hit_events(page)[0]
    _hover_split_chart(page, entry["x"], entry["beads"][0]["y"])
    page.wait_for_selector("#split-chart-tip:not([hidden])")
    assert _split_tip_lines(page)[0].startswith("Threshold ")
    assert page.pageerrors == []


# ---------------------------------------------------------------------------
# Split chart axis ticks and chart export
# ---------------------------------------------------------------------------


def _x_tick_model(page):
    return page.evaluate(
        """() => splitChartLayout(splitCanvas.width, splitCanvas.height).xTicks
            .map(tick => ({value: tick.value, label: tick.label, x: tick.x}))"""
    )


def test_split_chart_tick_labels_name_the_value_they_sit_at(dense_page):
    """The reported bug: ticks were placed at fixed fractions of the domain, so one sat
    at 0.8736 and printed "0.87" through formatValue's two decimals -- which then looked
    misaligned beside a lollipop whose own readout also said 0.87. Every tick label must
    parse back to the exact value the tick is drawn at."""
    page = dense_page
    ticks = _x_tick_model(page)
    assert len(ticks) >= 3
    for tick in ticks:
        assert float(tick["label"].replace(",", "")) == pytest.approx(tick["value"], abs=1e-12), tick

    # Round numbers, not arbitrary fractions of the domain: the step between
    # consecutive ticks is constant and the values are multiples of it.
    steps = [
        round(right["value"] - left["value"], 12)
        for left, right in zip(ticks, ticks[1:])
    ]
    assert len(set(steps)) == 1, steps
    step = steps[0]
    for tick in ticks:
        assert abs(tick["value"] / step - round(tick["value"] / step)) < 1e-6, tick

    # And the axis mapping is honored, so label, value and pixel all agree.
    domain = page.evaluate(
        """() => {
            const layout = splitChartLayout(splitCanvas.width, splitCanvas.height);
            return {lo: layout.minThreshold, span: layout.thresholdSpan,
                    left: layout.margin.left, width: layout.width};
        }"""
    )
    for tick in ticks:
        expected = domain["left"] + ((tick["value"] - domain["lo"]) / domain["span"]) * domain["width"]
        assert tick["x"] == pytest.approx(expected, abs=1e-9)
    assert page.pageerrors == []


def test_split_chart_tick_labels_stay_distinct_on_a_narrow_axis(page):
    """An extraction can span a few thousandths, where two decimals print every tick
    identically. The label precision follows the tick step instead."""
    labels = page.evaluate(
        "() => splitAxisTicks(0.86612, 0.87104, 8).map(tick => tick.label)"
    )
    assert len(labels) == len(set(labels)), labels
    assert all(len(label.split(".")[1]) == 3 for label in labels), labels

    # A wide axis does not pay for that precision with trailing zeros.
    wide = page.evaluate("() => splitAxisTicks(0.4992, 0.9983, 8).map(tick => tick.label)")
    assert wide == ["0.5", "0.6", "0.7", "0.8", "0.9"]
    assert page.pageerrors == []


def test_split_chart_count_axes_use_whole_numbers(dense_page):
    """Both vertical axes count nodes, so a tick at 2.5 names nothing. Large counts keep
    the thousands separators the rest of the UI uses."""
    page = dense_page
    ticks = page.evaluate(
        """() => {
            const layout = splitChartLayout(splitCanvas.width, splitCanvas.height);
            return {y: layout.yTicks.map(t => t.label), y2: layout.y2Ticks.map(t => t.label)};
        }"""
    )
    for label in ticks["y"] + ticks["y2"]:
        assert "." not in label, label
        assert float(label.replace(",", "")).is_integer()

    assert page.evaluate(
        "() => splitAxisTicks(0, 3, 6, {integer: true}).map(tick => tick.label)"
    ) == ["0", "1", "2", "3"]
    assert page.evaluate(
        "() => splitAxisTicks(0, 55328, 6, {integer: true}).map(tick => tick.label)"
    ) == ["0", "10,000", "20,000", "30,000", "40,000", "50,000"]
    assert page.pageerrors == []


def test_split_chart_svg_draws_what_the_canvas_drew(dense_page):
    """The export shares splitChartLayout with the canvas painter and the hit-test, so it
    has one mark per mark and one tick label per tick."""
    page = dense_page
    svg = page.evaluate("() => buildSplitChartSVG()")
    model = page.evaluate(
        """() => {
            const layout = splitChartLayout(splitCanvas.width, splitCanvas.height);
            return {
                marks: layout.marks.length,
                beads: layout.marks.reduce((sum, mark) => sum + mark.beads.length, 0),
                tickLabels: layout.xTicks.concat(layout.yTicks, layout.y2Ticks).map(t => t.label),
                titles: [layout.titles.x, layout.titles.y, layout.titles.y2],
                hasMarker: layout.markerX !== null,
                movingSumPoints: layout.movingSum.points.length,
            };
        }"""
    )
    assert svg.startswith("<svg xmlns=")
    # One stem path per event, plus the axis frame, the right axis, the tick groups'
    # paths and the moving-sum trace.
    assert svg.count('<path d="M') >= model["marks"]
    # Exactly one bead circle per bead: the only other circle the chart has is the
    # threshold marker's dot, which an export leaves out.
    assert svg.count("<circle") == model["beads"]
    for label in model["tickLabels"]:
        assert ">" + label + "<" in svg, label
    for title in model["titles"]:
        assert title in svg
    assert page.pageerrors == []


def test_split_chart_exports_leave_out_the_threshold_cursor(dense_page):
    """The dashed line and orange dot say where the slider is, not anything the chart
    measures, so a figure taken from it should not carry them.

    The PNG path is checked by rendering the export twice at two different thresholds:
    the marks do not move with the threshold, so if the cursor were drawn the two images
    would differ, and if it is not they are byte-identical.
    """
    page = dense_page
    page.wait_for_function("() => splitChartLayout(1100, 320).markerX !== null")
    svg = page.evaluate("() => buildSplitChartSVG()")
    assert "stroke-dasharray" not in svg
    assert "#e29b4b" not in svg   # SPLIT_CHART_COLORS.markerDot

    # ... while the on-screen chart still draws it.
    assert page.evaluate(
        """() => {
            const canvas = document.createElement('canvas');
            canvas.width = splitCanvas.width;
            canvas.height = splitCanvas.height;
            const context = canvas.getContext('2d');
            drawSplitChart(context, splitCanvas.width, splitCanvas.height);
            const onScreen = canvas.toDataURL();
            context.clearRect(0, 0, canvas.width, canvas.height);
            drawSplitChart(context, splitCanvas.width, splitCanvas.height, {onScreen: false});
            return onScreen !== canvas.toDataURL();
        }"""
    ), "the live chart must still draw the threshold cursor"

    def exported_at(stop_index):
        return page.evaluate(
            """index => {
                snapSliderToStop(state.sliderModel.stops[index]);
                const canvas = document.createElement('canvas');
                canvas.width = splitCanvas.width;
                canvas.height = splitCanvas.height;
                drawSplitChart(canvas.getContext('2d'), splitCanvas.width, splitCanvas.height,
                    {onScreen: false});
                return canvas.toDataURL();
            }""",
            stop_index,
        )

    stop_count = page.evaluate("() => state.sliderModel.stops.length")
    first = exported_at(0)
    later = exported_at(stop_count // 2)
    assert first == later
    assert page.pageerrors == []


def test_split_chart_svg_export_is_well_formed_xml(dense_page, tmp_path):
    """Downloaded, then parsed: a malformed attribute or an unescaped label would make the
    file unopenable, which no assertion on the string would necessarily catch."""
    import xml.etree.ElementTree as ElementTree

    page = dense_page
    with page.expect_download() as download_info:
        page.click("#export-split-svg")
    download = download_info.value
    assert download.suggested_filename.endswith("_split_events.svg")
    target = tmp_path / download.suggested_filename
    download.save_as(target)

    root = ElementTree.parse(target).getroot()
    assert root.tag == "{http://www.w3.org/2000/svg}svg"
    assert root.get("width") == "1100"
    assert root.get("height") == "320"
    assert page.pageerrors == []


def test_split_chart_png_export_follows_the_resolution_selector(dense_page, tmp_path):
    """The PNG is re-rasterized at the chosen density rather than upscaled, so its pixel
    dimensions are the canvas's times the scale. Read out of the file's IHDR."""
    page = dense_page

    def exported_png(scale):
        page.select_option("#export-png-scale", scale)
        with page.expect_download() as download_info:
            page.click("#export-split-png")
        download = download_info.value
        target = tmp_path / (scale + "_" + download.suggested_filename)
        download.save_as(target)
        header = target.read_bytes()[:24]
        assert header[:8] == b"\x89PNG\r\n\x1a\n"
        return download.suggested_filename, (
            int.from_bytes(header[16:20], "big"),
            int.from_bytes(header[20:24], "big"),
        )

    name1, size1 = exported_png("1")
    name4, size4 = exported_png("4")
    assert name1.endswith("_split_events.png")
    assert name4.endswith("_split_events@4x.png")
    assert size1 == (1100, 320)
    assert size4 == (4400, 1280)
    assert page.pageerrors == []


def test_split_chart_export_leaves_the_hover_geometry_alone(dense_page):
    """The export paints into an offscreen canvas, and only the on-screen pass may record
    the hit-test geometry -- otherwise an export at 4x would leave the hover reading a
    chart nobody is looking at."""
    page = dense_page
    page.wait_for_function("() => state.splitChartHit !== null")
    before = page.evaluate("() => JSON.stringify(state.splitChartHit)")
    with page.expect_download() as download_info:
        page.click("#export-split-png")
    download_info.value  # settle the download before reading state back
    assert page.evaluate("() => JSON.stringify(state.splitChartHit)") == before

    entry = _split_hit_events(page)[0]
    _hover_split_chart(page, entry["x"], entry["beads"][0]["y"])
    page.wait_for_selector("#split-chart-tip:not([hidden])")
    assert _split_tip_lines(page)[0].startswith("Threshold ")
    assert page.pageerrors == []


# ---------------------------------------------------------------------------
# Label level-of-detail: gated on the mark's size, not on the zoom
# ---------------------------------------------------------------------------


def test_cluster_label_shows_whenever_the_bubble_has_room_for_it(page):
    """A big cluster keeps its label at any zoom, which a scale gate got wrong.

    A 160k-node bubble is thousands of world units across, so the whole layout only
    fits at a small zoom -- and the old `scale >= 0.11` rule then hid a label with
    hundreds of pixels of room. What matters is the bubble's size on screen.
    """
    shown = page.evaluate("""() => {
        const out = [];
        for (const scale of [0.02, 0.05, 0.11, 0.5]) {
            state.viewTransform.scale = scale;
            out.push({
                scale,
                bigWidth: Math.round(itemScreenExtent({radius: 3000}).width),
                big: clusterCountLabelFits('161,764', {radius: 3000}),
                small: clusterCountLabelFits('12', {radius: 40}),
            });
        }
        return out;
    }""")
    for row in shown:
        # The 3000-unit bubble is 120px wide even at scale 0.02, so it is labeled
        # at every zoom; the 40-unit one only once it is wide enough for "12".
        assert row["big"] is True, row
        assert row["bigWidth"] >= 100, row
    assert [row["small"] for row in shown] == [False, False, False, True], shown
    assert page.pageerrors == []


def test_cluster_label_is_dropped_when_the_text_does_not_fit(page):
    """Per-mark, so a crowded layout still drops the labels with no room."""
    verdicts = page.evaluate("""() => {
        state.viewTransform.scale = 1;
        return {
            roomy: clusterCountLabelFits('8', {radius: 30}),
            tooNarrow: clusterCountLabelFits('161,764', {radius: 12}),
            tooShort: clusterCountLabelFits('8', {radius: 3}),
        };
    }""")
    assert verdicts == {"roomy": True, "tooNarrow": False, "tooShort": False}
    assert page.pageerrors == []


def test_a_box_shaped_cluster_is_measured_by_its_width_not_its_diagonal(page):
    """`radius` is a half-diagonal for lattice/rect items and overstates the room.

    A tall, narrow treemap cluster has a long diagonal and almost no width, so the
    diagonal would let a label spill straight out of the box.
    """
    result = page.evaluate("""() => {
        state.viewTransform.scale = 1;
        const tall = {shape: 'lattice', x0: 0, x1: 18, y0: 0, y1: 400,
                      radius: Math.hypot(9, 200)};
        const wide = {shape: 'lattice', x0: 0, x1: 400, y0: 0, y1: 18,
                      radius: Math.hypot(200, 9)};
        return {
            tallExtent: itemScreenExtent(tall),
            tallLabel: clusterCountLabelFits('1,234', tall),
            wideLabel: clusterCountLabelFits('1,234', wide),
            diagonalWouldSay: Math.round(tall.radius * 2),
        };
    }""")
    assert result["tallExtent"] == {"width": 18, "height": 400}
    assert result["diagonalWouldSay"] > 300      # the trap being avoided
    assert result["tallLabel"] is False
    assert result["wideLabel"] is True
    assert page.pageerrors == []


def test_edge_score_label_is_gated_on_the_link_length(page):
    """The badge is drawn across the link's midpoint, so the link must be longer."""
    verdicts = page.evaluate("""() => ({
        long: edgeScoreLabelFits('9.87', 200),
        exact: edgeScoreLabelFits('9.87', 40),
        short: edgeScoreLabelFits('9.87', 20),
        longTextNeedsMoreRoom: edgeScoreLabelFits('123456.78', 40),
    })""")
    assert verdicts["long"] is True
    assert verdicts["short"] is False
    # A wider number needs a wider badge, so the same link no longer qualifies.
    assert verdicts["longTextNeedsMoreRoom"] is False
    assert page.pageerrors == []


def test_labels_are_drawn_at_low_zoom_when_the_marks_are_large(page):
    """End to end on the canvas: zoom out hard, keep the bubbles big, count the SVG.

    buildClusterViewSVG mirrors renderClusterView, so its <text> elements are a
    readable proxy for what the canvas just painted.
    """
    page.check("#show-node-counts")
    # Center the view as well as zooming it, so the off-screen cull is held constant
    # and the fit rule is the only thing under test.
    before = page.evaluate("""() => {
        state.viewTransform.scale = 0.04;
        state.viewTransform.offsetX = document.getElementById('cluster-view').width / 2;
        state.viewTransform.offsetY = document.getElementById('cluster-view').height / 2;
        renderClusterView();
        return (buildClusterViewSVG().match(/<text/g) || []).length;
    }""")
    # At that zoom this fixture's bubbles are a few pixels wide, so none is labeled...
    assert before == 0

    after = page.evaluate("""() => {
        // Same zoom, bubbles 100x bigger in world units: now they have room.
        state.visibleLayout.forEach(item => { item.radius *= 100; });
        renderClusterView();
        return (buildClusterViewSVG().match(/<text/g) || []).length;
    }""")
    assert after == page.evaluate("() => state.visibleLayout.length")
    assert page.pageerrors == []


def test_dot_labels_are_dropped_for_sub_pixel_dots(page):
    """A label beside an invisible dot points at nothing."""
    page.select_option("#label-by", "__node_id__")
    labeled = page.evaluate("""() => {
        const counts = {};
        for (const scale of [0.01, 1]) {
            state.viewTransform.scale = scale;
            state.viewTransform.offsetX = 400;
            state.viewTransform.offsetY = 300;
            renderClusterView();
            counts[scale] = (buildClusterViewSVG().match(/<text/g) || []).length;
        }
        return counts;
    }""")
    assert labeled["0.01"] == 0
    assert labeled["1"] > 0
    assert page.pageerrors == []


# ---------------------------------------------------------------------------
# "Collapse long paths": contracting chains of pass-through clusters
# ---------------------------------------------------------------------------


def _build_embedded_viewer_spindles(out_dir):
    """A network built to exercise path collapsing, with three shapes on purpose.

    Five 6-node blobs (internal weight 10, so they hold together at any threshold
    below that) plus:

    * ``b1 -- s0 -- s1 -- s2 -- b2``: a chain of three singletons whose four links
      are 3.0 / 3.5 / 2.5 / 3.2, so the chain's weakest link is **2.5**. This is
      the spindle the option exists to contract.
    * ``b1 -- t0 -- t1``: a dangling tail, which leaf pruning removes on its own.
    * ``h``: a singleton joined to ``b3``, ``b4`` and ``b5``. It is below any
      interesting minimum size but has three links, so it is a branch point rather
      than a pass-through and must survive.

    Every link outside the blobs is <= 3.5, so at the 3.5 stop the blobs are whole
    and every small cluster is a singleton.
    """
    names = []

    def blob(prefix, size=6):
        ids = [f"{prefix}{i}" for i in range(size)]
        names.extend(ids)
        return ids

    b1 = blob("b1_")
    names.extend(["s0", "s1", "s2"])
    b2 = blob("b2_")
    b3 = blob("b3_")
    names.append("h")
    b4 = blob("b4_")
    b5 = blob("b5_")
    names.extend(["t0", "t1"])

    index = {name: position for position, name in enumerate(names)}
    data = np.zeros((len(names), len(names)), dtype=float)

    def link(left, right, weight):
        data[index[left], index[right]] = weight
        data[index[right], index[left]] = weight

    for group in (b1, b2, b3, b4, b5):
        for position, left in enumerate(group):
            for right in group[position + 1:]:
                link(left, right, 10.0)

    link("b1_0", "s0", 3.0)
    link("s0", "s1", 3.5)
    link("s1", "s2", 2.5)
    link("s2", "b2_0", 3.2)
    link("b1_1", "t0", 3.05)
    link("t0", "t1", 2.9)
    link("h", "b3_0", 3.1)
    link("h", "b4_0", 3.3)
    link("h", "b5_0", 3.4)

    input_file = out_dir / "spindles.hdf5"
    html_file = out_dir / "viewer_spindles.html"
    DenseDataMatrix(data, names, names).write(str(input_file), output_type="dense")
    build_ssn_viewer.main([
        "-i", str(input_file),
        "--html", str(html_file),
        "--embed_data",
        "--name", "Spindle Test Viewer",
    ])
    return html_file


@pytest.fixture(scope="module")
def spindle_viewer_html(tmp_path_factory):
    return _build_embedded_viewer_spindles(
        tmp_path_factory.mktemp("ssn_viewer_spindles")
    )


@pytest.fixture
def spindle_page(spindle_viewer_html):
    yield from _yield_loaded_page(spindle_viewer_html)


def _spindle_state(page, layout="tree", min_cluster_size=2, expected_visible=9):
    """Put the spindle viewer at the 3.5 stop with leaf pruning on.

    ``tree`` is chosen because it draws edges and is computed synchronously, so
    ``state.splitLinks`` is populated without waiting on the layout worker.

    ``expected_visible`` is how many clusters are left standing once every setting
    here has landed: the five blobs plus the four pass-through singletons that leaf
    pruning keeps (s0, s1, s2, h), with the dangling tail (t0, t1) pruned away.
    """
    page.select_option("#layout-algorithm", layout)
    page.evaluate(
        """() => {
            snapSliderToStop(nearestStopForThreshold(3.5));
            scheduleThresholdUI(true);
        }"""
    )
    page.fill("#min-cluster-size", str(min_cluster_size))
    page.check("#leaf-pruning-only")
    # Wait on the outcome, not on layoutComputing: the re-render that the fill and
    # the check schedule is one animation frame away, and layoutComputing is false
    # both before it starts and after it finishes, so waiting on it alone can return
    # while the page still shows the pre-pruning state. No state on the way here has
    # this count -- at the 3.5 stop it is 11 clusters before the fill (minimum size 1,
    # everything visible) and 5 after it but before the check (pruning off, so only
    # the blobs clear the minimum) -- so this settles only once both have been applied.
    page.wait_for_function(
        "expected => selectedThresholdValue() === 3.5 && !state.layoutComputing"
        " && state.visibleClusters.length === expected",
        arg=expected_visible,
    )


def _visible_sizes(page):
    return page.evaluate(
        """() => state.visibleClusters
            .map(id => state.bundle.graph.hierarchy.nodes[id].size)
            .sort((left, right) => left - right)"""
    )


def _link_summary(page):
    return page.evaluate(
        """() => state.splitLinks
            .map(link => ({weight: link.threshold, collapsed: link.collapsed}))
            .sort((left, right) => left.weight - right.weight)"""
    )


def test_collapse_long_paths_is_a_sub_option_of_leaf_pruning(spindle_page):
    """The checkbox is only available while leaf pruning is on: with leaf pruning off
    every below-minimum cluster is dropped outright, so no chain is left to contract."""
    page = spindle_page
    assert page.is_disabled("#collapse-long-paths")
    page.check("#leaf-pruning-only")
    assert page.is_enabled("#collapse-long-paths")
    page.uncheck("#leaf-pruning-only")
    assert page.is_disabled("#collapse-long-paths")
    assert page.pageerrors == []


def test_collapse_long_paths_contracts_a_chain_into_one_edge(spindle_page):
    """The headline behavior: three pass-through singletons and their four links become
    one link carrying the chain's weakest weight."""
    page = spindle_page
    _spindle_state(page)

    # Leaf pruning alone keeps the whole spindle: five blobs, the three chain
    # singletons and the branch point, joined by seven links.
    assert _visible_sizes(page) == [1, 1, 1, 1, 6, 6, 6, 6, 6]
    assert [link["weight"] for link in _link_summary(page)] == [2.5, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5]
    assert all(link["collapsed"] == 0 for link in _link_summary(page))

    page.check("#collapse-long-paths")
    page.wait_for_function("() => state.visibleClusters.length === 6")

    # The chain is gone; the branch point (size 1, three links) is not.
    assert _visible_sizes(page) == [1, 6, 6, 6, 6, 6]
    links = _link_summary(page)
    assert [link["weight"] for link in links] == [2.5, 3.1, 3.3, 3.4]
    # The one collapsed link stands for all three contracted clusters and carries the
    # chain's minimum weight (2.5), not the weight of either end of the chain.
    collapsed = [link for link in links if link["collapsed"]]
    assert collapsed == [{"weight": 2.5, "collapsed": 3}]
    assert page.pageerrors == []


def test_collapse_long_paths_counts_contracted_clusters_as_hidden(spindle_page):
    """A contracted cluster is not drawn, so it is hidden like any other below-minimum
    cluster -- and the summary says how many paths went away."""
    page = spindle_page
    _spindle_state(page)

    # The dangling tail (t0, t1) is pruned by leaf pruning whether or not paths collapse.
    assert page.locator("#stat-hidden-nodes").inner_text() == "2"
    assert (
        page.locator("#hidden-summary").inner_text()
        == "2 nodes hidden by minimum cluster size"
    )

    page.check("#collapse-long-paths")
    page.wait_for_function("() => state.visibleClusters.length === 6")
    assert page.locator("#stat-hidden-nodes").inner_text() == "5"
    assert (
        page.locator("#hidden-summary").inner_text()
        == "5 nodes hidden by minimum cluster size, 1 path collapsed"
    )
    assert page.pageerrors == []


def test_collapse_long_paths_does_nothing_while_leaf_pruning_is_off(spindle_page):
    """Checked but inert: with leaf pruning off the plain minimum-size rule applies, so
    the result must be identical either way."""
    page = spindle_page
    _spindle_state(page)
    page.uncheck("#leaf-pruning-only")
    # Waiting on the outcome, not on layoutComputing: the re-render the uncheck
    # schedules is one animation frame away, and layoutComputing is false both
    # before it starts and after it finishes. Once the plain minimum-size rule is
    # in force no cluster below it is left standing, which is the state to compare.
    page.wait_for_function(
        "() => !state.layoutComputing && state.visibleClusters.every("
        "  id => state.bundle.graph.hierarchy.nodes[id].size >= 2)"
    )
    before = {"sizes": _visible_sizes(page), "links": _link_summary(page)}

    # The control is disabled, so drive it the way a restored session would.
    page.evaluate(
        """() => {
            const input = document.getElementById('collapse-long-paths');
            input.checked = true;
            input.dispatchEvent(new Event('change', {bubbles: true}));
        }"""
    )
    page.wait_for_function("() => !state.layoutComputing")

    assert {"sizes": _visible_sizes(page), "links": _link_summary(page)} == before
    assert page.locator("#hidden-summary").inner_text() == "6 nodes hidden by minimum cluster size"
    assert page.pageerrors == []


def test_collapsed_edge_weight_is_the_path_minimum(spindle_page):
    """Checked against the uncollapsed graph rather than against a hand-written number:
    for every collapsed link, a plain BFS over the pre-collapse links must find a path
    whose interior is all pass-through clusters and whose minimum weight is the link's."""
    page = spindle_page
    _spindle_state(page)

    problems = page.evaluate(
        """() => {
            const nodes = state.bundle.graph.hierarchy.nodes;
            const active = activeClustersAtThreshold(selectedThresholdValue());
            const minSize = Number(document.getElementById('min-cluster-size').value);
            const pruned = mstLinksForActiveClusters(active, minSize, true, false);
            const collapsed = mstLinksForActiveClusters(active, minSize, true, true);

            const visible = new Set(pruned.visibleIds);
            const adjacency = new Map();
            pruned.links.forEach(link => {
                if (!adjacency.has(link.sourceId)) { adjacency.set(link.sourceId, []); }
                if (!adjacency.has(link.targetId)) { adjacency.set(link.targetId, []); }
                adjacency.get(link.sourceId).push({other: link.targetId, weight: link.weight});
                adjacency.get(link.targetId).push({other: link.sourceId, weight: link.weight});
            });
            const passThrough = id => visible.has(id)
                && (adjacency.get(id) || []).length === 2
                && nodes[id].size < minSize;

            const problems = [];
            const collapsedLinks = collapsed.links.filter(link => link.collapsed);
            if (collapsedLinks.length === 0) { problems.push('nothing collapsed'); }
            collapsedLinks.forEach(link => {
                const previous = new Map([[link.sourceId, null]]);
                const queue = [link.sourceId];
                for (let i = 0; i < queue.length; i++) {
                    (adjacency.get(queue[i]) || []).forEach(edge => {
                        if (previous.has(edge.other)) { return; }
                        previous.set(edge.other, {from: queue[i], weight: edge.weight});
                        queue.push(edge.other);
                    });
                }
                if (!previous.has(link.targetId)) {
                    problems.push('no pre-collapse path for ' + link.sourceId + '-' + link.targetId);
                    return;
                }
                const interior = [];
                let current = link.targetId;
                let minWeight = Infinity;
                while (previous.get(current)) {
                    const step = previous.get(current);
                    minWeight = Math.min(minWeight, step.weight);
                    if (current !== link.targetId) { interior.push(current); }
                    current = step.from;
                }
                if (Math.abs(minWeight - link.weight) > 1e-12) {
                    problems.push('weight ' + link.weight + ' is not the path minimum ' + minWeight);
                }
                if (interior.length !== link.collapsed) {
                    problems.push('collapsed ' + link.collapsed + ' but interior is ' + interior.length);
                }
                if (!interior.every(passThrough)) { problems.push('interior is not all pass-through'); }
                if (passThrough(link.sourceId) || passThrough(link.targetId)) {
                    problems.push('an endpoint is itself pass-through');
                }
            });
            return problems;
        }"""
    )
    assert problems == []
    assert page.pageerrors == []


def test_collapse_long_paths_preserves_which_clusters_are_connected(spindle_page):
    """Contracting a path must not change the relationships it stands for: two surviving
    clusters are connected after the collapse exactly when they were before."""
    page = spindle_page
    _spindle_state(page)

    changed = page.evaluate(
        """() => {
            const active = activeClustersAtThreshold(selectedThresholdValue());
            const minSize = Number(document.getElementById('min-cluster-size').value);
            const pruned = mstLinksForActiveClusters(active, minSize, true, false);
            const collapsed = mstLinksForActiveClusters(active, minSize, true, true);

            const componentOf = (ids, links) => {
                const adjacency = new Map(ids.map(id => [id, []]));
                links.forEach(link => {
                    adjacency.get(link.sourceId)?.push(link.targetId);
                    adjacency.get(link.targetId)?.push(link.sourceId);
                });
                const label = new Map();
                let next = 0;
                ids.forEach(id => {
                    if (label.has(id)) { return; }
                    const queue = [id];
                    label.set(id, next);
                    for (let i = 0; i < queue.length; i++) {
                        (adjacency.get(queue[i]) || []).forEach(other => {
                            if (label.has(other)) { return; }
                            label.set(other, next);
                            queue.push(other);
                        });
                    }
                    next += 1;
                });
                return label;
            };

            const before = componentOf(pruned.visibleIds, pruned.links);
            const after = componentOf(collapsed.visibleIds, collapsed.links);
            const survivors = collapsed.visibleIds;
            const changed = [];
            for (let i = 0; i < survivors.length; i++) {
                for (let j = i + 1; j < survivors.length; j++) {
                    const sameBefore = before.get(survivors[i]) === before.get(survivors[j]);
                    const sameAfter = after.get(survivors[i]) === after.get(survivors[j]);
                    if (sameBefore !== sameAfter) { changed.push([survivors[i], survivors[j]]); }
                }
            }
            return changed;
        }"""
    )
    assert changed == []
    assert page.pageerrors == []


def test_collapsed_edges_export_dashed(spindle_page):
    """A collapsed path is the weakest link of a contracted chain, not a measured edge
    between the two clusters it joins, so it is drawn dashed -- in the export too."""
    page = spindle_page
    _spindle_state(page)
    solid_only = page.evaluate("() => buildClusterViewSVG()")
    assert 'stroke-dasharray' not in solid_only

    page.check("#collapse-long-paths")
    page.wait_for_function("() => state.visibleClusters.length === 6")
    svg = page.evaluate("() => buildClusterViewSVG()")

    # One dashed group holding exactly the collapsed links, and the surviving real
    # edges still in an undashed group of their own.
    assert svg.count('stroke-dasharray="6 4"') == 1
    dashed_group = svg.split('stroke-dasharray="6 4">')[1].split("</g>")[0]
    assert dashed_group.count("<path") == 1
    assert page.evaluate("() => state.splitLinks.filter(link => link.collapsed).length") == 1
    assert page.pageerrors == []


def test_collapse_long_paths_round_trips_through_a_session(spindle_page):
    """The toggle is part of the saved view state, and restoring it also re-derives
    whether the control is available."""
    page = spindle_page
    _spindle_state(page)
    page.check("#collapse-long-paths")
    page.wait_for_function("() => state.visibleClusters.length === 6")

    saved = page.evaluate("() => collectSessionState()")
    assert saved["view"]["collapse_long_paths"] is True
    assert saved["view"]["leaf_pruning_only"] is True

    # A session saved with leaf pruning off keeps the sub-option's value but must
    # restore it as unavailable.
    note = page.evaluate(
        """() => applySessionState({
            view: {leaf_pruning_only: false, collapse_long_paths: true},
        })"""
    )
    assert note == ""
    assert page.evaluate("() => document.getElementById('collapse-long-paths').checked") is True
    assert page.is_disabled("#collapse-long-paths")
    assert page.pageerrors == []


def test_collapse_long_paths_is_stable_across_layouts_and_thresholds(spindle_page):
    """Every threshold stop x minimum size x layout, with the toggle on: no JS errors and
    no link pointing at a cluster that is not laid out."""
    page = spindle_page
    page.check("#leaf-pruning-only")
    page.check("#collapse-long-paths")
    for layout in ("tree", "packed", "treemap"):
        page.select_option("#layout-algorithm", layout)
        for stop_index in range(len(page.evaluate("() => state.sliderModel.stops"))):
            page.evaluate(
                "index => { snapSliderToStop(state.sliderModel.stops[index]);"
                " scheduleThresholdUI(false); }",
                stop_index,
            )
            for size in (1, 2, 4, 8, 40):
                page.fill("#min-cluster-size", str(size))
                page.wait_for_function("() => !state.layoutComputing")
                assert page.evaluate(
                    """() => {
                        const laidOut = new Set(state.visibleLayout.map(item => item.componentId));
                        return state.splitLinks.every(link => Number.isFinite(link.threshold)
                            && laidOut.has(link.left.componentId)
                            && laidOut.has(link.right.componentId));
                    }"""
                )
    assert page.pageerrors == []


# ---------------------------------------------------------------------------
# Crowded thresholds: stepping between stops, and zooming the split chart
# ---------------------------------------------------------------------------


def _build_embedded_viewer_crowded(out_dir, blob_count=10, blob_size=20):
    """A network whose merges crowd into a narrow band of scores.

    Ten 20-node blobs, internally scored 0.90-0.99, chained to each other by
    single links scored 0.80-0.89. Every score is distinct, so every merge gets
    its own slider stop -- but they are packed into a fifth of the axis, so many
    stops round to the same integer slider position. That collision is what made
    the threshold arrows appear to stick, and it is also what makes this network
    worth zooming into: at full extent its lollipops pile into one column.
    """
    rng = np.random.default_rng(11)
    node_count = blob_count * blob_size
    data = np.zeros((node_count, node_count), dtype=float)
    for blob in range(blob_count):
        low, high = blob * blob_size, (blob + 1) * blob_size
        block = rng.uniform(0.90, 0.99, size=(blob_size, blob_size))
        block = np.triu(block, 1)
        data[low:high, low:high] = block + block.T
    for blob in range(blob_count - 1):
        left, right = (blob * blob_size) + 3, ((blob + 1) * blob_size) + 7
        data[left, right] = data[right, left] = rng.uniform(0.80, 0.89)

    names = [f"n{i:03d}" for i in range(node_count)]
    input_file = out_dir / "crowded.hdf5"
    html_file = out_dir / "viewer_crowded.html"
    DenseDataMatrix(data, names, names).write(str(input_file), output_type="dense")
    build_ssn_viewer.main([
        "-i", str(input_file),
        "--html", str(html_file),
        "--embed_data",
        "--name", "Crowded Test Viewer",
    ])
    return html_file


@pytest.fixture(scope="module")
def crowded_viewer_html(tmp_path_factory):
    return _build_embedded_viewer_crowded(
        tmp_path_factory.mktemp("ssn_viewer_crowded")
    )


@pytest.fixture
def crowded_page(crowded_viewer_html):
    """Every merge gets a stop, which is the crowding this fixture exists to show.

    The viewer's default cap selects from the chart's *window*, so at full extent it
    deliberately thins the track; uncapping is how this fixture gets the dense track its
    tests are about.
    """
    for page in _yield_loaded_page(crowded_viewer_html):
        _set_split_event_cap(page, 0)
        yield page


def _split_window(page):
    """The x window the last paint used, as the hover geometry recorded it."""
    return page.evaluate("""() => ({
        min: state.splitChartHit.minThreshold,
        span: state.splitChartHit.thresholdSpan,
        marks: state.splitChartHit.events.length,
    })""")


def _split_chart_fraction_point(page, fraction=0.5):
    """Client coordinates of a point `fraction` of the way across the chart."""
    box = page.locator("#split-chart").bounding_box()
    return box["x"] + (box["width"] * fraction), box["y"] + (box["height"] * 0.4)


def _threshold_under(page, client_x):
    return page.evaluate(
        """clientX => {
            const rect = splitCanvas.getBoundingClientRect();
            return splitChartThresholdAt((clientX - rect.left) * (splitCanvas.width / rect.width));
        }""",
        client_x,
    )


def test_crowded_stops_each_get_their_own_slider_position(crowded_page):
    """The premise of the stepping tests, and the thing the warped map answers.

    These stops crowd: scaled straight onto the track by value, most of them land
    on a position some other stop already has, and only the first stop at each
    position can ever be selected by dragging. The map reserves a position per stop
    and spends what is left over on saying where they are, so the crowding stays
    visible in the spacing without costing anyone their own slot.
    """
    page = crowded_page
    measured = page.evaluate("""() => {
        const stops = state.sliderModel.stops;
        const finite = stops.filter(stop => stop.threshold_value !== null);
        const low = finite[0].threshold_value;
        const high = finite[finite.length - 1].threshold_value;
        // What a straight value-linear scale would have done.
        const naive = finite.map(stop => Math.round(
            ((stop.threshold_value - low) / (high - low)) * 920));
        // What a drag can actually select: every position the slider can take,
        // resolved the way a drag resolves it -- by position, with no remembered stop.
        const slider = document.getElementById('threshold-slider');
        const keep = slider.value;
        const reachable = new Set();
        for (let position = 0; position <= state.sliderModel.maxPosition; position++) {
            slider.value = String(position);
            state.selectedStop = null;
            reachable.add(stops.indexOf(currentSliderStop()));
        }
        slider.value = keep;
        return {
            stops: stops.length,
            naivePositions: new Set(naive).size,
            actualPositions: new Set(stops.map(stop => stop.sliderPosition)).size,
            monotone: finite.every((stop, index) =>
                index === 0 || stop.sliderPosition > finite[index - 1].sliderPosition),
            ends: [finite[0].sliderPosition, finite[finite.length - 1].sliderPosition],
            infinityAt: stops[stops.length - 1].sliderPosition,
            reachableByDragging: reachable.size,
        };
    }""")
    # The crowding is real: a straight scale would lose most of these stops.
    assert measured["naivePositions"] < measured["stops"] / 2
    # And none of them are lost.
    assert measured["actualPositions"] == measured["stops"]
    assert measured["reachableByDragging"] == measured["stops"]
    # Still ordered by threshold, and still spanning the whole track.
    assert measured["monotone"]
    assert measured["ends"] == [0, 920]
    assert measured["infinityAt"] == 1000
    assert page.pageerrors == []


def test_slider_positions_track_threshold_where_there_is_room(dense_page):
    """Uncrowded stops keep the plain proportional placement they always had.

    The map only borrows track where stops would otherwise collide, so on a network
    whose merges are spread out it should be indistinguishable from a straight scale
    -- which is what keeps the track reading like the chart above it.
    """
    page = dense_page
    drift = page.evaluate("""() => {
        const finite = state.sliderModel.stops.filter(s => s.threshold_value !== null);
        const low = finite[0].threshold_value;
        const high = finite[finite.length - 1].threshold_value;
        return finite.map(stop => Math.abs(stop.sliderPosition - Math.round(
            ((stop.threshold_value - low) / (high - low)) * 920)));
    }""")
    # Within 5% of the track of where a straight value-linear scale would put them.
    assert max(drift) < 46
    assert page.pageerrors == []


def test_slider_positions_degrade_gracefully_past_the_track(page):
    """The two branches no real bundle here reaches, driven directly.

    More stops than the track has positions (only with the split-event cap raised
    past ~900) and every stop at one threshold both have to stay ordered and inside
    the track; what they cannot keep is a position each, which is the one thing a
    thousand positions cannot give a thousand stops.
    """
    overflowing = page.evaluate("""() => {
        const stops = Array.from({length: 1000}, (_, index) => ({threshold_value: 0.5 + (index * 1e-4)}));
        positionSliderStops(stops);
        const positions = stops.map(stop => stop.sliderPosition);
        return {
            nonDecreasing: positions.every((p, i) => i === 0 || p >= positions[i - 1]),
            withinTrack: Math.min(...positions) >= 0 && Math.max(...positions) <= 920,
            ends: [positions[0], positions[positions.length - 1]],
        };
    }""")
    assert overflowing["nonDecreasing"]
    assert overflowing["withinTrack"]
    assert overflowing["ends"] == [0, 920]

    flat = page.evaluate("""() => {
        const stops = Array.from({length: 5}, () => ({threshold_value: 7.0}));
        positionSliderStops(stops);
        return stops.map(stop => stop.sliderPosition);
    }""")
    # Nothing to be proportional to, so they spread by rank instead.
    assert flat == [0, 230, 460, 690, 920]
    assert page.pageerrors == []


def test_zooming_the_chart_magnifies_that_band_of_the_slider(crowded_page):
    """The track is warped by the chart's window, which is what pays for the room."""
    page = crowded_page

    def band_track():
        return page.evaluate("""() => {
            const finite = state.sliderModel.stops.filter(s => s.threshold_value !== null);
            const window = state.splitChartZoom;
            const band = window
                ? finite.filter(s => s.threshold_value >= window.min && s.threshold_value <= window.max)
                : finite;
            return {
                count: band.length,
                track: band[band.length - 1].sliderPosition - band[0].sliderPosition,
                distinct: new Set(finite.map(s => s.sliderPosition)).size,
                total: finite.length,
            };
        }""")

    threshold_before = page.evaluate("() => selectedThresholdValue()")
    x, y = _split_chart_fraction_point(page, fraction=0.85)
    page.mouse.move(x, y)
    page.mouse.wheel(0, -1500)
    page.wait_for_function("() => Boolean(state.splitChartZoom)")
    zoomed_window = page.evaluate("() => state.splitChartZoom")

    # Measure the same stops before and after, by pinning the comparison to the
    # window the zoom settled on rather than to a band chosen in advance.
    zoomed = band_track()
    page.click("#split-chart-reset-zoom")
    page.wait_for_function("() => state.splitChartZoom === null")
    flat = page.evaluate("""window => {
        const finite = state.sliderModel.stops.filter(s => s.threshold_value !== null);
        const band = finite.filter(s => s.threshold_value >= window.min && s.threshold_value <= window.max);
        return {count: band.length,
                track: band[band.length - 1].sliderPosition - band[0].sliderPosition};
    }""", zoomed_window)

    assert zoomed["count"] == flat["count"]
    assert zoomed["track"] > flat["track"]
    # Zooming buys the band room; it never costs another stop its own position.
    assert zoomed["distinct"] == zoomed["total"]
    # And the threshold does not move just because the map under it did.
    assert page.evaluate("() => selectedThresholdValue()") == threshold_before
    assert page.pageerrors == []


def test_threshold_arrows_land_on_a_new_threshold_every_press(crowded_page):
    """→ walks every stop to ∞, and ← walks back, without ever standing still.

    Resolving the current stop from the slider's position alone returns the first
    stop sharing that position, so pressing → from any of the others put the
    slider back where it already was. The walk below is the regression.
    """
    page = crowded_page
    walk = page.evaluate("""() => {
        const stops = state.sliderModel.stops;
        snapSliderToStop(stops[0]);
        const climbed = [];
        for (let press = 0; press < stops.length + 5; press++) {
            const before = selectedThresholdValue();
            stepThreshold(1);
            const after = selectedThresholdValue();
            if (after === before) { break; }
            climbed.push(after);
        }
        const descended = [];
        for (let press = 0; press < stops.length + 5; press++) {
            const before = selectedThresholdValue();
            stepThreshold(-1);
            const after = selectedThresholdValue();
            if (after === before) { break; }
            descended.push(after);
        }
        return {
            stops: stops.length,
            climbed: climbed.length,
            descended: descended.length,
            strictlyRising: climbed.every((value, i) => i === 0 || value > climbed[i - 1]),
            endedAtInfinity: !Number.isFinite(climbed[climbed.length - 1]),
            endedAtFloor: descended[descended.length - 1] === stops[0].threshold_value,
        };
    }""")
    assert walk["climbed"] == walk["stops"] - 1
    assert walk["descended"] == walk["stops"] - 1
    assert walk["strictlyRising"]
    assert walk["endedAtInfinity"]
    assert walk["endedAtFloor"]
    assert page.pageerrors == []


def test_threshold_arrows_grey_out_only_at_the_ends(crowded_page):
    """An arrow is enabled exactly when pressing it would move the threshold."""
    page = crowded_page
    ends = page.evaluate("""() => {
        const stops = state.sliderModel.stops;
        const read = () => ({
            down: document.getElementById('threshold-step-down').disabled,
            up: document.getElementById('threshold-step-up').disabled,
        });
        snapSliderToStop(stops[0]);
        updateThresholdStepButtons();
        const bottom = read();
        snapSliderToStop(stops[Math.floor(stops.length / 2)]);
        updateThresholdStepButtons();
        const middle = read();
        snapSliderToStop(stops[stops.length - 1]);
        updateThresholdStepButtons();
        return {bottom, middle, top: read()};
    }""")
    assert ends["bottom"] == {"down": True, "up": False}
    assert ends["middle"] == {"down": False, "up": False}
    assert ends["top"] == {"down": False, "up": True}
    assert page.pageerrors == []


def test_split_chart_wheel_zooms_about_the_pointer(crowded_page):
    """Scrolling narrows the window and leaves the threshold under the cursor put."""
    page = crowded_page
    assert page.is_disabled("#split-chart-reset-zoom")
    assert page.text_content("#split-chart-zoom-hint") == "Scroll to zoom · drag to pan"

    x, y = _split_chart_fraction_point(page, fraction=0.6)
    page.mouse.move(x, y)
    before = _split_window(page)
    anchor_before = _threshold_under(page, x)

    page.mouse.wheel(0, -600)
    page.wait_for_function("() => Boolean(state.splitChartZoom)")
    after = _split_window(page)
    anchor_after = _threshold_under(page, x)

    assert after["span"] < before["span"] / 1.5
    # The lens moves, not the chart: a tenth of a percent of the window's width is
    # the pointer's own rounding into canvas pixels.
    assert abs(anchor_after - anchor_before) < after["span"] * 0.01
    assert not page.is_disabled("#split-chart-reset-zoom")
    # The readout names the window, at enough precision for its two ends to differ --
    # a fixed two decimals would print a tight window as "0.87 - 0.87".
    hint = page.text_content("#split-chart-zoom-hint")
    assert "double-click to reset" in hint
    low, high = (float(part) for part in hint.split("Showing ")[1].split(" ·")[0].split(" – "))
    assert low < high
    assert low == pytest.approx(after["min"], abs=after["span"] / 50)
    assert high == pytest.approx(after["min"] + after["span"], abs=after["span"] / 50)
    assert page.pageerrors == []


def test_split_chart_zoom_stops_at_the_whole_series(crowded_page):
    """Scrolling out past the ends leaves no window at all, rather than a wider one."""
    page = crowded_page
    x, y = _split_chart_fraction_point(page)
    page.mouse.move(x, y)
    page.mouse.wheel(0, -600)
    page.wait_for_function("() => Boolean(state.splitChartZoom)")
    full_span = page.evaluate("""() => {
        const hit = state.splitChartHit;
        return hit.dataMax - hit.dataMin;
    }""")

    page.mouse.wheel(0, 4000)
    page.wait_for_function("() => state.splitChartZoom === null")
    assert _split_window(page)["span"] == pytest.approx(full_span)
    assert page.is_disabled("#split-chart-reset-zoom")
    assert page.pageerrors == []


def test_split_chart_drag_pans_without_setting_the_threshold(crowded_page):
    """Dragging slides the window; the release is a pan, not a click-to-jump."""
    page = crowded_page
    x, y = _split_chart_fraction_point(page)
    page.mouse.move(x, y)
    page.mouse.wheel(0, -600)
    page.wait_for_function("() => Boolean(state.splitChartZoom)")
    before = _split_window(page)
    threshold_before = page.evaluate("() => selectedThresholdValue()")

    page.mouse.move(x, y)
    page.mouse.down()
    page.mouse.move(x - 150, y, steps=8)
    page.mouse.up()
    after = _split_window(page)

    # Dragging left moves the window to higher thresholds, and only moves it.
    assert after["min"] > before["min"]
    assert after["span"] == pytest.approx(before["span"])
    assert page.evaluate("() => selectedThresholdValue()") == threshold_before
    assert page.pageerrors == []


def test_split_chart_click_still_jumps_after_a_drag(crowded_page):
    """The click suppression lasts for the drag's own release, not beyond it."""
    page = crowded_page
    x, y = _split_chart_fraction_point(page)
    page.mouse.move(x, y)
    page.mouse.down()
    page.mouse.move(x - 120, y, steps=6)
    page.mouse.up()
    threshold_after_drag = page.evaluate("() => selectedThresholdValue()")

    page.mouse.click(x, y)
    page.wait_for_function(
        "before => selectedThresholdValue() !== before", arg=threshold_after_drag
    )
    assert page.pageerrors == []


def test_split_chart_zoom_drops_the_marks_outside_the_window(crowded_page):
    """Hit-testing must not offer an event the chart is not drawing."""
    page = crowded_page
    before = _split_window(page)
    x, y = _split_chart_fraction_point(page)
    page.mouse.move(x, y)
    page.mouse.wheel(0, -900)
    page.wait_for_function("() => Boolean(state.splitChartZoom)")
    after = _split_window(page)

    assert after["marks"] < before["marks"]
    assert page.evaluate("""() => {
        const hit = state.splitChartHit;
        const slack = hit.thresholdSpan * 0.02;
        return hit.events.every(entry =>
            entry.event.threshold_value >= hit.minThreshold - slack &&
            entry.event.threshold_value <= hit.minThreshold + hit.thresholdSpan + slack);
    }""")
    assert page.pageerrors == []


def test_split_chart_double_click_and_button_reset_the_zoom(crowded_page):
    page = crowded_page
    x, y = _split_chart_fraction_point(page)
    for reset in ("dblclick", "button"):
        page.mouse.move(x, y)
        page.mouse.wheel(0, -600)
        page.wait_for_function("() => Boolean(state.splitChartZoom)")
        if reset == "dblclick":
            page.mouse.dblclick(x, y)
        else:
            page.click("#split-chart-reset-zoom")
        page.wait_for_function("() => state.splitChartZoom === null")
        assert page.is_disabled("#split-chart-reset-zoom")
        assert page.text_content("#split-chart-zoom-hint") == "Scroll to zoom · drag to pan"
    assert page.pageerrors == []


def test_split_chart_double_click_leaves_the_threshold_alone(crowded_page):
    """Its two clicks jump the threshold; the reset gesture has to undo them.

    A click cannot know a second one is coming, so the alternative would be making
    every single click wait out the double-click interval.
    """
    page = crowded_page
    x, y = _split_chart_fraction_point(page, fraction=0.3)
    page.mouse.move(x, y)
    page.mouse.wheel(0, -600)
    page.wait_for_function("() => Boolean(state.splitChartZoom)")
    before = page.evaluate("() => selectedThresholdValue()")

    page.mouse.dblclick(x, y)
    page.wait_for_function("() => state.splitChartZoom === null")
    assert page.evaluate("() => selectedThresholdValue()") == before

    # A single click on the same spot does move it, which is the whole point of
    # having to undo it above.
    page.mouse.click(x, y)
    page.wait_for_function("before => selectedThresholdValue() !== before", arg=before)
    assert page.pageerrors == []


def test_split_chart_threshold_marker_leaves_the_window(crowded_page):
    """A marker pinned to the edge would claim the threshold is there."""
    page = crowded_page
    markers = page.evaluate("""() => {
        const hit = state.splitChartHit;
        const span = hit.dataMax - hit.dataMin;
        const layout = () => splitChartLayout(splitCanvas.width, splitCanvas.height).markerX;
        const results = {};
        snapSliderToStop(state.sliderModel.stops[state.sliderModel.stops.length - 1]);
        results.infinityFullRange = layout();
        // The floor stop sits below every plotted event, so its marker belongs at the
        // left edge of the plot box rather than outside it.
        snapSliderToStop(state.sliderModel.stops[0]);
        results.floorFullRange = layout();
        results.plotLeft = hit.plotLeft;
        snapSliderToStop(state.sliderModel.stops[state.sliderModel.stops.length - 1]);
        state.splitChartZoom = {min: hit.dataMin, max: hit.dataMin + (span * 0.05)};
        results.infinityZoomedAway = layout();
        snapSliderToStop(nearestStopForThreshold(hit.dataMin + (span * 0.02)));
        results.insideWindow = layout();
        state.splitChartZoom = {min: hit.dataMax - (span * 0.05), max: hit.dataMax};
        results.outsideWindow = layout();
        state.splitChartZoom = null;
        drawSplitChart();
        return results;
    }""")
    assert markers["infinityFullRange"] is not None
    assert markers["floorFullRange"] == pytest.approx(markers["plotLeft"])
    assert markers["infinityZoomedAway"] is None
    assert markers["insideWindow"] is not None
    assert markers["outsideWindow"] is None
    assert page.pageerrors == []


def test_split_chart_export_follows_the_zoom(crowded_page):
    """The export is the chart in front of you, clipped the same way."""
    page = crowded_page
    x, y = _split_chart_fraction_point(page)
    page.mouse.move(x, y)
    page.mouse.wheel(0, -900)
    page.wait_for_function("() => Boolean(state.splitChartZoom)")

    model = page.evaluate("""() => {
        const layout = splitChartLayout(splitCanvas.width, splitCanvas.height);
        return {
            stems: layout.marks.length,
            beads: layout.marks.reduce((sum, mark) => sum + mark.beads.length, 0),
            firstTick: layout.xTicks[0].label,
        };
    }""")
    svg = page.evaluate("() => buildSplitChartSVG()")
    assert svg.count("<circle") == model["beads"]
    assert svg.count("<path d=\"M") >= model["stems"]
    assert model["firstTick"] in svg
    # One clip path, applied to each group of data marks: stems, beads, the line.
    assert svg.count('<clipPath id="split-plot-clip">') == 1
    assert svg.count('clip-path="url(#split-plot-clip)"') == 3
    assert page.pageerrors == []


def test_split_chart_zoom_round_trips_through_a_session(crowded_page):
    """Saved like the canvas's own pan and zoom, and clamped on the way back in."""
    page = crowded_page
    x, y = _split_chart_fraction_point(page)
    page.mouse.move(x, y)
    page.mouse.wheel(0, -600)
    page.wait_for_function("() => Boolean(state.splitChartZoom)")
    saved = page.evaluate("() => collectSessionState().view.split_chart_zoom")
    assert saved["max"] > saved["min"]

    page.click("#split-chart-reset-zoom")
    assert page.evaluate("() => collectSessionState().view.split_chart_zoom") is None

    note = page.evaluate("value => applySessionState({view: {split_chart_zoom: value}})", saved)
    assert note == ""
    assert page.evaluate("() => state.splitChartZoom") == pytest.approx(saved)

    # A window wider than this bundle's own range is re-fitted, not obeyed.
    page.evaluate("() => applySessionState({view: {split_chart_zoom: {min: -50, max: 50}}})")
    page.evaluate("() => drawSplitChart()")
    assert _split_window(page)["span"] == pytest.approx(
        page.evaluate("() => state.splitChartHit.dataMax - state.splitChartHit.dataMin")
    )
    assert page.pageerrors == []



# ---------------------------------------------------------------------------
# Layout compactness, and the edge score badge sitting on the line it labels
# ---------------------------------------------------------------------------


# Everything the compactness tests need, measured off the same two arrays the
# canvas draws from. state.splitLinks[i].left/.right are the *same objects* as
# entries in state.visibleLayout (applyComputedLayout builds both from one map),
# so identity comparison is enough to skip a link's own endpoints. The geometry
# comes from the viewer's own pointSegmentDistance rather than a reimplementation.
_LAYOUT_METRICS = """() => {
    const med = values => {
        const sorted = [...values].sort((left, right) => left - right);
        return sorted.length ? sorted[Math.floor(sorted.length / 2)] : 0;
    };
    const items = state.visibleLayout;
    const links = state.splitLinks;
    const drawn = links.map(link =>
        Math.hypot(link.right.x - link.left.x, link.right.y - link.left.y)
        - link.left.radius - link.right.radius);

    let minBubbleSlack = Infinity;
    for (let a = 0; a < items.length; a++) {
        for (let b = a + 1; b < items.length; b++) {
            const left = items[a], right = items[b];
            minBubbleSlack = Math.min(minBubbleSlack,
                Math.hypot(right.x - left.x, right.y - left.y) - left.radius - right.radius);
        }
    }

    let minEdgeClearance = Infinity;
    links.forEach(link => items.forEach(item => {
        if (item === link.left || item === link.right) { return; }
        const hit = pointSegmentDistance(item, link.left, link.right);
        if (hit.t <= 0.03 || hit.t >= 0.97) { return; }
        minEdgeClearance = Math.min(minEdgeClearance, hit.distance - item.radius);
    }));

    const side = (p, q, r) => ((q.y - p.y) * (r.x - q.x)) - ((q.x - p.x) * (r.y - q.y));
    const crosses = (a, b, c, d) =>
        side(a, b, c) * side(a, b, d) < 0 && side(c, d, a) * side(c, d, b) < 0;
    let crossings = 0;
    for (let i = 0; i < links.length; i++) {
        for (let j = i + 1; j < links.length; j++) {
            const first = links[i], second = links[j];
            if (first.left === second.left || first.left === second.right
                || first.right === second.left || first.right === second.right) { continue; }
            if (crosses(first.left, first.right, second.left, second.right)) { crossings++; }
        }
    }

    return {
        ratio: med(drawn) / med(items.map(item => item.radius * 2)),
        minDrawnEdge: Math.min(...drawn),
        minBubbleSlack,
        minEdgeClearance: minEdgeClearance === Infinity ? null : minEdgeClearance,
        crossings,
        links: links.length,
    };
}"""


@pytest.mark.parametrize("layout", ["tree", "force"])
def test_layout_is_compact_without_collapsing(spindle_page, layout):
    """Edges should be sized by the bubbles they join, not by the biggest bubble around.

    The measurable is median drawn edge length over median bubble diameter. Both
    layouts used to derive their spacing from a component-wide maximum radius, which
    on this fixture (six-node blobs beside singletons) put it at 8.5 for the tree and
    5.8 for the force layout; they now sit at 3.6 and 2.4. The ceiling is set between
    the two so it fails loudly if that spacing regresses.

    The other four assertions are the half that matters just as much: without them
    "compact" would also be satisfied by collapsing the drawing into a pile.
    """
    page = spindle_page
    _spindle_state(page, layout=layout)
    page.wait_for_function(
        "() => !state.layoutComputing && state.visibleLayout.length > 0"
    )
    metrics = page.evaluate(_LAYOUT_METRICS)

    assert metrics["links"] > 0
    assert metrics["ratio"] < 5.0
    # Edges stay visible rather than shrinking to stubs between touching bubbles.
    assert metrics["minDrawnEdge"] >= 8
    # No bubble overlaps and no edge driven through a third bubble -- the two
    # guarantees refineLayoutGeometry exists to provide.
    assert metrics["minBubbleSlack"] >= -0.5
    assert metrics["minEdgeClearance"] is None or metrics["minEdgeClearance"] >= -0.5
    # This fixture lays out cleanly today; compaction must not buy itself crossings.
    assert metrics["crossings"] == 0
    assert page.pageerrors == []


def test_radial_seed_rings_stay_concentric_around_a_huge_cluster(spindle_page):
    """One radius per depth, even when a single bubble dwarfs the rest.

    radialTreeSeed lays subtrees out as disjoint annular sectors, which is what keeps
    sibling subtrees from interleaving. Sizing each parent->child step individually
    breaks that: nodes on the same depth land at different radii, the sectors stop
    being disjoint, and subtrees tangle. On a real network with one 14k-member cluster
    and a few hundred singletons hanging off it that produced 216 edge crossings and 11
    overlapping bubbles, which on screen was the cluster's neighbours piled up in
    stacked rows off to one side.

    The seed is a pure function of the topology, so this drives it directly with a
    synthetic hub rather than needing a fixture big enough to grow a 14k cluster.
    """
    result = spindle_page.evaluate(
        """() => {
            // One huge hub, a tail so the tree centre is not the hub itself, and a
            // few hundred singletons hanging off the hub.
            const hierarchyNodes = {};
            const adjacency = new Map();
            const ids = [];
            const add = (id, size) => {
                ids.push(id);
                hierarchyNodes[id] = {size, leaf_start: id};
                adjacency.set(id, []);
            };
            const link = (a, b) => { adjacency.get(a).push(b); adjacency.get(b).push(a); };
            add(0, 14000);
            for (let i = 1; i <= 300; i++) { add(i, 1); link(0, i); }
            let prev = 0;
            for (let t = 0; t < 8; t++) { const id = 400 + t; add(id, 1); link(prev, id); prev = id; }

            const seed = radialTreeSeed(ids, adjacency, hierarchyNodes, 90);
            // Recover each node's depth from the same rooting the seed used.
            const rootId = treeCenter(ids, adjacency, hierarchyNodes);
            const tree = rootedTree(rootId, adjacency, hierarchyNodes);
            const depth = new Map([[rootId, 0]]);
            tree.order.forEach(n => (tree.children.get(n) || []).forEach(
                c => depth.set(c, (depth.get(n) || 0) + 1)));

            // Radius from the root, bucketed by depth: every bucket must be a single value.
            const byDepth = new Map();
            ids.forEach(id => {
                const p = seed.positionById.get(id);
                const r = Math.hypot(p.x, p.y);
                const d = depth.get(id) || 0;
                if (!byDepth.has(d)) { byDepth.set(d, []); }
                byDepth.get(d).push(r);
            });
            let worstSpread = 0;
            byDepth.forEach(radii => {
                worstSpread = Math.max(worstSpread, Math.max(...radii) - Math.min(...radii));
            });

            // And the ring beyond the hub has to clear the hub's own rim.
            const hub = seed.positionById.get(0);
            const kids = (tree.children.get(0) || []);
            const clearance = Math.min(...kids.map(k => {
                const p = seed.positionById.get(k);
                return Math.hypot(p.x - hub.x, p.y - hub.y) - hub.radius - p.radius;
            }));
            return {worstSpread, clearance, hubRadius: hub.radius, kids: kids.length};
        }"""
    )
    assert result["kids"] > 250
    assert result["hubRadius"] > 400          # the hub really does dwarf the singletons
    # Same depth, same radius (the jitter the seed adds is sub-pixel).
    assert result["worstSpread"] < 2
    # Nothing on the next ring is buried inside the hub.
    assert result["clearance"] > 0
    assert spindle_page.pageerrors == []


def test_radial_seed_ring_has_room_for_everything_on_it(spindle_page):
    """A ring's circumference has to fit the bubbles standing on it.

    Sizing rings by radius alone is not enough: a cluster with hundreds of
    neighbours got them seeded shoulder to shoulder on a circle with no room, and
    the collision forces then blew the ring apart. On the reported network -- main
    cluster of 13,841 with 351 neighbours -- every one of those 351 ended up in a
    single beam off one side, with the median neighbour 6,800px from the cluster's
    rim instead of ~70.
    """
    result = spindle_page.evaluate(
        """() => {
            const hierarchyNodes = {};
            const adjacency = new Map();
            const ids = [];
            const add = (id, size) => {
                ids.push(id);
                hierarchyNodes[id] = {size, leaf_start: id};
                adjacency.set(id, []);
            };
            const link = (a, b) => { adjacency.get(a).push(b); adjacency.get(b).push(a); };
            // The reported shape: one dominant cluster carrying a large fan of singletons.
            add(0, 13841);
            for (let i = 1; i <= 351; i++) { add(i, 1); link(0, i); }

            const seed = radialTreeSeed(ids, adjacency, hierarchyNodes, 90);
            const hub = seed.positionById.get(0);
            const kids = ids.filter(id => id !== 0);
            // Arc available per neighbour on their shared ring, against what they occupy.
            const ringRadius = Math.hypot(
                seed.positionById.get(kids[0]).x - hub.x,
                seed.positionById.get(kids[0]).y - hub.y);
            const circumference = 2 * Math.PI * ringRadius;
            const occupied = kids.reduce(
                (sum, k) => sum + (2 * seed.positionById.get(k).radius), 0);
            return {ringRadius, circumference, occupied, kids: kids.length,
                    hubRadius: hub.radius};
        }"""
    )
    assert result["kids"] == 351
    # The ring is measured from the hub centre, so it must clear the hub itself...
    assert result["ringRadius"] > result["hubRadius"]
    # ...and still leave room for every bubble standing on it, with clearance between.
    assert result["circumference"] > result["occupied"] * 1.5
    assert spindle_page.pageerrors == []


def test_worker_and_main_thread_force_layouts_agree(crowded_page):
    """The Worker and the fallback must draw the same network the same way.

    They are one source (`_layout_core_js`) sharing one tuning factory, but that is
    exactly the invariant worth pinning: the layout used to exist as two hand-kept
    copies, and nothing caught it when they drifted. A viewer whose picture depends
    on whether `new Worker(...)` succeeded is a bug nobody would think to look for.

    This is an end-to-end check that the fallback path runs and lands in the same
    place, not a tripwire for the two copies drifting again: the drift that actually
    happened was a floating-point difference in how the anchor and gravity terms were
    summed, and neither this fixture nor the spindle one is large enough for that to
    amplify past 1e-6 -- both paths agreed here even while the sources differed. The
    guard for "there is only one copy" is structural and lives in
    test_build_ssn_viewer.py::test_layout_core_is_shared_by_the_worker_and_the_page.
    """
    page = crowded_page
    page.select_option("#layout-algorithm", "force")
    page.wait_for_function(
        "() => !state.layoutComputing && state.visibleLayout.length > 0"
    )
    positions = """() => state.visibleLayout
        .map(item => [item.componentId, item.x, item.y])
        .sort((left, right) => left[0] - right[0])"""
    with_worker = page.evaluate(positions)
    assert page.evaluate("() => state.layoutWorker !== null")

    page.evaluate("() => { state.layoutWorker = null; state.layoutCache.clear(); }")
    page.evaluate("() => drawClusterView(false)")
    page.wait_for_function(
        "() => !state.layoutComputing && state.visibleLayout.length > 0"
    )
    without_worker = page.evaluate(positions)

    assert len(with_worker) == len(without_worker)
    for (left_id, left_x, left_y), (right_id, right_x, right_y) in zip(
        with_worker, without_worker
    ):
        assert left_id == right_id
        assert left_x == pytest.approx(right_x, abs=1e-6)
        assert left_y == pytest.approx(right_y, abs=1e-6)
    assert page.pageerrors == []


def test_force_layout_is_deterministic(spindle_page):
    """Same network, same settings, same picture -- twice."""
    page = spindle_page
    _spindle_state(page, layout="force")
    page.wait_for_function(
        "() => !state.layoutComputing && state.visibleLayout.length > 0"
    )
    positions = "() => state.visibleLayout.map(item => [item.componentId, item.x, item.y])"
    first = page.evaluate(positions)
    page.evaluate("() => { state.layoutCache.clear(); drawClusterView(false); }")
    page.wait_for_function(
        "() => !state.layoutComputing && state.visibleLayout.length > 0"
    )
    assert page.evaluate(positions) == first
    assert page.pageerrors == []


@pytest.mark.parametrize("layout", ["tree", "force"])
def test_edge_score_badge_sits_on_the_drawn_link(spindle_page, layout):
    """The badge anchor has to land on the polyline the viewer actually strokes.

    It used to be the midpoint of the two bubble *centers*, which is off the drawn
    line by (rightRadius - leftRadius) / 2 once each end is trimmed to its own
    bubble -- and in the tree layout, where the link is drawn as a three-segment
    elbow, generally nowhere near it.
    """
    page = spindle_page
    _spindle_state(page, layout=layout)
    page.check("#show-edge-scores")
    page.wait_for_function(
        "() => !state.layoutComputing && state.visibleLayout.length > 0"
    )
    worst = page.evaluate(
        """() => {
            let worstDistance = 0;
            let worstLengthError = 0;
            state.splitLinks.forEach(link => {
                const anchor = linkLabelAnchor(link);
                const segments = renderedLinkSegments(link);
                const nearest = Math.min(...segments.map(segment => pointSegmentDistance(
                    anchor,
                    {x: segment.startX, y: segment.startY},
                    {x: segment.endX, y: segment.endY}).distance));
                worstDistance = Math.max(worstDistance, nearest);
                // `length` must be the drawn length, since that is what gates the badge.
                const drawn = segments.reduce((sum, segment) => sum + Math.hypot(
                    segment.endX - segment.startX, segment.endY - segment.startY), 0);
                worstLengthError = Math.max(worstLengthError, Math.abs(anchor.length - drawn));
            });
            return {worstDistance, worstLengthError, links: state.splitLinks.length};
        }"""
    )
    assert worst["links"] > 0
    assert worst["worstDistance"] < 0.5
    assert worst["worstLengthError"] < 1e-6
    assert page.pageerrors == []


def test_edge_score_label_is_gated_on_the_drawn_length_not_the_centers(spindle_page):
    """Two big bubbles nearly touching are far apart center to center, but show
    almost no edge. Gating on the center distance let the badge overflow onto them."""
    page = spindle_page
    _spindle_state(page, layout="force")
    verdict = page.evaluate(
        """() => {
            const link = {left: {x: 0, y: 0, radius: 140},
                          right: {x: 300, y: 0, radius: 140}};
            const anchor = linkLabelAnchor(link);
            return {
                centerDistance: 300,
                drawnLength: anchor.length,
                fitsOnDrawn: edgeScoreLabelFits('9.87', anchor.length),
                fitsOnCenters: edgeScoreLabelFits('9.87', 300),
            };
        }"""
    )
    assert verdict["drawnLength"] == pytest.approx(20)
    assert verdict["fitsOnCenters"] is True   # what the old rule saw
    assert verdict["fitsOnDrawn"] is False    # what is actually on screen
    assert page.pageerrors == []
