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
    # The moving sum is computed in Python and shipped in the bundle, not derived here.
    assert page.evaluate("() => Array.isArray(state.bundle.graph.merge_moving_sum.y)")
    # Stems are drawn at the largest single merge, so every event carries the new keys.
    assert page.evaluate(
        "() => state.bundle.graph.merge_event_series.every("
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
    ):
        page.click(f"#{checkbox_id}")
        page.click(f"#{checkbox_id}")
    assert page.pageerrors == []


# --- Custom color palette: picker, color-table TSV, and legend export ---


def _open_color_picker(page):
    page.click("#customize-colors")
    page.wait_for_selector("#color-picker-overlay:not([hidden])")


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
    """The "Domainator distinct" palette assigns the colors build_ssn.py would."""
    from domainator.utils import get_palette

    page = meta_page
    _open_color_picker(page)
    page.wait_for_selector("#color-picker-discrete:not([hidden])")
    before = _canvas_snapshot(page)
    page.select_option("#color-palette", "domainator")
    _wait_for_canvas_change(page, before)

    expected = get_palette(pd.Series(["alpha", "alpha", "beta", "beta", "gamma", "gamma"]))
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
    """A short palette cycles across many values, and the default palette comes back."""
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

    page.select_option("#color-palette", "__default__")
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


def _save_session(page, out_dir, filename="session.ssnv"):
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
    """Decompress a saved .ssnv session and return the parsed JSON."""
    import gzip
    import json

    raw = open(path, "rb").read()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return json.loads(raw.decode("utf-8"))


def test_saved_session_is_a_valid_v4_bundle(meta_page, tmp_path):
    """Saving produces a bundle Python readers accept, carrying an app_state."""
    from domainator.ssn_bundle import (
        SSN_VIEWER_BUNDLE_VERSION,
        load_bundle,
    )

    page = meta_page
    saved = _save_session(page, tmp_path)

    bundle = load_bundle(saved)  # would raise on a bad format/version
    assert bundle["version"] == SSN_VIEWER_BUNDLE_VERSION == 4
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
    edited = tmp_path / "edited.ssnv"
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
    edited = tmp_path / "missing.ssnv"
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
    edited = tmp_path / "v99.ssnv"
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

    saved = _save_session(page, tmp_path, "presets.ssnv")
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
    """An edit must invalidate the color cache, not just the table."""
    page = meta_page
    before = page.evaluate("() => state.nodeColorCache.slice()")
    _edit_cell(page, 0, "family", "a_brand_new_family")
    after = page.evaluate("() => state.nodeColorCache.slice()")
    assert after[0] != before[0]
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
    saved = _save_session(page, tmp_path, "panels.ssnv")
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

    saved = _save_session(page, tmp_path, "edits.ssnv")
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

    bundle = load_bundle(_save_session(page, tmp_path, "annotated.ssnv"))
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

    saved = _save_session(page, tmp_path, "colors.ssnv")
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
    saved = _save_session(page, tmp_path, "renamed.ssnv")
    bundle = _read_session_bundle(saved)
    bundle["name"] = "Some Other Network"
    renamed = tmp_path / "renamed_bundle.ssnv"
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

    saved = _save_session(page, tmp_path, "renamed_col.ssnv")
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
    # Cluster ids are labels, not magnitudes, so they colour as discrete categories.
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


def _save_extraction(page, out_dir, filename="extraction.ssnv"):
    with page.expect_download() as download_info:
        page.click("#save-extraction")
    target = out_dir / filename
    download_info.value.save_as(str(target))
    return target


def test_extraction_of_everything_reproduces_the_python_bundle(meta_page, tmp_path):
    """The JS hierarchy/merge-series port must agree with ssn_hierarchy.py.

    Selecting every node makes the extraction a rebuild of the whole network, so
    it can be compared field-for-field against what build_ssn_viewer.py wrote.
    """
    page = meta_page
    original = page.evaluate("() => state.bundle.graph")
    _select_nodes(page, list(range(6)))

    extracted = _read_session_bundle(_save_extraction(page, tmp_path))["graph"]

    assert extracted["nodes"] == original["nodes"]
    assert extracted["mst_edges"] == original["mst_edges"]
    assert extracted["hierarchy"] == original["hierarchy"]
    assert extracted["cluster_count_by_threshold"] == original["cluster_count_by_threshold"]
    assert extracted["merge_impact_metric"] == original["merge_impact_metric"]

    assert len(extracted["merge_event_series"]) == len(original["merge_event_series"])
    for got, want in zip(extracted["merge_event_series"], original["merge_event_series"]):
        assert got == pytest.approx(want, rel=1e-9) if not isinstance(want, dict) else True
        for field in ("edge_index", "threshold_value", "threshold_to", "threshold_from",
                      "merge_impact", "largest_merge", "merge_count", "merge_size_counts",
                      "summary_row_from", "summary_row_to"):
            assert got[field] == want[field], field
        for field in ("delta_largest", "delta_avg_non_singleton", "threshold_from_value"):
            assert got[field] == pytest.approx(want[field]), field

    assert extracted["slider_stops"] == original["slider_stops"]
    assert extracted["merge_moving_sum"]["x"] == pytest.approx(
        original["merge_moving_sum"]["x"], rel=1e-9
    )
    assert extracted["merge_moving_sum"]["y"] == pytest.approx(
        original["merge_moving_sum"]["y"], rel=1e-9
    )
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
    saved = _save_extraction(page, tmp_path, "abc.ssnv")

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
    bundle = _read_session_bundle(_save_extraction(page, tmp_path, "def.ssnv"))

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
    app_state = _read_session_bundle(_save_extraction(page, tmp_path, "state.ssnv"))["app_state"]

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
    saved = _save_extraction(page, tmp_path, "clean.ssnv")

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

    extracted = _read_session_bundle(_save_extraction(page, tmp_path, "all120.ssnv"))["graph"]

    assert extracted["hierarchy"] == original["hierarchy"]
    assert extracted["mst_edges"] == original["mst_edges"]
    assert extracted["cluster_count_by_threshold"] == original["cluster_count_by_threshold"]
    assert extracted["slider_stops"] == original["slider_stops"]
    assert extracted["merge_event_series"] == original["merge_event_series"]
    assert extracted["merge_moving_sum"] == original["merge_moving_sum"]
    # One tie group covering all 119 merges.
    assert len(original["merge_event_series"]) == 1
    assert original["merge_event_series"][0]["merge_count"] == 119
    assert page.pageerrors == []


def test_extraction_of_a_contiguous_run_is_connected(many_cat_page, tmp_path):
    """A contiguous stretch of the path is MST-connected; a gapped one is not."""
    page = many_cat_page
    _select_nodes(page, list(range(10, 30)))
    bundle = _read_session_bundle(_save_extraction(page, tmp_path, "run.ssnv"))
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
    assert len(original["merge_event_series"]) > 40      # many distinct thresholds
    _select_nodes(page, list(range(60)))

    extracted = _read_session_bundle(_save_extraction(page, tmp_path, "all60.ssnv"))["graph"]

    assert extracted["hierarchy"] == original["hierarchy"]
    assert extracted["cluster_count_by_threshold"] == original["cluster_count_by_threshold"]
    assert extracted["slider_stops"] == original["slider_stops"]
    assert len(extracted["merge_event_series"]) == len(original["merge_event_series"])
    for got, want in zip(extracted["merge_event_series"], original["merge_event_series"]):
        assert got["edge_index"] == want["edge_index"]
        assert got["merge_size_counts"] == want["merge_size_counts"]
        assert got["merge_count"] == want["merge_count"]
        assert got["threshold_value"] == pytest.approx(want["threshold_value"])
        assert got["merge_impact"] == pytest.approx(want["merge_impact"])
        assert got["largest_merge"] == pytest.approx(want["largest_merge"])
        assert got["delta_largest"] == pytest.approx(want["delta_largest"])
        assert got["delta_avg_non_singleton"] == pytest.approx(want["delta_avg_non_singleton"])
    assert extracted["merge_moving_sum"]["window"] == pytest.approx(
        original["merge_moving_sum"]["window"]
    )
    assert extracted["merge_moving_sum"]["x"] == pytest.approx(original["merge_moving_sum"]["x"])
    assert extracted["merge_moving_sum"]["y"] == pytest.approx(original["merge_moving_sum"]["y"])
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
    bundle = _read_session_bundle(_save_extraction(page, tmp_path, "multi.ssnv"))
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
    saved = _save_extraction(page, tmp_path, "floor.ssnv")

    graph = _read_session_bundle(saved)["graph"]
    stops = graph["slider_stops"]
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
