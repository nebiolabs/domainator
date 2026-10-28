"""Browser tests of the HTML report pages (report_html.py), in headless Chromium via Playwright.

test_report_html.py checks the data that goes into the pages; these check that the pages work
with it. Like test_ssn_viewer_browser.py they need the optional `browser` extra:

    pip install -e ".[test,browser]"
    playwright install chromium

Every page is loaded with all network requests refused: the reports must work offline.
"""
import json

import pytest

from domainator import table_writers
from helpers import report_data

sync_api = pytest.importorskip("playwright.sync_api")


def _table_page(path, rows=60, max_height=None):
    columns, types = ["name", "length", "score"], ["str", "int", "float"]
    with open(path, "w") as handle:
        writer = table_writers.HTMLTableWriter(columns, types, handle, max_height)
        writer.write_header()
        for i in range(rows):
            # Every fifth row has no score, so missing values are exercised.
            writer.write_row([f"seq_{i:02d}", 100 + (i * 37) % 60, None if i % 5 == 4 else i * 1.25])
        writer.write_footer()
    return path


@pytest.fixture
def open_page(tmp_path):
    """open_page(html_path) -> a loaded Playwright page, with page.pageerrors collecting
    uncaught JS errors and page.requests every request that was refused."""
    with sync_api.sync_playwright() as p:
        try:
            browser = p.chromium.launch()
        except sync_api.Error as exc:
            pytest.skip(f"Chromium not available for Playwright; run 'playwright install chromium'. ({exc})")
        context = browser.new_context(accept_downloads=True)
        context.grant_permissions(["clipboard-read", "clipboard-write"])

        def refuse(route):
            if not route.request.url.startswith("file:"):
                page.requests.append(route.request.url)
                route.abort()
            else:
                route.continue_()

        def opener(html_path):
            nonlocal page
            page = context.new_page()
            page.pageerrors = []
            page.requests = []
            page.on("pageerror", lambda exc: page.pageerrors.append(str(exc)))
            page.route("**/*", refuse)
            page.goto(html_path.as_uri())
            return page

        page = None
        try:
            yield opener
        finally:
            browser.close()


def _column(page, index):
    return page.eval_on_selector_all(f".report-table tbody tr td:nth-child({index})", "cells => cells.map(c => c.textContent)")


def _header(page):
    return page.eval_on_selector_all(".report-table thead th .column-title span:first-child", "spans => spans.map(s => s.textContent)")


def test_table_loads_offline_and_pages(open_page, tmp_path):
    page = open_page(_table_page(tmp_path / "table.html"))
    page.wait_for_selector(".report-table tbody tr")
    assert page.pageerrors == [] and page.requests == []
    assert _header(page) == ["name", "length", "score"]
    assert page.text_content("#table-counter") == "Showing 1-25 of 60 rows"
    # Floats are shown to one decimal; missing values are blank.
    assert _column(page, 3)[:5] == ["0.0", "1.2", "2.5", "3.8", ""]

    page.click("text=Next ›")
    assert page.text_content("#table-counter") == "Showing 26-50 of 60 rows"
    assert _column(page, 1)[0] == "seq_25"
    page.select_option("select", "100")
    assert page.text_content("#table-counter") == "Showing 1-60 of 60 rows"
    assert page.is_disabled("text=Next ›")


def test_table_sorts_and_filters(open_page, tmp_path):
    page = open_page(_table_page(tmp_path / "table.html"))
    page.wait_for_selector(".report-table tbody tr")

    page.click(".report-table thead th:nth-child(3) .column-title")
    scores = _column(page, 3)
    assert scores[:3] == ["0.0", "1.2", "2.5"]
    page.select_option("select", "100")
    assert _column(page, 3)[-12:] == [""] * 12  # missing values sort last
    page.click(".report-table thead th:nth-child(3) .column-title")
    descending = _column(page, 3)
    assert descending[0] == "72.5" and descending[-1] == ""

    page.fill("input[aria-label='Filter score']", ">= 70")
    assert page.text_content("#table-counter") == "Showing 1-3 of 3 rows (filtered from 60)"
    page.fill("input[aria-label='Filter score']", "")
    page.fill("input[aria-label='Filter name']", "SEQ_1")
    # Still sorted by score, descending; seq_14 and seq_19 have none, so they come last, in order.
    assert _column(page, 1) == [f"seq_{i}" for i in (18, 17, 16, 15, 13, 12, 11, 10, 14, 19)]
    assert page.pageerrors == []


def test_table_columns_move_and_downloads_keep_full_values(open_page, tmp_path):
    page = open_page(_table_page(tmp_path / "table.html", rows=6))
    page.wait_for_selector(".report-table tbody tr")
    page.drag_and_drop(".report-table thead th:nth-child(3)", ".report-table thead th:nth-child(1)")
    assert _header(page) == ["score", "name", "length"]

    page.fill("input[aria-label='Filter length']", "< 200")
    with page.expect_download() as download:
        page.click("#download-json")
    saved = json.loads(open(download.value.path()).read())
    assert saved[1] == {"score": 1.25, "name": "seq_01", "length": 137}
    assert saved[4] == {"score": None, "name": "seq_04", "length": 128}
    assert list(saved[0]) == ["score", "name", "length"]

    with page.expect_download() as download:
        page.click("#download-csv")
    assert download.value.suggested_filename == "enum_report.csv"
    lines = open(download.value.path()).read().splitlines()
    assert lines[:3] == ["score,name,length", "0.0,seq_00,100", "1.25,seq_01,137"]
    assert lines[5] == ",seq_04,128"
    assert page.pageerrors == []


def test_table_max_height(open_page, tmp_path):
    page = open_page(_table_page(tmp_path / "table.html", max_height=120))
    page.wait_for_selector(".report-table tbody tr")
    assert page.eval_on_selector(".report-table-wrap", "el => el.getBoundingClientRect().height") <= 120
    # The header stays in view while the rows scroll under it.
    page.eval_on_selector(".report-table-wrap", "el => { el.scrollTop = 200; }")
    top = page.eval_on_selector(".report-table-wrap", "el => el.getBoundingClientRect().top")
    header_top = page.eval_on_selector(".report-table thead th", "el => el.getBoundingClientRect().top")
    assert abs(header_top - top) < 2


def _width(page, selector):
    return page.eval_on_selector(selector, "el => el.getBoundingClientRect().width")


def test_table_pager_goes_to_first_last_and_typed_pages(open_page, tmp_path):
    page = open_page(_table_page(tmp_path / "table.html", rows=130))
    page.wait_for_selector(".report-table tbody tr")
    assert page.input_value("#page-input") == "1"
    assert page.text_content("#page-count") == "of 6"
    assert page.is_disabled("#page-first") and page.is_disabled("#page-prev")

    page.click("#page-last")
    assert page.text_content("#table-counter") == "Showing 126-130 of 130 rows"
    assert page.input_value("#page-input") == "6"
    assert page.is_disabled("#page-next") and page.is_disabled("#page-last")
    page.click("#page-first")
    assert page.text_content("#table-counter") == "Showing 1-25 of 130 rows"

    page.fill("#page-input", "3")
    page.press("#page-input", "Enter")
    assert page.text_content("#table-counter") == "Showing 51-75 of 130 rows"
    assert _column(page, 1)[0] == "seq_50"
    # Out of range goes to the nearest page; something that is not a number changes nothing.
    page.fill("#page-input", "99")
    page.press("#page-input", "Enter")
    assert page.input_value("#page-input") == "6"
    page.fill("#page-input", "abc")
    page.press("#page-input", "Enter")
    assert page.input_value("#page-input") == "6"
    # A filter that leaves fewer pages pulls the page back to one that exists.
    page.fill("input[aria-label='Filter name']", "seq_1")
    assert page.text_content("#page-count") == "of 2"
    page.select_option("select", "50")
    assert page.text_content("#page-count") == "of 1"
    assert page.pageerrors == []


def test_table_columns_resize(open_page, tmp_path):
    page = open_page(_table_page(tmp_path / "table.html"))
    page.wait_for_selector(".report-table tbody tr")
    name_header = ".report-table thead th:nth-child(1)"
    assert _width(page, name_header) == pytest.approx(150, abs=1)

    handle = page.locator(f"{name_header} .resize-handle")
    box = handle.bounding_box()
    page.mouse.move(box["x"] + box["width"] / 2, box["y"] + box["height"] / 2)
    page.mouse.down()
    page.mouse.move(box["x"] + 100, box["y"] + box["height"] / 2, steps=5)
    page.mouse.up()
    assert _width(page, name_header) == pytest.approx(250, abs=4)
    # The cells follow the header, and the other columns keep their widths.
    assert _width(page, ".report-table tbody tr td:nth-child(1)") == pytest.approx(_width(page, name_header), abs=1)
    assert _width(page, ".report-table thead th:nth-child(2)") == pytest.approx(150, abs=1)
    # Resizing neither sorted nor moved the column.
    assert _header(page) == ["name", "length", "score"]
    assert page.text_content(f"{name_header} .sort-mark") == ""

    # Double-click fits the column to what it shows; the arrow keys nudge it.
    handle.dblclick()
    fitted = _width(page, name_header)
    assert 50 < fitted < 120
    assert page.eval_on_selector(".report-table tbody tr td:nth-child(1)", "el => el.scrollWidth <= el.clientWidth")
    handle.focus()
    page.keyboard.press("ArrowRight")
    assert _width(page, name_header) == pytest.approx(fitted + 10, abs=1)

    # Widths survive sorting, which rebuilds the header.
    page.click(".report-table thead th:nth-child(3) .column-title")
    assert _width(page, name_header) == pytest.approx(fitted + 10, abs=1)
    assert page.pageerrors == []


def test_table_context_menu_copies_the_cell_value(open_page, tmp_path):
    page = open_page(_table_page(tmp_path / "table.html"))
    page.wait_for_selector(".report-table tbody tr")
    # Row 1's score is 1.25, shown as "1.2"; the copy is the value itself.
    page.click(".report-table tbody tr:nth-child(2) td:nth-child(3)", button="right")
    assert page.is_visible(".context-menu")
    page.click("text=Copy cell value")
    assert not page.is_visible(".context-menu")
    assert page.evaluate("navigator.clipboard.readText()") == "1.25"
    assert page.text_content("#table-copy-status") == "Copied."

    page.click(".report-table tbody tr:nth-child(1) td:nth-child(1)", button="right")
    page.click("text=Copy cell value")
    assert page.evaluate("navigator.clipboard.readText()") == "seq_00"

    # After sorting, the copy comes from the row now shown there; a missing value copies as "".
    page.click(".report-table thead th:nth-child(3) .column-title")
    page.click(".report-table thead th:nth-child(3) .column-title")
    page.click(".report-table tbody tr:nth-child(1) td:nth-child(3)", button="right")
    page.click("text=Copy cell value")
    assert page.evaluate("navigator.clipboard.readText()") == "72.5"

    # Escape and a click elsewhere close the menu without copying.
    page.click(".report-table tbody tr:nth-child(1) td:nth-child(1)", button="right")
    page.keyboard.press("Escape")
    assert not page.is_visible(".context-menu")
    page.click(".report-table tbody tr:nth-child(1) td:nth-child(1)", button="right")
    page.mouse.click(5, 5)
    assert not page.is_visible(".context-menu")
    assert page.evaluate("navigator.clipboard.readText()") == "72.5"
    assert page.pageerrors == []


def _summary_page(path, *extra):
    from domainator import summary_report
    summary_report.main(["-i", *extra, "--html_out", str(path)])
    return path


def test_summary_page_draws_every_section(open_page, tmp_path, shared_datadir):
    page = open_page(_summary_page(tmp_path / "summary.html", str(shared_datadir / "FeSOD_20_pfam.gb"), "--domains", "Sod_Fe_C", "Sod_Fe_N"))
    page.wait_for_selector("#domain-table .report-table tbody tr")
    assert page.pageerrors == [] and page.requests == []
    stats = page.eval_on_selector_all(".stat-grid div", "items => items.map(i => [i.querySelector('dt').textContent, i.querySelector('dd').textContent])")
    assert ["contigs", "20"] in stats
    # One bar per bin, and the bins hold every contig.
    bins = page.evaluate("window.lengthHistogram.bins()")
    bars = page.eval_on_selector_all("#length-histogram svg rect[data-tip]", "bars => bars.map(b => b.dataset.tip)")
    assert len(bars) == len(bins["counts"]) and sum(bins["counts"]) == 20
    page.hover("#length-histogram svg rect[data-tip]:not([fill-opacity='0'])")
    assert "contigs" in page.text_content(".chart-tooltip")
    # The two tables keep their controls apart.
    assert page.text_content("#domains-table-counter").startswith("Showing 1-")
    assert page.text_content("#cooccurrence-table-counter") == "Showing 1-2 of 2 rows"
    cooccurrence = page.eval_on_selector_all("#cooccurrence-table thead th .column-title span:first-child", "s => s.map(x => x.textContent)")
    assert cooccurrence[:3] == ["domain", "Sod_Fe_C", "Sod_Fe_N"] or cooccurrence[:3] == ["domain", "Sod_Fe_N", "Sod_Fe_C"]
    with page.expect_download() as download:
        page.click("#lengths-download-svg")
    assert download.value.suggested_filename == "summary_report_lengths.svg"
    assert open(download.value.path()).read().count("<rect") == len(bars)
    assert page.query_selector("#taxonomy-chart") is None


def _selected_text(page):
    return page.evaluate("window.getSelection().toString()")


def test_summary_page_histogram_bins_by_the_integer_rule(open_page, tmp_path, shared_datadir):
    page = open_page(_summary_page(tmp_path / "summary.html", str(shared_datadir / "FeSOD_20_pfam.gb")))
    page.wait_for_selector("#length-histogram svg rect")

    def bins(values, counts=None):
        counts = counts or [1] * len(values)
        return page.evaluate("""([values, counts]) => {
            const n = counts.reduce((a, b) => a + b, 0);
            const lo = values[0], hi = values[values.length - 1];
            return binValues(values, counts, lo, hi + 1, autoBinWidth(lo, hi, n));
        }""", [values, counts])

    assert bins([7]) == {"edges": [7, 8], "counts": [1]}
    # 20 lengths aim for max(6, 4) = 6 bins: ceil(1000 / 6) = 167, rounded up to 200.
    assert bins([1000 + 50 * i for i in range(20)]) == {"edges": [1000, 1200, 1400, 1600, 1800, 2000], "counts": [4, 4, 4, 4, 4]}
    # The first edge is the minimum rounded down to a multiple of the width. 4 lengths aim for
    # max(3, 2) = 3 bins: ceil(46 / 3) = 16, rounded up to 20.
    assert bins([95, 113, 140], [2, 1, 1]) == {"edges": [80, 100, 120, 140, 160], "counts": [2, 1, 0, 1]}
    # Many lengths: at most 50 bins.
    many = bins(list(range(100000)))
    assert len(many["counts"]) <= 50 and sum(many["counts"]) == 100000
    assert page.pageerrors == []


def test_summary_page_histogram_zooms_pans_and_rebins(open_page, tmp_path, shared_datadir):
    page = open_page(_summary_page(tmp_path / "summary.html", str(shared_datadir / "FeSOD_20_pfam.gb")))
    page.wait_for_selector("#length-histogram svg rect")
    full = page.evaluate("window.lengthHistogram.bins()")
    assert page.evaluate("window.lengthHistogram.view()") is None
    assert page.is_disabled("#lengths-reset-zoom")
    svg = page.locator("#length-histogram svg")
    svg.scroll_into_view_if_needed()
    box = svg.bounding_box()
    page.mouse.move(box["x"] + box["width"] / 2, box["y"] + box["height"] / 2)
    page.mouse.wheel(0, -300)
    view = page.evaluate("window.lengthHistogram.view()")
    assert view is not None and view[1] - view[0] < full["edges"][-1] - full["edges"][0]
    assert page.is_enabled("#lengths-reset-zoom")

    # A drag pans, and selects no text on the way.
    x, y = box["x"] + box["width"] / 2, box["y"] + box["height"] / 2
    page.mouse.move(x, y)
    page.mouse.down()
    page.mouse.move(x - 150, y + 60, steps=5)
    page.mouse.move(x - 150, box["y"] + box["height"] + 80, steps=5)
    page.mouse.up()
    panned = page.evaluate("window.lengthHistogram.view()")
    assert panned[0] > view[0] and panned[1] - panned[0] == pytest.approx(view[1] - view[0])
    assert _selected_text(page) == ""

    # A typed bin width holds until it is cleared.
    page.click("#lengths-reset-zoom")
    page.fill("#lengths-bin-width", "500")
    page.press("#lengths-bin-width", "Enter")
    edges = page.evaluate("window.lengthHistogram.bins().edges")
    assert {b - a for a, b in zip(edges, edges[1:])} == {500}
    assert page.text_content("#lengths-bin-status").endswith(", 500 wide")
    assert sum(page.evaluate("window.lengthHistogram.bins().counts")) == 20
    page.fill("#lengths-bin-width", "")
    page.press("#lengths-bin-width", "Enter")
    assert page.evaluate("window.lengthHistogram.bins()") == full
    assert page.pageerrors == []


def test_summary_page_taxonomy_zooms(open_page, tmp_path, shared_datadir):
    page = open_page(_summary_page(tmp_path / "summary.html", str(shared_datadir / "swissprot_CuSOD_subset.fasta"), "--taxonomy", "--ncbi_taxonomy_path", str(shared_datadir / "taxdmp")))
    page.wait_for_selector("#taxonomy-chart svg")
    assert page.pageerrors == [] and page.requests == []
    labels = page.eval_on_selector_all("#taxonomy-chart svg text[visibility='visible']", "t => t.map(x => x.textContent)")
    assert any(label.startswith("All") for label in labels)
    assert any(label.startswith("Bacteria") for label in labels)
    before = page.get_attribute("#taxonomy-chart svg g:nth-of-type(2)", "transform")
    page.evaluate("window.taxonomyChart.zoom('Bacteria')")
    page.wait_for_timeout(800)
    # Bacteria now fills the chart, at the left edge.
    assert page.get_attribute("#taxonomy-chart svg g:nth-of-type(2)", "transform") != before
    assert page.get_attribute("#taxonomy-chart svg g:nth-of-type(2)", "transform") == "translate(0,0)"
    assert page.text_content("#taxa-table-counter").startswith("Showing 1-")
    assert page.pageerrors == []


def test_summary_page_with_no_contigs(open_page, tmp_path, shared_datadir):
    page = open_page(_summary_page(tmp_path / "summary.html", str(shared_datadir / "empty.gb")))
    page.wait_for_selector("#domain-table .report-table")
    assert page.text_content("#length-histogram") == "No contigs."
    assert page.text_content("#domains-table-counter") == "No rows"
    assert page.pageerrors == [] and page.requests == []


def _matrix_page(path, matrix_path, *extra):
    from domainator import matrix_report
    matrix_report.main(["-i", str(matrix_path), "--html_out", str(path), *extra])
    return path


def _matrix_data(path):
    from helpers import report_data
    return report_data(path.read_text())["data"]


def _figures(page):
    return dict(page.eval_on_selector_all(".stat-grid div", "items => items.map(i => [i.querySelector('dt').textContent, i.querySelector('dd').textContent])"))


def test_matrix_page_walks_the_slider_through_every_stop(open_page, tmp_path, shared_datadir):
    path = _matrix_page(tmp_path / "matrix.html", shared_datadir / "scorefull.tsv", "--knn_projection")
    data = _matrix_data(path)
    page = open_page(path)
    page.wait_for_selector("#cluster-bubbles svg circle")
    assert page.pageerrors == [] and page.requests == []
    figures = _figures(page)
    assert figures["nodes"] == "20" and figures["non-zero edges"] == "190"
    # One stem per split event, and the caption counts them.
    events = data["summary"]["split_events"]
    assert len(page.query_selector_all("#split-chart svg g[stroke='#a8d2cf'] line")) == len(events)
    assert page.text_content("#split-event-count").startswith(f"All {len(events)} split events plotted")

    # Each stop shows what its row of the threshold table says: the page's union-find over the
    # MST agrees with Python's cluster counts.
    table = data["threshold_table"]
    knn = data["knn_projection"]
    page.focus("#threshold-slider")
    for i, stop in enumerate(data["slider_stops"]):
        if i > 0:
            page.keyboard.press("ArrowLeft")
        row = stop["threshold_index"]
        figures = _figures(page)
        assert page.text_content("#threshold-number") == stop["threshold_label"]
        assert figures["clusters"] == str(table["clusters"][row] if row >= 0 else 20)
        assert figures["edges kept"] == f'{table["edges"][row] if row >= 0 else 0:,}'
        assert figures["projected MST + Knn edges"] == f'{knn["counts"][row][knn["default_k"] - 1] if row >= 0 else 0:,}'
        assert len(page.query_selector_all("#cluster-bubbles svg circle")) == min(100, int(figures["clusters"]))
    assert page.text_content("#threshold-number") == data["slider_stops"][-1]["threshold_label"]
    assert _figures(page)["clusters"] == "1"
    page.keyboard.press("Home")
    assert page.text_content("#threshold-number") == "∞"

    # The k slider changes the projected MST + kNN edge count.
    page.keyboard.press("End")
    page.eval_on_selector("#mst-knn-k-slider", "s => { s.value = '2'; s.dispatchEvent(new Event('input')); }")
    floor_row = data["slider_stops"][-1]["threshold_index"]
    assert _figures(page)["projected MST + Knn edges"] == f'{knn["counts"][floor_row][1]:,}'
    assert page.text_content("#mst-knn-k-number") == "2"
    assert page.pageerrors == []


def test_matrix_page_split_chart_hover_click_and_download(open_page, tmp_path, shared_datadir):
    path = _matrix_page(tmp_path / "matrix.html", shared_datadir / "scorefull.tsv")
    page = open_page(path)
    page.wait_for_selector("#split-chart svg circle")
    element = page.query_selector("#split-chart svg g[fill='#72b7b2'] circle")
    element.scroll_into_view_if_needed()
    bead = element.bounding_box()
    page.mouse.move(bead["x"] + bead["width"] / 2, bead["y"] + bead["height"] / 2)
    tip = page.text_content(".chart-tooltip")
    assert "split(s) of" in tip and "nodes" in tip and "From:" in tip
    # Clicking an event moves the threshold to it.
    page.mouse.click(bead["x"] + bead["width"] / 2, bead["y"] + bead["height"] / 2)
    assert page.text_content("#threshold-number") in tip
    with page.expect_download() as download:
        page.click("#split-chart-download-svg")
    assert download.value.suggested_filename == "split_events.svg"
    with page.expect_download() as download:
        page.click("#split-chart-download-png")
    assert download.value.suggested_filename == "split_events.png"
    assert page.pageerrors == []


def test_matrix_page_spaces_slider_stops_by_threshold(open_page, tmp_path):
    import numpy as np
    from domainator.data_matrix import DataMatrix
    data = np.array([
        [0, 10, 0, 0, 0],
        [10, 0, 5, 0, 0],
        [0, 5, 0, 4, 0],
        [0, 0, 4, 0, 1],
        [0, 0, 0, 1, 0],
    ], dtype=float)
    DataMatrix(data, list("ABCDE"), list("ABCDE")).write(str(tmp_path / "m.hdf5"), output_type="hdf")
    path = _matrix_page(tmp_path / "matrix.html", tmp_path / "m.hdf5")
    stops = _matrix_data(path)["slider_stops"]
    assert [stop["threshold_label"] for stop in stops] == ["∞", "10.00", "5.00", "4.00", "1.00", "0.91"]
    page = open_page(path)
    page.wait_for_selector("#cluster-bubbles svg circle")
    positions = page.evaluate("window.matrixReport.threshold.positions")
    # Infinity at one end, then each stop as far along as its threshold is down the axis.
    assert positions[0] == 0
    finite = [(stop["threshold_value"], position) for stop, position in zip(stops[1:], positions[1:])]
    (t0, p0), (t1, p1) = finite[0], finite[-1]
    for threshold, position in finite:
        assert position == pytest.approx(p0 + (t0 - threshold) / (t0 - t1) * (p1 - p0), abs=1)
    assert 0 < p0 < p1 < 10000
    assert page.pageerrors == []


def test_matrix_page_threshold_chart_zooms(open_page, tmp_path, shared_datadir):
    page = open_page(_matrix_page(tmp_path / "matrix.html", shared_datadir / "scorefull.tsv"))
    page.wait_for_selector("#clusters-by-threshold svg path")
    ticks = lambda: page.eval_on_selector_all("#clusters-by-threshold svg text[text-anchor='middle']", "t => t.map(x => x.textContent)")
    before = ticks()
    assert page.is_disabled("#threshold-chart-reset-zoom")
    overlay = page.query_selector("#clusters-by-threshold svg rect[data-overlay]")
    overlay.scroll_into_view_if_needed()
    box = overlay.bounding_box()
    y = box["y"] + box["height"] / 2
    page.mouse.move(box["x"] + box["width"] * 0.5, y)
    page.mouse.down()
    page.mouse.move(box["x"] + box["width"] * 0.7, y, steps=5)
    # Wandering off the chart while dragging, as a hand does.
    page.mouse.move(box["x"] + box["width"] * 0.7, box["y"] + box["height"] + 120, steps=5)
    page.mouse.up()
    assert page.is_enabled("#threshold-chart-reset-zoom")
    assert ticks() != before
    # Neither the drag nor the double-click selects the text in and around the chart.
    assert _selected_text(page) == ""
    page.dblclick("#clusters-by-threshold svg rect[data-overlay]")
    assert page.is_disabled("#threshold-chart-reset-zoom")
    assert ticks() == before
    assert _selected_text(page) == ""
    assert page.pageerrors == []


@pytest.mark.parametrize("values", [[[0.0, 0.0, 0.0]] * 3, [[0.0]], [[0.0, 5.0], [5.0, 0.0]]])
def test_matrix_page_with_too_few_edges(open_page, tmp_path, values):
    import numpy as np
    from domainator.data_matrix import DataMatrix
    names = [f"n{i}" for i in range(len(values))]
    DataMatrix(np.array(values), names, names).write(str(tmp_path / "m.hdf5"), output_type="hdf")
    page = open_page(_matrix_page(tmp_path / "matrix.html", tmp_path / "m.hdf5", "--knn_projection"))
    page.wait_for_selector("#cluster-bubbles svg circle")
    assert page.pageerrors == [] and page.requests == []
    assert _figures(page)["clusters"] == str(len(values))
    assert page.is_disabled("#threshold-slider") == (len(values) < 2 or values[0][1] == 0)


def _contig_page(path, *extra):
    from domainator import plot_contigs
    plot_contigs.main(["--html_out", str(path), *extra])
    return path


def _view(page):
    return page.evaluate("window.contigPlot.view.view()")


def test_contig_plot_draws_every_contig_and_feature(open_page, tmp_path, shared_datadir):
    path = _contig_page(tmp_path / "contigs.html", "-i", str(shared_datadir / "MT_nbs.gb"))
    page = open_page(path)
    page.wait_for_selector("svg.contig-plot")
    assert page.pageerrors == [] and page.requests == []
    data = report_data(path.read_text())["data"]
    labels = page.eval_on_selector_all("svg.contig-plot [data-contig] > text", "t => t.map(x => x.textContent)")
    assert labels == [contig["name"] for contig in data["contigs"]]
    parts = sum(len(feature["parts"]) for contig in data["contigs"] for feature in contig["features"])
    assert page.locator("svg.contig-plot [data-feature] path").count() == parts
    # Each contig is below the one before, and its CDSs share the level next to it.
    offsets = page.evaluate("window.contigPlot.contigs.map(c => c.offset)")
    assert offsets == sorted(offsets) and len(set(offsets)) == len(offsets)
    levels = page.evaluate("window.contigPlot.contigs.flatMap(c => c.features.filter(f => f.type === 'CDS').map(f => f.level))")
    assert set(levels) == {0}
    # The widest contig is shrunk to fit, but to no less than a tenth.
    widest = max(contig["length"] for contig in data["contigs"])
    assert _view(page)["k"] == max(0.1, min(1, 1200 / widest))


def test_contig_plot_stacks_overlapping_features(open_page, tmp_path, shared_datadir):
    page = open_page(_contig_page(tmp_path / "contigs.html", "-i", str(shared_datadir / "FeSOD_20_pfam.gb")))
    page.wait_for_selector("svg.contig-plot")
    features = [
        {"name": "a", "type": "Domainator", "parts": [[0, 100, "+"]]},
        {"name": "c", "type": "Domainator", "parts": [[120, 200, "+"]]},
        {"name": "b", "type": "Domainator", "parts": [[50, 150, "+"]]},
        {"name": "y", "type": "CDS", "parts": [[400, 500, "-"]]},
        {"name": "x", "type": "CDS", "parts": [[0, 300, "+"]]},
        {"name": "split", "type": "misc_feature", "parts": [[350, 360, "+"], [0, 10, "+"]]},
    ]
    stacked = page.evaluate("features => stackFeatures(features, ['CDS', 'Domainator']).map(f => [f.name, f.level])", features)
    # By type, then start; CDSs share a level, and every other feature takes the lowest
    # level free of the features it overlaps.
    assert stacked == [["x", 0], ["y", 0], ["a", 1], ["b", 2], ["c", 1], ["split", 2]]


def test_contig_plot_hover_click_zoom_pan_and_download(open_page, tmp_path, shared_datadir):
    page = open_page(_contig_page(tmp_path / "contigs.html", "-i", str(shared_datadir / "FeSOD_20_pfam.gb"),
                                  "--metadata", str(shared_datadir / "metadata_FeSOD_20.tsv"), "--color_by", "group"))
    page.wait_for_selector("svg.contig-plot")
    assert page.pageerrors == [] and page.requests == []
    feature = page.locator("svg.contig-plot [data-contig='0'] [data-feature]").first
    feature.hover()
    assert page.text_content(".chart-tooltip").startswith("Domainator\nstart: ")
    feature.click()
    assert page.text_content("#details .details-heading") == "Domainator"
    assert "evalue" in page.eval_on_selector_all("#details dt", "terms => terms.map(t => t.textContent)")

    page.click("svg.contig-plot [data-contig='0'] > text")
    assert page.text_content("#details .details-heading") == "Contig"
    terms = page.eval_on_selector_all("#details dt", "terms => terms.map(t => t.textContent)")
    assert terms[:3] == ["Name", "Description", "Length"] and terms[-1] == "group"
    # The --color_by value is marked with the contig's color.
    line_color = page.get_attribute("svg.contig-plot [data-contig='0'] > line", "stroke")
    swatch = page.eval_on_selector("#details .color-swatch", "s => getComputedStyle(s).backgroundColor")
    assert swatch == "rgb({}, {}, {})".format(*(int(line_color[i:i + 2], 16) for i in (1, 3, 5)))

    page.click("#zoom-in")
    assert _view(page)["k"] == pytest.approx(1.5)
    page.click("#reset-view")
    assert _view(page) == {"k": 1, "x": 0, "y": 0}
    box = page.locator("svg.contig-plot").bounding_box()
    page.mouse.move(box["x"] + 300, box["y"] + 300)
    page.mouse.wheel(0, -200)
    assert _view(page)["k"] > 1
    page.click("#reset-view")

    # A drag pans, and does not select what it started on.
    page.click("svg.contig-plot [data-contig='0'] > text")
    start = feature.bounding_box()
    page.mouse.move(start["x"] + 5, start["y"] + 5)
    page.mouse.down()
    page.mouse.move(start["x"] + 105, start["y"] + 45, steps=5)
    page.mouse.up()
    view = _view(page)
    assert view["x"] > 0 and view["y"] > 0
    assert page.text_content("#details .details-heading") == "Contig"
    # Nor does a drag that starts on a feature's label select the labels it crosses.
    label = page.locator("svg.contig-plot [data-contig='0'] [data-feature] text").first.bounding_box()
    page.mouse.move(label["x"] + 2, label["y"] + label["height"] / 2)
    page.mouse.down()
    page.mouse.move(label["x"] + 200, label["y"] + 120, steps=8)
    page.mouse.up()
    assert _selected_text(page) == ""
    # Double-clicking a label zooms in without selecting its word.
    label = page.locator("svg.contig-plot [data-contig='0'] [data-feature] text").first.bounding_box()
    page.mouse.dblclick(label["x"] + 4, label["y"] + label["height"] / 2)
    assert _selected_text(page) == ""
    page.click("#reset-view")

    with page.expect_download() as download:
        page.click("#download-svg")
    assert download.value.suggested_filename == "contigs.svg"
    assert open(download.value.path()).read().count("<path") == page.locator("svg.contig-plot path").count()
    assert page.pageerrors == []


def test_contig_plot_with_no_contigs(open_page, tmp_path, shared_datadir):
    page = open_page(_contig_page(tmp_path / "contigs.html", "-i", str(shared_datadir / "empty.gb")))
    page.wait_for_selector(".panel .status")
    assert page.text_content(".panel .status") == "No contigs to plot."
    assert page.pageerrors == [] and page.requests == []


def _projection_page(path, *extra):
    from domainator import build_projection
    build_projection.main(["--html_out", str(path), *extra])
    return path


def _projection(page, expression):
    return page.evaluate(f"window.projectionView.{expression}")


def test_projection_draws_every_point_and_hides_legend_groups(open_page, tmp_path, shared_datadir):
    path = _projection_page(tmp_path / "projection.html", "-i", str(shared_datadir / "scorefull.tsv"),
                            "--metadata", str(shared_datadir / "metadata_FeSOD_20.tsv"), "--color_by", "group")
    page = open_page(path)
    page.wait_for_function("window.projectionView !== undefined")
    assert page.pageerrors == [] and page.requests == []
    data = report_data(path.read_text())["data"]
    assert _projection(page, "drawnPoints()") == 20
    labels = page.eval_on_selector_all(".projection-plot [data-group] text", "t => t.map(x => x.textContent)")
    assert labels == [f"{group['value']} ({group['count']})" for group in data["groups"]]

    page.click(".projection-plot [data-group='1']")
    page.wait_for_function(f"window.projectionView.drawnPoints() === {20 - data['groups'][1]['count']}")
    assert page.get_attribute(".projection-plot [data-group='1']", "aria-pressed") == "false"
    # A hidden group is left out of the SVG download, points and legend.
    with page.expect_download() as download:
        page.click("#download-svg")
    assert download.value.suggested_filename == "projection.svg"
    svg = open(download.value.path()).read()
    assert svg.count("<circle") == 20 - data["groups"][1]["count"]
    assert data["groups"][1]["color"] not in svg
    page.focus(".projection-plot [data-group='1']")
    page.keyboard.press("Enter")
    page.wait_for_function("window.projectionView.drawnPoints() === 20")

    with page.expect_download() as download:
        page.click("#download-png")
    assert download.value.suggested_filename == "projection.png"
    assert open(download.value.path(), "rb").read(8) == b"\x89PNG\r\n\x1a\n"
    assert page.pageerrors == []


def test_projection_hover_zoom_pan_and_find(open_page, tmp_path, shared_datadir):
    page = open_page(_projection_page(tmp_path / "projection.html", "-i", str(shared_datadir / "scorefull.tsv"),
                                      "--metadata", str(shared_datadir / "FeSOD_metadata.tsv")))
    page.wait_for_function("window.projectionView !== undefined")
    name = "FeSOD_A0A067LT26|unreviewed|Superoxide"
    point = _projection(page, f"pointPosition({json.dumps(name)})")
    page.mouse.move(point["x"], point["y"])
    tip = page.text_content(".chart-tooltip").split("\n")
    assert tip[0] == name and "Accession: A0A067LT26" in tip and "expression level: 0.828329655" in tip

    fitted = _projection(page, "domain()")
    # Equal units on both axes.
    assert fitted["x1"] - fitted["x0"] == pytest.approx(fitted["y1"] - fitted["y0"])
    page.mouse.wheel(0, -200)
    page.wait_for_function(f"window.projectionView.domain().x1 - window.projectionView.domain().x0 < {fitted['x1'] - fitted['x0']}")
    # The point under the pointer stays there.
    moved = _projection(page, f"pointPosition({json.dumps(name)})")
    assert moved["x"] == pytest.approx(point["x"], abs=1) and moved["y"] == pytest.approx(point["y"], abs=1)
    page.click("#reset-view")
    page.wait_for_function(f"window.projectionView.domain().x0 === {fitted['x0']}")

    page.mouse.move(point["x"], point["y"])
    page.mouse.down()
    page.mouse.move(point["x"] + 100, point["y"] + 50, steps=5)
    page.mouse.up()
    page.wait_for_function(f"window.projectionView.domain().x0 < {fitted['x0']}")
    assert _projection(page, "domain()")["y0"] > fitted["y0"]

    page.fill("#point-search", "a0a067lt26")
    page.click("#find-point")
    assert page.text_content("#search-status") == name
    page.fill("#point-search", "no such point")
    page.keyboard.press("Enter")
    assert page.text_content("#search-status") == "No point matches."
    assert page.pageerrors == [] and page.requests == []


def test_projection_3d_switches_axis_pairs(open_page, tmp_path, shared_datadir):
    page = open_page(_projection_page(tmp_path / "projection.html", "-i", str(shared_datadir / "scorefull.tsv"), "--three_d"))
    page.wait_for_function("window.projectionView !== undefined")
    assert _projection(page, "axes()") == ["projection_X", "projection_Y"]
    page.click("#axes-yz")
    assert _projection(page, "axes()") == ["projection_Y", "projection_Z"]
    assert page.get_attribute("#axes-yz", "aria-pressed") == "true"
    assert page.get_attribute("#axes-xy", "aria-pressed") == "false"
    page.wait_for_function("[...document.querySelectorAll('.projection-plot svg text')].some(t => t.textContent === 'projection_Z')")
    with page.expect_download() as download:
        page.click("#download-svg")
    assert download.value.suggested_filename == "projection_YZ.svg"
    assert page.pageerrors == [] and page.requests == []
