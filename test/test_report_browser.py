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
    # One bar per bin, and the bars hold every contig.
    data = page.evaluate("JSON.parse(document.getElementById('report-data').textContent).data")
    bars = page.eval_on_selector_all("#length-histogram svg rect[data-tip]", "bars => bars.map(b => b.dataset.tip)")
    assert len(bars) == len(data["length_histogram"]["counts"])
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
