"""The HTML report pages (report_html.py): how data goes into them, and that their scripts parse."""
import json
import re
import shutil
import subprocess
from io import StringIO

import pytest

from domainator import report_html, table_writers
from helpers import report_data, report_rows

REPORTS = ["table_report", "summary_report"]
# Pages that stream rows into a report-rows block.
TABLE_PAGES = {"table_report"}


def test_table_writer_streams_the_ndjson_rows():
    html, ndjson = StringIO(), StringIO()
    columns, types = ["name", "length", "score"], ["str", "int", "float"]
    for writer in (table_writers.HTMLTableWriter(columns, types, html, 300, title="HMMER Report", report="hmmer_report"),
                   table_writers.JSONTableWriter(columns, types, ndjson)):
        writer.write_header()
        writer.write_row(["a", 3, 1.25])
        writer.write_row([None, None, None])
        writer.write_footer()

    page = html.getvalue()
    assert "<title>HMMER Report</title>" in page
    envelope = report_data(page)
    assert envelope["report"] == "hmmer_report"
    assert envelope["payload_version"] == report_html.PAYLOAD_VERSION
    assert envelope["data"] == {
        "columns": [{"name": "name", "type": "str"}, {"name": "length", "type": "int"}, {"name": "score", "type": "float"}],
        "max_height": 300,
        "file_stem": "hmmer_report",
    }
    # The same rows as --json_out, full precision and null for missing values.
    import json
    assert report_rows(page) == [json.loads(line) for line in ndjson.getvalue().splitlines()]
    assert report_rows(page) == [{"name": "a", "length": 3, "score": 1.25}, {"name": None, "length": None, "score": None}]


def test_no_value_can_close_a_script_element():
    html = StringIO()
    writer = table_writers.HTMLTableWriter(["name"], ["str"], html, None, title="</title><script>x</script>")
    writer.write_header()
    writer.write_row(["</script><script>alert(1)</script>"])
    writer.write_row(["{{rows}} <!-- x"])
    writer.write_footer()
    page = html.getvalue()
    assert page.count("</script>") == 3  # report-data, report-rows and the page script
    assert "<title>&lt;/title&gt;&lt;script&gt;x&lt;/script&gt;</title>" in page
    assert [row["name"] for row in report_rows(page)] == ["</script><script>alert(1)</script>", "{{rows}} <!-- x"]


def test_compressed_data_round_trips():
    data = {"values": list(range(100)), "label": "<b>"}
    envelope = report_html.report_envelope("matrix_report", data, encoding="gzip+base64")
    page = f'<script type="application/json" id="report-data">{envelope}</script>'
    assert report_data(page)["data"] == data
    with pytest.raises(ValueError):
        report_html.report_envelope("matrix_report", data, encoding="zip")


@pytest.mark.parametrize("name", REPORTS)
def test_shell_has_each_token_where_the_page_has_its_value(name):
    shell = report_html.report_shell(name)
    envelope = report_html.report_envelope(name, {"columns": []})
    filled = (shell.replace(report_html.TITLE_TOKEN, "T &amp; &lt;T&gt;")
                   .replace(report_html.DATA_TOKEN, envelope))
    if name in TABLE_PAGES:
        head, tail = report_html.report_parts(name, "T & <T>", name, {"columns": []})
        assert filled == head + report_html.ROWS_TOKEN + tail
        assert shell.count(report_html.ROWS_TOKEN) == 1
    else:
        assert filled == report_html.render_report(name, "T & <T>", name, {"columns": []})
        assert report_html.ROWS_TOKEN not in shell
        with pytest.raises(ValueError):
            report_html.report_parts(name, "T", name, {})
    assert shell.count(report_html.TITLE_TOKEN) == 2
    assert shell.count(report_html.DATA_TOKEN) == 1
    assert "{{" not in shell


@pytest.mark.parametrize("name", REPORTS)
def test_pages_use_no_external_resources(name):
    shell = report_html.report_shell(name)
    assert not re.search(r'<(script|link)[^>]+(src|href)=', shell)
    # The SVG namespace is a name, never fetched.
    shell = shell.replace('"http://www.w3.org/2000/svg"', "")
    assert "http://" not in shell and "https://" not in shell


@pytest.mark.skipif(shutil.which("node") is None, reason="node is not installed")
@pytest.mark.parametrize("name", REPORTS)
def test_report_script_parses(tmp_path, name):
    (tmp_path / "page.js").write_text(report_html.report_script(name))
    subprocess.run(["node", "--check", "page.js"], capture_output=True, text=True, check=True, cwd=tmp_path)


@pytest.mark.skipif(shutil.which("node") is None, reason="node is not installed")
def test_float_cells_match_the_tsv_formatting(tmp_path):
    values = [1.25, 1.35, -1.25, 2.75, 0.05, 0.15, -0.04, -0.0, 0.0, 12345.65, 1e22, 3.999, -2.5, 1e-300, 7.0]
    (tmp_path / "table.js").write_text(report_html._asset("common/table.js"))
    (tmp_path / "run.js").write_text(
        'const fs = require("fs"); eval(fs.readFileSync("table.js", "utf8") + ";globalThis.f = formatFixed1;");'
        f'process.stdout.write(JSON.stringify({json.dumps(values)}.map(v => f(v))));'
    )
    shown = json.loads(subprocess.run(["node", "run.js"], capture_output=True, text=True, check=True, cwd=tmp_path).stdout)
    assert shown == ["%0.1f" % v for v in values]
