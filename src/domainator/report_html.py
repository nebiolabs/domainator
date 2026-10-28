"""The HTML reports' pages: self-contained files with no third-party JavaScript.

The page sources are package files in ``domainator/web/``: ``reports/<name>.html`` (the page),
``reports/<name>.js`` (its one script) and ``common/`` (the stylesheet and the scripts the
reports share). Each ``// @include <path>`` line in a report's script is replaced by that file
from ``web/``, so the included files share its scope.

A report's data goes in one ``<script type="application/json" id="report-data">`` block, an
envelope ``{"report", "payload_version", "encoding", "data"}`` that ``common/payload.js``
reads. A table report also streams its rows, one JSON object per line, into a
``<script type="application/x-ndjson" id="report-rows">`` block, so that it never holds a
whole table in memory.

The pages and their data are not part of the 1.x interface (see docs/stability.md); the
JSON outputs (``--json_out``) are. ``payload_version`` changes whenever the data a page reads
changes shape.
"""

import base64
import gzip
import json
import re
from importlib.resources import files

PAYLOAD_VERSION = 2

# The tokens report_shell() leaves in a page, for programs that fill pages without this module.
TITLE_TOKEN = "@@REPORT_TITLE@@"
DATA_TOKEN = "@@REPORT_DATA@@"
ROWS_TOKEN = "@@REPORT_ROWS@@"

_ASSETS = files("domainator") / "web"
_INCLUDE = re.compile(r"^// @include (\S+)$", re.MULTILINE)
_PLACEHOLDER = re.compile(r"\{\{(\w+)\}\}")
_ROWS = "{{rows}}"


def _asset(path: str) -> str:
    return (_ASSETS / path).read_text(encoding="utf-8")


def report_script(name: str) -> str:
    """reports/<name>.js with its ``// @include`` lines replaced by the included files."""
    return _INCLUDE.sub(lambda match: _asset(match.group(1)), _asset(f"reports/{name}.js"))


def escape_html(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def script_json(value) -> str:
    """``value`` as compact JSON that is safe inside a <script> element.

    Only a string can hold "<", and inside a string "\\u003c" means the same thing, so no value
    can close the element or open a comment.
    """
    return json.dumps(value, separators=(",", ":")).replace("<", "\\u003c")


def report_envelope(report: str, data, encoding: str = "json") -> str:
    """The report-data block's contents. ``encoding="gzip+base64"`` compresses ``data``."""
    if encoding == "gzip+base64":
        compressed = gzip.compress(json.dumps(data, separators=(",", ":")).encode("utf-8"), mtime=0)
        data = base64.b64encode(compressed).decode("ascii")
    elif encoding != "json":
        raise ValueError(f"Unknown report data encoding: {encoding!r}")
    return script_json({"report": report, "payload_version": PAYLOAD_VERSION, "encoding": encoding, "data": data})


def _render(name: str, title: str, data: str) -> list[str]:
    """reports/<name>.html filled in, as the text before each {{rows}} and after the last.

    The page is split before it is filled, so a value that happens to contain "{{rows}}" stays
    a value.
    """
    values = {
        "title": title,
        "css": _asset("common/report.css").rstrip("\n"),
        "data": data,
        "script": report_script(name).rstrip("\n"),
    }
    return [
        _PLACEHOLDER.sub(lambda match: values[match.group(1)], part)
        for part in _asset(f"reports/{name}.html").split(_ROWS)
    ]


def render_report(name: str, title: str, report: str, data, encoding: str = "json") -> str:
    """A page with no row block: reports/<name>.html filled with ``data``."""
    parts = _render(name, escape_html(title), report_envelope(report, data, encoding))
    if len(parts) != 1:
        raise ValueError(f"reports/{name}.html streams rows; use report_parts")
    return parts[0]


def report_parts(name: str, title: str, report: str, data, encoding: str = "json") -> tuple[str, str]:
    """A page that streams rows, as the text before the rows and the text after them.

    Each row written between the two is report_row(row).
    """
    head, tail = _render(name, escape_html(title), report_envelope(report, data, encoding))
    return head, tail


def report_row(row) -> str:
    """One line of a report-rows block."""
    return script_json(row) + "\n"


def report_shell(name: str) -> str:
    """The page with its title, data and rows left as tokens.

    For programs that write these reports without this module (domainatoredux vendors the
    shells): replace TITLE_TOKEN with the HTML-escaped title, DATA_TOKEN with
    report_envelope()'s text, and ROWS_TOKEN (table reports only) with the report_row() lines.
    """
    return ROWS_TOKEN.join(_render(name, TITLE_TOKEN, DATA_TOKEN))
