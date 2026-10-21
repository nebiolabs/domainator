"""The standalone Domainator Similarity Network Viewer: one self-contained HTML page.

The page's sources are package files in ``domainator/ssn_viewer/``: ``viewer.html`` (the
shell), ``viewer.css``, and the scripts in ``js/``. ``js/viewer.js`` is the page's one script;
each ``// @include <file>`` line in it is replaced by that file from ``js/``, so the included
files share its scope. The values the page is built with (versions, palettes, an embedded
bundle) go in one ``<script type="application/json" id="dsnv-config">`` block, which viewer.js
reads into ``DSNV_CONFIG``.
"""

import base64
import json
import re
from importlib.resources import files
from pathlib import Path

from domainator import __version__
from domainator.output_guardrails import make_temporary_output_path
from domainator.ssn_bundle import (
    SSN_VIEWER_BUNDLE_FORMAT,
    SSN_VIEWER_BUNDLE_VERSION,
    SUPPORTED_SSN_VIEWER_BUNDLE_VERSIONS,
)
from domainator.ssn_hierarchy import (
    MERGE_IMPACT_CHOICES,
    MERGE_IMPACT_MIN_CHILD,
    merge_impact_axis_labels,
)
from domainator.colors import NAMED_CATEGORICAL_PALETTES, OTHER_COLOR


VIEWER_APP_NAME = "Domainator Similarity Network Viewer"
# Titles that name no particular network, so the heading shows the app name alone
# rather than "<app name>: <app name>". The old app name is kept here so bundles
# written before the rename still get the bare heading.
_GENERIC_TITLES = frozenset({"", VIEWER_APP_NAME, "Domainator SSN Viewer"})


def viewer_heading(network_name: str | None) -> str:
    """The page's <h1>: the app name, plus the network name when there is one."""
    name = (network_name or "").strip()
    return VIEWER_APP_NAME if name in _GENERIC_TITLES else f"{VIEWER_APP_NAME}: {name}"


_ASSETS = files("domainator") / "ssn_viewer"
_INCLUDE = re.compile(r"^// @include (\S+)$", re.MULTILINE)
_PLACEHOLDER = re.compile(r"\{\{(\w+)\}\}")


def _asset(name: str) -> str:
    return (_ASSETS / name).read_text(encoding="utf-8")


def _layout_core_js() -> str:
    """The layout math, shared by the layout Web Worker and the page (see js/layout_core.js)."""
    return _asset("js/layout_core.js")


def _layout_worker_js() -> str:
    """The layout Web Worker's source: the shared core plus its entry point."""
    return _layout_core_js() + _asset("js/layout_worker.js")


def _viewer_script() -> str:
    """js/viewer.js with its ``// @include`` lines replaced by the included files."""
    return _INCLUDE.sub(lambda match: _asset("js/" + match.group(1)), _asset("js/viewer.js"))


def _escape_html(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def ssn_viewer_html(
    title: str = VIEWER_APP_NAME,
    embedded_bundle_json: bytes | None = None,
) -> str:
    """The viewer page, optionally with a bundle embedded so that it opens with it."""
    config = {
        "embeddedBundleBase64": None if embedded_bundle_json is None else base64.b64encode(embedded_bundle_json).decode("ascii"),
        "layoutWorkerCode": _layout_worker_js(),
        # Kept in lockstep with ssn_bundle.py so the browser and the Python readers
        # agree on which bundle revisions they understand.
        "bundleFormat": SSN_VIEWER_BUNDLE_FORMAT,
        "bundleVersion": SSN_VIEWER_BUNDLE_VERSION,
        "supportedBundleVersions": list(SUPPORTED_SSN_VIEWER_BUNDLE_VERSIONS),
        "domainatorVersion": __version__,
        "viewerAppName": VIEWER_APP_NAME,
        # Shared with get_palette (domainator build_ssn and friends) so "Domainator distinct"
        # in the viewer means the same colors those tools write.
        "categoricalPalettes": [
            {"name": palette["name"], "label": palette["label"], "note": palette["note"], "colors": list(palette["colors"])}
            for palette in NAMED_CATEGORICAL_PALETTES
        ],
        "otherColor": OTHER_COLOR,
        # Generated from the same helper matrix_report uses, keyed by metric so a bundle built
        # with --merge_impact_metric product gets titles that match what its numbers mean.
        "splitAxisLabels": {
            metric: {
                "largest": merge_impact_axis_labels(metric)["largest"],
                "movingSum": merge_impact_axis_labels(metric)["moving_sum_short"],
                # One-placeholder templates for the hover readout, so it names an impact
                # the same way matrix_report's hovertemplates do.
                "impactAmount": merge_impact_axis_labels(metric)["impact_amount"],
                "movingSumHover": merge_impact_axis_labels(metric)["moving_sum_hover"],
            }
            for metric in MERGE_IMPACT_CHOICES
        },
        "defaultMergeImpactMetric": MERGE_IMPACT_MIN_CHILD,
    }
    values = {
        # The browser tab keeps the bare network name: it is the part that
        # distinguishes one open viewer from another once the tab is truncated.
        "title": _escape_html(title),
        "heading": _escape_html(viewer_heading(title)),
        "css": _asset("viewer.css").rstrip("\n"),
        # "<" escaped so that no value can close the <script> element.
        "config": json.dumps(config).replace("<", "\\u003c"),
        "script": _viewer_script().rstrip("\n"),
    }
    return _PLACEHOLDER.sub(lambda match: values[match.group(1)], _asset("viewer.html"))


def write_ssn_viewer_html(
    out_path: str,
    title: str = VIEWER_APP_NAME,
    embedded_bundle_json: bytes | None = None,
) -> None:
    temp_path = make_temporary_output_path(out_path)
    try:
        with open(temp_path, "w", encoding="utf-8") as out_handle:
            out_handle.write(ssn_viewer_html(title=title, embedded_bundle_json=embedded_bundle_json))
        Path(temp_path).replace(out_path)
        temp_path = None
    finally:
        if temp_path is not None and Path(temp_path).exists():
            Path(temp_path).unlink()