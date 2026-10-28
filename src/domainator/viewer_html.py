"""The Domainator Viewer: one self-contained HTML page that opens the network, tree and
projection bundles (.dvb, viewer_bundle.py) that build_ssn, build_tree and build_projection
write, and the network bundles (.dsnv) of earlier versions.

The page's sources are package files in ``domainator/viewer_app/``: ``viewer.html`` (the shell),
``viewer.css``, and the scripts in ``js/``. ``js/viewer.js`` is the page's one script; each
``// @include <path>`` line in a script is replaced by the file at ``<path>`` (relative to
``js/``), recursively, so every included file shares one scope. Code every view uses is in
``js/common/``, and each view's own code in ``js/views/<view>/``. The values the page is built
with (versions, palettes, an embedded bundle) go in one
``<script type="application/json" id="viewer-config">`` block, which viewer.js reads into
``VIEWER_CONFIG``.
"""

import base64
import json
import re
from importlib.resources import files
from pathlib import Path

from domainator import __version__
from domainator.colors import NAMED_CATEGORICAL_PALETTES, OTHER_COLOR
from domainator.output_guardrails import enforce_output_limit, make_temporary_output_path
from domainator.ssn_bundle import SSN_VIEWER_BUNDLE_FORMAT, SUPPORTED_SSN_VIEWER_BUNDLE_VERSIONS
from domainator.ssn_hierarchy import (
    MERGE_IMPACT_CHOICES,
    MERGE_IMPACT_MIN_CHILD,
    merge_impact_axis_labels,
)
from domainator.viewer_bundle import (
    SUPPORTED_VIEWER_BUNDLE_VERSIONS,
    VIEWER_BUNDLE_FORMAT,
    VIEWER_BUNDLE_VERSION,
    serialize_bundle,
    write_bundle_bytes,
)


VIEWER_APP_NAME = "Domainator Viewer"
# Titles that name nothing in particular, so the heading shows the app name alone rather
# than "<app name>: <app name>". The older app names are kept here so bundles written
# under them still get the bare heading.
_GENERIC_TITLES = frozenset({"", VIEWER_APP_NAME, "Domainator Similarity Network Viewer", "Domainator SSN Viewer"})

_ASSETS = files("domainator") / "viewer_app"
_INCLUDE = re.compile(r"^// @include (\S+)$", re.MULTILINE)
_PLACEHOLDER = re.compile(r"\{\{(\w+)\}\}")


def viewer_heading(name: str | None) -> str:
    """The page's <h1>: the app name, plus the bundle's name when there is one."""
    name = (name or "").strip()
    return VIEWER_APP_NAME if name in _GENERIC_TITLES else f"{VIEWER_APP_NAME}: {name}"


def _asset(name: str) -> str:
    return (_ASSETS / name).read_text(encoding="utf-8")


def _script(name: str, _including: tuple = ()) -> str:
    """js/<name> with its ``// @include`` lines replaced by the included files."""
    if name in _including:
        raise ValueError(f"js/{name} includes itself: {' -> '.join(_including + (name,))}")
    return _INCLUDE.sub(lambda match: _script(match.group(1), _including + (name,)), _asset("js/" + name))


def _layout_core_js() -> str:
    """The SSN view's layout math, shared by its layout Web Worker and the page (see
    js/views/ssn/layout_core.js)."""
    return _asset("js/views/ssn/layout_core.js")


def _layout_worker_js() -> str:
    """The layout Web Worker's source: the shared core plus its entry point."""
    return _layout_core_js() + _asset("js/views/ssn/layout_worker.js")


def _escape_html(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def viewer_html(title: str = VIEWER_APP_NAME, embedded_bundle: bytes | None = None) -> str:
    """The viewer page, optionally with a bundle file's bytes (gzip-compressed JSON) embedded
    so that it opens with it."""
    config = {
        "embeddedBundleBase64": None if embedded_bundle is None else base64.b64encode(embedded_bundle).decode("ascii"),
        "layoutWorkerCode": _layout_worker_js(),
        # Kept in lockstep with ssn_bundle.py and viewer_bundle.py so the browser and the
        # Python readers agree on which bundle revisions they understand.
        "ssnBundleFormat": SSN_VIEWER_BUNDLE_FORMAT,
        "supportedSsnBundleVersions": list(SUPPORTED_SSN_VIEWER_BUNDLE_VERSIONS),
        "viewerBundleFormat": VIEWER_BUNDLE_FORMAT,
        "viewerBundleVersion": VIEWER_BUNDLE_VERSION,
        "supportedViewerBundleVersions": list(SUPPORTED_VIEWER_BUNDLE_VERSIONS),
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
        # The browser tab keeps the bare bundle name: it is the part that distinguishes one
        # open viewer from another once the tab is truncated.
        "title": _escape_html(title),
        "heading": _escape_html(viewer_heading(title)),
        "css": _asset("viewer.css").rstrip("\n"),
        # "<" escaped so that no value can close the <script> element.
        "config": json.dumps(config).replace("<", "\\u003c"),
        "script": _script("viewer.js").rstrip("\n"),
    }
    return _PLACEHOLDER.sub(lambda match: values[match.group(1)], _asset("viewer.html"))


def write_viewer_html(out_path, title: str = VIEWER_APP_NAME, embedded_bundle: bytes | None = None) -> None:
    """Write the viewer page to ``out_path`` (a path or a text handle). A path is written through
    a temporary file, so a failed write leaves no partial page."""
    if hasattr(out_path, "write"):
        out_path.write(viewer_html(title=title, embedded_bundle=embedded_bundle))
        return
    temp_path = make_temporary_output_path(out_path)
    try:
        with open(temp_path, "w", encoding="utf-8") as out_handle:
            out_handle.write(viewer_html(title=title, embedded_bundle=embedded_bundle))
        Path(temp_path).replace(out_path)
        temp_path = None
    finally:
        if temp_path is not None and Path(temp_path).exists():
            Path(temp_path).unlink()


def write_viewer_outputs(bundle: dict, *, html_out=None, dvb_out=None, max_output_bytes: int | None = None,
                         mitigation_options=()) -> None:
    """Write ``bundle`` (viewer_bundle.viewer_bundle) to a .dvb file at ``dvb_out`` and as a
    viewer page that opens it at ``html_out``, each optional, failing first if either would
    exceed ``max_output_bytes``."""
    compressed = serialize_bundle(bundle)
    if dvb_out is not None:
        write_bundle_bytes(dvb_out, compressed, description="Domainator Viewer bundle output",
                           max_output_bytes=max_output_bytes, mitigation_options=mitigation_options)
    if html_out is not None:
        # The page embeds the bundle's bytes in base64.
        enforce_output_limit(projected_bytes=len(compressed) * 4 // 3, max_output_bytes=max_output_bytes,
                             output_description=f"Domainator Viewer page '{html_out}'",
                             mitigation_options=list(mitigation_options))
        write_viewer_html(html_out, title=bundle["name"] or VIEWER_APP_NAME, embedded_bundle=compressed)
