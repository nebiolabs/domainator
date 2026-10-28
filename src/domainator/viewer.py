"""Write the Domainator Viewer, empty or opening a viewer bundle.

The viewer is one self-contained HTML page that opens the network, tree and projection
bundles (.dvb) that build_ssn, build_tree and build_projection write, and the sessions it
saves itself. Without -i/--input the page opens bundles you choose in it; with one, it opens
that bundle on load.
"""

import gzip
import json
from os import PathLike
from pathlib import Path
from typing import Optional, Union

from jsonargparse import ActionConfigFile

from domainator.cli.common import add_input, add_output, check_outputs, declare_outputs, tool_parser
from domainator.errors import FormatVersionError, InputError
from domainator.ssn_bundle import SSN_VIEWER_BUNDLE_FORMAT, SUPPORTED_SSN_VIEWER_BUNDLE_VERSIONS
from domainator.viewer_bundle import SUPPORTED_VIEWER_BUNDLE_VERSIONS, VIEWER_BUNDLE_FORMAT
from domainator.viewer_html import VIEWER_APP_NAME, write_viewer_html


def _read_bundle(path: Union[str, PathLike]) -> tuple:
    """The bundle file at ``path`` as (gzip-compressed bytes, parsed JSON)."""
    raw = Path(path).read_bytes()
    compressed = raw if raw[:2] == b"\x1f\x8b" else gzip.compress(raw, mtime=0)
    try:
        bundle = json.loads(gzip.decompress(compressed))
    except (OSError, EOFError, ValueError) as exc:
        raise InputError(f"'{path}' is not a Domainator Viewer bundle: {exc}") from exc
    fmt = bundle.get("format") if isinstance(bundle, dict) else None
    supported = {VIEWER_BUNDLE_FORMAT: SUPPORTED_VIEWER_BUNDLE_VERSIONS,
                 SSN_VIEWER_BUNDLE_FORMAT: SUPPORTED_SSN_VIEWER_BUNDLE_VERSIONS}
    if fmt not in supported:
        raise InputError(f"'{path}' is not a Domainator Viewer bundle (format='{fmt}', expected '{VIEWER_BUNDLE_FORMAT}').")
    if bundle.get("version") not in supported[fmt]:
        raise FormatVersionError(f"Unsupported {fmt} version {bundle.get('version')}; this build understands "
                                 f"version(s) {', '.join(str(v) for v in supported[fmt])}.")
    return compressed, bundle


def viewer(input: Optional[Union[str, PathLike]] = None, *, html_out: Union[str, PathLike],
           name: Optional[str] = None) -> None:
    """Write the Domainator Viewer page.

    Args:
        input: a viewer bundle (.dvb, or an older network bundle, .dsnv) for the page to open
            on load. None: the page opens bundles chosen in it.
        html_out: write the page to this HTML file.
        name: the page's title. Defaults to the bundle's name.
    """
    if input is None:
        write_viewer_html(html_out, title=name or VIEWER_APP_NAME)
        return
    compressed, bundle = _read_bundle(input)
    write_viewer_html(html_out, title=name or bundle.get("name") or VIEWER_APP_NAME, embedded_bundle=compressed)


def build_parser():
    parser = tool_parser(__doc__)
    add_input(parser, "A Domainator Viewer bundle (.dvb) for the page to open on load: from build_ssn, build_tree or build_projection --dvb_out, or the viewer's Save session. Older network bundles (.dsnv) work too. Without it, the page opens the bundles you choose in it.",
              multiple=False, required=False, stdin=False)
    add_output(parser, "html", "Write the viewer page to this HTML file.", primary=True, stdout=False)
    parser.add_argument("--name", type=str, default=None,
                        help="The page's title. Defaults to the bundle's name.")
    parser.add_argument("--config", action=ActionConfigFile)
    declare_outputs(parser, ["html_out"])
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)
    check_outputs(parser, params)
    viewer(params.input, html_out=params.html_out, name=params.name)
