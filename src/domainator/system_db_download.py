"""Download system models, or convert a local copy, into a directory that other commands read

Converts a foreign system model database into Domainator's own system model format, and writes
it to --output_dir:

    systems.yaml     the models (see docs/file_formats.md)
    padloc.hmm       every profile the models use, named after the collection's database
    LICENSE          the source's license
    conversion.tsv   every approximation the conversion made: model, feature, what was done

The directory can then be given to --system_models (find_systems, domainate, select_by_system
...). Its profile file can be searched with domainate -r, and the hits are labeled the same way
as with domainate --system_models:

    domainator system_db_download --source padloc --output_dir padloc_db
    domainator domainate -i genomes.gb -r padloc_db/padloc.hmm -e 1e-3 -o annotated.gb
    domainator find_systems -i annotated.gb --system_models padloc_db -o systems.gb

Sources:
    padloc: PADLOC's padloc-db (https://github.com/padlocbio/padloc-db, MIT license). Without
        --input, a release is downloaded from GitHub (--release, default: the latest). With
        --input, a local copy (a directory with sys/, hmm/, hmm_meta.txt and sys_meta.txt) is
        converted. CRISPR arrays can't be detected, so they are dropped from the models.
"""
import os
from pathlib import Path
import sys
import tempfile
from typing import List, Optional

from jsonargparse import ActionConfigFile

from domainator.cli.common import tool_parser
from domainator.errors import InputError
from domainator.systems.formats import Note, summarize, write_conversion
from domainator.systems.formats import padloc

SOURCES = ("padloc",)


def system_db_download(input: Optional[str] = None, *, source: str, output_dir: str,
                       release: Optional[str] = None) -> List[Note]:
    """Convert a foreign system model database into a native model directory.

    Args:
        input: a local copy of the database. None: download it.
        source: the database's format: "padloc".
        output_dir: the directory to write (see the module docstring).
        release: the release to download (default: the latest). Only without input.

    Returns:
        The notes on what the conversion approximated or dropped, also written to
        conversion.tsv in output_dir.
    """
    if source not in SOURCES:
        raise ValueError(f"source must be one of {', '.join(SOURCES)}, not {source!r}")
    if input is not None and release is not None:
        raise InputError("--release selects a release to download, so it can't be combined with --input.")
    if input is not None:
        conversion = padloc.convert(input, release=release)
        write_conversion(conversion, output_dir)
        return conversion.notes
    with tempfile.TemporaryDirectory() as directory:
        path, release, url = padloc.download(release, directory)
        conversion = padloc.convert(path, release=release, url=url)
        write_conversion(conversion, output_dir)
    return conversion.notes


def build_parser():
    parser = tool_parser(__doc__)
    parser.add_argument("--source", type=str, required=True, choices=SOURCES,
                        help="The format of the database to download or convert.")
    parser.add_argument("-i", "--input", type=str, default=None,
                        help="A local copy of the database to convert, instead of downloading it.")
    parser.add_argument("--release", type=str, default=None,
                        help="The release to download, e.g. v2.0.0. Default: the latest, which is recorded in "
                             "systems.yaml. Not with --input.")
    parser.add_argument("--output_dir", type=str, required=True,
                        help="The directory to write the models, profiles, license and conversion report to.")
    parser.add_argument("--config", action=ActionConfigFile)
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)
    if params.input is not None and params.release is not None:
        parser.error("--release selects a release to download, so it can't be combined with --input.")
    if params.input is not None and not os.path.isdir(params.input):
        parser.error(f"--input must be a directory: {params.input}")
    notes = system_db_download(params.input, source=params.source, output_dir=params.output_dir,
                               release=params.release)
    report = Path(params.output_dir) / "conversion.tsv"
    if notes:
        print(f"The conversion made {len(notes)} changes, listed in {report}:\n{summarize(notes)}", file=sys.stderr)
