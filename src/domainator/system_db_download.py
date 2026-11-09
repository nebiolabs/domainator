"""Download system models, or convert a local copy, into a directory that other commands read

Converts a foreign system model database into Domainator's own system model format, and writes
it to --output_dir:

    systems.yaml     the models (see docs/file_formats.md)
    <database>.hmm   every profile the models use, named after the collection's database
                     (padloc, defensefinder, or a MacSyFinder package's name)
    LICENSE          the source's license
    conversion.tsv   every approximation the conversion made: model, feature, what was done

The directory can then be given to --system_models (find_systems, domainate, system_search
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
    defensefinder: DefenseFinder's defense-finder-models (https://github.com/mdmparis/defense-finder-models;
        its metadata gives the licenses GPL-3.0 and CC BY-NC-SA 4.0). Without --input, a release
        is downloaded from GitHub. With --input, a local copy (a MacSyFinder package directory,
        with metadata.yml, definitions/ and profiles/, as installed by defense-finder update or
        macsydata) is converted. Hits must cover --coverage_profile of a profile (default 0.4,
        as DefenseFinder uses) and pass its gathering threshold.
    macsyfinder: any other MacSyFinder model package (models format 2.0), given with --input. It
        is named after its directory. --coverage_profile defaults to the package's
        model_conf.xml, else MacSyFinder's 0.5.

MacSyFinder loners and multi-locus models are found as single loci, and listed in
conversion.tsv.
"""
import os
from pathlib import Path
import sys
import tempfile
from typing import List, Optional

from jsonargparse import ActionConfigFile

from domainator.cli.common import tool_parser
from domainator.errors import InputError
from domainator.systems.formats import Conversion, Note, download_github_release, summarize, write_conversion
from domainator.systems.formats import macsyfinder, padloc

SOURCES = ("padloc", "defensefinder", "macsyfinder")
DOWNLOADS = {"padloc": (padloc.REPOSITORY, "padloc-db", None),
             "defensefinder": (macsyfinder.DEFENSEFINDER_REPOSITORY, "defense-finder-models", macsyfinder.defensefinder_asset)}


def _convert(source: str, path, *, release: Optional[str], url: Optional[str],
             coverage_profile: Optional[float]) -> Conversion:
    if source == "padloc":
        return padloc.convert(path, release=release, url=url)
    if source == "defensefinder":
        return macsyfinder.convert(path, database=macsyfinder.DEFENSEFINDER_DATABASE, release=release, url=url,
                                   coverage_profile=coverage_profile,
                                   default_coverage=macsyfinder.DEFENSEFINDER_COVERAGE_PROFILE,
                                   source_name="defense-finder-models")
    return macsyfinder.convert(path, database=macsyfinder.database_name(path), coverage_profile=coverage_profile)


def system_db_download(input: Optional[str] = None, *, source: str, output_dir: str,
                       release: Optional[str] = None, coverage_profile: Optional[float] = None) -> List[Note]:
    """Convert a foreign system model database into a native model directory.

    Args:
        input: a local copy of the database. None: download it (not for "macsyfinder").
        source: the database's format: "padloc", "defensefinder" or "macsyfinder".
        output_dir: the directory to write (see the module docstring).
        release: the release to download (default: the latest). Only without input.
        coverage_profile: for MacSyFinder packages, the fraction of a profile a hit must cover.
            None: DefenseFinder's 0.4, or a package's model_conf.xml, else 0.5.

    Returns:
        The notes on what the conversion approximated or dropped, also written to
        conversion.tsv in output_dir.
    """
    if source not in SOURCES:
        raise ValueError(f"source must be one of {', '.join(SOURCES)}, not {source!r}")
    if input is not None and release is not None:
        raise InputError("--release selects a release to download, so it can't be combined with --input.")
    if coverage_profile is not None and source == "padloc":
        raise InputError("--coverage_profile is for MacSyFinder packages: padloc-db gives each profile's coverage.")
    if input is None and source not in DOWNLOADS:
        raise InputError(f"--source {source} has nothing to download: give the package directory with --input.")
    if input is not None:
        conversion = _convert(source, input, release=None, url=None, coverage_profile=coverage_profile)
        write_conversion(conversion, output_dir)
        return conversion.notes
    repository, label, asset = DOWNLOADS[source]
    with tempfile.TemporaryDirectory() as directory:
        path, release, url = download_github_release(repository, release, directory, label=label, asset=asset)
        conversion = _convert(source, path, release=release, url=url, coverage_profile=coverage_profile)
        write_conversion(conversion, output_dir)
    return conversion.notes


def build_parser():
    parser = tool_parser(__doc__)
    parser.add_argument("--source", type=str, required=True, choices=SOURCES,
                        help="The format of the database to download or convert.")
    parser.add_argument("-i", "--input", type=str, default=None,
                        help="A local copy of the database to convert, instead of downloading it.")
    parser.add_argument("--release", type=str, default=None,
                        help="The release to download, e.g. v2.0.0 (padloc) or 3.1.0 (defensefinder). Default: the "
                             "latest, which is recorded in systems.yaml. Not with --input.")
    parser.add_argument("--coverage_profile", type=float, default=None,
                        help="For defensefinder and macsyfinder: the fraction of a profile a hit must cover. Default: "
                             "0.4 for defensefinder (as DefenseFinder runs MacSyFinder); for macsyfinder, the "
                             "package's model_conf.xml, else 0.5.")
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
    if params.input is None and params.source not in DOWNLOADS:
        parser.error(f"--source {params.source} needs --input: the package directory to convert.")
    if params.coverage_profile is not None and params.source == "padloc":
        parser.error("--coverage_profile is for MacSyFinder packages (defensefinder, macsyfinder).")
    if params.coverage_profile is not None and not 0 <= params.coverage_profile <= 1:
        parser.error("--coverage_profile must be between 0 and 1.")
    notes = system_db_download(params.input, source=params.source, output_dir=params.output_dir,
                               release=params.release, coverage_profile=params.coverage_profile)
    report = Path(params.output_dir) / "conversion.tsv"
    if notes:
        print(f"The conversion made {len(notes)} changes, listed in {report}:\n{summarize(notes)}", file=sys.stderr)
