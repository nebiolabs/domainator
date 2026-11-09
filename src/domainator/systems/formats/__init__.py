"""Converters from foreign system model formats to the native format, used by system_db_download.

A converter reads a local copy of a foreign database and returns a Conversion: the models as a
ModelCollection, the profile files they need, and a note for everything it approximated or
dropped. write_conversion writes it as a native model directory:

    <output_dir>/
      systems.yaml     the models, with database:, semantics:, source: and reference_files:
      <database>.hmm   every profile the models use (and <database>.cm for covariance models)
      LICENSE          copied from the source
      conversion.tsv   one row per approximation: model, feature, what was done

Profile files are named after the collection's database, because a hit's database qualifier
comes from its reference file's name (seqio.db_name_from_path). Then annotating with
``domainate -r <output_dir>/<database>.hmm`` and with ``domainate --system_models <output_dir>``
label hits the same way, so both find the same systems.
"""
from dataclasses import dataclass, field
import gzip
import os
from pathlib import Path
import shutil
import tarfile
import tempfile
from typing import Callable, Dict, List, NamedTuple, Optional, Sequence, Tuple

from domainator import __version__
from domainator.errors import InputError

from domainator.systems.model import ModelCollection
from domainator.systems.native import DIRECTORY_FILE, write_system_models

CONVERSION_REPORT = "conversion.tsv"


class Note(NamedTuple):
    """One thing a converter approximated or dropped."""
    model: str
    feature: str
    action: str


@dataclass
class Conversion:
    collection: ModelCollection
    profile_files: Dict[str, List[Path]] = field(default_factory=dict)
    """{output file name: the source files to concatenate into it}. The names are the
    collection's reference_files."""
    license_file: Optional[Path] = None
    notes: List[Note] = field(default_factory=list)


def write_conversion(conversion: Conversion, output_dir: os.PathLike) -> None:
    """Write ``conversion`` as a native model directory (see the module docstring)."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, sources in conversion.profile_files.items():
        _concatenate(sources, output_dir / name)
    if conversion.license_file is not None:
        shutil.copyfile(conversion.license_file, output_dir / "LICENSE")
    with open(output_dir / CONVERSION_REPORT, "w") as handle:
        print("model\tfeature\taction", file=handle)
        for note in conversion.notes:
            print("\t".join(value.replace("\t", " ") for value in note), file=handle)
    # Last, so a directory with a systems.yaml is complete.
    write_system_models(conversion.collection, output_dir / DIRECTORY_FILE)


def _concatenate(sources: Sequence[Path], destination: Path) -> None:
    """Concatenate text profile files (HMMER3 or Infernal, or gzipped), each ending in a newline."""
    with open(destination, "wb") as out:
        for source in sources:
            with (gzip.open if str(source).endswith(".gz") else open)(source, "rb") as handle:
                data = handle.read()
            out.write(data)
            if data and not data.endswith(b"\n"):
                out.write(b"\n")


def summarize(notes: Sequence[Note]) -> str:
    """A short summary of the notes, for the terminal: how often each kind of change was made."""
    counts: Dict[str, int] = dict()
    for note in notes:
        kind = note.action.split(":", 1)[0]
        counts[kind] = counts.get(kind, 0) + 1
    return "\n".join(f"  {count} x {kind}" for kind, count in sorted(counts.items(), key=lambda item: -item[1]))


def download_github_release(repository: str, release: Optional[str], directory: os.PathLike, *, label: str,
                            asset: Optional[Callable[[str], str]] = None) -> Tuple[Path, str, str]:
    """Download a release of a GitHub repository (default: the latest) into directory, and
    extract it.

    asset: a function of the release tag giving the name of the release's tarball asset to
    download. None: download the source archive of the tag.

    Returns:
        (the single top-level directory of the tarball, the release tag, the URL it came from)
    """
    import requests

    session = requests.Session()
    session.headers["User-Agent"] = f"domainator/{__version__}"
    if release is None:
        response = session.get(f"https://api.github.com/repos/{repository}/releases/latest", timeout=60)
        response.raise_for_status()
        release = response.json()["tag_name"]
    if asset is None:
        url = f"https://github.com/{repository}/archive/refs/tags/{release}.tar.gz"
    else:
        url = f"https://github.com/{repository}/releases/download/{release}/{asset(release)}"
    archive = Path(directory) / f"{label}-{release}.tar.gz"
    with session.get(url, stream=True, timeout=60) as response:
        if response.status_code == 404:
            raise InputError(f"{label} has no release {release!r} (see https://github.com/{repository}/releases).")
        response.raise_for_status()
        with open(archive, "wb") as handle:
            for chunk in response.iter_content(chunk_size=1 << 20):
                handle.write(chunk)
    extract_to = Path(tempfile.mkdtemp(dir=directory))
    with tarfile.open(archive) as tar:
        tar.extractall(extract_to, filter="data")
    roots = [child for child in extract_to.iterdir() if child.is_dir()]
    if len(roots) != 1:
        raise InputError(f"Unexpected layout in {url}: expected one top-level directory.")
    return roots[0], release, url
