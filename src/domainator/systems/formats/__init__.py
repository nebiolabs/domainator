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
import os
from pathlib import Path
import shutil
from typing import Dict, List, NamedTuple, Optional, Sequence

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
    """Concatenate text profile files (HMMER3 or Infernal), each ending in a newline."""
    with open(destination, "wb") as out:
        for source in sources:
            with open(source, "rb") as handle:
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
