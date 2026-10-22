"""Extract contig regions using explicit coordinate ranges.

Reads a labeled coordinate table into memory, then iterates through the input
contigs and extracts one region for each matching row.

Required table columns:
    contig, start, end

Optional table columns:
    strand, file

Coordinates are 1-based and inclusive.

If end < start, the region is assumed to wrap around the origin of a circular
contig. The same coordinates on a linear contig raise an exception.

Extracted regions are named like those of the other Domainator tools
(see utils.region_id): contig_first:last, where first and last are the
coordinates of the first and last bases of the region as written, followed by
'rc' for reverse-strand regions. For example, start 10, end 20, strand - is
named contig_20:10rc. Whole contigs keep their names (contig_rc if reverse).

If a file column is supplied, it must contain basenames only (with extension).
Rows without a file value match any input file.
"""

from typing import Iterable, Iterator
from collections import defaultdict
from dataclasses import dataclass
import csv
import os

from jsonargparse import ActionConfigFile

from domainator.cli.common import add_input, add_seq_outputs, check_outputs, declare_outputs, tool_parser
from domainator.outputs import write_seq_outputs
from Bio.SeqFeature import CompoundLocation, FeatureLocation
from domainator.seqio import parse_seqfiles
from domainator.model import region_id
from domainator.locations import slice_record_from_location
from Bio.SeqRecord import SeqRecord


ALLOWED_COLUMNS = {"contig", "start", "end", "strand", "file"}
REQUIRED_COLUMNS = {"contig", "start", "end"}


@dataclass(frozen=True)
class CoordinateSelection:
    contig: str
    start: int
    end: int
    strand: int = 1
    file: str | None = None


def _parse_strand(value, row_number):
    if value is None:
        return 1

    normalized = str(value).strip().lower()
    if normalized == "":
        return 1
    if normalized in {"+", "1", "f", "forward"}:
        return 1
    if normalized in {"-", "-1", "r", "reverse"}:
        return -1
    raise ValueError(f"Invalid strand value on row {row_number}: {value}")


def _parse_coordinate(value, column_name, row_number):
    try:
        parsed = int(str(value).strip())
    except ValueError as exc:
        raise ValueError(f"Column '{column_name}' must contain integers. Row {row_number}: {value}") from exc

    if parsed < 1:
        raise ValueError(f"Column '{column_name}' must be >= 1. Row {row_number}: {parsed}")
    return parsed


def load_coordinate_table(table_path):
    with open(table_path, newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        if sample.strip() == "":
            raise ValueError("Coordinate table is empty.")

        try:
            dialect = csv.Sniffer().sniff(sample, delimiters="\t,;")
        except csv.Error:
            dialect = csv.excel_tab

        reader = csv.DictReader(handle, dialect=dialect)
        if reader.fieldnames is None:
            raise ValueError("Coordinate table must include a header row.")

        fieldnames = [name.lstrip("\ufeff").strip() for name in reader.fieldnames]
        if len(fieldnames) != len(set(fieldnames)):
            raise ValueError("Coordinate table contains duplicate column names.")

        unexpected_columns = set(fieldnames) - ALLOWED_COLUMNS
        if unexpected_columns:
            unexpected = ", ".join(sorted(unexpected_columns))
            raise ValueError(f"Unexpected coordinate table columns: {unexpected}")

        missing_columns = REQUIRED_COLUMNS - set(fieldnames)
        if missing_columns:
            missing = ", ".join(sorted(missing_columns))
            raise ValueError(f"Missing required coordinate table columns: {missing}")

        reader.fieldnames = fieldnames
        selections = []
        for row_number, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(f"Row {row_number} has more fields than the header defines.")

            if all(value is None or str(value).strip() == "" for value in row.values()):
                continue

            contig = row["contig"].strip()
            if contig == "":
                raise ValueError(f"Column 'contig' must not be blank. Row {row_number}")

            file_name = row.get("file")
            if file_name is not None:
                file_name = file_name.strip() or None
            if file_name is not None and os.path.basename(file_name) != file_name:
                raise ValueError(f"Column 'file' must contain basenames only. Row {row_number}: {file_name}")

            selections.append(
                CoordinateSelection(
                    contig=contig,
                    start=_parse_coordinate(row["start"], "start", row_number),
                    end=_parse_coordinate(row["end"], "end", row_number),
                    strand=_parse_strand(row.get("strand"), row_number),
                    file=file_name,
                )
            )

    return selections


def parse_seqfiles_with_basename(seqfiles, contigs=None):
    for file in seqfiles:
        if isinstance(file, (str, os.PathLike)):
            basename = os.path.basename(file)
        else:
            basename = os.path.basename(getattr(file, "name", "<stdin>"))

        for rec in parse_seqfiles([file], contigs=contigs):
            yield basename, rec


def _build_location(contig, selection):
    contig_length = len(contig)
    if selection.start > contig_length or selection.end > contig_length:
        raise ValueError(
            f"Coordinates out of bounds for contig '{contig.id}': "
            f"{selection.start}-{selection.end} on length {contig_length}"
        )

    start = selection.start - 1
    end = selection.end

    if selection.end >= selection.start:
        return FeatureLocation(start, end, selection.strand)

    if contig.annotations.get("topology", "linear") != "circular":
        raise ValueError(
            f"end < start is only valid for circular contigs. "
            f"Contig '{contig.id}' is not circular."
        )

    if selection.strand == -1:
        parts = [
            FeatureLocation(0, end, selection.strand),
            FeatureLocation(start, contig_length, selection.strand),
        ]
    else:
        parts = [
            FeatureLocation(start, contig_length, selection.strand),
            FeatureLocation(0, end, selection.strand),
        ]
    return CompoundLocation(parts, "join")


def extract_coordinate_region(contig, selection):
    location = _build_location(contig, selection)
    record = slice_record_from_location(contig, location)
    record.id = region_id(contig, location)
    return record


def select_by_coord(records: Iterable[SeqRecord], *, coordinate_selections: Iterable[CoordinateSelection]) -> Iterator[SeqRecord]:
    """Yield the regions of ``records`` that ``coordinate_selections`` lists.

    coordinate_selections: CoordinateSelection rows, as load_coordinate_table reads them from a
        table of contig, start, end (1-based, inclusive) and optionally strand and file. Each row
        that names a record gives one region, named like the regions of the other tools
        (model.region_id); a reverse-strand region is reverse complemented.
    """
    selections_by_contig = defaultdict(list)
    for selection in coordinate_selections:
        selections_by_contig[selection.contig].append(selection)

    for basename, contig in records:
        for selection in selections_by_contig.get(contig.id, []):
            if selection.file is not None and selection.file != basename:
                continue
            yield extract_coordinate_region(contig, selection)


def build_parser():
    parser = tool_parser(__doc__)

    add_input(parser, 'Names of input genbank files.')
    add_seq_outputs(parser)
    parser.add_argument("--coords", required=True, type=str,
                        help="labeled coordinate table with columns contig, start, end, strand (optional), file (optional).")
    parser.add_argument("--config", action=ActionConfigFile)
    declare_outputs(parser, ["genbank_out", "fasta_out"])
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)

    coordinate_selections = load_coordinate_table(params.coords)
    target_contigs = {selection.contig for selection in coordinate_selections}

    genbanks = params.input

    check_outputs(parser, params)

    extracted_regions = select_by_coord(
        parse_seqfiles_with_basename(genbanks, contigs=target_contigs),
        coordinate_selections=coordinate_selections,
    )
    write_seq_outputs(extracted_regions, genbank=params.genbank_out, fasta=params.fasta_out)
