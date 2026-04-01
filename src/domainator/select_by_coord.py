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

If a file column is supplied, it must contain basenames only (with extension).
Rows without a file value match any input file.
"""

from collections import defaultdict
from dataclasses import dataclass
import csv
import os
import sys

from jsonargparse import ActionConfigFile, ArgumentParser

from domainator import __version__, RawAndDefaultsFormatter
from domainator.Bio.SeqFeature import CompoundLocation, FeatureLocation
from domainator.utils import parse_seqfiles, slice_record_from_location, write_genbank


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


def parse_seqfiles_with_basename(seqfiles, contigs=None, filetype_override=None):
    for file in seqfiles:
        if isinstance(file, (str, os.PathLike)):
            basename = os.path.basename(file)
        else:
            basename = os.path.basename(getattr(file, "name", "<stdin>"))

        for rec in parse_seqfiles([file], contigs=contigs, filetype_override=filetype_override):
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
    record.id = f"{record.id}_{selection.start}:{selection.end}"
    if selection.strand == -1:
        record.id += "rc"
    record.name = record.id
    return record


def select_by_coord(records, coordinate_selections):
    selections_by_contig = defaultdict(list)
    for selection in coordinate_selections:
        selections_by_contig[selection.contig].append(selection)

    for basename, contig in records:
        for selection in selections_by_contig.get(contig.id, []):
            if selection.file is not None and selection.file != basename:
                continue
            yield extract_coordinate_region(contig, selection)


def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument("-i", "--input", default=None, nargs="+", type=str, required=False,
                        help="names of input genbank files. If not supplied, reads from stdin.")
    parser.add_argument("-o", "--output", default=None, required=False,
                        help="genbank output file name. If not supplied, writes to stdout.")
    parser.add_argument("-c", "--coords", required=True, type=str,
                        help="labeled coordinate table with columns contig, start, end, strand (optional), file (optional).")
    parser.add_argument("--config", action=ActionConfigFile)

    params = parser.parse_args(argv)

    coordinate_selections = load_coordinate_table(params.coords)
    target_contigs = {selection.contig for selection in coordinate_selections}

    filetype_override = None
    if params.input is None:
        genbanks = [sys.stdin]
        filetype_override = "genbank"
    else:
        genbanks = params.input

    if params.output is None:
        output_handle = sys.stdout
    else:
        output_handle = open(params.output, "w")

    extracted_regions = select_by_coord(
        parse_seqfiles_with_basename(genbanks, contigs=target_contigs, filetype_override=filetype_override),
        coordinate_selections,
    )
    write_genbank(extracted_regions, output_handle)

    if params.output is not None:
        output_handle.close()


def _entrypoint():
    main(sys.argv[1:])


if __name__ == "__main__":
    main(sys.argv[1:])