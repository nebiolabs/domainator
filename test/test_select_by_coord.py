import tempfile

from domainator import select_by_coord
from domainator.Bio import SeqIO
from domainator.utils import parse_seqfiles
import pytest


def _write_table(path, header, rows):
    with open(path, "w") as handle:
        print("\t".join(header), file=handle)
        for row in rows:
            print("\t".join(str(value) for value in row), file=handle)


def test_select_by_coord_linear_region_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/extraction.gb"
        coords = output_dir + "/coords.tsv"
        _write_table(coords, ["contig", "start", "end"], [("pDONR201_3", 2, 6)])

        select_by_coord.main(["-i", str(shared_datadir / "simple_genpept.gb"), "-o", out, "-c", coords])

        source = {record.id: record for record in SeqIO.parse(str(shared_datadir / "simple_genpept.gb"), "genbank")}
        records = list(parse_seqfiles([out]))
        assert len(records) == 1
        assert records[0].id == "pDONR201_3_2:6"
        assert str(records[0].seq) == str(source["pDONR201_3"].seq[1:6])


def test_select_by_coord_repeated_contig_rows_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/extraction.gb"
        coords = output_dir + "/coords.tsv"
        _write_table(coords, ["contig", "start", "end"], [("pDONR201_3", 1, 3), ("pDONR201_3", 4, 6)])

        select_by_coord.main(["-i", str(shared_datadir / "simple_genpept.gb"), "-o", out, "-c", coords])

        source = {record.id: record for record in SeqIO.parse(str(shared_datadir / "simple_genpept.gb"), "genbank")}
        records = list(parse_seqfiles([out]))
        assert len(records) == 2
        assert records[0].id == "pDONR201_3_1:3"
        assert str(records[0].seq) == str(source["pDONR201_3"].seq[:3])
        assert records[1].id == "pDONR201_3_4:6"
        assert str(records[1].seq) == str(source["pDONR201_3"].seq[3:6])


def test_select_by_coord_reverse_strand_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/extraction.gb"
        coords = output_dir + "/coords.tsv"
        _write_table(coords, ["contig", "start", "end", "strand"], [("pDONR201_3", 1, 6, "-")])

        select_by_coord.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "-c", coords])

        source = next(record for record in SeqIO.parse(str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "genbank") if record.id == "pDONR201_3")
        records = list(parse_seqfiles([out]))
        assert len(records) == 1
        assert records[0].id == "pDONR201_3_1:6rc"
        assert str(records[0].seq) == str(source.seq[:6].reverse_complement())


def test_select_by_coord_circular_origin_span_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/extraction.gb"
        coords = output_dir + "/coords.tsv"
        _write_table(coords, ["contig", "start", "end"], [("pDONR201", 4468, 3)])

        select_by_coord.main(["-i", str(shared_datadir / "pDONR201_domainator_circular.gb"), "-o", out, "-c", coords])

        source = next(SeqIO.parse(str(shared_datadir / "pDONR201_domainator_circular.gb"), "genbank"))
        expected = str(source.seq[4467:] + source.seq[:3])
        records = list(parse_seqfiles([out]))
        assert len(records) == 1
        assert records[0].id == "pDONR201_4468:3"
        assert str(records[0].seq) == expected
        assert records[0].features[0].type == "source"


def test_select_by_coord_circular_origin_span_reverse_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/extraction.gb"
        coords = output_dir + "/coords.tsv"
        _write_table(coords, ["contig", "start", "end", "strand"], [("pDONR201", 4468, 3, "reverse")])

        select_by_coord.main(["-i", str(shared_datadir / "pDONR201_domainator_circular.gb"), "-o", out, "-c", coords])

        source = next(SeqIO.parse(str(shared_datadir / "pDONR201_domainator_circular.gb"), "genbank"))
        expected = str((source.seq[4467:] + source.seq[:3]).reverse_complement())
        records = list(parse_seqfiles([out]))
        assert len(records) == 1
        assert records[0].id == "pDONR201_4468:3rc"
        assert str(records[0].seq) == expected


def test_select_by_coord_linear_origin_span_raises_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/extraction.gb"
        coords = output_dir + "/coords.tsv"
        _write_table(coords, ["contig", "start", "end"], [("pDONR201_3", 10, 5)])

        with pytest.raises(ValueError, match="only valid for circular contigs"):
            select_by_coord.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "-c", coords])


def test_select_by_coord_file_basename_filter_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/extraction.gb"
        coords = output_dir + "/coords.tsv"
        _write_table(
            coords,
            ["contig", "start", "end", "file"],
            [("pDONR201_3", 1, 3, "simple_genpept.gb")],
        )

        select_by_coord.main([
            "-i",
            str(shared_datadir / "simple_genpept.gb"),
            str(shared_datadir / "pdonr_peptides.fasta"),
            "-o",
            out,
            "-c",
            coords,
        ])

        records = list(parse_seqfiles([out]))
        assert len(records) == 1
        assert records[0].id == "pDONR201_3_1:3"


def test_select_by_coord_missing_file_matches_all_files_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/extraction.gb"
        coords = output_dir + "/coords.tsv"
        _write_table(coords, ["contig", "start", "end"], [("pDONR201_3", 1, 3)])

        select_by_coord.main([
            "-i",
            str(shared_datadir / "simple_genpept.gb"),
            str(shared_datadir / "pdonr_peptides.fasta"),
            "-o",
            out,
            "-c",
            coords,
        ])

        records = list(parse_seqfiles([out]))
        assert len(records) == 2


def test_select_by_coord_rejects_extra_columns_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/extraction.gb"
        coords = output_dir + "/coords.tsv"
        _write_table(coords, ["contig", "start", "end", "extra"], [("pDONR201_3", 1, 3, "x")])

        with pytest.raises(ValueError, match="Unexpected coordinate table columns"):
            select_by_coord.main(["-i", str(shared_datadir / "simple_genpept.gb"), "-o", out, "-c", coords])