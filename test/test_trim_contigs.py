import tempfile
from Bio import SeqIO
from domainator import trim_contigs
from domainator.seqio import parse_seqfiles
from domainator.model import DomainatorCDS
import pytest
from io import StringIO
import sys
import subprocess


def test_trim_contigs_cds_both_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        trim_contigs.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--contigs", "pDONR201_1", "--cds_range", "2",])

        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        for record in new_file:
            assert len(record) == 1000
            assert record.name == "pDONR201_1_1266:2265"

def test_trim_contigs_cds_both_2(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        trim_contigs.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--contigs", "pDONR201_1", "--cds_range", "3",])

        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 0

def test_trim_contigs_cds_both_domains_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        trim_contigs.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--domains", "APH", "--contigs", "pDONR201_1", "--cds_range", "1",])

        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        for record in new_file:
            assert len(record) == 1307
            assert record.name == "pDONR201_1_959:2265"


def test_trim_contigs_domain_expr_2(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        trim_contigs.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--domain_expr", "(APH & TCAD9) | (CcdB)", "--contigs", "pDONR201_1", "--cds_range", "1",])

        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        for record in new_file:
            assert len(record) == 1000
            assert record.name == "pDONR201_1_1266:2265"

def test_trim_contigs_no_domain_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        trim_contigs.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--domain_expr", "(APH & TCAD9) | (CcdB)", "--contigs", "pDONR201_1", "--unannotated"])

        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        for record in new_file:
            assert len(record) == 1000
            assert record.name == "pDONR201_1_1266:2265"


def test_trim_contigs_kb_both_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        trim_contigs.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--contigs", "pDONR201_1", "--kb_range", "2",])

        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        for record in new_file:
            assert len(record) == 470
            assert record.name == "pDONR201_1_2001:2470"

def test_trim_contigs_kb_both_2(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        trim_contigs.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--contigs", "pDONR201_1", "--kb_range", "4",])

        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 0


def test_trim_contigs_kb_both_peptides_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        trim_contigs.main(["-i", str(shared_datadir / "pdonr_peptides.fasta"), "-o", out, "--contigs", "pDONR201_2", "--kb_range", "0.010",])

        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        assert len(new_file[0]) == 82




@pytest.mark.parametrize("filters,expected_name", [
    ([], "pDONR201_1_959:2265"),
    (["--databases", "pdonr_hmms"], "pDONR201_1_959:2265"),
    (["--databases", "other_db"], "pDONR201_1_959:3677"),  # APH is not from other_db, so it is not stripped
    (["-e", "1e-20"], "pDONR201_1_959:2265"),
    (["-e", "1e-30"], "pDONR201_1_959:3677"),  # APH has E-value 1.9e-28, so it is not stripped
])
def test_trim_contigs_domain_filters(shared_datadir, tmp_path, filters, expected_name):
    out = str(tmp_path / "out.gb")
    trim_contigs.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out,
                       "--domains", "APH", "--contigs", "pDONR201_1", "--cds_range", "1"] + filters)
    assert [record.name for record in SeqIO.parse(out, "genbank")] == [expected_name]


def test_trim_contigs_unannotated_alone(shared_datadir, tmp_path):
    """--unannotated on its own strips the unannotated CDSs at both ends."""
    out = str(tmp_path / "out.gb")
    trim_contigs.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out,
                       "--contigs", "pDONR201_1", "--unannotated"])
    assert [record.name for record in SeqIO.parse(out, "genbank")] == ["pDONR201_1_959:3677"]


def test_trim_contigs_cds_range_without_cdss(tmp_path):
    """A contig with no CDS annotations is left alone by --cds_range."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from domainator.seqio import write_genbank
    in_gb = str(tmp_path / "in.gb")
    write_genbank([SeqRecord(Seq("ACGT" * 50), id="no_cds", annotations={"molecule_type": "DNA"})], in_gb)
    out = str(tmp_path / "out.gb")
    trim_contigs.main(["-i", in_gb, "-o", out, "--cds_range", "1"])
    records = list(SeqIO.parse(out, "genbank"))
    assert len(records) == 1 and len(records[0]) == 200
