import tempfile
from domainator.Bio import SeqIO
from domainator import trim_contigs
from domainator.utils import parse_seqfiles, DomainatorCDS
import pytest
from io import StringIO
import sys
import subprocess


def test_trim_contigs_cds_both_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        trim_contigs.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--contigs", "pDONR201_1", "--cds_both", "2",])

        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        for record in new_file:
            assert len(record) == 1000
            assert record.name == "pDONR201_1_1266:2265"

def test_trim_contigs_cds_both_2(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        trim_contigs.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--contigs", "pDONR201_1", "--cds_both", "3",])

        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 0

def test_trim_contigs_cds_both_domains_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        trim_contigs.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--domains", "APH", "--contigs", "pDONR201_1", "--cds_both", "1",])

        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        for record in new_file:
            assert len(record) == 1307
            assert record.name == "pDONR201_1_959:2265"


def test_trim_contigs_domain_expr_2(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        trim_contigs.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--domain_expr", "(APH & TCAD9) | (CcdB)", "--contigs", "pDONR201_1", "--cds_both", "1",])

        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        for record in new_file:
            assert len(record) == 1000
            assert record.name == "pDONR201_1_1266:2265"

def test_trim_contigs_no_domain_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        trim_contigs.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--domain_expr", "(APH & TCAD9) | (CcdB)", "--contigs", "pDONR201_1", "--no_domain"])

        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        for record in new_file:
            assert len(record) == 1000
            assert record.name == "pDONR201_1_1266:2265"


def test_trim_contigs_kb_both_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        trim_contigs.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--contigs", "pDONR201_1", "--kb_both", "2",])

        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        for record in new_file:
            assert len(record) == 470
            assert record.name == "pDONR201_1_2001:2470"

def test_trim_contigs_kb_both_2(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        trim_contigs.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--contigs", "pDONR201_1", "--kb_both", "4",])

        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 0


def test_trim_contigs_kb_both_peptides_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        trim_contigs.main(["-i", str(shared_datadir / "pdonr_peptides.fasta"), "-o", out, "--contigs", "pDONR201_2", "--kb_both", "0.010",])

        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        assert len(new_file[0]) == 82


