import tempfile
from domainator.Bio import SeqIO
from domainator import extract_unannotated
from domainator.utils import parse_seqfiles, DomainatorCDS
import pytest
from io import StringIO
import sys
import subprocess


def test_extract_unannotated_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        extract_unannotated.main(["-i", str(shared_datadir / "simple_genpept.gb"), "-o", out])
        # assert 0
        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 8
        [str(seq) for seq in seqs] == ["FPALSPDSVDNRITASQEEFVETQKGHPSGWPSA*","M","I*","*","M","DEWQGGA*","MDADLYGYKWARDNVGQSGATIYRLYGKPDAPELFLKHGKGSVANDV","DMNKLQFHLMLDEFF*"]

def test_extract_unannotated_largest_keep_name_2(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        extract_unannotated.main(["-i", str(shared_datadir / "simple_genpept.gb"), "-o", out, "--largest", "--keep_name"])
        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 5
        [str(seq) for seq in seqs] == ["FPALSPDSVDNRITASQEEFVETQKGHPSGWPSA*","I*","*","DEWQGGA*","MDADLYGYKWARDNVGQSGATIYRLYGKPDAPELFLKHGKGSVANDV"]

def test_extract_unannotated_lb_2(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        extract_unannotated.main(["-i", str(shared_datadir / "simple_genpept.gb"), "-o", out, "--lb", "10"])
        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 3
        [str(seq) for seq in seqs] == ["FPALSPDSVDNRITASQEEFVETQKGHPSGWPSA*", "MDADLYGYKWARDNVGQSGATIYRLYGKPDAPELFLKHGKGSVANDV", "DMNKLQFHLMLDEFF*"]