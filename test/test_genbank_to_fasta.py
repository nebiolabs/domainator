"""Tests for domainator genbank_to_fasta."""
from Bio import SeqIO

from domainator import genbank_to_fasta
from domainator.seqio import parse_seqfiles


def test_function_and_command_write_the_same_fasta(shared_datadir, tmp_path):
    gb = str(shared_datadir / "pDONR201_multi_genemark_domainator.gb")
    from_function, from_command = tmp_path / "f.fasta", tmp_path / "c.fasta"
    genbank_to_fasta.genbank_to_fasta([gb], fasta_out=str(from_function))
    genbank_to_fasta.main(["-i", gb, "-o", str(from_command)])
    assert from_function.read_text() == from_command.read_text()
    records = list(SeqIO.parse(from_function, "fasta"))
    expected = list(parse_seqfiles([gb]))
    assert [(r.id, str(r.seq)) for r in records] == [(r.id, str(r.seq)) for r in expected]


def test_standard_output(shared_datadir, capsys):
    genbank_to_fasta.main(["-i", str(shared_datadir / "simple_genpept.gb"), "-o", "-"])
    assert capsys.readouterr().out.startswith(">")
