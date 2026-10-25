"""Tests for output files: temporary names and permissions."""
import os
import stat
import sys

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from domainator.output_guardrails import make_temporary_output_path
from domainator.outputs import output_file, write_seq_outputs

posix_only = pytest.mark.skipif(sys.platform == "win32", reason="POSIX file modes")


@pytest.fixture
def umask_022():
    old = os.umask(0o022)
    yield
    os.umask(old)


def _mode(path):
    return stat.S_IMODE(os.stat(path).st_mode)


def test_make_temporary_output_path(tmp_path):
    target = tmp_path / "out.gb"
    first = make_temporary_output_path(str(target))
    second = make_temporary_output_path(str(target))
    assert first != second
    for path in (first, second):
        assert os.path.exists(path)
        assert os.path.dirname(path) == str(tmp_path)
        assert os.path.basename(path).startswith(".out.gb.")
        assert path.endswith(".tmp")


@posix_only
def test_output_file_mode_follows_umask(tmp_path, umask_022):
    path = tmp_path / "report.tsv"
    with output_file(str(path)) as handle:
        handle.write("a\tb\n")
    assert _mode(path) == 0o644


@posix_only
def test_seq_outputs_mode_follows_umask(tmp_path, umask_022):
    record = SeqRecord(Seq("ACGT"), id="r1", name="r1", description="")
    genbank, fasta = tmp_path / "out.gb", tmp_path / "out.fasta"
    write_seq_outputs([record], genbank=str(genbank), fasta=str(fasta))
    assert _mode(genbank) == 0o644
    assert _mode(fasta) == 0o644
