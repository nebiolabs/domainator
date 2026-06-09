"""Tests for domainator.foldseek.search, including CPU and GPU (CUDA) modes.

The end-to-end domainate path converts proteins to 3Di with a 3B parameter
ESM2->3Di model, which is too heavy to run in CI. These tests instead exercise
foldseek.search directly, using the precomputed 3Di sequences in
test/data/foldseek/FeSOD_20.3di.fasta, so the model is not required.

GPU mode requires a recent foldseek, a CUDA capable GPU, and a GPU-compatible
(padded) target database. The padded database is rebuilt from the committed
FeSOD foldseek database with `foldseek makepaddedseqdb` (so it always matches the
installed foldseek version), and the GPU test is skipped when no GPU is present.
"""
import os
import shutil
import subprocess

import pytest
import pyhmmer

from domainator import foldseek as foldseek_lib


def _has_foldseek():
    return shutil.which("foldseek") is not None


def _has_cuda_gpu():
    try:
        import torch
        return torch.cuda.is_available()
    except Exception:
        return False


pytestmark = pytest.mark.skipif(not _has_foldseek(), reason="foldseek is not installed")


def _read_fasta(path):
    recs = []
    header = None
    seq = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    recs.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line)
    if header is not None:
        recs.append((header, "".join(seq)))
    return recs


def _build_search_inputs(prot_fasta, threedi_fasta):
    """Build the (proteins, foldseek) inputs that foldseek.search expects.

    Mirrors what domainate's foldseekBuilder produces at runtime, but reads
    precomputed 3Di sequences (index-aligned with the protein fasta) so the
    heavy ESM2->3Di model is not needed. fasta2foldseek requires the protein and
    3Di records to share an identical header and have equal length, so we assert
    that here too.
    """
    prot_recs = _read_fasta(prot_fasta)
    threedi_recs = _read_fasta(threedi_fasta)
    assert len(prot_recs) == len(threedi_recs)

    proteins = []
    foldseek = []
    for (p_header, p_seq), (t_header, t_seq) in zip(prot_recs, threedi_recs):
        assert p_header == t_header, f"{p_header!r} != {t_header!r}"
        assert len(p_seq) == len(t_seq)
        name, _, desc = p_header.partition(" ")
        proteins.append(
            pyhmmer.easel.TextSequence(
                name=name.encode(), description=desc.encode(), sequence=p_seq
            ).digitize(pyhmmer.easel.Alphabet.amino())
        )
        foldseek.append(f">{t_header}\n{t_seq}")
    return proteins, foldseek


def _make_padded_db(src_db, dest_db):
    """Build a GPU-compatible (padded) foldseek database from src_db.

    Note: foldseek's makepaddedseqdb runs a generated shell script whose cleanup
    step can exit non-zero even on success, so we verify the GPU mapping file was
    produced rather than trusting the return code.
    """
    proc = subprocess.run(
        ["foldseek", "makepaddedseqdb", str(src_db), str(dest_db)],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    gpu_mapping = str(dest_db) + "_ss.gpu_mapping1"
    assert os.path.exists(gpu_mapping), (
        f"padded database was not built (missing {gpu_mapping}):\n"
        + proc.stdout.decode("utf-8", errors="replace")
    )
    return dest_db


@pytest.fixture
def fesod_inputs(shared_datadir):
    return _build_search_inputs(
        shared_datadir / "foldseek" / "FeSOD_20.fasta",
        shared_datadir / "foldseek" / "FeSOD_20.3di.fasta",
    )


def _name(field):
    """foldseek query/target fields may carry a description; take the accession."""
    return field.split()[0]


def _check_self_hits(hits):
    """Each of the 20 FeSOD sequences should find itself in the database."""
    queries = {_name(h.query) for h in hits}
    assert len(queries) == 20
    self_hits = {_name(h.query) for h in hits if _name(h.query) == _name(h.target)}
    assert len(self_hits) == 20, f"expected 20 self hits, got {len(self_hits)}"


def test_foldseek_search_cpu(fesod_inputs, shared_datadir):
    proteins, foldseek = fesod_inputs
    db = shared_datadir / "foldseek" / "FeSOD"

    hits = list(
        foldseek_lib.search(str(db), proteins=proteins, foldseek=foldseek, cpu=1, E=0.001, device="cpu")
    )

    assert len(hits) >= 20
    _check_self_hits(hits)
    # every reported hit must respect the E-value threshold
    assert all(float(h.evalue) <= 0.001 for h in hits)


def test_foldseek_search_device_none_matches_cpu(fesod_inputs, shared_datadir):
    """device=None (the default) must behave the same as an explicit "cpu"."""
    proteins, foldseek = fesod_inputs
    db = shared_datadir / "foldseek" / "FeSOD"

    hits_none = list(
        foldseek_lib.search(str(db), proteins=proteins, foldseek=foldseek, cpu=1, E=0.001, device=None)
    )
    hits_cpu = list(
        foldseek_lib.search(str(db), proteins=proteins, foldseek=foldseek, cpu=1, E=0.001, device="cpu")
    )

    pairs_none = sorted((_name(h.query), _name(h.target)) for h in hits_none)
    pairs_cpu = sorted((_name(h.query), _name(h.target)) for h in hits_cpu)
    assert pairs_none == pairs_cpu


@pytest.mark.skipif(not _has_cuda_gpu(), reason="no CUDA GPU available for foldseek GPU search")
def test_foldseek_search_gpu_matches_cpu(fesod_inputs, shared_datadir, tmp_path):
    proteins, foldseek = fesod_inputs
    cpu_db = shared_datadir / "foldseek" / "FeSOD"
    padded_db = _make_padded_db(cpu_db, tmp_path / "FeSOD_padded")

    cpu_hits = list(
        foldseek_lib.search(str(cpu_db), proteins=proteins, foldseek=foldseek, cpu=1, E=0.001, device="cpu")
    )
    gpu_hits = list(
        foldseek_lib.search(str(padded_db), proteins=proteins, foldseek=foldseek, cpu=1, E=0.001, device="cuda:0")
    )

    _check_self_hits(gpu_hits)
    # the GPU search must find the same query/target pairs as the CPU search
    cpu_pairs = sorted((_name(h.query), _name(h.target)) for h in cpu_hits)
    gpu_pairs = sorted((_name(h.query), _name(h.target)) for h in gpu_hits)
    assert gpu_pairs == cpu_pairs
