"""Tests for the .didx offset-index format, staleness handling, and shard
resolution (db_index.py), plus the get_offsets/i_get_offsets index fast-path."""
import os
import shutil
import warnings
from array import array

import pytest

from domainator import utils, db_index
from domainator.Bio import bgzf

GB_FILES = ["pDONR201.gb", "MT_nbs.gb", "206.gb"]
FASTA_FILES = ["FeSOD_20.fasta"]


def _bgzip(src, dst):
    with open(src, "rb") as fh, bgzf.BgzfWriter(str(dst)) as w:
        w.write(fh.read())
    return str(dst)


def test_index_roundtrip_genbank(shared_datadir, tmp_path):
    for fn in GB_FILES:
        f = tmp_path / fn
        shutil.copy(shared_datadir / fn, f)
        off, npr = utils.get_genbank_offsets(str(f))
        db_index.write_index(str(f), off, npr, filetype="genbank", compression=None)
        roff, rnpr = db_index.read_index(str(f), filetype="genbank", compression=None)
        assert isinstance(roff, array) and roff.typecode == "Q"
        assert isinstance(rnpr, array) and rnpr.typecode == "Q"
        assert list(roff) == list(off) and list(rnpr) == list(npr)
        streamed = list(db_index.i_read_index(str(f), filetype="genbank", compression=None))
        assert streamed == list(zip(off, npr))


def test_index_roundtrip_fasta(shared_datadir, tmp_path):
    for fn in FASTA_FILES:
        f = tmp_path / fn
        shutil.copy(shared_datadir / fn, f)
        off, npr = utils.get_fasta_offsets(str(f))
        db_index.write_index(str(f), off, npr, filetype="fasta", compression=None)
        roff, rnpr = db_index.read_index(str(f), filetype="fasta", compression=None)
        assert list(roff) == list(off) and list(rnpr) == list(npr)


def test_read_total_cds_from_header(shared_datadir, tmp_path):
    """The total CDS count is available from the header in O(1), matches the body
    sum, and is None when there is no usable/fresh index."""
    f = tmp_path / "MT_nbs.gb"
    shutil.copy(shared_datadir / "MT_nbs.gb", f)
    off, npr = utils.get_offsets(str(f))
    db_index.write_index(str(f), off, npr, filetype="genbank", compression=None)
    assert db_index.read_total_cds(str(f), filetype="genbank", compression=None) == sum(npr)
    assert utils.index_total_cds(str(f)) == sum(npr)
    # stale source -> None (no count trusted)
    with open(f, "ab") as fh:
        fh.write(b"\n")
    db_index._warned_index_paths.clear()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert db_index.read_total_cds(str(f), filetype="genbank", compression=None) is None
    # no index at all
    assert utils.index_total_cds(str(tmp_path / "absent.gb")) is None


def test_get_offsets_uses_index(shared_datadir, tmp_path):
    """With a valid index present, get_offsets/i_get_offsets must return exactly the
    same arrays/tuples as a fresh scan."""
    f = tmp_path / "pDONR201.gb"
    shutil.copy(shared_datadir / "pDONR201.gb", f)
    off, npr = utils.get_offsets(str(f))               # fresh scan
    db_index.write_index(str(f), off, npr, filetype="genbank", compression=None)
    off2, npr2 = utils.get_offsets(str(f))             # via index
    assert list(off2) == list(off) and list(npr2) == list(npr)
    assert list(utils.i_get_offsets(str(f))) == list(zip(off, npr))


def test_index_bgzf_virtual_offsets(shared_datadir, tmp_path):
    """BGZF index offsets are virtual offsets that seek to record boundaries."""
    bgz = _bgzip(shared_datadir / "MT_nbs.gb", tmp_path / "MT_nbs.gb.bgz")
    off, npr = utils.get_offsets(bgz)                  # BGZF scan (Rust or Python)
    db_index.write_index(bgz, off, npr, filetype="genbank", compression="bgzf")
    roff, rnpr = db_index.read_index(bgz, filetype="genbank", compression="bgzf")
    assert list(roff) == list(off) and list(rnpr) == list(npr)
    reader = bgzf.BgzfReader(bgz, "rb")
    try:
        for vo in roff:
            reader.seek(int(vo))
            assert reader.readline().startswith(b"LOCUS "), vo
    finally:
        reader.close()


def test_index_stale_on_append(shared_datadir, tmp_path):
    f = tmp_path / "stale_append.gb"
    shutil.copy(shared_datadir / "MT_nbs.gb", f)
    off, npr = utils.get_offsets(str(f))
    db_index.write_index(str(f), off, npr, filetype="genbank", compression=None)
    with open(f, "ab") as fh:                          # changes size + mtime
        fh.write(b"\n")
    db_index._warned_index_paths.clear()
    with pytest.warns(RuntimeWarning, match="stale index"):
        assert db_index.read_index(str(f), filetype="genbank", compression=None) is None
    # get_offsets must transparently recompute correct (current) offsets.
    off2, npr2 = utils.get_offsets(str(f))
    assert len(off2) == len(off)


def test_index_corrupt_returns_none(shared_datadir, tmp_path):
    f = tmp_path / "corrupt.gb"
    shutil.copy(shared_datadir / "pDONR201.gb", f)
    off, npr = utils.get_offsets(str(f))
    idx = db_index.write_index(str(f), off, npr, filetype="genbank", compression=None)

    # truncated body
    with open(idx, "r+b") as fh:
        fh.truncate(db_index._HEADER_STRUCT.size + 4)
    db_index._warned_index_paths.clear()
    with pytest.warns(RuntimeWarning):
        assert db_index.read_index(str(f), filetype="genbank", compression=None) is None

    # bad magic
    db_index.write_index(str(f), off, npr, filetype="genbank", compression=None)
    with open(idx, "r+b") as fh:
        fh.seek(0); fh.write(b"XXXX")
    db_index._warned_index_paths.clear()
    with pytest.warns(RuntimeWarning):
        assert db_index.read_index(str(f), filetype="genbank", compression=None) is None

    # wrong file-kind flag (built genbank, queried as fasta) -> None, no raise
    db_index.write_index(str(f), off, npr, filetype="genbank", compression=None)
    db_index._warned_index_paths.clear()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert db_index.read_index(str(f), filetype="fasta", compression=None) is None


def test_index_missing_returns_none_silently(shared_datadir, tmp_path):
    f = tmp_path / "noidx.gb"
    shutil.copy(shared_datadir / "pDONR201.gb", f)
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # absence must NOT warn
        assert db_index.read_index(str(f), filetype="genbank", compression=None) is None


def test_write_index_rejects_gzip(shared_datadir, tmp_path):
    f = tmp_path / "x.gb"
    shutil.copy(shared_datadir / "pDONR201.gb", f)
    off, npr = utils.get_offsets(str(f))
    with pytest.raises(ValueError):
        db_index.write_index(str(f), off, npr, filetype="genbank", compression="gzip")


# --- shard resolution ----------------------------------------------------------

def _touch(path):
    with open(path, "w") as f:
        f.write("x")


def test_resolve_only_unsharded(tmp_path):
    p = tmp_path / "mydb.gb"; _touch(p)
    assert db_index.resolve_database_shards(str(p)) == [str(p)]


@pytest.mark.filterwarnings("ignore:Database shards")  # gap (10) warning tested separately
def test_resolve_only_shards_numeric_order(tmp_path):
    for i in (0, 1, 2, 10):                            # 10 must sort after 2, not before
        _touch(tmp_path / f"mydb.{i}.gb")
    resolved = db_index.resolve_database_shards(str(tmp_path / "mydb.gb"))
    assert [os.path.basename(p) for p in resolved] == ["mydb.0.gb", "mydb.1.gb", "mydb.2.gb", "mydb.10.gb"]


def test_resolve_both_present_raises(tmp_path):
    _touch(tmp_path / "mydb.gb")
    _touch(tmp_path / "mydb.0.gb")
    with pytest.raises(ValueError, match="[Bb]oth"):
        db_index.resolve_database_shards(str(tmp_path / "mydb.gb"))


def test_resolve_neither_returns_path(tmp_path):
    p = tmp_path / "missing.gb"
    assert db_index.resolve_database_shards(str(p)) == [str(p)]


def test_resolve_leading_zero_collision_raises(tmp_path):
    _touch(tmp_path / "mydb.0.gb")
    _touch(tmp_path / "mydb.00.gb")
    with pytest.raises(ValueError, match="collision"):
        db_index.resolve_database_shards(str(tmp_path / "mydb.gb"))


def test_resolve_gap_warns(tmp_path):
    for i in (0, 1, 3):
        _touch(tmp_path / f"mydb.{i}.gb")
    with pytest.warns(RuntimeWarning, match="contiguous"):
        resolved = db_index.resolve_database_shards(str(tmp_path / "mydb.gb"))
    assert len(resolved) == 3


def test_resolve_mixed_compression(tmp_path):
    _touch(tmp_path / "mydb.0.gb")
    _touch(tmp_path / "mydb.1.gb.bgz")
    resolved = db_index.resolve_database_shards(str(tmp_path / "mydb.gb"))
    assert [os.path.basename(p) for p in resolved] == ["mydb.0.gb", "mydb.1.gb.bgz"]
