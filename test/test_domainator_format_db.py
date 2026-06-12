"""Tests for domainator_format_db: BGZF compression, sharding, indexing, the
compose pipeline, idempotency, and end-to-end use by domain_search."""
import glob
import gzip
import os
import shutil
import warnings

import pytest

from domainator import domainator_format_db, utils, db_index, bgzf_compress
from domainator.domain_search import main as domain_search_main

native_compress = pytest.mark.skipif(
    not bgzf_compress.native_available(), reason="native BGZF compressor not built"
)


def _records(path):
    """(id, seq) for every record in a database file, via the faithful parser."""
    recs = list(utils.parse_seqfiles([path], genbank_parser="biopython"))
    return [(r.id, str(r.seq)) for r in recs]


def _shard_files(prefix_glob):
    files = glob.glob(prefix_glob)
    return sorted(files, key=lambda p: int(os.path.basename(p).split(".")[1]))


# --- BGZF compressor -----------------------------------------------------------

@native_compress
def test_native_compress_roundtrip_and_ordering(shared_datadir, tmp_path):
    """Native parallel compression must round-trip exactly (block order preserved)
    at multiple thread counts on a multi-block input."""
    data = (shared_datadir / "Staph_phages.gb").read_bytes()
    big = tmp_path / "big.gb"
    with open(big, "wb") as f:
        for _ in range(8):              # several MB -> many 64 KiB BGZF blocks
            f.write(data)
    orig = big.read_bytes()
    for threads in (1, 4):
        out = tmp_path / f"big_{threads}.gb.bgz"
        bgzf_compress.compress_to_bgzf(str(big), str(out), threads=threads)
        assert utils.detect_compression(str(out)) == "bgzf"
        with gzip.open(str(out), "rb") as fh:
            assert fh.read() == orig


def test_compress_fallback_matches_source(shared_datadir, tmp_path, monkeypatch):
    """The Python BgzfWriter fallback must also produce valid BGZF of the source."""
    data = (shared_datadir / "MT_nbs.gb").read_bytes()
    src = tmp_path / "f.gb"
    src.write_bytes(data)
    monkeypatch.setattr(bgzf_compress, "_gbfast", None)   # force fallback
    out = tmp_path / "fb.gb.bgz"
    bgzf_compress.compress_to_bgzf(str(src), str(out))
    assert utils.detect_compression(str(out)) == "bgzf"
    with gzip.open(str(out), "rb") as fh:
        assert fh.read() == data


# --- sharding ------------------------------------------------------------------

def test_shard_reassembles_to_original(shared_datadir, tmp_path):
    """Byte-exact shards concatenate back to the original record set, in order."""
    src = tmp_path / "r.gb"
    shutil.copy(shared_datadir / "MT_nbs.gb", src)
    original = _records(str(src))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")   # "wrote shards alongside original"
        domainator_format_db.format_db([str(src)], shards=4, compress=False, index=False, cpu=1)
    shards = _shard_files(str(tmp_path / "r.*.gb"))
    assert len(shards) == 4
    reassembled = []
    for s in shards:
        reassembled += _records(s)
    assert reassembled == original


def test_compose_shard_compress_index(shared_datadir, tmp_path):
    """plain .gb -> --shards 3 --compress --index produces 3 valid BGZF shards each
    with a usable .didx, covering all records."""
    src = tmp_path / "db.gb"
    shutil.copy(shared_datadir / "MT_nbs.gb", src)
    original = _records(str(src))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        domainator_format_db.format_db([str(src)], shards=3, compress=True, index=True, cpu=2)
    os.remove(src)                         # so the shard resolver is unambiguous

    shards = db_index.resolve_database_shards(str(src))
    assert len(shards) == 3
    reassembled = []
    for s in shards:
        assert utils.detect_compression(s) == "bgzf"
        assert os.path.exists(db_index.index_path_for(s))
        # get_offsets here goes through the freshly written .didx
        assert len(utils.get_offsets(s)[0]) > 0
        reassembled += _records(s)
    assert sorted(reassembled) == sorted(original)


def test_idempotency(shared_datadir, tmp_path):
    src = tmp_path / "idem.gb"
    shutil.copy(shared_datadir / "MT_nbs.gb", src)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        domainator_format_db.format_db([str(src)], shards=2, compress=True, index=True, cpu=1)
        shards = _shard_files(str(tmp_path / "idem.*.gb.bgz"))
        assert len(shards) == 2

        for s in shards:                   # mark as old
            os.utime(s, ns=(0, 0))
        domainator_format_db.format_db([str(src)], shards=2, compress=True, index=True, cpu=1, force=False)
        assert all(os.stat(s).st_mtime_ns == 0 for s in shards), "outputs rewritten without --force"

        domainator_format_db.format_db([str(src)], shards=2, compress=True, index=True, cpu=1, force=True)
        assert all(os.stat(s).st_mtime_ns != 0 for s in shards), "--force did not rewrite"


def test_index_only_in_place(shared_datadir, tmp_path):
    """No shard, no compress: the input is indexed in place (not rewritten)."""
    src = tmp_path / "inplace.gb"
    shutil.copy(shared_datadir / "MT_nbs.gb", src)
    before = os.stat(src).st_mtime_ns
    os.utime(src, ns=(123456789, 123456789))
    marked = os.stat(src).st_mtime_ns
    domainator_format_db.format_db([str(src)], shards=None, compress=False, index=True, cpu=1)
    assert os.path.exists(db_index.index_path_for(str(src)))
    assert os.stat(src).st_mtime_ns == marked, "input file was rewritten for index-only run"
    # and the index is valid/fresh
    assert utils.get_offsets(str(src))[0] is not None


def test_nothing_to_do_raises(shared_datadir, tmp_path):
    """No operation requested (no shard/compress/index) is an error, not a no-op."""
    src = tmp_path / "noop.gb"
    shutil.copy(shared_datadir / "MT_nbs.gb", src)
    with pytest.raises(ValueError, match="[Nn]othing to do"):
        domainator_format_db.format_db([str(src)], shards=None, compress=False, index=False)


def test_name_renames_sharded_output(shared_datadir, tmp_path):
    """--name renames the output database; format extension is preserved, and the
    renamed shards reassemble to the original records."""
    src = tmp_path / "foo.gb"
    shutil.copy(shared_datadir / "MT_nbs.gb", src)
    original = _records(str(src))
    domainator_format_db.format_db([str(src)], shards=2, compress=True, index=True,
                                   name="mydb", cpu=1)
    shards = _shard_files(str(tmp_path / "mydb.*.gb.bgz"))
    assert [os.path.basename(s) for s in shards] == ["mydb.0.gb.bgz", "mydb.1.gb.bgz"]
    assert all(os.path.exists(db_index.index_path_for(s)) for s in shards)
    # The original (foo.gb) is untouched and not shadowed (different name).
    assert os.path.exists(src)
    reassembled = []
    for s in shards:
        reassembled += _records(s)
    assert sorted(reassembled) == sorted(original)


def test_name_index_only_warns_and_copies(shared_datadir, tmp_path):
    """--name on an index-only run warns, then copies the input to the new name and
    indexes the copy (an index alone cannot rename a database)."""
    src = tmp_path / "foo.gb"
    shutil.copy(shared_datadir / "MT_nbs.gb", src)
    with pytest.warns(RuntimeWarning, match="index-only"):
        domainator_format_db.format_db([str(src)], index=True, name="renamed", cpu=1)
    out = tmp_path / "renamed.gb"
    assert out.exists() and os.path.exists(db_index.index_path_for(str(out)))
    assert _records(str(out)) == _records(str(src))


def test_name_multiple_inputs_raises(shared_datadir, tmp_path):
    a = tmp_path / "a.gb"; b = tmp_path / "b.gb"
    shutil.copy(shared_datadir / "MT_nbs.gb", a)
    shutil.copy(shared_datadir / "MT_nbs.gb", b)
    with pytest.raises(ValueError, match="multiple inputs"):
        domainator_format_db.format_db([str(a), str(b)], compress=True, name="x", cpu=1)


# --- end-to-end with domain_search --------------------------------------------

def test_domain_search_over_shards_matches_unsharded(shared_datadir, tmp_path):
    """Searching a logical (sharded) database must return the same hits as searching
    the equivalent unsharded file."""
    one = (shared_datadir / "pDONR201.gb").read_bytes()
    hmm = str(shared_datadir / "CcdB.hmm")

    unsharded = tmp_path / "multi.gb"
    with open(unsharded, "wb") as f:
        for _ in range(4):                 # 4 records, each with a CcdB hit
            f.write(one)

    sharded = tmp_path / "sdb.gb"
    shutil.copy(unsharded, sharded)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        domainator_format_db.format_db([str(sharded)], shards=2, compress=True, index=True, cpu=1)
    os.remove(sharded)                     # leave only the shards

    common = ["-r", hmm, "--evalue", "0.1", "--max_overlap", "1", "-Z", "1000", "--add_annotations"]
    out_u = str(tmp_path / "out_u.gb")
    out_s = str(tmp_path / "out_s.gb")
    domain_search_main(["--input", str(unsharded), "-o", out_u] + common)
    domain_search_main(["--input", str(sharded), "-o", out_s] + common)   # resolves to sdb.0/.1

    recs_u = sorted(_records(out_u))
    recs_s = sorted(_records(out_s))
    assert len(recs_u) > 0
    assert recs_u == recs_s


def test_domain_search_Z0_uses_index_count(shared_datadir, tmp_path):
    """-Z 0 over an indexed sharded database (true target count read from the .didx
    headers) must match -Z 0 over the equivalent unsharded file."""
    one = (shared_datadir / "pDONR201.gb").read_bytes()
    hmm = str(shared_datadir / "CcdB.hmm")

    unsharded = tmp_path / "z0.gb"
    with open(unsharded, "wb") as f:
        for _ in range(4):
            f.write(one)
    sharded = tmp_path / "z0s.gb"
    shutil.copy(unsharded, sharded)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        domainator_format_db.format_db([str(sharded)], shards=2, index=True, cpu=1)
    os.remove(sharded)

    common = ["-r", hmm, "--evalue", "0.1", "--max_overlap", "1", "-Z", "0", "--add_annotations"]
    out_u = str(tmp_path / "z0_u.gb")
    out_s = str(tmp_path / "z0_s.gb")
    domain_search_main(["--input", str(unsharded), "-o", out_u] + common)
    domain_search_main(["--input", str(sharded), "-o", out_s] + common)
    recs_u = sorted(_records(out_u))
    assert len(recs_u) > 0
    assert recs_u == sorted(_records(out_s))
