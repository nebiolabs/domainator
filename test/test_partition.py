import os
import pytest
from domainator.offsets import i_partition_seqfiles, partition_seqfile
from domainator.seqio import parse_seqfiles

# (file, cdss_per_partition, total CDSs, partition offsets, records per partition)
PARTITION_CASES = [
    ("pDONR201.gb", 1, 3, [0], [1]),
    ("pDONR201_multigenemark_partition.gb", 12, 24, [0, 16382], [2, 2]),
    ("pDONR201_multigenemark_partition.gb", 24, 24, [0], [4]),
    ("pDONR201_multigenemark_partition.gb", 2, 24, [0, 8191, 16382, 24573], [1, 1, 1, 1]),
    ("pDONR201_multigenemark_partition.gb", 6, 24, [0, 8191, 16382, 24573], [1, 1, 1, 1]),
    ("pDONR201_multigenemark_partition.gb", 7, 24, [0, 16382], [2, 2]),
]


@pytest.mark.parametrize("file,cdss_per_partition,num_proteins,offsets,recs_to_read", PARTITION_CASES)
def test_partition_seqfile(file, cdss_per_partition, num_proteins, offsets, recs_to_read, shared_datadir):
    path = str(shared_datadir / file)
    cds_count, partitions = partition_seqfile(path, cdss_per_partition)
    assert cds_count == num_proteins
    assert [p[0] for p in partitions] == [path] * len(offsets)
    assert [p[1] for p in partitions] == offsets
    assert [p[2] for p in partitions] == recs_to_read


@pytest.mark.parametrize("file,cdss_per_partition,num_proteins,offsets,recs_to_read", PARTITION_CASES)
def test_i_partition_seqfiles_matches_partition_seqfile(file, cdss_per_partition, num_proteins, offsets, recs_to_read, shared_datadir):
    path = str(shared_datadir / file)
    _, expected = partition_seqfile(path, cdss_per_partition)
    assert list(i_partition_seqfiles([path, path], cdss_per_partition)) == expected + expected


def test_partitions_read_every_record_once(shared_datadir):
    path = str(shared_datadir / "pDONR201_multigenemark_partition.gb")
    _, partitions = partition_seqfile(path, 7)
    ids = [rec.id for _, offset, n in partitions for rec in parse_seqfiles([path], seek_to=offset, max_recs=n)]
    assert ids == [rec.id for rec in parse_seqfiles([path])]


def test_partition_seqfile_rejects_empty_file(tmp_path):
    empty = tmp_path / "empty.fasta"
    empty.write_text("")
    with pytest.raises(ValueError, match="No sequence records"):
        partition_seqfile(str(empty), 1)


def test_gbfast_offsets_match_python(shared_datadir):
    """If the native _gbfast extension is built, its offsets must exactly match
    the pure-Python scanner (so partitions are identical)."""
    import glob
    from domainator.offsets import get_fasta_offsets, get_genbank_offsets
    from domainator.offsets import _gbfast as gbfast
    if gbfast is None:
        import pytest
        pytest.skip("native _gbfast extension not built")
    checked = 0
    for f in sorted(glob.glob(str(shared_datadir / "*.gb"))):
        try:
            off, n = get_genbank_offsets(f)
        except Exception:
            continue
        assert list(zip(off, n)) == gbfast.genbank_offsets(f), f
        checked += 1
    for f in sorted(glob.glob(str(shared_datadir / "*.fasta"))) + sorted(glob.glob(str(shared_datadir / "*.fna"))):
        try:
            off, n = get_fasta_offsets(f)
        except Exception:
            continue
        assert list(zip(off, n)) == gbfast.fasta_offsets(f), f
        checked += 1
    assert checked > 0


def test_gbfast_bgzf_offsets_match_python(shared_datadir, tmp_path):
    """The native BGZF offset scanner must produce the exact same virtual offsets
    (block_start<<16 | within) as the Python BgzfReader scanner, so partitions of
    a BGZF database are identical to the uncompressed ones."""
    import glob
    from domainator.offsets import get_fasta_offsets, get_genbank_offsets, get_offsets
    from domainator.seqio import detect_compression
    from domainator.Bio import bgzf
    from domainator.offsets import _gbfast as gbfast
    if gbfast is None or not hasattr(gbfast, "genbank_offsets_bgzf"):
        pytest.skip("native _gbfast BGZF scanner not built")

    def bgzip(src):
        dst = tmp_path / (os.path.basename(src) + ".bgz")
        with open(src, "rb") as fh, bgzf.BgzfWriter(str(dst)) as w:
            w.write(fh.read())
        return str(dst)

    checked = 0
    for f in sorted(glob.glob(str(shared_datadir / "*.gb"))):
        try:
            py = list(zip(*get_genbank_offsets(bgz := bgzip(f))))
        except Exception:
            continue
        assert detect_compression(bgz) == "bgzf"
        assert py == gbfast.genbank_offsets_bgzf(bgz), f
        # get_offsets must dispatch to the native scanner and agree.
        off, n = get_offsets(bgz)
        assert list(zip(off, n)) == py, f
        checked += 1
    for f in sorted(glob.glob(str(shared_datadir / "*.fasta"))) + sorted(glob.glob(str(shared_datadir / "*.fna"))):
        try:
            py = list(zip(*get_fasta_offsets(bgz := bgzip(f))))
        except Exception:
            continue
        assert py == gbfast.fasta_offsets_bgzf(bgz), f
        checked += 1
    assert checked > 0


def test_gbfast_bgzf_offsets_multiblock(shared_datadir, tmp_path):
    """Records spanning many BGZF blocks: offsets must still match the Python
    scanner exactly, and every offset must seek to a record boundary."""
    from domainator.offsets import get_genbank_offsets
    from domainator.Bio import bgzf
    from domainator.offsets import _gbfast as gbfast
    if gbfast is None or not hasattr(gbfast, "genbank_offsets_bgzf"):
        pytest.skip("native _gbfast BGZF scanner not built")

    src = shared_datadir / "MT_nbs.gb"  # 20 records; repeat to span many blocks
    data = src.read_bytes()
    bgz = tmp_path / "multiblock.gb.bgz"
    with bgzf.BgzfWriter(str(bgz)) as w:
        for _ in range(50):
            w.write(data)

    py = list(zip(*get_genbank_offsets(str(bgz))))
    rust = gbfast.genbank_offsets_bgzf(str(bgz))
    assert len({o >> 16 for o, _ in py}) > 1, "test file did not span multiple BGZF blocks"
    assert py == rust

    reader = bgzf.BgzfReader(str(bgz), "rb")
    try:
        for vo, _ in rust:
            reader.seek(int(vo))
            assert reader.readline().startswith(b"LOCUS "), vo
    finally:
        reader.close()
