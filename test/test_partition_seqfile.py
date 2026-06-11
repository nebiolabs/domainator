from math import prod
from domainator.partition_seqfile import main
import os
import tempfile
from glob import glob
from io import StringIO
import pytest

@pytest.mark.parametrize("file,option,value,num_proteins,offsets,recs_to_read",
[("pDONR201.gb","--partitions",2,3,[0],[1]),
("pDONR201.gb","--cdss_per_partition",1,3,[0],[1]),
("pDONR201_multigenemark_partition.gb","--partitions",2,24,[0,16382],[2,2]),
("pDONR201_multigenemark_partition.gb","--partitions",1,24,[0],[4]),
("pDONR201_multigenemark_partition.gb","--partitions",10,24,[0,8191,16382,24573],[1,1,1,1]),
("pDONR201_multigenemark_partition.gb","--partitions",4,24,[0,8191,16382,24573],[1,1,1,1]),
("pDONR201_multigenemark_partition.gb","--cdss_per_partition",6,24,[0,8191,16382,24573],[1,1,1,1]),
("pDONR201_multigenemark_partition.gb","--cdss_per_partition",7,24,[0,16382],[2,2]),
])
def test_partition_seqfile(file, num_proteins, option, value, offsets, recs_to_read, shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        outfile = output_dir + "/" + "outfile.txt"
        args = ["-i", str(shared_datadir / file), "-o", outfile, option, str(value)]
        main(args)
        with open(outfile) as f:
            proteins = int(f.readline().strip())
            assert proteins == num_proteins
            produced_offsets = list()
            produced_recs_to_read = list()
            for line in f:
                parts = line.strip().split()
                produced_offsets.append(int(parts[0]))
                produced_recs_to_read.append(int(parts[1]))
            assert produced_offsets == offsets
            assert produced_recs_to_read == recs_to_read


def test_gbfast_offsets_match_python(shared_datadir):
    """If the native _gbfast extension is built, its offsets must exactly match
    the pure-Python scanner (so partitions are identical)."""
    import glob
    from domainator import utils
    gbfast = getattr(utils, "_gbfast", None)
    if gbfast is None:
        import pytest
        pytest.skip("native _gbfast extension not built")
    checked = 0
    for f in sorted(glob.glob(str(shared_datadir / "*.gb"))):
        try:
            off, n = utils.get_genbank_offsets(f)
        except Exception:
            continue
        assert list(zip(off, n)) == gbfast.genbank_offsets(f), f
        checked += 1
    for f in sorted(glob.glob(str(shared_datadir / "*.fasta"))) + sorted(glob.glob(str(shared_datadir / "*.fna"))):
        try:
            off, n = utils.get_fasta_offsets(f)
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
    from domainator import utils
    from domainator.Bio import bgzf
    gbfast = getattr(utils, "_gbfast", None)
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
            py = list(zip(*utils.get_genbank_offsets(bgz := bgzip(f))))
        except Exception:
            continue
        assert utils.detect_compression(bgz) == "bgzf"
        assert py == gbfast.genbank_offsets_bgzf(bgz), f
        # get_offsets must dispatch to the native scanner and agree.
        off, n = utils.get_offsets(bgz)
        assert list(zip(off, n)) == py, f
        checked += 1
    for f in sorted(glob.glob(str(shared_datadir / "*.fasta"))) + sorted(glob.glob(str(shared_datadir / "*.fna"))):
        try:
            py = list(zip(*utils.get_fasta_offsets(bgz := bgzip(f))))
        except Exception:
            continue
        assert py == gbfast.fasta_offsets_bgzf(bgz), f
        checked += 1
    assert checked > 0


def test_gbfast_bgzf_offsets_multiblock(shared_datadir, tmp_path):
    """Records spanning many BGZF blocks: offsets must still match the Python
    scanner exactly, and every offset must seek to a record boundary."""
    from domainator import utils
    from domainator.Bio import bgzf
    gbfast = getattr(utils, "_gbfast", None)
    if gbfast is None or not hasattr(gbfast, "genbank_offsets_bgzf"):
        pytest.skip("native _gbfast BGZF scanner not built")

    src = shared_datadir / "MT_nbs.gb"  # 20 records; repeat to span many blocks
    data = src.read_bytes()
    bgz = tmp_path / "multiblock.gb.bgz"
    with bgzf.BgzfWriter(str(bgz)) as w:
        for _ in range(50):
            w.write(data)

    py = list(zip(*utils.get_genbank_offsets(str(bgz))))
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
