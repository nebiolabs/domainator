from domainator import domainator_db_download
import tempfile
import gzip
import shutil
import pytest
from domainator import utils
from pathlib import Path


def test_gbff_is_valid(shared_datadir):
    src_gb = shared_datadir / "bacillus_phage_SPR.gb"
    with tempfile.TemporaryDirectory() as output_dir:
        # valid gzip (direct-backend style: <accession>_<name>_genomic.gbff.gz)
        good = Path(output_dir) / "good.gb.gz"
        with open(src_gb, "rb") as src, gzip.open(good, "wb") as dst:
            shutil.copyfileobj(src, dst)
        assert domainator_db_download.gbff_is_valid(good) is True

        # valid uncompressed (NCBI Datasets style: genomic.gbff)
        plain = Path(output_dir) / "genomic.gbff"
        shutil.copyfile(src_gb, plain)
        assert domainator_db_download.gbff_is_valid(plain) is True

        # truncated / corrupt gzip
        truncated = Path(output_dir) / "bad.gb.gz"
        truncated.write_bytes(good.read_bytes()[: good.stat().st_size // 2])
        assert domainator_db_download.gbff_is_valid(truncated) is False


def test_is_gzip_and_open_maybe_gzip(shared_datadir):
    src_gb = shared_datadir / "bacillus_phage_SPR.gb"
    expected = src_gb.read_text()
    with tempfile.TemporaryDirectory() as output_dir:
        gz = Path(output_dir) / "genome.gbff.gz"
        with open(src_gb, "rb") as src, gzip.open(gz, "wb") as dst:
            shutil.copyfileobj(src, dst)
        plain = Path(output_dir) / "genomic.gbff"
        shutil.copyfile(src_gb, plain)

        assert domainator_db_download._is_gzip(gz) is True
        assert domainator_db_download._is_gzip(plain) is False

        # both open transparently to the same text content
        with domainator_db_download._open_maybe_gzip(gz, "rt") as fh:
            assert fh.read() == expected
        with domainator_db_download._open_maybe_gzip(plain, "rt") as fh:
            assert fh.read() == expected


def test_process_local_gbff_no_gene_call(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        gbff_gz = Path(output_dir) / "genome.gbff.gz"
        with open(shared_datadir / "bacillus_phage_SPR.gb", "rb") as src, gzip.open(gbff_gz, "wb") as dst:
            shutil.copyfileobj(src, dst)
        outfile = Path(output_dir) / "out.gb"
        outfile.touch()

        ok = domainator_db_download.process_local_gbff(str(gbff_gz), str(outfile), gene_call=None)
        assert ok is True
        recs = list(utils.parse_seqfiles([str(outfile)]))
        assert len(recs) == 1
        assert any(f.type == "CDS" for f in recs[0].features)


def test_process_local_gbff_skips_corrupt_gzip(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        bad = Path(output_dir) / "genome.gbff.gz"
        bad.write_bytes(b"not really gzip")
        outfile = Path(output_dir) / "out.gb"
        outfile.touch()
        ok = domainator_db_download.process_local_gbff(str(bad), str(outfile), gene_call=None)
        assert ok is False
        assert outfile.stat().st_size == 0



def test_uniprot_download_fasta(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        outfile = Path(output_dir) / "uniprot_sprot.fasta"
        domainator_db_download.main(["--db", "swissprot", "--num_recs", "2", "--output", str(outfile)])
        assert outfile.exists()
        assert outfile.stat().st_size > 0
        recs = list(utils.parse_seqfiles([str(outfile)]))
        assert len(recs) == 2

def test_uniprot_download_genbank(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        outfile = Path(output_dir) / "uniprot_sprot.gb"
        domainator_db_download.main(["--db", "swissprot_gb", "--num_recs", "2", "--output", str(outfile)])
        assert outfile.exists()
        assert outfile.stat().st_size > 0
        recs = list(utils.parse_seqfiles([str(outfile)]))
        assert len(recs) == 2



def test_genbank_download_genbank_1(shared_datadir):
    small_genbanks=['https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/766/775/GCA_008766775.1_ASM876677v1','https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/031/580/395/GCA_031580395.1_ASM3158039v1']
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        outfile = Path(output_dir) / "gb.gb"
        domainator_db_download.process_genbank_accessions([{'ftp_path':small_genbank} for small_genbank in small_genbanks], outfile, gene_call=None, num_recs=1, download_workers=3)
        assert outfile.exists()
        # read output file
        recs = list(utils.parse_seqfiles([str(outfile)]))
        assert len(recs) == 1

def test_genbank_download_genbank_2(shared_datadir):
    small_genbanks=['https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/766/775/GCA_008766775.1_ASM876677v1','https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/031/580/395/GCA_031580395.1_ASM3158039v1']
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        outfile = Path(output_dir) / "gb.gb"
        domainator_db_download.process_genbank_accessions([{'ftp_path':small_genbank} for small_genbank in small_genbanks], outfile, gene_call=None, num_recs=None, download_workers=3)
        assert outfile.exists()
        # read output file
        recs = list(utils.parse_seqfiles([str(outfile)]))
        assert len(recs) == 2

def test_genbank_download_genbank_3(shared_datadir):
    small_genbanks=['https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/766/775/GCA_008766775.1_ASM876677v1', 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/987/885/GCA_002987885.1_ASM298788v1']
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        outfile = Path(output_dir) / "gb.gb"
        domainator_db_download.process_genbank_accessions([{'ftp_path':small_genbank} for small_genbank in small_genbanks], outfile, gene_call="all", num_recs=None, download_workers=2)
        assert outfile.exists()
        # read output file
        recs = list(utils.parse_seqfiles([str(outfile)]))
        assert len(recs) == 2
        outfile_text = outfile.read_text()
        assert "CDS" in outfile_text
        assert '/gene_id="AM260465_1"' in outfile_text

# def test_genbank_download_genbank_skipped_record_log(shared_datadir):
#     small_genbanks=['https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/766/775/GCA_008766775.1_ASM87667v1', 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/987/885/GCA_002987885.1_ASM298788v1']
#     with tempfile.TemporaryDirectory() as output_dir:
#         #output_dir = "test_out"
#         outfile = Path(output_dir) / "gb.gb"
#         skipped_record_log = Path(output_dir) / "skipped_record_log.txt"
#         domainator_db_download.process_genbank_accessions([{'ftp_path':small_genbank} for small_genbank in small_genbanks], outfile, gene_call="all", num_recs=None, cpus=2, skipped_record_log=skipped_record_log)
#         assert outfile.exists()
#         # read output file
#         recs = list(utils.parse_seqfiles([str(outfile)]))
#         assert len(recs) == 1
#         outfile_text = outfile.read_text()
#         assert "CDS" in outfile_text
#         assert '/gene_id="AM260465_1"' in outfile_text
#         assert skipped_record_log.exists()
#         skipped_record_log_text = skipped_record_log.read_text()
#         assert "GCA_008766775.1_ASM87667v1" in skipped_record_log_text
#         assert "GCA_002987885.1_ASM298788v1" not in skipped_record_log_text

def test_process_local_gbff_uncompressed(shared_datadir):
    """NCBI Datasets writes an uncompressed genomic.gbff; process_local_gbff must handle it
    (regression: the datasets backend previously assumed gzip and produced no output)."""
    src_gb = shared_datadir / "bacillus_phage_SPR.gb"
    with tempfile.TemporaryDirectory() as output_dir:
        plain = Path(output_dir) / "genomic.gbff"   # uncompressed, datasets-style name
        shutil.copyfile(src_gb, plain)
        outfile = Path(output_dir) / "out.gb"
        outfile.touch()

        # fast path (no gene calling): records copied through verbatim
        ok = domainator_db_download.process_local_gbff(str(plain), str(outfile), gene_call=None)
        assert ok is True
        recs = list(utils.parse_seqfiles([str(outfile)]))
        assert len(recs) == 1
        assert any(f.type == "CDS" for f in recs[0].features)


def test_process_local_gbff_uncompressed_gene_call(shared_datadir):
    """gene_call='unannotated' over an uncompressed gbff (the mode the update scripts use)."""
    src_gb = shared_datadir / "bacillus_phage_SPR.gb"
    with tempfile.TemporaryDirectory() as output_dir:
        plain = Path(output_dir) / "genomic.gbff"
        shutil.copyfile(src_gb, plain)
        outfile = Path(output_dir) / "out.gb"
        outfile.touch()
        ok = domainator_db_download.process_local_gbff(str(plain), str(outfile), gene_call="unannotated")
        assert ok is True
        recs = list(utils.parse_seqfiles([str(outfile)]))
        assert len(recs) == 1
        assert any(f.type == "CDS" for f in recs[0].features)


def _make_fake_rehydrated_tree(workdir, src_gb, accessions):
    """Create workdir/extracted/ncbi_dataset/data/<accession>/genomic.gbff for each accession,
    mimicking what `datasets rehydrate` produces (uncompressed genomic.gbff per assembly)."""
    data_dir = Path(workdir) / "extracted" / "ncbi_dataset" / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    for acc in accessions:
        d = data_dir / acc
        d.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src_gb, d / "genomic.gbff")


def _mock_datasets_cli(monkeypatch):
    """Make the datasets/unzip subprocess calls no-op successes so process_via_datasets can be
    exercised offline against a pre-populated rehydrated tree."""
    monkeypatch.setattr(domainator_db_download.shutil, "which", lambda name: "/usr/bin/" + name)

    class _Result:
        returncode = 0

    monkeypatch.setattr(domainator_db_download.subprocess, "run", lambda *a, **k: _Result())


def test_process_via_datasets_offline(shared_datadir, monkeypatch):
    """End-to-end of the datasets backend with the CLI mocked: the glob must match the real
    datasets layout (<accession>/genomic.gbff), the uncompressed files must be appended, and
    the success log must be consistent with what was written (no appended-but-unlogged genome,
    which would duplicate on resume)."""
    src_gb = shared_datadir / "bacillus_phage_SPR.gb"
    accs = ["GCA_000000011.1", "GCA_000000022.1"]
    _mock_datasets_cli(monkeypatch)

    with tempfile.TemporaryDirectory() as td:
        workdir = Path(td) / "workdir"
        workdir.mkdir(parents=True, exist_ok=True)
        _make_fake_rehydrated_tree(workdir, src_gb, accs)  # what `datasets rehydrate` would produce

        outfile = Path(td) / "out.gb"
        success_log = Path(td) / "success_new.txt"
        genbank_accessions = [{"assembly_accession": a} for a in accs]

        summary = domainator_db_download.process_via_datasets(
            genbank_accessions, str(outfile), gene_call=None, num_recs=None, cpus=1,
            success_rec_log=str(success_log), datasets_workdir=str(workdir),
            datasets_include="gbff", datasets_max_workers=1, api_key=None,
            datasets_path="datasets",
        )

        assert summary["appended"] == 2
        assert summary["failed"] == 0
        assert summary["missing"] == 0

        recs = list(utils.parse_seqfiles([str(outfile)]))
        assert len(recs) == 2

        logged = set(success_log.read_text().split())
        assert logged == {"GCA_000000011", "GCA_000000022"}  # logged == appended (consistent)


def test_process_via_datasets_num_recs_slices_download(shared_datadir, monkeypatch):
    """--num_recs caps the *download*: only num_recs accessions are written to the datasets
    input file, and only those are appended/logged."""
    src_gb = shared_datadir / "bacillus_phage_SPR.gb"
    accs = ["GCA_000000011.1", "GCA_000000022.1", "GCA_000000033.1"]
    _mock_datasets_cli(monkeypatch)

    with tempfile.TemporaryDirectory() as td:
        workdir = Path(td) / "workdir"
        workdir.mkdir(parents=True, exist_ok=True)
        _make_fake_rehydrated_tree(workdir, src_gb, accs[:1])  # rehydrate would fetch only the first

        outfile = Path(td) / "out.gb"
        success_log = Path(td) / "success_new.txt"
        genbank_accessions = [{"assembly_accession": a} for a in accs]

        summary = domainator_db_download.process_via_datasets(
            genbank_accessions, str(outfile), gene_call=None, num_recs=1, cpus=1,
            success_rec_log=str(success_log), datasets_workdir=str(workdir),
            datasets_include="gbff", datasets_max_workers=1, api_key=None,
            datasets_path="datasets",
        )

        requested = (workdir / "accessions.txt").read_text().split()
        assert requested == ["GCA_000000011.1"]   # download capped to num_recs
        assert summary["appended"] == 1
        assert summary["missing"] == 0
        assert success_log.read_text().split() == ["GCA_000000011"]


def test_success_log_equals_exclude_file_raises():
    """--success_rec_log (truncated each run) must not be the same file as
    --exclude_accessions_file, or the canonical exclusion log would be erased."""
    with tempfile.TemporaryDirectory() as td:
        out = Path(td) / "out.fasta"
        log = Path(td) / "log.txt"
        with pytest.raises(Exception, match="must be different files"):
            domainator_db_download.main([
                "--db", "swissprot", "--output", str(out),
                "--success_rec_log", str(log),
                "--exclude_accessions_file", str(log),
            ])
        assert not out.exists()  # guard fires before the output is created/blanked
