from domainator import db_download
import tempfile
import gzip
import shutil
import pytest
from domainator.partial import get_fragment_status
from domainator.seqio import parse_seqfiles
from domainator.errors import InputError
from pathlib import Path


def test_gbff_is_valid(shared_datadir):
    src_gb = shared_datadir / "bacillus_phage_SPR.gb"
    with tempfile.TemporaryDirectory() as output_dir:
        # valid gzip (direct-backend style: <accession>_<name>_genomic.gbff.gz)
        good = Path(output_dir) / "good.gb.gz"
        with open(src_gb, "rb") as src, gzip.open(good, "wb") as dst:
            shutil.copyfileobj(src, dst)
        assert db_download.gbff_is_valid(good) is True

        # valid uncompressed (NCBI Datasets style: genomic.gbff)
        plain = Path(output_dir) / "genomic.gbff"
        shutil.copyfile(src_gb, plain)
        assert db_download.gbff_is_valid(plain) is True

        # truncated / corrupt gzip
        truncated = Path(output_dir) / "bad.gb.gz"
        truncated.write_bytes(good.read_bytes()[: good.stat().st_size // 2])
        assert db_download.gbff_is_valid(truncated) is False


def test_is_gzip_and_open_maybe_gzip(shared_datadir):
    src_gb = shared_datadir / "bacillus_phage_SPR.gb"
    expected = src_gb.read_text()
    with tempfile.TemporaryDirectory() as output_dir:
        gz = Path(output_dir) / "genome.gbff.gz"
        with open(src_gb, "rb") as src, gzip.open(gz, "wb") as dst:
            shutil.copyfileobj(src, dst)
        plain = Path(output_dir) / "genomic.gbff"
        shutil.copyfile(src_gb, plain)

        assert db_download._is_gzip(gz) is True
        assert db_download._is_gzip(plain) is False

        # both open transparently to the same text content
        with db_download._open_maybe_gzip(gz, "rt") as fh:
            assert fh.read() == expected
        with db_download._open_maybe_gzip(plain, "rt") as fh:
            assert fh.read() == expected


def test_process_local_gbff_no_gene_call(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        gbff_gz = Path(output_dir) / "genome.gbff.gz"
        with open(shared_datadir / "bacillus_phage_SPR.gb", "rb") as src, gzip.open(gbff_gz, "wb") as dst:
            shutil.copyfileobj(src, dst)
        outfile = Path(output_dir) / "out.gb"
        outfile.touch()

        ok = db_download.process_local_gbff(str(gbff_gz), str(outfile), gene_call=None)
        assert ok is True
        recs = list(parse_seqfiles([str(outfile)]))
        assert len(recs) == 1
        assert any(f.type == "CDS" for f in recs[0].features)


def test_process_local_gbff_skips_corrupt_gzip(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        bad = Path(output_dir) / "genome.gbff.gz"
        bad.write_bytes(b"not really gzip")
        outfile = Path(output_dir) / "out.gb"
        outfile.touch()
        ok = db_download.process_local_gbff(str(bad), str(outfile), gene_call=None)
        assert ok is False
        assert outfile.stat().st_size == 0



@pytest.mark.network
def test_uniprot_download_fasta(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        outfile = Path(output_dir) / "uniprot_sprot.fasta"
        db_download.main(["--db", "swissprot", "--max_records", "2", "--fasta_out", str(outfile)])
        assert outfile.exists()
        assert outfile.stat().st_size > 0
        recs = list(parse_seqfiles([str(outfile)]))
        assert len(recs) == 2

@pytest.mark.network
def test_uniprot_download_genbank(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        outfile = Path(output_dir) / "uniprot_sprot.gb"
        db_download.main(["--db", "swissprot_gb", "--max_records", "2", "--genbank_out", str(outfile)])
        assert outfile.exists()
        assert outfile.stat().st_size > 0
        recs = list(parse_seqfiles([str(outfile)]))
        assert len(recs) == 2



@pytest.mark.network
def test_genbank_download_genbank_1(shared_datadir):
    small_genbanks=['https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/766/775/GCA_008766775.1_ASM876677v1','https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/031/580/395/GCA_031580395.1_ASM3158039v1']
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        outfile = Path(output_dir) / "gb.gb"
        db_download.process_genbank_accessions([{'ftp_path':small_genbank} for small_genbank in small_genbanks], outfile, gene_call=None, num_recs=1, download_workers=3)
        assert outfile.exists()
        # read output file
        recs = list(parse_seqfiles([str(outfile)]))
        assert len(recs) == 1

@pytest.mark.network
def test_genbank_download_genbank_2(shared_datadir):
    small_genbanks=['https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/766/775/GCA_008766775.1_ASM876677v1','https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/031/580/395/GCA_031580395.1_ASM3158039v1']
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        outfile = Path(output_dir) / "gb.gb"
        db_download.process_genbank_accessions([{'ftp_path':small_genbank} for small_genbank in small_genbanks], outfile, gene_call=None, num_recs=None, download_workers=3)
        assert outfile.exists()
        # read output file
        recs = list(parse_seqfiles([str(outfile)]))
        assert len(recs) == 2

@pytest.mark.network
def test_genbank_download_genbank_3(shared_datadir):
    small_genbanks=['https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/766/775/GCA_008766775.1_ASM876677v1', 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/987/885/GCA_002987885.1_ASM298788v1']
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        outfile = Path(output_dir) / "gb.gb"
        db_download.process_genbank_accessions([{'ftp_path':small_genbank} for small_genbank in small_genbanks], outfile, gene_call="all", num_recs=None, download_workers=2)
        assert outfile.exists()
        # read output file
        recs = list(parse_seqfiles([str(outfile)]))
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
#         db_download.process_genbank_accessions([{'ftp_path':small_genbank} for small_genbank in small_genbanks], outfile, gene_call="all", num_recs=None, cpus=2, skipped_record_log=skipped_record_log)
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
        ok = db_download.process_local_gbff(str(plain), str(outfile), gene_call=None)
        assert ok is True
        recs = list(parse_seqfiles([str(outfile)]))
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
        ok = db_download.process_local_gbff(str(plain), str(outfile), gene_call="unannotated")
        assert ok is True
        recs = list(parse_seqfiles([str(outfile)]))
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
    monkeypatch.setattr(db_download, "resolve_binary", lambda name, path=None, **kwargs: path or "/usr/bin/" + name)

    class _Result:
        returncode = 0

    monkeypatch.setattr(db_download.subprocess, "run", lambda *a, **k: _Result())


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

        summary = db_download.process_via_datasets(
            genbank_accessions, str(outfile), gene_call=None, num_recs=None, cpus=1,
            success_rec_log=str(success_log), datasets_workdir=str(workdir),
            datasets_include="gbff", datasets_max_workers=1, api_key=None,
            datasets_path="datasets",
        )

        assert summary["appended"] == 2
        assert summary["failed"] == 0
        assert summary["missing"] == 0

        recs = list(parse_seqfiles([str(outfile)]))
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

        summary = db_download.process_via_datasets(
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


def test_success_log_equals_exclude_file_raises(capsys):
    """--success_log_out (truncated each run) must not be the same file as
    --exclude_accessions_file, or the canonical exclusion log would be erased."""
    with tempfile.TemporaryDirectory() as td:
        out = Path(td) / "out.fasta"
        log = Path(td) / "log.txt"
        with pytest.raises(SystemExit):
            db_download.main([
                "--db", "swissprot", "--fasta_out", str(out),
                "--success_log_out", str(log),
                "--exclude_accessions_file", str(log),
            ])
        assert "must be different files" in capsys.readouterr().err
        assert not out.exists()  # guard fires before the output is created/blanked


def test_process_via_datasets_gzip_flag(shared_datadir, monkeypatch):
    """datasets_gzip=True must add `--gzip` to the `datasets rehydrate` command (and not otherwise)."""
    src_gb = shared_datadir / "bacillus_phage_SPR.gb"
    accs = ["GCA_000000011.1"]

    calls = []

    class _Result:
        returncode = 0

    def _record(cmd, *a, **k):
        calls.append(list(cmd))
        return _Result()

    monkeypatch.setattr(db_download, "resolve_binary", lambda name, path=None, **kwargs: path or "/usr/bin/" + name)
    monkeypatch.setattr(db_download.subprocess, "run", _record)

    def _run(datasets_gzip):
        calls.clear()
        with tempfile.TemporaryDirectory() as td:
            workdir = Path(td) / "workdir"
            workdir.mkdir(parents=True, exist_ok=True)
            _make_fake_rehydrated_tree(workdir, src_gb, accs)
            db_download.process_via_datasets(
                [{"assembly_accession": a} for a in accs], str(Path(td) / "out.gb"),
                gene_call=None, num_recs=None, cpus=1, success_rec_log=None,
                datasets_workdir=str(workdir), datasets_include="gbff", datasets_max_workers=10,
                api_key=None, datasets_path="datasets", datasets_gzip=datasets_gzip)
        return [c for c in calls if "rehydrate" in c]

    on = _run(True)
    assert on and all("--gzip" in c for c in on)

    off = _run(False)
    assert off and not any("--gzip" in c for c in off)


# --- offline downloads from local mirrors (test/data/db_download, see its README) ---

SPROT = ["C0HLM9", "P0C6K8", "Q01457", "P80438", "P0A7Y4"]
SPROT_FRAGMENTS = ["C0HLM9", "P80438"]
TREMBL = ["A0A0U3GSY3", "Q50831", "A0A0K2RVI7"]


def _accessions(path):
    """The UniProt accessions of the records in path, in file order (the id of a GenBank record,
    the middle field of a UniProt FASTA header)."""
    return [r.id.split("|")[1] if "|" in r.id else r.id for r in parse_seqfiles([str(path)])]


@pytest.mark.parametrize("db,out_flag,out_name", [("swissprot", "--fasta_out", "out.fasta"), ("swissprot_gb", "--genbank_out", "out.gb")])
@pytest.mark.parametrize("partial", ["include", "exclude", "only"])
def test_uniprot_download_partial(shared_datadir, tmp_path, db, out_flag, out_name, partial):
    outfile = tmp_path / out_name
    db_download.main(["--db", db, out_flag, str(outfile), "--partial", partial, "--uniprot_dir", str(shared_datadir / "db_download" / "uniprot")])
    expected = {"include": SPROT,
                "exclude": [a for a in SPROT if a not in SPROT_FRAGMENTS],
                "only": SPROT_FRAGMENTS}[partial]
    assert _accessions(outfile) == expected
    if db == "swissprot_gb" and partial == "only":
        assert get_fragment_status(next(iter(parse_seqfiles([str(outfile)])))) is not None


@pytest.mark.parametrize("db,out_flag,out_name", [("uniprot", "--fasta_out", "out.fasta"), ("uniprot_gb", "--genbank_out", "out.gb")])
@pytest.mark.parametrize("max_records", ["1", "5", "6", "100"])
def test_uniprot_download_includes_swissprot_and_trembl(shared_datadir, tmp_path, db, out_flag, out_name, max_records):
    outfile = tmp_path / out_name
    db_download.main(["--db", db, out_flag, str(outfile), "--max_records", max_records, "--uniprot_dir", str(shared_datadir / "db_download" / "uniprot")])
    assert _accessions(outfile) == (SPROT + TREMBL)[:int(max_records)]


@pytest.mark.parametrize("db,out_flag,out_name", [("trembl", "--fasta_out", "out.fasta"), ("trembl_gb", "--genbank_out", "out.gb")])
def test_uniprot_download_trembl(shared_datadir, tmp_path, db, out_flag, out_name):
    outfile = tmp_path / out_name
    db_download.main(["--db", db, out_flag, str(outfile), "--uniprot_dir", str(shared_datadir / "db_download" / "uniprot")])
    assert _accessions(outfile) == TREMBL


@pytest.mark.parametrize("db,out_flag,out_name", [("uniprot", "--fasta_out", "out.fasta"), ("uniprot_gb", "--genbank_out", "out.gb")])
@pytest.mark.parametrize("taxonomy_args,expected", [
    (["--include_taxids", "2"], ["Q01457", "P80438", "P0A7Y4"]),
    (["--exclude_taxids", "2", "10239"], ["C0HLM9", "A0A0U3GSY3", "Q50831"]),
    (["--taxonomy_expr", "2157 | (10239 & ~136187)"], ["P0C6K8", "A0A0U3GSY3", "Q50831"]),
])
def test_uniprot_download_taxonomy(shared_datadir, tmp_path, db, out_flag, out_name, taxonomy_args, expected):
    outfile = tmp_path / out_name
    db_download.main(["--db", db, out_flag, str(outfile), "--uniprot_dir", str(shared_datadir / "db_download" / "uniprot"),
                      "--ncbi_taxonomy_path", str(shared_datadir / "db_download" / "taxdmp"), "--no_taxonomy_update"] + taxonomy_args)
    assert _accessions(outfile) == expected


def test_uniprot_download_append(shared_datadir, tmp_path):
    outfile = tmp_path / "out.fasta"
    uniprot_dir = shared_datadir / "db_download" / "uniprot"
    db_download.main(["--db", "trembl", "--fasta_out", str(outfile), "--uniprot_dir", str(uniprot_dir)])
    db_download.main(["--db", "swissprot", "--fasta_out", str(outfile), "--uniprot_dir", str(uniprot_dir), "--append"])
    assert _accessions(outfile) == TREMBL + SPROT


def test_uniprot_dir_missing_release_file(shared_datadir, tmp_path):
    with pytest.raises(InputError, match="uniprot_trembl.dat.gz"):
        db_download.main(["--db", "trembl_gb", "--genbank_out", str(tmp_path / "out.gb"), "--uniprot_dir", str(tmp_path)])


def test_uniprot_dir_rejects_ncbi(shared_datadir, tmp_path):
    with pytest.raises(ValueError, match="UniProt"):
        db_download.db_download(None, output=str(tmp_path / "out.gb"), db="ncbi_all", uniprot_dir=str(shared_datadir / "db_download" / "uniprot"))


# The mirror's assemblies, in summary order: GCA_000000088 is a truncated download, GCA_000000099
# has no genome file and GCA_000000101 has no ftp_path. GCA_000000055/66 are metagenomes and
# GCA_000000077/88 environmental samples, which the nonredundant databases keep however many there are.
@pytest.mark.parametrize("db,expected", [
    ("ncbi_all", ["GCA_000000011", "GCA_000000022", "GCA_000000033", "GCA_000000044", "GCA_000000055", "GCA_000000066", "GCA_000000077", "GCA_000000112"]),
    ("ncbi_representative_all", ["GCA_000000022", "GCA_000000033", "GCA_000000044", "GCA_000000112"]),
    ("ncbi_representative_proks", ["GCA_000000022", "GCA_000000033", "GCA_000000112"]),
    ("ncbi_nonredundant_all", ["GCA_000000033", "GCA_000000044", "GCA_000000022", "GCA_000000112", "GCA_000000055", "GCA_000000066", "GCA_000000077"]),
    ("ncbi_nonredundant_proks", ["GCA_000000033", "GCA_000000022", "GCA_000000112", "GCA_000000055", "GCA_000000066", "GCA_000000077"]),
    ("ncbi_complete_genome_proks", ["GCA_000000011", "GCA_000000022", "GCA_000000077", "GCA_000000112"]),
    ("ncbi_complete_genome_nonredundant_proks", ["GCA_000000022", "GCA_000000112", "GCA_000000077"]),
])
def test_ncbi_download_from_mirror(shared_datadir, tmp_path, db, expected):
    data = shared_datadir / "db_download"  # the NCBI mirror: its assembly summary names local genome directories
    ncbi = ["--download_backend", "direct", "--ncbi_summary", str(data / "assembly_summary.txt"),
            "--ncbi_taxonomy_path", str(data / "taxdmp"), "--no_taxonomy_update"]
    outfile = tmp_path / "out.gb"
    success_log = tmp_path / "success_log.txt"
    db_download.main(["--db", db, "-o", str(outfile), "--success_log_out", str(success_log)] + ncbi)
    assert success_log.read_text().split() == expected
    assert len(list(parse_seqfiles([str(outfile)]))) > 0


@pytest.mark.parametrize("filter_args,expected", [
    (["--include_taxids", "562"], ["GCA_000000011", "GCA_000000022", "GCA_000000112"]),
    (["--exclude_taxids", "562", "2759"], ["GCA_000000033", "GCA_000000055", "GCA_000000066", "GCA_000000077"]),
    (["--taxonomy_expr", "2 & ~81490"], ["GCA_000000011", "GCA_000000022", "GCA_000000033", "GCA_000000112"]),
    (["--before", "2021-06-01"], ["GCA_000000011", "GCA_000000033", "GCA_000000077"]),
    (["--after", "2021-06-01"], ["GCA_000000044", "GCA_000000055", "GCA_000000066", "GCA_000000112"]),
    (["--exclude_accessions", "GCA_000000011", "GCA_000000044"], ["GCA_000000022", "GCA_000000033", "GCA_000000055", "GCA_000000066", "GCA_000000077", "GCA_000000112"]),
    (["--max_records", "3"], ["GCA_000000011", "GCA_000000022", "GCA_000000033"]),
])
def test_ncbi_download_filters(shared_datadir, tmp_path, filter_args, expected):
    data = shared_datadir / "db_download"  # the NCBI mirror: its assembly summary names local genome directories
    ncbi = ["--download_backend", "direct", "--ncbi_summary", str(data / "assembly_summary.txt"),
            "--ncbi_taxonomy_path", str(data / "taxdmp"), "--no_taxonomy_update"]
    outfile = tmp_path / "out.gb"
    success_log = tmp_path / "success_log.txt"
    db_download.main(["--db", "ncbi_all", "-o", str(outfile), "--success_log_out", str(success_log)] + ncbi + filter_args)
    assert success_log.read_text().split() == expected


def test_ncbi_download_exclude_accessions_file(shared_datadir, tmp_path):
    data = shared_datadir / "db_download"  # the NCBI mirror: its assembly summary names local genome directories
    ncbi = ["--download_backend", "direct", "--ncbi_summary", str(data / "assembly_summary.txt"),
            "--ncbi_taxonomy_path", str(data / "taxdmp"), "--no_taxonomy_update"]
    outfile = tmp_path / "out.gb"
    success_log = tmp_path / "success_log.txt"
    exclude = shared_datadir / "db_download" / "exclude_accessions.txt"
    db_download.main(["--db", "ncbi_complete_genome_proks", "-o", str(outfile), "--success_log_out", str(success_log),
                      "--exclude_accessions_file", str(exclude)] + ncbi)
    assert success_log.read_text().split() == ["GCA_000000022", "GCA_000000077"]


@pytest.mark.parametrize("backend", ["direct", "datasets"])
def test_ncbi_download_empty_selection_succeeds(shared_datadir, tmp_path, monkeypatch, backend):
    """An incremental update with nothing new (every selected accession excluded) leaves an
    empty output and an empty success log and succeeds, without running `datasets`, which
    rejects an empty accession file."""
    data = shared_datadir / "db_download"
    def no_datasets(*args, **kwargs):
        raise AssertionError(f"datasets must not run on an empty selection: {args}")
    monkeypatch.setattr(db_download.subprocess, "run", no_datasets)
    outfile = tmp_path / "out.gb"
    success_log = tmp_path / "success_log.txt"
    db_download.main(["--db", "ncbi_all", "-o", str(outfile), "--success_log_out", str(success_log),
                      "--download_backend", backend, "--datasets_workdir", str(tmp_path / "wd"),
                      "--ncbi_summary", str(data / "assembly_summary.txt"),
                      "--ncbi_taxonomy_path", str(data / "taxdmp"), "--no_taxonomy_update",
                      "--include_taxids", "562", "--exclude_accessions",
                      "GCA_000000011", "GCA_000000022", "GCA_000000099", "GCA_000000112"])
    assert outfile.read_text() == ""
    assert success_log.read_text() == ""


@pytest.mark.parametrize("download_workers", ["1", "4"])
def test_ncbi_download_keeps_summary_order(shared_datadir, tmp_path, download_workers):
    """Genomes are appended in summary order whatever order the workers finish in, so the
    output does not depend on timing."""
    data = shared_datadir / "db_download"  # the NCBI mirror: its assembly summary names local genome directories
    ncbi = ["--download_backend", "direct", "--ncbi_summary", str(data / "assembly_summary.txt"),
            "--ncbi_taxonomy_path", str(data / "taxdmp"), "--no_taxonomy_update"]
    outfile = tmp_path / "out.gb"
    db_download.main(["--db", "ncbi_all", "-o", str(outfile), "--download_workers", download_workers] + ncbi)
    expected = []
    for source in ["pDONR201.gb", "pDONR201_no_CDSs.gb", "pDONR201_partly_CDSs.gb", "saccharomyces_extraction.gb",
                   "pDONR201_multi_genemark.gb", "bacillus_phage_SPR.gb", "ccdb.gb", "pDONR201_genemark.gb"]:
        expected += [r.id for r in parse_seqfiles([str(shared_datadir / source)])]
    assert [r.id for r in parse_seqfiles([str(outfile)])] == expected


@pytest.mark.parametrize("gene_call", ["unannotated", "all"])
def test_ncbi_download_gene_call(shared_datadir, tmp_path, gene_call):
    data = shared_datadir / "db_download"  # the NCBI mirror: its assembly summary names local genome directories
    ncbi = ["--download_backend", "direct", "--ncbi_summary", str(data / "assembly_summary.txt"),
            "--ncbi_taxonomy_path", str(data / "taxdmp"), "--no_taxonomy_update"]
    outfile = tmp_path / "out.gb"
    db_download.main(["--db", "ncbi_all", "-o", str(outfile), "--include_taxids", "562", "--gene_call", gene_call] + ncbi)
    records = list(parse_seqfiles([str(outfile)]))
    assert len(records) == 3
    assert all(any(f.type == "CDS" for f in r.features) for r in records)  # pDONR201_no_CDSs.gb gets genes called


def test_ncbi_nonredundant_all_needs_no_taxonomy_filter(shared_datadir, tmp_path):
    """The nonredundant databases load the taxonomy themselves, to compare lineages."""
    data = shared_datadir / "db_download"  # the NCBI mirror: its assembly summary names local genome directories
    ncbi = ["--download_backend", "direct", "--ncbi_summary", str(data / "assembly_summary.txt"),
            "--ncbi_taxonomy_path", str(data / "taxdmp"), "--no_taxonomy_update"]
    outfile = tmp_path / "out.gb"
    success_log = tmp_path / "success_log.txt"
    db_download.main(["--db", "ncbi_nonredundant_all", "-o", str(outfile), "--success_log_out", str(success_log), "--max_records", "2"] + ncbi)
    assert success_log.read_text().split() == ["GCA_000000033", "GCA_000000044"]


def test_uniprot_download_partial_rejects_ncbi(tmp_path):
    with pytest.raises(ValueError, match="UniProt"):
        db_download.db_download(None, output=str(tmp_path / "out.gb"), include_taxids=None, exclude_taxids=None, db="ncbi_all", gene_call=None, max_records=None, ncbi_summary=None, before=None, after=None, partial="exclude")
