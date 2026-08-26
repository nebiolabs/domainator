"""Tests for structure_search.py, the structure-input counterpart of domain_search.py."""
import shutil

import pytest

from domainator import structure_search, utils


pytestmark = pytest.mark.skipif(shutil.which("foldseek") is None, reason="foldseek is not installed")

# 1UBQ is a copy of the reference; 1NDD (NEDD8) is a genuine remote homolog of it.
EXPECTED_HITS = {"1UBQ_A", "1NDD_A"}


def _read(path):
    return list(utils.parse_seqfiles((str(path),), default_molecule_type="protein"))


def _features(record, feature_type):
    return [f for f in record.features if f.type == feature_type]


@pytest.fixture
def database(shared_datadir):
    return str(shared_datadir / "structures" / "inputs")


@pytest.fixture
def references(shared_datadir):
    return str(shared_datadir / "structures" / "refs")


@pytest.fixture
def hits(database, references, tmp_path):
    out = tmp_path / "hits.gb"
    structure_search.main(["-i", database, "-r", references, "-o", str(out), "--cpu", "2"])
    return _read(out)


def test_only_hit_records_are_written(hits):
    assert {record.id for record in hits} == EXPECTED_HITS


def test_hits_carry_domain_search_features(hits):
    for record in hits:
        assert len(_features(record, "Domain_Search")) == 1
        assert not _features(record, "Domainator"), (
            "Domainator features are only written with --add_annotations")


def test_records_are_built_from_backend_sequences(hits):
    """The searched database is never read end to end: sequences come from the hit rows."""
    lengths = {record.id: len(record) for record in hits}
    assert lengths == {"1UBQ_A": 76, "1NDD_A": 74}
    for record in hits:
        assert record.annotations["molecule_type"] == "protein"
        assert str(record.seq)


def test_hit_coordinates_are_on_the_database_record(hits):
    """The query/target inversion test: the feature must describe the database record and
    /rlen the reference, which differ in length for the NEDD8 fixture."""
    by_id = {record.id: record for record in hits}
    for record in hits:
        feature = _features(record, "Domain_Search")[0]
        start, end = int(feature.location.start), int(feature.location.end)
        assert 0 <= start < end <= len(record)
        assert int(feature.qualifiers["rlen"][0]) == 76

    nedd8 = by_id["1NDD_A"]
    assert len(nedd8) == 74
    assert int(_features(nedd8, "Domain_Search")[0].qualifiers["rlen"][0]) == 76


def test_identity_is_a_percentage(hits):
    by_id = {record.id: record for record in hits}
    identity = float(_features(by_id["1UBQ_A"], "Domain_Search")[0].qualifiers["identity"][0])
    assert identity == pytest.approx(100.0, abs=0.1)


def test_add_annotations(database, references, tmp_path):
    out = tmp_path / "annotated.gb"
    structure_search.main(["-i", database, "-r", references, "-o", str(out),
                           "--add_annotations", "--cpu", "2"])
    for record in _read(out):
        assert len(_features(record, "Domain_Search")) == 1
        assert len(_features(record, "Domainator")) >= 1


def test_max_hits_keeps_the_best_scoring_record(database, references, tmp_path):
    out = tmp_path / "one.gb"
    structure_search.main(["-i", database, "-r", references, "-o", str(out),
                           "--max_hits", "1", "--cpu", "2"])
    records = _read(out)
    assert len(records) == 1
    # the exact copy of the reference must outscore the remote homolog
    assert records[0].id == "1UBQ_A"


def test_max_hits_sorts_best_first(database, references, tmp_path):
    out = tmp_path / "sorted.gb"
    structure_search.main(["-i", database, "-r", references, "-o", str(out),
                           "-e", "10", "--max_hits", "10", "--cpu", "2"])
    records = _read(out)
    scores = [max(float(f.qualifiers["score"][0]) for f in _features(r, "Domain_Search"))
              for r in records]
    assert scores == sorted(scores, reverse=True)


def test_min_evalue_excludes_the_identical_match(database, references, tmp_path):
    out = tmp_path / "no_self.gb"
    structure_search.main(["-i", database, "-r", references, "-o", str(out),
                           "--min_evalue", "1e-15", "--cpu", "2"])
    assert {record.id for record in _read(out)} == {"1NDD_A"}


def test_evalue_threshold_changes_hit_count(database, references, tmp_path):
    strict = tmp_path / "strict.gb"
    loose = tmp_path / "loose.gb"
    structure_search.main(["-i", database, "-r", references, "-o", str(strict),
                           "-e", "1e-15", "--cpu", "2"])
    structure_search.main(["-i", database, "-r", references, "-o", str(loose),
                           "-e", "10", "--cpu", "2"])
    assert len(_read(strict)) < len(_read(loose))


def test_prebuilt_database_input(shared_datadir, references, tmp_path):
    """The AlphaFold-DB shape: a prebuilt database searched with reference structures."""
    out = tmp_path / "fesod_hits.gb"
    structure_search.main(["-i", str(shared_datadir / "foldseek" / "FeSOD"),
                           "-r", references, "-o", str(out), "-e", "10", "--cpu", "2"])
    records = _read(out)
    assert records, "expected some hits at -e 10"
    for record in records:
        assert str(record.seq), "sequences must come from the backend, not the database"
        assert len(_features(record, "Domain_Search")) == 1


def test_database_prefix_input_matches_structure_file_input(database, references, tmp_path):
    from_files = tmp_path / "files.gb"
    prefix = tmp_path / "indb"
    structure_search.main(["-i", database, "-r", references, "-o", str(from_files),
                           "--keep_db", str(prefix), "--cpu", "2"])
    from_db = tmp_path / "db.gb"
    structure_search.main(["-i", str(prefix), "-r", references, "-o", str(from_db), "--cpu", "2"])
    assert ([(r.id, str(r.seq)) for r in _read(from_files)]
            == [(r.id, str(r.seq)) for r in _read(from_db)])


def test_metrics_and_alignment_type_1(database, references, tmp_path):
    out = tmp_path / "metrics.gb"
    structure_search.main(["-i", database, "-r", references, "-o", str(out),
                           "--alignment_type", "1", "-e", "10",
                           "--metrics", "tmscore", "lddt", "--cpu", "2"])
    by_id = {record.id: record for record in _read(out)}
    feature = _features(by_id["1UBQ_A"], "Domain_Search")[0]
    assert float(feature.qualifiers["tmscore"][0]) > 0.95
    assert float(feature.qualifiers["lddt"][0]) == pytest.approx(1.0, abs=0.01)


def test_hits_tsv(database, references, tmp_path):
    """The raw table exists so results can be cross-checked against a hand-run search."""
    out = tmp_path / "hits.gb"
    tsv = tmp_path / "raw.tsv"
    structure_search.main(["-i", database, "-r", references, "-o", str(out),
                           "--hits_tsv", str(tsv), "--cpu", "2"])
    lines = [line for line in tsv.read_text().split("\n") if line]
    assert len(lines) == 2
    targets = {line.split("\t")[1] for line in lines}
    assert targets == EXPECTED_HITS


def test_max_seqs_saturation_warns(database, references, tmp_path):
    """Silent truncation would read as 'this reference has few homologs'."""
    out = tmp_path / "capped.gb"
    with pytest.warns(RuntimeWarning, match=r"--max_seqs \(1\) was reached"):
        structure_search.main(["-i", database, "-r", references, "-o", str(out),
                               "-e", "10", "--max_seqs", "1", "--cpu", "2"])


def test_max_output_gb(database, references, tmp_path):
    out = tmp_path / "toobig.gb"
    with pytest.raises(SystemExit, match="max_output_gb"):
        structure_search.main(["-i", database, "-r", references, "-o", str(out),
                               "--max_output_gb", "0.000000001", "--cpu", "2"])
    assert not out.exists()


def test_max_hits_must_be_positive(database, references, tmp_path):
    with pytest.raises(RuntimeError, match="--max_hits must be at least 1"):
        structure_search.main(["-i", database, "-r", references,
                               "-o", str(tmp_path / "x.gb"), "--max_hits", "0"])


def test_downstream_reports(database, references, tmp_path):
    """Hit records must be consumable by the rest of the suite unchanged."""
    from domainator import enum_report

    out = tmp_path / "hits.gb"
    structure_search.main(["-i", database, "-r", references, "-o", str(out),
                           "--add_annotations", "--cpu", "2"])
    report = tmp_path / "report.tsv"
    enum_report.main(["-i", str(out), "--domains", "-o", str(report)])
    text = report.read_text()
    assert "UBQref_A" in text
