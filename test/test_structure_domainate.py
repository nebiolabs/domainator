"""Tests for structure_domainate.py, the structure-input counterpart of domainate.py."""
import shutil

import pytest

from domainator import structure_domainate, utils


pytestmark = pytest.mark.skipif(shutil.which("foldseek") is None, reason="foldseek is not installed")

ALL_NAMES = ["1CRN_A", "1IGD_A", "1NDD_A", "1UBQ_A",
             "1ZNI_A", "1ZNI_B", "1ZNI_C", "1ZNI_D", "some.long.name_A"]
# 1UBQ is a copy of the reference; 1NDD (NEDD8) is a genuine remote homolog of it.
EXPECTED_HITS = {"1UBQ_A", "1NDD_A"}


def _read(path):
    return list(utils.parse_seqfiles((str(path),), default_molecule_type="protein"))


def _domains(record):
    return [f for f in record.features if f.type == "Domainator"]


@pytest.fixture
def inputs(shared_datadir):
    return str(shared_datadir / "structures" / "inputs")


@pytest.fixture
def references(shared_datadir):
    return str(shared_datadir / "structures" / "refs")


@pytest.fixture
def annotated(inputs, references, tmp_path):
    out = tmp_path / "annotated.gb"
    structure_domainate.main(["-i", inputs, "-r", references, "-o", str(out), "--cpu", "2"])
    return _read(out)


def test_every_input_record_is_written(annotated):
    assert [record.id for record in annotated] == ALL_NAMES


def test_only_homologs_are_annotated(annotated):
    annotated_ids = {record.id for record in annotated if _domains(record)}
    assert annotated_ids == EXPECTED_HITS


def test_annotation_coordinates_are_on_the_input(annotated):
    """The load-bearing test for the query/target inversion.

    References are the query and inputs the target, so the feature location must describe
    the *input* record and /rstart, /rend, /rlen must describe the *reference*. A swap
    would leave the location referring to reference coordinates, which this catches
    because the reference (76 aa) is longer than the NEDD8 input record (74 aa).
    """
    by_id = {record.id: record for record in annotated}

    for record_id in EXPECTED_HITS:
        record = by_id[record_id]
        feature = _domains(record)[0]
        start, end = int(feature.location.start), int(feature.location.end)
        assert 0 <= start < end <= len(record), (
            f"{record_id}: feature {start}-{end} is not inside a {len(record)} aa record")
        # these fixtures align over essentially their whole length
        assert end - start >= len(record) - 3
        # the reference length is the reference's, not the record's
        assert int(feature.qualifiers["rlen"][0]) == 76

    nedd8 = by_id["1NDD_A"]
    assert len(nedd8) == 74
    assert int(_domains(nedd8)[0].qualifiers["rlen"][0]) == 76, (
        "rlen must be the reference length (76), not the input length (74)")


def test_identity_is_a_percentage(annotated):
    """foldseek reports fident as a fraction; every other domainator producer uses percent."""
    by_id = {record.id: record for record in annotated}
    self_identity = float(_domains(by_id["1UBQ_A"])[0].qualifiers["identity"][0])
    assert self_identity == pytest.approx(100.0, abs=0.1)
    homolog_identity = float(_domains(by_id["1NDD_A"])[0].qualifiers["identity"][0])
    assert 40.0 < homolog_identity < 90.0


def test_annotation_qualifiers(annotated):
    feature = _domains({r.id: r for r in annotated}["1UBQ_A"])[0]
    assert feature.qualifiers["program"] == ["foldseek"]
    assert feature.qualifiers["name"] == ["UBQref_A"]
    assert feature.qualifiers["database"] == ["refs"]
    # structural metrics are absent unless asked for
    for metric in ("tmscore", "lddt", "rmsd", "prob"):
        assert metric not in feature.qualifiers


def test_hits_only(inputs, references, tmp_path):
    out = tmp_path / "hits.gb"
    structure_domainate.main(["-i", inputs, "-r", references, "-o", str(out),
                              "--hits_only", "--cpu", "2"])
    assert {record.id for record in _read(out)} == EXPECTED_HITS


def test_database_prefix_input_matches_structure_file_input(inputs, references, tmp_path):
    from_files = tmp_path / "files.gb"
    prefix = tmp_path / "indb"
    structure_domainate.main(["-i", inputs, "-r", references, "-o", str(from_files),
                              "--keep_db", str(prefix), "--cpu", "2"])
    from_db = tmp_path / "db.gb"
    structure_domainate.main(["-i", str(prefix), "-r", references, "-o", str(from_db), "--cpu", "2"])

    def summarize(path):
        return [(r.id, str(r.seq),
                 [(f.qualifiers["name"][0], int(f.location.start), int(f.location.end))
                  for f in _domains(r)])
                for r in _read(path)]

    assert summarize(from_files) == summarize(from_db)


def test_metrics_and_alignment_type_1(inputs, references, tmp_path):
    """What the standalone structure tools buy over annotating through domainate.py:
    the databases carry coordinates, so TM-align and coordinate metrics are available."""
    out = tmp_path / "metrics.gb"
    structure_domainate.main(["-i", inputs, "-r", references, "-o", str(out),
                              "--alignment_type", "1", "-e", "10",
                              "--metrics", "tmscore", "lddt", "rmsd", "prob", "--cpu", "2"])
    by_id = {record.id: record for record in _read(out)}
    feature = _domains(by_id["1UBQ_A"])[0]
    assert float(feature.qualifiers["tmscore"][0]) > 0.95
    assert float(feature.qualifiers["lddt"][0]) == pytest.approx(1.0, abs=0.01)
    assert float(feature.qualifiers["rmsd"][0]) < 0.1
    assert float(feature.qualifiers["prob"][0]) > 0.9


def test_metrics_rejected_without_coordinates(shared_datadir, references, tmp_path):
    with pytest.raises(RuntimeError, match="C-alpha coordinates"):
        structure_domainate.main(["-i", str(shared_datadir / "foldseek" / "FeSOD"),
                                  "-r", references, "-o", str(tmp_path / "x.gb"),
                                  "--metrics", "tmscore"])


def test_evalue_threshold_changes_hit_count(inputs, references, tmp_path):
    strict = tmp_path / "strict.gb"
    loose = tmp_path / "loose.gb"
    structure_domainate.main(["-i", inputs, "-r", references, "-o", str(strict),
                              "-e", "1e-15", "--hits_only", "--cpu", "2"])
    structure_domainate.main(["-i", inputs, "-r", references, "-o", str(loose),
                              "-e", "10", "--hits_only", "--cpu", "2"])
    assert len(_read(strict)) < len(_read(loose))


def test_min_evalue_excludes_the_identical_match(inputs, references, tmp_path):
    """--min_evalue has no foldseek equivalent, so it is enforced by domainator itself."""
    out = tmp_path / "no_self.gb"
    structure_domainate.main(["-i", inputs, "-r", references, "-o", str(out),
                              "--min_evalue", "1e-15", "--hits_only", "--cpu", "2"])
    assert {record.id for record in _read(out)} == {"1NDD_A"}


def test_database_name_override(inputs, references, tmp_path):
    out = tmp_path / "named.gb"
    structure_domainate.main(["-i", inputs, "-r", references, "-o", str(out),
                              "--database_name", "mystructures", "--cpu", "2"])
    databases = {f.qualifiers["database"][0] for r in _read(out) for f in _domains(r)}
    assert databases == {"mystructures"}


def test_no_annotations_requires_hits_only(inputs, references, tmp_path):
    with pytest.raises(RuntimeError, match="only makes sense together with --hits_only"):
        structure_domainate.main(["-i", inputs, "-r", references,
                                  "-o", str(tmp_path / "x.gb"), "--no_annotations"])


def test_max_output_gb(inputs, references, tmp_path):
    out = tmp_path / "toobig.gb"
    with pytest.raises(SystemExit, match="max_output_gb"):
        structure_domainate.main(["-i", inputs, "-r", references, "-o", str(out),
                                  "--max_output_gb", "0.000000001", "--cpu", "2"])
    assert not out.exists(), "a size-limited run must not leave a partial output file"


def test_missing_binary(inputs, references, tmp_path, monkeypatch):
    from domainator import structure_lib
    monkeypatch.setattr(structure_lib.shutil, "which", lambda name: None)
    with pytest.raises(RuntimeError, match=r"--foldseek_path"):
        structure_domainate.main(["-i", inputs, "-r", references, "-o", str(tmp_path / "x.gb")])
