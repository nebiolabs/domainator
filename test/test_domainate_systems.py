"""domainate --system_models: search the models' profiles and detect systems in one pass."""
import shutil

import pytest

from domainator import DOMAIN_FEATURE_NAME, SYSTEM_FEATURE_NAME
from domainator import domainate as domainate_module
from domainator.errors import InputError
from domainator.find_systems import find_systems
from domainator.seqio import parse_seqfiles

MODELS = """\
domainator_systems: 1
database: ccd
reference_files: [pdonr_hmms_2.hmm]
systems:
  - name: CcdAB
    category: toxin_antitoxin
    condition: "CcdA & CcdB"
    window: {max_cds: 2}
  - name: CcdB_Condensation
    window: {max_cds: 5}
    components:
      - {name: toxin, role: core, domains: [CcdB]}
      - {name: condensation, role: core, domains: [{name: Condensation, evalue: CONDENSATION_EVALUE}]}
"""


@pytest.fixture
def models(shared_datadir, tmp_path):
    def write(condensation_evalue="1.0e-3"):
        shutil.copy(shared_datadir / "pdonr_hmms_2.hmm", tmp_path / "pdonr_hmms_2.hmm")
        path = tmp_path / "ccd.yaml"
        path.write_text(MODELS.replace("CONDENSATION_EVALUE", condensation_evalue))
        return path
    return write


def run(shared_datadir, tmp_path, *args):
    out = tmp_path / "out.gb"
    domainate_module.main(["-i", str(shared_datadir / "pDONR201_multi_genemark.gb"), "-o", str(out), "-Z", "0", *args])
    return list(parse_seqfiles([str(out)]))


def systems_of(record):
    return sorted(f.qualifiers["name"][0] for f in record.features if f.type == SYSTEM_FEATURE_NAME)


def hits(record, cds_id, database=None):
    return sorted((f.qualifiers["name"][0], f.qualifiers["evalue"][0]) for f in record.features
                  if f.type == DOMAIN_FEATURE_NAME and f.qualifiers["cds_id"] == [cds_id]
                  and (database is None or f.qualifiers["database"] == [database]))


def test_models_alone(shared_datadir, tmp_path, models):
    records = run(shared_datadir, tmp_path, "--system_models", str(models()), "-e", "1e-3")
    assert [systems_of(record) for record in records] == [["CcdAB", "CcdB_Condensation"]] * len(records)
    system = next(f for f in records[0].features if f.type == SYSTEM_FEATURE_NAME and f.qualifiers["name"] == ["CcdAB"])
    assert system.qualifiers["component"] == ["core|CcdB|1264_-1_959|pdonr_hmms_2|CcdB",
                                              "core|CcdA|1391_-1_1266|pdonr_hmms_2|CcdA"]
    assert {f.qualifiers["database"][0] for f in records[0].features if f.type == DOMAIN_FEATURE_NAME} == {"pdonr_hmms_2"}


def test_model_profiles_are_exempt_from_hit_caps(shared_datadir, tmp_path, models):
    records = run(shared_datadir, tmp_path, "--system_models", str(models()), "-r", str(shared_datadir / "pdonr_hmms_1.hmm"),
                  "-e", "10", "--max_hits_per_cds", "1")
    # The -r profiles are capped at one hit per CDS; the model profiles are not.
    assert hits(records[0], "2265_-1_1606", "pdonr_hmms_1") == [("CAT", "9.6e-102")]
    assert hits(records[0], "2265_-1_1606", "pdonr_hmms_2") == [("Condensation", "1.4e-04"), ("Condensation", "7.9e-01")]


def test_model_profiles_are_searched_with_their_own_looser_threshold(shared_datadir, tmp_path, models):
    tight = run(shared_datadir, tmp_path, "--system_models", str(models("1.0e-6")), "-e", "1e-5")
    assert hits(tight[0], "2265_-1_1606") == []
    loose = run(shared_datadir, tmp_path, "--system_models", str(models("1.0")), "-e", "1e-5")
    assert hits(loose[0], "2265_-1_1606") == [("Condensation", "1.4e-04"), ("Condensation", "7.9e-01")]
    # -e still filters the evidence, as in find_systems, so the looser hits make no system.
    assert systems_of(loose[0]) == ["CcdAB"]
    # And the two routes agree.
    again, = list(find_systems([loose[0]], system_models=[str(models("1.0"))], evalue=1e-5))
    assert systems_of(again) == ["CcdAB"]


def test_a_reference_that_is_a_model_file_is_searched_once(shared_datadir, tmp_path, models):
    path = models()
    records = run(shared_datadir, tmp_path, "--system_models", str(path), "-r", str(path.parent / "pdonr_hmms_2.hmm"),
                  "-e", "1e-3")
    assert hits(records[0], "1264_-1_959") == [("CcdB", "1.1e-32")]


def test_systems_restricts_models(shared_datadir, tmp_path, models):
    records = run(shared_datadir, tmp_path, "--system_models", str(models()), "-e", "1e-3", "--systems", "CcdAB")
    assert systems_of(records[0]) == ["CcdAB"]


def test_models_need_something_to_search(shared_datadir, tmp_path, capsys):
    path = tmp_path / "bare.yaml"
    path.write_text("domainator_systems: 1\nsystems:\n  - {name: AB, condition: 'A & B', window: {max_cds: 2}}\n")
    with pytest.raises(SystemExit):
        run(shared_datadir, tmp_path, "--system_models", str(path))
    assert "list no reference_files" in capsys.readouterr().err
    records = parse_seqfiles([str(shared_datadir / "pDONR201_multi_genemark.gb")])
    with pytest.raises(InputError, match="Nothing to search"):
        list(domainate_module.domainate(records, references=None, z=1000, system_models=[str(path)]))


def test_library_with_collections_and_references(shared_datadir, models):
    from domainator.systems.native import read_system_models
    collections = read_system_models([models()])
    records = list(domainate_module.domainate(parse_seqfiles([str(shared_datadir / "pDONR201_multi_genemark.gb")]),
                                              references=[str(shared_datadir / "pdonr_hmms_1.hmm")], z=1000,
                                              system_models=collections, max_systems_per_contig=1))
    assert [systems_of(record) for record in records] == [["CcdAB"]] * len(records)
    assert {f.qualifiers["database"][0] for f in records[0].features if f.type == DOMAIN_FEATURE_NAME} == {
        "pdonr_hmms_1", "pdonr_hmms_2"}
