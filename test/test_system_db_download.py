"""system_db_download and the padloc-db converter (test/data/systems/padloc_db is a miniature padloc-db)."""
import pytest

from domainator import SYSTEM_FEATURE_NAME
from domainator import domainate as domainate_module
from domainator import system_db_download as system_db_download_module
from domainator.errors import InputError
from domainator.hmmio import iter_hmm_names
from domainator.seqio import parse_seqfiles
from domainator.systems.formats import padloc
from domainator.systems.model import ReferenceRule
from domainator.systems.native import read_system_models


@pytest.fixture
def padloc_db(shared_datadir):
    return shared_datadir / "systems" / "padloc_db"


def models_by_name(collection):
    return {model.name: model for model in collection.models}


def test_read_hmm_meta(padloc_db):
    rows = {row.name: row for row in padloc.read_hmm_meta(padloc_db / "hmm_meta.txt")}
    assert rows["GajB_00001"].proteins == ("GajB", "GajB2")  # split on '|', spaces removed
    assert rows["Cas1_00001"].groups == ("cas_adaptation",)
    # NA thresholds get padloc's defaults.
    assert (rows["Cas1_00001"].evalue, rows["Cas1_00001"].reference_coverage, rows["Cas1_00001"].target_coverage) == (1e-5, 0.3, 0.3)
    assert (rows["Cas2_00001"].evalue, rows["Cas2_00001"].reference_coverage, rows["Cas2_00001"].target_coverage) == (1e-3, 0.2, 0.1)


def test_convert(padloc_db):
    conversion = padloc.convert(padloc_db, release="v0.test")
    collection = conversion.collection
    assert (collection.database, collection.semantics, collection.source.release) == ("padloc", "padloc", "v0.test")
    assert collection.reference_files == ("padloc.hmm", "padloc.cm")
    models = models_by_name(collection)
    assert sorted(models) == ["cas_test", "gabija", "gabija_other", "retron_test"]

    gabija = models["gabija"]
    assert (gabija.description, gabija.accession, gabija.category, gabija.program) == ("Gabija", "gabija", "defense", "padloc")
    assert gabija.window.describe() == "chain:cds=0"
    assert (gabija.min_core, gabija.min_total, gabija.same_strand) == (2, 2, False)
    assert [(c.name, c.role) for c in gabija.components] == [("GajA", "core"), ("GajB", "core"), ("Tnp", "forbidden")]
    assert gabija.components[0].references == (
        ReferenceRule("GajA_00001", database="padloc", evalue=1e-5, min_reference_coverage=0.3, min_target_coverage=0.3),)
    assert gabija.literature == ("doron2018:10/gcwwm5",)

    # The catch-all yields to the canonical models it shares GajA with.
    assert models["gabija_other"].superiors == (("padloc", "gabija"), ("padloc", "retron_test"))

    # The group expands to its proteins; Cas1 is then core, not secondary, and Cas2 is not neutral.
    cas = models["cas_test"]
    assert [(c.name, c.role) for c in cas.components] == [("Cas1", "core"), ("Cas2", "core")]
    assert (cas.min_core, cas.min_total) == (2, 2)  # minimum_total 4 can't be reached without the CRISPR array

    retron = models["retron_test"]
    assert retron.same_strand
    ncrna, = [c for c in retron.components if c.name == "ncRNA"]
    assert (ncrna.role, ncrna.references) == ("accessory", (ReferenceRule("CopA", database="padloc"),))

    assert {(note.model, note.feature, note.action.split(":")[0]) for note in conversion.notes} == {
        ("no_model_file", "model", "dropped"),
        ("gabija_other", "superiors", "added"),
        ("cas_test", "CRISPR_array", "dropped"),
        ("cas_test", "Missing_gene", "dropped"),
        ("cas_test", "minimum_total", "lowered"),
    }


def test_system_db_download_with_input(padloc_db, tmp_path):
    output_dir = tmp_path / "padloc_out"
    system_db_download_module.main(["--source", "padloc", "--input", str(padloc_db), "--output_dir", str(output_dir)])
    assert sorted(path.name for path in output_dir.iterdir()) == ["LICENSE", "conversion.tsv", "padloc.cm", "padloc.hmm", "systems.yaml"]
    # Only the profiles the models use.
    assert sorted(iter_hmm_names(output_dir / "padloc.hmm")) == ["Cas1_00001", "Cas2_00001", "GajA_00001", "GajB_00001", "Tnp_00001"]
    collection, = read_system_models([output_dir])
    assert collection.models == padloc.convert(padloc_db).collection.models
    report = (output_dir / "conversion.tsv").read_text().splitlines()
    assert report[0] == "model\tfeature\taction" and len(report) == 6


def test_converted_models_find_systems(padloc_db, shared_datadir, tmp_path):
    """GajA and GajB are copies of the CcdB and CcdA profiles, so pDONR201's ccdB-ccdA is a Gabija."""
    output_dir = tmp_path / "padloc_out"
    system_db_download_module.system_db_download(str(padloc_db), source="padloc", output_dir=str(output_dir))
    out = tmp_path / "out.gb"
    domainate_module.main(["-i", str(shared_datadir / "pDONR201_multi_genemark.gb"), "--system_models", str(output_dir),
                           "-o", str(out)])
    record = next(iter(parse_seqfiles([str(out)])))
    systems = [f for f in record.features if f.type == SYSTEM_FEATURE_NAME]
    assert [(f.qualifiers["name"][0], f.qualifiers["database"][0]) for f in systems] == [("gabija", "padloc")]
    assert systems[0].qualifiers["component"] == ["core|GajA|1264_-1_959|padloc|GajA_00001",
                                                  "core|GajB|1391_-1_1266|padloc|GajB_00001"]


def test_input_and_release_conflict(padloc_db, tmp_path, capsys):
    with pytest.raises(SystemExit):
        system_db_download_module.main(["--source", "padloc", "--input", str(padloc_db), "--release", "v2.0.0",
                                        "--output_dir", str(tmp_path / "out")])
    assert "can't be combined with --input" in capsys.readouterr().err


def test_not_a_padloc_db(tmp_path):
    with pytest.raises(InputError, match="not a padloc-db directory"):
        padloc.convert(tmp_path)


def test_padloc_db_given_to_system_models_names_the_converter(padloc_db):
    with pytest.raises(InputError, match="system_db_download --source padloc --input"):
        read_system_models([padloc_db])


@pytest.mark.network
def test_download_padloc(tmp_path):
    output_dir = tmp_path / "padloc_db"
    system_db_download_module.main(["--source", "padloc", "--release", "v2.0.0", "--output_dir", str(output_dir)])
    collection, = read_system_models([output_dir])
    assert collection.source.release == "v2.0.0"
    assert len(collection.models) > 300
    assert "gabija" in models_by_name(collection)
