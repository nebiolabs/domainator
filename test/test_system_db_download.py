"""system_db_download and the padloc-db converter (test/data/systems/padloc_db is a miniature padloc-db)."""
import pytest

from domainator import SYSTEM_FEATURE_NAME
from domainator import domainate as domainate_module
from domainator import system_db_download as system_db_download_module
from domainator.errors import InputError
from domainator.hmmio import iter_hmm_names
from domainator.seqio import parse_seqfiles
from domainator.systems.formats import antismash, macsyfinder, padloc
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
    assert [(c.name, c.role) for c in cas.components] == [("Cas1", "core"), ("Cas2", "core"), ("CRISPR_array", "accessory")]
    assert cas.components[2].references == (ReferenceRule("repeat_region", feature=True, qualifiers=(("rpt_family", "CRISPR"),)),)
    assert (cas.min_core, cas.min_total) == (2, 3)  # minimum_total 4 can't be reached without Missing_gene

    retron = models["retron_test"]
    assert retron.same_strand
    ncrna, = [c for c in retron.components if c.name == "ncRNA"]
    assert (ncrna.role, ncrna.references) == ("accessory", (ReferenceRule("CopA", database="padloc"),))

    assert {(note.model, note.feature, note.action.split(":")[0]) for note in conversion.notes} == {
        ("no_model_file", "model", "dropped"),
        ("gabija_other", "superiors", "added"),
        ("cas_test", "CRISPR_array", "feature component"),
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


# MacSyFinder packages. test/data/systems/macsyfinder_pkg is a synthetic package (Pfam profiles,
# renamed): ccd (CcdB, CcdA with an exchangeable, CAT accessory) and kan (APH, gzipped).

@pytest.fixture
def macsyfinder_pkg(shared_datadir):
    return shared_datadir / "systems" / "macsyfinder_pkg"


def copy_package(source, destination):
    import shutil
    shutil.copytree(source, destination)
    return destination


def test_convert_macsyfinder(macsyfinder_pkg):
    conversion = macsyfinder.convert(macsyfinder_pkg, database="toxins")
    collection = conversion.collection
    assert (collection.database, collection.semantics, collection.reference_files) == ("toxins", "macsyfinder", ("toxins.hmm",))
    assert collection.source.release == "0.1.0"  # metadata.yml's vers
    assert collection.source.license == "CC0-1.0; CC0-1.0 (second key, as in defense-finder-models)"
    models = models_by_name(collection)
    ccd, kan = models["ccd"], models["kan"]
    assert (ccd.accession, ccd.category, ccd.program, ccd.window.describe()) == ("TA/ccd", "TA", "macsyfinder", "chain:cds=1")
    assert (ccd.min_core, ccd.min_total, ccd.max_components) == (2, 2, 3)
    assert [(c.name, c.role, c.loner, c.max_gap_cds) for c in ccd.components] == [
        ("CcdB", "core", False, None), ("CcdA", "core", True, None), ("CAT", "accessory", False, 0)]
    # Profiles are named by their NAME line; thresholds from model_conf.xml and each profile's domain GA.
    assert ccd.components[1].references == (
        ReferenceRule("CcdA_renamed", database="toxins", evalue=0.01, min_score=22.4, min_reference_coverage=0.3),
        ReferenceRule("TA__CcdA_2", database="toxins", evalue=0.01, min_score=22.4, min_reference_coverage=0.3,
                      exchangeable=True))
    # Defaults: the number of mandatory genes.
    assert (kan.min_core, kan.min_total, kan.multi_loci) == (1, 1, True)
    assert [(c.name, c.role, c.multi_system, c.multi_model) for c in kan.components] == [
        ("APH", "core", True, True), ("Tnp", "forbidden", False, False)]
    assert {(note.model, note.feature, note.action.split(":")[0]) for note in conversion.notes} == {
        ("kan", "multi_loci", "approximated"), ("ccd", "CcdA", "approximated"), ("ccd", "TA__Missing", "dropped")}
    assert "Kan__APH.hmm.gz" in {path.name for path in conversion.profile_files["toxins.hmm"]}
    # --coverage_profile overrides model_conf.xml.
    ccd = models_by_name(macsyfinder.convert(macsyfinder_pkg, database="toxins", coverage_profile=0.6).collection)["ccd"]
    assert {rule.min_reference_coverage for c in ccd.components for rule in c.references} == {0.6}


def test_macsyfinder_package_quirks(macsyfinder_pkg, tmp_path):
    """Genes listed twice are merged, colliding model names are qualified, and without
    model_conf.xml the source's default coverage applies."""
    package = copy_package(macsyfinder_pkg, tmp_path / "toxins-2.1")
    (package / "model_conf.xml").unlink()
    (package / "definitions" / "Kan" / "ccd.xml").write_text(
        '<model inter_gene_max_space="0" vers="2.0">\n<gene name="Kan__APH" presence="mandatory">'
        '<exchangeables><gene name="Kan__APH"/></exchangeables></gene>\n'
        '<gene name="Kan__Tnp" presence="forbidden"/>\n<gene name="Kan__Tnp" presence="forbidden"/>\n</model>\n')
    assert macsyfinder.database_name(package) == "toxins"
    conversion = macsyfinder.convert(package, database="toxins", default_coverage=0.4)
    models = models_by_name(conversion.collection)
    assert sorted(models) == ["Kan/ccd", "TA/ccd", "kan"]
    assert [len(c.references) for c in models["Kan/ccd"].components] == [1, 1]
    assert models["kan"].components[0].references[0].min_reference_coverage == 0.4
    assert {(note.model, note.feature, note.action.split(":")[0]) for note in conversion.notes} >= {
        ("Kan/ccd", "name", "qualified"), ("TA/ccd", "name", "qualified"),
        ("Kan/ccd", "Kan__APH", "merged"), ("Kan/ccd", "Kan__Tnp", "merged")}


def test_macsyfinder_version_1_refused(macsyfinder_pkg, tmp_path):
    from domainator.errors import FormatVersionError
    package = copy_package(macsyfinder_pkg, tmp_path / "old")
    (package / "definitions" / "TA" / "ccd.xml").write_text('<system inter_gene_max_space="1"><gene name="TA__CcdB" presence="mandatory"/></system>')
    with pytest.raises(FormatVersionError, match="format version 1.0"):
        macsyfinder.convert(package, database="old")


@pytest.mark.parametrize("source, database, coverage", [("macsyfinder", "macsyfinder_pkg", 0.3),
                                                        ("defensefinder", "defensefinder", 0.3)])
def test_system_db_download_macsyfinder(macsyfinder_pkg, shared_datadir, tmp_path, source, database, coverage):
    output_dir = tmp_path / "out"
    system_db_download_module.main(["--source", source, "--input", str(macsyfinder_pkg), "--output_dir", str(output_dir)])
    assert sorted(path.name for path in output_dir.iterdir()) == ["LICENSE", "conversion.tsv", f"{database}.hmm", "systems.yaml"]
    assert sorted(iter_hmm_names(output_dir / f"{database}.hmm")) == ["CcdA_renamed", "Kan__APH", "Kan__Tnp", "TA__CAT", "TA__CcdA_2", "TA__CcdB"]
    collection, = read_system_models([output_dir])
    assert collection.database == database
    assert {rule.min_reference_coverage for m in collection.models for c in m.components for rule in c.references} == {coverage}

    out = tmp_path / "out.gb"
    domainate_module.main(["-i", str(shared_datadir / "pDONR201_multi_genemark.gb"), "--system_models", str(output_dir),
                           "-o", str(out)])
    record = next(iter(parse_seqfiles([str(out)])))
    systems = {f.qualifiers["name"][0]: f for f in record.features if f.type == SYSTEM_FEATURE_NAME}
    assert sorted(systems) == ["ccd", "kan"]
    assert [value.split("|")[1] for value in systems["ccd"].qualifiers["component"]] == ["CcdB", "CcdA", "CAT"]


def test_system_db_download_errors(macsyfinder_pkg, padloc_db, tmp_path, capsys):
    with pytest.raises(SystemExit):
        system_db_download_module.main(["--source", "macsyfinder", "--output_dir", str(tmp_path / "out")])
    assert "needs --input" in capsys.readouterr().err
    with pytest.raises(SystemExit):
        system_db_download_module.main(["--source", "padloc", "--input", str(padloc_db), "--coverage_profile", "0.5",
                                        "--output_dir", str(tmp_path / "out")])
    assert "--coverage_profile is for MacSyFinder" in capsys.readouterr().err
    with pytest.raises(InputError, match="system_db_download --source defensefinder"):
        read_system_models([macsyfinder_pkg])


@pytest.mark.network
def test_download_defensefinder(tmp_path):
    output_dir = tmp_path / "defensefinder"
    system_db_download_module.main(["--source", "defensefinder", "--release", "3.1.0", "--output_dir", str(output_dir)])
    collection, = read_system_models([output_dir])
    assert (collection.database, collection.source.release) == ("defensefinder", "3.1.0")
    assert len(collection.models) > 600
    gabija = models_by_name(collection)["Gabija"]
    assert gabija.accession == "DefenseFinder/Gabija/Gabija"
    assert gabija.components[0].references[0].min_reference_coverage == 0.4


@pytest.fixture
def antismash_rules(shared_datadir):
    """A synthetic antiSMASH checkout: rules written for these tests, in antiSMASH's language, with Pfam profiles."""
    return shared_datadir / "systems" / "antismash_rules"


def test_convert_antismash(antismash_rules):
    conversion = antismash.convert(antismash_rules, strictness="loose")
    models = models_by_name(conversion.collection)
    assert list(models) == ["ccd_cat", "ccdb_alone", "resistance", "needs_code"]
    ccd = models["ccd_cat"]
    assert ccd.condition == "ANY_TOX & CAT & ~cds(CcdB & APH)"
    assert ccd.define == (("ANY_TOX", "(CcdB | TCAD9)"),)
    assert ccd.extenders == "APH"
    assert (ccd.window.mode, ccd.window.max_kb, ccd.neighborhood_kb) == ("anchor", 0.999, 2.0)
    assert (ccd.category, ccd.description, ccd.accession) == ("toxin", "A toxin near a chloramphenicol acetyltransferase", "antismash:ccd_cat")
    assert [(e.accession, e.start, e.end, e.note) for e in ccd.examples] == [("pDONR201.1", 0, 4470, "synthetic")]
    assert ccd.related == ("CcdA", "Condensation")
    assert [(c.name, c.role) for c in ccd.components] == [("CcdB", "core"), ("TCAD9", "core"), ("CAT", "core"), ("APH", "neutral")]
    assert ccd.components[2].references == (ReferenceRule("CAT", database="antismash", min_score=200.0),)
    assert models["resistance"].condition == "minimum(2, [CAT, APH, Condensation]) | repeatable(3, [CAT, APH])"
    assert models["ccdb_alone"].superiors == (("antismash", "ccd_cat"),)
    assert models["needs_code"].superiors == ()  # its superior isn't a rule
    assert [(c.name, c.role, c.references) for c in models["needs_code"].components][1:] == [
        ("fake_precursor", "core", ()), ("Unknown_profile", "forbidden", ())]
    assert set(conversion.collection.competing) == {frozenset({"CcdB", "TCAD9"}), frozenset({"CAT", "Condensation"})}
    assert {(note.model, note.feature, note.action.split(":")[0]) for note in conversion.notes} == {
        ("needs_code", "fake_precursor", "no references"), ("needs_code", "Unknown_profile", "no references"),
        ("needs_code", "superiors", "dropped"), (".", "license", "AGPL-3.0")}
    assert [len(antismash.convert(antismash_rules, strictness=level).collection.models) for level in antismash.STRICTNESS] == [2, 3, 4]


def test_antismash_substitutes(antismash_rules, tmp_path):
    substitutes = tmp_path / "substitutes.yaml"
    substitutes.write_text("domainator_substitutes: 1\ncomponents:\n"
                           "  - {name: fake_precursor, domains: [{name: W_W_K, database: MOTIF, max_target_length: 80}]}\n"
                           "  - {name: unused, domains: [X]}\n")
    conversion = antismash.convert(antismash_rules, strictness="loose", substitutes_file=str(substitutes))
    precursor = models_by_name(conversion.collection)["needs_code"].components[1]
    assert precursor.references == (ReferenceRule("W_W_K", database="MOTIF", max_target_length=80),)
    assert {(note.feature, note.action.split(":")[0]) for note in conversion.notes if note.model == "."} == {
        ("fake_precursor", "substituted"), ("unused", "unused substitute"), ("license", "AGPL-3.0")}
    substitutes.write_text("components: []\n")
    with pytest.raises(InputError, match="not a substitutes file"):
        antismash.convert(antismash_rules, substitutes_file=str(substitutes))


def test_system_db_download_antismash(antismash_rules, shared_datadir, tmp_path, capsys):
    output_dir = tmp_path / "out"
    system_db_download_module.main(["--source", "antismash", "--input", str(antismash_rules / "antismash" / "detection" / "hmm_detection"),
                                    "--output_dir", str(output_dir)])
    assert "for your own use" in capsys.readouterr().err
    assert sorted(path.name for path in output_dir.iterdir()) == ["LICENSE", "NOTICE", "antismash.hmm", "conversion.tsv", "systems.yaml"]
    assert "Synthetic test rules" in (output_dir / "LICENSE").read_text()  # found above the hmm_detection directory
    assert "AGPL-3.0" in (output_dir / "NOTICE").read_text()
    assert sorted(iter_hmm_names(output_dir / "antismash.hmm")) == ["APH", "CAT", "CcdB", "Condensation", "TCAD9"]
    collection, = read_system_models([output_dir])
    assert [model.name for model in collection.models] == ["ccd_cat", "ccdb_alone", "resistance"]

    out = tmp_path / "out.gb"
    domainate_module.main(["-i", str(shared_datadir / "pDONR201.gb"), "--system_models", str(output_dir), "-o", str(out)])
    record = next(iter(parse_seqfiles([str(out)])))
    systems = {f.qualifiers["name"][0]: f for f in record.features if f.type == SYSTEM_FEATURE_NAME}
    # ccdb_alone yields to ccd_cat, its superior; ccd_cat grows over APH, its extender
    assert sorted(systems) == ["ccd_cat", "resistance"]
    assert systems["ccd_cat"].qualifiers["component"] == [
        "core|CcdB|1264_-1_959|antismash|CcdB", "core|CAT|2265_-1_1606|antismash|CAT", "extender|APH|2868_1_3677|antismash|APH"]
    assert [value.split("|")[1] for value in systems["resistance"].qualifiers["component"]] == ["CAT", "APH"]
    assert systems["ccd_cat"].qualifiers["neighborhood_kb"] == ["2"]


def test_antismash_errors(antismash_rules, tmp_path, capsys):
    with pytest.raises(SystemExit):
        system_db_download_module.main(["--source", "antismash", "--output_dir", str(tmp_path / "out")])
    assert "needs --input" in capsys.readouterr().err
    with pytest.raises(InputError, match="system_db_download --source antismash"):
        read_system_models([antismash_rules / "antismash" / "detection" / "hmm_detection"])
    with pytest.raises(InputError, match="not an antiSMASH hmm_detection directory"):
        antismash.convert(tmp_path)
    rules = tmp_path / "rules"
    (rules / "cluster_rules").mkdir(parents=True)
    (rules / "data").mkdir()
    (rules / "data" / "hmmdetails.txt").write_text("")
    for text, message in (("CATEGORY x\n", "outside a RULE"), ("RULE a\nCATEGORY x\nCUTOFF 1\nNEIGHBOURHOOD 1\n", "has no CONDITIONS"),
                          ("DEFINE X Y\n", "DEFINE <label> AS"), ("RULE a\nCATEGORY x\nCUTOFF one\nNEIGHBOURHOOD 1\nCONDITIONS a\n", "must be numbers")):
        (rules / "cluster_rules" / "strict.txt").write_text(text)
        with pytest.raises(InputError, match=message):
            antismash.convert(rules, strictness="strict")


@pytest.mark.network
def test_convert_antismash_release(tmp_path):
    """Every rule of a pinned antiSMASH release converts, and the result reads back unchanged."""
    from domainator.systems.formats import download_github_release, write_conversion
    root, release, _ = download_github_release("antismash/antismash", "8-0-4", tmp_path, label="antismash")
    conversion = antismash.convert(root, strictness="loose")
    assert len(conversion.collection.models) >= 100
    assert conversion.collection.source.release == "8.0.4"
    assert sum(1 for note in conversion.notes if "dynamic profile" in note.action) == 6
    write_conversion(conversion, tmp_path / "out")
    collection, = read_system_models([tmp_path / "out"])
    assert collection.models == conversion.collection.models
