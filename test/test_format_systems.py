"""The format_systems command."""
import shutil

import pytest

from domainator import format_systems as format_systems_module
from domainator.errors import InputError
from domainator.format_systems import TSV_COLUMNS, format_systems
from domainator.hmmio import iter_hmm_names
from domainator.systems.native import read_system_models


SOD_MODELS = """\
domainator_systems: 1
database: sod
reference_files: [FeSOD_pfam.hmm, notes.fasta]
defaults: {window: {max_cds: 3}}
systems:
  - {name: FeSOD, condition: "Sod_Fe_N & Sod_Fe_C"}
  - name: CcdB_only
    category: toxin
    condition: "CcdB & Missing"
    window: {max_cds: 1}
"""


@pytest.fixture
def sod_models(shared_datadir, tmp_path):
    shutil.copy(shared_datadir / "FeSOD_pfam.hmm", tmp_path / "FeSOD_pfam.hmm")
    (tmp_path / "notes.fasta").write_text(">x\nMA\n")
    path = tmp_path / "sod.yaml"
    path.write_text(SOD_MODELS)
    return path


def test_command_line_model_to_yaml(tmp_path, shared_datadir):
    out = tmp_path / "xyz.yaml"
    format_systems_module.main(["--domain_expr", "X & Y & Z & (Q | R)", "--window_cds", "10", "--system_name", "XYZ_QR",
                                "--yaml_out", str(out)])
    written, = read_system_models([out])
    assert written.database == "xyz"
    model, = written.models
    assert (model.database, model.name, model.condition) == ("xyz", "XYZ_QR", "X & Y & Z & (Q | R)")
    assert model.window.max_cds == 10


def test_round_trip_selects_models(sod_models, tmp_path):
    out = tmp_path / "one.yaml"
    format_systems([str(sod_models)], systems={"FeSOD"}, yaml_out=str(out))
    original, = read_system_models([sod_models])
    written, = read_system_models([out])
    assert written.models == original.models[:1]


def test_tsv(sod_models, tmp_path):
    out = tmp_path / "models.tsv"
    format_systems([str(sod_models)], tsv_out=str(out))
    header, *rows = [line.split("\t") for line in out.read_text().splitlines()]
    assert tuple(header) == TSV_COLUMNS
    column = {name: i for i, name in enumerate(TSV_COLUMNS)}
    assert [tuple(row[column[name]] for name in ("system", "component", "role", "reference")) for row in rows] == [
        ("FeSOD", "Sod_Fe_N", "core", "Sod_Fe_N"), ("FeSOD", "Sod_Fe_C", "core", "Sod_Fe_C"),
        ("CcdB_only", "CcdB", "core", "CcdB"), ("CcdB_only", "Missing", "core", "Missing")]
    assert rows[2][2] == "toxin" and rows[2][3] == "span:cds=1"


def test_hmm_out(sod_models, tmp_path):
    out = tmp_path / "sod.hmm"
    with pytest.warns(UserWarning) as warned:
        format_systems([str(sod_models)], hmm_out=str(out))
    assert sorted(iter_hmm_names(out)) == ["Sod_Fe_C", "Sod_Fe_N"]
    messages = " ".join(str(w.message) for w in warned)
    assert "Skipping reference file" in messages and "notes.fasta" in messages
    assert "No profile found in the reference files for: CcdB, Missing." in messages


def test_categories(sod_models, tmp_path):
    out = tmp_path / "models.tsv"
    format_systems([str(sod_models)], system_categories={"toxin"}, tsv_out=str(out))
    assert {line.split("\t")[1] for line in out.read_text().splitlines()[1:]} == {"CcdB_only"}
    with pytest.raises(InputError, match="No models match"):
        format_systems([str(sod_models)], system_categories={"PKS"}, tsv_out=str(out))


def test_yaml_out_needs_one_collection(sod_models, tmp_path):
    with pytest.raises(InputError, match="writes one collection"):
        format_systems([str(sod_models)], domain_expr="A & B", window_cds=3, yaml_out=str(tmp_path / "out.yaml"))


def test_unknown_system_name(sod_models, tmp_path):
    with pytest.raises(InputError, match="Nope"):
        format_systems([str(sod_models)], systems={"Nope"}, tsv_out=str(tmp_path / "out.tsv"))


def test_hmm_out_needs_reference_files(tmp_path):
    with pytest.raises(InputError, match="no reference_files"):
        format_systems(domain_expr="A & B", window_cds=3, hmm_out=str(tmp_path / "out.hmm"))


def test_needs_an_output(capsys):
    with pytest.raises(SystemExit):
        format_systems_module.main(["--domain_expr", "A", "--window_cds", "3"])


def test_check_examples(shared_datadir, tmp_path, monkeypatch):
    """--check_examples_out downloads each example (here, a stand-in for NCBI) and looks for its model there."""
    from domainator.seqio import parse_seqfiles
    from domainator.system_db_download import system_db_download
    from domainator.systems import examples

    models = tmp_path / "antismash_rules"
    system_db_download(str(shared_datadir / "systems" / "antismash_rules"), source="antismash", output_dir=str(models))
    fetched = list()

    def fetch(self, example):
        fetched.append((example.accession, example.start, example.end))
        if example.accession == "missing.1":
            raise RuntimeError("no such record")
        return list(parse_seqfiles([str(shared_datadir / "pDONR201.gb")]))

    monkeypatch.setattr(examples.NcbiFetcher, "__init__", lambda self: None)
    monkeypatch.setattr(examples.NcbiFetcher, "__call__", fetch)
    collection, = read_system_models([models])
    from dataclasses import replace
    from domainator.systems.model import Example
    ccd, alone, resistance = collection.models
    collection = replace(collection, models=(ccd, replace(alone, examples=(Example("missing.1"),)),
                                             replace(resistance, examples=(Example("pDONR201.1", 1, 4470),)),
                                             replace(ccd, name="ccd_cat_copy")))
    out = tmp_path / "examples.tsv"
    with pytest.warns(UserWarning, match="0 examples were not found by their models, and 1 couldn't be checked"):
        format_systems([collection], check_examples_out=str(out))
    header, *rows = [line.split("\t") for line in out.read_text().splitlines()]
    assert header == list(examples.REPORT_COLUMNS)
    assert [(row[1], row[2], row[6], row[7], row[9][:12]) for row in rows] == [
        ("ccd_cat", "pDONR201.1", "yes", "1", ""), ("ccdb_alone", "missing.1", "no", "0", "RuntimeError"),
        ("resistance", "pDONR201.1", "yes", "1", ""), ("ccd_cat_copy", "pDONR201.1", "yes", "1", "")]
    assert fetched[0] == ("pDONR201.1", 0, 4470)


@pytest.mark.network
def test_ncbi_fetcher():
    from domainator.systems.examples import NcbiFetcher
    from domainator.systems.model import Example
    record, = NcbiFetcher()(Example("NC_001422.1", 100, 600))  # phage phiX174
    assert len(record) == 501
