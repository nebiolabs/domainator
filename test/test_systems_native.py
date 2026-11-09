"""Reading and writing native system model files (domainator.systems.native)."""
import os
import tempfile

import pytest
import yaml

from domainator.errors import FormatVersionError, InputError
from domainator.systems.detect import detect_systems
from domainator.systems.model import ReferenceRule
from domainator.systems.native import collection_from_data, parse_version, read_system_models, window_model, write_system_models

from system_records import contig, member_ids


def read_text(text, name="models.yaml"):
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, name)
        with open(path, "w") as handle:
            handle.write(text)
        return read_system_models([path])[0]


def test_short_form_is_the_full_form(shared_datadir):
    full, = read_system_models([shared_datadir / "systems" / "xyz_full.yaml"])
    short, = read_system_models([shared_datadir / "systems" / "xyz_short.yaml"])
    assert full == short
    assert full.models[0].components[0].references == (ReferenceRule("X", evalue=1e-5),)


def test_three_spellings_find_the_same_systems(shared_datadir):
    full, = read_system_models([shared_datadir / "systems" / "xyz_full.yaml"])
    short, = read_system_models([shared_datadir / "systems" / "xyz_short.yaml"])
    command_line = window_model("X & Y & Z & (Q | R)", window_cds=10)
    record = contig(["", "X", "", "Y", "Z", "Q", "", "", "", "", "", "", "", "X", "Y", "", "R", "Z"])
    calls = [[(str(hit.location), member_ids(hit), hit.score) for hit in detect_systems(record, [collection])]
             for collection in (full, short, command_line)]
    assert calls[0] == calls[1] == calls[2]
    assert [ids for _, ids, _ in calls[0]] == [["cds1", "cds3", "cds4", "cds5"], ["cds13", "cds14", "cds16", "cds17"]]


def test_command_line_model():
    collection = window_model("HicB < RelE & ~Tnp", window_kb=2.5, window_mode="chain", system_name="TA", same_strand=True)
    model, = collection.models
    assert (collection.database, model.database, model.name, model.description) == ("command_line", "command_line", "TA", "HicB < RelE & ~Tnp")
    assert {c.name: c.role for c in model.components} == {"HicB": "core", "RelE": "core", "Tnp": "forbidden"}
    assert (model.window.mode, model.window.max_kb, model.window.max_cds, model.same_strand) == ("chain", 2.5, None, True)
    with pytest.raises(InputError, match="needs a window"):
        window_model("A & B")


def test_database_prefix_limits_a_name():
    model, = window_model("pfam:GajA & GajB", window_cds=3).models
    assert model.components[0].references == (ReferenceRule("GajA", database="pfam"),)
    assert model.components[1].references == (ReferenceRule("GajB"),)


@pytest.mark.parametrize("fixture", ["xyz_full.yaml", "xyz_short.yaml", "defense.yaml"])
def test_round_trip(shared_datadir, fixture):
    collection, = read_system_models([shared_datadir / "systems" / fixture])
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "out.yaml")
        write_system_models(collection, path)
        again, = read_system_models([path])
        with open(path) as handle:
            data = yaml.safe_load(handle)
    assert again == collection
    if fixture != "defense.yaml":
        assert data["systems"] == {"XYZ_QR": "X & Y & Z & (Q | R)"}  # the writer chooses the short form


def test_round_trip_everything():
    text = """
domainator_systems: 1
database: everything
semantics: padloc
source: {name: padloc, url: "https://example.org", release: v2.0.0, date: 2026-11-09, converter: domainator 1.1.0}
reference_files: [profiles/defense.hmm]
competing:
  - [B_1, B_2]
competing_overlap_aa: 30
systems:
  - name: First
    category: defense
    description: a model
    accession: everything/First
    window: {mode: chain, max_cds: 3, max_span_kb: 40}
    min_core: 1
    min_total: 2
    components:
      - {name: A, role: core, domains: [A, {name: A2, database: other, evalue: 1e-5, min_score: 20, min_reference_coverage: 0.4, min_target_coverage: 0.3, min_target_length: 10, max_target_length: 80, exchangeable: true}], max_gap_cds: 1}
      - {name: B, role: accessory, domains: [B_1, B_2], loner: true, multi_system: true, multi_model: true}
      - {name: C, role: neutral, domains: [{name: nc, attach_to_cds: true}]}
    same_strand: true
    neighborhood_kb: 5
    neighborhood_cds: 3
    max_components: 4
    multi_loci: true
    orientation_from: B
    examples: [{accession: NZ_CP000000.1, start: 1, end: 900, note: an example}]
    literature: ["10.1000/x"]
    reference_files: [extra.hmm]
  - name: Second
    window: {mode: anchor, max_kb: 19.999}
    condition: "A & ~B"
    superiors: [First]
"""
    collection = read_text(text)
    assert collection.models[1].superiors == (("everything", "First"),)
    assert collection.source.date == "2026-11-09"
    assert collection.models[0].components[0].references[1].evalue == 1e-5  # 1e-5 is text to YAML 1.1
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "out.yaml")
        write_system_models(collection, path)
        again, = read_system_models([path])
    assert again == collection


@pytest.mark.parametrize("value,version", [(1, (1, 0)), (1.0, (1, 0)), ("1.1", (1, 1)), ("1.10", (1, 10))])
def test_parse_version(value, version):
    assert parse_version(value) == version


@pytest.mark.parametrize("version", ["2", "\"1.1\""])
def test_newer_version_refused(version):
    with pytest.raises(FormatVersionError, match="format"):
        read_text(f"domainator_systems: {version}\nsystems: {{A: \"A\"}}\ndefaults: {{window: {{max_cds: 1}}}}\n")


@pytest.mark.parametrize("text,message", [
    ("systems: {A: A}", "no top-level 'domainator_systems:'"),
    ("domainator_systems: 1\ndefaults: {window: {max_cds: 1}}\nsystem: {A: A}\n", r"line 1: unknown key 'system'"),
    ("domainator_systems: 1\nsystems:\n  - name: A\n    window: {max_cds: 1}\n    condition: A\n    extenders: B\n", "format 1.1"),
    ("domainator_systems: 1\ndefine: {X: A}\nsystems: {A: A}\n", "format 1.1"),
    ("domainator_systems: 1\nsystems: {A: A}\n", "no window"),
    ("domainator_systems: 1\ndefaults: {window: {mode: span}}\nsystems: {A: A}\n", "needs max_cds, max_kb or both"),
    ("domainator_systems: 1\ndefaults: {window: {max_cds: 2, max_span_cds: 3}}\nsystems: {A: A}\n", "in span mode"),
    ("domainator_systems: 1\ndefaults: {window: {max_cds: 1}}\nsystems:\n  - name: A\n    components: [{name: X, domains: [X]}]\n    condition: X & Y\n", "Y, which are not components"),
    ("domainator_systems: 1\ndefaults: {window: {max_cds: 1}}\nsystems:\n  - name: A\n    components: [{name: X, role: accessory, domains: [X]}]\n", "no requirement"),
    ("domainator_systems: 1\ndefaults: {window: {max_cds: 1}}\nsystems: {A: \"~B\"}\n", "no positive term"),
    ("domainator_systems: 1\ndefaults: {window: {max_cds: 1}}\nsystems: {\"A B\": A}\n", "can't contain"),
    ("domainator_systems: 1\ndefaults: {window: {max_cds: 1}}\nsystems:\n  - {name: A, condition: A}\n  - {name: A, condition: B}\n", "duplicate model names"),
    ("domainator_systems: 1\ndefaults: {window: {max_cds: 1}}\nsystems:\n  - {name: A, condition: A, superiors: [B]}\n", "superior 'B'"),
    ("domainator_systems: 1\ndefaults: {window: {max_cds: 1}}\nsystems:\n  - {name: A, condition: A, superiors: [\"other:B\"]}\n", "superior other:B"),
    ("domainator_systems: 1\ndefaults: {window: {max_cds: one}}\nsystems: {A: A}\n", "must be a number"),
    ("domainator_systems: 1\ndefaults: {window: {max_cds: 1}}\nsystems: {yes: A}\n", "must be text"),
])
def test_malformed(text, message):
    with pytest.raises(InputError, match=message):
        read_text(text)


def test_duplicate_across_files(shared_datadir):
    path = shared_datadir / "systems" / "xyz_full.yaml"
    with pytest.raises(InputError, match="defined twice"):
        read_system_models([path, shared_datadir / "systems" / "xyz_short.yaml"])


def test_directory(shared_datadir):
    with tempfile.TemporaryDirectory() as directory:
        collection, = read_system_models([shared_datadir / "systems" / "xyz_short.yaml"])
        write_system_models(collection, os.path.join(directory, "systems.yaml"))
        assert read_system_models([directory])[0] == collection
        os.remove(os.path.join(directory, "systems.yaml"))
        with pytest.raises(InputError, match="without a systems.yaml"):
            read_system_models([directory])


@pytest.mark.parametrize("files,label", [
    ({"sys/Gabija.yaml": "maximum_separation: 3\n", "hmm_meta.txt": ""}, "padloc-db directory"),
    ({"metadata.yml": "maintainer: x\n", "definitions/Gabija.xml": "<model/>"}, "MacSyFinder model package"),
    ({"cluster_rules/relaxed.txt": "RULE x\n"}, "antiSMASH"),
])
def test_foreign_directory_refused(files, label):
    with tempfile.TemporaryDirectory() as directory:
        for name, content in files.items():
            os.makedirs(os.path.dirname(os.path.join(directory, name)), exist_ok=True)
            with open(os.path.join(directory, name), "w") as handle:
                handle.write(content)
        with pytest.raises(InputError, match=f"looks like a {label}.*must be converted"):
            read_system_models([directory])


@pytest.mark.parametrize("name,text,label", [
    ("Gabija.yaml", "maximum_separation: 3\ncore_genes: [GajA]\n", "PADLOC system model"),
    ("Gabija.xml", "<model inter_gene_max_space='3'/>", "MacSyFinder model definition"),
    ("relaxed.txt", "# rules\nRULE T1PKS\n  CUTOFF 20\n", "antiSMASH cluster rule file"),
])
def test_foreign_file_refused(name, text, label):
    with pytest.raises(InputError, match=f"looks like a {label}"):
        read_text(text, name)


def test_collection_from_data_default_database():
    collection = collection_from_data({"domainator_systems": 1, "defaults": {"window": {"max_cds": 1}}, "systems": {"A": "A"}},
                                      default_database="mine")
    assert collection.database == collection.models[0].database == "mine"
