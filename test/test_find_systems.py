"""The find_systems command."""
import os
import tempfile

import pytest

from domainator import SYSTEM_FEATURE_NAME
from domainator import find_systems as find_systems_module
from domainator.errors import InputError
from domainator.find_systems import find_systems
from domainator.select_by_cds import select_by_cds
from domainator.seqio import parse_seqfiles
from domainator.systems.features import component_values

from system_records import contig


def systems_of(record):
    return [(f.qualifiers["database"][0], f.qualifiers["name"][0], str(f.location)) for f in record.features if f.type == SYSTEM_FEATURE_NAME]


def run(shared_datadir, *args):
    with tempfile.TemporaryDirectory() as directory:
        out = os.path.join(directory, "out.gb")
        find_systems_module.main(["-i", str(shared_datadir / "systems" / "contigs.gb"), "-o", out, *args])
        return {record.id: record for record in parse_seqfiles([out])}


def test_models_from_files(shared_datadir):
    records = run(shared_datadir, "--system_models", str(shared_datadir / "systems" / "defense.yaml"),
                  str(shared_datadir / "systems" / "xyz_short.yaml"))
    assert {rid: systems_of(r) for rid, r in records.items()} == {
        "xyzq": [("xyz", "XYZ_QR", "[1000:5900](+)")],
        "xyz_only": [],
        "gabija_ta": [("defense", "Gabija", "[0:2900](+)"), ("defense", "TA", "[6000:7900](+)")],
        "vetoed": [],
        "circular": [("defense", "Gabija", "join{[6000:7000](+), [0:900](+)}")],
    }
    gabija = next(f for f in records["gabija_ta"].features if f.type == SYSTEM_FEATURE_NAME)
    assert gabija.qualifiers["category"] == ["defense"]
    assert gabija.qualifiers["window"] == ["chain:cds=2"]
    assert gabija.qualifiers["component"] == ["core|GajA|cds0|db|GajA", "core|GajB|cds2|db|GajB_1"]


def test_command_line_model(shared_datadir):
    records = run(shared_datadir, "--domain_expr", "X & Y & Z & (Q | R)", "--window_cds", "10", "--system_name", "XYZ_QR")
    feature, = [f for f in records["xyzq"].features if f.type == SYSTEM_FEATURE_NAME]
    assert (feature.qualifiers["database"], feature.qualifiers["name"], feature.qualifiers["description"]) == (
        ["command_line"], ["XYZ_QR"], ["X & Y & Z & (Q | R)"])
    assert feature.qualifiers["completeness"] == ["."]


def test_stdin_stdout(shared_datadir, monkeypatch, capsys):
    with open(shared_datadir / "systems" / "contigs.gb") as handle:
        monkeypatch.setattr("sys.stdin", handle)
        find_systems_module.main(["-i", "-", "--domain_expr", "GajA & GajB_1", "--window_cds", "3", "-o", "-"])
    out = capsys.readouterr().out
    assert out.count("Domain_System   ") == 2


def test_evidence_filter(shared_datadir):
    defense = str(shared_datadir / "systems" / "defense.yaml")
    assert systems_of(run(shared_datadir, "--system_models", defense, "-e", "1e-10")["gabija_ta"]) != []
    assert systems_of(run(shared_datadir, "--system_models", defense, "-e", "1e-11")["gabija_ta"]) == []
    assert systems_of(run(shared_datadir, "--system_models", defense, "--databases", "other")["gabija_ta"]) == []


def test_systems_option(shared_datadir):
    defense = str(shared_datadir / "systems" / "defense.yaml")
    assert [name for _, name, _ in systems_of(run(shared_datadir, "--system_models", defense, "--systems", "defense:TA")["gabija_ta"])] == ["TA"]
    with pytest.raises(InputError, match="Lamassu"):
        list(find_systems([contig(["A"])], system_models=[defense], systems={"Lamassu"}))


def test_score_filters():
    assert len(systems_of(next(find_systems([contig(["A", "B"])], domain_expr="A & B", window_cds=2, min_score=100)))) == 1
    assert systems_of(next(find_systems([contig(["A", "B"])], domain_expr="A & B", window_cds=2, min_score=100.1))) == []
    # command-line models have no completeness, so min_completeness keeps them
    assert len(systems_of(next(find_systems([contig(["A", "B"])], domain_expr="A & B", window_cds=2, min_completeness=1)))) == 1


def test_clear(shared_datadir):
    defense = shared_datadir / "systems" / "defense.yaml"
    record = contig(["GajA", "GajB_1"])
    record = next(find_systems([record], system_models=[defense]))
    record = next(find_systems([record], domain_expr="GajA & GajB_1", window_cds=2, system_name="first"))
    record = next(find_systems([record], domain_expr="GajA & GajB_1", window_cds=2, system_name="second"))
    assert [name for _, name, _ in systems_of(record)] == ["Gabija", "first", "second"]
    record = next(find_systems([record], domain_expr="GajA & GajB_1", window_cds=2, system_name="first"))
    record = next(find_systems([record], system_models=[defense]))
    assert sorted(name for _, name, _ in systems_of(record)) == ["Gabija", "first", "second"]
    record = next(find_systems([record], system_models=[defense], clear=False))
    assert sorted(name for _, name, _ in systems_of(record)) == ["Gabija", "Gabija", "first", "second"]


def test_needs_models():
    with pytest.raises(InputError, match="No system models"):
        list(find_systems([contig(["A"])]))
    with pytest.raises(SystemExit):
        find_systems_module.main(["-i", "x.gb", "-o", "y.gb"])
    with pytest.raises(SystemExit):
        find_systems_module.main(["-i", "x.gb", "-o", "y.gb", "--system_models", "m.yaml", "--window_cds", "3"])


def test_genbank_round_trip_with_pipe_in_a_reference_name():
    record = contig(["A", "B"])
    record.features[1].qualifiers["name"] = ["A|v2"]
    with tempfile.TemporaryDirectory() as directory:
        models = os.path.join(directory, "models.yaml")
        with open(models, "w") as handle:
            handle.write("domainator_systems: 1\nsystems:\n  - name: S\n    window: {max_cds: 2}\n"
                         "    components: [{name: A, domains: [\"A|v2\"]}, {name: B, domains: [B]}]\n")
        out = os.path.join(directory, "out.gb")
        from domainator.seqio import write_genbank
        write_genbank(find_systems([record], system_models=[models]), out)
        again, = parse_seqfiles([out])
    feature, = [f for f in again.features if f.type == SYSTEM_FEATURE_NAME]
    assert [tuple(value) for value in component_values(feature)] == [("core", "A", "cds0", "db", "A|v2"), ("core", "B", "cds1", "db", "B")]


def test_cds_ids_survive_slicing(shared_datadir):
    records = list(find_systems(parse_seqfiles([shared_datadir / "systems" / "contigs.gb"]),
                                system_models=[shared_datadir / "systems" / "xyz_short.yaml"]))
    xyzq = [r for r in records if r.id == "xyzq"]
    region, = select_by_cds(xyzq, domains={"Y"}, cds_range=(2, 2))
    feature, = [f for f in region.features if f.type == SYSTEM_FEATURE_NAME]
    assert str(feature.location) == "[0:4900](+)"
    region_cds_ids = {f.qualifiers["cds_id"][0] for f in region.features if f.type == "CDS"}
    assert {value.cds_id for value in component_values(feature)} <= region_cds_ids
    cut, = select_by_cds(xyzq, domains={"Y"}, cds_range=(1, 1))
    assert [f for f in cut.features if f.type == SYSTEM_FEATURE_NAME] == []


def test_truncated_systems(shared_datadir):
    expr = ("--domain_expr", "X & Y & Z & (Q | R)", "--window_cds", "10")
    assert systems_of(run(shared_datadir, *expr)["xyz_only"]) == []
    records = run(shared_datadir, *expr, "--truncated_systems", "include")
    feature, = [f for f in records["xyz_only"].features if f.type == SYSTEM_FEATURE_NAME]
    # Q or R could lie beyond either end of this 5-CDS contig
    assert (str(feature.location), feature.qualifiers["system_partial"]) == ("[0:2900](+)", ["<>"])
    assert systems_of(records["xyzq"]) == [("command_line", "domain_expr", "[1000:5900](+)")]
