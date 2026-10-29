"""Shared list/file normalization without changing tool-specific selection roles."""
import importlib
import json

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from domainator.cli.common import list_and_file_to_dict_keys
from domainator.seqio import parse_seqfiles, write_genbank


TOOLS = ("extract_domains", "extract_peptides", "extract_unannotated",
         "select_by_cds", "select_by_contig")
DOMAIN_TOOLS = tuple(tool for tool in TOOLS if tool != "extract_unannotated")


def test_list_file_normalization(tmp_path):
    names = tmp_path / "names.txt"
    names.write_text(" B \n\nA\nB\n C\t\n")
    assert list(list_and_file_to_dict_keys(["A", "A"], names)) == ["A", "B", "C"]
    assert list_and_file_to_dict_keys(["A"], names, as_set=True) == {"A", "B", "C"}
    assert list_and_file_to_dict_keys([" A ", "a"], as_set=True) == {" A ", "a"}
    assert list_and_file_to_dict_keys() is None
    assert list_and_file_to_dict_keys([], as_set=True) is None
    names.write_text(" \n\t\n")
    assert list_and_file_to_dict_keys(input_file=names) == {}
    assert list_and_file_to_dict_keys(input_file=names, as_set=True) == set()
    assert list_and_file_to_dict_keys(["A"], names, as_set=True) == {"A"}


@pytest.mark.parametrize("suffix", ["fasta", "gb"])
def test_sequence_file_names(tmp_path, suffix):
    path = tmp_path / f"names.{suffix}"
    if suffix == "fasta":
        path.write_text(">A description\nMA\n>B\nMG\n>A duplicate\nMA\n")
    else:
        records = [SeqRecord(Seq("MA"), id=name, annotations={"molecule_type": "protein"})
                   for name in ("A", "B", "A")]
        with path.open("w") as handle:
            write_genbank(records, handle)
    assert list_and_file_to_dict_keys(["C"], path, as_set=True) == {"A", "B", "C"}


def test_hmm_file_names(shared_datadir):
    assert list_and_file_to_dict_keys(
        ["other", "CcdB"], shared_datadir / "CcdB.hmm", as_set=True) == {"other", "CcdB"}


@pytest.fixture
def inputs(tmp_path):
    paths = {}
    for protein in (False, True):
        records = []
        for label, sequence in (("A", "ATG"), ("B", "GCT")):
            name = f"record_{label}"
            translation = "M" * 10 if label == "A" else "A" * 10
            record = SeqRecord(Seq(translation if protein else sequence * 10), id=name,
                               annotations={"molecule_type": "protein" if protein else "DNA"})
            if not protein:
                record.features.append(SeqFeature(
                    SimpleLocation(0, 30, strand=1), type="CDS",
                    qualifiers={"gene": [name], "cds_id": [name], "translation": [translation]}))
            record.features.append(SeqFeature(
                SimpleLocation(3, 6 if protein else 12, strand=1), type="Domainator",
                qualifiers={"name": [f"domain_{label}"], "database": ["test"],
                            "cds_id": ["0_1_10" if protein else name], "score": ["50"],
                            "evalue": ["1e-10"], "description": ["test"], "program": ["hmmsearch"]}))
            records.append(record)
        path = tmp_path / ("proteins.gb" if protein else "contigs.gb")
        with path.open("w") as handle:
            write_genbank(records, handle)
        paths[protein] = path
    return paths


def run_tool(tool, inputs, tmp_path, options, *, config=False):
    module = importlib.import_module(f"domainator.{tool}")
    output = tmp_path / "output.gb"
    args = ["-i", str(inputs[tool == "extract_unannotated"]), "-o", str(output)]
    if config:
        path = tmp_path / "config.json"
        path.write_text(json.dumps(options))
        args += ["--config", str(path)]
    else:
        for name, value in options.items():
            args.append(f"--{name}")
            if value is not True:
                args.extend(str(item) for item in (value if isinstance(value, list) else [value]))
    module.main(args)
    return [(record.id, str(record.seq)) for record in parse_seqfiles([output])]


@pytest.mark.parametrize("tool", TOOLS)
def test_existing_list_defaults(tool):
    params = importlib.import_module(f"domainator.{tool}").build_parser().parse_args(
        ["-i", "input.gb", "-o", "output.gb"])
    assert params.contigs == ([] if tool.startswith("extract_") else None)
    assert params.contigs_file is None
    if tool in DOMAIN_TOOLS:
        assert params.domains is None
        assert params.domains_file is None


@pytest.mark.parametrize("tool", TOOLS)
@pytest.mark.parametrize("config", [False, True])
def test_contig_union_and_empty_restriction(tool, inputs, tmp_path, config):
    names = tmp_path / "contigs.txt"
    names.write_text("record_B\n record_A \nrecord_B\n\n")
    baseline = run_tool(tool, inputs, tmp_path, {})
    assert baseline
    assert run_tool(tool, inputs, tmp_path,
                    {"contigs": ["record_A"], "contigs_file": str(names)},
                    config=config) == baseline
    names.write_text("")
    assert run_tool(tool, inputs, tmp_path, {"contigs_file": str(names)}, config=config) == []
    if tool in ("select_by_contig", "select_by_cds", "extract_peptides"):
        expected = baseline if tool == "select_by_contig" else []
        inverted = {"contigs_file": str(names), "invert": True}
        if tool == "select_by_contig":
            assert run_tool(tool, inputs, tmp_path, inverted, config=config) == expected
        else:
            with pytest.warns(UserWarning, match="Invert supplied with no selection criteria"):
                assert run_tool(tool, inputs, tmp_path, inverted, config=config) == expected
    direct = run_tool(tool, inputs, tmp_path, {"contigs": ["record_A"]})
    assert direct
    assert run_tool(tool, inputs, tmp_path,
                    {"contigs": ["record_A"], "contigs_file": str(names)},
                    config=config) == direct


@pytest.mark.parametrize("tool", DOMAIN_TOOLS)
@pytest.mark.parametrize("config", [False, True])
def test_domain_union_and_empty_selector(tool, inputs, tmp_path, config):
    names = tmp_path / "domains.txt"
    names.write_text(" domain_B \n\ndomain_A\ndomain_B\n")
    baseline = run_tool(tool, inputs, tmp_path, {})
    assert baseline
    assert run_tool(tool, inputs, tmp_path,
                    {"domains": ["domain_A"], "domains_file": str(names)},
                    config=config) == baseline
    names.write_text("")
    assert run_tool(tool, inputs, tmp_path, {"domains_file": str(names)}, config=config) == []
    direct = run_tool(tool, inputs, tmp_path, {"domains": ["domain_A"]})
    assert direct
    assert run_tool(tool, inputs, tmp_path,
                    {"domains": ["domain_A"], "domains_file": str(names)},
                    config=config) == direct
    if tool != "extract_domains":
        assert run_tool(tool, inputs, tmp_path,
                        {"domains_file": str(names), "invert": True}, config=config) == baseline


@pytest.mark.parametrize("tool,invert", [
    (tool, invert) for tool in DOMAIN_TOOLS for invert in (False, True)
    if not (tool == "extract_domains" and invert)
])
def test_record_selector_versus_input_restriction(tool, invert, inputs, tmp_path):
    options = {"contigs": ["record_A"], "domains": ["domain_B"]}
    if invert:
        options["invert"] = True
    records = run_tool(tool, inputs, tmp_path, options)
    if tool == "select_by_contig":
        expected = [] if invert else [
            ("record_A", "ATG" * 10), ("record_B", "GCT" * 10)]
    elif invert:
        expected = [("record_A", "M" * 10 if tool == "extract_peptides" else "ATG" * 10)]
    else:
        expected = []
    assert records == expected
