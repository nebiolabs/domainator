"""Shared selection contracts for peptide and CDS-neighborhood extraction."""
import importlib
import json

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import BeforePosition, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from domainator import DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME
from domainator.errors import InputError
from domainator.seqio import parse_seqfiles, write_genbank


TOOLS = ("extract_peptides", "select_by_cds")
NAMES = ("f_A", "f_B", "r_A", "r_B")


@pytest.fixture
def record():
    features = []
    for i, name in enumerate(NAMES):
        start = i * 30
        strand = 1 if i < 2 else -1
        features.append(SeqFeature(
            SimpleLocation(start, start + 18, strand=strand), type="CDS",
            qualifiers={"gene": [name], "cds_id": [name], "translation": ["M" * 6]},
        ))
        names = (("A", "B"), ("B",), ("A",), ())[i]
        for domain in names:
            features.append(SeqFeature(
                SimpleLocation(start + 3, start + 12, strand=strand),
                type=DOMAIN_FEATURE_NAME,
                qualifiers={"name": [domain], "cds_id": [name],
                            "database": ["dbA" if domain == "A" else "dbB"],
                            "evalue": ["1e-10" if domain == "A" else "1e-2"],
                            "score": ["50"], "program": ["hmmsearch"],
                            "description": ["test domain"]},
            ))
    features.append(SeqFeature(
        SimpleLocation(33, 42, strand=1), type=DOMAIN_SEARCH_BEST_HIT_NAME,
        qualifiers={"name": ["hit"], "cds_id": ["f_B"], "database": ["search"],
                    "evalue": ["1e-10"], "score": ["80"], "program": ["hmmsearch"],
                    "description": ["test search hit"]},
    ))
    return SeqRecord(Seq("ATG" * 40), id="contig", description="selection fixture",
                     annotations={"molecule_type": "DNA"}, features=features)


def selected_names(records):
    return [next(f.qualifiers["gene"][0] for f in rec.features if f.type == "CDS")
            for rec in records]


@pytest.mark.parametrize("tool", TOOLS)
@pytest.mark.parametrize("strand", [None, "f", "r"])
@pytest.mark.parametrize("invert", [False, True])
def test_selector_strand_invert_matrix(record, tool, strand, invert):
    function = getattr(importlib.import_module(f"domainator.{tool}"), tool)
    selected = list(function([record], cds_ids={"f_A", "r_A"}, strand=strand, invert=invert))
    expected = [name for name in NAMES
                if (name.endswith("A") != invert)
                and (strand is None or name.startswith(strand))]
    assert selected_names(selected) == expected


@pytest.mark.parametrize("tool", TOOLS)
@pytest.mark.parametrize("invert", [False, True])
def test_no_selectors_library_defaults(record, tool, invert):
    function = getattr(importlib.import_module(f"domainator.{tool}"), tool)
    if invert:
        with pytest.warns(UserWarning, match="Invert supplied with no selection criteria") as caught:
            selected = list(function([record, record], strand="f", invert=True))
        assert len(caught) == 1
    else:
        selected = list(function([record, record], strand="f"))
    assert selected_names(selected) == ["f_A", "f_B"] * 2
    if tool == "extract_peptides":
        assert [rec.id for rec in selected] == ["f_A", "f_B"] * 2
        assert all(str(rec.seq) == "MMMMMM" for rec in selected)


@pytest.mark.parametrize("tool", TOOLS)
@pytest.mark.parametrize("kwargs,expected", [
    ({"domain_expr": "A & B"}, ["f_A"]),
    ({"domain_expr": "A & ~B"}, ["r_A"]),
    ({"domain_expr": "(A & B) | ~A"}, ["f_A", "f_B", "r_B"]),
    ({"domain_expr": "A & B", "cds_ids": {"r_B"}}, ["f_A", "r_B"]),
    ({"domain_expr": "A & B", "domains": {"B"}}, ["f_A", "f_B"]),
    ({"domain_expr": "A & B", "unannotated": True}, ["f_A", "r_B"]),
    ({"domain_expr": "A & B", "has_search_hit": True}, ["f_A", "f_B"]),
    ({"domain_expr": "A & ~B", "databases": {"dbA"}}, ["f_A", "r_A"]),
    ({"domain_expr": "A & ~B", "evalue": 1e-5}, ["f_A", "r_A"]),
    ({"domain_expr": "hit", "annotation_type": "domain_search"}, ["f_B"]),
    ({"domain_expr": "B & hit", "annotation_type": "all"}, ["f_B"]),
    ({"domain_expr": "A & B", "invert": True, "strand": "f"}, ["f_B"]),
    ({"domain_expr": "missing", "has_search_hit": True, "databases": {"absent"}}, ["f_B"]),
])
def test_domain_expression_parity(record, tool, kwargs, expected):
    function = getattr(importlib.import_module(f"domainator.{tool}"), tool)
    assert selected_names(function([record], **kwargs)) == expected


@pytest.mark.parametrize("tool", TOOLS)
@pytest.mark.parametrize("expression", ["", " ", "A &", "(A", "A B"])
def test_invalid_domain_expression(record, tool, expression):
    function = getattr(importlib.import_module(f"domainator.{tool}"), tool)
    with pytest.raises(InputError, match="boolean expression"):
        list(function([record], domain_expr=expression))


@pytest.mark.parametrize("tool", TOOLS)
@pytest.mark.parametrize("partial,expected", [("exclude", ["f_B", "r_A", "r_B"]),
                                            ("only", ["f_A"])])
def test_partial_is_not_inverted(record, tool, partial, expected):
    cds = record.features[0]
    cds.location = SimpleLocation(BeforePosition(0), 18, strand=1)
    function = getattr(importlib.import_module(f"domainator.{tool}"), tool)
    assert selected_names(function([record], cds_ids={"absent"}, invert=True,
                                   partial=partial)) == expected


@pytest.mark.parametrize("tool", TOOLS)
def test_name_field_matching_fallback_and_defaults(record, tool):
    module = importlib.import_module(f"domainator.{tool}")
    function = getattr(module, tool)
    cds = record.features[0]
    cds.qualifiers.update(protein_id=["WP_123"], locus_tag=["LOC_456"], custom=["CUSTOM"])
    default_name = "WP_123" if tool == "extract_peptides" else "LOC_456"
    assert selected_names(function([record], cds_ids={default_name})) == ["f_A"]
    selected = list(function([record], cds_ids={"CUSTOM", "f_B"}, name_field="custom"))
    assert selected_names(selected) == ["f_A", "f_B"]
    fallback = list(function([record], cds_ids={default_name}, name_field="missing"))
    assert selected_names(fallback) == ["f_A"]
    if tool == "extract_peptides":
        assert [rec.id for rec in selected] == ["CUSTOM", "f_B"]
        assert [rec.id for rec in function([record], cds_ids={"CUSTOM"},
                                           name_field="custom", keep_name=True)] == ["contig"]
    else:
        original = list(function([record], cds_ids={default_name}))
        assert selected[0].id == original[0].id


@pytest.mark.parametrize("tool", TOOLS)
@pytest.mark.parametrize("invert", [False, True])
@pytest.mark.parametrize("use_config", [False, True])
def test_cli_library_parity_and_contig_restriction(record, tool, invert, use_config, tmp_path):
    module = importlib.import_module(f"domainator.{tool}")
    function = getattr(module, tool)
    other = record[:]
    other.id = "other"
    infile = tmp_path / "input.gb"
    outfile = tmp_path / "output.gb"
    write_genbank([record, other], infile)
    config = {"input": [str(infile)], "genbank_out": str(outfile),
              "contigs": ["contig"], "domain_expr": "A & B",
              "name_field": "gene", "invert": invert, "strand": "f"}
    if use_config:
        configfile = tmp_path / "params.json"
        configfile.write_text(json.dumps(config))
        argv = ["--config", str(configfile)]
    else:
        argv = ["-i", str(infile), "-o", str(outfile), "--contigs", "contig",
                "--domain_expr", "A & B", "--name_field", "gene", "--strand", "f"]
        if invert:
            argv.append("--invert")
    module.main(argv)
    actual = list(parse_seqfiles([outfile]))
    expected = list(function([record], domain_expr="A & B", name_field="gene",
                             invert=invert, strand="f"))
    assert selected_names(actual) == (["f_B"] if invert else ["f_A"])
    assert [(r.id, str(r.seq)) for r in actual] == [(r.id, str(r.seq)) for r in expected]


@pytest.mark.parametrize("tool", TOOLS)
def test_cli_no_selectors_warns_once_and_matches_library(record, tool, tmp_path):
    module = importlib.import_module(f"domainator.{tool}")
    infile = tmp_path / "input.gb"
    outfile = tmp_path / "output.gb"
    write_genbank([record], infile)
    with pytest.warns(UserWarning, match="Invert supplied with no selection criteria") as caught:
        module.main(["-i", str(infile), "-o", str(outfile), "--strand", "f", "--invert"])
    assert len(caught) == 1
    actual = list(parse_seqfiles([outfile]))
    expected = list(getattr(module, tool)([record], strand="f"))
    assert [(r.id, str(r.seq)) for r in actual] == [(r.id, str(r.seq)) for r in expected]


def test_neighborhood_restrictions_apply_only_to_focus(record):
    from domainator.select_by_cds import select_by_cds

    selected = list(select_by_cds([record], cds_ids={"f_B"}, strand="f",
                                  cds_range=(1, 1)))
    assert len(selected) == 1
    assert [f.qualifiers["gene"][0] for f in selected[0].features if f.type == "CDS"] == [
        "f_A", "f_B", "r_A",
    ]
    assert next(f for f in selected[0].features
                if f.type == "CDS" and f.qualifiers["gene"] == ["r_A"]).location.strand == -1


@pytest.mark.parametrize("tool", TOOLS)
def test_explicit_empty_selector_is_not_absent(record, tool):
    function = getattr(importlib.import_module(f"domainator.{tool}"), tool)
    assert list(function([record], cds_ids=set())) == []
    assert selected_names(function([record], cds_ids=set(), invert=True,
                                   strand="r")) == ["r_A", "r_B"]


def test_peptide_expression_preserves_annotations(record):
    from domainator.extract_peptides import extract_peptides

    peptide, = extract_peptides([record], domain_expr="A & B")
    assert peptide.id == "f_A"
    assert str(peptide.seq) == "MMMMMM"
    cds = next(f for f in peptide.features if f.type == "CDS")
    assert cds.qualifiers["source_contig"] == ["contig"]
    assert "translation" not in cds.qualifiers
    domains = [f for f in peptide.features if f.type == DOMAIN_FEATURE_NAME]
    assert [f.qualifiers["name"][0] for f in domains] == ["A", "B"]
    assert [(int(f.location.start), int(f.location.end)) for f in domains] == [(1, 4), (1, 4)]
