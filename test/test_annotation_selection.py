"""Annotation sources and the independent legacy-compatible search-hit selector."""
import importlib

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from domainator import DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME
from domainator.extract_domains import extract_domains
from domainator.extract_unannotated import extract_unannotated
from domainator.filter_domains import filter_domains
from domainator.model import DomainatorCDS
from domainator.rotate_contigs import rotate_contigs
from domainator.seqio import parse_seqfiles, write_genbank


def _annotation(start, end, kind, name, cds_id, *, evalue="0.001", score="100"):
    return SeqFeature(FeatureLocation(start, end, strand=1), type=kind, qualifiers={
        "name": [name], "cds_id": [cds_id], "evalue": [evalue], "score": [score],
        "database": ["domains" if kind == DOMAIN_FEATURE_NAME else "searches"],
        "description": [name],
    })


@pytest.fixture
def mixed():
    record = SeqRecord(Seq("ATG" * 50), id="mixed",
                       annotations={"molecule_type": "DNA", "topology": "circular"})
    for i, name in enumerate(["domain_only", "search_only", "both", "neither", "weak"]):
        start = i * 30
        record.features.append(SeqFeature(FeatureLocation(start, start + 30, strand=1),
                                         type="CDS", qualifiers={
                                             "cds_id": [name], "protein_id": [name],
                                             "translation": ["M" * 10],
                                         }))
        if name in ("domain_only", "both"):
            record.features.append(_annotation(start + 3, start + 9, DOMAIN_FEATURE_NAME, "A", name))
        if name in ("search_only", "both", "weak"):
            record.features.append(_annotation(start + 15, start + 21, DOMAIN_SEARCH_BEST_HIT_NAME,
                                               "H", name, evalue="0.1" if name == "weak" else "0.001"))
    return record


@pytest.fixture
def protein():
    record = SeqRecord(Seq("M" * 30), id="protein", annotations={"molecule_type": "protein"})
    record.features = [
        _annotation(3, 9, DOMAIN_FEATURE_NAME, "A", "."),
        _annotation(15, 21, DOMAIN_SEARCH_BEST_HIT_NAME, "H", "."),
    ]
    return record


def _selected(tool, record, **kwargs):
    module = importlib.import_module(f"domainator.{tool}")
    records = [record]
    if tool == "select_by_contig":
        records = [record[i:i + 30] for i in range(0, len(record), 30)]
        for rec in records:
            rec.id = next(f.qualifiers["protein_id"][0] for f in rec.features if f.type == "CDS")
    output = list(getattr(module, tool)(records, **kwargs))
    if tool == "select_by_cds":
        return {next(f.qualifiers["protein_id"][0] for f in rec.features if f.type == "CDS")
                for rec in output}
    return {rec.id for rec in output}


@pytest.mark.parametrize("tool", ["select_by_cds", "extract_peptides", "select_by_contig"])
@pytest.mark.parametrize("annotation_type,domains,expected", [
    ("domainator", {"A"}, {"domain_only", "both"}),
    ("domainator", {"H"}, set()),
    ("domain_search", {"A"}, set()),
    ("domain_search", {"H"}, {"search_only", "both"}),
    ("all", {"A", "H"}, {"domain_only", "search_only", "both"}),
])
def test_annotation_sources(mixed, tool, annotation_type, domains, expected):
    assert _selected(tool, mixed, annotation_type=annotation_type, domains=domains,
                     evalue=0.001) == expected


@pytest.mark.parametrize("tool", ["select_by_cds", "extract_peptides", "select_by_contig"])
@pytest.mark.parametrize("annotation_type", ["domainator", "domain_search", "all"])
@pytest.mark.parametrize("invert", [False, True])
def test_search_hit_bypasses_names_and_databases(mixed, tool, annotation_type, invert):
    expected = {"domain_only", "neither", "weak"} if invert else {"search_only", "both"}
    assert _selected(tool, mixed, annotation_type=annotation_type, domains={"absent"},
                     databases={"absent"}, has_search_hit=True, invert=invert,
                     evalue=0.001) == expected


@pytest.mark.parametrize("tool", ["select_by_cds", "extract_peptides", "select_by_contig"])
@pytest.mark.parametrize("invert", [False, True])
def test_search_hit_alone_is_a_selector(mixed, tool, invert):
    expected = {"domain_only", "neither", "weak"} if invert else {"search_only", "both"}
    assert _selected(tool, mixed, has_search_hit=True, evalue=0.001, invert=invert) == expected


@pytest.mark.parametrize("tool", ["select_by_cds", "extract_peptides", "select_by_contig"])
def test_search_hit_unions_with_normal_selection(mixed, tool):
    assert _selected(tool, mixed, domains={"A"}, has_search_hit=True,
                     evalue=0.001) == {"domain_only", "search_only", "both"}


@pytest.mark.parametrize("tool", ["select_by_cds", "extract_peptides", "select_by_contig"])
@pytest.mark.parametrize("annotation_type", ["domain_search", "all"])
def test_database_filter_applies_without_search_hit_selector(mixed, tool, annotation_type):
    assert _selected(tool, mixed, domains={"H"}, databases={"domains"},
                     annotation_type=annotation_type) == set()


@pytest.mark.parametrize("tool", ["select_by_cds", "extract_peptides", "select_by_contig"])
@pytest.mark.parametrize("annotation_type,expected", [
    ("domainator", {"search_only", "neither", "weak"}),
    ("domain_search", {"domain_only", "neither", "weak"}),
    ("all", {"neither", "weak"}),
])
def test_unannotated_uses_selected_sources(mixed, tool, annotation_type, expected):
    assert _selected(tool, mixed, unannotated=True, annotation_type=annotation_type,
                     evalue=0.001) == expected


@pytest.mark.parametrize("tool", ["select_by_cds", "select_by_contig"])
def test_expression_uses_union_before_evaluation(mixed, tool):
    assert _selected(tool, mixed, domain_expr="A & H", annotation_type="all") == {"both"}
    assert _selected(tool, mixed, domain_expr="A & H", annotation_type="domain_search") == set()
    assert _selected(tool, mixed, domain_expr="absent", has_search_hit=True,
                     databases={"absent"}, evalue=0.001) == {"search_only", "both"}


@pytest.mark.parametrize("tool", ["select_by_cds", "select_by_contig"])
def test_source_alone_does_not_select_hits(mixed, tool):
    assert _selected(tool, mixed, annotation_type="domain_search") == {
        "domain_only", "search_only", "both", "neither", "weak",
    }


@pytest.mark.parametrize("tool", ["select_by_cds", "extract_peptides", "select_by_contig"])
@pytest.mark.parametrize("invert", [False, True])
def test_search_hit_still_respects_partial(mixed, tool, invert):
    assert _selected(tool, mixed, has_search_hit=True, partial="only", invert=invert) == set()


def test_nucleotide_search_hit_sources_and_selector():
    record = SeqRecord(Seq("A" * 30), id="nucleotide",
                       annotations={"molecule_type": "DNA"})
    record.features = [
        _annotation(3, 9, DOMAIN_FEATURE_NAME, "A", "."),
        _annotation(15, 21, DOMAIN_SEARCH_BEST_HIT_NAME, "H", "."),
    ]
    from domainator.select_by_cds import select_by_cds

    assert [r.id for r in select_by_cds([record], has_search_hit=True)] == ["nucleotide_16:21"]
    assert [r.id for r in select_by_cds([record], annotation_type="all", domains={"H"},
                                       include_nucleic_acids=True)] == ["nucleotide_16:21"]
    assert list(select_by_cds([record], has_search_hit=True, invert=True)) == []


def test_search_hit_contig_score_cutoff(mixed):
    assert _selected("select_by_contig", mixed, has_search_hit=True, min_search_score=101,
                     databases={"absent"}) == set()


@pytest.mark.parametrize("annotation_type,expected", [
    ("domainator", {"protein_4:9"}),
    ("domain_search", {"protein_16:21"}),
    ("all", {"protein_4:9", "protein_16:21"}),
])
def test_extract_domain_sources(protein, annotation_type, expected):
    assert {rec.id for rec in extract_domains([protein], annotation_type=annotation_type)} == expected
    expected_named = {"protein_16:21"} if annotation_type != "domainator" else set()
    assert {rec.id for rec in extract_domains([protein], annotation_type=annotation_type,
                                              domains={"H"})} == expected_named


def test_extract_domain_search_hit_selector(protein):
    assert [r.id for r in extract_domains([protein], has_search_hit=True)] == ["protein_16:21"]
    assert [r.id for r in extract_domains([protein], has_search_hit=True, domains={"absent"},
                                          databases={"absent"})] == ["protein_16:21"]
    assert list(extract_domains([protein], has_search_hit=True, min_score=101)) == []
    assert list(extract_domains([protein], has_search_hit=True, evalue=0.0001)) == []
    assert {r.id for r in extract_domains([protein], has_search_hit=True,
                                         domains={"A"})} == {"protein_4:9", "protein_16:21"}


@pytest.mark.parametrize("annotation_type,expected", [
    ("domainator", {"protein_1:3", "protein_10:30"}),
    ("domain_search", {"protein_1:15", "protein_22:30"}),
    ("all", {"protein_1:3", "protein_10:15", "protein_22:30"}),
])
def test_unannotated_region_sources(protein, annotation_type, expected):
    assert {rec.id for rec in extract_unannotated([protein], annotation_type=annotation_type)} == expected


@pytest.mark.parametrize("annotation_type,expected", [
    ("domainator", {DOMAIN_SEARCH_BEST_HIT_NAME}),
    ("domain_search", {DOMAIN_FEATURE_NAME}),
    ("all", set()),
])
def test_filter_domain_sources(protein, annotation_type, expected):
    filtered = next(filter_domains([protein], annotation_type=annotation_type, domains={"absent"}))
    assert {f.type for f in filtered.features} == expected


@pytest.mark.parametrize("annotation_type", ["domain_search", "all"])
def test_rotation_source_filters(mixed, annotation_type):
    result = next(rotate_contigs([mixed], annotation_type=annotation_type, domains={"H"},
                                databases={"searches"}, evalue=0.001))
    assert result.id == "mixed_46:45"
    with pytest.warns(UserWarning, match="No qualifying domain"):
        assert list(rotate_contigs([mixed], annotation_type=annotation_type, domains={"H"},
                                  databases={"domains"}, evalue=0.001)) == []


def test_rotation_search_hit_bypasses_names_and_databases(mixed):
    result = next(rotate_contigs([mixed], has_search_hit=True, domains={"absent"},
                                databases={"absent"}, evalue=0.001))
    assert result.id == "mixed_46:45"
    with pytest.warns(UserWarning, match="No qualifying domain"):
        assert list(rotate_contigs([mixed], has_search_hit=True, evalue=0.0001)) == []


@pytest.mark.parametrize("tool", ["select_by_cds", "extract_peptides", "select_by_contig",
                                 "extract_domains", "rotate_contigs"])
def test_cli_search_hit_forwarding(mixed, tool, tmp_path):
    module = importlib.import_module(f"domainator.{tool}")
    source, output = tmp_path / "input.gb", tmp_path / "output.gb"
    write_genbank([mixed], source)
    module.main(["-i", str(source), "-o", str(output), "--has_search_hit",
                 "--domains", "absent", "--databases", "absent", "-e", "0.001"])
    assert len(list(parse_seqfiles([output]))) == {
        "select_by_cds": 2, "extract_peptides": 2, "select_by_contig": 1,
        "extract_domains": 2, "rotate_contigs": 1,
    }[tool]


@pytest.mark.parametrize("tool", ["select_by_cds", "extract_peptides", "select_by_contig"])
def test_cli_search_hit_inversion(mixed, tool, tmp_path):
    module = importlib.import_module(f"domainator.{tool}")
    source, output = tmp_path / "input.gb", tmp_path / "output.gb"
    write_genbank([mixed], source)
    module.main(["-i", str(source), "-o", str(output), "--has_search_hit", "--invert",
                 "-e", "0.001"])
    assert len(list(parse_seqfiles([output]))) == (0 if tool == "select_by_contig" else 3)


@pytest.mark.parametrize("tool", ["select_by_cds", "extract_peptides", "select_by_contig",
                                 "extract_domains", "rotate_contigs"])
@pytest.mark.parametrize("annotation_type", ["domain_search", "all"])
def test_cli_source_and_database_forwarding(mixed, tool, annotation_type, tmp_path):
    module = importlib.import_module(f"domainator.{tool}")
    source, output = tmp_path / "input.gb", tmp_path / "output.gb"
    write_genbank([mixed], source)
    module.main(["-i", str(source), "-o", str(output), "--annotation_type", annotation_type,
                 "--domains", "H", "--databases", "searches", "-e", "0.001"])
    assert len(list(parse_seqfiles([output]))) == {
        "select_by_cds": 2, "extract_peptides": 2, "select_by_contig": 1,
        "extract_domains": 2, "rotate_contigs": 1,
    }[tool]


def test_cds_names_keep_default_and_filter_selected_sources(mixed):
    cds = next(c for c in DomainatorCDS.list_from_contig(mixed) if c.name == "both")
    assert cds.get_domain_names() == {"A"}
    assert cds.get_domain_names(annotation_type="all") == {"A", "H"}
    assert cds.get_domain_names({"domains"}, annotation_type="all") == {"A"}
    assert cds.get_domain_names(evalue=0.0001, annotation_type="all") == set()


@pytest.mark.parametrize("tool", ["select_by_cds", "select_by_contig", "extract_peptides",
                                 "extract_domains", "extract_unannotated", "rotate_contigs",
                                 "filter_domains"])
def test_all_tools_accept_and_validate_sources(mixed, tool):
    module = importlib.import_module(f"domainator.{tool}")
    action = next(a for a in module.build_parser()._actions if a.dest == "annotation_type")
    assert action.choices == ["domainator", "domain_search", "all"]
    with pytest.raises(ValueError, match="annotation.type"):
        list(getattr(module, tool)([mixed], annotation_type="invalid"))
