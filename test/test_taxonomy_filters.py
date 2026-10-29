"""Shared taxonomy contract across extraction, selection, and search."""
import importlib
import inspect
import json
import pickle
import shutil
from types import SimpleNamespace
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from domainator import DOMAIN_FEATURE_NAME
from domainator.errors import InputError
from domainator.seqio import parse_seqfiles, write_genbank
from domainator.taxonomy import NCBITaxonomy
from domainator.taxonomy.records import compile_taxonomy_filter


TOOLS = ("extract_domains", "extract_peptides", "extract_unannotated", "select_by_cds", "select_by_contig")
FEATURE_TOOLS = ("extract_domains", "extract_peptides", "select_by_cds")


@pytest.fixture(scope="module")
def taxonomy_path(tmp_path_factory):
    directory = tmp_path_factory.mktemp("filter-taxonomy") / "taxdmp"
    shutil.copytree(Path(__file__).parent / "data" / "taxdmp", directory)
    return directory


@pytest.fixture(scope="module")
def taxonomy(taxonomy_path):
    return NCBITaxonomy(taxonomy_path, overwrite=False)


def source(start, end, taxid):
    qualifiers = {} if taxid is None else {"db_xref": [f"taxon:{taxid}"]}
    return SeqFeature(SimpleLocation(start, end), type="source", qualifiers=qualifiers)


def domain(start, end, cds_id="cds"):
    return SeqFeature(SimpleLocation(start, end, strand=1), type=DOMAIN_FEATURE_NAME,
                      qualifiers={"name": ["hit"], "database": ["db"], "evalue": ["1e-10"],
                                  "score": ["50"], "cds_id": [cds_id], "description": ["test hit"],
                                  "program": ["hmmsearch"]})


def protein(taxid=562, id="protein"):
    features = [source(0, 100, taxid), domain(10, 30)]
    return SeqRecord(Seq("M" * 100), id=id, description=id,
                     annotations={"molecule_type": "protein"}, features=features)


def mixed():
    host = source(0, 90, 562)
    host.location = CompoundLocation([SimpleLocation(0, 90), SimpleLocation(180, 300)], operator="order")
    features = [host, source(90, 180, 1423)]
    for i, (start, end) in enumerate(((0, 60), (90, 150), (165, 225))):
        cds_id = f"cds{i}"
        features.extend([
            SeqFeature(SimpleLocation(start, end, strand=1), type="CDS",
                       qualifiers={"protein_id": [cds_id], "cds_id": [cds_id], "translation": ["M" * 20]}),
            domain(start + 3, end - 3, cds_id),
        ])
    return SeqRecord(Seq("ATG" * 100), id="mixed", description="mixed",
                     annotations={"molecule_type": "DNA"}, features=features)


@pytest.mark.parametrize("policy", ["exclude", "include", "evaluate", "error"])
@pytest.mark.parametrize("criteria,empty", [
    ({"taxonomy_expr": "~2"}, True),
    ({"exclude_taxids": [2]}, True),
    ({"include_taxids": [2]}, False),
    ({"taxonomy_expr": "2 & ~1224"}, False),
])
def test_unknown_policies(taxonomy, policy, criteria, empty):
    gate = compile_taxonomy_filter(taxonomy, taxonomy_unknown=policy, **criteria)
    gate = pickle.loads(pickle.dumps(gate))
    for taxid in (None, 32644, 999999999):
        if policy == "error":
            with pytest.raises(InputError, match="record 'unknown'"):
                gate.allows(taxid, "record 'unknown'")
        else:
            assert gate.allows(taxid) == (policy == "include" or (policy == "evaluate" and empty))


def test_known_merged_deleted_and_criteria(taxonomy):
    gate = compile_taxonomy_filter(taxonomy, include_taxids=[2], exclude_taxids=[1224])
    assert gate.allows(1423)
    assert not gate.allows(562)
    for old, row in taxonomy._merged_tab.items():
        replacement = int(row["new_tax_id"])
        assert gate.allows(old) == gate.allows(replacement)
    for deleted in taxonomy._deleted_tab:
        assert not gate.allows(int(deleted))
        assert compile_taxonomy_filter(taxonomy, taxonomy_expr="~2",
                                       taxonomy_unknown="evaluate").allows(int(deleted))
    old = next((t for t, row in taxonomy._merged_tab.items()
                if int(row["new_tax_id"]) in taxonomy._nodes_tab), None)
    if old is not None:
        assert compile_taxonomy_filter(taxonomy, taxonomy_expr=str(old)) == compile_taxonomy_filter(
            taxonomy, include_taxids=[int(taxonomy._merged_tab[old]["new_tax_id"])])
    for kwargs in ({"include_taxids": [999999999]}, {"taxonomy_expr": "Bacteria"}):
        with pytest.raises(InputError, match="criterion"):
            compile_taxonomy_filter(taxonomy, **kwargs)


@pytest.mark.parametrize("tool", TOOLS)
def test_no_filter_is_noop_without_database(tool):
    function = getattr(importlib.import_module(f"domainator.{tool}"), tool)
    records = [mixed()] if tool in FEATURE_TOOLS else [protein()]
    assert [(r.id, str(r.seq)) for r in function(records, taxonomy_expr="  ")] == [
        (r.id, str(r.seq)) for r in function(records)]
    with pytest.raises(InputError, match="requires"):
        list(function(records, taxonomy_expr="~2"))
    with pytest.raises(InputError, match="Invalid taxonomy_unknown"):
        list(function(records, taxonomy_unknown="typo"))
    with pytest.raises(InputError, match="mutually exclusive"):
        list(function(records, include_taxids=[2], taxonomy_expr="~2"))


@pytest.mark.parametrize("tool", FEATURE_TOOLS)
def test_feature_taxonomy_precedes_extraction(taxonomy, tool):
    function = getattr(importlib.import_module(f"domainator.{tool}"), tool)
    kwargs = {"extract_all": True} if tool == "extract_peptides" else {}
    records = list(function([mixed()], ncbi_taxonomy=taxonomy, include_taxids=[1423], **kwargs))
    assert len(records) == 1
    assert records[0].id == ("cds1" if tool == "extract_peptides" else
                             "mixed_94:147" if tool == "extract_domains" else "mixed_91:150")
    assert len(records[0]) == (20 if tool == "extract_peptides" else
                               54 if tool == "extract_domains" else 60)
    assert records[0].annotations["molecule_type"] == ("protein" if tool == "extract_peptides" else "DNA")


@pytest.mark.parametrize("tool", ["extract_peptides", "select_by_cds"])
@pytest.mark.parametrize("invert", [False, True])
def test_taxonomy_only_selection_uses_core_defaults(taxonomy, tool, invert):
    function = getattr(importlib.import_module(f"domainator.{tool}"), tool)
    kwargs = {"ncbi_taxonomy": taxonomy, "include_taxids": [1423], "invert": invert}
    if invert:
        with pytest.warns(UserWarning, match="Invert supplied with no selection criteria") as caught:
            records = list(function([mixed()], **kwargs))
        assert len(caught) == 1
    else:
        records = list(function([mixed()], **kwargs))
    assert [(r.id, len(r)) for r in records] == [
        ("cds1", 20) if tool == "extract_peptides" else ("mixed_91:150", 60)
    ]


@pytest.mark.parametrize("tool", FEATURE_TOOLS)
def test_boundary_unknown_and_inversion(taxonomy, tool):
    function = getattr(importlib.import_module(f"domainator.{tool}"), tool)
    kwargs = {"extract_all": True} if tool == "extract_peptides" else {}
    assert list(function([mixed()], ncbi_taxonomy=taxonomy, taxonomy_expr="~2", **kwargs)) == []
    records = list(function([mixed()], ncbi_taxonomy=taxonomy, taxonomy_expr="~2",
                            taxonomy_unknown="evaluate", **kwargs))
    assert len(records) == 1
    with pytest.raises(InputError, match="mixed.*(CDS|Domainator)"):
        list(function([mixed()], ncbi_taxonomy=taxonomy, taxonomy_expr="~2",
                      taxonomy_unknown="error", **kwargs))
    if tool != "extract_domains":
        assert not list(function([mixed()], ncbi_taxonomy=taxonomy, taxonomy_expr="~2",
                                 cds_ids={"nonexistent"}, invert=True))


@pytest.mark.parametrize("mode", ["combine", "splice"])
def test_only_allowed_domains_are_combined(taxonomy, mode):
    from domainator.extract_domains import extract_domains
    records = list(extract_domains([mixed()], ncbi_taxonomy=taxonomy, include_taxids=[1423],
                                   extend_range=40, **{mode: True}))
    assert len(records) == 1
    assert len(records[0]) == 134


def test_neighborhood_taxonomy_is_anchor_not_output(taxonomy):
    from domainator.select_by_cds import select_by_cds
    records = list(select_by_cds([mixed()], ncbi_taxonomy=taxonomy, include_taxids=[1423],
                                 whole_contig=True))
    assert len(records) == 1
    assert str(records[0].seq) == str(mixed().seq)


@pytest.mark.parametrize("tool", ["extract_domains", "extract_unannotated", "select_by_contig"])
def test_protein_taxonomy_policy_and_malformed(taxonomy, tool):
    function = getattr(importlib.import_module(f"domainator.{tool}"), tool)
    assert not list(function([protein(None)], ncbi_taxonomy=taxonomy, taxonomy_expr="~2"))
    assert list(function([protein(None)], ncbi_taxonomy=taxonomy, taxonomy_expr="~2",
                         taxonomy_unknown="evaluate"))
    with pytest.raises(InputError, match="Invalid taxonomy.*protein"):
        list(function([protein("bad")], ncbi_taxonomy=taxonomy, taxonomy_expr="~2"))


@pytest.mark.parametrize("tool", TOOLS)
def test_cli_library_taxonomy_parity(taxonomy, taxonomy_path, tmp_path, tool):
    module = importlib.import_module(f"domainator.{tool}")
    input_path, output_path = tmp_path / "input.gb", tmp_path / "out.gb"
    record = mixed() if tool in FEATURE_TOOLS else protein()
    with input_path.open("w") as handle:
        write_genbank([record], handle)
    config = {"input": [str(input_path)], "genbank_out": str(output_path),
              "taxonomy_expr": "2 & ~1224", "ncbi_taxonomy_path": str(taxonomy_path)}
    config_path = tmp_path / "params.json"
    config_path.write_text(json.dumps(config))
    module.main(["--config", str(config_path)])
    kwargs = {"extract_all": True} if tool == "extract_peptides" else {}
    expected = list(getattr(module, tool)([record], ncbi_taxonomy=taxonomy,
                                         taxonomy_expr="2 & ~1224", **kwargs))
    actual = list(parse_seqfiles([output_path]))
    assert [(r.id, str(r.seq)) for r in actual] == [(r.id, str(r.seq)) for r in expected]


@pytest.mark.parametrize("tool", ["domainate", "structure_domainate"])
@pytest.mark.parametrize("keyword", ["hits_only", "no_annotations"])
def test_annotation_api_removes_search_options(tool, keyword):
    function = getattr(importlib.import_module(f"domainator.{tool}"), tool)
    assert keyword not in inspect.signature(function).parameters
    kwargs = {"z": 1000} if tool == "domainate" else {"aligner": None}
    with pytest.raises(TypeError, match=keyword):
        function([], references=[], **kwargs, **{keyword: True})


@pytest.mark.parametrize("flag", ["--hits_only", "--no_annotations", "--include_taxids",
                                 "--exclude_taxids", "--taxonomy_expr", "--taxonomy_update",
                                 "--ncbi_taxonomy_path"])
def test_domainate_rejects_removed_flags(flag):
    from domainator.domainate import build_parser
    with pytest.raises(SystemExit):
        build_parser().parse_args(["-i", "input.gb", "-r", "reference.hmm", "-o", "out.gb", flag])


def test_lazy_cli_taxonomy_setup(monkeypatch):
    from domainator.cli.common import taxonomy_from_args

    def unexpected_load(*args, **kwargs):
        pytest.fail("Inactive taxonomy criteria must not initialize a database")

    monkeypatch.setattr("domainator.taxonomy.NCBITaxonomy", unexpected_load)
    assert taxonomy_from_args(SimpleNamespace(include_taxids=None, exclude_taxids=None,
                                               taxonomy_expr="  ", taxonomy_update=True)) is None


@pytest.mark.parametrize("tool", TOOLS)
def test_compiles_once_per_operation(taxonomy, monkeypatch, tool):
    module = importlib.import_module(f"domainator.{tool}")
    calls = []

    def compile_once(*args):
        calls.append(True)
        return compile_taxonomy_filter(*args)

    monkeypatch.setattr(module, "compile_taxonomy_filter", compile_once)
    kwargs = {"extract_all": True} if tool == "extract_peptides" else {}
    records = [mixed(), mixed()] if tool in FEATURE_TOOLS else [protein(), protein()]
    list(getattr(module, tool)(records, ncbi_taxonomy=taxonomy, include_taxids=[2], **kwargs))
    assert len(calls) == 1


@pytest.mark.parametrize("policy", ["exclude", "include", "evaluate", "error"])
def test_native_and_python_feature_policy_parity(taxonomy, tmp_path, policy):
    from domainator.domainate import get_prot_list
    from domainator.lean_record import iter_lean_search_native
    from domainator.seqio import native_parser_available

    if not native_parser_available():
        pytest.skip("native parser not built")
    record = mixed()
    path = tmp_path / "mixed.gb"
    with path.open("w") as handle:
        write_genbank([record], handle)
    native = next(iter_lean_search_native(str(path), 0, -1, "nucleotide"))
    gate = compile_taxonomy_filter(taxonomy, taxonomy_expr="~2", taxonomy_unknown=policy)
    if policy == "error":
        for candidate in (record, native):
            with pytest.raises(InputError, match="mixed.*CDS"):
                list(get_prot_list(candidate, 0, allowed_taxids=gate))
    else:
        python_peptides = list(get_prot_list(record, 0, allowed_taxids=gate))
        assert list(get_prot_list(native, 0, allowed_taxids=gate)) == python_peptides
        assert len(python_peptides) == (0 if policy == "exclude" else 1)


def test_nucleotide_error_prefilter_does_not_invent_unknown_candidates(taxonomy):
    from domainator.domainate import _nucleotide_contig_searchable
    record = mixed()
    record.features = record.features[:2]
    gate = compile_taxonomy_filter(taxonomy, taxonomy_expr="~2", taxonomy_unknown="error")
    # No whole-contig source exists, but both sources are known. Check actual hits later.
    assert _nucleotide_contig_searchable(record, gate)


@pytest.mark.parametrize("method", ["threads", "spawn"])
@pytest.mark.parametrize("policy", ["exclude", "evaluate", "error"])
def test_search_worker_unknown_policy(taxonomy, shared_datadir, monkeypatch, method, policy):
    from domainator.domain_search import domain_search
    monkeypatch.setenv("DOMAINATOR_MP_START_METHOD", method)
    reference = shared_datadir / "swissprot_CuSOD_subset.fasta"
    record = next(parse_seqfiles([reference]))
    # Same sequence, but no taxid: a negated criterion must not admit it by default.
    record.id = "unlabeled"
    record.description = "unlabeled"
    record.annotations.pop("ncbi_taxid", None)
    record.features = []
    iterator = domain_search([record], references=[str(reference)], z=1000, cpu=2,
                             ncbi_taxonomy=taxonomy, taxonomy_expr="~2", taxonomy_unknown=policy)
    if policy == "error":
        with pytest.raises(InputError, match="unlabeled"):
            list(iterator)
    else:
        assert [r.id for r in iterator] == ([] if policy == "exclude" else ["unlabeled"])


def test_nucleotide_taxonomy_before_hit_limit(taxonomy, monkeypatch):
    from domainator import domainate
    from domainator.model import get_hit_best_score

    record = mixed()
    hits = [domainate.SearchResult("unknown", "unknown", ".", 1e-10, 1000, 170, 210, "db", 1, 0, 40),
            domainate.SearchResult("coli", "coli", ".", 1e-10, 900, 10, 50, "db", 1, 0, 40),
            domainate.SearchResult("subtilis", "subtilis", ".", 1e-10, 50, 100, 140, "db", 1, 0, 40)]
    monkeypatch.setattr(domainate, "run_search", lambda *args, **kwargs: {0: {-1: hits}})
    domainate.domainator_inner(
        [record], [], [], [], [], {}, evalue=0.001, cpu=1, z=1000,
        hits_only=True, no_annotations=True, max_hits=1, max_overlap=1,
        best_annotation=True, max_hits_per_contig=1,
        allowed_taxids=compile_taxonomy_filter(taxonomy, include_taxids=[1423]),
    )
    selected = [f for f in record.features if f.type == "Domain_Search"]
    assert len(selected) == 1
    assert selected[0].qualifiers["name"] == ["subtilis"]
    assert (int(selected[0].location.start), int(selected[0].location.end)) == (100, 140)
    assert get_hit_best_score(record) == 50


@pytest.mark.parametrize("strand", [1, -1])
def test_circular_anchor_source_precedence(taxonomy, strand):
    record = mixed()
    record.annotations["topology"] = "circular"
    wrapped = SeqFeature(CompoundLocation([SimpleLocation(270, 300, strand=strand),
                                           SimpleLocation(0, 30, strand=strand)]),
                         type="CDS")
    gate = compile_taxonomy_filter(taxonomy, include_taxids=[562])
    assert gate.matches(record, wrapped)
    # A full-contig source takes precedence over nested/gappy sources.
    record.features.append(source(0, 300, 1423))
    assert not gate.matches(record, wrapped)
    assert compile_taxonomy_filter(taxonomy, include_taxids=[1423]).matches(record, wrapped)


@pytest.mark.parametrize("tool", ["domainate", "structure_domainate"])
@pytest.mark.parametrize("keyword", ["hits_only", "no_annotations"])
def test_removed_annotation_config_keys(tmp_path, tool, keyword):
    module = importlib.import_module(f"domainator.{tool}")
    path = tmp_path / "params.json"
    path.write_text(json.dumps({"input": ["in.gb"], "reference": ["ref.hmm"],
                                "genbank_out": "out.gb", keyword: True}))
    with pytest.raises(SystemExit):
        module.build_parser().parse_args(["--config", str(path)])
