"""Tests for taxonomy.records: the taxid of a record, and filtering records by taxonomy."""
import shutil
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from domainator.taxonomy import NCBITaxonomy
from domainator.taxonomy.records import (TaxonomyData, compile_taxonomy_allowed, filter_by_taxonomy, get_taxid,
                                         source_taxids)

E_COLI, B_SUBTILIS, BACTERIA, UNIDENTIFIED = 562, 1423, 2, 32644


@pytest.fixture(scope="module")
def ncbi_taxonomy(tmp_path_factory):
    taxdmp = tmp_path_factory.mktemp("taxonomy") / "taxdmp"
    shutil.copytree(Path(__file__).parent / "data" / "taxdmp", taxdmp)
    return NCBITaxonomy(taxdmp, overwrite=False)


def _source(start, end, *xrefs):
    return SeqFeature(SimpleLocation(start, end, 1), type="source", qualifiers={"db_xref": list(xrefs)})


def _record(description="", annotations=None, features=(), length=100, id="rec"):
    return SeqRecord(Seq("A" * length), id=id, description=description, annotations=annotations or {},
                     features=list(features))


@pytest.mark.parametrize("record,taxid", [
    (_record(annotations={"ncbi_taxid": ["562"]}), E_COLI),  # SwissProt text
    (_record(description="sp|P0AGD3|SODF_ECOLI Superoxide dismutase OS=Escherichia coli OX=562 GN=sodB PE=1"), E_COLI),
    (_record(description="UniRef50_P0AGD3 Superoxide dismutase n=1 Tax=Escherichia coli TaxID=562 RepID=SODF_ECOLI"), E_COLI),
    # GenBank: the longest source feature's taxon
    (_record(features=[_source(0, 10, "taxon:1423"), _source(0, 100, "GOA:x", "taxon:562")]), E_COLI),
    (_record(features=[_source(0, 100, "taxon:unknown")]), UNIDENTIFIED),
    (_record(features=[_source(0, 100)]), UNIDENTIFIED),
    (_record(), UNIDENTIFIED),
])
def test_get_taxid(record, taxid):
    assert get_taxid(record) == taxid


def test_source_taxids_skips_sources_without_a_taxon():
    record = _record(features=[_source(0, 10, "taxon:1423"), _source(10, 100), _source(10, 50, "taxon:562")])
    assert source_taxids(record) == [B_SUBTILIS, E_COLI]


def test_taxonomy_data(ncbi_taxonomy):
    from_taxid = TaxonomyData(ncbi_taxonomy, taxid=E_COLI)
    from_record = TaxonomyData(ncbi_taxonomy, record=_record(annotations={"ncbi_taxid": ["562"]}))
    assert str(from_taxid) == str(from_record) == repr(from_record)
    assert from_taxid.lineage[0] == 1 and from_taxid.lineage[-1] == E_COLI
    assert from_taxid.rank_to_name["species"] == "Escherichia coli"
    assert from_taxid.rank_to_taxid["superkingdom" if "superkingdom" in from_taxid.ranks else "domain"] == BACTERIA
    with pytest.raises(RuntimeError, match="either a SeqRecord or a taxid"):
        TaxonomyData(ncbi_taxonomy)
    with pytest.raises(RuntimeError, match="not both"):
        TaxonomyData(ncbi_taxonomy, record=_record(), taxid=E_COLI)


def test_filter_by_taxonomy(ncbi_taxonomy):
    records = [_record(id="coli", annotations={"ncbi_taxid": ["562"]}),
               _record(id="subtilis", annotations={"ncbi_taxid": ["1423"]}),
               _record(id="unknown")]

    def ids(**kwargs):
        args = {"include_taxids": None, "exclude_taxids": None, **kwargs}
        return [r.id for r in filter_by_taxonomy(records, ncbi_taxonomy=ncbi_taxonomy, **args)]

    assert ids() == ["coli", "subtilis", "unknown"]
    assert ids(include_taxids=[BACTERIA]) == ["coli", "subtilis"]
    assert ids(include_taxids=[BACTERIA], exclude_taxids=[B_SUBTILIS]) == ["coli"]
    assert ids(exclude_taxids=[E_COLI]) == ["subtilis"]  # the test taxonomy has no 32644 ("unidentified")
    assert ids(taxonomy_expr=f"{BACTERIA} & ~{E_COLI}") == ["subtilis"]


def test_compile_taxonomy_allowed_rejects_an_expression_with_include_or_exclude(ncbi_taxonomy):
    assert compile_taxonomy_allowed(ncbi_taxonomy) is None
    with pytest.raises(ValueError, match="mutually exclusive"):
        compile_taxonomy_allowed(ncbi_taxonomy, include_taxids=[BACTERIA], taxonomy_expr=str(E_COLI))
