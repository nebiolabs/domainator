"""Taxonomy of records and regions, and filtering records by taxonomy.
"""
from typing import Optional

from Bio.SeqRecord import SeqRecord
from domainator.locations import feature_intervals, location_covers
from domainator.model import get_sources
from domainator.taxonomy.ncbi import NCBITaxonomy


class TaxonomyData:
    """
    Dataclass for storing pre_computed taxonomy data
    """
    # lineage: List[str]
    # rank: str
    def __init__(self, ncbi_taxonomy:NCBITaxonomy, record:SeqRecord=None, taxid:int=None):
        if record is None and taxid is None:
            raise RuntimeError("Must provide either a SeqRecord or a taxid")
        if record is not None and taxid is not None:
            raise RuntimeError("Must provide either a SeqRecord or a taxid, not both")
        if record is not None:
            self.taxid = get_taxid(record)
        else:
            self.taxid = taxid

        self.taxid = ncbi_taxonomy.normalized_taxid(self.taxid)
        self.lineage = ncbi_taxonomy.lineage(self.taxid)
        self.lineage.reverse() # reverse so that self has the highest index
        self.ranks = [ ncbi_taxonomy.rank(x).lower() for x in self.lineage ]
        self.names = [ ncbi_taxonomy.name(x) for x in self.lineage ]
        self.rank_to_name = dict(zip(self.ranks, self.names))
        self.rank_to_taxid = dict(zip(self.ranks, self.lineage))

    def __str__(self):
        return f"TaxonomyData(taxid={self.taxid}, lineage={self.lineage}, ranks={self.ranks}, names={self.names})"
    def __repr__(self):
        return str(self)

def get_taxid(record:SeqRecord) -> Optional[int]:
    """
        returns 32644, "unidentified" if no taxid is found
    """
    taxid = 32644
    if "ncbi_taxid" in record.annotations: # from SwissProt-style text file
        return int(record.annotations["ncbi_taxid"][0])
    elif " OX=" in record.description: # from SwissProt-style fasta file
        return int(record.description.split(" OX=")[-1].split(" ")[0])
    elif " TaxID=" in record.description: # from Uniref-style fasta file
        return int(record.description.split(" TaxID=")[-1].split(" ")[0])
    else: # from GenBank-style genbank file, take taxid from the longest source feature
        sources = get_sources(record)
        sources.sort(key=len, reverse=True)
        if len(sources) > 0:
            if "db_xref" in sources[0].qualifiers:
                for xref in sources[0].qualifiers["db_xref"]:
                    if xref.startswith("taxon:"):
                        try:
                            return int(xref.split(":")[-1])
                        except ValueError:
                            pass
    return taxid

_UNIDENTIFIED_TAXID = 32644  # "unidentified", get_taxid's default for records lacking a taxid


def _feature_taxon(feature) -> Optional[int]:
    for xref in feature.qualifiers.get("db_xref", []):
        if xref.startswith("taxon:"):
            try:
                return int(xref.split(":")[-1])
            except ValueError:
                pass
    return None


def location_taxid(region_intervals, source_features) -> Optional[int]:
    """Taxid of the longest `source` feature whose location *covers* `region_intervals`.

    `region_intervals` is a list of (start, end) half-open intervals (a CDS, a whole contig
    ``[(0, len)]``, or a hit region); `source_features` are the record's "source" features
    (Biopython or LeanFeature, e.g. from get_sources). Returns the covering source's taxon,
    _UNIDENTIFIED_TAXID if a source covers it but carries no taxon, or None if no source
    covers it. None is distinguishable so a whole-contig coverage test can tell a
    single-source contig (non-None) from a gappy multi-source one (None). Mirrors
    get_taxid's longest-source rule, evaluated against a specific region.
    """
    ranked = sorted(
        source_features,
        key=lambda f: sum(e - s for s, e in feature_intervals(f)),
        reverse=True,
    )
    for src in ranked:
        if location_covers(feature_intervals(src), region_intervals):
            taxon = _feature_taxon(src)
            return taxon if taxon is not None else _UNIDENTIFIED_TAXID
    return None


def taxid_allowed(taxid, allowed_taxids) -> bool:
    """True if `taxid` (None -> unidentified 32644) is in the allowed-taxid set."""
    return (_UNIDENTIFIED_TAXID if taxid is None else taxid) in allowed_taxids


def source_taxids(record):
    """Taxids of the `source` features in a record (Biopython SeqRecord or LeanContig);
    sources without a taxon db_xref are omitted. Mirrors the Rust LeanSearchContig.source_taxids
    accessor (used to decide whether a gappy multi-source contig could yield an allowed hit)."""
    return [t for t in (_feature_taxon(s) for s in get_sources(record)) if t is not None]


def record_matches_taxonomy(record: SeqRecord, include_taxids, exclude_taxids, ncbi_taxonomy) -> bool:
    taxid = get_taxid(record)
    if taxid is None:
        return False

    lineage = set(ncbi_taxonomy.lineage(taxid))
    if include_taxids and not lineage.intersection(include_taxids):
        return False
    if exclude_taxids and lineage.intersection(exclude_taxids):
        return False
    return True


def compile_taxonomy_allowed(ncbi_taxonomy, include_taxids=None, exclude_taxids=None, taxonomy_expr=None):
    """Compile taxonomy filtering criteria into a frozenset of allowed taxids.

    taxonomy_expr is a boolean expression over taxids (operators & | ~ and
    parentheses; see NCBITaxonomy.compile_taxonomy_expression) and is mutually
    exclusive with include_taxids/exclude_taxids. The verdict for every taxid is
    computed in a single O(N) pass, so callers can filter records by an O(1)
    membership test. Returns None when no filtering criteria are given.
    """
    if taxonomy_expr and (include_taxids or exclude_taxids):
        raise ValueError("taxonomy_expr is mutually exclusive with include_taxids/exclude_taxids.")
    if taxonomy_expr:
        return ncbi_taxonomy.compile_taxonomy_expression(taxonomy_expr)
    if include_taxids or exclude_taxids:
        return ncbi_taxonomy.compile_taxonomy_filter(include_taxids, exclude_taxids)
    return None


def filter_by_taxonomy(records, include_taxids, exclude_taxids, ncbi_taxonomy, taxonomy_expr=None):
    allowed = compile_taxonomy_allowed(ncbi_taxonomy, include_taxids, exclude_taxids, taxonomy_expr)
    if allowed is None:
        yield from records
        return
    for record in records:
        taxid = get_taxid(record)
        if taxid is not None and taxid in allowed:
            yield record
