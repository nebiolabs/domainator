"""Taxonomy of records and regions, and filtering records by taxonomy.
"""
from dataclasses import dataclass
import re
from typing import FrozenSet, Iterable, Optional

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from domainator.locations import feature_intervals, location_covers
from domainator.model import get_sources
from domainator.taxonomy.ncbi import NCBITaxonomy
from domainator.errors import InputError
from domainator.expr import BooleanEvaluator


_UNIDENTIFIED_TAXID = 32644  # "unidentified", get_taxid's default for records lacking a taxid
TAXONOMY_UNKNOWN_CHOICES = ("exclude", "include", "evaluate", "error")


@dataclass(frozen=True)
class TaxonomyFilter:
    """A compiled, picklable taxonomy gate; no database or lineage walks in workers."""

    allowed: FrozenSet[int]
    resolved: FrozenSet[int]
    unknown: str
    empty_matches: bool

    def allows(self, taxid: Optional[int], context: str = "candidate") -> bool:
        if taxid in self.resolved:
            return taxid in self.allowed
        if self.unknown == "error":
            raise InputError(f"Unresolved taxonomy for {context}: taxid {taxid!r} "
                             "(missing, unidentified, deleted, or absent from the taxonomy database).")
        return self.unknown == "include" or (self.unknown == "evaluate" and self.empty_matches)

    def __contains__(self, taxid: int) -> bool:
        return self.allows(taxid)

    def matches(self, record: SeqRecord, feature: Optional[SeqFeature] = None) -> bool:
        context = f"record {record.id!r}"
        if feature is not None:
            context += f", {feature.type} {feature.location}"
        try:
            if feature is None or record.annotations.get("molecule_type") == "protein":
                taxid = get_taxid(record)
                sources = get_sources(record)
                if ("ncbi_taxid" in record.annotations or " OX=" in record.description
                        or " TaxID=" in record.description):
                    sources = []
            else:
                sources = get_sources(record)
                taxid = location_taxid(feature_intervals(feature), sources)
            # The older resolver treats malformed source taxids as missing. New filtering
            # must report malformed labels rather than applying the unknown policy to them.
            if feature is not None and record.annotations.get("molecule_type") != "protein":
                sources = [source for source in sources
                           if location_covers(feature_intervals(source), feature_intervals(feature))]
            for source in sorted(sources, key=lambda source: sum(e - s for s, e in feature_intervals(source)),
                                 reverse=True)[:1]:
                for xref in source.qualifiers.get("db_xref", []):
                    if xref.startswith("taxon:"):
                        int(xref.split(":")[-1])
        except (ValueError, IndexError) as error:
            raise InputError(f"Invalid taxonomy for {context}: {error}") from error
        return self.allows(taxid, context)


def taxonomy_active(include_taxids=None, exclude_taxids=None, taxonomy_expr=None) -> bool:
    """Whether nonempty taxonomy criteria were supplied."""
    return bool(include_taxids or exclude_taxids or (taxonomy_expr and taxonomy_expr.strip()))


def compile_taxonomy_filter(
    ncbi_taxonomy: Optional[NCBITaxonomy], include_taxids: Optional[Iterable[int]] = None,
    exclude_taxids: Optional[Iterable[int]] = None, taxonomy_expr: Optional[str] = None,
    taxonomy_unknown: str = "exclude",
) -> Optional[TaxonomyFilter]:
    """Compile a taxonomy gate with explicit handling of unresolved candidates."""
    if taxonomy_unknown not in TAXONOMY_UNKNOWN_CHOICES:
        raise InputError(f"Invalid taxonomy_unknown {taxonomy_unknown!r}; choose from {TAXONOMY_UNKNOWN_CHOICES}.")
    include = set(include_taxids) if include_taxids is not None else set()
    exclude = set(exclude_taxids) if exclude_taxids is not None else set()
    expression = taxonomy_expr.strip() if taxonomy_expr else None
    if expression and (include or exclude):
        raise InputError("taxonomy_expr is mutually exclusive with include_taxids/exclude_taxids.")
    if not taxonomy_active(include, exclude, expression):
        return None
    if ncbi_taxonomy is None:
        raise InputError("Taxonomy filtering requires an ncbi_taxonomy database.")

    def normalize(taxid):
        try:
            taxid = int(taxid)
        except (ValueError, TypeError) as error:
            raise InputError(f"Invalid taxonomy criterion: {taxid!r}. Expected an NCBI taxid.") from error
        seen = set()
        while taxid in ncbi_taxonomy._merged_tab:
            if taxid in seen:
                raise InputError(f"Cyclic merged taxonomy criterion: {taxid}.")
            seen.add(taxid)
            taxid = int(ncbi_taxonomy._merged_tab[taxid]["new_tax_id"])
        if taxid not in ncbi_taxonomy._nodes_tab:
            raise InputError(f"Taxonomy criterion {taxid} is deleted or absent from the database; "
                             "check the taxid or use --taxonomy_update.")
        return taxid

    include = {normalize(t) for t in include}
    exclude = {normalize(t) for t in exclude}
    if expression:
        evaluator = BooleanEvaluator(expression)
        replacements = {name: str(normalize(name)) for name in evaluator.names}
        expression = re.sub(r"[^&|~()\s]+", lambda match: replacements[match.group()], expression)
        empty_matches = BooleanEvaluator(expression).check_expression(set())
    else:
        empty_matches = not include
    allowed = set(compile_taxonomy_allowed(ncbi_taxonomy, include, exclude, expression))
    resolved = set(ncbi_taxonomy._nodes_tab) - {_UNIDENTIFIED_TAXID}
    for old in ncbi_taxonomy._merged_tab:
        replacement = old
        seen = set()
        while replacement in ncbi_taxonomy._merged_tab and replacement not in seen:
            seen.add(replacement)
            replacement = int(ncbi_taxonomy._merged_tab[replacement]["new_tax_id"])
        if replacement in resolved:
            resolved.add(old)
            if replacement in allowed:
                allowed.add(old)
    return TaxonomyFilter(frozenset(allowed), frozenset(resolved), taxonomy_unknown, empty_matches)


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

def get_taxid(record:SeqRecord) -> int:
    """The record's NCBI taxid, or 32644 ("unidentified") if it has none."""
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
    return _UNIDENTIFIED_TAXID

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


def taxid_allowed(taxid, allowed_taxids, context="candidate") -> bool:
    """True if `taxid` (None -> unidentified 32644) is in the allowed-taxid set."""
    if isinstance(allowed_taxids, TaxonomyFilter):
        return allowed_taxids.allows(taxid, context)
    return (_UNIDENTIFIED_TAXID if taxid is None else taxid) in allowed_taxids


def source_taxids(record):
    """Taxids of the `source` features in a record (Biopython SeqRecord or LeanContig);
    sources without a taxon db_xref are omitted. Mirrors the Rust LeanSearchContig.source_taxids
    accessor (used to decide whether a gappy multi-source contig could yield an allowed hit)."""
    return [t for t in (_feature_taxon(s) for s in get_sources(record)) if t is not None]


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
        if taxid in allowed:
            yield record
