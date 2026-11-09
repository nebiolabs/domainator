"""Systems already annotated on records, for the commands that select by them."""
from typing import Dict, Iterable, List, NamedTuple, Optional, Sequence, Set, Tuple

from Bio.SeqFeature import SeqFeature

from domainator import SYSTEM_SEARCH_BEST_HIT_NAME
from domainator.location_ops import overlap_size
from domainator.model import DomainatorCDS, copy_feature, system_feature_types
from domainator.systems.features import component_values

MEMBERSHIP_QUALIFIER = "system_membership"
"""The qualifier that extract_peptides puts on the CDS feature of a protein from a system member."""
NUCLEOTIDE_MEMBER_PREFIX = "nuc_"
"""The cds_id prefix of a contig-level nucleotide member: its coordinates when it was detected."""


def system_features(record, system_type: str = "domain_system") -> List[SeqFeature]:
    """The record's systems of ``system_type`` (see model.SYSTEM_TYPES), in feature order."""
    types = system_feature_types(system_type)
    return [feature for feature in record.features if feature.type in types]


def system_strand(feature: SeqFeature) -> int:
    """The orientation of a system: the strand of its location, or 1 if it has none."""
    return feature.location.parts[0].strand or 1


def member_indices(feature: SeqFeature, cds_list: Sequence[DomainatorCDS]) -> List[int]:
    """The indices into cds_list (from DomainatorCDS.list_from_contig with
    include_nucleic_acid_annotations) of a system's members, sorted.

    CDS members are found by cds_id. A nucleotide member's id is its coordinates when it was
    detected, which slicing changes, so it is found by its hit's database and name instead,
    among the nucleotide annotations that overlap the system.
    """
    values = component_values(feature)
    cds_ids = {value.cds_id for value in values if not value.cds_id.startswith(NUCLEOTIDE_MEMBER_PREFIX)}
    nucleotide_hits = {(value.database, value.reference) for value in values
                       if value.cds_id.startswith(NUCLEOTIDE_MEMBER_PREFIX)}
    indices = list()
    for index, cds in enumerate(cds_list):
        if cds.is_nucleic_acid:
            if nucleotide_hits and overlap_size(cds.feature.location, feature.location) > 0 and any(
                    (hit.qualifiers["database"][0], hit.qualifiers["name"][0]) in nucleotide_hits
                    for hit in _hits(cds)):
                indices.append(index)
        elif cds.num in cds_ids:
            indices.append(index)
    return indices


def _hits(cds: DomainatorCDS):
    yield from cds.domain_features
    if cds.domain_search_feature is not None:
        yield cds.domain_search_feature


class MemberHits:
    """The domain hits that are members' evidence in some systems: ``hit in member_hits``.

    A hit on a CDS is evidence when the system's component qualifiers name its cds_id,
    database and name. A contig-level nucleotide hit is evidence when they name its database
    and name for a nucleotide member, and it overlaps the system.
    """

    def __init__(self, features: Iterable[SeqFeature]):
        self.cds_hits: Set[Tuple[str, str, str]] = set()
        self.nucleotide_hits: List[Tuple[str, str, SeqFeature]] = list()
        for feature in features:
            for value in component_values(feature):
                if value.cds_id.startswith(NUCLEOTIDE_MEMBER_PREFIX):
                    self.nucleotide_hits.append((value.database, value.reference, feature))
                else:
                    self.cds_hits.add((value.cds_id, value.database, value.reference))

    def __contains__(self, hit: SeqFeature) -> bool:
        database, name = hit.qualifiers["database"][0], hit.qualifiers["name"][0]
        if (hit.qualifiers.get("cds_id", ["."])[0], database, name) in self.cds_hits:
            return True
        return any((database, name) == (system_database, system_name) and overlap_size(hit.location, system.location) > 0
                   for system_database, system_name, system in self.nucleotide_hits)


def system_completeness(feature: SeqFeature) -> Optional[float]:
    """The completeness qualifier, or None for systems without one (models without a quorum)."""
    value = feature.qualifiers.get("completeness", ["."])[0]
    return None if value == "." else float(value)


def as_search_hit(feature: SeqFeature) -> SeqFeature:
    """A System_Search copy of a system feature, to mark it as the focus of a region."""
    if feature.type == SYSTEM_SEARCH_BEST_HIT_NAME:
        return feature
    search_hit = copy_feature(feature)
    search_hit.type = SYSTEM_SEARCH_BEST_HIT_NAME
    return search_hit



class Membership(NamedTuple):
    """One system a CDS is a member of, and the component it fills."""
    role: str
    component: str
    database: str
    name: str

    def format(self) -> str:
        """The system_membership qualifier value, ``role|component|database|name``."""
        return "|".join(self)

    @property
    def names(self) -> Set[str]:
        """The names a selector matches this system by: ``name`` and ``database:name``."""
        return {self.name, f"{self.database}:{self.name}"}


def memberships(record, system_type: str = "domain_system") -> Dict[str, List[Membership]]:
    """The systems each CDS of record is a member of, by cds_id, from its system features."""
    out: Dict[str, List[Membership]] = dict()
    for feature in system_features(record, system_type):
        database, name = feature.qualifiers["database"][0], feature.qualifiers["name"][0]
        for value in component_values(feature):
            if not value.cds_id.startswith(NUCLEOTIDE_MEMBER_PREFIX):
                membership = Membership(value.role, value.component, database, name)
                if membership not in out.setdefault(value.cds_id, []):
                    out[value.cds_id].append(membership)
    return out
