"""Systems already annotated on records, for the commands that select by them."""
from typing import List, Optional, Sequence

from Bio.SeqFeature import SeqFeature

from domainator import SYSTEM_SEARCH_BEST_HIT_NAME
from domainator.location_ops import overlap_size
from domainator.model import DomainatorCDS, copy_feature, system_feature_types
from domainator.systems.features import component_values

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

