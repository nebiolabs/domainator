"""Detected systems as GenBank features (Domain_System and System_Search).

See docs/design/systems/data_model.md and the "Systems" section of docs/file_formats.md.
"""
from typing import Dict, List, NamedTuple

from Bio.SeqFeature import SeqFeature

from domainator import SYSTEM_FEATURE_NAME
from domainator.errors import InputError
from domainator.systems.model import SystemHit

COMPONENT_QUALIFIER = "component"


class ComponentValue(NamedTuple):
    """One value of the repeated ``component`` qualifier: a (member, component) pair."""
    role: str
    component: str
    cds_id: str
    database: str
    reference: str

    def format(self) -> str:
        for field in ("role", "component", "cds_id", "database"):
            if "|" in getattr(self, field):
                raise InputError(f"Can't write system component {self!r}: its {field} contains '|'.")
        return "|".join(self)


def parse_component(value: str) -> ComponentValue:
    """A ``component`` qualifier value, ``role|component|cds_id|database|reference``."""
    fields = value.split("|", 4)
    if len(fields) != 5:
        raise InputError(f"Malformed system component qualifier {value!r}: expected role|component|cds_id|database|reference.")
    return ComponentValue(*fields)


def system_feature(hit: SystemHit, feature_type: str = SYSTEM_FEATURE_NAME) -> SeqFeature:
    """``hit`` as a Domain_System (or, with feature_type, System_Search) feature."""
    model = hit.model
    qualifiers: Dict[str, List[str]] = {
        "name": [model.name],
        "database": [model.database],
        "program": [model.program],
        "category": [model.category],
        "description": [model.description],
        "accession": [model.accession],
        "window": [model.window.describe()],
        "score": [f"{hit.score:.1f}"],
        "evalue": [f"{hit.evalue:.1e}"],
        "completeness": ["." if hit.completeness is None else f"{hit.completeness:.3f}"],
        "core_count": [str(hit.core_count)],
        "total_count": [str(hit.total_count)],
        "member_count": [str(hit.member_count)],
    }
    if model.neighborhood_kb is not None:
        qualifiers["neighborhood_kb"] = [f"{model.neighborhood_kb:g}"]
    if model.neighborhood_cds is not None:
        qualifiers["neighborhood_cds"] = [str(model.neighborhood_cds)]
    if hit.partial:
        qualifiers["system_partial"] = [hit.partial]
    qualifiers[COMPONENT_QUALIFIER] = [
        ComponentValue(member.role, member.component.name, member.cds_id, member.hit.qualifiers["database"][0],
                       member.hit.qualifiers["name"][0]).format()
        for member in hit.members]
    return SeqFeature(hit.location, type=feature_type, qualifiers=qualifiers)


def component_values(feature: SeqFeature) -> List[ComponentValue]:
    """The members of a Domain_System or System_Search feature."""
    return [parse_component(value) for value in feature.qualifiers.get(COMPONENT_QUALIFIER, [])]
