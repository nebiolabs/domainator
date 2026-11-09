"""Turn antiSMASH's cluster calls (in its GenBank output) into Domain_System features.

Checked against antiSMASH 8.0.4 output. antiSMASH writes, per record:

- ``protocluster`` features (``/protocluster_number``, ``/product``, ``/category``, ``/cutoff``,
  ``/neighbourhood``, ``/detection_rule``, ``/contig_edge``), each with a ``proto_core`` feature
  (the core, with the same ``/protocluster_number``);
- ``cand_cluster`` features (``/protoclusters``: protocluster numbers, ``/kind``) and
  ``region`` features (``/candidate_cluster_numbers``);
- on core CDSs, ``/gene_functions="biosynthetic (rule-based-clusters) <product>: <profile>"``
  naming the profiles that triggered the product, and
  ``/sec_met_domain="<profile> (E-value: <e>, bitscore: <s>, ...)"`` with their scores.

A system's members are the CDSs in its protocluster cores that triggered its products. Each
fills one component per triggering profile, with the role core. Profiles are referenced by name
only (database ``antismash``): antiSMASH gives no positions within the protein, so there are
no matching Domainator features. See docs/design/systems/antismash.md.
"""
import re
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from domainator import SYSTEM_FEATURE_NAME
from domainator.location_ops import overlap_size
from domainator.model import DomainatorCDS, get_cds_unique_name
from domainator.systems.features import system_feature
from domainator.systems.model import Component, ReferenceRule, SystemHit, SystemMember, SystemModel, Window

DATABASE = "antismash"
LEVELS = ("protocluster", "candidate_cluster", "region")
SOURCE_FEATURE_TYPES = ("region", "cand_cluster", "protocluster", "proto_core")
"""The antiSMASH features that import_systems --clear_source removes."""

_SEC_MET_DOMAIN = re.compile(r"^(\S+) \(E-value: ([^,]+), bitscore: ([^,)]+)")
_CORE_FUNCTION = re.compile(r"^biosynthetic \(rule-based-clusters\) (.+): (\S+)$")


def _first(feature: SeqFeature, key: str, default: Optional[str] = None) -> Optional[str]:
    return feature.qualifiers.get(key, [default])[0]


def _number(value: Optional[str]) -> Optional[float]:
    try:
        return None if value is None else float(value)
    except ValueError:
        return None


class _Protocluster:
    """A protocluster, its core, and the members of its core."""

    def __init__(self, record: SeqRecord, protocluster: SeqFeature, core: SeqFeature, cdss: Sequence[DomainatorCDS]):
        self.feature = protocluster
        self.core = core
        self.number = _first(protocluster, "protocluster_number")
        self.product = _first(protocluster, "product", ".")
        self.category = _first(protocluster, "category", ".")
        self.cutoff = _number(_first(protocluster, "cutoff"))
        self.neighbourhood = _number(_first(protocluster, "neighbourhood"))
        self.rule = " ".join(_first(protocluster, "detection_rule", ".").split())
        self.contig_edge = _first(protocluster, "contig_edge", "False") == "True"
        # (cds, profile, (evalue, score)) for each profile a core CDS contributes to the product
        self.members: List[Tuple[DomainatorCDS, str, Tuple[float, float]]] = list()
        for cds in cdss:
            if not overlap_size(core.location, cds.feature.location):
                continue
            scores = _sec_met_domains(cds.feature)
            for function in cds.feature.qualifiers.get("gene_functions", []):
                match = _CORE_FUNCTION.match(" ".join(function.split()))
                if match and match.group(1) == self.product:
                    profile = match.group(2)
                    self.members.append((cds, profile, scores.get(profile, (float("nan"), 0.0))))


def _sec_met_domains(feature: SeqFeature) -> Dict[str, Tuple[float, float]]:
    """{profile: (e-value, bit score)} from a CDS's sec_met_domain qualifiers."""
    out = dict()
    for value in feature.qualifiers.get("sec_met_domain", []):
        match = _SEC_MET_DOMAIN.match(" ".join(value.split()))
        if match:
            evalue, score = _number(match.group(2)), _number(match.group(3))
            out[match.group(1)] = (float("nan") if evalue is None else evalue, 0.0 if score is None else score)
    return out


def _partial(record: SeqRecord, location, contig_edge: bool) -> str:
    """system_partial for a feature antiSMASH marks as touching a contig edge."""
    if not contig_edge or record.annotations.get("topology") == "circular":
        return ""
    left = int(location.start) == 0
    right = int(location.end) >= len(record)
    return ("<" if left else "") + (">" if right else "")


def _hit(record: SeqRecord, protoclusters: Sequence[_Protocluster], location, *, name: str) -> Optional[SystemHit]:
    members = list()
    components: Dict[str, Component] = dict()
    seen = set()
    best: Dict[str, Tuple[float, float]] = dict()  # by cds_id: the best (score, evalue) of its profiles
    for protocluster in protoclusters:
        for cds, profile, (evalue, score) in protocluster.members:
            if profile not in components:
                components[profile] = Component(profile, "core", (ReferenceRule(profile, database=DATABASE),))
            if (cds.num, profile) in seen:
                continue
            seen.add((cds.num, profile))
            hit = SeqFeature(cds.feature.location, type="sec_met_domain",
                             qualifiers={"database": [DATABASE], "name": [profile]})
            members.append(SystemMember(cds, cds.num, components[profile], "core", hit))
            if score > best.get(cds.num, (float("-inf"), 0.0))[0]:
                best[cds.num] = (score, evalue)
    if not members:
        return None
    members.sort(key=lambda member: int(member.cds.feature.location.start))
    cutoffs = [p.cutoff for p in protoclusters if p.cutoff is not None]
    neighbourhoods = [p.neighbourhood for p in protoclusters if p.neighbourhood is not None]
    categories = list(dict.fromkeys(p.category for p in protoclusters))
    model = SystemModel(
        name=name, database=DATABASE, program="antismash", components=tuple(components.values()),
        window=Window("anchor", max_kb=max(cutoffs) / 1000 if cutoffs else None),
        category="+".join(categories), neighborhood_kb=max(neighbourhoods) / 1000 if neighbourhoods else None,
        description="; ".join(dict.fromkeys(p.rule for p in protoclusters)))
    evalues = [evalue for _, evalue in best.values() if evalue == evalue]  # not nan
    return SystemHit(model=model, members=tuple(members), location=location,
                     score=sum(score for score, _ in best.values()), evalue=max(evalues) if evalues else float("nan"),
                     completeness=None, core_count=len(components), total_count=len(components),
                     partial=_partial(record, location, any(p.contig_edge for p in protoclusters)))


def import_antismash(record: SeqRecord, level: str = "protocluster", clear_source: bool = False) -> SeqRecord:
    """Add a Domain_System feature for each antiSMASH protocluster (or candidate cluster, or
    region) of record, replacing systems imported before. CDSs get cds_id qualifiers, which the
    systems' component qualifiers refer to."""
    if level not in LEVELS:
        raise ValueError(f"level must be one of {', '.join(LEVELS)}, not {level!r}")
    for feature in record.features:
        if feature.type == "CDS" and "cds_id" not in feature.qualifiers:
            feature.qualifiers["cds_id"] = [get_cds_unique_name(feature)]
    cdss = DomainatorCDS.list_from_contig(record)
    cores = {_first(f, "protocluster_number"): f for f in record.features if f.type == "proto_core"}
    protoclusters = {_first(f, "protocluster_number"): _Protocluster(record, f, cores[_first(f, "protocluster_number")], cdss)
                     for f in record.features
                     if f.type == "protocluster" and _first(f, "protocluster_number") in cores}
    candidates = {_first(f, "candidate_cluster_number"): f for f in record.features if f.type == "cand_cluster"}

    def products(features: Iterable[_Protocluster]) -> str:
        return "+".join(dict.fromkeys(p.product for p in features))

    systems = list()
    if level == "protocluster":
        for protocluster in protoclusters.values():
            hit = _hit(record, [protocluster], protocluster.core.location, name=protocluster.product)
            if hit is not None:
                systems.append(system_feature(hit))
    elif level == "candidate_cluster":
        for candidate in candidates.values():
            members = [protoclusters[n] for n in candidate.qualifiers.get("protoclusters", []) if n in protoclusters]
            hit = _hit(record, members, candidate.location, name=products(members))
            if hit is not None:
                feature = system_feature(hit)
                feature.qualifiers["kind"] = [_first(candidate, "kind", ".")]
                systems.append(feature)
    else:
        for region in (f for f in record.features if f.type == "region"):
            numbers = [n for c in region.qualifiers.get("candidate_cluster_numbers", []) if c in candidates
                       for n in candidates[c].qualifiers.get("protoclusters", [])]
            members = [protoclusters[n] for n in dict.fromkeys(numbers) if n in protoclusters]
            hit = _hit(record, members, region.location, name=products(members))
            if hit is not None:
                systems.append(system_feature(hit))

    record.features = [feature for feature in record.features
                       if not (feature.type == SYSTEM_FEATURE_NAME and _first(feature, "database") == DATABASE)
                       and not (clear_source and feature.type in SOURCE_FEATURE_TYPES)]
    record.features.extend(systems)
    return record
