"""System detection: a pure function of one annotated record and the model collections.

Follows docs/design/systems/detection.md. It reads Domainator (and, with annotation_type,
Domain_Search) features and never searches.
"""
from dataclasses import dataclass, replace
import math
from typing import Collection, Dict, FrozenSet, List, Optional, Sequence, Set, Tuple
import warnings

from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from domainator.location_ops import circular_envelope, overlap_size
from domainator.model import (DomainatorCDS, annotation_feature_types, check_annotation_type, evalue_passes, feature_evalue,
                              feature_score, get_cds_unique_name)
from domainator import DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME
from domainator.partial import PARTIAL_CHOICES, cds_fragment_status, fragment_status_allowed
from domainator.systems import windows
from domainator.systems.conditions import Condition
from domainator.systems.model import (COUNTED_ROLES, FEATURE_DATABASE, ROLES, STRUCTURAL_THRESHOLDS, Component,
                                      ModelCollection, ReferenceRule, SystemHit, SystemMember, SystemModel)

_ROLE_PRIORITY = {role: i for i, role in enumerate(ROLES)}
TRUNCATED_CHOICES = ("exclude", "include")
_FEATURE_HIT = "system_feature_component"
"""The type of the stand-in hit a record feature gets when a feature component can match it."""
_SCORE_TIE = 1e-8  # hits within this many bits of the best in a competing group are all kept


@dataclass
class _Gene:
    """A CDS, or the pseudo-CDS of a contig-level nucleotide hit, in contig order."""
    index: int
    cds: DomainatorCDS
    start: int
    end: int
    strand: int
    is_cds: bool
    rank: int
    cds_id: str
    length: int  # bp of the CDS (or hit)
    hits: List[SeqFeature]


@dataclass(frozen=True)
class _Match:
    component: Component
    hit: SeqFeature
    score: float
    evalue: float
    positions: Tuple[int, ...] = ()
    """Where every passing hit to the component starts along the protein, N- to C-terminal."""


@dataclass
class _Member:
    """A gene as one model sees it: the components it fills. ``start``, ``end`` and ``rank`` are
    on the unwrapped axis, so copies shifted around a circular contig differ in them."""
    gene: int
    record_gene: _Gene
    filled: Dict[str, _Match]
    primary: str
    role: str
    start: int
    end: int
    rank: int
    contig_length: int
    cds_count: int

    @property
    def is_cds(self) -> bool:
        return self.record_gene.is_cds

    @property
    def virtual(self) -> bool:
        """A stand-in beyond a contig end for the components a truncated system lacks."""
        return self.record_gene.cds is None

    @property
    def strand(self) -> int:
        return self.record_gene.strand

    @property
    def names(self) -> FrozenSet[str]:
        return frozenset(self.filled)

    def score(self, name: str) -> float:
        return math.inf if self.virtual else self.filled[name].score

    def positions(self, name: str) -> Tuple[int, ...]:
        return self.filled[name].positions

    def shifted(self, k: int) -> "_Member":
        return _Member(self.gene, self.record_gene, self.filled, self.primary, self.role,
                       self.start + k * self.contig_length, self.end + k * self.contig_length,
                       self.rank + k * self.cds_count, self.contig_length, self.cds_count)


@dataclass
class _Evaluation:
    passed: bool
    strand: int
    counted: Dict[str, Tuple[_Member, _Match]]
    """counted component name -> (member, match) with the best score for it"""
    contributors: FrozenSet[str]
    core_count: int
    total_count: int


@dataclass
class _Prepared:
    model: SystemModel
    condition: Optional[Condition]
    extenders: Optional[Condition]
    components: Dict[str, Component]
    positive: FrozenSet[str]
    """Components whose presence counts: core and accessory for quorum models, the condition's positive names otherwise."""


@dataclass
class _Found:
    """A system found on one contig, before overlap resolution."""
    prepared: _Prepared
    members: List[_Member]
    evaluation: _Evaluation
    radius_cds: Optional[float]
    radius_bp: Optional[float]
    truncated: str = ""
    """For a system found only with stand-ins beyond the contig ends: which ends, "<", ">" or "<>" (contig order)."""


def _gene_bounds(location, contig_length: int, circular: bool) -> Tuple[int, int]:
    """(start, end) of a location on the unwrapped axis: end passes the contig length when it wraps the origin."""
    if circular and len(location.parts) > 1:
        envelope = circular_envelope(location, contig_length)
        if len(envelope.parts) == 2:
            first = max(envelope.parts, key=lambda part: int(part.start))
            other = min(envelope.parts, key=lambda part: int(part.start))
            return int(first.start), contig_length + int(other.end)
    return int(location.start), int(location.end)


def _location_strand(location) -> int:
    strand = location.parts[0].strand
    return -1 if strand == -1 else 1


class SystemDetector:
    """Detects the systems of ``collections`` on records, one record at a time.

    evalue: only hits with an e-value at most this count as evidence (in addition to each
        model's own thresholds).
    databases: only hits from these databases count as evidence.
    annotation_type: which hits count: domainator, domain_search or all.
    partial: which CDSs can be members, by fragment status: include, exclude or only.
    max_system_overlap: reject a system when more than this fraction of its counted members
        are already used by better systems (1: keep everything).
    max_systems_per_contig: keep at most this many systems per record, best first.
    truncated_systems: "include" also reports systems that a linear contig end may have cut:
        windows at an end that fail the requirement, but would pass if every component they lack
        lay just beyond the end. They are flagged with ``partial``. "exclude" (the default)
        reports only systems that pass.
    """

    def __init__(self, collections: Sequence[ModelCollection], *, evalue: Optional[float] = None,
                 databases: Optional[Collection[str]] = None, annotation_type: str = "domainator",
                 partial: str = "include", max_system_overlap: float = 1.0,
                 max_systems_per_contig: Optional[int] = None, truncated_systems: str = "exclude"):
        if partial not in PARTIAL_CHOICES:
            raise ValueError(f"partial must be one of {', '.join(PARTIAL_CHOICES)}, not {partial}")
        if truncated_systems not in TRUNCATED_CHOICES:
            raise ValueError(f"truncated_systems must be one of {', '.join(TRUNCATED_CHOICES)}, not {truncated_systems}")
        self.truncated_systems = truncated_systems
        if not 0 <= max_system_overlap <= 1:
            raise ValueError("max_system_overlap must be between 0 and 1.")
        self.collections = list(collections)
        self.evalue = evalue
        self.databases = None if databases is None else set(databases)
        self.feature_types = annotation_feature_types(check_annotation_type(annotation_type))
        self.partial = partial
        self.max_system_overlap = max_system_overlap
        self.max_systems_per_contig = max_systems_per_contig
        self._warned: Set[Tuple[str, str]] = set()
        self._prepared: List[List[_Prepared]] = list()
        self._lookup: List[Dict[str, List[Tuple[int, Component, ReferenceRule]]]] = list()
        self._feature_types: Set[str] = set()
        for collection in self.collections:
            prepared = list()
            lookup: Dict[str, List[Tuple[int, Component, ReferenceRule]]] = dict()
            for model_index, model in enumerate(collection.models):
                define = dict(model.define)
                condition = None if model.condition is None else Condition(model.condition, define)
                extenders = None if model.extenders is None else Condition(model.extenders, define)
                if model.has_quorum:
                    positive = frozenset(c.name for c in model.components if c.role in COUNTED_ROLES)
                else:
                    positive = condition.positive_names
                prepared.append(_Prepared(model, condition, extenders, {c.name: c for c in model.components}, positive))
                for component in model.components:
                    for rule in component.references:
                        lookup.setdefault(rule.name, list()).append((model_index, component, rule))
                        if rule.feature:
                            self._feature_types.add(rule.name)
            self._prepared.append(prepared)
            self._lookup.append(lookup)

    # ---- step 1: order the contig

    def _genes(self, record: SeqRecord) -> List[_Gene]:
        circular = record.annotations.get("topology") == "circular"
        length = len(record)
        cdss = DomainatorCDS.list_from_contig(record, domain_evalue=self.evalue, include_nucleic_acid_annotations=True)
        genes = list()
        for cds in cdss:
            hits = list()
            if DOMAIN_FEATURE_NAME in self.feature_types:
                hits.extend(cds.domain_features)
            if DOMAIN_SEARCH_BEST_HIT_NAME in self.feature_types and cds.domain_search_feature is not None:
                hits.append(cds.domain_search_feature)
            if self.databases is not None:
                hits = [hit for hit in hits if hit.qualifiers["database"][0] in self.databases]
            location = cds.feature.location
            start, end = _gene_bounds(location, length, circular)
            if cds.is_nucleic_acid:
                cds_id = "nuc_" + get_cds_unique_name(SeqFeature(location, type="misc_feature"))
            else:
                cds_id = cds.num
                if self.partial != "include" and not fragment_status_allowed(cds_fragment_status(record, cds.feature), self.partial):
                    hits = []
            genes.append(_Gene(0, cds, start, end, _location_strand(location), not cds.is_nucleic_acid, 0, cds_id,
                               len(location), hits))
        genes.extend(self._feature_genes(record, circular))
        genes.sort(key=lambda gene: (gene.start, gene.end))
        rank = 0
        for index, gene in enumerate(genes):
            gene.index = index
            gene.rank = rank
            if gene.is_cds:
                rank += 1
        return genes

    def _feature_genes(self, record: SeqRecord, circular: bool) -> List[_Gene]:
        """Pseudo-CDSs for the record features that feature components can match, each with a stand-in hit
        (database ``feature``, name the feature type, score 0, e-value 0) carrying the feature's qualifiers."""
        genes = list()
        for feature in record.features:
            if feature.type not in self._feature_types:
                continue
            location = feature.location
            hit = SeqFeature(location, type=_FEATURE_HIT, qualifiers={
                **feature.qualifiers, "database": [FEATURE_DATABASE], "name": [feature.type], "score": ["0"], "evalue": ["0"]})
            cds_id = "nuc_" + get_cds_unique_name(SeqFeature(location, type="misc_feature"))
            cds = DomainatorCDS(feature.type, cds_id, -1, feature, [hit], {feature.type}, None, True)
            start, end = _gene_bounds(location, len(record), circular)
            genes.append(_Gene(0, cds, start, end, _location_strand(location), False, 0, cds_id, len(location), [hit]))
        return genes

    # ---- steps 2-4: competing hits, matching, primary components

    def _competing(self, collection: ModelCollection, gene: _Gene) -> List[SeqFeature]:
        if not collection.competing:
            return gene.hits
        hits = sorted(gene.hits, key=feature_score, reverse=True)
        dropped = set()
        divisor = 3 if gene.is_cds else 1
        for group in collection.competing:
            kept: List[SeqFeature] = list()
            for hit in (h for h in hits if h.qualifiers["name"][0] in group):
                beaten = any(feature_score(other) > feature_score(hit) + _SCORE_TIE
                             and overlap_size(other.location, hit.location) / divisor > collection.competing_overlap_aa
                             for other in kept)
                if beaten:
                    dropped.add(id(hit))
                else:
                    kept.append(hit)
        return [hit for hit in gene.hits if id(hit) not in dropped]

    @staticmethod
    def _rule_passes(rule: ReferenceRule, hit: SeqFeature, gene: _Gene) -> bool:
        qualifiers = hit.qualifiers
        if rule.feature != (hit.type == _FEATURE_HIT):
            return False
        if rule.feature:
            return all(any(str(value).lower() == wanted.lower() for value in qualifiers.get(key, ()))
                       for key, wanted in rule.qualifiers)
        if rule.database is not None and qualifiers["database"][0] != rule.database:
            return False
        if not evalue_passes(hit, rule.evalue):
            return False
        if rule.min_score is not None and feature_score(hit) < rule.min_score:
            return False
        if rule.min_reference_coverage is not None:
            try:
                coverage = (int(qualifiers["rend"][0]) - int(qualifiers["rstart"][0]) + 1) / int(qualifiers["rlen"][0])
            except (KeyError, ValueError, ZeroDivisionError):
                return False
            if coverage < rule.min_reference_coverage:
                return False
        if rule.min_target_coverage is not None:
            coverage = len(hit.location) / gene.length if gene.is_cds else 1.0
            if coverage < rule.min_target_coverage:
                return False
        if rule.min_target_length is not None or rule.max_target_length is not None:
            if not gene.is_cds:
                return False
            translation = gene.cds.feature.qualifiers.get("translation")
            length = len(translation[0]) if translation else gene.length // 3
            if rule.min_target_length is not None and length < rule.min_target_length:
                return False
            if rule.max_target_length is not None and length > rule.max_target_length:
                return False
        for key, qualifier in STRUCTURAL_THRESHOLDS.items():
            threshold = getattr(rule, key)
            if threshold is not None:
                try:
                    if float(qualifiers[qualifier][0]) < threshold:
                        return False
                except (KeyError, IndexError, ValueError):
                    return False
        return True

    @staticmethod
    def _attached_gene(hit: SeqFeature, gene: _Gene, genes: Sequence[_Gene]) -> _Gene:
        """The CDS a nucleotide hit overlaps on its strand by at least half its length, or its own pseudo-CDS."""
        best, best_overlap = gene, len(hit.location) / 2
        for other in genes:
            if not other.is_cds or other.strand != gene.strand:
                continue
            overlap = overlap_size(hit.location, other.cds.feature.location)
            if overlap >= best_overlap:
                best, best_overlap = other, overlap
        return best

    def _matches(self, collection_index: int, genes: Sequence[_Gene]) -> Dict[int, Dict[int, Dict[str, _Match]]]:
        """model index -> gene index -> component name -> the best match."""
        collection = self.collections[collection_index]
        lookup = self._lookup[collection_index]
        out: Dict[int, Dict[int, Dict[str, _Match]]] = dict()
        for gene in genes:
            for hit in self._competing(collection, gene):
                for model_index, component, rule in lookup.get(hit.qualifiers["name"][0], ()):
                    if not self._rule_passes(rule, hit, gene):
                        continue
                    target = gene
                    if rule.attach_to_cds and not gene.is_cds:
                        target = self._attached_gene(hit, gene, genes)
                    position = int(hit.location.start) if target.strand == 1 else -int(hit.location.end)
                    match = _Match(component, hit, feature_score(hit), feature_evalue(hit), (position,))
                    filled = out.setdefault(model_index, dict()).setdefault(target.index, dict())
                    current = filled.get(component.name)
                    if current is None:
                        filled[component.name] = match
                    else:
                        best = match if (match.score, -match.evalue) > (current.score, -current.evalue) else current
                        filled[component.name] = replace(best, positions=tuple(sorted(current.positions + match.positions)))
        return out

    # ---- the requirement

    def _orientation(self, prepared: _Prepared, members: Sequence[_Member]) -> int:
        model = prepared.model
        members = [m for m in members if not m.virtual]
        if model.orientation_from is not None:
            filling = [(m.filled[model.orientation_from].score, m.strand) for m in members if model.orientation_from in m.filled]
            if filling:
                return max(filling)[1]
        counted = [m for m in members if m.is_cds and m.role in COUNTED_ROLES]
        if counted:
            plus = sum(1 for m in counted if m.strand == 1)
            minus = len(counted) - plus
            if plus != minus:
                return 1 if plus > minus else -1
            core = [m for m in counted if m.role == "core"] or counted
            return max(core, key=lambda m: m.filled[m.primary].score).strand
        if members:
            return max(members, key=lambda m: m.filled[m.primary].score).strand
        return 1

    def _evaluate(self, prepared: _Prepared, run: Sequence[_Member], reach: Sequence[_Member]) -> _Evaluation:
        model = prepared.model
        strand = self._orientation(prepared, run)
        contributors: FrozenSet[str] = frozenset()
        condition_ok = True
        if prepared.condition is not None:
            reverse = strand == -1
            ordered_run = sorted(run, key=lambda m: m.start, reverse=reverse)
            ordered_reach = sorted(reach, key=lambda m: m.start, reverse=reverse)
            condition_ok, contributors = prepared.condition.evaluate(ordered_run, ordered_reach)
        counted, core_count, total_count = self._counted(prepared, run, contributors)
        forbidden = any(m.role == "forbidden" for m in reach)
        passed = (condition_ok and not forbidden and core_count >= model.min_core and total_count >= model.min_total
                  and bool(counted))
        return _Evaluation(passed, strand, counted, contributors, core_count, total_count)

    @staticmethod
    def _counted(prepared: _Prepared, run: Sequence[_Member], contributors: FrozenSet[str]):
        """The counted components of ``run`` (with the member and match that fill each best), and the
        numbers of core and of counted components among them."""
        counted: Dict[str, Tuple[_Member, _Match]] = dict()
        for member in run:
            if member.virtual:  # a stand-in fills every component
                names = [name for name in member.filled if prepared.components[name].role in COUNTED_ROLES]
            elif prepared.model.has_quorum:
                names = [member.primary] if member.role in COUNTED_ROLES else []
            else:
                names = contributors & member.names
            for name in names:
                match = member.filled[name]
                current = counted.get(name)
                if current is None or match.score > current[1].score:
                    counted[name] = (member, match)
        core_count = sum(1 for name in counted if prepared.components[name].role == "core")
        total_count = sum(1 for name in counted if prepared.components[name].role in COUNTED_ROLES)
        return counted, core_count, total_count

    def _passing(self, prepared):
        def evaluate(run, reach):
            evaluation = self._evaluate(prepared, run, reach)
            return evaluation if evaluation.passed else None
        return evaluate

    # ---- step 5: windows

    def _members(self, prepared: _Prepared, matches: Dict[int, Dict[str, _Match]], genes: Sequence[_Gene],
                 contig_length: int, cds_count: int) -> List[_Member]:
        members = list()
        for gene_index in sorted(matches):
            filled = matches[gene_index]
            primary = min(filled, key=lambda name: (_ROLE_PRIORITY[filled[name].component.role], -filled[name].score,
                                                    filled[name].evalue))
            gene = genes[gene_index]
            members.append(_Member(gene.index, gene, filled, primary, filled[primary].component.role,
                                   gene.start, gene.end, gene.rank, contig_length, cds_count))
        return members

    def _chain(self, prepared: _Prepared, members: List[_Member], circular: bool) -> List[_Found]:
        window = prepared.model.window

        def gap_limit(a: _Member, b: _Member) -> Optional[int]:
            overrides = [prepared.components[m.primary].max_gap_cds for m in (a, b)
                         if prepared.components[m.primary].max_gap_cds is not None]
            return min(overrides) if overrides else window.max_cds

        def linked(a: _Member, b: _Member) -> bool:
            limit = gap_limit(a, b)
            if limit is not None and windows.cds_between(a, b) > limit:
                return False
            if window.max_kb is not None and windows.gap_bp(a, b) > window.max_kb * 1000:
                return False
            return True

        found = list()
        for cluster in windows.chain(windows.unwrap(members, linked, circular), linked):
            found.extend(self._capped(prepared, cluster, window.max_cds, window.max_kb))
        return found

    def _capped(self, prepared: _Prepared, cluster: List[_Member], radius_cds, radius_kb) -> List[_Found]:
        """A chain or anchor cluster as systems: itself, or, when it is longer than the model's
        max_span caps, the span-mode systems inside it."""
        window = prepared.model.window
        if (window.max_span_cds is not None or window.max_span_kb is not None) and not windows.fits_span(
                cluster, window.max_span_cds, window.max_span_kb):
            return self._span(prepared, cluster, window.max_span_cds, window.max_span_kb, circular=False)
        evaluation = self._evaluate(prepared, cluster, cluster)
        if not evaluation.passed:
            return []
        members = self._trim(cluster)
        return [_Found(prepared, members, self._evaluate(prepared, members, members), radius_cds,
                       None if radius_kb is None else radius_kb * 1000)]

    @staticmethod
    def _trim(members: List[_Member]) -> List[_Member]:
        start, end = 0, len(members)
        while start < end and members[start].role == "neutral":
            start += 1
        while end > start and members[end - 1].role == "neutral":
            end -= 1
        return members[start:end]

    def _span(self, prepared: _Prepared, members: List[_Member], max_cds, max_kb, *, circular: bool) -> List[_Found]:
        if not members:
            return []
        n = len(members)
        if circular:
            items = [m.shifted(-1) for m in members] + members + [m.shifted(1) for m in members]
            starts = range(n, 2 * n)
        else:
            items = members
            starts = None
        evaluate = self._passing(prepared)
        candidates = [(items[i:j + 1], result) for i, j, result in windows.span_candidates(
            items, lambda run: windows.fits_span(run, max_cds, max_kb), evaluate, starts)]
        picked = windows.select_greedily(
            candidates, key=lambda run, result: (-sum(match.score for _, match in result.counted.values()),
                                                 windows.cds_inclusive(run), run[0].start))
        found = list()
        for run, result in picked:
            # a missing member past an end takes a CDS of its own
            radius_cds = None if max_cds is None else max_cds - windows.cds_inclusive(run) - 1
            radius_bp = None if max_kb is None else max_kb * 1000 - windows.span_bp(run)
            found.append(_Found(prepared, list(run), result, radius_cds, radius_bp))
        return found

    def _anchor(self, prepared: _Prepared, members: List[_Member], circular: bool) -> List[_Found]:
        window = prepared.model.window
        n = len(members)
        items = [m.shifted(-1) for m in members] + members + [m.shifted(1) for m in members] if circular else members
        core_genes: Dict[int, _Member] = dict()
        for anchor in members:
            if anchor.role not in COUNTED_ROLES and not (anchor.names & prepared.positive):
                continue
            region = windows.anchor_window(items, anchor, window.max_cds, window.max_kb)
            evaluation = self._evaluate(prepared, region, region)
            if not evaluation.passed:
                continue
            if prepared.model.has_quorum:
                contributes = anchor.role in COUNTED_ROLES
                helpers = [m for m in region if m.role in COUNTED_ROLES]
            else:
                contributes = bool(anchor.names & evaluation.contributors)
                helpers = [m for m in region if m.names & evaluation.contributors]
            if not contributes:
                continue
            for member in [anchor, *helpers]:
                core_genes.setdefault(member.gene, next(m for m in members if m.gene == member.gene))
        core = sorted(core_genes.values(), key=lambda m: m.start)

        def linked(a, b):
            return windows.within_radius(a, b, window.max_cds, window.max_kb)

        found = list()
        for cluster in windows.chain(windows.unwrap(core, linked, circular), linked):
            if (window.max_span_cds is not None or window.max_span_kb is not None) and not windows.fits_span(
                    cluster, window.max_span_cds, window.max_span_kb):
                found.extend(self._span(prepared, cluster, window.max_span_cds, window.max_span_kb, circular=False))
                continue
            found.append(_Found(prepared, cluster, self._evaluate(prepared, cluster, cluster), window.max_cds,
                                None if window.max_kb is None else window.max_kb * 1000))
        return found

    def _windows(self, prepared: _Prepared, members: List[_Member], circular: bool) -> List[_Found]:
        window = prepared.model.window
        if window.mode == "chain":
            return self._chain(prepared, members, circular)
        if window.mode == "anchor":
            return self._anchor(prepared, members, circular)
        return self._span(prepared, members, window.max_cds, window.max_kb, circular=circular)

    def _truncated(self, prepared: _Prepared, members: List[_Member], systems: Sequence[_Found], contig_length: int,
                   cds_count: int) -> List[_Found]:
        """The systems that fail, but pass with a stand-in member just beyond a contig end that fills
        every component, and share no member with ``systems``."""
        names = prepared.positive

        def stand_in(gene, start, end, rank):
            record_gene = _Gene(gene, None, start, end, 0, True, rank, "", 0, [])
            filled = {name: _Match(prepared.components[name], None, -math.inf, math.inf) for name in names}
            primary = min(filled, key=lambda name: (_ROLE_PRIORITY[prepared.components[name].role], name))
            return _Member(gene, record_gene, filled, primary, prepared.components[primary].role, start, end, rank,
                           contig_length, cds_count)

        left, right = stand_in(-1, -1, 0, -1), stand_in(-2, contig_length, contig_length + 1, cds_count)
        window = prepared.model.window
        if window.mode == "span":
            # the window at each end: a stand-in passes on its own, so the tightest runs would be the stand-ins
            candidates = list()
            for end, ordered in ((left, members), (right, members[::-1])):
                run = [end]
                for member in ordered:
                    if not windows.fits_span(sorted(run + [member], key=lambda m: m.start), window.max_cds, window.max_kb):
                        break
                    run.append(member)
                run.sort(key=lambda m: m.start)
                evaluation = self._evaluate(prepared, run, run)
                if evaluation.passed:
                    needed = [i for i, m in enumerate(run) if not m.virtual and self._counted(prepared, [m], evaluation.contributors)[0]]
                    if needed:
                        first, last = (0, max(needed)) if end is left else (min(needed), len(run) - 1)
                        candidates.append(_Found(prepared, run[first:last + 1], evaluation, None, None))
        else:
            candidates = self._windows(prepared, [left, *members, right], circular=False)
        used = {m.gene for system in systems for m in system.members}
        out = list()
        for system in candidates:
            ends = {m.gene for m in system.members if m.virtual}
            real = self._trim([m for m in system.members if not m.virtual])
            if not ends or not real or used & {m.gene for m in real}:
                continue
            evaluation = self._evaluate(prepared, real, real)
            if evaluation.passed:
                continue
            # count what the stand-ins helped satisfy, as far as the real members fill it
            counted, core_count, total_count = self._counted(prepared, real, system.evaluation.contributors)
            evaluation = replace(evaluation, counted=counted, contributors=system.evaluation.contributors,
                                 core_count=core_count, total_count=total_count)
            used |= {m.gene for m in real}
            if window.mode == "span":
                radius_cds = None if window.max_cds is None else window.max_cds - windows.cds_inclusive(real) - 1
                radius_bp = None if window.max_kb is None else window.max_kb * 1000 - windows.span_bp(real)
            else:
                radius_cds, radius_bp = system.radius_cds, system.radius_bp
            out.append(replace(system, members=real, evaluation=evaluation, radius_cds=radius_cds, radius_bp=radius_bp,
                               truncated=("<" if -1 in ends else "") + (">" if -2 in ends else "")))
        return out

    # ---- step 7: extend

    @staticmethod
    def _extend(found: _Found, members: Sequence[_Member], circular: bool) -> _Found:
        """``found`` grown outward over the members next to it that meet the model's extenders, each
        within the window's radius of the last member on its side."""
        prepared = found.prepared
        window = prepared.model.window
        if not found.members:
            return found
        pool = [m.shifted(k) for k in (-1, 0, 1) for m in members] if circular else list(members)
        system = sorted(found.members, key=lambda m: m.start)
        used = {m.gene for m in system}
        extended: List[_Member] = list()
        sides = ((system[-1], sorted((m for m in pool if m.start > system[-1].start), key=lambda m: m.start)),
                 (system[0], sorted((m for m in pool if m.start < system[0].start), key=lambda m: m.start, reverse=True)))
        for edge, candidates in sides:
            for candidate in candidates:
                if not windows.within_radius(edge, candidate, window.max_cds, window.max_kb):
                    break
                if candidate.gene in used:
                    continue
                passed, contributors = prepared.extenders.evaluate([candidate])
                filled = {name: candidate.filled[name] for name in contributors & candidate.names}
                if not passed or not filled:
                    continue
                primary = min(filled, key=lambda name: (-filled[name].score, filled[name].evalue, name))
                extended.append(replace(candidate, filled=filled, primary=primary, role="extender"))
                used.add(candidate.gene)
                edge = candidate
        if not extended:
            return found
        grown = sorted(system + extended, key=lambda m: m.start)
        if grown[0].start < 0:  # extended back across the origin: keep starts on the contig
            grown = [m.shifted(1) for m in grown]
        return replace(found, members=grown)

    # ---- detection

    def detect(self, record: SeqRecord) -> List[SystemHit]:
        """The systems on ``record``, in contig order. Protein records have none."""
        if record.annotations.get("molecule_type") == "protein":
            return []
        circular = record.annotations.get("topology") == "circular"
        contig_length = len(record)
        genes = self._genes(record)
        cds_count = sum(1 for gene in genes if gene.is_cds)
        found: List[_Found] = list()
        for collection_index, collection in enumerate(self.collections):
            matches = self._matches(collection_index, genes)
            for model_index, prepared in enumerate(self._prepared[collection_index]):
                model_matches = matches.get(model_index)
                if not model_matches:
                    continue
                if not any(set(filled) & prepared.positive for filled in model_matches.values()):
                    continue
                window = prepared.model.window
                if cds_count == 0 and window.max_kb is None:
                    if prepared.model.key not in self._warned:
                        self._warned.add(prepared.model.key)
                        warnings.warn(f"System model {prepared.model.database}:{prepared.model.name} measures its window only "
                                      "in CDSs, so it is skipped on records without CDSs.")
                    continue
                members = self._members(prepared, model_matches, genes, contig_length, cds_count)
                strands = [[m for m in members if m.strand == s] for s in (1, -1)] if prepared.model.same_strand else [members]
                for group in strands:
                    systems = self._windows(prepared, group, circular)
                    if self.truncated_systems == "include" and not circular and group:
                        systems.extend(self._truncated(prepared, group, systems, contig_length, cds_count))
                    if prepared.extenders is not None:
                        systems = [self._extend(system, group, circular) for system in systems]
                    found.extend(systems)
        found = [f for f in found if f.members and f.evaluation.counted]
        hits = [self._hit(f, record, genes, circular, cds_count) for f in found]
        accepted = self._resolve(hits, found)
        accepted.sort(key=lambda i: (min(m.start for m in found[i].members), hits[i].model.database, hits[i].model.name))
        return [hits[i] for i in accepted]

    # ---- steps 7-9: build, resolve

    def _hit(self, found: _Found, record: SeqRecord, genes: Sequence[_Gene], circular: bool, cds_count: int) -> SystemHit:
        prepared, evaluation = found.prepared, found.evaluation
        model = prepared.model
        members = sorted(found.members, key=lambda m: m.start)
        system_members = list()
        for member in members:
            if model.has_quorum:
                names = [member.primary] if member.role != "forbidden" else []
            else:
                names = [name for name in member.filled if prepared.components[name].role != "forbidden"]
                names.sort(key=lambda name: (_ROLE_PRIORITY[prepared.components[name].role], name))
            if member.role == "extender":
                names = sorted(member.filled)
            for name in names:
                match = member.filled[name]
                component = prepared.components[name]
                role = "extender" if member.role == "extender" else component.role
                system_members.append(SystemMember(member.record_gene.cds, member.record_gene.cds_id, component,
                                                   role, match.hit))
        contig_length = len(record)
        strand = evaluation.strand
        start = min(m.start for m in members)
        end = max(m.end for m in members)
        if end - start >= contig_length:
            location = FeatureLocation(0, contig_length, strand)
        elif circular and end > contig_length:
            after = FeatureLocation(start, contig_length, strand)
            before = FeatureLocation(0, end - contig_length, strand)
            location = CompoundLocation([after, before] if strand == 1 else [before, after])
        else:
            location = FeatureLocation(start, end, strand)
        counted_matches = [match for _, match in evaluation.counted.values()]
        score = sum(match.score for match in counted_matches)
        evalue = max(match.evalue for match in counted_matches)
        completeness = None
        if model.has_quorum:
            total = model.max_components or sum(1 for c in model.components if c.role in COUNTED_ROLES)
            completeness = evaluation.total_count / total if total else None
        partial = self._partial(found, members, record, genes, circular, cds_count)
        return SystemHit(model, tuple(system_members), location, score, evalue, completeness,
                         evaluation.core_count, evaluation.total_count, partial)

    @staticmethod
    def _partial(found: _Found, members: Sequence[_Member], record, genes, circular, cds_count) -> str:
        if circular:
            return ""
        found_names = set().union(*(m.names for m in members))
        if not (found.prepared.positive - found_names):
            return ""
        first, last = members[0], members[-1]
        cds_left = first.rank
        cds_right = cds_count - (last.rank + (1 if last.is_cds else 0))
        bp_left = first.start
        bp_right = len(record) - max(m.end for m in members)

        def reachable(cds, bp):
            return ((found.radius_cds is None or cds <= found.radius_cds)
                    and (found.radius_bp is None or bp <= found.radius_bp))

        left = reachable(cds_left, bp_left) or "<" in found.truncated
        right = reachable(cds_right, bp_right) or ">" in found.truncated
        if found.evaluation.strand == -1:
            left, right = right, left
        return ("<" if left else "") + (">" if right else "")

    def _resolve(self, hits: List[SystemHit], found: List[_Found]) -> List[int]:
        """The indices of the hits kept by superiors, greedy overlap and the per-contig cap."""
        genes_of = [frozenset(m.gene for m in f.members) for f in found]
        by_key: Dict[Tuple[str, str], List[int]] = dict()
        for i, hit in enumerate(hits):
            by_key.setdefault(hit.model.key, list()).append(i)
        # superiors
        keep = list()
        for i, hit in enumerate(hits):
            dropped = False
            for superior in hit.model.superiors:
                for j in by_key.get(superior, ()):
                    if genes_of[i] & genes_of[j] or self._contains(hits[j].location, hit.location):
                        dropped = True
                        break
                if dropped:
                    break
            if not dropped:
                keep.append(i)
        # greedy overlap
        order = sorted(keep, key=lambda i: (-hits[i].score, hits[i].model.database, hits[i].model.name, min(genes_of[i])))
        used: Dict[int, List[Tuple[str, str]]] = dict()
        accepted = list()
        for i in order:
            hit = hits[i]
            counted_genes = {m.gene: m for m in found[i].members if m.role in COUNTED_ROLES}
            if counted_genes and self.max_system_overlap < 1:
                conflicts = 0
                for gene, member in counted_genes.items():
                    component = found[i].prepared.components[member.primary]
                    for key in used.get(gene, ()):
                        if (key != hit.model.key and not component.multi_model) or (key == hit.model.key and not component.multi_system):
                            conflicts += 1
                            break
                if conflicts / len(counted_genes) > self.max_system_overlap:
                    continue
            for gene in counted_genes:
                used.setdefault(gene, list()).append(hit.model.key)
            accepted.append(i)
        if self.max_systems_per_contig is not None:
            accepted = accepted[:self.max_systems_per_contig]
        return accepted

    @staticmethod
    def _contains(outer, inner) -> bool:
        outer_parts = [(int(p.start), int(p.end)) for p in outer.parts]
        return all(any(s <= int(p.start) and int(p.end) <= e for s, e in outer_parts) for p in inner.parts)


def detect_systems(record: SeqRecord, collections: Sequence[ModelCollection], **options) -> List[SystemHit]:
    """The systems of ``collections`` on ``record``. ``options`` are those of SystemDetector."""
    return SystemDetector(collections, **options).detect(record)
