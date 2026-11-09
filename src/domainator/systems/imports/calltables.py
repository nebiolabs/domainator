"""System calls from DefenseFinder's and PADLOC's output tables, as Domain_System features.

Checked against DefenseFinder's test outputs (defense-finder's tests/data) and padloc's
bin/padloc.R. Read:

- DefenseFinder: ``<name>_defense_finder_genes.tsv``, one row per member: ``replicon``,
  ``hit_id`` (the protein), ``hit_pos`` (its position on the replicon, 1-based), ``gene_name``
  (the profile that hit), ``hit_gene_ref`` (the model gene it fills), ``hit_status`` (mandatory,
  accessory, ...), ``hit_i_eval``, ``hit_score``, ``sys_id``, ``model_fqn``, ``sys_score``,
  ``sys_wholeness``, ``type`` and ``subtype``.
- PADLOC: ``<name>_padloc.csv``, one row per member: ``system.number``, ``seqid``, ``system``
  (the model's yaml.name), ``target.name`` (the protein), ``hmm.name``, ``protein.name`` (the
  model gene it fills), ``domain.iE.value``, and the gene's ``start``, ``end`` (GFF, 1-based) and
  ``strand``. It doesn't say which genes are core, so every member is imported as core.

Members are matched to the record's CDSs: PADLOC's by coordinates, then by name; DefenseFinder's
by name (protein_id, locus_tag, gene_id, gene or cds_id), then, for DefenseFinder's own gene
calls (``<replicon>_<n>``), as the record's n-th CDS. Members with no CDS (PADLOC's CRISPR
arrays and ncRNAs) are skipped.
"""
import csv
from dataclasses import dataclass, field, replace
import math
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
import warnings

from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from domainator import SYSTEM_FEATURE_NAME
from domainator.errors import InputError
from domainator.location_ops import stranded_start
from domainator.model import DomainatorCDS, get_cds_unique_name
from domainator.systems.features import system_feature
from domainator.systems.model import Component, ReferenceRule, SystemHit, SystemMember, SystemModel, Window

DEFENSEFINDER_DATABASE = "defensefinder"
PADLOC_DATABASE = "padloc"
_DEFENSEFINDER_COLUMNS = {"replicon", "hit_id", "hit_pos", "gene_name", "hit_gene_ref", "hit_status", "sys_id",
                          "model_fqn", "subtype"}
_PADLOC_COLUMNS = {"system.number", "seqid", "system", "target.name", "hmm.name", "protein.name", "start", "end", "strand"}
_ROLES = {"mandatory": "core", "accessory": "accessory", "neutral": "neutral", "forbidden": "forbidden"}
_NAME_QUALIFIERS = ("protein_id", "locus_tag", "gene_id", "gene", "cds_id")


@dataclass
class CallMember:
    protein: str
    component: str
    role: str
    reference: str
    evalue: float
    score: float
    position: Optional[int] = None
    """DefenseFinder: the gene's position on its replicon (1-based)."""
    coordinates: Optional[Tuple[int, int, int]] = None
    """PADLOC: (start, end, strand) of the gene, 0-based half-open."""


@dataclass
class Call:
    """One system another tool called on one contig."""
    contig: str
    model: SystemModel
    score: float
    completeness: Optional[float]
    members: List[CallMember] = field(default_factory=list)


def _number(value: Optional[str], default: float = math.nan) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _rows(path: Path, delimiter: str, required: set, label: str) -> List[Dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise InputError(f"{path} is not {label}: it has no {', '.join(sorted(missing))} column.")
        return list(reader)


def read_defensefinder(paths: Iterable[str]) -> List[Call]:
    """The systems in DefenseFinder ``*_defense_finder_genes.tsv`` tables."""
    calls: Dict[Tuple[str, str], Call] = dict()
    for path in paths:
        for row in _rows(Path(path), "\t", _DEFENSEFINDER_COLUMNS,
                         "a DefenseFinder genes table (*_defense_finder_genes.tsv)"):
            key = (row["replicon"], row["sys_id"])
            if key not in calls:
                accession = row["model_fqn"].split("/", 1)[-1]  # without the package name
                parts = accession.split("/")
                model = SystemModel(name=row["subtype"], database=DEFENSEFINDER_DATABASE, program="macsyfinder",
                                    components=(), window=Window("chain"), accession=accession,
                                    category=parts[0] if len(parts) > 1 else ".", description=row.get("type") or ".")
                calls[key] = Call(row["replicon"], model, _number(row.get("sys_score"), 0.0),
                                  _number(row.get("sys_wholeness"), math.nan))
            calls[key].members.append(CallMember(
                protein=row["hit_id"], component=row["hit_gene_ref"].split("__", 1)[-1],
                role=_ROLES.get(row["hit_status"], "core"), reference=row["gene_name"],
                evalue=_number(row.get("hit_i_eval")), score=_number(row.get("hit_score"), 0.0),
                position=int(row["hit_pos"]) if row["hit_pos"].isdigit() else None))
    for call in calls.values():
        if math.isnan(call.completeness):
            call.completeness = None
    return list(calls.values())


def read_padloc(paths: Iterable[str]) -> List[Call]:
    """The systems in PADLOC ``*_padloc.csv`` tables."""
    calls: Dict[Tuple[str, str], Call] = dict()
    for path in paths:
        for row in _rows(Path(path), ",", _PADLOC_COLUMNS, "a PADLOC table (*_padloc.csv)"):
            key = (row["seqid"], row["system.number"])
            if key not in calls:
                model = SystemModel(name=row["system"], database=PADLOC_DATABASE, program="padloc", components=(),
                                    window=Window("chain"), accession=row["system"], category="defense")
                calls[key] = Call(row["seqid"], model, 0.0, None)
            try:
                start, end = int(float(row["start"])) - 1, int(float(row["end"]))
            except ValueError:
                coordinates = None
            else:
                coordinates = (start, end, -1 if row["strand"] == "-" else 1)
            calls[key].members.append(CallMember(
                protein=row["target.name"], component=row["protein.name"], role="core", reference=row["hmm.name"],
                evalue=_number(row.get("domain.iE.value")), score=0.0, coordinates=coordinates))
    return list(calls.values())


def _names(cds: DomainatorCDS) -> List[str]:
    return [value for key in _NAME_QUALIFIERS for value in cds.feature.qualifiers.get(key, [])]


def _envelope(locations: Sequence, contig_length: int, circular: bool, strand: int):
    """The smallest location holding every member: across the origin of a circular contig
    when the biggest gap between members is not the one across the origin."""
    spans = sorted((int(location.start), int(location.end)) for location in locations)
    start, end = spans[0][0], max(span_end for _, span_end in spans)
    if circular and len(spans) > 1:
        gaps = [(spans[i + 1][0] - spans[i][1], i) for i in range(len(spans) - 1)]
        widest, index = max(gaps)
        if widest > contig_length - end + start:  # skipping that gap is shorter than going around
            after = FeatureLocation(spans[index + 1][0], contig_length, strand)
            before = FeatureLocation(0, spans[index][1], strand)
            return CompoundLocation([after, before] if strand == 1 else [before, after])
    return FeatureLocation(start, end, strand)


def import_calls(record: SeqRecord, calls: Sequence[Call], database: str) -> SeqRecord:
    """Add a Domain_System feature to record for each of its calls (calls on other contigs are
    ignored), replacing the systems of database already on it."""
    for feature in record.features:
        if feature.type == "CDS" and "cds_id" not in feature.qualifiers:
            feature.qualifiers["cds_id"] = [get_cds_unique_name(feature)]
    cdss = DomainatorCDS.list_from_contig(record)
    cdss.sort(key=lambda cds: stranded_start(cds.feature.location))
    by_name: Dict[str, DomainatorCDS] = dict()
    for cds in cdss:
        for name in _names(cds):
            by_name.setdefault(name, cds)
    by_coordinates = {(int(cds.feature.location.start), int(cds.feature.location.end),
                       -1 if cds.feature.location.parts[0].strand == -1 else 1): cds for cds in cdss}
    circular = record.annotations.get("topology") == "circular"

    def find(member: CallMember, contig: str) -> Optional[DomainatorCDS]:
        if member.coordinates is not None and member.coordinates in by_coordinates:
            return by_coordinates[member.coordinates]
        if member.protein in by_name:
            return by_name[member.protein]
        if member.position is not None and member.protein == f"{contig}_{member.position}" and member.position <= len(cdss):
            return cdss[member.position - 1]  # DefenseFinder's own gene calls, in contig order
        return None

    systems = list()
    unmatched = 0
    for call in calls:
        if call.contig not in (record.id, record.name):
            continue
        components: Dict[str, Component] = dict()
        members = list()
        for member in call.members:
            cds = find(member, call.contig)
            if cds is None:
                unmatched += 1
                continue
            if member.component not in components:
                components[member.component] = Component(member.component, member.role,
                                                         (ReferenceRule(member.reference, database=database),))
            hit = SeqFeature(cds.feature.location, type="system_call",
                             qualifiers={"database": [database], "name": [member.reference]})
            members.append((cds, member, SystemMember(cds, cds.num, components[member.component], member.role, hit)))
        if not members:
            continue
        members.sort(key=lambda item: int(item[0].feature.location.start))
        strands = [-1 if cds.feature.location.parts[0].strand == -1 else 1 for cds, _, _ in members]
        strand = 1 if strands.count(1) >= strands.count(-1) else -1
        location = _envelope([cds.feature.location for cds, _, _ in members], len(record), circular, strand)
        counted = {member.component for _, member, _ in members if member.role in ("core", "accessory")}
        core = {member.component for _, member, _ in members if member.role == "core"}
        evalues = [member.evalue for _, member, _ in members if not math.isnan(member.evalue)]
        model = replace(call.model, components=tuple(components.values()))
        systems.append(system_feature(SystemHit(
            model=model, members=tuple(system_member for _, _, system_member in members), location=location,
            score=call.score, evalue=max(evalues) if evalues else math.nan, completeness=call.completeness,
            core_count=len(core), total_count=len(counted))))
    if unmatched:
        warnings.warn(f"{record.id}: {unmatched} system members match no CDS of the record, and were skipped "
                      f"(CRISPR arrays and ncRNAs never do).")
    record.features = [feature for feature in record.features
                       if not (feature.type == SYSTEM_FEATURE_NAME and feature.qualifiers.get("database", [None])[0] == database)]
    record.features.extend(systems)
    return record
