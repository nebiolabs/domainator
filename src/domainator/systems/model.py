"""The internal model of system definitions, and of detected systems.

Frozen dataclasses, so model collections can be shared read-only between threads and
processes. See docs/design/systems/definition_formats.md.
"""
from dataclasses import dataclass, field
from typing import FrozenSet, Optional, Tuple

from Bio.SeqFeature import SeqFeature

from domainator.model import DomainatorCDS

ROLES = ("core", "accessory", "neutral", "forbidden")
"""Component roles, in the priority order used to choose a CDS's primary component."""
COUNTED_ROLES = frozenset({"core", "accessory"})
"""Roles that count toward a quorum, the score and completeness."""
WINDOW_MODES = ("span", "chain", "anchor")
"""How "close together" is measured. The first is the default."""
COMMAND_LINE_DATABASE = "command_line"
"""The database of a model defined with --domain_expr."""


@dataclass(frozen=True)
class ReferenceRule:
    """One way to satisfy a component: a hit to this reference that passes these thresholds."""
    name: str
    """The hit's ``name`` qualifier (HMM NAME, sequence id, structure name)."""
    database: Optional[str] = None
    """The hit's ``database`` qualifier; None matches any database."""
    evalue: Optional[float] = None
    min_score: Optional[float] = None
    min_reference_coverage: Optional[float] = None
    min_target_coverage: Optional[float] = None
    min_target_length: Optional[int] = None
    max_target_length: Optional[int] = None
    min_tmscore: Optional[float] = None
    """Structural hits (foldseek): the least TM-score, from the hit's tmscore qualifier."""
    min_lddt: Optional[float] = None
    """Structural hits: the least lDDT, from the lddt qualifier."""
    min_prob: Optional[float] = None
    """Structural hits: the least probability of homology, from the prob qualifier."""
    exchangeable: bool = False
    attach_to_cds: bool = False


THRESHOLD_KEYS = ("evalue", "min_score", "min_reference_coverage", "min_target_coverage",
                  "min_target_length", "max_target_length", "min_tmscore", "min_lddt", "min_prob")
"""The ReferenceRule fields that are thresholds, which a collection's ``defaults`` can set."""
STRUCTURAL_THRESHOLDS = {"min_tmscore": "tmscore", "min_lddt": "lddt", "min_prob": "prob"}
"""The structural thresholds, and the hit qualifiers they read. A hit without the qualifier fails them."""


@dataclass(frozen=True)
class Component:
    name: str
    role: str
    references: Tuple[ReferenceRule, ...]
    loner: bool = False
    multi_system: bool = False
    multi_model: bool = False
    max_gap_cds: Optional[int] = None


@dataclass(frozen=True)
class Window:
    mode: str = "span"
    max_cds: Optional[int] = None
    max_kb: Optional[float] = None
    max_span_cds: Optional[int] = None
    max_span_kb: Optional[float] = None

    def describe(self) -> str:
        """The window as written in the ``window`` qualifier, e.g. ``chain:cds=5`` or ``span:kb=2.5``."""
        limits = [(key, value) for key, value in (("cds", self.max_cds), ("kb", self.max_kb),
                                                  ("span_cds", self.max_span_cds), ("span_kb", self.max_span_kb))
                  if value is not None]
        return f"{self.mode}:" + ",".join(f"{key}={value:g}" for key, value in limits)


@dataclass(frozen=True)
class Example:
    """A known occurrence of a system."""
    accession: str
    start: Optional[int] = None
    end: Optional[int] = None
    note: Optional[str] = None


@dataclass(frozen=True)
class SystemModel:
    name: str
    database: str
    program: str
    components: Tuple[Component, ...]
    window: Window
    min_core: int = 0
    min_total: int = 0
    condition: Optional[str] = None
    extenders: Optional[str] = None
    """A condition, evaluated on one CDS at a time, that CDSs next to the system must meet to extend it (format 1.1)."""
    define: Tuple[Tuple[str, str], ...] = ()
    """The aliases (name, expression) the condition and extenders use (format 1.1)."""
    superiors: Tuple[Tuple[str, str], ...] = ()
    """(database, name) of the models that win when they overlap this one."""
    category: str = "."
    neighborhood_kb: Optional[float] = None
    neighborhood_cds: Optional[int] = None
    description: str = "."
    accession: str = "."
    same_strand: bool = False
    max_components: Optional[int] = None
    multi_loci: bool = False
    orientation_from: Optional[str] = None
    examples: Tuple[Example, ...] = ()
    related: Tuple[str, ...] = ()
    """References that are relevant but not required (antiSMASH RELATED). Informational (format 1.1)."""
    literature: Tuple[str, ...] = ()
    reference_files: Tuple[str, ...] = ()

    @property
    def key(self) -> Tuple[str, str]:
        return (self.database, self.name)

    @property
    def has_quorum(self) -> bool:
        return self.min_core > 0 or self.min_total > 0


@dataclass(frozen=True)
class Source:
    """Where a converted collection came from. Informational."""
    name: Optional[str] = None
    url: Optional[str] = None
    release: Optional[str] = None
    date: Optional[str] = None
    converter: Optional[str] = None
    license: Optional[str] = None


@dataclass(frozen=True)
class ModelCollection:
    """What one --system_models input yields."""
    database: str
    models: Tuple[SystemModel, ...]
    semantics: str = "domainator"
    source: Optional[Source] = None
    competing: Tuple[FrozenSet[str], ...] = ()
    competing_overlap_aa: int = 20
    reference_files: Tuple[str, ...] = ()
    path: Optional[str] = field(default=None, compare=False)
    """The file the collection was read from; relative reference_files resolve against its directory."""


@dataclass(frozen=True)
class SystemMember:
    """One (CDS, component) pair of a detected system."""
    cds: DomainatorCDS
    cds_id: str
    """The CDS's cds_id, or ``nuc_<start>_<strand>_<end>`` for a contig-level nucleotide hit."""
    component: Component
    role: str
    """core, accessory, neutral, or extender for a CDS added by the model's extenders."""
    hit: SeqFeature
    """The Domainator feature that satisfied the component."""


@dataclass(frozen=True)
class SystemHit:
    model: SystemModel
    members: Tuple[SystemMember, ...]
    """In contig order."""
    location: object
    """The envelope of the member CDSs (a FeatureLocation, or a CompoundLocation across the origin)."""
    score: float
    evalue: float
    completeness: Optional[float]
    """None for models without a quorum."""
    core_count: int
    total_count: int
    partial: str = ""
    """"", "<", ">" or "<>": which ends (in system orientation) a missing component could lie beyond."""

    @property
    def member_count(self) -> int:
        return len({member.cds_id for member in self.members})
