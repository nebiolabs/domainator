"""Extract contig regions around systems already annotated in GenBank records

Takes GenBank records with systems annotated by find_systems (Domain_System features) and
extracts a region around each selected system, the way select_by_cds extracts regions around
CDSs. Padding (--cds_range, --kb_range, ...) is counted from the system's outermost members, in
system orientation, and --model_neighborhood pads each system by its model's own neighborhood.

Each region gets a System_Search feature marking the system it was selected for, so it can be
handled like system_search output.

--systems/--systems_file and --system_categories select systems (both must match when both
are given). --invert negates that selection only. The other criteria (--databases, --min_score,
--min_completeness, --partial_systems, --strand and --min_context_*) restrict which systems can
be selected, and are not inverted.

--min_context_* skips systems without enough context beyond them, for example systems cut by
the end of a linear contig: --cds_range 10 --min_context_cds 10 selects only systems with 10
CDSs on each side.
"""
from typing import Iterable, Iterator, Optional, Set, Tuple
import warnings

from Bio.SeqRecord import SeqRecord
from jsonargparse import ActionConfigFile
from jsonargparse.typing import NonNegativeFloat, NonNegativeInt

from domainator import SYSTEM_SEARCH_BEST_HIT_NAME
from domainator.cli.common import (add_dedup, add_input, add_list_and_file, add_min_context, add_model_neighborhood,
                                   add_seq_outputs, add_system_selection, add_system_type, check_outputs,
                                   declare_outputs, list_and_file_to_dict_keys, tool_parser)
from domainator.errors import InputError
from domainator.location_ops import overlap_size, stranded_start
from domainator.locations import WHERE_CHOICES
from domainator.model import SYSTEM_TYPES, DomainatorCDS, feature_score, set_hit_best_score
from domainator.outputs import SeqOutputs
from domainator.partial import PARTIAL_CHOICES
from domainator.select_by_cds import (deduplicate_and_pad, get_neighborhood, has_context, neighborhood_feature_index,
                                      up_down_pair, validate_range)
from domainator.seqio import parse_seqfiles
from domainator.systems.records import (as_search_hit, member_indices, system_completeness, system_features,
                                        system_strand)


def _model_neighborhood(feature) -> Tuple[Optional[Tuple[float, float]], Optional[Tuple[int, int]]]:
    """(kb_range, cds_range) from a system's neighborhood_kb, else neighborhood_cds, qualifier."""
    if "neighborhood_kb" in feature.qualifiers:
        kb = float(feature.qualifiers["neighborhood_kb"][0])
        return (kb, kb), None
    if "neighborhood_cds" in feature.qualifiers:
        cds = int(feature.qualifiers["neighborhood_cds"][0])
        return None, (cds, cds)
    return None, None


def _partial_allowed(feature, partial_systems: str) -> bool:
    is_partial = bool(feature.qualifiers.get("system_partial", [""])[0])
    if partial_systems == "exclude":
        return not is_partial
    if partial_systems == "only":
        return is_partial
    return True


def select_by_system(records: Iterable[SeqRecord], *, systems: Optional[Set[str]] = None,
                     system_categories: Optional[Set[str]] = None, databases: Optional[Set[str]] = None,
                     system_type: str = "domain_system", min_score: Optional[float] = None,
                     min_completeness: Optional[float] = None, partial_systems: str = "include",
                     invert: bool = False, strand: Optional[str] = None,
                     cds_range: Optional[Tuple[int, int]] = None, kb_range: Optional[Tuple[float, float]] = None,
                     whole_contig: bool = False, model_neighborhood: bool = False,
                     min_context_cds: Optional[Tuple[int, int]] = None,
                     min_context_kb: Optional[Tuple[float, float]] = None,
                     normalize_direction: bool = True, max_region_overlap: float = 1.0) -> Iterator[SeqRecord]:
    """Yield a region of ``records`` around each selected system.

    A system is selected if it matches ``systems`` (names, or ``database:name``) and
    ``system_categories``, or, with neither, always; ``invert`` negates that. It must also be
    a ``system_type`` system from one of ``databases`` (None: any), score at least
    ``min_score``, have a completeness of at least ``min_completeness`` (systems without a
    quorum have none and pass), pass ``partial_systems`` ("include", "exclude" or "only", by
    the system_partial qualifier), lie on ``strand`` ("f" or "r"; None: both), and have the
    context that ``min_context_cds`` and ``min_context_kb`` require.

    Args:
        records: the contigs (SeqRecords), with Domain_System or System_Search features.
        cds_range, kb_range: (upstream, downstream) CDSs or kb to extract beyond the system's
            outermost members, in system orientation.
        whole_contig: extract the whole contig around each system.
        model_neighborhood: pad each system by its model's neighborhood_kb, else
            neighborhood_cds, instead of by cds_range or kb_range.
        min_context_cds, min_context_kb: (upstream, downstream) CDSs or kb the contig must
            have beyond the system, in system orientation. On a circular contig, the CDSs or kb
            outside the system must add up to at least upstream + downstream. CDSs are counted
            as cds_range counts them. Independent of the padding.
        normalize_direction: reverse-complement regions whose system is on the reverse strand.
        max_region_overlap: the largest fraction of overlap allowed between two output regions
            (1 or more: no filter).

    Yields:
        SeqRecords of the regions, each with a System_Search feature for its system.
    """
    if system_type not in SYSTEM_TYPES:
        raise ValueError(f"system_type must be one of {', '.join(SYSTEM_TYPES)}, not {system_type}")
    if partial_systems not in PARTIAL_CHOICES:
        raise ValueError(f"partial_systems must be one of {', '.join(PARTIAL_CHOICES)}, not {partial_systems}")
    if strand not in (None, "f", "r"):
        raise ValueError(f"strand must be f, r or None, not {strand}")
    if model_neighborhood and (cds_range is not None or kb_range is not None or whole_contig):
        raise InputError("model_neighborhood can't be combined with cds_range, kb_range or whole_contig.")
    systems = set(systems) if systems is not None else None
    system_categories = set(system_categories) if system_categories is not None else None
    databases = set(databases) if databases is not None else None
    if systems is None and system_categories is None and invert:
        warnings.warn("Invert supplied with no selection criteria. Invert will be ignored.")
        invert = False

    def selected(feature) -> bool:
        name = feature.qualifiers["name"][0]
        database = feature.qualifiers["database"][0]
        keep = ((systems is None or name in systems or f"{database}:{name}" in systems)
                and (system_categories is None or feature.qualifiers.get("category", ["."])[0] in system_categories))
        if invert:
            keep = not keep
        if databases is not None and database not in databases:
            return False
        if min_score is not None and feature_score(feature) < min_score:
            return False
        completeness = system_completeness(feature)
        if min_completeness is not None and completeness is not None and completeness < min_completeness:
            return False
        if strand is not None and (system_strand(feature) == -1) != (strand == "r"):
            return False
        return keep and _partial_allowed(feature, partial_systems)

    for record in records:
        focuses = [feature for feature in system_features(record, system_type) if selected(feature)]
        if not focuses:
            continue
        cds_list = DomainatorCDS.list_from_contig(record, include_nucleic_acid_annotations=True)
        cds_list.sort(key=lambda cds: stranded_start(cds.feature.location))
        feature_index = neighborhood_feature_index(record, cds_list)
        selected_locations = list()
        for feature in focuses:
            members = member_indices(feature, cds_list)
            if not members:
                warnings.warn(f"Skipping system {feature.qualifiers['name'][0]} on {record.id}: none of its member CDSs "
                              f"are in the record.")
                continue
            orientation = system_strand(feature)
            if not has_context(record, cds_list, members, feature.location, orientation, min_context_cds, min_context_kb):
                continue
            system_kb_range, system_cds_range = _model_neighborhood(feature) if model_neighborhood else (kb_range, cds_range)
            region, location = get_neighborhood(
                record, cds_list, members, feature.location, orientation, cds_range=system_cds_range,
                kb_range=system_kb_range, whole_contig=whole_contig, normalize_direction=normalize_direction,
                keep_feature=lambda f: f.type != SYSTEM_SEARCH_BEST_HIT_NAME, extra_features=[as_search_hit(feature)],
                _feature_index=feature_index, describe=f"system {feature.qualifiers['name'][0]}")
            set_hit_best_score(region, feature_score(feature))
            if max_region_overlap < 1.0:
                if any(overlap_size(location, other) / len(location) > max_region_overlap for other in selected_locations):
                    continue
                selected_locations.append(location)
            yield region


def build_parser():
    parser = tool_parser(__doc__)
    add_input(parser, "GenBank files with system annotations (from find_systems).")
    add_seq_outputs(parser)

    add_system_selection(parser)
    parser.add_argument("--databases", type=str, nargs="+", default=None,
                        help="Only select systems from these model collections (the database qualifier). Default: all.")
    add_system_type(parser, help="Which systems to select.")
    parser.add_argument("--min_score", type=float, default=None,
                        help="Only select systems with at least this score.")
    parser.add_argument("--min_completeness", type=NonNegativeFloat, default=None,
                        help="Only select systems with at least this completeness (0 to 1). Systems of models without "
                             "a quorum have no completeness, and are not filtered by it.")
    parser.add_argument("--partial_systems", type=str, default="include", choices=PARTIAL_CHOICES,
                        help="Which systems to select by their system_partial flag, set when a system touches a linear "
                             "contig end and a missing component could lie beyond it. include: all. exclude: skip "
                             "flagged systems. only: only flagged systems.")
    parser.add_argument("--invert", action="store_true", default=False,
                        help="Invert --systems/--systems_file and --system_categories, not the other criteria. "
                             "With neither, warn and ignore.")
    parser.add_argument("--strand", type=str, default=None, choices=["f", "r"],
                        help="Only select systems in this orientation.")

    parser.add_argument("--kb_range", type=NonNegativeFloat, default=None,
                        help="How many kb upstream and downstream of the system to extract, from its outermost members.")
    parser.add_argument("--kb_up", type=NonNegativeFloat, default=None,
                        help="How many kb upstream of the system to extract, in system orientation.")
    parser.add_argument("--kb_down", type=NonNegativeFloat, default=None,
                        help="How many kb downstream of the system to extract, in system orientation.")
    parser.add_argument("--cds_range", type=NonNegativeInt, default=None,
                        help="How many CDSs upstream and downstream of the system's outermost members to extract. "
                             "Only protein-coding CDSs are counted, never hits on a nucleotide reference.")
    parser.add_argument("--cds_up", type=NonNegativeInt, default=None,
                        help="How many CDSs upstream of the system to extract, in system orientation.")
    parser.add_argument("--cds_down", type=NonNegativeInt, default=None,
                        help="How many CDSs downstream of the system to extract, in system orientation.")
    parser.add_argument("--whole_contig", action="store_true", default=False,
                        help="Extract the whole contig around each system.")
    add_model_neighborhood(parser)
    add_min_context(parser)

    parser.add_argument("--pad", action="store_true", default=False,
                        help="Pad the ends of the regions so that their systems align. See --anchor.")
    parser.add_argument("--anchor", type=str, default="start", choices=WHERE_CHOICES,
                        help="Where to align the systems when --pad is set, in system orientation. start: their first "
                             "base. end: their last base. center: their midpoints.")
    parser.add_argument("--keep_direction", action="store_true", default=False,
                        help="By default, regions are flipped so that their system is on the forward strand. Set this to "
                             "keep each region in the direction of its contig.")
    add_dedup(parser, default="id")
    parser.add_argument("--max_region_overlap", type=float, default=1.0,
                        help="The largest fraction of a region that may overlap a region already output from the same "
                             "contig. If >= 1, then no overlap filtering will be done.")
    add_list_and_file(parser, "contigs", help="Only read input records with these IDs; not inverted by --invert.")

    parser.add_argument("--config", action=ActionConfigFile)
    declare_outputs(parser, ["genbank_out", "fasta_out"])
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)
    check_outputs(parser, params)
    try:
        kb_range, cds_range = validate_range(params.kb_range, params.kb_up, params.kb_down,
                                             params.cds_range, params.cds_up, params.cds_down)
        min_context_cds = up_down_pair(params.min_context_cds, params.min_context_cds_up, params.min_context_cds_down,
                                       "--min_context_cds", "--min_context_cds_up", "--min_context_cds_down")
        min_context_kb = up_down_pair(params.min_context_kb, params.min_context_kb_up, params.min_context_kb_down,
                                      "--min_context_kb", "--min_context_kb_up", "--min_context_kb_down")
    except ValueError as error:
        parser.error(str(error))
    if params.model_neighborhood and (kb_range is not None or cds_range is not None or params.whole_contig):
        parser.error("--model_neighborhood can't be combined with --kb_*, --cds_* or --whole_contig.")
    if params.whole_contig and (kb_range is not None or cds_range is not None):
        parser.error("--whole_contig can't be combined with --kb_* or --cds_*.")

    target_contigs = list_and_file_to_dict_keys(params.contigs, params.contigs_file, as_set=True)
    selected = select_by_system(
        parse_seqfiles(params.input, target_contigs),
        systems=list_and_file_to_dict_keys(params.systems, params.systems_file, as_set=True),
        system_categories=params.system_categories,
        databases=params.databases,
        system_type=params.system_type,
        min_score=params.min_score,
        min_completeness=params.min_completeness,
        partial_systems=params.partial_systems,
        invert=params.invert,
        strand=params.strand,
        cds_range=cds_range,
        kb_range=kb_range,
        whole_contig=params.whole_contig,
        model_neighborhood=params.model_neighborhood,
        min_context_cds=min_context_cds,
        min_context_kb=min_context_kb,
        normalize_direction=not params.keep_direction,
        max_region_overlap=params.max_region_overlap,
    )
    with SeqOutputs(params.genbank_out, params.fasta_out) as outputs:
        outputs.write_all(deduplicate_and_pad(selected, dedup=params.dedup, pad=params.pad, anchor=params.anchor))
