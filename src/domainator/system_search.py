"""Find systems in sequence databases and select the regions around them.

A system is a set of domains found together in one genomic neighborhood (see find_systems).
system_search annotates the contigs of -i/--input with the profiles the system models use (and
any -r/--reference), finds the systems, and returns a region around each one, as domain_search
does for single hits:

    domainator system_search -i genomes.gb --system_models padloc_db --systems gabija --cds_range 5 -o gabija.gb
    domainator system_search -i genomes.gb -r toxins.hmm --domain_expr "TA_toxin & TA_antitoxin" --window_cds 2 -o ta.gb

Systems are defined in model files (--system_models), or on the command line with --domain_expr
and a window, whose names must be references given with -r. Inputs must be nucleotide contigs
with CDSs, or nucleotide sequences with --gene_call.

Padding (--cds_range, --kb_range, ...) is counted from the system's outermost members, in system
orientation, and --model_neighborhood pads each system by its model's own neighborhood.
--min_context_* skips systems without enough context beyond them, such as systems cut by the
end of a linear contig.

Each region has one System_Search feature, marking the system it was selected for, and the
Domainator hits of the system's members. --add_annotations keeps every hit in the region, and
the region's other systems (Domain_System features).

--systems/--systems_file, --system_expr and --system_categories choose which models to look
for. --system_expr is evaluated over the names of the systems found on each contig, and only
systems of the models it names positively become regions: "gabija & ~thoeris_I" returns Gabija
systems on contigs without a Thoeris I system.

Every model profile is searched with the larger of -e and its model's own threshold, and a hit
counts as evidence for a system only if it also passes -e (as in find_systems).

--max_hits applies to the whole search, keeping the best-scoring systems. When --max_hits or
--pad is given, regions are written after the search completes, best first; otherwise they are
written as they are found, in input order.
"""
import functools
import os
from typing import Iterable, Iterator, List, Optional, Sequence, Set, Tuple, Union

from Bio.SeqRecord import SeqRecord
from jsonargparse import ActionConfigFile
from jsonargparse.typing import NonNegativeFloat, NonNegativeInt, PositiveFloat

from domainator import DOMAIN_FEATURE_NAME, SYSTEM_FEATURE_NAME
from domainator import domainate
from domainator._mp import Shared, imap_bounded, make_pool, resolve_cpu
from domainator.cli.common import (add_dedup, add_input, add_min_context, add_model_neighborhood, add_seq_outputs,
                                   add_system_detection, add_system_models, add_system_selection, add_taxonomy_filters,
                                   add_window_model, check_outputs, declare_outputs, list_and_file_to_dict_keys,
                                   partial_help, reject_stdin, taxonomy_from_args, tool_parser)
from domainator.dedup import RecordDeduplicator, TopHits
from domainator.domain_search import _search_tasks, get_input_molecule_types, search_partitions
from domainator.errors import InputError
from domainator.expr import BooleanEvaluator
from domainator.extract_peptides import extract_peptides
from domainator.find_systems import annotate_systems
from domainator.locations import WHERE_CHOICES
from domainator.location_ops import overlap_size
from domainator.model import DomainatorCDS, get_hit_best_score, set_hit_best_score
from domainator import db_index
from domainator.output_guardrails import add_max_output_gb_argument, max_output_gb_to_bytes
from domainator.outputs import SeqOutputs
from domainator.partial import PARTIAL_CHOICES
from domainator.select_by_cds import deduplicate_and_pad, up_down_pair, validate_range
from domainator.select_by_system import select_by_system
from domainator.seqio import parse_seqfiles
from domainator.systems.features import component_values
from domainator.systems.model import ModelCollection
from domainator.systems.native import read_system_models, select_models, window_model
from domainator.systems.records import NUCLEOTIDE_MEMBER_PREFIX, system_features
from domainator.systems.search import model_profile_evalues, model_reference_files
from domainator.taxonomy import NCBITaxonomy
from domainator.taxonomy.records import compile_taxonomy_filter


def keep_member_hits(region: SeqRecord, databases: Set[str], system_databases: Set[str]) -> SeqRecord:
    """Remove from region the Domainator hits from ``databases`` that aren't evidence for its
    System_Search system, and the Domain_System features from ``system_databases``.

    A hit is evidence when it is on a member CDS (or, for a nucleotide member, overlaps the
    system) and its database and name are those the system's component qualifiers give for that
    member.
    """
    focus, = system_features(region, "system_search")
    values = component_values(focus)
    cds_hits = {(value.cds_id, value.database, value.reference) for value in values
                if not value.cds_id.startswith(NUCLEOTIDE_MEMBER_PREFIX)}
    nucleotide_hits = {(value.database, value.reference) for value in values
                       if value.cds_id.startswith(NUCLEOTIDE_MEMBER_PREFIX)}
    evidence = set()
    for cds in DomainatorCDS.list_from_contig(region, include_nucleic_acid_annotations=True):
        for hit in cds.domain_features:
            key = (hit.qualifiers["database"][0], hit.qualifiers["name"][0])
            if cds.is_nucleic_acid:
                if key in nucleotide_hits and overlap_size(hit.location, focus.location) > 0:
                    evidence.add(id(hit))
            elif (cds.num, *key) in cds_hits:
                evidence.add(id(hit))

    def kept(feature) -> bool:
        if feature.type == DOMAIN_FEATURE_NAME:
            return id(feature) in evidence or feature.qualifiers.get("database", [None])[0] not in databases
        if feature.type == SYSTEM_FEATURE_NAME:
            return feature.qualifiers.get("database", [None])[0] not in system_databases
        return True

    region.features = [feature for feature in region.features if kept(feature)]
    return region


class _system_search_worker():
    """Pool task: search one partition or batch of records, and return the regions around the
    systems found.

    ``references`` (the parsed -r references and model profiles) and ``allowed_taxids`` are
    Shared (see domainator._mp), so each worker gets them once.
    """

    def __init__(self, *, references: Shared, allowed_taxids: Shared, collections: Sequence[ModelCollection],
                 profile_evalues, focus: Set[str], system_expr: Optional[str], z: int, evalue: float,
                 min_evalue: float, max_mode: bool, batch_size: int, gene_call: Optional[str], partial: str,
                 detection: dict, selection: dict, add_annotations: bool, max_hits: Optional[int], dedup: str):
        self.references = references
        self.allowed_taxids = allowed_taxids
        self.collections = collections
        self.profile_evalues = profile_evalues
        self.focus = focus
        self.system_expr = system_expr
        self.z = z
        self.evalue = evalue
        self.min_evalue = min_evalue
        self.max_mode = max_mode
        self.batch_size = batch_size
        self.gene_call = gene_call
        self.partial = partial
        self.detection = detection
        self.selection = selection
        self.add_annotations = add_annotations
        self.max_hits = max_hits  # see _domain_search_worker: a task's top max_hits suffice
        self.dedup = dedup

    def __call__(self, task) -> List[SeqRecord]:
        if isinstance(task, list):
            records = task
        else:
            path, offset, count = task
            records = parse_seqfiles((path,), None, filetype_override=None, seek_to=offset, max_recs=count,
                                     default_molecule_type="nucleotide", genbank_parser="lean")
        references = self.references.get()
        databases = {database for groups in references.values() for database in groups}
        system_databases = {collection.database for collection in self.collections}
        evaluator = BooleanEvaluator(self.system_expr) if self.system_expr is not None else None
        annotated = domainate._domainate(
            records, references=None, pre_parsed_references=references, z=self.z, evalue=self.evalue,
            max_overlap=1, cpu=1, batch_size=self.batch_size, hits_only=True, gene_call=self.gene_call,
            min_evalue=self.min_evalue, max_mode=self.max_mode, allowed_taxids=self.allowed_taxids.get(),
            partial=self.partial, profile_evalues=self.profile_evalues)
        with_systems = annotate_systems(annotated, self.collections, evalue=self.evalue, partial=self.partial,
                                        **self.detection)
        out = TopHits(self.max_hits, self.dedup) if self.max_hits is not None else list()
        for record in with_systems:
            try:
                if evaluator is not None:
                    names = set()
                    for feature in system_features(record):
                        database, name = feature.qualifiers["database"][0], feature.qualifiers["name"][0]
                        if database in system_databases:
                            names.update((name, f"{database}:{name}"))
                    if not evaluator.check_expression(names):
                        continue
                for region in select_by_system((record,), systems=self.focus, **self.selection):
                    if not self.add_annotations:
                        keep_member_hits(region, databases, system_databases)
                    if isinstance(out, list):
                        out.append(region)
                    else:
                        out.add(region)
            except Exception as e:
                e.add_note(f"while processing {record.id}")
                raise
        return out if isinstance(out, list) else out.records()


def _translate(regions: Iterable[SeqRecord]) -> Iterator[SeqRecord]:
    """One protein per member CDS of each region's System_Search system, scored as the system."""
    for region in regions:
        focus, = system_features(region, "system_search")
        name = {focus.qualifiers["database"][0] + ":" + focus.qualifiers["name"][0]}
        score = get_hit_best_score(region)
        for peptide in extract_peptides((region,), systems=name, system_type="system_search", keep_cds_feature=True):
            set_hit_best_score(peptide, score)
            yield peptide


def system_search(targets: Iterable[Union[SeqRecord, Tuple[str, int, int]]], *, z: int,
                  references: Optional[List[str]] = None,
                  system_models: Optional[Sequence[Union[str, os.PathLike, ModelCollection]]] = None,
                  domain_expr: Optional[str] = None, window_cds: Optional[int] = None, window_kb: Optional[float] = None,
                  window_mode: str = "span", system_name: str = "domain_expr", same_strand: bool = False,
                  systems: Optional[Set[str]] = None, system_expr: Optional[str] = None,
                  system_categories: Optional[Set[str]] = None, evalue: float = 0.001, min_evalue: float = 0.0,
                  max_mode: bool = False, min_score: Optional[float] = None, min_completeness: Optional[float] = None,
                  max_system_overlap: float = 1.0, max_systems_per_contig: Optional[int] = None,
                  max_hits: Optional[int] = None, cpu: int = 0, batch_size: int = 10000, add_annotations: bool = False,
                  cds_range: Optional[Tuple[int, int]] = None, kb_range: Optional[Tuple[float, float]] = None,
                  whole_contig: bool = False, model_neighborhood: bool = False,
                  min_context_cds: Optional[Tuple[int, int]] = None, min_context_kb: Optional[Tuple[float, float]] = None,
                  normalize_direction: bool = True, strand: Optional[str] = None, max_region_overlap: float = 1.0,
                  translate: bool = False, gene_call: Optional[str] = None,
                  ncbi_taxonomy: Optional[NCBITaxonomy] = None, include_taxids: Optional[Iterable[int]] = None,
                  exclude_taxids: Optional[Iterable[int]] = None, taxonomy_expr: Optional[str] = None,
                  taxonomy_unknown: str = "exclude", partial: str = "include",
                  dedup: str = "none") -> Iterator[SeqRecord]:
    """Search ``targets`` for systems in parallel and yield the regions around them.

    The targets are annotated in batches with the references and the models' profiles (as by
    domainate --system_models), the systems are found (as by find_systems), and the region
    around each is extracted (as by select_by_system).

    Args:
        targets: nucleotide SeqRecords, or (sequence file path, file offset, number of records)
            partitions of sequence files, but not both (see domain_search.domain_search).
        z: the database size for e-values (hmmsearch's Z), as ``-Z``.
        references: references to search besides the models' profiles: the names a
            command-line model uses, and other annotations.
        system_models: model files, directories holding systems.yaml, or ModelCollections.
        domain_expr, window_cds, window_kb, window_mode, system_name, same_strand: a model
            defined on the command line, as for find_systems. At least one of system_models and
            domain_expr is required.
        systems: look only for the models with these names (``name`` or ``database:name``).
        system_expr: a boolean expression over the names of the systems found on a contig. Its
            names are looked for too, and only systems of models it names positively are
            selected, on contigs where it holds.
        system_categories: look only for models in these categories.
        evalue: search with this e-value (model profiles: or their model's looser threshold),
            and count only hits passing it as evidence.
        min_evalue, max_mode: as for domainate.
        min_score, min_completeness, max_system_overlap, max_systems_per_contig: as for
            find_systems.
        max_hits: keep only this many systems, the best-scoring. None: all.
        cpu: number of workers (0: one per physical core); values below 2 are raised to 2.
        batch_size: about how many CDSs a worker searches at once.
        add_annotations: keep every hit in each region, and its other systems. Otherwise a region
            keeps only its system's member hits from the searched databases.
        cds_range, kb_range, whole_contig, model_neighborhood, min_context_cds, min_context_kb,
            normalize_direction, max_region_overlap: as for select_by_system.
        strand: "f" or "r" to select only systems in that orientation. None: both.
        translate: yield one protein per member CDS of each system, instead of regions.
        gene_call: as for domainate.
        ncbi_taxonomy, include_taxids, exclude_taxids, taxonomy_expr, taxonomy_unknown: as for
            domain_search.
        partial: which CDSs to search and count as members, by fragment status.
        dedup: output uniqueness policy, as for domain_search.

    Yields:
        The regions, best first with max_hits, else in input order; or with translate, the
        proteins of their members.
    """
    if max_hits is not None and max_hits < 1:
        raise InputError("max_hits must be at least 1.")
    if dedup == "sequence_rc" and translate:
        raise InputError("--dedup sequence_rc is only supported for nucleotide output.")
    collections = read_system_models(system_models or [])
    if domain_expr is not None:
        collections.append(window_model(domain_expr, window_cds=window_cds, window_kb=window_kb, window_mode=window_mode,
                                        system_name=system_name, same_strand=same_strand))
    if not collections:
        raise InputError("No system models: give system_models (--system_models), domain_expr (--domain_expr), or both.")
    evaluator = BooleanEvaluator(system_expr) if system_expr is not None else None
    named = None
    if systems is not None or evaluator is not None:
        named = set(systems or ()) | (evaluator.names if evaluator is not None else set())
    collections = select_models(collections, named, system_categories)
    if not collections:
        raise InputError("No system models are in the --system_categories given.")
    focus_names = None
    if systems is not None or evaluator is not None:
        focus_names = set(systems or ()) | (evaluator.positive_names if evaluator is not None else set())
        if not focus_names:
            raise InputError(f"--system_expr {system_expr!r} names no system positively, so no system would be selected.")
    focus = {f"{model.database}:{model.name}" for collection in collections for model in collection.models
             if focus_names is None or {model.name, f"{model.database}:{model.name}"} & focus_names}

    model_files = model_reference_files(collections)
    model_paths = {os.path.realpath(path) for path in model_files}
    # A reference that is also a model's reference file is searched once, as a model profile.
    search_files = [path for path in references or [] if os.path.realpath(path) not in model_paths] + model_files
    if not search_files:
        raise InputError("Nothing to search: give references (-r), or system models that list reference_files.")
    profile_evalues = model_profile_evalues(collections, model_files, evalue)
    cpu = max(2, resolve_cpu(cpu))
    allowed_taxids = compile_taxonomy_filter(ncbi_taxonomy, include_taxids, exclude_taxids, taxonomy_expr,
                                            taxonomy_unknown)
    shared_references = Shared(build=functools.partial(domainate.read_references, search_files, foldseek=None))
    shared_taxids = Shared(allowed_taxids)
    detection = dict(min_score=min_score, min_completeness=min_completeness, max_system_overlap=max_system_overlap,
                     max_systems_per_contig=max_systems_per_contig)
    selection = dict(strand=strand, cds_range=cds_range, kb_range=kb_range, whole_contig=whole_contig,
                     model_neighborhood=model_neighborhood, min_context_cds=min_context_cds,
                     min_context_kb=min_context_kb, normalize_direction=normalize_direction,
                     max_region_overlap=max_region_overlap)
    # Translated output is deduplicated as proteins, below, so the regions aren't.
    region_dedup = "none" if translate else dedup
    deduplicator = RecordDeduplicator(region_dedup)
    worker = _system_search_worker(
        references=shared_references, allowed_taxids=shared_taxids, collections=collections,
        profile_evalues=profile_evalues, focus=focus, system_expr=system_expr, z=z, evalue=evalue,
        min_evalue=min_evalue, max_mode=max_mode, batch_size=batch_size, gene_call=gene_call, partial=partial,
        detection=detection, selection=selection, add_annotations=add_annotations, max_hits=max_hits,
        dedup=region_dedup)
    top_hits = TopHits(max_hits, region_dedup) if max_hits is not None else None

    def regions() -> Iterator[SeqRecord]:
        workers = cpu - 1  # the parent merges and writes (see domain_search)
        with make_pool(processes=workers, shared=[shared_references, shared_taxids]) as pool:
            for hits in imap_bounded(pool, worker, _search_tasks(targets, batch_size, "nucleotide"),
                                     max_pending=2 * workers):
                for region in hits:
                    if top_hits is not None:
                        top_hits.add(region)
                    elif deduplicator.keep(region):
                        yield region
        if top_hits is not None:
            yield from top_hits.records()

    if not translate:
        yield from regions()
        return
    peptide_deduplicator = RecordDeduplicator(dedup)
    for peptide in _translate(regions()):
        if peptide_deduplicator.keep(peptide):
            yield peptide


def build_parser():
    parser = tool_parser(__doc__)
    add_input(parser, "The GenBank or FASTA files to search: nucleotide contigs with CDSs, or nucleotide sequences "
                      "with --gene_call.", stdin=False)
    parser.add_argument("-r", "--reference", type=str, nargs="+", default=None,
                        help="References to search besides the models' profiles (HMM, FASTA or CM files): the names "
                             "--domain_expr uses, and other annotations, which are kept with --add_annotations.")
    add_system_models(parser)
    add_window_model(parser)
    add_seq_outputs(parser, genbank_help="Write the regions around the systems to this GenBank file.",
                    fasta_help="Write the regions around the systems to this FASTA file.")

    add_system_selection(parser, expr=True)
    parser.add_argument("-Z", default=1000, type=int,
                        help="The database size for e-values (hmmsearch's -Z). 0: the number of CDSs in the inputs, "
                             "which takes an extra pass to count unless they have .didx indexes.")
    parser.add_argument("-e", "--evalue", type=PositiveFloat, default=0.001,
                        help="Count only hits with an e-value at most this as evidence. Model profiles are searched "
                             "with this or their model's own threshold, whichever is larger.")
    parser.add_argument("--min_evalue", type=NonNegativeFloat, default=0.0,
                        help="Hits with an e-value below this are not reported.")
    parser.add_argument("--max_mode", action="store_true", default=False,
                        help="Turn off hmmsearch's heuristic filters: slower, but finds weak hits.")
    parser.add_argument("--min_score", type=float, default=None,
                        help="Only select systems with at least this score.")
    parser.add_argument("--min_completeness", type=NonNegativeFloat, default=None,
                        help="Only select systems with at least this completeness (0 to 1). Systems of models without "
                             "a quorum have no completeness, and are not filtered by it.")
    add_system_detection(parser)
    parser.add_argument("--max_hits", type=int, default=None,
                        help="The most systems to return, best score first. Default: all.")
    add_taxonomy_filters(parser)
    parser.add_argument("--partial", type=str, default="include", choices=PARTIAL_CHOICES,
                        help=partial_help("Which CDSs to search and count as system members"))
    parser.add_argument("--cpu", type=int, default=0,
                        help="The number of processes to use (0: one per physical core; at least 2).")
    parser.add_argument("--add_annotations", action="store_true", default=False,
                        help="Keep every hit in each region, and its other systems. Default: keep only the hits of "
                             "the system's members.")
    add_dedup(parser)
    parser.add_argument("--max_region_overlap", type=float, default=1.0,
                        help="The largest fraction of a region that may overlap a region already output from the same "
                             "contig. If >= 1, then no overlap filtering will be done.")
    parser.add_argument("--translate", action="store_true", default=False,
                        help="Output one protein per member CDS of each system, instead of the regions. Each protein's "
                             "CDS feature has a system_membership qualifier.")

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
    parser.add_argument("--strand", type=str, default=None, choices=["f", "r"],
                        help="Only select systems in this orientation.")
    parser.add_argument("--pad", action="store_true", default=False,
                        help="Pad the ends of the regions so that their systems align. See --anchor.")
    parser.add_argument("--anchor", type=str, default="start", choices=WHERE_CHOICES,
                        help="Where to align the systems when --pad is set, in system orientation. start: their first "
                             "base. end: their last base. center: their midpoints.")
    parser.add_argument("--keep_direction", action="store_true", default=False,
                        help="By default, regions are flipped so that their system is on the forward strand. Set this to "
                             "keep each region in the direction of its contig.")
    parser.add_argument("--batch_size", type=int, default=10000,
                        help="About how many CDSs each process searches at once.")
    parser.add_argument("--gene_call", type=str, default=None, choices=["all", "unannotated"],
                        help="Call genes with Prodigal (metagenomic mode) before searching. all: replace every CDS. "
                             "unannotated: only contigs without CDSs.")
    add_max_output_gb_argument(parser)
    parser.add_argument("--config", action=ActionConfigFile)
    declare_outputs(parser, ["genbank_out", "fasta_out"])
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)
    check_outputs(parser, params)
    reject_stdin(parser, params.input, "it reads its input more than once (to count and partition it)")
    if params.system_models is None and params.domain_expr is None:
        parser.error("Give --system_models, --domain_expr, or both.")
    if params.evalue > 10:
        parser.error("-e/--evalue must be <= 10.")
    if params.Z == 0 and params.gene_call is not None:
        parser.error("-Z 0 can't be used with --gene_call: give the database size explicitly.")
    if params.batch_size < 1:
        parser.error("--batch_size must be at least 1.")
    if params.max_hits is not None and params.max_hits < 1:
        parser.error("--max_hits must be at least 1. Leave it out to return all systems.")
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
    if params.translate and (kb_range is not None or cds_range is not None or params.whole_contig or params.pad
                             or params.model_neighborhood):
        parser.error("--translate outputs proteins, so it can't be combined with neighborhood extraction or --pad.")
    if params.dedup == "sequence_rc" and params.translate:
        parser.error("--dedup sequence_rc is only supported for nucleotide output.")

    cpus = resolve_cpu(params.cpu)
    input_files = [shard for db in params.input for shard in db_index.resolve_database_shards(db)]
    if "protein" in get_input_molecule_types(input_files, default_molecule_type="nucleotide"):
        parser.error("system_search needs nucleotide inputs: systems are found on contigs.")
    partitions, z = search_partitions(input_files, params.Z, params.batch_size, cpus)

    mitigation_options = ["--max_hits", "--cds_range", "--kb_range"]
    extra_guidance = "Avoid --whole_contig on broad searches unless you expect very few systems." if params.whole_contig else None
    outputs = SeqOutputs(params.genbank_out, params.fasta_out, max_output_bytes=max_output_gb_to_bytes(params.max_output_gb),
                         mitigation_options=mitigation_options, extra_guidance=extra_guidance)
    with outputs:
        found = system_search(
            partitions,
            z=z,
            references=params.reference,
            system_models=params.system_models,
            domain_expr=params.domain_expr,
            window_cds=params.window_cds,
            window_kb=params.window_kb,
            window_mode=params.window_mode,
            system_name=params.system_name,
            same_strand=params.same_strand,
            systems=list_and_file_to_dict_keys(params.systems, params.systems_file, as_set=True),
            system_expr=params.system_expr,
            system_categories=params.system_categories,
            evalue=params.evalue,
            min_evalue=params.min_evalue,
            max_mode=params.max_mode,
            min_score=params.min_score,
            min_completeness=params.min_completeness,
            max_system_overlap=params.max_system_overlap,
            max_systems_per_contig=params.max_systems_per_contig,
            max_hits=params.max_hits,
            cpu=cpus,
            batch_size=params.batch_size,
            add_annotations=params.add_annotations,
            cds_range=cds_range,
            kb_range=kb_range,
            whole_contig=params.whole_contig,
            model_neighborhood=params.model_neighborhood,
            min_context_cds=min_context_cds,
            min_context_kb=min_context_kb,
            normalize_direction=not params.keep_direction,
            strand=params.strand,
            max_region_overlap=params.max_region_overlap,
            translate=params.translate,
            gene_call=params.gene_call,
            ncbi_taxonomy=taxonomy_from_args(params),
            include_taxids=params.include_taxids,
            exclude_taxids=params.exclude_taxids,
            taxonomy_expr=params.taxonomy_expr,
            taxonomy_unknown=params.taxonomy_unknown,
            partial=params.partial,
            dedup=params.dedup,
        )
        outputs.write_all(deduplicate_and_pad(found, pad=params.pad, anchor=params.anchor, best_first=True))
