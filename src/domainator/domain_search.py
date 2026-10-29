"""Select sequence records or hit neighborhoods matching protein or nucleic-acid references.

Only matching records or regions are returned. Use domainate to annotate and retain
all input records, including those without qualifying hits.
This tool is optimized for large sequence databases with relatively few references.
The -i/--input files are the database being searched; -r/--reference supplies the
references to search for.
Protein references are searched against CDSs or protein records.
Nucleic-acid references are searched against whole nucleotide contigs.
Supported references include protein HMMs, protein FASTA, nucleotide FASTA, nucleotide HMMs, and infernal CM files.
Protein inputs cannot be searched with nucleotide references.

Protein input yields matching protein records. Nucleotide input yields hit-centered
regions or, with --whole_contig, whole contigs. For protein-reference searches of
CDSs, --translate yields peptides from matching CDSs. The output sequence need not
be the same length as the aligned region.

Search best-hit annotations (Domain_Search features) mark the selected hits.
--add_annotations also adds domain annotations (Domainator features) for retained
hits; it does not enable the best-hit annotations, which are already present.

When --max_hits is specified, it applies to the entire search, not individual references. Thus if you supply 5 reference profiles and --max_hits 100,
at most 100 hits will be returned (not 500). 

When max_hits or pad are specified, hits will be sorted by best domain score and written after the entire search completes.
If neither is set, then hits will be written on the fly and not sorted.

--dedup selects one output uniqueness policy: none (the default), id,
sequence, or sequence_rc. With --max_hits, deduplication keeps the best-scoring
representatives before applying the limit, retaining at most --max_hits winners
per worker and in the global merge. Domainator assumes that the same record ID
always identifies the same sequence, so sequence uniqueness also ensures unique
IDs. Hit scores and annotations can still differ between copies.
Without --max_hits, the first arriving representative is kept; parallel arrival
order is arbitrary.
"""

from jsonargparse import ActionConfigFile
from jsonargparse.typing import NonNegativeFloat, PositiveFloat
from domainator.cli.common import add_alphabet, add_dedup, add_input, add_seq_outputs, check_outputs, declare_outputs, reject_stdin, tool_parser
from domainator.dedup import RecordDeduplicator, TopHits
from domainator.errors import InputError
from domainator.outputs import SeqOutputs
from domainator.seqio import parse_seqfiles, native_parser_available
from domainator.offsets import count_peptides_in_record, partition_seqfile, i_partition_seqfiles, index_total_cds
from domainator.cli.common import list_and_file_to_dict_keys
from domainator._mp import Shared, imap_unordered_bounded, make_pool, resolve_cpu
from domainator.taxonomy.records import compile_taxonomy_filter
from domainator.cli.common import add_taxonomy_filters, taxonomy_from_args
from domainator.partial import PARTIAL_CHOICES
from domainator.locations import WHERE_CHOICES
from domainator.model import get_hit_names, set_dist_to_start
from domainator import select_by_cds
from domainator import db_index
from domainator import domainate
import functools
from typing import List, Tuple, Set, Optional, Iterable, Iterator, Union
from domainator import extract_peptides
from pathlib import Path
from domainator.taxonomy import NCBITaxonomy
from domainator.output_guardrails import add_max_output_gb_argument, max_output_gb_to_bytes
from Bio.SeqRecord import SeqRecord


class _domain_search_worker():
    """Pool task: search one partition or batch of records, and return its hit regions.

    ``references`` and ``allowed_taxids`` are Shared (see domainator._mp), so each worker gets them
    once: the parsed references, and the compiled taxonomy filter (None: no filter).
    """

    def __init__(self, *, references: Shared, allowed_taxids: Shared, z: int, evalue: float, max_overlap: float, add_annotations: bool, cds_range: Tuple, kb_range: Tuple, whole_contig: bool, normalize_direction: bool, translate: bool, gene_call:str = None, min_evalue:float = 0.0, batch_size: int = 10000, fasta_type: str = "protein", max_mode: bool = False, max_region_overlap=1.0, strand: Optional[str] = None, decoy_names: Optional[Set[str]] = None, max_hits_per_contig: Optional[int] = None, max_hits: Optional[int] = None, partial: str = "include", dedup: str = "none"):

        self.references = references
        self.allowed_taxids = allowed_taxids
        self.z = z
        self.evalue = evalue
        self.max_overlap = max_overlap
        self.add_annotations = add_annotations
        self.cds_range = cds_range
        self.kb_range = kb_range
        self.whole_contig = whole_contig
        self.normalize_direction = normalize_direction
        self.batch_size = batch_size
        self.translate = translate
        self.gene_call = gene_call
        self.min_evalue = min_evalue
        self.fasta_type = fasta_type
        self.max_mode = max_mode
        self.max_region_overlap = max_region_overlap
        self.strand = strand
        self.decoy_names = decoy_names
        self.max_hits_per_contig = max_hits_per_contig
        self.partial = partial
        # A task's top-N unique keys suffice for the global top-N: a discarded
        # key already has N better keys in that task, and merging cannot remove
        # those distinct keys under the ID-to-sequence invariant.
        self.max_hits = max_hits
        self.dedup = dedup

    def __call__(self, task):
        """Search one task: a (sequence file path, file offset, number of records) partition, or a
        list of SeqRecords.

        Returns:
            list of SeqRecords. When max_hits is set, at most max_hits records (the highest
            scoring of this task) are returned.
        """
        if isinstance(task, list):
            records = task
        else:
            path, offset, count = task
            records = parse_seqfiles((path,), None, filetype_override=None, seek_to=offset, max_recs=count, default_molecule_type=self.fasta_type, genbank_parser="lean")
        out = list()
        top_hits = TopHits(self.max_hits, self.dedup) if self.max_hits is not None else None

        def _emit(record):
            if top_hits is None:
                out.append(record)
            else:
                top_hits.add(record)

        for rec in domainate._domainate(
                records, #TODO: clear best_domain_hit ?
                references=None,
                pre_parsed_references=self.references.get(),
                z=self.z,
                evalue=self.evalue,
                max_overlap=self.max_overlap,
                no_annotations=not self.add_annotations,
                cpu=1,
                batch_size=self.batch_size,
                hits_only=True,
                best_annotation=True,
                gene_call=self.gene_call,
                min_evalue=self.min_evalue,
                max_mode=self.max_mode,
                max_hits_per_contig=self.max_hits_per_contig,
                allowed_taxids=self.allowed_taxids.get(),  # per-CDS / per-region pre-hmmer taxonomy filter
                partial=self.partial,
            ):
            try:
                if rec.annotations['molecule_type'] == "protein":
                    if self.decoy_names is not None and get_hit_names(rec)[0] in self.decoy_names:
                        continue
                    set_dist_to_start(rec, 0)
                    _emit(rec)
                else: # it's a nucleotide, so pass it through select_by_cds
                    if not self.translate: 
                        for region in select_by_cds.select_by_cds(
                            (rec,),
                            cds_ids=None,
                            domains=None,
                            domain_expr=None,
                            cds_range=self.cds_range,
                            kb_range=self.kb_range,
                            whole_contig=self.whole_contig,
                            normalize_direction=self.normalize_direction,
                            max_region_overlap=self.max_region_overlap,
                            strand=self.strand,
                            _from_domain_search=True,
                            _domain_search_negatives=self.decoy_names
                            ):
                                _emit(region)
                    else: # self.translate == True, so we need to translate it
                        for peptide in extract_peptides.extract_peptides(
                            (rec,),
                            evalue=100000000,
                            domains=None,
                            cds_ids=None,
                            keep_cds_feature=True,
                            strand=self.strand,
                            _from_domain_search=True,
                            _domain_search_negatives=self.decoy_names
                        ):
                            set_dist_to_start(rec, 0)
                            _emit(peptide)
            except Exception as e:
                e.add_note(f"while processing {rec.id}")
                raise
        if top_hits is None:
            return out
        return top_hits.records()
 
class _partition_seqfile_worker():
    def __init__(self, cdss_per_partition):
        self.cdss_per_partition = cdss_per_partition
    
    def __call__(self, input_file):
        return partition_seqfile(input_file, self.cdss_per_partition)
    

def _search_tasks(targets, batch_size, alphabet):
    """The pool tasks for ``targets``: each partition as it is, and records in lists of about
    ``batch_size`` CDSs. A record without a molecule type gets ``alphabet``, as parse_seqfiles
    gives FASTA records."""
    kind = None
    batch, batch_cdss = [], 0
    for target in targets:
        target_kind = "partitions" if isinstance(target, tuple) else "records"
        if kind is None:
            kind = target_kind
        elif target_kind != kind:
            raise TypeError("domain_search targets must be all SeqRecords or all (path, offset, count) partitions, not both.")
        if kind == "partitions":
            yield target
            continue
        target.annotations.setdefault("molecule_type", alphabet)
        batch.append(target)
        batch_cdss += max(1, count_peptides_in_record(target))  # a contig without CDSs still counts
        if batch_cdss >= batch_size:
            yield batch
            batch, batch_cdss = [], 0
    if batch:
        yield batch


def domain_search(targets: Iterable[Union[SeqRecord, Tuple[str, int, int]]], *, references: List[str], z: int, evalue: float = 0.001,
                  max_hits: Optional[int] = None, max_overlap: float = 1, cpu: int = 0, batch_size: int = 10000, add_annotations: bool = False,
                  cds_range: Optional[Tuple[int, int]] = None, kb_range: Optional[Tuple[float, float]] = None,
                  whole_contig: bool = False, normalize_direction: bool = True, translate: bool = False, gene_call: Optional[str] = None, min_evalue: float = 0.0,
                  ncbi_taxonomy: Optional[NCBITaxonomy] = None, include_taxids: Optional[Iterable[int]] = None,
                  exclude_taxids: Optional[Iterable[int]] = None, taxonomy_expr: Optional[str] = None,
                  taxonomy_unknown: str = "exclude",
                  alphabet: str = "protein", max_mode: bool = False, max_region_overlap: float = 1.0,
                  strand: Optional[str] = None, decoys: Optional[Set[str]] = None, max_hits_per_contig: Optional[int] = None,
                  partial: str = "include", dedup: str = "none") -> Iterator[SeqRecord]:
    """Search ``references`` against ``targets`` in parallel and yield the hit regions.

    The targets are searched in batches, each annotated as by domainate, and the region around
    each hit is extracted as by select_by_cds.

    Args:
        targets: the sequences to search: SeqRecords (as parse_seqfiles reads them), or
            (sequence file path, file offset, number of records) partitions of sequence files
            (as offsets.i_partition_seqfiles makes them), but not both. Records are grouped into
            batches of about ``batch_size`` CDSs, and read only a few batches ahead of the
            search, so they can be a stream; a partition is one batch, read by the worker that
            searches it.
        references: references to search for: protein HMMs/FASTA, nucleotide HMMs/FASTA,
            or Infernal covariance models.
        z: the database size for e-values (hmmsearch's Z), as ``-Z``.
        evalue: keep hits with an e-value at most this.
        max_hits: keep only this many hits, the best-scoring. None: all.
        max_overlap: the largest fraction of a hit that may overlap a better hit on the same CDS.
        cpu: number of workers (0: one per physical core); values below 2 are raised to 2 (one
            merges the results, the others search).
        batch_size: about how many CDSs (or proteins) a worker searches at once.
        add_annotations: also add domain annotations (Domainator features) for retained
            hits, in addition to the search best-hit annotations (Domain_Search features).
        cds_range: (upstream, downstream) numbers of CDSs to extract around each hit CDS.
        kb_range: (upstream, downstream) kb to extract around each hit CDS. Partially enclosed
            features are truncated, with '<' or '>' marking the cut ends.
        whole_contig: extract the whole contig of each hit (once, if it has several hits).
        normalize_direction: reverse-complement regions whose hit CDS is on the reverse strand.
        translate: emit peptides from matching CDSs in protein-reference searches.
        gene_call: as for domainate.
        min_evalue: as for domainate.
        max_mode: as for domainate.
        ncbi_taxonomy: an NCBITaxonomy, needed by the taxonomy filters.
        include_taxids: search proteins/CDSs or nucleotide hits in these taxids or descendants.
        exclude_taxids: exclude these taxids or descendants, overriding inclusion.
        taxonomy_expr: a boolean expression over taxids (operators & | ~ and parentheses), e.g.
            "2 & ~1224". Mutually exclusive with include_taxids/exclude_taxids.
        taxonomy_unknown: unresolved taxonomy policy: exclude (default), include, evaluate
            against an empty lineage, or error. Protein/CDS targets are filtered before
            searching, nucleotide hits before ranking; expanded output is not re-filtered.
        alphabet: "protein" or "nucleotide", the molecule type of FASTA databases.
        max_region_overlap: the largest fraction of overlap allowed between two output regions
            (1 or more: no filter). Regions are kept greedily in order of CDS start, if less than
            this fraction of them overlaps a region already kept.
        strand: "f" or "r" to extract only regions around CDSs on that strand. None: both.
        decoys: names of decoy domains: a CDS or protein whose best hit is a decoy is not a hit.
        max_hits_per_contig: the most independent hits per contig. None: all. For contig-level
            nucleotide hits, each non-overlapping hit region is an independent hit.
        partial: which proteins/CDSs to search, by fragment status: "include" (all, the
            default), "exclude" (skip fragments), or "only" (only fragments). See
            domainate.domainate.
        dedup: output uniqueness policy: "none" (the default), "id", "sequence",
            or "sequence_rc". Records with the same ID must have the same sequence,
            so either sequence policy also ensures unique IDs.
            Sequence comparisons ignore case and terminal '*'; sequence_rc also matches
            nucleotide reverse complements. With max_hits, each worker and the parent
            retain at most max_hits winners, replacing copies with better-scoring hits,
            and sorting only the final winners. Equal scores favor earlier arrivals.
            Without max_hits, the first arriving representative is kept; parallel
            arrival order is arbitrary.

    Yields:
        Matching protein records, nucleotide hit regions/whole contigs, or translated
        CDSs, according to the input and output options. Search best-hit feature
        locations mark alignments, not necessarily the extent of the output sequence.
    """

    deduplicator = RecordDeduplicator(dedup)
    if dedup == "sequence_rc" and translate:
        raise InputError("--dedup sequence_rc is only supported for nucleotide output.")
    if max_hits is not None and max_hits < 1:
        raise InputError("max_hits must be at least 1.")
    cpu = max(2, resolve_cpu(cpu))

    allowed_taxids = compile_taxonomy_filter(ncbi_taxonomy, include_taxids, exclude_taxids,
                                            taxonomy_expr, taxonomy_unknown)
    ncbi_taxonomy = None

    # Each worker process parses the references itself, which is cheaper than
    # pickling the parsed pyhmmer/pyinfernal objects; threads share one parsed copy.
    shared_references = Shared(build=functools.partial(domainate.read_references, references, foldseek=None))
    shared_taxids = Shared(allowed_taxids)

    top_hits = TopHits(max_hits, dedup) if max_hits is not None else None
    worker = _domain_search_worker(references=shared_references, allowed_taxids=shared_taxids, z=z, evalue=evalue, max_overlap=max_overlap, add_annotations=add_annotations, cds_range=cds_range, kb_range=kb_range, whole_contig=whole_contig, normalize_direction=normalize_direction, translate=translate, gene_call=gene_call, min_evalue=min_evalue, batch_size=batch_size, fasta_type=alphabet, max_mode=max_mode, max_region_overlap=max_region_overlap, strand=strand, decoy_names=decoys, max_hits_per_contig=max_hits_per_contig, max_hits=max_hits, partial=partial, dedup=dedup)

    # Pool tuning notes:
    #   - Workers persist for the whole search, so the references parsed once per
    #     worker stay loaded.
    #   - Tasks are large (~batch_size CDSs each), so per-task dispatch overhead is
    #     negligible. Two tasks per worker are in flight, which keeps the workers busy
    #     while reading records only that far ahead.
    #   - cpu - 1 workers: the parent does the heap merge and GenBank writing, so it
    #     is reserved one core rather than oversubscribing.
    workers = cpu - 1
    with make_pool(processes=workers, shared=[shared_references, shared_taxids]) as pool:
        # When max_hits is set, each list is already deduplicated and
        # trimmed to that task's own top-max_hits (see _domain_search_worker), so
        # the parent merges O(tasks * max_hits) records rather than every hit; the
        # central heap below still selects the global top-max_hits.
        for hits in imap_unordered_bounded(pool, worker, _search_tasks(targets, batch_size, alphabet), max_pending=2 * workers):

            for rec in hits:
                if top_hits is None:
                    if deduplicator.keep(rec):
                        yield rec
                else:
                    top_hits.add(rec)

    if top_hits is not None:
        yield from top_hits.records()

def get_input_molecule_types(input_files: List[str], default_molecule_type="protein") -> Set[str]:
    """reads the first record from each input file and gets the 'molecule_type' annotation

    Args:
        input_files (List[str]): a list of input file paths
        default_molecule_type (str, optional): the molecule type to use if the input file does not have a molecule_type annotation. Defaults to "protein".
        
    
    Returns:
        a set of molecule type names.
    """

    # Only the first record of each file is read (max_recs=1) just to learn its
    # molecule type; use the fast native parser (it falls back to Biopython for
    # records it cannot parse) so this peek is cheap on large GenBank inputs.
    return {x.annotations['molecule_type'] for x in parse_seqfiles(input_files, default_molecule_type=default_molecule_type, max_recs=1, genbank_parser="lean")}



def build_parser():
    if native_parser_available():
        parser_status = "Native (Rust) parser: installed (fast GenBank/FASTA parsing and offset scanning)."
    else:
        parser_status = ("Native (Rust) parser: NOT installed; using the slower Biopython parser. "
                         "See the README 'Native (Rust) acceleration' section to build it.")
    parser = tool_parser(f"{parser_status}\n\n" + __doc__)

    add_input(parser, "the genbank or fasta files to search. Genbank files can be nucleotide (with CDS annotations) or peptide. Fasta files can be protein or nucleotide; use --alphabet to specify which.", stdin=False)
    add_alphabet(parser)
    parser.add_argument('-r', '--reference', required=True, type=str,
                        default=None, nargs='+',
                        help="References to search for in the --input records. Supported types are protein HMMs, protein FASTA, nucleotide FASTA, nucleotide HMMs, and infernal CM files.")

    add_seq_outputs(parser, genbank_help="Write the hits to this GenBank file. "
                    "If the filename ends in .gz/.bgz, the output is written as BGZF (block-gzip): "
                    "it is gunzip-readable and ~3-4x smaller, and can be used directly as a "
                    "random-access, partitionable compressed input to a later domainator step.",
                    fasta_help="Write the hits to this FASTA file.")


    parser.add_argument('-Z', default=1000, type=int, 
                        help="Passed as the -Z parameter to hmmsearch: Assert that the total number of targets in your searches is <x>, for the purposes of per-sequence E-value calculations, rather than the actual number of targets seen. "
                        "\n"
                        "Set to 0 to use the number of target sequences in the input, counted before any taxonomy filter, so filtering does not change E-values. "
                        "Supplying -Z is recommended for most use cases of domainate, as it will speed up the search and make comparisons between searches more meaningful.")

    parser.add_argument('-e', '--evalue', type=PositiveFloat, default=0.001,
                        help="threshold E value for the domain hit. Must be <=10.")
    
    parser.add_argument('--min_evalue', type=NonNegativeFloat, default=0.0,
                        help="hits with E value lower than this will be filtered out. Not frequently used. Use only if you want to eliminate close matches. Must be >=0. [default 0]")

    parser.add_argument('--max_mode', action="store_true", default=False,
                        help="Run hmmsearch, phmmer, or nhmmer in maximum sensitivity mode, which is much slower, but more sensitive.")

    parser.add_argument('--max_hits', type=int, default=None,
                        help="the maximum number of CDSs returned by the search, keeping the best-scoring ones by the bitscore of their best scoring profile. If not supplied, return all hits passing the evalue threshold.")
    
    parser.add_argument('--max_hits_per_contig', type=int, default=None,
                        help="the maximum number of independent hits per contig. For nucleotide contig-level hits, each non-overlapping hit region is treated as an independent hit. [default: unlimited]")


    add_taxonomy_filters(parser)


    parser.add_argument("--partial", type=str, default="include", choices=PARTIAL_CHOICES, help=domainate.PARTIAL_HELP)
    parser.add_argument('--decoys', type=str, default=None, nargs='+',
                        help="Names of decoy domains. A decoy domain is a domain that is not expected to be found in the input sequences. If the best hit for a CDS/protein is a domain from the decoys list, it will be not be returned as a hit.")
    parser.add_argument('--decoys_file', type=str, default=None,
                        help="text file containing names of decoy domains.")


    #TODO: add max_domains option

    parser.add_argument('--cpu', type=int, default=0,
                        help="the number of cores of the cpu which are used at a time to run the search (0: one thread per physical core)")
    parser.add_argument('--max_overlap', type=float, default=1,
                        help="the maximum fractional overlap between domains to be included in the annotated genbank. If >= 1, then no overlap filtering will be done.")
    parser.add_argument('--add_annotations', action='store_true', default=False,
                        help="Add domain annotations (Domainator features) for retained hits, in addition to the search best-hit annotations (Domain_Search features) already present.")
    add_dedup(parser)
    parser.add_argument('--max_region_overlap', type=float, default=1.0,
                       help="the maximum fractional of overlap between any two output regions. If >= 1, then no overlap filtering will be done. Regions are output in a greedy fashion based on CDS start site. New regions are output if less than this fraction of them overlaps with any previously output region.")


    parser.add_argument('--translate', action='store_true', default=False,
                            help="by default, nucleotide inputs return nucleotide regions. When --translate is set, CDS hits are translated before writing. Whole-contig nucleic-acid hits are not translated. Note that --translate is incompatble with neighborhood extraction.")


    ### Select by CDS options ###

    parser.add_argument('--kb_range', type=float, default=None,
                        help="How many kb upstream and downstream of the selected CDS to extract. Starting from the ends of the CDS.")
    parser.add_argument('--kb_up', type=float, default=None,
                        help="How many kb upstream of the selected CDS to extract. Starting from the end of the CDS.")
    parser.add_argument('--kb_down', type=float, default=None,
                        help="How many kb downstream of the selected CDS to extract. Starting from the end of the CDS.")

    parser.add_argument('--cds_range', type=int, default=None,
                        help="How many CDSs upstream and downstream of the selected CDS to extract. Example, 1 would mean to return contigs including the selected CDSs and the immediate upstream and downstream CDSs.")
    parser.add_argument('--cds_up', type=int, default=None,
                            help="How many CDSs upstream of the selected CDS to extract. Example, 1 would mean to return contigs including the selected CDSs and the immediate upstream CDS.")
    parser.add_argument('--cds_down', type=int, default=None,
                        help="How many CDSs downstream of the selected CDS to extract. Example, 1 would mean to return contigs including the selected CDSs and the immediate downstream CDS.")

    parser.add_argument('--whole_contig', default=False, action="store_true",
                        help="If set, then return the entire matching contigs. Circular contigs stay circular, and are rotated so that the hit is in the center.")

    parser.add_argument('--strand', type=str, default=None, choices=["f", "r"], help="Only extract regions around CDSs on the specified strand.")

    parser.add_argument('--pad', action='store_true', default=False,
                       help="If set, then ends of sequences will be padded so that the focus CDS aligns. See --anchor.")
    parser.add_argument('--anchor', type=str, default="start", choices=WHERE_CHOICES,
                       help="Where to align the focus CDSs when --pad is set, relative to the focus CDS in its own orientation. "
                       "start: align the first base of the focus CDSs. end: align their last base. center: align their midpoints.")

    parser.add_argument('--keep_direction', action='store_true', default=False,
                       help="by default extracted regions will be flipped so that the focus cds is on the forward strand. Setting this option will keep the focus cds on whatever strand it started on.")

    parser.add_argument("--batch_size", type=int, default=10000, required=False,
                        help="Approximately how many target sequences to search at one time in a batch. "
                        "This also controls partition size: each parallel task processes about this many CDSs/sequences. "
                        "Larger values mean fewer, bigger partitions (lower per-task overhead, but coarser load balancing); "
                        "smaller values give finer parallelism at the cost of more tasks. The input is scanned once to find "
                        "partition boundaries regardless of this value.")

    parser.add_argument('--gene_call', type=str, default=None, choices=["all", "unannotated"], required=False,
                        help="When activated, new CDS annotations will be added with Prodigal in Metagenomic mode. If 'all', then any existing CDS annotations will be deleted and all contigs will be re-annotated. If 'unannotated', then only contigs without CDS annotations will be annotated. "
                        "-Z 0 can't be used with gene calling, because the number of target sequences isn't known until the genes are called.")

    add_max_output_gb_argument(parser)

    parser.add_argument('--config', action=ActionConfigFile)
    declare_outputs(parser, ["genbank_out", "fasta_out"])
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)

    # parameter validation
    cpus = resolve_cpu(params.cpu)
    check_outputs(parser, params)
    reject_stdin(parser, params.input, "it reads its input more than once (to count and partition it)")

    if params.evalue > 10:
        parser.error("-e/--evalue must be <= 10.")
    if params.Z == 0 and params.gene_call is not None:
        parser.error("-Z 0 can't be used with --gene_call: give the database size explicitly.")
    
    if params.batch_size < 1:
        raise ValueError("batch_size must be > 0")

    # Note: --max_output_gb is measured against UNCOMPRESSED GenBank size. For
    # BGZF-compressed outputs (-o foo.gb.gz) the file on disk will be smaller than
    # this limit (R1).
    max_output_bytes = max_output_gb_to_bytes(params.max_output_gb)
    mitigation_options = ["--max_hits", "--max_hits_per_contig", "--cds_range", "--kb_range"]
    extra_guidance = None
    if params.whole_contig:
        extra_guidance = "Avoid --whole_contig on broad searches unless you expect very few hits."

    if params.max_hits is not None and params.max_hits < 1:
        parser.error("--max_hits must be at least 1. Leave it out to return all hits.")
    max_hits = params.max_hits
        
    whole_contig = params.whole_contig
    normalize_direction = not params.keep_direction
    pad = params.pad
    translate = params.translate

    kb_range, cds_range = select_by_cds.validate_range(params.kb_range, params.kb_up, params.kb_down, params.cds_range, params.cds_up, params.cds_down)

    # Expand each input to its shards: a logical database name (mydb.gb) resolves
    # to its numbered shards (mydb.0.gb, mydb.1.gb, ...) if present, else to itself.
    # Done once here so every downstream consumer sees the real per-shard files.
    input_files = [shard for db in params.input for shard in db_index.resolve_database_shards(db)]

    molecule_types = get_input_molecule_types(input_files, default_molecule_type=params.alphabet)
    reference_groups = domainate.read_references(params.reference, foldseek=None)
    domainate.validate_input_molecule_types(molecule_types, reference_groups)
    if len(molecule_types) > 1 and "protein" in molecule_types and not translate:
        raise ValueError("Inputs are a combination of DNA and protein sequences, but --translate is not set. This would lead to an output genbank file with mixed protein and DNA sequences.")
    if ("protein" in molecule_types or translate) and (kb_range is not None or cds_range is not None or whole_contig or pad):
        raise ValueError("Protein sequence output is not compatible with neighborhood extraction or padding.")
    if params.dedup == "sequence_rc":
        if "protein" in molecule_types or translate:
            parser.error("--dedup sequence_rc is only supported for nucleotide output.")
    


    #partition the input files for parallelization later
    if params.Z == 0: #if Z is none, then we must know the total CDS count for per-sequence E-values.
        # Counted before the taxonomy filter, so filtering does not change E-values.
        # If every input has a valid .didx index, the total count is in the headers
        # (O(1) per file), so we can take the true count cheaply and still partition
        # lazily. Otherwise fall back to an eager parallel count over all files.
        header_counts = [index_total_cds(f) for f in input_files]
        if header_counts and all(c is not None for c in header_counts):
            Z = sum(header_counts)
            partitions = i_partition_seqfiles(input_files, params.batch_size)
        else:
            cds_count = 0
            partitions = list()
            with make_pool(processes=cpus) as pool:
                for file_cds_count, file_partitions in pool.imap_unordered(_partition_seqfile_worker(params.batch_size), input_files):
                    cds_count += file_cds_count
                    partitions.extend(file_partitions)
            Z = cds_count
    else: #Z is pre-set, so partitions can be a lazy iterator.
        partitions = i_partition_seqfiles(input_files, params.batch_size)
        Z = params.Z

    ncbi_taxonomy = taxonomy_from_args(params)

    decoy_names = list_and_file_to_dict_keys(params.decoys, params.decoys_file, as_set=True)

    ### Run
    outputs = SeqOutputs(params.genbank_out, params.fasta_out, max_output_bytes=max_output_bytes,
                         mitigation_options=mitigation_options, extra_guidance=extra_guidance)
    with outputs:
        hits = domain_search(
        partitions,
        references=params.reference,
        z=Z,
        evalue=params.evalue,
        max_hits=max_hits,
        max_overlap=params.max_overlap,
        cpu=cpus,
        batch_size=params.batch_size,
        add_annotations=params.add_annotations,
        cds_range=cds_range,
        kb_range=kb_range,
        whole_contig=whole_contig,
        normalize_direction=normalize_direction,
        translate=translate,
        gene_call=params.gene_call,
        min_evalue=params.min_evalue,
        ncbi_taxonomy=ncbi_taxonomy,
        include_taxids=params.include_taxids,
        exclude_taxids=params.exclude_taxids,
        taxonomy_expr=params.taxonomy_expr,
        taxonomy_unknown=params.taxonomy_unknown,
        alphabet=params.alphabet,
        max_mode=params.max_mode,
        max_region_overlap=params.max_region_overlap,
        strand=params.strand,
        decoys=decoy_names,
        max_hits_per_contig=params.max_hits_per_contig,
        partial=params.partial,
        dedup=params.dedup,
        )
        outputs.write_all(select_by_cds.deduplicate_and_pad(hits, pad=pad, anchor=params.anchor, best_first=True))
