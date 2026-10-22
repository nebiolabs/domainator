"""Search for matches to protein or nucleic-acid references.

For searching large databases of sequences (genbank files or fasta files) with small numbers (< ~100) of queries.
Protein references are searched against CDSs or protein records.
Nucleic-acid references are searched against whole nucleotide contigs.
Supported references include protein HMMs, protein FASTA, nucleotide FASTA, nucleotide HMMs, and infernal CM files.
Protein inputs cannot be searched with nucleotide references.

Returns the contigs with hits, possibly truncated to the neighborhood around the hits.

When --max_hits is specified, it applies to the entire search, not individual queries. Thus if you supply 5 query profiles and --max_hits 100, 
at most 100 hits will be returned (not 500). 

When max_hits or pad are specified, hits will be sorted by best domain score and written after the entire search completes.
If neither is set, then hits will be written on the fly and not sorted.
"""

from jsonargparse import ActionConfigFile
from jsonargparse.typing import NonNegativeFloat, PositiveFloat
from domainator.cli.common import add_alphabet, add_input, add_seq_outputs, add_switch, check_outputs, declare_outputs, reject_stdin, resolve_cpu, tool_parser
from domainator.outputs import SeqOutputs
from domainator.seqio import parse_seqfiles, native_parser_available
from domainator.offsets import partition_seqfile, i_partition_seqfiles, index_total_cds
from domainator.cli.common import list_and_file_to_dict_keys
from domainator._mp import make_pool
from domainator.taxonomy.records import compile_taxonomy_allowed
from domainator.partial import PARTIAL_CHOICES
from domainator.locations import WHERE_CHOICES
from domainator.model import sequence_hashes, get_hit_best_score, get_hit_names, set_dist_to_start
from domainator import select_by_cds
from domainator import db_index
from domainator import domainate
import heapq
from typing import List, Tuple, Set, Optional, Iterable, Iterator
from domainator import extract_peptides
from pathlib import Path
from domainator.taxonomy import NCBITaxonomy, default_ncbi_taxonomy_path
from domainator.output_guardrails import add_max_output_gb_argument, max_output_gb_to_bytes
from Bio.SeqRecord import SeqRecord


# Per-process taxonomy filter, populated by the pool initializer. The include/
# exclude filter is "compiled" once in the parent into a frozenset of allowed
# taxids (NCBITaxonomy.compile_taxonomy_filter); filtering is then a single O(1)
# membership test and the set is pickled once per process, not per partition.
_WORKER_ALLOWED_TAXIDS = None

# Per-process parsed references, populated by the pool initializer. References are
# parsed once per worker (from path strings passed through initargs) rather than
# re-parsed for every partition. Only path strings cross the spawn boundary, so no
# pyhmmer/pyinfernal object is ever pickled (see R4 in speed_up_plan.md).
_WORKER_REFERENCES = None


def _is_new_sequence(record, seen: Set[bytes], both_strands: bool) -> bool:
    """Exact-sequence dedup check: True (and record its hashes in seen) if record's sequence hasn't been seen yet."""
    hashes = sequence_hashes(record.seq, both_strands)
    if any(h in seen for h in hashes):
        return False
    seen.update(hashes)
    return True


def _init_search_worker(reference_paths, foldseek_paths, allowed_taxids):
    """Pool initializer: parse references once per worker, install the taxonomy filter.

    Args:
        reference_paths: list of reference file paths (parsed inside the worker).
        foldseek_paths: list of foldseek database paths, or None.
        allowed_taxids: frozenset of taxids that pass the taxonomy filter, or None
            if no taxonomy filtering is active.
    """
    global _WORKER_REFERENCES, _WORKER_ALLOWED_TAXIDS
    _WORKER_REFERENCES = domainate.read_references(reference_paths, foldseek=foldseek_paths)
    _WORKER_ALLOWED_TAXIDS = allowed_taxids


class _domain_search_worker():

    def __init__(self, references: List, z: int , evalue: float, max_overlap: float, add_annotations: bool, cds_range: Tuple, kb_range: Tuple, whole_contig: bool, normalize_direction: bool, translate: bool, gene_call:str = None, min_evalue:float = 0.0, batch_size: int = 10000, include_taxids: Optional[Set[int]] = None, exclude_taxids: Optional[Set[int]] = None, fasta_type: str = "protein", max_mode: bool = False, max_region_overlap=1.0, strand: Optional[str] = None, decoy_names: Optional[Set[str]] = None, max_hits_per_contig: Optional[int] = None, max_hits: Optional[int] = None, partial: str = "include", deduplicate_by_sequence: bool = False, dedup_both_strands: bool = False):

        self.z = z
        self.evalue = evalue
        self.max_overlap = max_overlap
        self.add_annotations = add_annotations
        self.cds_range = cds_range
        self.kb_range = kb_range
        self.whole_contig = whole_contig
        self.normalize_direction = normalize_direction
        self.batch_size = batch_size
        self.references = references
        self.translate = translate
        self.gene_call = gene_call
        self.min_evalue = min_evalue
        # Kept for the constructor contract; the actual filter is the compiled
        # frozenset of allowed taxids installed per worker process (_WORKER_ALLOWED_TAXIDS).
        self.include_taxids = include_taxids
        self.exclude_taxids = exclude_taxids
        self.fasta_type = fasta_type
        self.max_mode = max_mode
        self.max_region_overlap = max_region_overlap
        self.strand = strand
        self.decoy_names = decoy_names
        self.max_hits_per_contig = max_hits_per_contig
        self.partial = partial
        # When max_hits is set, each worker keeps only its own top-max_hits records
        # (by best domain score) before returning, so the parent receives at most
        # O(workers * max_hits) records to merge instead of every hit. The global
        # top-max_hits is unchanged: a single worker can contribute at most max_hits
        # records to a global top-max_hits, so trimming per worker is lossless.
        self.max_hits = max_hits
        # Sequence dedup is also done per partition so that, with max_hits set, the
        # per-worker top-max_hits isn't crowded out by copies of one sequence (which
        # would make the per-worker trim lossy). domain_search dedups again globally.
        self.deduplicate_by_sequence = deduplicate_by_sequence
        self.dedup_both_strands = dedup_both_strands
        # References are no longer parsed per-partition here; the pool initializer
        # (_init_search_worker) parses them once per worker into _WORKER_REFERENCES.

    def __call__(self, partition):
        """
            a single thread worker for domain_search

            input:
                partition: (sequence_file_path, file_offset, records_to_read)
                references: A list of hmm or protein fasta file paths, the input HMM databases to be searched against the sequences using hmmsearch or phmmer.
                z: Z parameter to pass to hmmsearch/phmmer
                evalue: The threshold E value for hmmer hit
                max_overlap: The maximum fractional of overlap to be allowed between domains
                add_annotations: When activated, new domainator annotations will be added to the file, not just the domain_search annotations. Useful if you want to see what the non-best hits score.
                cds_range: extract a contig region enclosing this many CDSs upstream and downstream of the CDS hits
                kb_range: extract a contig region enclosing this many kb upstream and downstream of the selected CDSs. Partially enclosed features are truncated, with '<' or '>' marking the cut ends.
                whole_contig: extract the whole contigs of containing the selected CDSs (if a single contig contains multiple selected CDSs, only one copy of the contig will be returned)
                normalize_direction: if True then if a target cds occurs on the reverse strand, reverse-complement the extracted region before returning
            returns:
                list of SeqRecords. When max_hits is set, at most max_hits records
                (the highest scoring for this partition) are returned.
        """
        out = list()
        max_hits = self.max_hits
        seen = set() if self.deduplicate_by_sequence else None
        arrivals = []

        def _emit(record):
            if seen is not None and not _is_new_sequence(record, seen, self.dedup_both_strands):
                return
            # Bounded min-heap over this worker's own hits: with max_hits set the
            # worker never holds (or returns) more than max_hits records. Entries are
            # (best score, -arrival, record), so the heap root is always this partition's
            # weakest kept hit, and of equal scores the earlier hit ranks higher.
            if max_hits is None:
                out.append(record)
                return
            entry = (get_hit_best_score(record), -len(arrivals), record)
            arrivals.append(None)
            if len(out) < max_hits:
                heapq.heappush(out, entry)
            else:
                heapq.heappushpop(out, entry)

        for rec in domainate.domainate(
                parse_seqfiles((partition[0],), None, filetype_override=None, seek_to=partition[1], max_recs=partition[2], default_molecule_type=self.fasta_type, genbank_parser="lean"), #TODO: clear best_domain_hit ?
                references=None,  # references are parsed once per worker by _init_search_worker (R4)
                pre_parsed_references=_WORKER_REFERENCES,
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
                allowed_taxids=_WORKER_ALLOWED_TAXIDS,  # per-CDS / per-region pre-hmmer taxonomy filter
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
        if max_hits is None:
            return out
        return [record for _, _, record in sorted(out, reverse=True)]  # best first, then in input order
 
class _partition_seqfile_worker():
    def __init__(self, cdss_per_partition):
        self.cdss_per_partition = cdss_per_partition
    
    def __call__(self, input_file):
        return partition_seqfile(input_file, self.cdss_per_partition)
    

def domain_search(partitions: Iterable[Tuple[str, int, int]], *, references: List[str], z: int, evalue: float,
                  max_hits: Optional[int], max_overlap: float, cpu: int, add_annotations: bool,
                  cds_range: Optional[Tuple[int, int]], kb_range: Optional[Tuple[float, float]], whole_contig: bool,
                  normalize_direction: bool, translate: bool, gene_call: Optional[str] = None, min_evalue: float = 0.0,
                  ncbi_taxonomy: Optional[NCBITaxonomy] = None, include_taxids: Optional[Iterable[int]] = None,
                  exclude_taxids: Optional[Iterable[int]] = None, taxonomy_expr: Optional[str] = None,
                  alphabet: str = "protein", max_mode: bool = False, max_region_overlap: float = 1.0,
                  strand: Optional[str] = None, decoys: Optional[Set[str]] = None, max_hits_per_contig: Optional[int] = None,
                  partial: str = "include", deduplicate_by_sequence: bool = False, both_strands: bool = False) -> Iterator[SeqRecord]:
    """
    runs hmmsearch in parallel on multiple sections of genbank or fasta files. 

    input:
        partitions: list of tuples (sequence_file_path, file_offset, records_to_read), each tuple will be passed to the threadpool and hmmsearch and select_by_cds will be run on it.
        references: A list of hmm or protein fasta file paths, the input HMM databases to be searched against the sequences using hmmsearch or phmmer.
        z: Z parameter to pass to hmmsearch/phmmer
        evalue: The threshold E value for hmmer hit
        max_hits: The maximum number of CDSs to be returned
        max_overlap: The maximum fractional of overlap to be allowed between domains
        cpu: number of threads to use. Must be at least 2.
        add_annotations: if True, then add new domainate annotations .
        cds_range: extract a contig region enclosing this many CDSs upstream and downstream of the CDS hits
        kb_range: extract a contig region enclosing this many kb upstream and downstream of the selected CDSs. Partially enclosed features are truncated, with '<' or '>' marking the cut ends.
        whole_contig: extract the whole contigs of containing the selected CDSs (if a single contig contains multiple selected CDSs, only one copy of the contig will be returned)
        normalize_direction: if True then if a target cds occurs on the reverse strand, reverse-complement the extracted region before returning
        translate: if True then translate nucleotide hits into proteins before returning.
        ncbi_taxonomy: an NCBITaxonomy object (needed for filtering by taxonomy, otherwise not needed)
        include_taxids: only keep contigs with a taxonomy id in this list
        exclude_taxids: only keep contigs with a taxonomy id not in this list
        taxonomy_expr: a boolean expression over taxids (operators & | ~ and parentheses) used to filter contigs by lineage, e.g. "2 & ~1224". Mutually exclusive with include_taxids/exclude_taxids.
        alphabet: "protein" or "nucleotide", the type of fasta file to be searched.
        max_mode: if True, then only the best hit for each hmm will be returned, otherwise all hits will be returned
        max_region_overlap: the maximum fractional of overlap between any two output regions. If >= 1, then no overlap filtering will be done. Regions are output in a greedy fashion based on CDS start site. New regions are output if less than this fraction of them overlaps with any previously output region. [default 1]
        strand: only extract regions around CDSs on the specified strand. If None, then extract regions around CDSs on both strands.
        decoys: a set of names of decoy domains. A decoy domain is a domain that is not expected to be found in the input sequences. If the best hit for a CDS/protein is a domain from the decoys list, it will be not be returned as a hit.
        max_hits_per_contig: the maximum number of independent hits per contig. If None, then all hits will be returned. For nucleotide contig-level hits, each non-overlapping hit region is treated as an independent hit.
        partial: which proteins/CDSs to search, by fragment status: "include" (all, the default), "exclude" (skip fragments), or "only" (only fragments). See domainate.domainate.
        deduplicate_by_sequence: if True, drop hits whose output sequence is identical (case-insensitive, ignoring terminal '*') to an earlier hit.
            The first copy to arrive is kept; partitions finish in arbitrary order, so which of several identical copies is kept is arbitrary.
            Only a set of 32-byte hashes is held, not the duplicate records.
        both_strands: with deduplicate_by_sequence, also treat complement and reverse-complement nucleotide sequences as duplicates.
    yields:
        SeqRecords of the selected regions
    """

    if cpu < 2: 
        cpu=2

    if include_taxids is not None:
        include_taxids = set(include_taxids)
    if exclude_taxids is not None:
        exclude_taxids = set(exclude_taxids)

    if taxonomy_expr and (include_taxids or exclude_taxids):
        raise ValueError("taxonomy_expr is mutually exclusive with include_taxids/exclude_taxids.")

    # "Compile" the taxonomy filter once in the parent into a frozenset of allowed
    # taxids (a single O(N) pass over the tree), then hand it to each worker process
    # via the pool initializer. Filtering in workers is then an O(1) membership test,
    # and the set -- not the full NCBITaxonomy -- is pickled once per process. The
    # full NCBITaxonomy is released afterward so the parent does not hold it.
    allowed_taxids = None
    if ncbi_taxonomy is not None:
        allowed_taxids = compile_taxonomy_allowed(ncbi_taxonomy, include_taxids, exclude_taxids, taxonomy_expr)
        ncbi_taxonomy = None

    # Always use the pool initializer: it parses references once per worker (R4) and,
    # when present, installs the compiled taxonomy filter.
    pool_kwargs = dict(initializer=_init_search_worker, initargs=(references, None, allowed_taxids))

    out_heap = []
    arrivals = 0
    worker = _domain_search_worker(references, z, evalue, max_overlap, add_annotations, cds_range, kb_range, whole_contig, normalize_direction, translate, gene_call, min_evalue, include_taxids=include_taxids, exclude_taxids=exclude_taxids, fasta_type=alphabet, max_mode=max_mode, max_region_overlap=max_region_overlap, strand=strand, decoy_names=decoys, max_hits_per_contig=max_hits_per_contig, max_hits=max_hits, partial=partial, deduplicate_by_sequence=deduplicate_by_sequence, dedup_both_strands=both_strands)
    seen = set() if deduplicate_by_sequence else None

    # Pool tuning notes (speed_up_plan.md Q4):
    #   - No maxtasksperchild: workers must persist so the references parsed once
    #     by _init_search_worker (R4) stay loaded. Recycling workers would re-parse
    #     references on every restart, defeating R4.
    #   - No chunksize on imap_unordered: partitions are already large (~batch_size
    #     CDSs each), so per-task dispatch overhead is negligible; chunking would
    #     only coarsen load balancing. (domainator hmmer_search reaches the same conclusion.)
    #   - cpu - 1 workers: the parent process does the heap merge and GenBank
    #     writing, so it is reserved one core rather than oversubscribing.
    with make_pool(processes=cpu - 1, **pool_kwargs) as pool:
        # hits are lists of SeqRecords. When max_hits is set, each list is already
        # trimmed to that worker's own top-max_hits (see _domain_search_worker), so
        # the parent merges O(workers * max_hits) records rather than every hit; the
        # central heap below still selects the global top-max_hits.
        for hits in pool.imap_unordered(worker, partitions):

            for rec in hits:
                if seen is not None and not _is_new_sequence(rec, seen, both_strands):
                    continue
                if max_hits is None:
                    yield rec
                else:
                    # (best score, -arrival, record): the root is the weakest kept hit
                    entry = (get_hit_best_score(rec), -arrivals, rec)
                    arrivals += 1
                    if len(out_heap) < max_hits:
                        heapq.heappush(out_heap, entry)
                    else:
                        heapq.heappushpop(out_heap, entry)

    out_heap.sort(reverse=True)
    for _, _, rec in out_heap:
        yield rec

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
                        help="reference files to search with. Supported types are protein HMMs, protein FASTA, nucleotide FASTA, nucleotide HMMs, and infernal CM files.") 

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


    parser.add_argument("--include_taxids", nargs='+', default=None, type=int, help="Space separated list of taxids to include")
    parser.add_argument("--exclude_taxids", nargs='+', default=None, type=int, help="Space separated list of taxids to exclude")
    parser.add_argument("--taxonomy_expr", type=str, default=None, help="A boolean expression over taxids using operators & (AND), | (OR), ~ (NOT), and parentheses, e.g. \"2 & ~1224\" (within Bacteria but not Proteobacteria). A taxid is true for a contig when it is in the contig's lineage. Mutually exclusive with --include_taxids/--exclude_taxids.")
    parser.add_argument("--ncbi_taxonomy_path", type=str,  default=default_ncbi_taxonomy_path(), help="Path to NCBI taxonomy database directory. Will be created and downloaded if it does not exist.")
    parser.add_argument("--taxonomy_update", action="store_true", help="If taxonomy database exists, check it against the version on the ncbi server and update if there is a newer version.")


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
                        help="When activated, new domainator annotations will be added to the file, not just the domain_search annotations. Useful if you want to see what the non-best hits score.")
    add_switch(parser, "deduplicate", default=False,
               help="if the same region is extracted for multiple hits, then keep only one copy of it (compared by record ID).",
               negative_help="keep every copy of a region extracted for multiple hits. This can result in duplicated contig names.")
    parser.add_argument('--deduplicate_by_sequence', action='store_true', default=False,
                        help="Drop hits whose output sequence is identical to an earlier hit (SHA-256 hash, as in domainator deduplicate_genbank --algorithm hash). "
                        "Useful for highly redundant databases. Independent of --deduplicate, which only drops repeated record IDs.")
    parser.add_argument('--both_strands', action='store_true', default=False,
                        help="With --deduplicate_by_sequence, also treat complement and reverse-complement nucleotide sequences as duplicates.")
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
    deduplicate = params.deduplicate
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
    if params.both_strands:
        if not params.deduplicate_by_sequence:
            parser.error("--both_strands requires --deduplicate_by_sequence.")
        if "protein" in molecule_types or translate:
            parser.error("--both_strands is only supported for nucleotide output.")
    


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

    if params.taxonomy_expr and (params.include_taxids or params.exclude_taxids):
        raise ValueError("--taxonomy_expr is mutually exclusive with --include_taxids/--exclude_taxids.")

    ncbi_taxonomy = None
    if params.include_taxids or params.exclude_taxids or params.taxonomy_expr:
        # create the path to the NCBI taxonomy database
        Path(params.ncbi_taxonomy_path).mkdir(parents=True, exist_ok=True)
        # load the NCBI taxonomy database
        ncbi_taxonomy = NCBITaxonomy(params.ncbi_taxonomy_path, params.taxonomy_update)

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
        alphabet=params.alphabet,
        max_mode=params.max_mode,
        max_region_overlap=params.max_region_overlap,
        strand=params.strand,
        decoys=decoy_names,
        max_hits_per_contig=params.max_hits_per_contig,
        partial=params.partial,
        deduplicate_by_sequence=params.deduplicate_by_sequence,
        both_strands=params.both_strands,
        )
        outputs.write_all(select_by_cds.deduplicate_and_pad(hits, deduplicate=deduplicate, pad=pad, anchor=params.anchor, best_first=True))
