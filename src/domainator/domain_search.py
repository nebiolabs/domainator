"""Search for matches to hmm profiles

For searching large databases of sequences (genbank files or fasta files) with small numbers (< ~100) of queries (profiles or protein sequences).

Returns the contigs with hits, possibly truncated to the neighborhood around the hits.

When --max_hits is specified, it applies to the entire search, not individual queries. Thus if you supply 5 query profiles and --max_hits 100, 
at most 100 hits will be returned (not 500). 

When max_hits or pad are specified, hits will be sorted by best domain score and written after the entire search completes.
If neither is set, then hits will be written on the fly and not sorted.
"""

import sys
from jsonargparse import ArgumentParser, ActionConfigFile
from domainator.utils import parse_seqfiles, write_genbank, list_and_file_to_dict_keys, make_pool
from domainator import __version__
from domainator import select_by_cds
import psutil
from domainator import partition_seqfile
from domainator import domainate, RawAndDefaultsFormatter
import heapq
from typing import List, Tuple, Set, Optional
from domainator import extract_peptides
import warnings
from pathlib import Path
from domainator.Taxonomy import NCBITaxonomy

class _domain_search_worker():
    
    def __init__(self, references: List, z: int , evalue: float, max_overlap: float, add_annotations: bool, cds_range: Tuple, kb_range: Tuple, whole_contig: bool, normalize_direction: bool, translate: bool, gene_call:str = None, min_evalue:float = 0.0, batch_size: int = 10000, ncbi_taxonomy: Optional[NCBITaxonomy] = None, include_taxids: Optional[Set[int]] = None, exclude_taxids: Optional[Set[int]] = None, fasta_type: str = "protein", max_mode: bool = False, max_region_overlap=1.0, strand: Optional[str] = None, decoy_names: Optional[Set[str]] = None):

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
        self.ncbi_taxonomy = ncbi_taxonomy
        self.include_taxids = include_taxids
        self.exclude_taxids = exclude_taxids
        self.fasta_type = fasta_type
        self.max_mode = max_mode
        self.max_region_overlap = max_region_overlap
        self.strand = strand
        self.decoy_names = decoy_names
        #self.pre_parsed_references = domainate.read_references(references) #TODO: figure out why pickling pyhmmer peptides causes multithreading issues, and try to patch pyhmmer.


    

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
                kb_range: extract a contig region enclosing this many kb upstream and downstream of the selected CDSs. Partially enclosed CDSs will not be annotated in the output.
                whole_contig: extract the whole contigs of containing the selected CDSs (if a single contig contains multiple selected CDSs, only one copy of the contig will be returned)
                normalize_direction: if True then if a target cds occurs on the reverse strand, reverse-complement the extracted region before returning
            returns: 
                list of SeqRecords
        """
        out = list()
        
        for rec in domainate.domainate(
                parse_seqfiles((partition[0],), None, filetype_override=None, seek_to=partition[1], max_recs=partition[2], default_molecule_type=self.fasta_type), #TODO: clear best_domain_hit ?
                references=self.references,
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
                ncbi_taxonomy=self.ncbi_taxonomy,
                include_taxids=self.include_taxids,
                exclude_taxids=self.exclude_taxids,
                max_mode=self.max_mode,
                #pre_parsed_references=self.pre_parsed_references
            ):
            try:
                if rec.annotations['molecule_type'] == "protein":
                    if self.decoy_names is not None and rec.get_hit_names()[0] in self.decoy_names:
                        continue
                    rec.set_dist_to_start(0)
                    out.append(rec)
                else: # it's a nucleotide, so pass it through select_by_cds
                    if not self.translate: 
                        for region in select_by_cds.select_by_cds(
                            (rec,),
                            target_cdss=None,
                            target_domains=None,
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
                                out.append(region)
                    else: # self.translate == True, so we need to translate it
                        for peptide in extract_peptides.extract_peptides(
                            (rec,),
                            evalue=100000000,
                            target_domains=None,
                            target_cdss=None,
                            keep_cds_feature=True,
                            strand=self.strand,
                            _from_domain_search=True,
                            _domain_search_negatives=self.decoy_names
                        ):
                            rec.set_dist_to_start(0)
                            out.append(peptide)
            except Exception as e:
                warnings.warn(f"Error processing {rec.id}.")
                raise e
        return out
 
class _partition_seqfile_worker():
    def __init__(self, cdss_per_partition):
        self.cdss_per_partition = cdss_per_partition
    
    def __call__(self, input_file):
        return partition_seqfile.partition_seqfile(input_file,cdss_per_partition=self.cdss_per_partition)
    

def domain_search(partitions, references, z, evalue, max_hits, max_overlap, cpu, add_annotations, cds_range, kb_range, whole_contig, normalize_direction, translate, gene_call=None, min_evalue=0.0, ncbi_taxonomy=None, include_taxids=None, exclude_taxids=None, fasta_type="protein", max_mode:bool=False, max_region_overlap=1.0, strand=None, decoy_names=None):
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
        kb_range: extract a contig region enclosing this many kb upstream and downstream of the selected CDSs. Partially enclosed CDSs will not be annotated in the output.
        whole_contig: extract the whole contigs of containing the selected CDSs (if a single contig contains multiple selected CDSs, only one copy of the contig will be returned)
        normalize_direction: if True then if a target cds occurs on the reverse strand, reverse-complement the extracted region before returning
        translate: if True then translate nucleotide hits into proteins before returning.
        ncbi_taxonomy: an NCBITaxonomy object (needed for filtering by taxonomy, otherwise not needed)
        include_taxids: only keep contigs with a taxonomy id in this list
        exclude_taxids: only keep contigs with a taxonomy id not in this list
        fasta_type: "protein" or "nucleotide", the type of fasta file to be searched.
        max_mode: if True, then only the best hit for each hmm will be returned, otherwise all hits will be returned
        max_region_overlap: the maximum fractional of overlap between any two output regions. If >= 1, then no overlap filtering will be done. Regions are output in a greedy fashion based on CDS start site. New regions are output if less than this fraction of them overlaps with any previously output region. [default 1]
        strand: only extract regions around CDSs on the specified strand. If None, then extract regions around CDSs on both strands.
        decoy_names: a set of names of decoy domains. A decoy domain is a domain that is not expected to be found in the input sequences. If the best hit for a CDS/protein is a domain from the decoys list, it will be not be returned as a hit.
    yields:
        SeqRecords of the selected regions
    """

    if cpu < 2: 
        cpu=2

    if include_taxids is not None:
        include_taxids = set(include_taxids)
    if exclude_taxids is not None:
        exclude_taxids = set(exclude_taxids)

    out_heap = []
    worker = _domain_search_worker(references, z, evalue, max_overlap, add_annotations, cds_range, kb_range, whole_contig, normalize_direction, translate, gene_call, min_evalue, ncbi_taxonomy=ncbi_taxonomy, include_taxids=include_taxids, exclude_taxids=exclude_taxids, fasta_type=fasta_type, max_mode=max_mode, max_region_overlap=max_region_overlap, strand=strand, decoy_names=decoy_names)

    with make_pool(processes=cpu - 1) as pool:
        # hits are lists of SeqRecords
        for hits in pool.imap_unordered(worker, partitions):
            
            for rec in hits:
                if max_hits is None:
                    yield rec
                else:
                    if len(out_heap) < max_hits:
                        heapq.heappush(out_heap, rec) # Domainator seqrecs compare based on DOMAINATOR_HIT_BEST_SCORE_ANNOTATION.
                    else:
                        heapq.heappushpop(out_heap, rec)
    
    if len(out_heap) > 0:
        out_heap.sort(reverse=True)
        for rec in out_heap:
            yield rec

def get_input_molecule_types(input_files: List[str], default_molecule_type="protein") -> Set[str]:
    """reads the first record from each input file and gets the 'molecule_type' annotation

    Args:
        input_files (List[str]): a list of input file paths
        default_molecule_type (str, optional): the molecule type to use if the input file does not have a molecule_type annotation. Defaults to "protein".
        
    
    Returns:
        a set of molecule type names.
    """

    return {x.annotations['molecule_type'] for x in parse_seqfiles(input_files, default_molecule_type=default_molecule_type, max_recs=1)}



def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument('-i', '--input', default=None, required=True,
                       nargs='+', type=str,
                       help="the genbank or fasta files to annotate. Genbank files can be nucleotide (with CDS annotations) or peptide. Fasta files must be peptide.")
    parser.add_argument("--fasta_type", type=str, default="protein", choices={"protein", "nucleotide"}, 
                        help="Whether the sequences in fasta files are protein or nucleotide sequences.")
    parser.add_argument('-r', '--references', required=True, type=str,
                        default=None, nargs='+',
                        help="the names of the HMM files with profiles to search. Or protein query files. ") 

    parser.add_argument('-o', '--output', default=None, type=str, required=False,
                        help="output genbank filename. If not supplied, writes to stdout.")


    parser.add_argument('-Z', default=1000, type=int, 
                        help="Passed as the -Z parameter to hmmsearch: Assert that the total number of targets in your searches is <x>, for the purposes of per-sequence E-value calculations, rather than the actual number of targets seen. "
                        "Default: 1000\n"
                        "Set to 0 to use the number actual target sequences (currently this does does not account for taxonomic filtering)." # TODO: account for taxonomy filtering
                        "Supplying -Z is recommended for most use cases of domainate, as it will speed up the search and make comparisons between searches more meaningful.")

    parser.add_argument('-e', '--evalue', type=float, default=0.001,
                        help="threshold E value for the domain hit. Must be <=10. [default 0.001]")
    
    parser.add_argument('--min_evalue', type=float, default=0.0,
                        help="hits with E value lower than this will be filtered out. Not frequently used. Use only if you want to eliminate close matches. Must be >=0. [default 0]")

    parser.add_argument('--max_mode', action="store_true", default=False,
                        help="Run hmmsearch/phmmer in maximum sensitivity mode, which is much slower, but more sensitive.")

    parser.add_argument('--max_hits', type=int, default=None,
                        help="the maximum number of CDSs returned by the search. Prioritized by bitscore of best scoring profile. [default: return all hits passing the evalue threshold]")
    


    parser.add_argument("--include_taxids", nargs='+', default=None, type=int, help="Space separated list of taxids to include")
    parser.add_argument("--exclude_taxids", nargs='+', default=None, type=int, help="Space separated list of taxids to exclude")
    parser.add_argument("--ncbi_taxonomy_path", type=str,  default="/tmp/ncbi_taxonomy", help="Path to NCBI taxonomy database directory. Will be created and downloaded if it does not exist.")
    parser.add_argument("--taxonomy_update", action="store_true", help="If taxonomy database exists, check it against the version on the ncbi server and update if there is a newer version.")


    parser.add_argument('--decoys', type=str, default=None, nargs='+',
                        help="Names of decoy domains. A decoy domain is a domain that is not expected to be found in the input sequences. If the best hit for a CDS/protein is a domain from the decoys list, it will be not be returned as a hit.")
    parser.add_argument('--decoys_file', type=str, default=None,
                        help="text file containing names of decoy domains.")


    #TODO: add max_domains option

    parser.add_argument('--cpu', type=int, default=0,
                        help="the number of cores of the cpu which are used at a time to run the search [default: use all available cores]")
    parser.add_argument('--max_overlap', type=float, default=1,
                        help="the maximum fractional overlap between domains to be included in the annotated genbank. If >= 1, then no overlap filtering will be done. [default 1]")
    parser.add_argument('--add_annotations', action='store_true', default=False,
                        help="When activated, new domainator annotations will be added to the file, not just the domain_search annotations. Useful if you want to see what the non-best hits score.")
    parser.add_argument('--deduplicate', action='store_true', default=False,
                        help="by default if the same region is extracted for multiple hits then both will be kept. Set this option to eliminate redundancies.")
    parser.add_argument('--max_region_overlap', type=float, default=1.0,
                       help="the maximum fractional of overlap between any two output regions. If >= 1, then no overlap filtering will be done. Regions are output in a greedy fashion based on CDS start site. New regions are output if less than this fraction of them overlaps with any previously output region.")


    parser.add_argument('--translate', action='store_true', default=False,
                            help="by default, nucleotide databases will return nucleotide hits. When --translate is set, CDS hits will be translated before writing. Note that --translate is incompatble with neighborhood extraction.")


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
                        help="If set, then return the entire matching contigs.")

    parser.add_argument('--strand', type=str, default=None, choices=["f", "r"], help="Only extract regions around CDSs on the specified strand.")

    parser.add_argument('--pad', action='store_true', default=False,
                       help="If set, then ends of sequences will be padded so that the focus CDS aligns.")

    parser.add_argument('--keep_direction', action='store_true', default=False,
                       help="by default extracted regions will be flipped so that the focus cds is on the forward strand. Setting this option will keep the focus cds on whatever strand it started on.")

    parser.add_argument("--batch_size", type=int, default=10000, required=False,
                        help="Approximately how many protein sequences to search at one time in a batch.")

    parser.add_argument('--gene_call', type=str, default=None, choices = {"all", "unannotated"}, required=False,
                        help="When activated, new CDS annotations will be added with Prodigal in Metagenomic mode. If 'all', then any existing CDS annotations will be deleted and all contigs will be re-annotated. If 'unannotated', then only contigs without CDS annotations will be annotated. [default: None] "
                        "Note that if you do gene calling, it is STRONGLY recommended that you also supply -Z, because database size is pre-calcuated at the beginning of the execution, whereas gene-calling is done on the fly. Not supplying Z may become an error in the future.")

    parser.add_argument('--config', action=ActionConfigFile)

    params = parser.parse_args(argv)

    # parameter validation
    if params.cpu <= 0:
        cpus = psutil.cpu_count(logical=False)
    else:
        cpus = params.cpu

    if params.evalue > 10:
        raise ValueError("evalue must be <= 10.")
    
    if params.batch_size < 1:
        raise ValueError("batch_size must be > 0")

    if params.output is None:
        out = sys.stdout
    else:
        out = open(params.output, "w")

    if params.max_hits == 0:
        max_hits = None
    else:
        max_hits = params.max_hits
        
    whole_contig = params.whole_contig
    normalize_direction = not params.keep_direction
    deduplicate = params.deduplicate
    pad = params.pad
    translate = params.translate

    kb_range, cds_range = select_by_cds.validate_range(params.kb_range, params.kb_up, params.kb_down, params.cds_range, params.cds_up, params.cds_down)

    molecule_types = get_input_molecule_types(params.input, default_molecule_type=params.fasta_type)
    if len(molecule_types) > 1 and "protein" in molecule_types and not translate:
        raise ValueError("Inputs are a combination of DNA and protein sequences, but --translate is not set. This would lead to an output genbank file with mixed protein and DNA sequences.")
    if ("protein" in molecule_types or translate) and (kb_range is not None or cds_range is not None or whole_contig or pad):
        raise ValueError("Protein sequence output is not compatible with neighborhood extraction or padding.")
    


    #partition the input files for parallelization later
    if params.Z == 0: #if Z is none, then we must read through the entire input to count CDSs and get file offsets. Luckily we can do this in parallel if there are multiple inputs.
        #TODO: this is broken for taxonomy filtering, because it doesn't filter for taxonomy. 
        cds_count = 0
        partitions = list()
        with make_pool(processes=cpus) as pool:
            for file_cds_count, file_partitions in pool.imap_unordered(_partition_seqfile_worker(params.batch_size), params.input):
                cds_count += file_cds_count
                partitions.extend(file_partitions)
        Z = cds_count
    else: #Z is pre-set, so partitions can be a lazy iterator.
        partitions = partition_seqfile.i_partition_seqfiles(params.input, cdss_per_partition=params.batch_size)
        Z = params.Z

    ncbi_taxonomy = None
    if params.include_taxids or params.exclude_taxids:
        # create the path to the NCBI taxonomy database
        Path(params.ncbi_taxonomy_path).mkdir(parents=True, exist_ok=True)
        # load the NCBI taxonomy database
        ncbi_taxonomy = NCBITaxonomy(params.ncbi_taxonomy_path, params.taxonomy_update)

    decoy_names = list_and_file_to_dict_keys(params.decoys, params.decoys_file, as_set=True)

    ### Run
    seen = set()
    records = list()
    for record in domain_search(
        partitions,
        references=params.references,
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
        fasta_type=params.fasta_type,
        max_mode=params.max_mode,
        max_region_overlap=params.max_region_overlap,
        strand=params.strand,
        decoy_names=decoy_names
        ):
        #skip pad
        if not pad:
            if not deduplicate or (record.id not in seen):
                write_genbank((record,), out)
        else:
            if not deduplicate or (record.id not in seen):
                records.append(record)
        if deduplicate:
            seen.add(record.id)
    if pad:
        records.sort(reverse=True)
        select_by_cds.pad_records(records)
        write_genbank(records, out)

    if params.output is not None:
        out.close()

def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])
