"""Create a new genbank file with a subset of contigs from the input

From a domainator-annotated genbank file, extracts a subset of contigs based on selection criteria.

Selection criteria are combined as follows:

(contig length within specified bounds [--min_length, --max_length])
AND 
(contig CDS count within specified bounds [--min_cds, --max_cds])
AND 
(contig matches specified taxonomy filters)
AND
(contig passes the fragment filter [partial])
AND # INVERT applies to this list of criteria
(
   (contig name in specified list [contigs, contigs_file])
OR (contig name matches specified regex [contigs_regex])
OR (contig definition matches specified regex [definition_regex]) 
OR (contig sequence matches specified regex [sequence_regex])
OR (contig contains at least one domain with name in specified list [target_domains, target_domains_file]) 
OR (contig contains a combination of domains that match the specified boolean expression [domain_expr])
OR (contig contains no domains [unannotated])
)
AND
(total number of selected contigs is less than or equal to specified maximum [first])

If --invert is specified, the selection criteria in the middle expression (the stuff linked by OR) are inverted. That is, length and first are still applied, 
but the contigs that would have been selected by the middle expression are instead discarded, and contigs that would have been discarded are selected.

"""
from jsonargparse import ActionConfigFile
from jsonargparse.typing import PositiveFloat
import sys

from domainator.cli.common import add_annotation_type, add_input, add_seq_outputs, check_outputs, declare_outputs, tool_parser
from domainator.outputs import write_seq_outputs
from domainator.seqio import parse_seqfiles
from domainator.cli.common import list_and_file_to_dict_keys
from domainator.expr import BooleanEvaluator
from domainator.taxonomy.records import filter_by_taxonomy
from domainator.partial import filter_by_partial, PARTIAL_CHOICES
from typing import Set, Optional, Iterable, Iterator
from domainator import DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME
import re
from pathlib import Path
from domainator.taxonomy import NCBITaxonomy, default_ncbi_taxonomy_path
from domainator.model import feature_evalue, feature_score
from Bio.SeqRecord import SeqRecord

# TODO: account for different domain databases
# TODO: account for domain order, like by making boolean searches binned. Try to emulate the behavior of the InterPro domain selector tool:
#       https://www.ebi.ac.uk/interpro/search/ida/
# TODO: add evalue or bitscore threshold for domains? If evalue is normalized, it should be ok.
# TODO: streaming for genbank files, converting genbank files to sqlite databases is quite inefficient because database building is presumably nlog(n), and it makes piping impossible (unless the db-building functions are re-written to not use seek operations)
#       possibly add streaming and piping in v 1.5, or something like that?
# TODO: when streaming have option to list non-selected contigs
# TODO: add option for random subset or top x contigs

# TODO: add option for selecting by seq.annotations, and seq.feature.qualifiers

def length_filter(contigs, length_lb=None, length_ub=None):
    if length_lb is None:
        length_lb = 0
    if length_ub is None:
        length_ub = sys.maxsize

    for contig in contigs:
        contig_size = len(contig)
        if contig_size < length_lb or contig_size > length_ub:
            continue
        yield contig


def count_cds(contig):
    count = 0
    for feature in contig.features:
        if feature.type == 'CDS':
            if "pseudo" in feature.qualifiers or "pseudogene" in feature.qualifiers:
                continue
            count += 1
    return count


def cds_count_filter(contigs, cds_count_lb=None, cds_count_ub=None):
    if cds_count_lb is None:
        cds_count_lb = 0
    if cds_count_ub is None:
        cds_count_ub = sys.maxsize

    for contig in contigs:
        cds_count = count_cds(contig)
        if cds_count < cds_count_lb or cds_count > cds_count_ub:
            continue
        yield contig



def select_by_contig(records: Iterable[SeqRecord], *, domains: Set[str] = None, domain_expr: str = None,
                     evalue: float = float("inf"), max_records: Optional[int] = None, min_length: Optional[int] = None,
                     max_length: Optional[int] = None, min_cds: Optional[int] = None, max_cds: Optional[int] = None,
                     definition_regex: Optional[str] = None, sequence_regex: Optional[str] = None,
                     contigs: Optional[Set[str]] = None, contigs_regex: Optional[str] = None, invert: bool = False,
                     ncbi_taxonomy: Optional[NCBITaxonomy] = None, include_taxids: Optional[Iterable[int]] = None,
                     exclude_taxids: Optional[Iterable[int]] = None, taxonomy_expr: Optional[str] = None,
                     databases: Optional[Set[str]] = None, unannotated: bool = False,
                     min_search_score: Optional[float] = None, max_search_score: Optional[float] = None,
                     annotation_type: str = "domainator", partial: str = "include") -> Iterator[SeqRecord]:
    """

        Args:
            records: an iterable/iterator of SeqRecords

            domains: Select all records containing domains with these names

            domain_expr: Use boolean expressions for domain names to select CDSs

            definition_regex: only keep records with a definition that matches this regex.

            sequence_regex: only keep records with a sequence that matches this regex.

            contigs: only keep records with these names

            contigs_regex: only keep records with names matching this regex.

            evalue: consider only domains with evalues less than this

            max_records: only keep the first max_records records

            min_length: only keep records with length greater than this

            max_length: only keep records with length less than this

            min_cds: only keep records with CDS count greater than or equal to this

            max_cds: only keep records with CDS count less than or equal to this

            invert: invert the selection criteria

            ncbi_taxonomy: an NCBITaxonomy object (needed for filtering by taxonomy, otherwise not needed)

            include_taxids: only keep records with a taxonomy id in this list

            exclude_taxids: only keep records with a taxonomy id not in this list

            taxonomy_expr: a boolean expression over taxids (operators & | ~ and parentheses), e.g. "2 & ~1224". Mutually exclusive with include_taxids/exclude_taxids.

            databases: only consider domains from these databases

            unannotated: keep records with no domain annotations

            min_search_score: only keep records with a domain_search annotation score greater than this

            max_search_score: only keep records with a domain_search annotation score less than this

            annotation_type: which annotations to consider: "domainator", "domain_search", or "all".
            partial: "include" (default) keeps all records, "exclude" drops fragments, "only" keeps only fragments.
                A protein record is a fragment if utils.get_fragment_status says so (e.g. a UniProt " (Fragment)" header);
                a nucleotide contig is a fragment if it contains a CDS with '<' or '>' in its location.

        Yields:
            SeqRecords of the selected records

    """

    databases = set(databases) if databases is not None else None
    returned = 0

    if domains is not None:
        domains = set(domains)

    if domain_expr is not None:
        expression_evaluator = BooleanEvaluator(domain_expr)
    
    if definition_regex is not None:
        definition_regex = re.compile(definition_regex)

    if sequence_regex is not None:
        sequence_regex = re.compile(sequence_regex)

    if contigs_regex is not None:
        contigs_regex = re.compile(contigs_regex)

    if include_taxids is not None:
        include_taxids = set(include_taxids)
    if exclude_taxids is not None:
        exclude_taxids = set(exclude_taxids)

    if include_taxids or exclude_taxids or taxonomy_expr:
        records = filter_by_taxonomy(records, include_taxids, exclude_taxids, ncbi_taxonomy, taxonomy_expr=taxonomy_expr)
    
    if partial != "include":
        records = filter_by_partial(records, partial)

    if min_length is not None or max_length is not None:
        records = length_filter(records, min_length, max_length)

    if min_cds is not None or max_cds is not None:
        records = cds_count_filter(records, min_cds, max_cds)

    if min_search_score is not None or max_search_score is not None:
        if min_search_score is None:
            min_search_score = float("-inf")
        if max_search_score is None:
            max_search_score = float("inf")

    for contig in records:
        selected = False
        at_least_one_test = False # TODO: make this more modular, so we don't need this variable.
        if contigs is not None:
            at_least_one_test = True
            if contig.id in contigs:
                selected = True

        if contigs_regex is not None:
            at_least_one_test = True
            if contigs_regex.search(contig.id):  # select the contig if the name matches the regex
                selected = True
                
        if definition_regex is not None:
            at_least_one_test = True
            if definition_regex.search(contig.description): # select the contig if the definition matches the regex
                selected = True
                
        if sequence_regex is not None:
            at_least_one_test = True
            if sequence_regex.search(str(contig.seq)): # select the contig if the sequence matches the regex
                selected = True
                
        if domains is not None or domain_expr is not None or unannotated or min_search_score is not None or max_search_score is not None:
            contig_domains = set()
            domain_search_feature = None
            at_least_one_test = True
            for feature in contig.features:
                if feature.type == DOMAIN_FEATURE_NAME and feature_evalue(feature) <= evalue and (databases is None or feature.qualifiers["database"][0] in databases):
                    if annotation_type in ("domainator", "all"):
                        contig_domains.add(feature.qualifiers["name"][0])
                if feature.type == DOMAIN_SEARCH_BEST_HIT_NAME and feature_evalue(feature) <= evalue and (databases is None or feature.qualifiers["database"][0] in databases):
                    domain_search_feature = feature
                    if annotation_type in ("domain_search", "all"):
                        contig_domains.add(feature.qualifiers["name"][0])

            if min_search_score is not None and domain_search_feature is not None:
                score = feature_score(domain_search_feature)
                if score >= min_search_score and score <= max_search_score:
                    selected = True
            
            if domains is not None and len(domains.intersection(contig_domains)) > 0: # select the contig if it contains any of the target domains
                selected = True
            
            if domain_expr is not None and expression_evaluator.check_expression(contig_domains): # select the contig if it contains the specified combination of domains
                selected = True
            
            if unannotated and len(contig_domains) == 0: # select the contig if it contains no domains
                selected = True
        if not at_least_one_test: # if no selection criteria are specified, select the contig
            selected = True
        
        if invert:
            selected = not selected
        
        if selected:
            if max_records is None:
                yield contig
            else:
                yield contig
                returned += 1
                if returned == max_records:
                    break

def build_parser():
    parser = tool_parser(__doc__)

    add_input(parser, 'Names of input genbank files.')

    add_seq_outputs(parser)

    parser.add_argument('--contigs', type=str, default=None, nargs='+',
                        help="Only consider contigs with these names.")
    parser.add_argument('--contigs_file', type=str, default=None,
                        help="text file containing names of contigs to select.")
    parser.add_argument('--contigs_regex', type=str, default=None,
                        help="Only consider contigs with names matching this regex.")
    
    parser.add_argument("--definition_regex", default=None, required=False, type=str,
                        help="only keep contigs with a definition that matches this regex.")

    parser.add_argument("--sequence_regex", default=None, required=False, type=str, #TODO: add a test for this.
                        help="only keep contigs with a sequence that matches this regex.")

    parser.add_argument('--domains', type=str, default=None, nargs='+',
                        help="Include domains with these names")  # TODO: add some way to specify the database of a domain
    parser.add_argument('--domains_file', type=str, default=None,
                        help="text file containing names of domains to include.")
    parser.add_argument('--domain_expr', type=str, default=None,
                        help="a boolean expression using operators & (AND), | (OR), and ~(NOT), to specify a desired combination of domains")

    parser.add_argument("--partial", type=str, default="include", choices=PARTIAL_CHOICES,
                        help="Select contigs by fragment status. include: all. exclude: drop fragments. only: keep only fragments. "
                        "A protein record is a fragment if it has a UniProt ' (Fragment)' or MGnify 'FL=0' header, a 'Flags: Fragment' definition, NON_TER features, "
                        "or a CDS/Protein feature with '<' or '>' spanning it. A nucleotide contig is a fragment if any of its CDSs has '<' or '>' in its location, or a /partial qualifier. Not affected by --invert.")
    add_annotation_type(parser, choices=("domainator", "domain_search", "all"), help="Which annotations to select contigs by.")    

    parser.add_argument("--databases", default=None, required=False, type=str, nargs="+",
                        help="Consider only domains from these databases when filtering by domain. default: all databases.")

    parser.add_argument('--unannotated', action='store_true', default=False, help="Select contigs with no domain annotations.")

    parser.add_argument('-e', '--evalue', default=None, type=PositiveFloat,
                        help="the evalue cutoff (max) for domain annotations to be considered when selecting by domain. Default: no cutoff.")

    parser.add_argument("--include_taxids", nargs='+', default=None, type=int, help="Space separated list of taxids to include")
    parser.add_argument("--exclude_taxids", nargs='+', default=None, type=int, help="Space separated list of taxids to exclude")
    parser.add_argument("--taxonomy_expr", type=str, default=None, help="A boolean expression over taxids using operators & (AND), | (OR), ~ (NOT), and parentheses, e.g. \"2 & ~1224\" (within Bacteria but not Proteobacteria). A taxid is true for a contig when it is in the contig's lineage. Mutually exclusive with --include_taxids/--exclude_taxids.")
    parser.add_argument("--ncbi_taxonomy_path", type=str,  default=default_ncbi_taxonomy_path(), help="Path to NCBI taxonomy database directory. Will be created and downloaded if it does not exist.")
    parser.add_argument("--taxonomy_update", action="store_true", help="If taxonomy database exists, check it against the version on the ncbi server and update if there is a newer version.")

    #

    parser.add_argument('--max_records', default=None, type=int,
                        help="If set then stop after returning this many contigs.")

    parser.add_argument('--min_length', default=None, type=int,
                        help="skip contigs smaller than this. In units of bp (aa for protein input).")

    parser.add_argument('--max_length', default=None, type=int,
                        help="skip contigs larger than this. In units of bp (aa for protein input).")

    parser.add_argument('--min_cds', default=None, type=int,
                        help="skip contigs with fewer CDSs than this. Pseudogenes are not counted.")

    parser.add_argument('--max_cds', default=None, type=int,
                        help="skip contigs with more CDSs than this. Pseudogenes are not counted.")
    
    parser.add_argument('--min_search_score', default=None, type=float,
                            help="Skip contigs with a domain_search annotation score smaller than this.")
    parser.add_argument('--max_search_score', default=None, type=float,
                            help="Skip contigs with a domain_search annotation score larger than this.")
    
    parser.add_argument('--invert', action='store_true',
                        help="Invert the selection criteria. i.e. return contigs that do not match the selection criteria.")

    parser.add_argument('--config', action=ActionConfigFile)  
    

    
    #TODO: error when evalue is defined but no domains are defined.
    declare_outputs(parser, ["genbank_out", "fasta_out"])
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)

    # figure out selection filters

    target_contigs = list_and_file_to_dict_keys(
        params.contigs, params.contigs_file, as_set=True)
    selected_domains = list_and_file_to_dict_keys(
        params.domains, params.domains_file, as_set=True)
    domain_expr = params.domain_expr

    genbanks = params.input

    check_outputs(parser, params)

    if params.taxonomy_expr and (params.include_taxids or params.exclude_taxids):
        raise ValueError("--taxonomy_expr is mutually exclusive with --include_taxids/--exclude_taxids.")

    ncbi_taxonomy = None
    if params.include_taxids or params.exclude_taxids or params.taxonomy_expr:
        # create the path to the NCBI taxonomy database
        Path(params.ncbi_taxonomy_path).mkdir(parents=True, exist_ok=True)
        # load the NCBI taxonomy database
        ncbi_taxonomy = NCBITaxonomy(params.ncbi_taxonomy_path, params.taxonomy_update)

    # Run
   
    extracted_contigs_iterator = select_by_contig(
                    parse_seqfiles(genbanks), 
                    domains=selected_domains, 
                    domain_expr=domain_expr, 
                    evalue=float("inf") if params.evalue is None else params.evalue,
                    max_records=params.max_records,
                    min_length=params.min_length,
                    max_length=params.max_length,
                    min_cds=params.min_cds,
                    max_cds=params.max_cds,
                    definition_regex=params.definition_regex,
                    sequence_regex=params.sequence_regex,
                    contigs=target_contigs,
                    contigs_regex=params.contigs_regex,
                    invert=params.invert,
                    ncbi_taxonomy=ncbi_taxonomy,
                    include_taxids=params.include_taxids,
                    exclude_taxids=params.exclude_taxids,
                    taxonomy_expr=params.taxonomy_expr,
                    databases=params.databases,
                    unannotated=params.unannotated,
                    min_search_score=params.min_search_score,
                    max_search_score=params.max_search_score,
                    annotation_type=params.annotation_type,
                    partial=params.partial
    )

    write_seq_outputs(extracted_contigs_iterator, genbank=params.genbank_out, fasta=params.fasta_out)
