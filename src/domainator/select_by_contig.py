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
OR (contig contains at least one domain with name in specified list [domains, domains_file])
OR (contig contains a combination of domains that match the specified boolean expression [domain_expr])
OR (contig contains no domains [unannotated])
OR (contig contains a qualifying search hit [has_search_hit])
OR (contig has a search-hit score within specified bounds [min_search_score, max_search_score])
)
AND
(total number of selected contigs is less than or equal to specified maximum [max_records])

If --invert is specified, the selection criteria in the middle expression (the stuff linked by OR) are inverted. Length, CDS count, taxonomy, fragment status, and max_records still apply,
but the contigs that would have been selected by the middle expression are instead discarded, and contigs that would have been discarded are selected.

"""
from jsonargparse import ActionConfigFile
from jsonargparse.typing import PositiveFloat
import sys

from domainator.cli.common import add_annotation_type, add_has_search_hit, add_input, add_seq_outputs, check_outputs, declare_outputs, tool_parser
from domainator.outputs import write_seq_outputs
from domainator.seqio import parse_seqfiles
from domainator.cli.common import add_list_and_file, list_and_file_to_dict_keys
from domainator.expr import BooleanEvaluator
from domainator.taxonomy.records import compile_taxonomy_filter
from domainator.cli.common import add_taxonomy_filters, taxonomy_from_args
from domainator.partial import filter_by_partial, PARTIAL_CHOICES
from typing import Set, Optional, Iterable, Iterator
from domainator import DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME
import re
from domainator.taxonomy import NCBITaxonomy
from domainator.model import annotation_feature_types, evalue_passes, feature_score
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



def select_by_contig(records: Iterable[SeqRecord], *, domains: Optional[Set[str]] = None, domain_expr: Optional[str] = None,
                     evalue: Optional[float] = None, max_records: Optional[int] = None, min_length: Optional[int] = None,
                     max_length: Optional[int] = None, min_cds: Optional[int] = None, max_cds: Optional[int] = None,
                     definition_regex: Optional[str] = None, sequence_regex: Optional[str] = None,
                     contigs: Optional[Set[str]] = None, contigs_regex: Optional[str] = None, invert: bool = False,
                     ncbi_taxonomy: Optional[NCBITaxonomy] = None, include_taxids: Optional[Iterable[int]] = None,
                     exclude_taxids: Optional[Iterable[int]] = None, taxonomy_expr: Optional[str] = None,
                     taxonomy_unknown: str = "exclude",
                     databases: Optional[Set[str]] = None, unannotated: bool = False,
                     min_search_score: Optional[float] = None, max_search_score: Optional[float] = None,
                     annotation_type: str = "domainator", partial: str = "include",
                     has_search_hit: bool = False) -> Iterator[SeqRecord]:
    """Yield the records of ``records`` that are selected.

    A record is selected if it matches any of the name, definition, sequence, domain and
    search-hit or search-score criteria (or with none of them, every record), and it must also pass the
    length, CDS-count, taxonomy and ``partial`` filters.

    Args:
        records: the contigs or proteins (SeqRecords).
        domains: select records with a domain of one of these names.
        domain_expr: select records whose domain names satisfy this boolean expression.
        evalue: consider only domains with an e-value at most this (None: no cutoff).
        max_records: yield at most this many records.
        min_length: keep only input records at least this long, in bp for nucleotide input
            or amino acids for protein input. Not affected by invert.
        max_length: keep only input records at most this long, in bp for nucleotide input
            or amino acids for protein input. Not affected by invert.
        min_cds: keep only records with at least this many CDSs.
        max_cds: keep only records with at most this many CDSs.
        definition_regex: select records whose definition matches this regular expression.
        sequence_regex: select records whose sequence matches this regular expression.
        contigs: select records with these names.
        contigs_regex: select records whose names match this regular expression.
        invert: invert the selection criteria (not the filters).
        ncbi_taxonomy: an NCBITaxonomy, needed by the taxonomy filters.
        include_taxids: keep only records with a taxid in this list (or below it).
        exclude_taxids: drop records with a taxid in this list (or below it).
        taxonomy_expr: a boolean expression over taxids (operators & | ~ and parentheses), e.g.
            "2 & ~1224". Mutually exclusive with include_taxids/exclude_taxids.
        taxonomy_unknown: unresolved record taxonomy policy: exclude (default), include,
            evaluate against an empty lineage, or error. Not affected by invert.
        databases: consider only domains from these databases. None: all.
        unannotated: select records with no domain annotations.
        min_search_score: also select records whose existing Domain_Search annotation
            passes evalue and databases and scores at least this. Combined with
            max_search_score as an inclusive range, OR with other selectors, and inverted
            by invert. Also restricts has_search_hit scoring, which ignores databases.
        max_search_score: upper bound of the inclusive search-score selector range.
            Not an independent restriction on records selected by other criteria.
        annotation_type: which annotations count as domains: "domainator", "domain_search", or
            "all".
        has_search_hit: also select records with search hits passing evalue and search-score
            cutoffs, ignoring domain names, expressions, and databases. Participates in invert.
        partial: "include" (the default) keeps all records, "exclude" drops fragments, "only"
            keeps only fragments. A protein record is a fragment if
            partial.get_fragment_status says so (e.g. a UniProt " (Fragment)" header); a
            nucleotide contig is a fragment if it has a CDS with '<' or '>' in its location.

    Yields:
        The selected records.
    """

    feature_types = annotation_feature_types(annotation_type)
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

    taxonomy_filter = compile_taxonomy_filter(ncbi_taxonomy, include_taxids, exclude_taxids,
                                              taxonomy_expr, taxonomy_unknown)
    
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
                
        if domains is not None or domain_expr is not None or unannotated or has_search_hit or min_search_score is not None or max_search_score is not None:
            contig_domains = set()
            domain_search_feature = None
            at_least_one_test = True
            for feature in contig.features:
                if has_search_hit and feature.type == DOMAIN_SEARCH_BEST_HIT_NAME and evalue_passes(feature, evalue):
                    if min_search_score is None or min_search_score <= feature_score(feature) <= max_search_score:
                        selected = True
                if feature.type == DOMAIN_FEATURE_NAME and evalue_passes(feature, evalue) and (databases is None or feature.qualifiers["database"][0] in databases):
                    if feature.type in feature_types:
                        contig_domains.add(feature.qualifiers["name"][0])
                if feature.type == DOMAIN_SEARCH_BEST_HIT_NAME and evalue_passes(feature, evalue) and (databases is None or feature.qualifiers["database"][0] in databases):
                    domain_search_feature = feature
                    if feature.type in feature_types:
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
            if taxonomy_filter is not None and not taxonomy_filter.matches(contig):
                continue
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

    add_list_and_file(parser, "contigs",
                      help="Select records with these IDs, OR with other record selectors; affected by --invert.")
    parser.add_argument('--contigs_regex', type=str, default=None,
                        help="Only consider contigs with names matching this regex.")
    
    parser.add_argument("--definition_regex", default=None, required=False, type=str,
                        help="only keep contigs with a definition that matches this regex.")

    parser.add_argument("--sequence_regex", default=None, required=False, type=str, #TODO: add a test for this.
                        help="only keep contigs with a sequence that matches this regex.")

    add_list_and_file(parser, "domains",
                      help="Select records containing these domain names, OR with other record selectors.")
    parser.add_argument('--domain_expr', type=str, default=None,
                        help="a boolean expression using operators & (AND), | (OR), and ~(NOT), to specify a desired combination of domains")

    parser.add_argument("--partial", type=str, default="include", choices=PARTIAL_CHOICES,
                        help="Select contigs by fragment status. include: all. exclude: drop fragments. only: keep only fragments. "
                        "A protein record is a fragment if it has a UniProt ' (Fragment)' or MGnify 'FL=0' header, a 'Flags: Fragment' definition, NON_TER features, "
                        "or a CDS/Protein feature with '<' or '>' spanning it. A nucleotide contig is a fragment if any of its CDSs has '<' or '>' in its location, or a /partial qualifier. Not affected by --invert.")
    add_annotation_type(parser, help="Which annotations to select contigs by.")
    add_has_search_hit(parser)

    parser.add_argument("--databases", default=None, required=False, type=str, nargs="+",
                        help="Consider only domains from these databases when filtering by domain. default: all databases.")

    parser.add_argument('--unannotated', action='store_true', default=False, help="Select contigs with no domain annotations.")

    parser.add_argument('-e', '--evalue', default=None, type=PositiveFloat,
                        help="the evalue cutoff (max) for domain annotations to be considered when selecting by domain. Default: no cutoff.")

    add_taxonomy_filters(parser)

    #

    parser.add_argument('--max_records', default=None, type=int,
                        help="If set then stop after returning this many contigs.")

    parser.add_argument('--min_length', default=None, type=int,
                        help="Keep only input records with length >= this, in bp for nucleotide input or amino acids for protein input. Not affected by --invert.")

    parser.add_argument('--max_length', default=None, type=int,
                        help="Keep only input records with length <= this, in bp for nucleotide input or amino acids for protein input. Not affected by --invert.")

    parser.add_argument('--min_cds', default=None, type=int,
                        help="skip contigs with fewer CDSs than this. Pseudogenes are not counted.")

    parser.add_argument('--max_cds', default=None, type=int,
                        help="skip contigs with more CDSs than this. Pseudogenes are not counted.")
    
    parser.add_argument('--min_search_score', default=None, type=float,
                            help="Also select records with an existing Domain_Search annotation scoring >= this and passing e-value/database filters. Paired with --max_search_score as an inclusive range; OR with other selectors and affected by --invert. Also restricts --has_search_hit scoring, which ignores databases.")
    parser.add_argument('--max_search_score', default=None, type=float,
                            help="Upper bound (inclusive) of the search-score selector range. OR with other selectors, not a restriction on records selected by other criteria. Also restricts --has_search_hit scoring, which ignores databases.")
    
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

    ncbi_taxonomy = taxonomy_from_args(params)

    # Run
   
    extracted_contigs_iterator = select_by_contig(
                    parse_seqfiles(genbanks), 
                    domains=selected_domains, 
                    domain_expr=domain_expr, 
                    evalue=params.evalue,
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
                    taxonomy_unknown=params.taxonomy_unknown,
                    databases=params.databases,
                    unannotated=params.unannotated,
                    min_search_score=params.min_search_score,
                    max_search_score=params.max_search_score,
                    annotation_type=params.annotation_type,
                    has_search_hit=params.has_search_hit,
                    partial=params.partial
    )

    write_seq_outputs(extracted_contigs_iterator, genbank=params.genbank_out, fasta=params.fasta_out)
