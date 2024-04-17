"""Create a new genbank file with a subset of contigs from the input

From a domainator-annotated genbank file, extracts a subset of contigs based on selection criteria.

Selection criteria are combined as follows:

(contig length within specified bounds [length_lb, length_ub])
AND 
(contig matches specified taxonomy filters)
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
from jsonargparse import ArgumentParser, ActionConfigFile
import sys

from domainator.utils import parse_seqfiles, write_genbank, list_and_file_to_dict_keys, BooleanEvaluator, filter_by_taxonomy
from typing import Set, Optional
from domainator import __version__, DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME, RawAndDefaultsFormatter
import re
from pathlib import Path
from domainator.Taxonomy import NCBITaxonomy
from domainator.Bio import SeqIO

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



def select_by_contig(contigs, target_domains: Set[str] = None, domain_expr: str =None, domain_evalue=float("inf"), first=None, length_lb=None, length_ub=None, definition_regex=None, sequence_regex=None, target_contigs=None, contigs_regex=None, invert=False, ncbi_taxonomy=None, include_taxids=None, exclude_taxids=None, databases:Optional[Set[str]]=None, unanntotated=False, search_score_lb=None, search_score_ub=None, domain_type="domain"):
    """

        Args:
            contigs: an iterable/iterator of SeqRecords

            target_domains: Select all contigs containing domains with these names

            domain_expr: Use boolean expressions for domain names to select CDSs

            definition_regex: only keep contigs with a definition that matches this regex.

            sequence_regex: only keep contigs with a sequence that matches this regex.

            target_contigs: only keep contigs with these names

            contigs_regex: only keep contigs with names matching this regex.

            domain_evalue: consider only domains with evalues less than this

            first: only keep the first n contigs

            length_lb: only keep contigs with length greater than this

            length_ub: only keep contigs with length less than this

            invert: invert the selection criteria

            ncbi_taxonomy: an NCBITaxonomy object (needed for filtering by taxonomy, otherwise not needed)

            include_taxids: only keep contigs with a taxonomy id in this list

            exclude_taxids: only keep contigs with a taxonomy id not in this list

            databases: only consider domains from these databases

            unanntotated: keep contigs with no domain annotations

            search_score_lb: only keep contigs with a domain_search annotation score greater than this

            search_score_ub: only keep contigs with a domain_search annotation score less than this

            domain_type: the type of domain to consider, domain, search, both. Default: domain

        Yields:
            SeqRecords of the selected contigs

    """

    databases = set(databases) if databases is not None else None
    returned = 0

    if target_domains is not None:
        target_domains = set(target_domains)

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

    if include_taxids or exclude_taxids:
        contigs = filter_by_taxonomy(contigs, include_taxids, exclude_taxids, ncbi_taxonomy)
    
    if length_lb is not None or length_ub is not None:
        contigs = length_filter(contigs, length_lb, length_ub)

    if search_score_lb is not None or search_score_ub is not None:
        if search_score_lb is None:
            search_score_lb = float("-inf")
        if search_score_ub is None:
            search_score_ub = float("inf")

    for contig in contigs:
        selected = False
        at_least_one_test = False # TODO: make this more modular, so we don't need this variable.
        if target_contigs is not None:
            at_least_one_test = True
            if contig.id in target_contigs:
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
                
        if target_domains is not None or domain_expr is not None or unanntotated or search_score_lb is not None or search_score_ub is not None:
            domains = set()
            domain_search_feature = None
            at_least_one_test = True
            for feature in contig.features:
                if feature.type == DOMAIN_FEATURE_NAME and float(feature.qualifiers['evalue'][0]) < domain_evalue and (databases is None or feature.qualifiers["database"][0] in databases):
                    if domain_type == "domain" or domain_type == "both":
                        domains.add(feature.qualifiers["name"][0])
                if feature.type == DOMAIN_SEARCH_BEST_HIT_NAME and float(feature.qualifiers['evalue'][0]) < domain_evalue and (databases is None or feature.qualifiers["database"][0] in databases):
                    domain_search_feature = feature
                    if domain_type == "search" or domain_type == "both":
                        domains.add(feature.qualifiers["name"][0])

            if search_score_lb is not None and domain_search_feature is not None:
                score = float(domain_search_feature.qualifiers["score"][0])
                if score >= search_score_lb and score <= search_score_ub:
                    selected = True
            
            if target_domains is not None and len(target_domains.intersection(domains)) > 0: # select the contig if it contains any of the target domains
                selected = True
            
            if domain_expr is not None and expression_evaluator.check_expression(domains): # select the contig if it contains the specified combination of domains
                selected = True
            
            if unanntotated and len(domains) == 0: # select the contig if it contains no domains
                selected = True
        if not at_least_one_test: # if no selection criteria are specified, select the contig
            selected = True
        
        if invert:
            selected = not selected
        
        if selected:
            if first is None:
                yield contig
            else:
                yield contig
                returned += 1
                if returned == first:
                    break

def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument("-i", "--input", default=None, required=False,
                          nargs="+", type=str,
                          help="names of input genbank files. If not supplied, reads from stdin.")

    parser.add_argument('-o', '--output', default=None, required=False,  
                        help="genbank output file name. If not supplied writes to stdout.")

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

    parser.add_argument("--domain_type", type=str, default="domain", choices={"domain", "search", "both"}, help="The type of domain to consider, domain, search, both. Default: domain")    

    parser.add_argument("--databases", default=None, required=False, type=str, nargs="+",
                        help="Consider only domains from these databases when filtering by domain. default: all databases.")

    parser.add_argument('--unannotated', action='store_true', default=False, help="Select contigs with no domain annotations.")

    parser.add_argument('-e', '--evalue', default=100000000.0, type=float,
                        help="the evalue cutoff for domain annotations to be considered when selecting by domain.")

    parser.add_argument('--fasta_out', action='store_true', default=False,
                        help="makes output a fasta file when activated")

    parser.add_argument("--include_taxids", nargs='+', default=None, type=int, help="Space separated list of taxids to include")
    parser.add_argument("--exclude_taxids", nargs='+', default=None, type=int, help="Space separated list of taxids to exclude")
    parser.add_argument("--ncbi_taxonomy_path", type=str,  default="/tmp/ncbi_taxonomy", help="Path to NCBI taxonomy database directory. Will be created and downloaded if it does not exist.")
    parser.add_argument("--taxonomy_update", action="store_true", help="If taxonomy database exists, check it against the version on the ncbi server and update if there is a newer version.")


    # parser.add_argument('--taxid', default=None, type=int, nargs="+",  #TODO: implement this. Probably need another argument to point to the ncbi taxonomy database files.
    #                     help="Exclude contigs based on the /db_xref='taxon:' qualifier of the source annotation. "
    #                     "A contig will be excluded iff the contig taxid is not a sub taxon of a supplied taxon, or if the contig lacks a taxon annotation. More than one taxid can be supplied.")
    #
    # parser.add_argument('--qualifier_regex', type=str, nargs=3, required=False, action=CustomAction, #TODO: implement this
    #                     help="Supply three strings, an feature type and a qualifier name and a regex. Keep only contigs that match.")

    parser.add_argument('--first', default=None, type=int,
                        help="If set then stop after returning this many contigs.")

    parser.add_argument('--length_lb', default=None, type=int,
                        help="skip contigs smaller than this. In units of bp (aa for protein input).")

    parser.add_argument('--length_ub', default=None, type=int,
                        help="skip contigs larger than this. In units of bp (aa for protein input).")
    
    parser.add_argument('--search_score_lb', default=None, type=float,
                            help="Skip contigs with a domain_search annotation score smaller than this.")
    parser.add_argument('--search_score_ub', default=None, type=float,
                            help="Skip contigs with a domain_search annotation score larger than this.")
    
    parser.add_argument('--invert', action='store_true',
                        help="Invert the selection criteria. i.e. return contigs that do not match the selection criteria.")


    parser.add_argument('--config', action=ActionConfigFile)  
    

    
    # parser.add_argument('-d', '--database', type=str, default=None,
    #                    help="the name of the database which you want to filter annotations from.") #TODO: see todo above.

    #TODO: error when evalue is defined but no domains are defined.

    params = parser.parse_args(argv)


    # figure out selection filters

    target_contigs = list_and_file_to_dict_keys(
        params.contigs, params.contigs_file, as_set=True)
    selected_domains = list_and_file_to_dict_keys(
        params.domains, params.domains_file, as_set=True)
    domain_expr = params.domain_expr

    filetype_override = None
    if params.input is None:
        genbanks = [sys.stdin]
        filetype_override = "genbank"
    else:
        genbanks = params.input

    if params.output is None:
        output_handle = sys.stdout
    else:
        output_handle = open(params.output, "w")

    ncbi_taxonomy = None
    if params.include_taxids or params.exclude_taxids:
        # create the path to the NCBI taxonomy database
        Path(params.ncbi_taxonomy_path).mkdir(parents=True, exist_ok=True)
        # load the NCBI taxonomy database
        ncbi_taxonomy = NCBITaxonomy(params.ncbi_taxonomy_path, params.taxonomy_update)

    # Run
   
    extracted_contigs_iterator = select_by_contig(
                    parse_seqfiles(genbanks, filetype_override=filetype_override), 
                    selected_domains, 
                    domain_expr, 
                    params.evalue,
                    params.first,
                    params.length_lb,
                    params.length_ub,
                    params.definition_regex,
                    params.sequence_regex,
                    target_contigs,
                    params.contigs_regex,
                    params.invert,
                    ncbi_taxonomy,
                    params.include_taxids,
                    params.exclude_taxids,
                    databases=params.databases,
                    unanntotated=params.unannotated,
                    search_score_lb=params.search_score_lb,
                    search_score_ub=params.search_score_ub,
                    domain_type=params.domain_type
    )

    if params.fasta_out:
        SeqIO.write(extracted_contigs_iterator, output_handle, "fasta")
    else:
        write_genbank(extracted_contigs_iterator, output_handle)

    if params.output is not None:
        output_handle.close()


def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])
