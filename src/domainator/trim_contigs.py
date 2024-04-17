"""Trim contigs based on CDS or domain names or presence of domains in CDSs.

Circular contigs are linearized before trimming.

Typically used for nucleotide sequences, but can be used for protein sequences as well, if you use the 'kb' options to mean 1000 amino acids. CDS and domain options are not supported for protein sequences.

"""

from jsonargparse import ArgumentParser, ActionConfigFile
import sys
from domainator.utils import parse_seqfiles, slice_record, write_genbank, list_and_file_to_dict_keys, BooleanEvaluator, DomainatorCDS, get_sources, copy_feature, pad_location, slice_record_from_location, FeatureLocation, CompoundLocation, get_non_domainator_features, circular_dist
from domainator.Bio.SeqRecord import SeqRecord
from typing import Tuple, List, Optional, Set, Iterable, Union
from domainator import __version__, RawAndDefaultsFormatter, DOMAIN_SEARCH_BEST_HIT_NAME
from domainator.select_by_cds import validate_range

class StripDomainsError(Exception):
    """
    Custom exception for when strip_domains_from_ends encounters a CDS that is entirely composed of domains to be stripped.
    """
    pass

def strip_domains_from_ends(cds_list:List[DomainatorCDS], left_cds, right_cds, strip_domains: Optional[Set[str]], domain_expr:Optional[Union[str, BooleanEvaluator]], no_domain:bool=False) -> SeqRecord:

    if strip_domains is None:
        strip_domains = set()



    if domain_expr is not None:
        if type(domain_expr) is not BooleanEvaluator:
            domain_expr = BooleanEvaluator(domain_expr)

    if left_cds > right_cds:
        raise ValueError("left_cds must be <= right_cds")

    if left_cds < 0 or right_cds >= len(cds_list):
        raise ValueError("left_cds and right_cds must be >= 0 and < len(cds_list)")

    if len(strip_domains) == 0 and domain_expr is None:
        return left_cds, right_cds
    
    while (left_cds < len(cds_list) and 
           (len(cds_list[left_cds].domain_names.intersection(strip_domains)) > 0  
            or (domain_expr is not None and domain_expr.check_expression(cds_list[left_cds].domain_names)) 
            or (no_domain and len(cds_list[left_cds].domain_names) == 0)
           )
        ):
        left_cds += 1
    
    while (right_cds >= 0 and 
            (len(cds_list[right_cds].domain_names.intersection(strip_domains)) > 0 
             or (domain_expr is not None and domain_expr.check_expression(cds_list[right_cds].domain_names)) 
             or (no_domain and len(cds_list[right_cds].domain_names) == 0)
            )
        ):
        right_cds -= 1

    if left_cds > right_cds:
        raise StripDomainsError("All CDSs in the contig are composed of domains to be stripped.")
    
    return left_cds, right_cds


def trim_contig(record:SeqRecord, kb_truncation:Optional[Tuple[float, float]], cds_truncation:Optional[Tuple[int, int]], strip_domains:Optional[Set[str]], domain_expr:Optional[str]=None, no_domain:bool=False) -> SeqRecord:
    if kb_truncation is None and cds_truncation is None and strip_domains is None and domain_expr is None:
        return record # TODO: should this be a copy?

    if domain_expr is not None:
        if type(domain_expr) is not BooleanEvaluator:
            domain_expr = BooleanEvaluator(domain_expr)
    
    left = 0
    right = len(record)
    if kb_truncation is not None:
        left += int(kb_truncation[0] * 1000)
        right -= int(kb_truncation[1] * 1000)
    
    cdss = DomainatorCDS.list_from_contig(record)
    cdss.sort(key=lambda x: x.feature.location.start)

    left_cds = 0
    right_cds = len(cdss) - 1    
    if len(cdss) > 0 and strip_domains is not None or domain_expr is not None or cds_truncation is not None or no_domain:
        if cds_truncation is not None:
            left_cds = cds_truncation[0]
            right_cds = (len(cdss) - 1) - cds_truncation[1]

        if strip_domains is not None or domain_expr is not None or no_domain:
            try:
                domain_left_cds, domain_right_cds = strip_domains_from_ends(cdss, left_cds, right_cds, strip_domains, domain_expr, no_domain)
                if domain_left_cds > left_cds:
                    left_cds = domain_left_cds
                if domain_right_cds < right_cds:
                    right_cds = domain_right_cds
            except StripDomainsError:
                print(f"WARNING: All CDSs in contig {record.id} are composed of domains to be stripped. Skipping contig.", file=sys.stderr)
                return None

        if cdss[left_cds].feature.location.start > cdss[right_cds].feature.location.end:
            print(f"WARNING: CDSs in contig {record.id} overlap after trimming. Skipping contig.", file=sys.stderr)
            return None
        
        if (left_cds > 0) and (cdss[left_cds].feature.location.start > left):
            left = cdss[left_cds].feature.location.start
        if (right_cds < len(cdss) - 1 ) and (cdss[right_cds].feature.location.end < right):
            right = cdss[right_cds].feature.location.end
    if left >= right:
        return None

    if left > 0 or right < len(record):
        record = slice_record(record, left, right)
        record.id = record.id + f"_{left+1}:{right}"
    else:
        record = record #TODO: should this be a copy?
    
    return record


def trim_contigs(records: Iterable[SeqRecord], kb_truncation: Optional[Tuple[float, float]], cds_truncation: Optional[Tuple[int, int]], strip_domains: Optional[Set[str]], domain_expr:Optional[str]=None, no_domain:bool=False) -> List[SeqRecord]:
    if strip_domains is not None:
        strip_domains = set(strip_domains)

    if domain_expr is not None:
        domain_expr = BooleanEvaluator(domain_expr)

    for record in records:
        trimmed = trim_contig(record, kb_truncation, cds_truncation, strip_domains, domain_expr, no_domain)
        if trimmed is not None:
            yield trimmed


def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument("-i", "--input", default=None, nargs="+", type=str, required=False,
                        help="names of input genbank files. If not supplied, reads from stdin.")

    parser.add_argument('-o', '--output', default=None, required=False, 
                        help="genbank output file name. If not supplied, writes to stdout.")


    parser.add_argument('--kb_both', type=float, default=None,
                        help="How many kb to trim from both sides of each contig.")
    parser.add_argument('--kb_left', type=float, default=None,
                        help="How many kb to trim from the left side of each contig.")
    parser.add_argument('--kb_right', type=float, default=None,
                        help="How many kb to trim from the right side of each contigs.")

    parser.add_argument('--cds_both', type=int, default=None,
                        help="How many CDSs to trim from both sides of each contig.")
    parser.add_argument('--cds_left', type=int, default=None,
                        help="How many CDSs to trim from the left side of each contig.")
    parser.add_argument('--cds_right', type=int, default=None,
                        help="How many CDSs to trim from the right side of each contig.")


    parser.add_argument('--contigs', type=str, default=None, nargs='+',
                        help="Only return contigs with these names.")
    parser.add_argument('--contigs_file', type=str, default=None,
                        help="text file containing names of contigs to return.")
    parser.add_argument('--domains', type=str, default=None, nargs='+',
                        help="Trim CDSs containing domains with these names from the ends of each contig.") #TODO: add some way to specify the database of a domain
    parser.add_argument('--no_domain', action='store_true', default=False,
                        help="Trim CDSs containing no domain annotations from the ends of each contig.")
    
    parser.add_argument('--domains_file', type=str, default=None, 
                        help="text file containing names of domains to trim contigs containing.")
    parser.add_argument('--domain_expr', type=str, default=None,
                        help="a boolean expression using operators & (AND), | (OR), and ~(NOT), to specify a desired combination of domains for CDSs to trim")
    
    parser.add_argument('--config', action=ActionConfigFile)
    
    params = parser.parse_args(argv)

    ### validate input


    filetype_override = None
    if params.input is None:
        genbanks = (sys.stdin,)
        filetype_override = "genbank"
    else:
        genbanks = params.input

    if params.output is None:
        output_handle = sys.stdout
    else:
        output_handle = open(params.output, "w")

    ### figure out selection filters
    kb_truncation, cds_truncation = validate_range(params.kb_both, params.kb_left, params.kb_right, params.cds_both, params.cds_left, params.cds_right)

    target_contigs = list_and_file_to_dict_keys(params.contigs, params.contigs_file, as_set=True)
    selected_domains = list_and_file_to_dict_keys(params.domains, params.domains_file, as_set=True)
    domain_expr = params.domain_expr


    ### Run
    #TODO: resolve code duplication with domain_search.py
    records = list()
    for record in trim_contigs(
        parse_seqfiles(genbanks, target_contigs, filetype_override=filetype_override),
        kb_truncation,
        cds_truncation,
        selected_domains,
        domain_expr,
        params.no_domain,
        ):
        write_genbank((record,), output_handle)
 
 
    if params.output is not None:
        output_handle.close()

def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])
  