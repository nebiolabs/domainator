"""Trim contigs based on CDS or domain names or presence of domains in CDSs

Circular contigs are linearized before trimming.

Typically used for nucleotide sequences, but can be used for protein sequences as well, if you use the 'kb' options to mean 1000 amino acids. CDS and domain options are not supported for protein sequences.

"""

from jsonargparse import ActionConfigFile
from jsonargparse.typing import PositiveFloat
import logging
from domainator.cli.common import add_input, add_seq_outputs, check_outputs, declare_outputs, tool_parser
from domainator.outputs import write_seq_outputs
from domainator.seqio import parse_seqfiles
from domainator.locations import slice_record
from domainator.cli.common import list_and_file_to_dict_keys
from domainator.expr import BooleanEvaluator
from domainator.model import DomainatorCDS, region_id
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord
from typing import Iterable, Iterator, List, Optional, Set, Tuple, Union
from domainator.select_by_cds import validate_range

logger = logging.getLogger(__name__)

class StripDomainsError(Exception):
    """
    Custom exception for when strip_domains_from_ends encounters a CDS that is entirely composed of domains to be stripped.
    """
    pass

def strip_domains_from_ends(cds_domain_names: List[Set[str]], left_cds, right_cds, strip_domains: Optional[Set[str]], domain_expr:Optional[Union[str, BooleanEvaluator]], no_domain:bool=False) -> Tuple[int, int]:
    """Move left_cds right and right_cds left past CDSs to strip, given each CDS's domain names in contig order."""
    if strip_domains is None:
        strip_domains = set()

    if domain_expr is not None:
        if type(domain_expr) is not BooleanEvaluator:
            domain_expr = BooleanEvaluator(domain_expr)

    if left_cds > right_cds:
        raise ValueError("left_cds must be <= right_cds")

    if left_cds < 0 or right_cds >= len(cds_domain_names):
        raise ValueError("left_cds and right_cds must be >= 0 and < len(cds_domain_names)")

    if len(strip_domains) == 0 and domain_expr is None and not no_domain:
        return left_cds, right_cds

    def strip(names):
        return (len(names.intersection(strip_domains)) > 0
                or (domain_expr is not None and domain_expr.check_expression(names))
                or (no_domain and len(names) == 0))

    while left_cds < len(cds_domain_names) and strip(cds_domain_names[left_cds]):
        left_cds += 1

    while right_cds >= 0 and strip(cds_domain_names[right_cds]):
        right_cds -= 1

    if left_cds > right_cds:
        raise StripDomainsError("All CDSs in the contig are composed of domains to be stripped.")

    return left_cds, right_cds


def trim_contig(record:SeqRecord, kb_truncation:Optional[Tuple[float, float]], cds_truncation:Optional[Tuple[int, int]], strip_domains:Optional[Set[str]], domain_expr:Optional[str]=None, no_domain:bool=False, evalue:Optional[float]=None, databases:Optional[Set[str]]=None) -> SeqRecord:
    """Trim one contig. Domain criteria only count domain annotations from databases (if given) with E-values at most evalue (if given)."""
    by_domain = strip_domains is not None or domain_expr is not None or no_domain
    if kb_truncation is None and cds_truncation is None and not by_domain:
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
    if len(cdss) > 0 and (by_domain or cds_truncation is not None):
        if cds_truncation is not None:
            left_cds = cds_truncation[0]
            right_cds = (len(cdss) - 1) - cds_truncation[1]
            if left_cds > right_cds:
                return None

        if by_domain:
            if evalue is None and databases is None:
                cds_domain_names = [cds.domain_names for cds in cdss]
            else:
                cds_domain_names = [cds.get_domain_names(databases, evalue) for cds in cdss]
            try:
                domain_left_cds, domain_right_cds = strip_domains_from_ends(cds_domain_names, left_cds, right_cds, strip_domains, domain_expr, no_domain)
                if domain_left_cds > left_cds:
                    left_cds = domain_left_cds
                if domain_right_cds < right_cds:
                    right_cds = domain_right_cds
            except StripDomainsError:
                logger.warning(f"All CDSs in contig {record.id} are composed of domains to be stripped. Skipping contig.")
                return None

        if cdss[left_cds].feature.location.start > cdss[right_cds].feature.location.end:
            logger.warning(f"CDSs in contig {record.id} overlap after trimming. Skipping contig.")
            return None

        if (left_cds > 0) and (cdss[left_cds].feature.location.start > left):
            left = cdss[left_cds].feature.location.start
        if (right_cds < len(cdss) - 1 ) and (cdss[right_cds].feature.location.end < right):
            right = cdss[right_cds].feature.location.end
    if left >= right:
        return None

    if left > 0 or right < len(record):
        new_id = region_id(record, FeatureLocation(left, right))
        record = slice_record(record, left, right)
        record.id = new_id

    return record


def trim_contigs(records: Iterable[SeqRecord], *, kb_truncation: Optional[Tuple[float, float]], cds_truncation: Optional[Tuple[int, int]], domains: Optional[Set[str]], domain_expr:Optional[str]=None, unannotated:bool=False, evalue:Optional[float]=None, databases:Optional[Set[str]]=None) -> Iterator[SeqRecord]:
    if domains is not None:
        domains = set(domains)

    if domain_expr is not None:
        domain_expr = BooleanEvaluator(domain_expr)

    for record in records:
        trimmed = trim_contig(record, kb_truncation, cds_truncation, domains, domain_expr, unannotated, evalue, databases)
        if trimmed is not None:
            yield trimmed


def build_parser():
    parser = tool_parser(__doc__)

    add_input(parser, "GenBank files.")

    add_seq_outputs(parser)


    parser.add_argument('--kb_range', type=float, default=None,
                        help="How many kb to trim from both sides of each contig.")
    parser.add_argument('--kb_left', type=float, default=None,
                        help="How many kb to trim from the left side of each contig.")
    parser.add_argument('--kb_right', type=float, default=None,
                        help="How many kb to trim from the right side of each contigs.")

    parser.add_argument('--cds_range', type=int, default=None,
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
    parser.add_argument('--unannotated', action='store_true', default=False,
                        help="Trim CDSs containing no domain annotations from the ends of each contig.")
    
    parser.add_argument('--domains_file', type=str, default=None, 
                        help="text file containing names of domains to trim contigs containing.")
    parser.add_argument('--domain_expr', type=str, default=None,
                        help="a boolean expression using operators & (AND), | (OR), and ~(NOT), to specify a desired combination of domains for CDSs to trim")
    parser.add_argument('-e', '--evalue', type=PositiveFloat, default=None,
                        help="the evalue cutoff (max) for domain annotations to be considered when trimming by domain. Default: no cutoff.")
    parser.add_argument('--databases', type=str, default=None, nargs='+',
                        help="Only consider domain annotations from these databases when trimming by domain.")
    
    parser.add_argument('--config', action=ActionConfigFile)
    declare_outputs(parser, ["genbank_out", "fasta_out"])
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)

    ### validate input


    check_outputs(parser, params)

    ### figure out selection filters
    kb_truncation, cds_truncation = validate_range(params.kb_range, params.kb_left, params.kb_right, params.cds_range, params.cds_left, params.cds_right)

    target_contigs = list_and_file_to_dict_keys(params.contigs, params.contigs_file, as_set=True)
    selected_domains = list_and_file_to_dict_keys(params.domains, params.domains_file, as_set=True)
    domain_expr = params.domain_expr


    ### Run
    #TODO: resolve code duplication with domainator domain_search
    write_seq_outputs(trim_contigs(
        parse_seqfiles(params.input, target_contigs),
        kb_truncation=kb_truncation,
        cds_truncation=cds_truncation,
        domains=selected_domains,
        domain_expr=domain_expr,
        unannotated=params.unannotated,
        evalue=params.evalue,
        databases=set(params.databases) if params.databases is not None else None,
        ), genbank=params.genbank_out, fasta=params.fasta_out)
