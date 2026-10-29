"""Extract unannotated regions

Takes a genbank file which has been annotated using the domainator program, looks for the specified domains, and extracts regions not containing them.

This tool supports only protein sequences.

can write a fasta or a genbank file
"""

from typing import Iterable, Iterator, Optional, Set
from jsonargparse import ActionConfigFile
from jsonargparse.typing import PositiveFloat
from domainator.cli.common import add_annotation_type, add_input, add_seq_outputs, check_outputs, declare_outputs, tool_parser
from domainator.outputs import write_seq_outputs
from domainator.seqio import parse_seqfiles
from domainator.cli.common import list_and_file_to_dict_keys
from domainator.locations import slice_record_from_location
from domainator.model import region_id, annotation_feature_types, feature_score, evalue_passes
from Bio.SeqFeature import FeatureLocation
from domainator.location_ops import inverted
from Bio.SeqRecord import SeqRecord
from domainator.cli.common import add_taxonomy_filters, taxonomy_from_args
from domainator.taxonomy import NCBITaxonomy
from domainator.taxonomy.records import compile_taxonomy_filter


#TODO: better filtering of output records when input is nucleic acids. Use functions from domainator select_by_cds to extract CDSs first, then extract the domains from the individual CDSs.
def extract_unannotated(records: Iterable[SeqRecord], *, evalue: Optional[float] = None, min_annotation_score: Optional[float] = None,
                        keep_direction: bool = False, keep_parent_id: bool = False, annotation_type: str = "domainator",
                        databases: Optional[Set[str]] = None, largest: bool = False, min_region_length: Optional[int] = None,
                        ncbi_taxonomy: Optional[NCBITaxonomy] = None,
                        include_taxids: Optional[Iterable[int]] = None, exclude_taxids: Optional[Iterable[int]] = None,
                        taxonomy_expr: Optional[str] = None, taxonomy_unknown: str = "exclude") -> Iterator[SeqRecord]:
    """Yield the regions of the proteins in ``records`` that have no domain annotation.

    Args:
        records: protein SeqRecords (only proteins are supported).
        evalue: count only domain/search-hit annotations with an e-value at most this
            when finding unannotated gaps. None: no cutoff.
        min_annotation_score: count only domain/search-hit annotations with a score at least this
            when finding unannotated gaps. Raising this can enlarge the gaps.
            None: no cutoff.
        keep_direction: has no effect, since only proteins are supported; kept for
            compatibility with extract_domains.
        keep_parent_id: use the parent protein ID instead of a coordinate-based region ID.
            This can produce duplicate IDs for different extracted sequences.
        annotation_type: which annotations count as annotated: "domainator", "domain_search",
            or "all".
        databases: count only domains from these databases as annotations. None: all.
        largest: yield only the largest unannotated region of each protein that is at least
            ``min_region_length`` long.
        min_region_length: yield only unannotated gaps at least this many amino acids long,
            before choosing the largest gap. Not a cutoff on the input protein length.
        ncbi_taxonomy: database required when taxonomy criteria are supplied.
        include_taxids: include these taxids and descendants (OR).
        exclude_taxids: exclude these taxids and descendants, overriding inclusion.
        taxonomy_expr: boolean lineage expression, exclusive with include/exclude_taxids.
        taxonomy_unknown: unresolved taxonomy policy: exclude (default), include,
            evaluate against an empty lineage, or error. Filters original protein records
            before extracting gaps.

    Yields:
        The unannotated regions.
    """
    feature_types = annotation_feature_types(annotation_type)
    taxonomy_filter = compile_taxonomy_filter(ncbi_taxonomy, include_taxids, exclude_taxids,
                                              taxonomy_expr, taxonomy_unknown)
    databases = set(databases) if databases is not None else None
    for rec in records:
        hit_domains = list() #list of SeqFeatures
        if rec.annotations['molecule_type'] != "protein":
            raise ValueError("extract_unannotated currently only supports protein sequences.")
        if taxonomy_filter is not None and not taxonomy_filter.matches(rec):
            continue
        
        for feature in rec.features:
            if feature.type in feature_types:
                if databases is not None and feature.qualifiers["database"][0] not in databases:
                    continue
                keep = True
                if not evalue_passes(feature, evalue):
                    keep = False
                elif min_annotation_score is not None and feature_score(feature) < min_annotation_score:
                    keep = False
                if keep:
                    hit_domains.append(feature)
        
        
        domains_envelope = None
        for feature in hit_domains:
            if domains_envelope is None:
                domains_envelope = feature.location
            else:
                domains_envelope = domains_envelope + feature.location
        
        if domains_envelope is None:
            unannotated_regions = FeatureLocation(0, len(rec))
        else:
            unannotated_regions = inverted(domains_envelope, length=len(rec))
        
        if unannotated_regions is None:
            continue
        
        unannotated_regions = unannotated_regions.parts
        # filter by size
        if min_region_length is not None:
            unannotated_regions = [region for region in unannotated_regions if len(region) >= min_region_length]

        if len(unannotated_regions) > 0:
            if largest:
                unannotated_regions.sort(key=lambda x: len(x), reverse=True)
                unannotated_regions = [unannotated_regions[0]]

            contig_sections = list()
            for location in unannotated_regions:
                contig_sections.append(slice_record_from_location(rec, location)) # slice_record_from_location handles maintaining source and taxonomy annotations.
                if not keep_parent_id:
                    contig_sections[-1].id = region_id(rec, location)

            for out_rec in contig_sections:
                yield out_rec 

def build_parser():
    parser = tool_parser(__doc__)

    add_input(parser, 'Genbank filenames.')

    add_seq_outputs(parser)
    add_taxonomy_filters(parser)

    parser.add_argument('--contigs', default=[], nargs='+',
                        help="only extract regions from contigs with ids in this list. Additive with --contigs_file. default: all contigs.")

    parser.add_argument('--contigs_file', default=None, 
                        help="only extract regions from contigs with ids listed in this file (one per line). Additive with --contigs. default: all contigs.")

    parser.add_argument('-e', '--evalue', default=None, type=PositiveFloat,
                        help="Count only domain/search-hit annotations with e-value <= this when finding unannotated gaps. Default: no cutoff.")
    
    parser.add_argument("--min_annotation_score", type=float, default=None, required=False,
                        help="Count only domain/search-hit annotations with score >= this when finding unannotated gaps. Raising this can enlarge the gaps. Default: no cutoff.")

    parser.add_argument("--databases", default=None, required=False, type=str, nargs="+",
                        help="Only consider domains from these databases. default: all databases.")

    add_annotation_type(parser, help="Which annotations count as annotated.")

    parser.add_argument('--keep_direction', action='store_true', default=False,
                       help="Has no effect: this tool supports only protein sequences and never reverse-complements them. Retained for compatibility with extract_domains.")

    parser.add_argument('--keep_parent_id', action='store_true', default=False,
                       help="Use the parent protein ID instead of a coordinate-based region ID. This can produce duplicate IDs for different extracted sequences.")
    
    parser.add_argument('--largest', action='store_true', default=False,
                        help="Extract only the largest unannotated region from each protein, providing it passes the length cutoff specified by --min_region_length")
    
    parser.add_argument('--min_region_length', type=int, default=None,
                        help="Extract only unannotated gaps of at least this many amino acids, before choosing --largest. Not a cutoff on the input protein length. Default: no cutoff.")

    parser.add_argument('--config', action=ActionConfigFile)
    declare_outputs(parser, ["genbank_out", "fasta_out"])
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)
   
    ### Figure out what input and output files ####

    genbanks = params.input

    check_outputs(parser, params)

    target_contigs = list_and_file_to_dict_keys(
        params.contigs, params.contigs_file, as_set=True)

    
    # Run
    extracted_domains_iterator = extract_unannotated(
        parse_seqfiles(genbanks, target_contigs), 
        evalue = params.evalue,
        min_annotation_score = params.min_annotation_score,
        keep_direction = params.keep_direction, 
        keep_parent_id = params.keep_parent_id,
        annotation_type=params.annotation_type,
        databases=params.databases,
        ncbi_taxonomy=taxonomy_from_args(params),
        include_taxids=params.include_taxids,
        exclude_taxids=params.exclude_taxids,
        taxonomy_expr=params.taxonomy_expr,
        taxonomy_unknown=params.taxonomy_unknown,
        largest=params.largest,
        min_region_length=params.min_region_length,
        )

    write_seq_outputs(extracted_domains_iterator, genbank=params.genbank_out, fasta=params.fasta_out)
