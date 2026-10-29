"""Extract domain-annotated regions

Takes a genbank file which has been annotated using the domainator program, looks for the specified domains, and extracts them from the contigs.

This tool is typically used with protein sequences, but can also be used with nucleic acid sequences.
Handling of nucleic acids is less stable and will may change in future updates, so use with caution.

can write a fasta or a genbank file
"""

#TODO: better filtering of output records when input is nucleic acids. Use functions from domainator select_by_cds to extract CDSs first, then extract the domains from the individual CDSs.

from typing import Iterable, Iterator, Optional, Set
from jsonargparse import ActionConfigFile
from jsonargparse.typing import PositiveFloat
from domainator.cli.common import add_annotation_type, add_has_search_hit, add_input, add_seq_outputs, check_outputs, declare_outputs, tool_parser
from domainator.outputs import write_seq_outputs
from domainator.seqio import parse_seqfiles
from domainator.cli.common import list_and_file_to_dict_keys, partial_help
from domainator.locations import FeatureIndex, slice_record_from_location, pad_location
from domainator.model import region_id, get_cds_unique_name, annotation_feature_types, feature_score, evalue_passes
from domainator.partial import cds_fragment_status, fragment_status_allowed, PARTIAL_CHOICES
from domainator import DOMAIN_SEARCH_BEST_HIT_NAME
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from domainator.location_ops import merge_parts
from Bio.SeqRecord import SeqRecord
from domainator.cli.common import add_taxonomy_filters, taxonomy_from_args
from domainator.taxonomy import NCBITaxonomy
from domainator.taxonomy.records import compile_taxonomy_filter

def extract_domains(records: Iterable[SeqRecord], *, evalue: Optional[float] = None, min_score: Optional[float] = None,
                    domains: Optional[Set[str]] = None, extend_range: Optional[int] = None,
                    extend_up: Optional[int] = None, extend_down: Optional[int] = None, combine: bool = False,
                    keep_direction: bool = False, keep_parent_id: bool = False, annotation_type: str = "domainator",
                    databases: Optional[Set[str]] = None, splice: bool = False, partial: str = "include",
                    has_search_hit: bool = False, ncbi_taxonomy: Optional[NCBITaxonomy] = None,
                    include_taxids: Optional[Iterable[int]] = None, exclude_taxids: Optional[Iterable[int]] = None,
                    taxonomy_expr: Optional[str] = None, taxonomy_unknown: str = "exclude") -> Iterator[SeqRecord]:
    """Yield the regions of ``records`` covered by their domain annotations.

    Args:
        records: the contigs or proteins (SeqRecords).
        evalue: extract only domain/search-hit annotations with an e-value at most this,
            before extension, combination, or splicing. None: no cutoff.
        min_score: extract only domain/search-hit annotations with a score at least this,
            before extension, combination, or splicing. None: no cutoff.
        domains: select only domains with these names. None: all.
        extend_range: the default of extend_up and extend_down (None: 0).
        extend_up: also extract this many residues (or bases) upstream of each domain. None:
            extend_range.
        extend_down: also extract this many residues (or bases) downstream of each domain. None:
            extend_range.
        combine: when a contig has several selected domains, extract one region enclosing them
            all. Raises an error for a circular contig, whose enclosing region is ambiguous.
        keep_direction: don't reverse-complement regions. By default, regions are flipped so
            that the domain (with ``combine`` or ``splice``, the one with the lowest start) is on
            the forward strand.
        keep_parent_id: use the parent record ID instead of a coordinate-based region ID.
            Applies to nucleotide and protein inputs and can produce duplicate IDs for
            different extracted sequences.
        annotation_type: which annotations to extract: "domainator", "domain_search", or "all".
        has_search_hit: also extract search hits, ignoring domain names and databases but
            honoring evalue, min_score, and partial. Alone, extracts only search hits.
        databases: select only domains from these databases. None: all.
        ncbi_taxonomy: database required when taxonomy criteria are supplied.
        include_taxids: include these taxids and descendants (OR).
        exclude_taxids: exclude these taxids and descendants, overriding inclusion.
        taxonomy_expr: boolean lineage expression, exclusive with include/exclude_taxids.
        taxonomy_unknown: unresolved taxonomy policy: exclude (default), include, evaluate
            against an empty lineage, or error. Applied to the original domain location
            before extension/combine/splice; protein inputs use record taxonomy.
        splice: when a contig has several selected domains, delete the sequence between them and
            join the domains, in the order they are read on the strand of the domain with the
            lowest start.
        partial: select domains by the fragment status of their CDS (or protein record):
            "include" (all, the default), "exclude" (skip domains on fragments), or "only" (only
            domains on fragments). Contig-level nucleotide hits are never on fragments.
            See partial.get_fragment_status.

    Yields:
        The extracted regions.
    """
    feature_types = annotation_feature_types(annotation_type)
    taxonomy_filter = compile_taxonomy_filter(ncbi_taxonomy, include_taxids, exclude_taxids,
                                              taxonomy_expr, taxonomy_unknown)

    if splice and combine:
        raise ValueError("Cannot use --splice and --combine at the same time.")

    databases = set(databases) if databases is not None else None
    extend_default = extend_range if extend_range is not None else 0
    extend_up = extend_up if extend_up is not None else extend_default
    extend_down = extend_down if extend_down is not None else extend_default
    #TODO: ensure naming consistency with select_by_cds.
    for rec in records:
        feature_index = FeatureIndex(rec.features)
        hit_domains = list() #list of SeqFeatures
        if partial != "include":
            cds_by_id = {get_cds_unique_name(f): f for f in rec.features if f.type == "CDS"}
        for feature in rec.features:
            search_hit = has_search_hit and feature.type == DOMAIN_SEARCH_BEST_HIT_NAME
            selected = (
                feature.type in feature_types
                and (domains is not None or not has_search_hit)
                and (domains is None or feature.qualifiers["name"][0] in domains)
                and (databases is None or feature.qualifiers["database"][0] in databases)
            )
            if selected or search_hit:
                keep = True
                if not evalue_passes(feature, evalue):
                    keep = False
                elif min_score is not None and feature_score(feature) < min_score:
                    keep = False
                elif partial != "include":
                    cds_id = feature.qualifiers.get("cds_id", ["."])[0]
                    if rec.annotations.get("molecule_type") == "protein":
                        status = cds_fragment_status(rec, None)
                    elif cds_id in cds_by_id:
                        status = cds_fragment_status(rec, cds_by_id[cds_id])
                    else: # contig-level nucleotide hit
                        status = None
                    keep = fragment_status_allowed(status, partial)
                if keep and taxonomy_filter is not None:
                    keep = taxonomy_filter.matches(rec, feature)
                if keep:
                    hit_domains.append(feature)
        
        if len(hit_domains) > 0:
            contig_sections = list()
            if combine:
                if rec.annotations.get("topology", None) == "circular": #TODO: this is kind of silly, because we should extract/linearize the CDS first
                    raise ValueError(f"Cannot combine domains on circular contigs because the direction to join on is ambiguous (you could go either way around the circle). Contig {rec.id} is circular.") #TODO: maybe warn and skip?
                start_i = min(list(range(len(hit_domains))), key = lambda x: int(hit_domains[x].location.start)) # find start coordinate of the domain with the lowest start coordinate
                end_i   = max(list(range(len(hit_domains))), key = lambda x: int(hit_domains[x].location.end)  ) # find end coordinate of the domain with the highest end coordinate
                if hit_domains[start_i].location.parts[0].strand == -1: # if the domain is on the reverse strand, then we pad down from start and up from end.
                    start = hit_domains[start_i].location.start - extend_down
                    end = hit_domains[end_i].location.end + extend_up
                else: # if the domain is on the forward strand, then we pad up from start and down from end.
                    start = hit_domains[start_i].location.start - extend_up
                    end = hit_domains[end_i].location.end + extend_down
    
                if keep_direction:
                    new_strand = 1
                else:
                    new_strand = hit_domains[start_i].location.parts[0].strand

                if start < 0:
                    start = 0
                if end > len(rec): 
                    end = len(rec)

                
                location = FeatureLocation(start, end, new_strand)
                contig_sections.append(slice_record_from_location(rec, location, feature_index.near(location))) # slice_record_from_location handles maintaining source and taxonomy annotations.
                if not keep_parent_id:
                    contig_sections[-1].id = region_id(rec, location)

            elif splice:
                # delete all sequence not covered by selected domains and glue parts together.
                # As for --combine, the domain with the lowest start coordinate decides the direction of the output.
                first_domain = min(hit_domains, key=lambda x: int(x.location.start))
                reverse = not keep_direction and first_domain.location.parts[0].strand == -1
                # Domains are glued in the order they are read: from right to left on the reverse strand.
                # The parts of each domain are already in reading order, which keeps domains that span the origin of circular contigs intact.
                if reverse:
                    ordered_domains = sorted(hit_domains, key=lambda x: int(x.location.parts[0].end), reverse=True)
                else:
                    ordered_domains = sorted(hit_domains, key=lambda x: int(x.location.parts[0].start))
                splice_parts = list()
                for domain in ordered_domains:
                    parts = pad_location(rec, domain.location, extend_up, extend_down).parts
                    if keep_direction:
                        if parts[0].strand == -1: # reverse strand parts are read from right to left, so read them from left to right on the forward strand.
                            parts = parts[::-1]
                        splice_parts.extend(FeatureLocation(int(part.start), int(part.end), strand=1) for part in parts)
                    else:
                        splice_parts.extend(FeatureLocation(int(part.start), int(part.end), strand=part.strand) for part in parts)

                if not splice_parts:
                    continue

                merged_parts = merge_parts(splice_parts)

                if len(merged_parts) == 1:
                    splice_location = merged_parts[0]
                else:
                    splice_location = CompoundLocation(merged_parts, operator="join")

                new_rec = slice_record_from_location(rec, splice_location, feature_index.near(splice_location)) # slice_record_from_location handles maintaining source and taxonomy annotations.
                if rec.annotations.get("topology", None) == "circular":
                    new_rec.annotations["topology"] = "circular"

                if not keep_parent_id:
                    new_rec.id = region_id(rec, splice_location)

                contig_sections.append(new_rec)
            
            else: # return single domains
                for domain in hit_domains:
                    strand = domain.location.parts[0].strand # strand of the first part.
                    
                    location = pad_location(rec, domain.location, extend_up, extend_down)

                    contig_sections.append(slice_record_from_location(rec, location, feature_index.near(location))) # slice_record_from_location handles maintaining source and taxonomy annotations.

                    # slice_from_location will flip the strand if the domain is on the reverse strand, so we need to flip it back if we want to keep the direction.
                    flip = rec.annotations['molecule_type'] != "protein" and keep_direction and strand == -1
                    if flip:
                        contig_sections[-1] = contig_sections[-1].reverse_complement(id=rec.id, name=True, description=True, annotations=True) 

                    if not keep_parent_id:
                        contig_sections[-1].id = region_id(rec, location, flip=flip)

            for out_rec in contig_sections:
                yield out_rec 

def build_parser():
    parser = tool_parser(__doc__)

    add_input(parser, 'Genbank filenames.')

    add_seq_outputs(parser)
    add_taxonomy_filters(parser)

    parser.add_argument('--contigs', default=[], nargs='+',
                        help="only extract domains from contigs with ids in this list. Additive with --contigs_file. default: all contigs.")

    parser.add_argument('--contigs_file', default=None, 
                        help="only extract domains from contigs with ids listed in this file (one per line). Additive with --contigs. default: all contigs.")

    parser.add_argument('-e', '--evalue', default=None, type=PositiveFloat,
                        help="Extract only domain/search-hit annotations with e-value <= this, before extension, combination, or splicing. Default: no cutoff.")
    parser.add_argument("--min_score", type=float, default=None, required=False, 
                        help="Extract only domain/search-hit annotations with score >= this, before extension, combination, or splicing. Default: no cutoff.")

    parser.add_argument('--domains', type=str, default=None, nargs='+',
                        help="Include domains with these names")  # TODO: add some way to specify the database of a domain
    parser.add_argument('--domains_file', type=str, default=None,
                        help="text file containing names of domains to include.")
        
    parser.add_argument("--databases", default=None, required=False, type=str, nargs="+",
                        help="Domain databases to use for domain extraction. default: all databases.")

    add_annotation_type(parser, help="Which annotations to extract.")
    add_has_search_hit(parser)
    parser.add_argument('--partial', type=str, default="include", choices=PARTIAL_CHOICES, help=partial_help("Which domains to extract, by the CDS or protein they are on") + " Contig-level nucleotide hits are never on fragments.")

    parser.add_argument("--extend_range", type=int, default=None, help="extract this many additional residues (or bases) on each side of the hit region. Sets the default of --extend_up and --extend_down.")
    parser.add_argument("--extend_up", type=int, default=None, help="extract this many additional residues (or bases) upstream of the hit region. Default: --extend_range, or 0.")
    parser.add_argument("--extend_down", type=int, default=None, help="extract this many additional residues (or bases) downstream of the hit region. Default: --extend_range, or 0.")

    parser.add_argument("--combine", action="store_true", default=False, help="if set, then when multiple target domains occur in the same contig, an envelope containing all of the hit domains will be extracted. "
                        "By default, each hit domain will be extracted separately. This option will raise an exception on circular contigs, because it would be ambiguous. Cannot be used with --splice.")
    
    parser.add_argument("--splice", action="store_true", default=False, help="if set, then when multiple target domains occur in the same contig, the intermediate sequence will be deleted and the extracted domains from the contig glued together, in the order they are read on the strand of the domain with the lowest start coordinate. "
                        "By default, each hit domain will be extracted separately. Cannot be used with --combine.")
    
    parser.add_argument('--keep_direction', action='store_true', default=False,
                       help="by default extracted regions will be flipped so that the domain (or the domain with the lowest start coordinate when --combine or --splice is active) is on the forward strand. Set this option to not reverse complement anything. In the rare case of a split domain with parts on different strands, the strand of the first part will be considered.")

    parser.add_argument('--keep_parent_id', action='store_true', default=False,
                       help="Use the parent nucleotide or protein record ID instead of a coordinate-based region ID. This can produce duplicate IDs for different extracted sequences.")

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
    selected_domains = list_and_file_to_dict_keys(
        params.domains, params.domains_file, as_set=True)
    
    # Run
    extracted_domains_iterator = extract_domains(
        parse_seqfiles(genbanks, target_contigs), 
        evalue=params.evalue, 
        min_score=params.min_score, 
        domains=selected_domains, 
        extend_range=params.extend_range,
        extend_up=params.extend_up,
        extend_down=params.extend_down,
        combine=params.combine, 
        keep_direction=params.keep_direction, 
        keep_parent_id=params.keep_parent_id,
        annotation_type=params.annotation_type,
        has_search_hit=params.has_search_hit,
        databases=params.databases,
        splice=params.splice,
        partial=params.partial,
        ncbi_taxonomy=taxonomy_from_args(params),
        include_taxids=params.include_taxids,
        exclude_taxids=params.exclude_taxids,
        taxonomy_expr=params.taxonomy_expr,
        taxonomy_unknown=params.taxonomy_unknown,
        )

    write_seq_outputs(extracted_domains_iterator, genbank=params.genbank_out, fasta=params.fasta_out)
