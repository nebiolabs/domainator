"""Extract unannotated regions

Takes a genbank file which has been annotated using the domainator program, looks for the specified domains, and extracts regions not containing them.

This tool is currently only supports protein sequences, support for nucleic acids may be added someday.

can write a fasta or a genbank file
"""

from jsonargparse import ActionConfigFile
from domainator.cli.common import add_annotation_type, add_input, add_seq_outputs, check_outputs, tool_parser
from domainator.outputs import write_seq_outputs
from domainator.utils import parse_seqfiles, list_and_file_to_dict_keys, slice_record_from_location, region_id
from domainator import DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME
from domainator.Bio.SeqFeature import FeatureLocation


#TODO: better filtering of output records when input is nucleic acids. Use functions from domainator select_by_cds to extract CDSs first, then extract the domains from the individual CDSs.
def extract_unannotated(records, evalue=None, score=None, keep_direction=False, keep_name=False, search_hits=False, databases=None, largest=False, lb=0):
    """

      input: 
        records: iterator of SeqRecords
        evalue: if not None then only select domains with evalue <= this parameter
        score: if not None then only select domains with score >= this parameter
        keep_direction: has no effect, because only protein sequences are supported. Kept for compatibility with extract_domains.
        keep_name: by default, extracted regions will be renamed to note their original coordinates (see utils.region_id). set this option to keep contigs with their original names. WARNING: if you aren't careful, this can easily result in duplicate names.
        search_hits: if True, then only select domains that are marked as search hits (from domainator domain_search).
        databases: if not None, then only select domains from these databases. default: all databases.
        largest: if True, then only select the largest unannotated region from each protein, providing it passes the length cutoff specified by --min_length
        lb: if largest is True, then only select unannotated regions of at least this length from each protein. default: no cutoff.

      yields the extracted regions as SeqRecord objects
    """
    databases = set(databases) if databases is not None else None
    for rec in records:
        hit_domains = list() #list of SeqFeatures
        if rec.annotations['molecule_type'] != "protein":
            raise ValueError("extract_unannotated currently only supports protein sequences.")
        
        for feature in rec.features:
            if (search_hits and feature.type == DOMAIN_SEARCH_BEST_HIT_NAME) or (feature.type == DOMAIN_FEATURE_NAME):
                if databases is not None and feature.qualifiers["database"][0] not in databases:
                    continue
                keep = True
                if evalue is not None and float(feature.qualifiers["evalue"][0]) > evalue:
                    keep = False
                elif score is not None and float(feature.qualifiers["score"][0]) < score:
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
            unannotated_regions = domains_envelope.inverted(length=len(rec))
        
        if unannotated_regions is None:
            continue
        
        unannotated_regions = unannotated_regions.parts
        # filter by size
        if lb is not None:
            unannotated_regions = [region for region in unannotated_regions if len(region) >= lb]

        if len(unannotated_regions) > 0:
            if largest:
                unannotated_regions.sort(key=lambda x: len(x), reverse=True)
                unannotated_regions = [unannotated_regions[0]]

            contig_sections = list()
            for location in unannotated_regions:
                contig_sections.append(slice_record_from_location(rec, location)) # slice_record_from_location handles maintaining source and taxonomy annotations.
                if not keep_name:
                    contig_sections[-1].id = region_id(rec, location)

            for out_rec in contig_sections:
                yield out_rec 

def build_parser():
    parser = tool_parser(__doc__)

    add_input(parser, 'Genbank filenames.')

    add_seq_outputs(parser)

    parser.add_argument('--contigs', default=[], nargs='+',
                        help="only extract regions from contigs with ids in this list. Additive with --contigs_file. default: all contigs.")

    parser.add_argument('--contigs_file', default=None, 
                        help="only extract regions from contigs with ids listed in this file (one per line). Additive with --contigs. default: all contigs.")

    parser.add_argument('-e', '--evalue', default=None, type=float,
                        help="the evalue cutoff (max) for domain annotations to be considered. default: no cutoff.")
    
    parser.add_argument("--min_score", type=float, default=None, required=False, 
                        help="the score cutoff (min) for domain annotations to be considered. default: no cutoff.")

    parser.add_argument("--databases", default=None, required=False, type=str, nargs="+",
                        help="Only consider domains from these databases. default: all databases.")

    add_annotation_type(parser, help="Which annotations count as annotated.")

    parser.add_argument('--keep_direction', action='store_true', default=False,
                       help="by default extracted regions will be flipped so that the source CDS will be on the forward strand. Set this option to not reverse complement anything. In the rare case of a split CDS with parts on different strands, the strand of the first part will be considered.")

    parser.add_argument('--keep_name', action='store_true', default=False,
                       help="by default, extracted regions will be renamed to note their original coordinates and whether they were reverse complemented. Set this option to keep contigs with their original names. WARNING: if you aren't careful, this can easily result in duplicate names.")
    
    parser.add_argument('--largest', action='store_true', default=False,
                        help="Extract only the largest unannotated region from each protein, providing it passes the length cutoff specified by --min_length")
    
    parser.add_argument('--min_length', type=int, default=None,
                        help=f"Extract unannotated regions of at least this length from each protein. default: no cutoff. ")

    parser.add_argument('--config', action=ActionConfigFile)
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)
   
    ### Figure out what input and output files ####

    genbanks = params.input

    check_outputs(parser, params, ["genbank_out", "fasta_out"])

    target_contigs = list_and_file_to_dict_keys(
        params.contigs, params.contigs_file, as_set=True)

    
    # Run
    extracted_domains_iterator = extract_unannotated(
        parse_seqfiles(genbanks, target_contigs), 
        evalue = params.evalue,
        score = params.min_score,
        keep_direction = params.keep_direction, 
        keep_name = params.keep_name,
        search_hits=params.annotation_type == "domain_search",
        databases=params.databases,
        largest=params.largest,
        lb=params.min_length,
        )

    write_seq_outputs(extracted_domains_iterator, genbank=params.genbank_out, fasta=params.fasta_out)
