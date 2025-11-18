"""Extract unannotated regions

Takes a genbank file which has been annotated using the domainator program, looks for the specified domains, and extracts regions not containing them.

This tool is currently only supports protein sequences, support for nucleic acids may be added someday.

can write a fasta or a genbank file
"""

import sys
from jsonargparse import ArgumentParser, ActionConfigFile
from domainator.Bio import SeqIO
from domainator.utils import parse_seqfiles, list_and_file_to_dict_keys, write_genbank, slice_record_from_location, slice_record, get_sources, pad_location, DomainatorCDS
from domainator import __version__, DOMAIN_FEATURE_NAME, RawAndDefaultsFormatter, DOMAIN_SEARCH_BEST_HIT_NAME
from domainator.Bio.SeqFeature import FeatureLocation


#TODO: better filtering of output records when input is nucleic acids. Use functions from select_by_cds.py to extract CDSs first, then extract the domains from the individual CDSs.
def extract_unannotated(records, evalue=None, score=None, keep_direction=False, keep_name=False, search_hits=False, databases=None, largest=False, lb=0):
    """

      input: 
        records: iterator of SeqRecords
        evalue: if not None then only select domains with evalue <= this parameter
        score: if not None then only select domains with score >= this parameter
        keep_direction: by default extracted regions will be flipped so that the domain (or the domain with the lowest start coordinate when --combine is active) is on the forawd strand. Set this option to not reverse complement anything.
        keep_name: by default, extracted regions will be renamed to note their original coordinates and whether they were reverse complemented. set this option to keep contigs with their original names. WARNING: if you aren't careful, this can easily result in duplicate names.
        search_hits: if True, then only select domains that are marked as search hits (from domain_search.py).
        databases: if not None, then only select domains from these databases. default: all databases.
        largest: if True, then only select the largest unannotated region from each protein, providing it passes the length cutoff specified by '--lb'
        lb: if largest is True, then only select unannotated regions of at least this length from each protein. default: no cutoff.

      yields the extracted regions as SeqRecord objects
    """
    databases = set(databases) if databases is not None else None
    #TODO: ensure naming consistency with select_by_cds.
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
                strand = location.strand # strand of the first part.
                
                contig_sections.append(slice_record_from_location(rec, location)) # slice_record_from_location handles maintaining source and taxonomy annotations.

                new_id = rec.id
                if len(contig_sections[-1]) != len(rec):
                    new_id = contig_sections[-1].id + f"_{location.stranded_start_human_readable}:{location.stranded_end_human_readable}"
                # slice_from_location will flip the strand if the domain is on the reverse strand, so we need to flip it back if we want to keep the direction.
                if rec.annotations['molecule_type'] != "protein" and keep_direction and strand == -1: 
                    contig_sections[-1] = contig_sections[-1].reverse_complement(id=rec.id, name=True, description=True, annotations=True) 
                elif rec.annotations['molecule_type'] != "protein" and strand == -1: # if strand of the location is -1, then slice_record_from_location will have flipped the strand, so we need to mark that in the name.
                    new_id += "rc"
                    
                if not keep_name:
                    contig_sections[-1].id = new_id

            for out_rec in contig_sections:
                yield out_rec 

def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument('-i', '--input', nargs='+', required=False,
                        default=None,
                        help="Genbank filenames. If not supplied, reads from stdin.")

    parser.add_argument('-o', '--output', default=None, required=False,
                        help="the name of the output file. If not supplied, writes to stdout.")

    parser.add_argument('--contigs', default=[], nargs='+',
                        help="only extract regions from contigs with ids in this list. Additive with --contigs_file. default: all contigs.")

    parser.add_argument('--contigs_file', default=None, 
                        help="only extract regions from contigs with ids listed in this file (one per line). Additive with --contigs. default: all contigs.")

    parser.add_argument('-e', '--evalue', default=None, type=float,
                        help="the evalue cutoff (max) for domain annotations to be considered. default: no cutoff.")
    
    parser.add_argument("--score", type=float, default=None, required=False, 
                        help="the score cutoff (min) for domain annotations to be considered. default: no cutoff.")

    parser.add_argument("--databases", default=None, required=False, type=str, nargs="+",
                        help="Only consider domains from these databases. default: all databases.")

    parser.add_argument('--search_hits', action='store_true', default=False, help="Consider domains that are marked as search hits (from domain_search.py).")

    parser.add_argument('--keep_direction', action='store_true', default=False,
                       help="by default extracted regions will be flipped so that the source CDS will be on the forward strand. Set this option to not reverse complement anything. In the rare case of a split CDS with parts on different strands, the strand of the first part will be considered.")

    parser.add_argument('--keep_name', action='store_true', default=False,
                       help="by default, extracted regions will be renamed to note their original coordinates and whether they were reverse complemented. Set this option to keep contigs with their original names. WARNING: if you aren't careful, this can easily result in duplicate names.")
    
    parser.add_argument('--largest', action='store_true', default=False,
                        help="Extract only the largest unannotated region from each protein, providing it passes the length cutoff specified by '--lb'")
    
    parser.add_argument('--lb', type=int, default=None,
                        help=f"Extract unannotated regions of at least this length from each protein. default: no cutoff. ")

    parser.add_argument('--fasta_out', action='store_true', default=False,
                        help="makes output a fasta file when activated")

    parser.add_argument('--config', action=ActionConfigFile)
        
    params = parser.parse_args(argv)
   
    ### Figure out what input and output files ####

    filetype_override = None
    if params.input is None:
        genbanks = [sys.stdin]
        filetype_override = "genbank"
    else:
        genbanks = params.input

    if params.output is None:
        out = sys.stdout
    else:
        out = open(params.output, "w")


    target_contigs = list_and_file_to_dict_keys(
        params.contigs, params.contigs_file, as_set=True)

    
    # Run
    extracted_domains_iterator = extract_unannotated(
        parse_seqfiles(genbanks, target_contigs, filetype_override=filetype_override), 
        evalue = params.evalue,
        score = params.score,
        keep_direction = params.keep_direction, 
        keep_name = params.keep_name,
        search_hits=params.search_hits,
        databases=params.databases,
        largest=params.largest,
        lb=params.lb,
        )

    if params.fasta_out:
        SeqIO.write(extracted_domains_iterator, out, "fasta")
    else:
        write_genbank(extracted_domains_iterator, out)

    if params.output is not None:
        out.close()

def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])
