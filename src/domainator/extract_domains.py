"""Extract domain-annotated regions

Takes a genbank file which has been annotated using the domainator program, looks for the specified domains, and extracts them from the contigs.

This tool is typically used with protein sequences, but can also be used with nucleic acid sequences.
Handling of nucleic acids is less stable and will may change in future updates, so use with caution.

can write a fasta or a genbank file
"""

#TODO: better filtering of output records when input is nucleic acids. Use functions from select_by_cds.py to extract CDSs first, then extract the domains from the individual CDSs.

import sys
from jsonargparse import ArgumentParser, ActionConfigFile
from domainator.Bio import SeqIO
from domainator.utils import parse_seqfiles, list_and_file_to_dict_keys, write_genbank, slice_record_from_location, slice_record, get_sources, pad_location
from domainator import __version__, DOMAIN_FEATURE_NAME, RawAndDefaultsFormatter, DOMAIN_SEARCH_BEST_HIT_NAME
from domainator.Bio.SeqFeature import FeatureLocation, CompoundLocation

def extract_domains(records, evalue=None, score=None, domains=None, pad_up=0, pad_down=0, combine=False, keep_direction=False, keep_name=False, search_hits=False, databases=None, splice=False):
    """

      input: 
        records: iterator of SeqRecords
        evalue: if not None then only select domains with evalue <= this parameter
        score: if not None then only select domains with score >= this parameter
        domains: collection if not None then only select domains with names in this list.
        pad_up: extract this many additional residues (or bases) upstream of the hit region.
        pad_down: extract this many additional residues (or bases) downstream of the hit region.
        combine: if True, then when multiple target domains occur in the same contig, an envelope containing all of the hit domains will be extracted. Will raise an error if the contig is circular, because the envelope cannot be determined.
        keep_direction: by default extracted regions will be flipped so that the domain (or the domain with the lowest start coordinate when --combine is active) is on the forawd strand. Set this option to not reverse complement anything.
        keep_name: by default, extracted regions will be renamed to note their original coordinates and whether they were reverse complemented. set this option to keep contigs with their original names. WARNING: if you aren't careful, this can easily result in duplicate names.
        search_hits: if True, then only select domains that are marked as search hits (from domain_search.py).
        databases: if not None, then only select domains from these databases. default: all databases.
        splice: if True, then when multiple target domains occur in the same contig, the intermediate sequence will be deleted and the extracted domains from the contig glued together.s
    
      yields the extracted regions as SeqRecord objects
    """

    if splice and combine:
        raise ValueError("Cannot use --splice and --combine at the same time.")



    databases = set(databases) if databases is not None else None
    #TODO: ensure naming consistency with select_by_cds.
    for rec in records:
        hit_domains = list() #list of SeqFeatures
        for feature in rec.features:
            if (search_hits and feature.type == DOMAIN_SEARCH_BEST_HIT_NAME) or (feature.type == DOMAIN_FEATURE_NAME):
                if databases is not None and feature.qualifiers["database"][0] not in databases:
                    continue
                keep = True
                if feature.type == DOMAIN_FEATURE_NAME and domains is not None and feature.qualifiers["name"][0] not in domains:
                    keep = False
                elif (search_hits and feature.type == DOMAIN_FEATURE_NAME) and domains is None:
                    keep = False
                elif evalue is not None and float(feature.qualifiers["evalue"][0]) > evalue:
                    keep = False
                elif score is not None and float(feature.qualifiers["score"][0]) < score:
                    keep = False
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
                    start = hit_domains[start_i].location.start - pad_down
                    end = hit_domains[end_i].location.end + pad_up
                else: # if the domain is on the forward strand, then we pad up from start and down from end.
                    start = hit_domains[start_i].location.start - pad_up
                    end = hit_domains[end_i].location.end + pad_down
    
                if keep_direction:
                    new_strand = 1
                else:
                    new_strand = hit_domains[start_i].location.parts[0].strand

                if start < 0:
                    start = 0
                if end > len(rec): 
                    end = len(rec)

                
                location = FeatureLocation(start, end, new_strand)
                contig_sections.append(slice_record_from_location(rec, location)) # slice_record_from_location handles maintaining source and taxonomy annotations.

                new_id = rec.id
                if len(contig_sections[-1]) != len(rec): # if the slice is not the same length as the original, then we need to add the coordinates to the name.
                    new_id = contig_sections[-1].id + f"_{location.stranded_start_human_readable}:{location.stranded_end_human_readable}"

                if rec.annotations['molecule_type'] != "protein" and new_strand == -1: # if new_strand is -1, then slice_record_from_location will have flipped the strand, so we need to mark that in the name.
                    # contig_sections[-1] = contig_sections[-1].reverse_complement(id=rec.id, name=True, description=True, annotations=True) 
                    new_id += "rc"
                if not keep_name:
                    contig_sections[-1].id = new_id

            elif splice:
                # delete all sequence not covered by selected domains and glue parts together.
                splice_parts = list()
                for domain in hit_domains:
                    padded_location = pad_location(rec, domain.location, pad_up, pad_down)
                    for part in padded_location.parts:
                        strand = part.strand
                        if keep_direction:
                            strand = 1 if strand in (-1, None) else strand
                        splice_parts.append(
                            FeatureLocation(int(part.start), int(part.end), strand=strand)
                        )

                if not splice_parts:
                    continue

                merged_parts = CompoundLocation.merge_parts(splice_parts)

                if len(merged_parts) == 1:
                    splice_location = merged_parts[0]
                else:
                    splice_location = CompoundLocation(merged_parts, operator="join")

                new_rec = slice_record_from_location(rec, splice_location, truncate_features=True) # slice_record_from_location handles maintaining source and taxonomy annotations.
                if rec.annotations.get("topology", None) == "circular":
                    new_rec.annotations["topology"] = "circular"

                if not keep_name:
                    new_id = new_rec.id
                    if len(new_rec) != len(rec):
                        if isinstance(splice_location, CompoundLocation):
                            parts = splice_location.parts
                            strand = splice_location.strand
                        else:
                            parts = [splice_location]
                            strand = splice_location.strand

                        min_start = min(int(part.start) for part in parts)
                        max_end = max(int(part.end) for part in parts)

                        if strand == -1:
                            start_label = max_end
                            end_label = min_start + 1
                        else:
                            start_label = min_start + 1
                            end_label = max_end

                        new_id = f"{new_id}_{start_label}:{end_label}"
                        if (
                            rec.annotations.get("molecule_type") != "protein"
                            and strand == -1
                            and not keep_direction
                        ):
                            new_id += "rc"
                    new_rec.id = new_id

                contig_sections.append(new_rec)
            
            else: # return single domains
                for domain in hit_domains:
                    strand = domain.location.parts[0].strand # strand of the first part.
                    
                    location = pad_location(rec, domain.location, pad_up, pad_down)

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
                        help="only extract domains from contigs with ids in this list. Additive with --contigs_file. default: all contigs.")

    parser.add_argument('--contigs_file', default=None, 
                        help="only extract domains from contigs with ids listed in this file (one per line). Additive with --contigs. default: all contigs.")


    parser.add_argument('-e', '--evalue', default=None, type=float,
                        help="the evalue cutoff (max) for domain annotations to be extracted")
    parser.add_argument("--score", type=float, default=None, required=False, 
                        help="the score cutoff (min) for domain annotations to be extracted")

    parser.add_argument('--domains', type=str, default=None, nargs='+',
                        help="Include domains with these names")  # TODO: add some way to specify the database of a domain
    parser.add_argument('--domains_file', type=str, default=None,
                        help="text file containing names of domains to include.")
        
    parser.add_argument("--databases", default=None, required=False, type=str, nargs="+",
                        help="Domain databases to use for domain extraction. default: all databases.")

    parser.add_argument('--search_hits', action='store_true', default=False, help="Select domains that are marked as search hits (from domain_search.py).")

    #parser.add_argument("--padding", type=int, default=0, help="extract this many additional residues (or bases) to each side of the hit region.")
    parser.add_argument("--pad_up", type=int, default=0, help="extract this many additional residues (or bases) upstream of the hit region.")
    
    parser.add_argument("--pad_down", type=int, default=0, help="extract this many additional residues (or bases) downstream of the hit region.")

    parser.add_argument("--combine", action="store_true", default=False, help="if set, then when multiple target domains occur in the same contig, an envelope containing all of the hit domains will be extracted. "
                        "By default, each hit domain will be extracted separately. This option will raise an exception on circular contigs, because it would be ambiguous. Cannot be used with --splice.")
    
    parser.add_argument("--splice", action="store_true", default=False, help="if set, then when multiple target domains occur in the same contig, the intermediate sequence will be deleted and the extracted domains from the contig glued together. "
                        "By default, each hit domain will be extracted separately. Cannot be used with --combine.")
    
    parser.add_argument('--keep_direction', action='store_true', default=False,
                       help="by default extracted regions will be flipped so that the domain (or the domain with the lowest start coordinate when --combine is active) is on the forward strand. Set this option to not reverse complement anything. In the rare case of a split domain with parts on different strands, the strand of the first part will be considered.")

    parser.add_argument('--keep_name', action='store_true', default=False,
                       help="by default, extracted regions will be renamed to note their original coordinates and whether they were reverse complemented. Set this option to keep contigs with their original names. WARNING: if you aren't careful, this can easily result in duplicate names.")

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
    selected_domains = list_and_file_to_dict_keys(
        params.domains, params.domains_file, as_set=True)
    
    # Run
    extracted_domains_iterator = extract_domains(
        parse_seqfiles(genbanks, target_contigs, filetype_override=filetype_override), 
        params.evalue, 
        params.score, 
        selected_domains, 
        params.pad_up,
        params.pad_down,
        params.combine, 
        params.keep_direction, 
        params.keep_name,
        search_hits=params.search_hits,
        databases=params.databases,
        splice=params.splice,
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
