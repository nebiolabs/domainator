"""Filter domain anntations

Takes a genbank file which has been annotated through the domainator program and filters domainator annotations,
based on percentage overlap between domains and evalues of domains.

"""

import argparse
import sys
from domainator.utils import list_and_file_to_dict_keys, write_genbank, parse_seqfiles
from domainator import __version__, DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME, RawAndDefaultsFormatter
from typing import List, Tuple, Optional, Set, Iterable, Dict
from domainator.Bio.SeqFeature import SeqFeature

#TODO: parallelize input parsing
#TODO: the logic of this program is a bit confusing, in terms of interplay between different filters.

def domain_overlap(hits:Dict[int, SeqFeature], allowed_fraction_overlap):
    """
    Removes domains that overlap another domain with a higher bitscore at a percentage
    higher than that of the allowed percent overlap
    Args:
        hits: 

        allowed_fraction_overlap: float
            The highest percentage overlap which will be allowed between domains
    Returns:
        a list of tuples of the hits of a coding sequence with the non-overlapping domains
        with the highest bitscores
    """
    no_overlap_hits = []
    hits = sorted(hits.values(), key=lambda x: float(x.qualifiers['score'][0]),reverse=True)  # sort by bitscore
    # Add the first hit as this is the one with the highest bitscore
    no_overlap_hits.append(hits[0])
    for i in range(1, len(hits)):
        # Number of non-overlapping hits the current hit does not overlap with
        for j in range(0, len(no_overlap_hits)):
            # Does this hit overlap the current kept hits
            overlap_fraction = hits[i].location.overlap_size(no_overlap_hits[j].location) / len(hits[i].location)
            if overlap_fraction > allowed_fraction_overlap:
                break
        else:
            no_overlap_hits.append(hits[i])
    out = {}
    for hit in no_overlap_hits:
        out[id(hit)] = hit
    return out


def find_kept_annotations(evalue, max_overlap, rec, databases_keep:Optional[Set[str]], databases_filter:Optional[Set[str]], selected_domains:Optional[Set[str]], feature_types:Optional[Set[str]]=None )-> Dict[str, Dict[int, SeqFeature]]: # 
    """
    Finds all the annotations in a contig of the genbank which meet 
    the new evalue and max overlap requirements
    Args:
        evalue: float
            the max evalue which will be accepted in the filtered genbank
        max_overlap: float
            the max overlap between domains which will be accepted
            in the filtered genbank
        rec: SeqRecord
            the current contig being parsed through
        databases_keep: str or None
            before applying filters, remove all domains not in this set of databases.        
        databases_filter: str or None
            the name of the databases which you want to filter annotations from, by default will filter all annoatations.
        selected_domains: set of strings or None. If not None, then only keep domains with names in this set.
        feature_types: set of strings or None. Whether to filter Domainator, or Domain_search annotations, or both. None is equivalent to {"Domainator"}.
        
    Returns:
        a dictionary mapping the number of the coding sequence to a list of tuples 
        with all of the hits of the coding sequence which meet the new filter requirements
    """
    kept_annotations = {} # {cds_id: {id: feature}}
    overlap_annotations = {} # {cds_id: {id: feature}}
    if feature_types is None:
        feature_types = {DOMAIN_FEATURE_NAME}
    for feature in rec.features:
        if feature.type not in feature_types:  # skip over non-domainator features
            continue
        # set the key for each CDS in the contig
        cds_id = feature.qualifiers['cds_id'][0]
        db = feature.qualifiers['database'][0]
        name = feature.qualifiers['name'][0]
        if databases_keep is not None and db not in databases_keep:
            continue

        annot_evalue = float(feature.qualifiers['evalue'][0])
        if cds_id not in kept_annotations:
            kept_annotations[cds_id] = {}
 

        if databases_filter is not None and db not in databases_filter:
            # we're not filtering this database, so just keep everything regardless of evalue and overlap.
            kept_annotations[cds_id][id(feature)] = feature
            continue
        
        if selected_domains is not None and name not in selected_domains:
            # we're filtering by domain name, and this domain is not in the list, so skip it.
            continue

        if annot_evalue > evalue:
            # skip over domainator features which have a higher evalue than the new evalue
            continue  

        # save domainator features which have a low enough evalue to corresponding cds
        if cds_id not in overlap_annotations:
            overlap_annotations[cds_id] = {}


        overlap_annotations[cds_id][id(feature)] = feature

    
    # check for overlap between domainator features of each CDS
    
    for cds_id in overlap_annotations:
        if len(overlap_annotations[cds_id]) > 0:
                if max_overlap < 1:
                    overlap_annotations[cds_id] = domain_overlap(
                        overlap_annotations[cds_id], max_overlap)
        for feature_id in overlap_annotations[cds_id]:
            kept_annotations[cds_id][feature_id] = overlap_annotations[cds_id][feature_id]
    return kept_annotations

def filter_domains(sequence_iterator, evalue, max_overlap, databases_keep:Optional[Set[str]]=None, databases_filter:Optional[Set[str]]=None, selected_domains=None, annotation_type="all"):
    """
    The main function which iterates through genbank files to find the domainator annotations
    to be removed
    Args:
        sequence_iterator: list of lists of SeqRecords
        evalue: float
            The new upper limit for evalues on domainator annotations
        max_overlap: float
            The new upper limit for max overlap between domainator annotations
        database: 
            database to filter from, if None, then filter all databases together
        selected_domains: set of strings or None. If not None, then only keep domains with names in this set.
    """
    feature_types = set()
    if annotation_type == "domain_search":
        feature_types.add(DOMAIN_SEARCH_BEST_HIT_NAME)
    elif annotation_type == "domainator":
        feature_types.add(DOMAIN_FEATURE_NAME)
    elif annotation_type == "all":
        feature_types.add(DOMAIN_FEATURE_NAME)
        feature_types.add(DOMAIN_SEARCH_BEST_HIT_NAME)
    else:
        raise ValueError(f"Invalid annotation type: {annotation_type}")



    for rec in sequence_iterator:
        keep = find_kept_annotations(
            evalue, max_overlap, rec, databases_keep, databases_filter, selected_domains, feature_types)

        new_features = []
        for feature in rec.features:
            if feature.type in feature_types:
                if feature.qualifiers['cds_id'][0] in keep and id(feature) in keep[feature.qualifiers['cds_id'][0]]:
                    new_features.append(feature)
            else:
                new_features.append(feature)
        
        rec.features = new_features
        yield rec


def main(argv):
    parser = argparse.ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument('-i', '--input', nargs='+', required=False,
                       default=None,
                       help="the genbank filenames to be filtered. If not supplied, reads from stdin.")
    
    parser.add_argument('-o', '--output', default=None, type=str, required=False,
                       help="output genbank filename. If not supplied writes to stdout.")

    parser.add_argument('-e', '--evalue', type=float, default=100000000,
                       help="threshold E value for the domain hit")
    parser.add_argument('-m', '--max_overlap', type=float, default=1.0,
                       help="the maximum fraction of overlap between domains to be included in the annotated genbank.")

    parser.add_argument('--contigs', default=[], nargs='+',
                        help="only annotate contigs with ids in this list. Additive with --contigs_file. default: annotate all contigs.")

    parser.add_argument('--contigs_file', default=None, 
                        help="only annotate contigs with ids listed in this file (one per line). Additive with --contigs. default: annotate all contigs.")
    
    parser.add_argument('--domains', type=str, default=None, nargs='+',
                        help="Remove all domains NOT matching these names") 
    parser.add_argument('--domains_file', type=str, default=None,
                        help="text file containing names of domains to include. (all domains not matching this list will be removed).")

    parser.add_argument('--databases_keep', type=str, default=None, nargs="+",
                       help="Remove annotations from all databases except these. Applied before other filters, such as domains, evalue, and max_overlap. By default will keep all annotations from all databases.")

    parser.add_argument('--databases_filter', type=str, default=None, nargs="+",
                       help="The name of the database which you want to filter annotations from, by default will filter all annoatations according to evalue, domains, and max_overlap. Domains from other databases will be left as-is.")
    
    parser.add_argument('--annotation_type', type=str, default="all", choices={"domainator", "domain_search", "all"},
                       help="Which type of annotations to filter (domainator, domain_search, all). Default: domainator")
    

    params = parser.parse_args(argv)

    target_contigs = list_and_file_to_dict_keys(params.contigs, params.contigs_file)    
    selected_domains = list_and_file_to_dict_keys(
        params.domains, params.domains_file, as_set=True)

    evalue = params.evalue
    max_overlap = params.max_overlap


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

    databases_filter = None
    if params.databases_filter is not None:
        databases_filter = set(params.databases_filter)

    databases_keep = None
    if params.databases_keep is not None:
        databases_keep = set(params.databases_keep)

    # Run
    write_genbank(filter_domains(
                parse_seqfiles(genbanks, target_contigs, filetype_override=filetype_override), 
                evalue=evalue,
                max_overlap=max_overlap,
                databases_keep=databases_keep,
                databases_filter=databases_filter,
                selected_domains=selected_domains,
                annotation_type=params.annotation_type
                ),
        out)

def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])
