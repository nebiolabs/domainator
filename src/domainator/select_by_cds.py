"""Extract contig regions around selected CDSs

Takes a domain-annotated genbank file and extracts regions of contigs based on contig or CDS name or presence of specified domains within individual CDSs.

Specify a range to extract around selected CDSs, if no range is specified, then only the selected CDSs will be returned and not any of their neighbors.

if target_cdss, target_domains or domain_expr is specified then selection logic is:

    contigs & (target_cdss or target_domains or domain_expr or unannotated)

otherwise:
    return every cds in contigs

"""
from jsonargparse import ArgumentParser, ActionConfigFile
import sys
from domainator.utils import parse_seqfiles, write_genbank, list_and_file_to_dict_keys, BooleanEvaluator, DomainatorCDS, pad_location, slice_record_from_location, FeatureLocation, get_non_domainator_features, circular_dist
from typing import Tuple, List, Optional, Set
from domainator import __version__, RawAndDefaultsFormatter, DOMAIN_SEARCH_BEST_HIT_NAME

#TODO: account for different domain databases
#TODO: account for domain order, like by making boolean searches binned. Try to emulate the behavior of the InterPro domain selector tool:
#       https://www.ebi.ac.uk/interpro/search/ida/
#TODO: add evalue or bitscore threshold for domains?
#TODO: domain description wildcard searches (using sqlite "LIKE" function) for example to pull all domains with "methyl" in the description
#TODO: when streaming have option to list non-selected contigs

#TODO: add support for selecting by qualifiers

def clear_best_hit_features(input_features: List, keep: str) -> List:
    """returns a copy of the input feature list, but with all features of type domainator.DOMAIN_SEARCH_BEST_HIT_NAME removed except the one with cds_id == keep

    Args:
        input_features (List): list[seqfeature]
        keep (str): cds_id 

    Returns:
        List: list of seqfeatures with all but one domainator.DOMAIN_SEARCH_BEST_HIT_NAME removed.
    """
    out_features = list()

    for feature in input_features:
        if feature.type != DOMAIN_SEARCH_BEST_HIT_NAME:
            out_features.append(feature)
        elif feature.qualifiers['cds_id'][0] == keep:
            out_features.append(feature)
    return out_features    

def pad_records(records, pad_char="N"):
    if len(records) == 0:
        return 
    offsets = [x.get_dist_to_start() for x in records]
    lengths = [len(x) for x in records]
    max_dist = max(offsets)
    shift_amount = [max_dist - x for x in offsets]
    ends = [shift_amount[i] + lengths[i] for i in range(len(records))]
    max_end = max(ends)
    right_pad_amount = [max_end - x for x in ends]
    for i in range(len(records)):
        records[i] = (pad_char * shift_amount[i]) + records[i] + (pad_char * right_pad_amount[i])

def pad_circular_cds(center:int, size:int, left_pad:int, right_pad:int) -> Tuple[int,int]:
    """
    center: the index of the focus CDS in the list of CDSs
    size: the number of CDSs in the contig
    left_pad: the number of CDSs to pad to the left of the focus CDS
    right_pad: the number of CDSs to pad to the right of the focus CDS

    returns: a tuple of (left_cds, right_cds) where left_cds is the index of the left-most CDS in the padded region, and right_cds is the index of the right-most CDS in the padded region.
    """
    if center < 0 or center >= size:
        raise(ValueError("center must be >= 0 and < size"))
    if left_pad < 0:
        raise(ValueError("Padding must be positive"))
    if right_pad < 0:
        raise(ValueError("Padding must be positive"))

    if size == 1: # cannot pad when there is only 1 CDS
        left_pad = 0
        right_pad = 0

    end_to_origin = size-(center + 1) # 
    min_coord = center-left_pad
    max_coord = center + right_pad
    
    if min_coord >= 0 and max_coord < size:
        # the bounds do not extend beyond the origin, on either side, so they can't overlap
        return min_coord, max_coord
    else:
        # If regions extending beyond the origin run into the center, truncate them there.
        # We want a maximum coverage of one for any CDS.
        
        # adjust min_coord and max_coord so that they don't wrap around the center.
        # we've already handled the case of size == 1, we should not be setting min_coord or max_coord to the center here.
        if min_coord < 0 and ((-1 * min_coord) > end_to_origin):
            min_coord = -1 * (end_to_origin)
        if max_coord >= size and (max_coord - size >= center):
            max_coord = size + (center -1)
        
        # If the adjusted regions overlap with each other, then divide the overlap equally between min and max.
        # shift the perspective so that min_coord is at the origin
        min_coord_shifted = 0 # min_coord_shifted = min_coord - min_coord
        max_coord_shifted = max_coord - min_coord
        
        # calculuate overlap size:
        overlap_size = (max_coord_shifted + 1) - size
        if overlap_size > 0:
            # divide overlap equally between min and max
            min_coord_shifted += overlap_size // 2
            max_coord_shifted -= overlap_size // 2
            # if overlap_size is odd, then we need to adjust one of the coordinates by 1
            if overlap_size % 2 == 1:
                    max_coord_shifted -= 1 
        # shift the perspective back to the original coordinate system
        max_coord = max_coord_shifted + min_coord
        min_coord = min_coord_shifted + min_coord
        if min_coord < 0:
            min_coord = size + min_coord

        return (min_coord, max_coord % size) # TODO: is % size necessary?

def get_cds_neighborhood(contig, cds_list, cds_idx, cds_range: Tuple[int, int]=None, kb_range:Tuple[float, float]=None, whole_contig=False, normalize_direction=True, _from_domain_search=False):
    """
        contig: a SeqRecord object

        cds_list: a list of DomainatorCDS object derived from the contig, MUST BE SORTED by stranded start coordinate

        cds_idx: the index of the focus cds within cds_list

        cds_range: extract a contig region enclosing this many CDSs upstream and downstream of the selected CDSs. A 2-tuple of (upstream-range, downstream-range). If None, then don't use a CDS range.

        kb_range: extract a contig region enclosing this many kb upstream and downstream of the selected CDSs. Partially enclosed CDSs will not be annotated in the output. A 2-tuple of (upstream-range, downstream-range). If None, then don't use a kb range.

        whole_contig: extract the whole contigs containing the selected CDSs (if a single contig contains multiple selected CDSs, only one copy of the contig will be returned)

        normalize_direction: if True then if a target cds occurs on the reverse strand, reverse-complement the extracted region before returning

        _from_domain_search: if True, then remove DOMAIN_SEARCH_BEST_HIT_NAME records from CDSs other than the focus CDS.
    """
   
    contig_id = contig.id
    focus_cds = cds_list[cds_idx]

    strand = focus_cds.feature.location.parts[0].strand

    if sum((cds_range is not None, kb_range is not None, whole_contig)) > 1:
        raise RuntimeError(
            "No more than one of: cds_range, kb_range, or whole_contig may be specified.")

    if "topology" not in contig.annotations or contig.annotations["topology"] == "linear":
        upstream_pad = 0
        downstream_pad = 0
        if whole_contig:
            upstream_pad = len(contig)
            downstream_pad = len(contig)
        elif kb_range is not None:
            upstream_pad = kb_range[0] * 1000
            downstream_pad = kb_range[1] * 1000
        elif cds_range is not None:
            if strand == -1: # reverse strand
                left_cds = max(0, cds_idx - cds_range[1])
                right_cds = min(len(cds_list) - 1, cds_idx + cds_range[0])
            else: #forward strand
                left_cds = max(0, cds_idx - cds_range[0])
                right_cds = min(len(cds_list) - 1, cds_idx + cds_range[1])
            left_coord = int(cds_list[left_cds].feature.location.start)
            right_coord = int(cds_list[right_cds].feature.location.end)
            left_pad = max(int(focus_cds.feature.location.start) - left_coord, 0) # min in case of overlapping CDSs
            right_pad = max(right_coord - int(focus_cds.feature.location.end), 0) # max in case of overlapping CDSs
            if strand == -1:
                upstream_pad = right_pad
                downstream_pad = left_pad
            else:
                upstream_pad = left_pad
                downstream_pad = right_pad
        else: # just keep this CDS
            pass

        envelope = FeatureLocation(focus_cds.feature.location.start, focus_cds.feature.location.end, strand)
        slice_location, (dist_to_start, dist_to_end) = pad_location(contig, envelope, int(upstream_pad), int(downstream_pad), return_extension_sizes=True)
    else: #Topology is circular, so we need to find the selection edges using modular arithmetic
        if whole_contig:
            upstream_pad = len(contig)
            downstream_pad = len(contig)
            slice_location, (dist_to_start, dist_to_end) = pad_location(contig, focus_cds.feature.location.circular_envelope(len(contig), strand), upstream_pad, downstream_pad, return_extension_sizes=True)
        elif kb_range is not None or cds_range is None:
            upstream_pad = 0
            downstream_pad = 0
            if kb_range is not None:
                upstream_pad = int(kb_range[0] * 1000)
                downstream_pad = int(kb_range[1] * 1000)
            slice_location, (dist_to_start, dist_to_end) = pad_location(contig, focus_cds.feature.location.circular_envelope(len(contig), strand), upstream_pad, downstream_pad, return_extension_sizes=True)
        elif cds_range is not None:
            if strand == -1:
                left_cds, right_cds = pad_circular_cds(cds_idx, len(cds_list), cds_range[1], cds_range[0])
            else:
                left_cds, right_cds = pad_circular_cds(cds_idx, len(cds_list), cds_range[0], cds_range[1])

            #TODO: account for overlapping CDSs
            left_envelope = cds_list[left_cds].feature.location.circular_envelope(len(contig), 1)
            focus_envelope = focus_cds.feature.location.circular_envelope(len(contig), 1)
            right_envelope = cds_list[right_cds].feature.location.circular_envelope(len(contig), 1)
            neighborhood = left_envelope + focus_envelope + right_envelope

            slice_location = neighborhood.circular_envelope(len(contig), strand)
            dist_location = slice_location
            if dist_location.parts[0].strand != 1:
                dist_location = dist_location.copy()
                for part in dist_location.parts:
                    part.strand = 1       

            dist_to_start = circular_dist(dist_location.stranded_start, focus_envelope.stranded_start, len(contig))
            dist_to_end = circular_dist(focus_envelope.stranded_end, dist_location.stranded_end, len(contig))
            if strand == -1: # reverse strand, so swap dist to start and end
                tmp = dist_to_start
                dist_to_start = dist_to_end
                dist_to_end = tmp
            
            if len(slice_location) == len(contig): # Account for the edge case where neighborhood spans the entire contig, which would make the neighborhood improperly centered.
                slice_location, (dist_to_start, dist_to_end) = pad_location(contig, focus_cds.feature.location.circular_envelope(len(contig), strand), len(contig), len(contig), return_extension_sizes=True)

    features = get_non_domainator_features(contig)
    
    for cds in cds_list:
        if slice_location.contains(cds.feature.location):
            # features.append(cds.feature)
            features.extend(cds.domain_features)
            if cds.domain_search_feature is not None:
                features.append(cds.domain_search_feature)
    try:
        record = slice_record_from_location(contig, slice_location, features)
        start = int(slice_location.stranded_start)
        end = int(slice_location.stranded_end)
    except Exception as e:
        print(f"contig: {contig.id}", file=sys.stderr)
        print(f"slice_location: {slice_location}", file=sys.stderr)
        print(f"focus_cds: {focus_cds}", file=sys.stderr)
        
        raise e


    # we need to add coordinates to the name any time we are making a subset of the record or changing the origin.
    if (start == 0) and end == len(contig) and strand != -1: #forward strand whole contig starting at the origin, so don't rename
        pass    
    elif (end == 0) and start == len(contig) and strand == -1: #reverse strand whole contig
        pass
    else:
        record.id = record.id + f"_{start+1}:{end}"

    # slice_from_location will flip the strand if the domain is on the reverse strand, so we need to flip it back if we want to keep the direction.
    if not normalize_direction and strand == -1: 
        record = record.reverse_complement(id=record.id, name=True, description=True, annotations=True)
        tmp = dist_to_start
        dist_to_start = dist_to_end
        dist_to_end = tmp
    elif strand == -1: # if strand of the location is -1, then slice_record_from_location will have flipped the strand, so we need to mark that in the name.
        record.id += "rc"
        
    if _from_domain_search:
        record.features = clear_best_hit_features(record.features, keep=focus_cds.num)
    record.set_dist_to_start(dist_to_start)
    return record, slice_location

def select_by_cds(contigs, target_cdss=None, target_domains=None, domain_expr=None, cds_range: Tuple[int, int]=None, kb_range:Tuple[float, float]=None, whole_contig=False, normalize_direction=True, evalue=float("inf"), invert=False, search_hits=False, max_region_overlap=1.0, strand=None, _from_domain_search=False, _domain_search_negatives:Optional[Set[str]]=None, databases:Optional[Set[str]]=None, unannotated=False):
    """

        Args:
            contigs: an interator over SeqRecord objects

            target_contigs: only extract regions from contigs with these names

            target_cdss: extract CDSs with these names

            target_domains: Select all CDSs containing domains with these names

            domain_expr: Use boolean expressions for domain names to select CDSs
            
            cds_range: extract a contig region enclosing this many CDSs upstream and downstream of the selected CDSs

            kb_range: extract a contig region enclosing this many kb upstream and downstream of the selected CDSs. Partially enclosed CDSs will not be annotated in the output.

            whole_contig: extract the whole contigs of containing the selected CDSs (if a single contig contains multiple selected CDSs, only one copy of the contig will be returned)

            normalize_direction: if True then if a target cds occurs on the reverse strand, reverse-complement the extracted region before returning

            evalue: ignore domains with evalue higher than this.

            invert: if True then invert the target_cdss, target_domains, and domain_expr selection criteria (contigs, and ranges are still selected normally, only CDS selection is inverted)

            search_hits: Select CDSs that are marked as search hits

            max_region_overlap: the maximum fractional of overlap between any two output regions. If >= 1, then no overlap filtering will be done.

            strand: only extract regions around CDSs on the specified strand. If None, then extract regions around CDSs on both strands. forward strand is "f", reverse strand is "r".

            _from_domain_search: if True then read target_cdss from rec.get_hit_scores

            _domain_search_negatives (Optional[Set[str]]): if _from_domain_search is True, then this is a set of domain names to negatively select on.

            databases: a set of database names to consider when filtering by domain. If None, then all databases will be considered.

        Yields:
            SeqRecords of the selected regions
    """ 
    databases = set(databases) if databases is not None else None

    if target_cdss is not None:
        target_cdss = set(target_cdss)

    if target_domains is not None:
        target_domains = set(target_domains)

    if domain_expr is not None:
        expression_evaluator = BooleanEvaluator(domain_expr)

    for rec in contigs:
        # First make a list of all the focus CDSs, by index
        cdss = DomainatorCDS.list_from_contig(rec, evalue)
        cdss.sort(key=lambda x: x.feature.location.stranded_start) #TODO: why is this stranded_start, not start?
        if max_region_overlap < 1.0:
            selected_locations = list()
        
        if _from_domain_search:
            hit_scores = rec.get_hit_scores()
            hit_names = rec.get_hit_names()

        for focus_index, cds in enumerate(cdss):
            keep = False
            #If no selectors then keep every cds
            if target_cdss is None and target_domains is None and domain_expr is None and not search_hits and (not _from_domain_search) and (not unannotated):
                keep = True
            elif _from_domain_search:
                if cds.index in hit_scores:
                    if _domain_search_negatives is None or hit_names[cds.index] not in _domain_search_negatives:
                        # if _domain_search_negatives is not None, then we are doing a negative selection, so skip any cds that has a negative domain
                        keep = True
                        best_score = hit_scores[cds.index]
            else:
                if target_cdss is not None: # if target_cdss then update based on cds name
                    if cds.name in target_cdss:
                        keep = True
                if target_domains is not None: # if target_domains, then update based on domain names
                    if len(cds.get_domain_names(databases).intersection(target_domains)):
                        keep = True
                if domain_expr is not None: # if target_domains, then update based on domain names
                    if expression_evaluator.check_expression(cds.get_domain_names(databases)):
                        keep = True
                if unannotated:
                    if len(cds.get_domain_names(databases)) == 0:
                        keep = True
                if search_hits and cds.domain_search_feature is not None:
                    keep = True
            
            if strand is not None:
                if strand == "f" and cds.feature.location.parts[0].strand == -1: # if strand is forward, but the cds is on the reverse strand, then don't keep it
                    keep = False
                elif strand == "r" and cds.feature.location.parts[0].strand != -1:
                    keep = False

            if invert:
                keep = not keep

            if keep:
                record, record_location = get_cds_neighborhood(rec, cdss, focus_index, cds_range=cds_range, kb_range=kb_range, whole_contig=whole_contig, normalize_direction=normalize_direction, _from_domain_search=_from_domain_search)
                if _from_domain_search:
                    record.set_hit_best_score(best_score)
                
                if max_region_overlap < 1.0:
                    loc_len = len(record_location)
                    for other_location in selected_locations:
                        if (record_location.overlap_size(other_location) / loc_len) > max_region_overlap:
                            break
                    else:
                        selected_locations.append(record_location)
                        yield record
                else:
                    yield record

def validate_range(kb_range, kb_up, kb_down, cds_range, cds_up, cds_down):
    """

    """

    if kb_range is not None and kb_range < 0:
        raise ValueError("kb_range must be >= 0")
    if kb_up is not None and kb_up < 0:
        raise ValueError("kb_range must be >= 0")
    if kb_down is not None and kb_down < 0:
        raise ValueError("kb_range must be >= 0")

    if cds_range is not None and cds_range < 0:
        raise ValueError("cds_range must be >= 0")
    if cds_up is not None and cds_up < 0:
        raise ValueError("cds_range must be >= 0")
    if cds_down is not None and cds_down < 0:
        raise ValueError("cds_range must be >= 0")

    cds_range_seen = False
    cds_range_out = None
    if cds_up is not None or cds_down is not None:
        cds_range_out = (0,0)

    if cds_up is not None:
        cds_range_out = (cds_up, cds_range_out[1])
        cds_range_seen = True
    if cds_down is not None:
        cds_range_out = (cds_range_out[0], cds_down)
        cds_range_seen = True

    if cds_range is not None:
        if cds_range_seen:
            raise ValueError("Cannot specify cds_range with cds_down or cds_up")
        cds_range_seen = True
        cds_range_out = (cds_range, cds_range)


    if cds_range_seen and ((kb_range is not None) or (kb_up is not None) or (kb_down is not None)):
        raise ValueError("Cannot specify kb_range and cds_range")



    kb_range_seen = False
    kb_range_out = None
    if kb_up is not None or kb_down is not None:
        kb_range_out = (0,0)

    if kb_up is not None:
        kb_range_out = (kb_up, kb_range_out[1])
        kb_range_seen = True
    if kb_down is not None:
        kb_range_out = (kb_range_out[0], kb_down)
        kb_range_seen = True

    if kb_range is not None:
        if kb_range_seen:
            raise ValueError("Cannot specify kb_range with kb_down or kb_up")
        kb_range_seen = True
        kb_range_out = (kb_range, kb_range)
   
    return (kb_range_out, cds_range_out)


def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument("-i", "--input", default=None, nargs="+", type=str, required=False,
                        help="names of input genbank files. If not supplied, reads from stdin.")

    parser.add_argument('-o', '--output', default=None, required=False, 
                        help="genbank output file name. If not supplied, writes to stdout.")


    parser.add_argument('--kb_range', type=float, default=None,
                        help="How many kb upstream and downstream of the selected CDS to extract. Starting from the ends of the CDS.")
    parser.add_argument('--kb_up', type=float, default=None,
                        help="How many kb upstream of the selected CDS to extract. Starting from the end of the CDS.")
    parser.add_argument('--kb_down', type=float, default=None,
                        help="How many kb downstream of the selected CDS to extract. Starting from the end of the CDS.")

    parser.add_argument('--cds_range', type=int, default=None,
                        help="How many CDSs upstream and downstream of the selected CDS to extract. Example, 1 would mean to return contigs including the selected CDSs and the immediate upstream and downstream CDSs.")
    parser.add_argument('--cds_up', type=int, default=None,
                        help="How many CDSs upstream of the selected CDS to extract. Example, 1 would mean to return contigs including the selected CDSs and the immediate upstream CDS.")
    parser.add_argument('--cds_down', type=int, default=None,
                        help="How many CDSs downstream of the selected CDS to extract. Example, 1 would mean to return contigs including the selected CDSs and the immediate downstream CDS.")

    parser.add_argument('--whole_contig', default=False, action="store_true",
                        help="If set, then return the entire matching contigs.")

           
    parser.add_argument('--contigs', type=str, default=None, nargs='+',
                        help="Only extract regions from contigs with these names.")
    parser.add_argument('--contigs_file', type=str, default=None,
                        help="text file containing names of contigs to extract from.")
    parser.add_argument('--cds', type=str, default=None, nargs='+',
                        help="Include CDSs with these names (additive with --domains)")
    parser.add_argument('--cds_file', type=str, default=None,
                        help="text file containing names of CDSs to include.")
    parser.add_argument('--domains', type=str, default=None, nargs='+',
                        help="Include domains with these names") #TODO: add some way to specify the database of a domain
    parser.add_argument('--domains_file', type=str, default=None,
                        help="text file containing names of domains to include.")
    parser.add_argument('--domain_expr', type=str, default=None,
                        help="a boolean expression using operators & (AND), | (OR), and ~(NOT), to specify a desired combination of domains")
    
    parser.add_argument("--databases", default=None, required=False, type=str, nargs="+",
                        help="Consider only domains from these databases when filtering by domain. default: all databases.")

    parser.add_argument('--strand', type=str, default=None, choices=["f", "r"], help="Only extract regions around CDSs on the specified strand.")

    parser.add_argument('--unannotated', action='store_true', default=False, help="Select contigs with no domain annotations.")

    parser.add_argument('--search_hits', action='store_true', default=False, help="Select CDSs that are marked as search hits (from domain_search.py).")
    #TODO: 
    # parser.add_argument('--domain_count_expr', type=str, default=None,
    #                     help="a comparator and a number, to specify a desired number of domains. Example: '> 2' would mean to return contigs with more than 2 domains.")

    parser.add_argument('--pad', action='store_true', default=False,
                        help="If set, then ends of sequences will be padded so that the focus CDS aligns.")

    parser.add_argument('--keep_direction', action='store_true', default=False,
                        help="by default extracted regions will be flipped so that the focus cds is on the forward strand. Setting this option will keep the focus cds on whatever strand it started on.")
    parser.add_argument('--skip_deduplicate', action='store_true', default=False,
                        help="by default if the same region is extracted for multiple hits, then only one will be kept. Set this option to retain redundancies (which might be useful if you are counting hits), but will result in duplicated contig names, which could be a problem for downstream domainator applications.")
    parser.add_argument('--max_region_overlap', type=float, default=1.0,
                       help="the maximum fractional of overlap between any two output regions. If >= 1, then no overlap filtering will be done. Regions are output in a greedy fashion based on CDS start site. New regions are output if less than this fraction of them overlaps with any previously output region.")


    parser.add_argument('-e', '--evalue', default=100000000, type=float,
                        help="the evalue cutoff for domain annotations to be put onto the new genbank file")

    parser.add_argument('--invert', action='store_true',
                        help="Invert the CDS selection criteria. i.e. return CDSs that don't match the CDS selection criteria. (only applies to CDS selection, not neighborhoods or contigs).")

    parser.add_argument('--config', action=ActionConfigFile)
    

    # parser.add_argument(, '--architecture', type=str, default=None,
    #                    help="space separated list of domains")
    # parser.add_argument(, '--architecture_file', type=str, default=None,
    #                    help="space separated lists of domains, one architecture per line")
    # parser.add_argument(, '--subarchitecture', type=str, default=None,
    #                    help="space separated list of domains")
    # parser.add_argument(, '--subarchitecture_file', type=str, default=None,
    #                    help="space separated lists of domains, one architecture per line")

    # parser.add_argument('-d', '--database', type=str, default=None,
    #                    help="the name of the database which you want to filter annotations from.") #TODO: see todo above.

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
    kb_range, cds_range = validate_range(params.kb_range, params.kb_up, params.kb_down, params.cds_range, params.cds_up, params.cds_down)

    whole_contig = params.whole_contig
    normalize_direction = not params.keep_direction
    deduplicate = not params.skip_deduplicate
    pad = params.pad

    target_contigs = list_and_file_to_dict_keys(params.contigs, params.contigs_file, as_set=True)
    selected_cdss = list_and_file_to_dict_keys(params.cds, params.cds_file, as_set=True)
    selected_domains = list_and_file_to_dict_keys(params.domains, params.domains_file, as_set=True)
    domain_expr = params.domain_expr


    ### Run
    #TODO: resolve code duplication with domain_search.py
    seen = set()
    records = list()
    for record in select_by_cds(
        parse_seqfiles(genbanks, target_contigs, filetype_override=filetype_override), 
        selected_cdss, 
        selected_domains, 
        domain_expr, 
        cds_range, 
        kb_range, 
        whole_contig, 
        normalize_direction,
        params.evalue,
        search_hits = params.search_hits,
        invert = params.invert,
        max_region_overlap = params.max_region_overlap,
        strand = params.strand,
        databases = params.databases,
        unannotated = params.unannotated
        ):
        #skip_deduplicate, pad
        if not pad:
            if not deduplicate or (record.id not in seen):
                write_genbank((record,), output_handle)
        else:
            if not deduplicate or (record.id not in seen):
                records.append(record)
        if deduplicate:
            seen.add(record.id)

    if pad:
        pad_records(records)
        write_genbank(records, output_handle)

    if params.output is not None:
        output_handle.close()

def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])
  