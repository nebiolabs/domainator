"""Extract contig regions around selected CDSs

Takes a domain-annotated genbank file and extracts regions of contigs based on contig or CDS name or presence of specified domains within individual CDSs.

Specify a range to extract around selected CDSs, if no range is specified, then only the selected CDSs will be returned and not any of their neighbors.

CDS names, domain names, domain expressions, unannotated, and search-hit presence
are OR selectors. With no selectors, all eligible CDSs are selected.
Contig IDs, strand, fragment status, and taxonomy are independent AND restrictions.
Invert negates only the combined selectors, never these restrictions. With no
selectors, invert warns and is ignored. Restrictions apply to the focus CDS,
not every feature in an expanded neighborhood.

"""
from jsonargparse import ActionConfigFile
from jsonargparse.typing import PositiveFloat
from domainator.cli.common import add_annotation_type, add_dedup, add_has_search_hit, add_input, add_seq_outputs, check_outputs, declare_outputs, tool_parser
from domainator.dedup import RecordDeduplicator
from domainator.outputs import SeqOutputs
from domainator.seqio import parse_seqfiles
from domainator.cli.common import list_and_file_to_dict_keys, partial_help
from domainator.expr import BooleanEvaluator
from domainator.model import DomainatorCDS, cds_name_precedence, get_non_domainator_features, region_id, check_annotation_type, get_dist_to_end, get_hit_best_score, get_dist_to_start, get_focus_strand, get_hit_names, get_hit_scores, set_dist_to_end, set_dist_to_start, set_focus_strand, set_hit_best_score
from domainator import DOMAIN_SEARCH_BEST_HIT_NAME
from domainator.locations import FeatureIndex, pad_location, slice_record_from_location, circular_dist, WHERE_CHOICES
from Bio.SeqFeature import FeatureLocation
from domainator.partial import fragment_status_allowed, cds_fragment_status, PARTIAL_CHOICES
from domainator.location_ops import circular_envelope, copy_location, overlap_size, stranded_end, stranded_start
from typing import Tuple, Optional, Set, Iterable, Iterator
from Bio.SeqRecord import SeqRecord
from domainator.cli.common import add_taxonomy_filters, taxonomy_from_args
from domainator.taxonomy import NCBITaxonomy
from domainator.taxonomy.records import compile_taxonomy_filter
import warnings
#TODO: account for domain order, like by making boolean searches binned. Try to emulate the behavior of the InterPro domain selector tool:
#       https://www.ebi.ac.uk/interpro/search/ida/
#TODO: add evalue or bitscore threshold for domains?
#TODO: domain description wildcard searches (using sqlite "LIKE" function) for example to pull all domains with "methyl" in the description
#TODO: when streaming have option to list non-selected contigs

#TODO: add support for selecting by qualifiers

def pad_anchor(record, where:str="start") -> int:
    """
        The coordinate in record to align on when padding.

        Args:
            record: a record from get_cds_neighborhood, with dist to start, dist to end, and focus strand set
            where: "start", "end", or "center", relative to the focus CDS in its own orientation

        Returns:
            a 0-based coordinate in record. For "end" on the forward strand, it is the position just after the last base of the focus CDS.
    """
    if where not in WHERE_CHOICES:
        raise ValueError(f"where must be one of {', '.join(WHERE_CHOICES)}, not {where}")
    left = get_dist_to_start(record)
    right = len(record) - get_dist_to_end(record)
    length = right - left
    if get_focus_strand(record) == -1:
        if where == "start":
            return right
        elif where == "end":
            return left
        return right - length // 2 # center
    if where == "start":
        return left
    elif where == "end":
        return right
    return left + length // 2 # center

def pad_records(records, pad_char="N", where:str="start"):
    """
        Pads the ends of records in place so that their focus CDSs align.
        Padded records are linear, because the padding is not part of the original molecule.

        Args:
            records: a list of records from get_cds_neighborhood
            pad_char: the character to pad with
            where: "start", "end", or "center". Where to align, relative to the focus CDS in its own orientation.
    """
    if where not in WHERE_CHOICES:
        raise ValueError(f"where must be one of {', '.join(WHERE_CHOICES)}, not {where}")
    if len(records) == 0:
        return 
    offsets = [pad_anchor(x, where) for x in records]
    lengths = [len(x) for x in records]
    max_dist = max(offsets)
    shift_amount = [max_dist - x for x in offsets]
    ends = [shift_amount[i] + lengths[i] for i in range(len(records))]
    max_end = max(ends)
    right_pad_amount = [max_end - x for x in ends]
    for i in range(len(records)):
        if shift_amount[i] == 0 and right_pad_amount[i] == 0:
            continue
        records[i] = (pad_char * shift_amount[i]) + records[i] + (pad_char * right_pad_amount[i])
        if records[i].annotations.get("topology") == "circular":
            records[i].annotations["topology"] = "linear"

def deduplicate_and_pad(regions, *, dedup: str = "none", pad=False, anchor="start", best_first=False):
    """Yield ``regions``, filtered by the ``dedup`` policy, and padded
    with N so that their focus CDSs line up (with ``pad``, which reads every region first).

    ``best_first`` sorts padded domain_search hits best first, as their unpadded order would be.
    """
    deduplicator = RecordDeduplicator(dedup)
    kept = list()
    for region in regions:
        if not deduplicator.keep(region):
            continue
        if pad:
            kept.append(region)
        else:
            yield region
    if pad:
        if best_first:
            kept.sort(key=get_hit_best_score, reverse=True)
        pad_records(kept, where=anchor)
        yield from kept

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
        overflow = (max_coord_shifted + 1) - size
        if overflow > 0:
            # divide overlap equally between min and max
            min_coord_shifted += overflow // 2
            max_coord_shifted -= overflow // 2
            # if overflow is odd, then we need to adjust one of the coordinates by 1
            if overflow % 2 == 1:
                    max_coord_shifted -= 1 
        # shift the perspective back to the original coordinate system
        max_coord = max_coord_shifted + min_coord
        min_coord = min_coord_shifted + min_coord
        if min_coord < 0:
            min_coord = size + min_coord

        return (min_coord, max_coord % size) # TODO: is % size necessary?

def circular_focus_distances(slice_location, focus_location, contig_length:int, strand:int) -> Tuple[int, int]:
    """
        Distances from the start of a slice of a circular contig to the start of a focus feature, and from the end of the focus feature to the end of the slice,
        in the direction of strand.

        slice_location: the slice, as returned by pad_location
        focus_location: the location of the focus feature, which must be inside the slice
        contig_length: the length of the contig
        strand: the strand of the slice and focus feature
    """
    focus_envelope = circular_envelope(focus_location, contig_length, strand)
    if strand == -1:
        dist_to_start = circular_dist(int(stranded_start(focus_envelope)), int(stranded_start(slice_location)), contig_length)
    else:
        dist_to_start = circular_dist(int(stranded_start(slice_location)), int(stranded_start(focus_envelope)), contig_length)
    dist_to_end = len(slice_location) - dist_to_start - len(focus_envelope)
    return dist_to_start, dist_to_end

def neighborhood_feature_index(contig, cds_list) -> FeatureIndex:
    """The features get_cds_neighborhood carries into a neighborhood: the contig's features other than
    domain annotations, then the domain annotations and search hit of each CDS in cds_list."""
    features = get_non_domainator_features(contig)
    for cds in cds_list:
        features.extend(cds.domain_features)
        if cds.domain_search_feature is not None:
            features.append(cds.domain_search_feature)
    return FeatureIndex(features)


def get_cds_neighborhood(contig, cds_list, cds_idx,cds_range: Tuple[int, int]=None, kb_range:Tuple[float, float]=None, whole_contig=False, normalize_direction=True, _from_domain_search=False, _feature_index=None):
    """
        contig: a SeqRecord object

        cds_list: a list of DomainatorCDS object derived from the contig, MUST BE SORTED by stranded start coordinate

        cds_idx: the index of the focus cds within cds_list

        cds_range: extract a contig region enclosing this many CDSs upstream and downstream of the selected CDSs. A 2-tuple of (upstream-range, downstream-range). If None, then don't use a CDS range.

        kb_range: extract a contig region enclosing this many kb upstream and downstream of the selected CDSs. Partially enclosed features are truncated, with '<' or '>' marking the cut ends; partially enclosed domain annotations are dropped. A 2-tuple of (upstream-range, downstream-range). If None, then don't use a kb range.

        whole_contig: extract the whole contigs containing the selected CDSs (if a single contig contains multiple selected CDSs, only one copy of the contig will be returned)

        normalize_direction: if True then if a target cds occurs on the reverse strand, reverse-complement the extracted region before returning

        _from_domain_search: if True, then remove DOMAIN_SEARCH_BEST_HIT_NAME records from CDSs other than the focus CDS.

        _feature_index: neighborhood_feature_index(contig, cds_list), to reuse when getting several neighborhoods of one contig.
    """
   
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
            if _from_domain_search and focus_cds.domain_search_feature is not None:
                # center the contig on the search hit, rather than on the CDS containing it. Distances are still to the focus CDS, so that padding aligns CDSs.
                slice_location = pad_location(contig, circular_envelope(focus_cds.domain_search_feature.location, len(contig), strand), upstream_pad, downstream_pad)
                dist_to_start, dist_to_end = circular_focus_distances(slice_location, focus_cds.feature.location, len(contig), strand)
            else:
                slice_location, (dist_to_start, dist_to_end) = pad_location(contig, circular_envelope(focus_cds.feature.location, len(contig), strand), upstream_pad, downstream_pad, return_extension_sizes=True)
        elif kb_range is not None or cds_range is None:
            upstream_pad = 0
            downstream_pad = 0
            if kb_range is not None:
                upstream_pad = int(kb_range[0] * 1000)
                downstream_pad = int(kb_range[1] * 1000)
            slice_location, (dist_to_start, dist_to_end) = pad_location(contig, circular_envelope(focus_cds.feature.location, len(contig), strand), upstream_pad, downstream_pad, return_extension_sizes=True)
        elif cds_range is not None:
            if strand == -1:
                left_cds, right_cds = pad_circular_cds(cds_idx, len(cds_list), cds_range[1], cds_range[0])
            else:
                left_cds, right_cds = pad_circular_cds(cds_idx, len(cds_list), cds_range[0], cds_range[1])

            #TODO: account for overlapping CDSs
            left_envelope = circular_envelope(cds_list[left_cds].feature.location, len(contig), 1)
            focus_envelope = circular_envelope(focus_cds.feature.location, len(contig), 1)
            right_envelope = circular_envelope(cds_list[right_cds].feature.location, len(contig), 1)
            neighborhood = left_envelope + focus_envelope + right_envelope

            slice_location = circular_envelope(neighborhood, len(contig), strand)
            dist_location = slice_location
            if dist_location.parts[0].strand != 1:
                dist_location = copy_location(dist_location)
                for part in dist_location.parts:
                    part.strand = 1       

            dist_to_start = circular_dist(stranded_start(dist_location), stranded_start(focus_envelope), len(contig))
            dist_to_end = circular_dist(stranded_end(focus_envelope), stranded_end(dist_location), len(contig))
            if strand == -1: # reverse strand, so swap dist to start and end
                tmp = dist_to_start
                dist_to_start = dist_to_end
                dist_to_end = tmp
            
            if len(slice_location) == len(contig): # Account for the edge case where neighborhood spans the entire contig, which would make the neighborhood improperly centered.
                slice_location, (dist_to_start, dist_to_end) = pad_location(contig, circular_envelope(focus_cds.feature.location, len(contig), strand), len(contig), len(contig), return_extension_sizes=True)

    if _feature_index is None:
        _feature_index = neighborhood_feature_index(contig, cds_list)
    features = _feature_index.near(slice_location)
    if _from_domain_search: # only the focus CDS keeps its search hit
        features = [f for f in features if f.type != DOMAIN_SEARCH_BEST_HIT_NAME or f is focus_cds.domain_search_feature]
    # Domain annotations cut by the edge of the region are dropped by slice_record_from_location.
    try:
        record = slice_record_from_location(contig, slice_location, features)
        if contig.annotations.get("topology") == "circular" and len(slice_location) == len(contig): # the whole circular contig, rotated
            record.annotations["topology"] = "circular"
    except Exception as e:
        e.add_note(f"while slicing contig {contig.id} at {slice_location} around CDS {focus_cds}")
        raise


    # slice_from_location will flip the strand if the domain is on the reverse strand, so we need to flip it back if we want to keep the direction.
    flip = not normalize_direction and strand == -1
    if flip:
        record = record.reverse_complement(id=record.id, name=True, description=True, annotations=True)
        tmp = dist_to_start
        dist_to_start = dist_to_end
        dist_to_end = tmp
    record.id = region_id(contig, slice_location, flip=flip)

    set_dist_to_start(record, dist_to_start)
    set_dist_to_end(record, dist_to_end)
    set_focus_strand(record, -1 if flip else 1)
    return record, slice_location

def select_by_cds(records: Iterable[SeqRecord], *, cds_ids: Optional[Set[str]] = None, domains: Optional[Set[str]] = None,
                  domain_expr: Optional[str] = None, cds_range: Optional[Tuple[int, int]] = None,
                  kb_range: Optional[Tuple[float, float]] = None, whole_contig: bool = False, normalize_direction: bool = True,
                  evalue: Optional[float] = None, invert: bool = False, annotation_type: str = "domainator",
                  max_region_overlap: float = 1.0, strand: Optional[str] = None, _from_domain_search: bool = False,
                  _domain_search_negatives: Optional[Set[str]] = None, databases: Optional[Set[str]] = None,
                  unannotated: bool = False, include_nucleic_acids: bool = False, partial: str = "include",
                  has_search_hit: bool = False, name_field: Optional[str] = None,
                  ncbi_taxonomy: Optional[NCBITaxonomy] = None,
                  include_taxids: Optional[Iterable[int]] = None, exclude_taxids: Optional[Iterable[int]] = None,
                  taxonomy_expr: Optional[str] = None, taxonomy_unknown: str = "exclude",
                  ) -> Iterator[SeqRecord]:
    """Yield the regions of ``records`` around the CDSs (or proteins) that are selected.

    A CDS is selected if it matches any of ``cds_ids``, ``domains``, ``domain_expr``,
    ``unannotated`` and ``has_search_hit`` (or with none of them, every CDS), and then the ``strand`` and ``partial``
    filters.

    Args:
        records: the contigs or proteins (SeqRecords).
        cds_ids: select the CDSs with these names.
        domains: select the CDSs with a domain of one of these names.
        domain_expr: select the CDSs whose domain names satisfy this boolean expression.
        cds_range: (upstream, downstream) numbers of CDSs to extract around each selected CDS.
        kb_range: (upstream, downstream) kb to extract around each selected CDS. Partially
            enclosed features are truncated, with '<' or '>' marking the cut ends; partially
            enclosed domain annotations are dropped.
        whole_contig: extract the whole contig of each selected CDS (once, if it has several).
        normalize_direction: reverse-complement regions whose selected CDS is on the reverse
            strand.
        evalue: ignore domains with an e-value above this (None: no cutoff).
        invert: invert the combined ``cds_ids``, ``domains``, ``domain_expr``, ``unannotated``
            and ``has_search_hit`` selection, not strand, partial, or taxonomy.
            With no selectors, warn and ignore inversion.
        annotation_type: which annotations supply domain names: "domainator", "domain_search",
            or "all". Does not itself select CDSs.
        has_search_hit: also select CDSs marked as domain_search hits passing evalue, ignoring
            domain names, expressions, and databases. Participates in invert and enables
            contig-level nucleotide search hits.
        max_region_overlap: the largest fraction of overlap allowed between two output regions
            (1 or more: no filter).
        strand: "f" or "r" to select only CDSs on that strand. None: both.
        name_field: prefer this CDS qualifier for name matching, falling back to
            gene_id, locus_tag, gene, then protein_id. Region IDs remain coordinate-based.
        _from_domain_search: internal (domain_search): select the CDSs recorded as hits on the
            record.
        _domain_search_negatives: internal (domain_search): with _from_domain_search, domain
            names to select against.
        databases: consider only domains from these databases. None: all.
        unannotated: select the CDSs with no domain annotations.
        include_nucleic_acids: treat nucleotide (contig-level) annotations as CDSs, so they can
            be selected too.
        partial: select CDSs by fragment status: "include" (all, the default), "exclude" (skip
            fragments), or "only" (only fragments). A CDS is a fragment if its location has '<'
            or '>' (see partial.get_fragment_status). Not affected by invert.
        ncbi_taxonomy: database required when taxonomy criteria are supplied.
        include_taxids: include these taxids and descendants (OR).
        exclude_taxids: exclude these taxids and descendants, overriding inclusion.
        taxonomy_expr: boolean lineage expression, exclusive with include/exclude_taxids.
        taxonomy_unknown: unresolved taxonomy policy: exclude (default), include, evaluate
            against an empty lineage, or error. Applies to the focus CDS/hit, not the
            expanded neighborhood; not affected by invert.

    Yields:
        SeqRecords of the selected regions.
    """ 
    check_annotation_type(annotation_type)
    taxonomy_filter = compile_taxonomy_filter(ncbi_taxonomy, include_taxids, exclude_taxids,
                                              taxonomy_expr, taxonomy_unknown)
    databases = set(databases) if databases is not None else None

    if cds_ids is not None:
        cds_ids = set(cds_ids)

    if domains is not None:
        domains = set(domains)

    if domain_expr is not None:
        expression_evaluator = BooleanEvaluator(domain_expr)
    name_precedence = cds_name_precedence(name_field)
    no_selectors = (cds_ids is None and domains is None and domain_expr is None
                    and not has_search_hit and not _from_domain_search and not unannotated)
    if no_selectors and invert:
        warnings.warn("Invert supplied with no selection criteria. Invert will be ignored.")
        invert = False

    for rec in records:
        # First make a list of all the focus CDSs, by index
        cdss = DomainatorCDS.list_from_contig(rec, evalue, name_precedence=name_precedence,
                                            include_nucleic_acid_annotations=True) # _from_domain_search
        cdss.sort(key=lambda x: stranded_start(x.feature.location)) #TODO: why is this stranded_start, not start?
        if max_region_overlap < 1.0:
            selected_locations = list()
        
        if _from_domain_search:
            hit_scores = get_hit_scores(rec)
            hit_names = get_hit_names(rec)

        feature_index = None
        for focus_index, cds in enumerate(cdss):
            keep = False
            
            if cds.is_nucleic_acid and not _from_domain_search and not include_nucleic_acids:
                if not has_search_hit or cds.domain_search_feature is None:
                    continue
            #If no selectors then keep every cds
            if no_selectors:
                keep = True
            elif _from_domain_search:
                if cds.index in hit_scores:
                    if _domain_search_negatives is None or hit_names[cds.index] not in _domain_search_negatives:
                        # if _domain_search_negatives is not None, then we are doing a negative selection, so skip any cds that has a negative domain
                        keep = True
                        best_score = hit_scores[cds.index]
            else:
                domain_names = cds.get_domain_names(databases, annotation_type=annotation_type)
                if cds_ids is not None: # if cds_ids then update based on cds name
                    if cds.name in cds_ids:
                        keep = True
                if domains is not None: # if target_domains, then update based on domain names
                    if len(domain_names.intersection(domains)):
                        keep = True
                if domain_expr is not None: # if target_domains, then update based on domain names
                    if expression_evaluator.check_expression(domain_names):
                        keep = True
                if unannotated:
                    if len(domain_names) == 0:
                        keep = True
                if has_search_hit and cds.domain_search_feature is not None:
                    keep = True
            
            if invert:
                keep = not keep

            if strand is not None:
                if strand == "f" and cds.feature.location.parts[0].strand == -1: # if strand is forward, but the cds is on the reverse strand, then don't keep it
                    keep = False
                elif strand == "r" and cds.feature.location.parts[0].strand != -1:
                    keep = False

            if keep and partial != "include" and not fragment_status_allowed(cds_fragment_status(rec, cds.feature), partial):
                keep = False
            if keep and taxonomy_filter is not None:
                keep = taxonomy_filter.matches(rec, cds.feature)

            if keep:
                if feature_index is None:
                    feature_index = neighborhood_feature_index(rec, cdss)
                record, record_location = get_cds_neighborhood(rec, cdss, focus_index, cds_range=cds_range, kb_range=kb_range, whole_contig=whole_contig, normalize_direction=normalize_direction, _from_domain_search=_from_domain_search, _feature_index=feature_index)
                if _from_domain_search:
                    set_hit_best_score(record, best_score)
                
                if max_region_overlap < 1.0:
                    loc_len = len(record_location)
                    for other_location in selected_locations:
                        if (overlap_size(record_location, other_location) / loc_len) > max_region_overlap:
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


def build_parser():
    parser = tool_parser(__doc__)

    add_input(parser, "GenBank files.")

    add_seq_outputs(parser)
    add_taxonomy_filters(parser)


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
    parser.add_argument('--name_field', type=str, default=None,
                        help="CDS qualifier to prefer for name matching, falling back to gene_id, locus_tag, gene, then protein_id. Does not change coordinate-based region IDs.")
    parser.add_argument('--domains', type=str, default=None, nargs='+',
                        help="Include domains with these names") #TODO: add some way to specify the database of a domain
    parser.add_argument('--domains_file', type=str, default=None,
                        help="text file containing names of domains to include.")
    parser.add_argument('--domain_expr', type=str, default=None,
                        help="a boolean expression using operators & (AND), | (OR), and ~(NOT), to specify a desired combination of domains")
    
    parser.add_argument("--databases", default=None, required=False, type=str, nargs="+",
                        help="Consider only domains from these databases when filtering by domain. default: all databases.")

    parser.add_argument('--strand', type=str, default=None, choices=["f", "r"], help="Only extract regions around CDSs on the specified strand.")
    parser.add_argument('--partial', type=str, default="include", choices=PARTIAL_CHOICES, help=partial_help("Which CDSs to extract regions around") + " Not affected by --invert.")

    parser.add_argument('--unannotated', action='store_true', default=False, help="Select CDSs with no domain annotations.")

    add_annotation_type(parser, help="Which annotations to select CDSs by.")
    add_has_search_hit(parser)

    parser.add_argument('--pad', action='store_true', default=False,
                        help="If set, then ends of sequences will be padded so that the focus CDS aligns. See --anchor.")
    parser.add_argument('--anchor', type=str, default="start", choices=WHERE_CHOICES,
                        help="Where to align the focus CDSs when --pad is set, relative to the focus CDS in its own orientation. "
                        "start: align the first base of the focus CDSs. end: align their last base. center: align their midpoints.")

    parser.add_argument('--keep_direction', action='store_true', default=False,
                        help="by default extracted regions will be flipped so that the focus cds is on the forward strand. Setting this option will keep the focus cds on whatever strand it started on.")
    add_dedup(parser, default="id")
    parser.add_argument('--max_region_overlap', type=float, default=1.0,
                       help="the maximum fractional of overlap between any two output regions. If >= 1, then no overlap filtering will be done. Regions are output in a greedy fashion based on CDS start site. New regions are output if less than this fraction of them overlaps with any previously output region.")


    parser.add_argument('-e', '--evalue', default=None, type=PositiveFloat,
                        help="the evalue cutoff (max) for domain annotations to be considered when selecting CDSs. Default: no cutoff.")

    parser.add_argument('--invert', action='store_true',
                        help="Invert the combined CDS selectors, not contigs, strand, partial, or taxonomy. With no selectors, warn and ignore.")

    parser.add_argument('--include_nucleic_acids', action='store_true',
                        help="Include nucleic-acid annotations with no associated CDS (those with cds_id='.') when applying CDS selection criteria.")


    parser.add_argument('--config', action=ActionConfigFile)
    

    declare_outputs(parser, ["genbank_out", "fasta_out"])
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)

    ### validate input


    check_outputs(parser, params)

    ### figure out selection filters
    kb_range, cds_range = validate_range(params.kb_range, params.kb_up, params.kb_down, params.cds_range, params.cds_up, params.cds_down)

    whole_contig = params.whole_contig
    normalize_direction = not params.keep_direction

    target_contigs = list_and_file_to_dict_keys(params.contigs, params.contigs_file, as_set=True)
    selected_cdss = list_and_file_to_dict_keys(params.cds, params.cds_file, as_set=True)
    selected_domains = list_and_file_to_dict_keys(params.domains, params.domains_file, as_set=True)
    domain_expr = params.domain_expr


    ### Run
    selected = select_by_cds(
        parse_seqfiles(params.input, target_contigs), 
        cds_ids=selected_cdss, 
        name_field=params.name_field,
        domains=selected_domains, 
        domain_expr=domain_expr, 
        cds_range=cds_range, 
        kb_range=kb_range, 
        whole_contig=whole_contig, 
        normalize_direction=normalize_direction,
        evalue=params.evalue,
        annotation_type = params.annotation_type,
        has_search_hit=params.has_search_hit,
        invert = params.invert,
        max_region_overlap = params.max_region_overlap,
        strand = params.strand,
        databases = params.databases,
        unannotated = params.unannotated,
        include_nucleic_acids = params.include_nucleic_acids,
        partial = params.partial,
        ncbi_taxonomy=taxonomy_from_args(params),
        include_taxids=params.include_taxids,
        exclude_taxids=params.exclude_taxids,
        taxonomy_expr=params.taxonomy_expr,
        taxonomy_unknown=params.taxonomy_unknown,
        )
    with SeqOutputs(params.genbank_out, params.fasta_out) as outputs:
        outputs.write_all(deduplicate_and_pad(selected, dedup=params.dedup, pad=params.pad, anchor=params.anchor))
