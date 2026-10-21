"""Slicing records and padding locations, with their features carried over.
"""
from typing import Dict, Iterable, List, Tuple, Union
from collections import defaultdict
import functools

from domainator import DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation, ExactPosition, BeforePosition, AfterPosition
from Bio.SeqRecord import SeqRecord
from domainator.partial import codon_offset, location_partial_ends, on_contig_translation
from domainator.location_ops import circular_envelope, merge_parts, stranded_end, stranded_start


def regions_overlap(region1, region2, min_overlap_fraction=0.0):
    """
        regions are tuples of start and stop coordinates
        returns true if a fraction of region2 >= min_overlap_fraction overlaps with region1
        coordinates within regions must be sorted low to high
    """
    if min_overlap_fraction >= 1:
        return False

    # TODO: does biopython coordinate system require us to adjust by 1 somewhere?
    region_1_size = region1[1] - region1[0]
    region_2_size = region2[1] - region2[0]
    if region_1_size == 0:
        region_1_size = 0.1
    if region_2_size == 0:
        region_2_size = 0.1

    r1_s, r1_e = region1
    r2_s, r2_e = region2
    
    overlap_size = max(0, min(r1_e, r2_e) - max(r1_s, r2_s))
    if overlap_size == 0:
        return False
    else:
        return (overlap_size / region_2_size) >= min_overlap_fraction

def overlap_fraction(region1, region2) -> float:
    """The fraction of region2 that overlaps region1. Regions are (start, end) tuples with start <= end."""
    overlap_size = max(0, min(region1[1], region2[1]) - max(region1[0], region2[0]))
    if overlap_size == 0:
        return 0.0
    return overlap_size / (region2[1] - region2[0])


def slice_record(r:SeqRecord, start:int, stop:int, features:Iterable[SeqFeature]=None, truncate_features=True) -> SeqRecord:
    """ Adapted from Biopython SeqRecord::SeqRecord::__getitem__
        
        Unlike __getitem__, this function keeps 'source' features (resizing them to the size of the slice),
            as well as other metadata such as dbxrefs, contig 'source', 'organism', and 'taxonomy'.
        
        It also allows you to specify a subset of features to slice, instead of all features.

        Args:
            r (SeqRecord): The record to slice
            start (int): The start coordinate of the slice
            stop (int): The stop coordinate of the slice
            features (Iterable[SeqFeature]): The features to slice. If None, then all features will be sliced.
            truncate_features (bool): If True, then features that are partially within the slice will be truncated to fit within the slice.
                If False, then features that are partially within the slice will be omitted from the output record.
                source features will always be truncated, regardless of this setting.

        #TODO: maybe make this a method of SeqRecord, like make a slice method and have an option to keep source.
    """

    #TODO: avoid code duplication with slice_record_from_location, by making this function a wrapper for that one, as it is a special case.

    return slice_record_from_location(r, FeatureLocation(start, stop), features, truncate_features=truncate_features)

def pad_location(r:SeqRecord, location:Union[FeatureLocation,CompoundLocation], pad_upstream:int = 0, pad_downstream:int = 0, return_extension_sizes:bool = False) -> Union[Union[FeatureLocation,CompoundLocation],Tuple[Union[FeatureLocation,CompoundLocation],Tuple[int,int]]]:
    """ Pad a location by a specified amount of bases on either side.
        If the location is a CompoundLocation, then the first and last subparts will be padded.
        Padding will take strand information into account.
            if the first subpart is on the reverse strand, then upstream padding will be added to its end coordinate, 
            and if the last subpart is on the reverse strand, then downstream padding will be added to its start coordinate.
        
        Padding will also take into account the fact that the location might extend beyond the end of the sequence.
            If the padded location extends beyond the end of the sequence, then the padding will be reduced to the maximum possible padding.
        
        Finally, padding will take into account circular sequences.
            If padding extends beyond the origin, then padding will be extended beyond the origin, with new parts added to the location as necessary.
            If upstream padding overlaps with downstream padding, then the overlapping region will be divided equally between the two parts.

        Args:
            r (SeqRecord): The record to which the location belongs
            location (Union[FeatureLocation,CompoundLocation]): The location to pad
            pad_upstream (int): The number of bases/amino acids to pad upstream
            pad_downstream (int): The number of bases/amino acids to pad downstream
        
        Returns:
            Union[FeatureLocation,CompoundLocation]: The padded location
    """
    if pad_upstream < 0 or pad_downstream < 0:
        raise ValueError("Padding cannot be negative")
    if pad_upstream == 0 and pad_downstream == 0:
        if return_extension_sizes:
            return location, (0, 0)
        else:
            return location
    
    rec_len = len(r)
    if r.annotations.get("topology", None) != "circular":
        if len(location.parts) == 1:
            location = location.parts[0]
            if location.strand == -1:
                min_coord = max(0, int(location.start)-pad_downstream)
                downstream_extension = int(location.start)-min_coord
                max_coord = min(rec_len, int(location.end)+pad_upstream)
                upstream_extension = max_coord-int(location.end)

            else: # location.strand == 1 or location.strand == 0
                min_coord = max(0, int(location.start)-pad_upstream)
                upstream_extension = int(location.start)-min_coord
                max_coord = min(rec_len, int(location.end)+pad_downstream)
                downstream_extension = max_coord-int(location.end)
            out_location = FeatureLocation(
                    min_coord,
                    max_coord,
                    strand=location.strand
                )
            if return_extension_sizes:
                return out_location, (upstream_extension, downstream_extension)
            else:
                return out_location
        else: # CompoundLocation with multiple parts
            # in the weird case where the first part and last part are on different strands or where a middle part is beyond the bounds of the outer parts, this code could result in covering the same regions multiple times.
            first_part, upstream_extensions = pad_location(r, location.parts[0], pad_upstream, 0, return_extension_sizes=True)
            last_part, downstream_extensions = pad_location(r, location.parts[-1], 0, pad_downstream, return_extension_sizes=True)
            new_parts = [first_part]
            for part in location.parts[1:-1]:
                new_parts.append(part)
            new_parts.append(last_part)
            out_location = CompoundLocation(new_parts, operator=location.operator)
            if return_extension_sizes:
                return out_location, (upstream_extensions[0], downstream_extensions[1])
            else:
                return out_location
    else: # contig is circular
        if len(location.parts) == 1: # simple location (one segment)
            location = location.parts[0]
            l_start = int(location.start)
            l_end = int(location.end)
            #start_to_origin = int(location.start)
            end_to_origin = rec_len-int(location.end)
            
            # find the bounds of the region to be extracted
            if location.strand == -1:
                min_coord = location.start-pad_downstream
                max_coord = location.end+pad_upstream
            else: # location.strand == 1 or location.strand == 0
                min_coord = location.start-pad_upstream
                max_coord = location.end+pad_downstream
            
            if min_coord >= 0 and max_coord <= rec_len:
                # the bounds to not extend beyond the origin, on either side, so  we don't care if they overlap
                out_location = FeatureLocation(
                    min_coord,
                    max_coord,
                    strand=location.strand
                )
                if return_extension_sizes:
                    return out_location, (pad_upstream, pad_downstream)
                else:
                    return out_location
            else:
                # If regions extending beyond the origin run into the FeatureLocation, truncate them there.
                # We prefer to not create a copy of any part of the contig.
                
                # adjust min_coord and max_coord so that they don't wrap around the FeatureLocation bounds.
                if min_coord < 0 and ((-1 * min_coord) > end_to_origin):
                    min_coord = -1 * end_to_origin
                if max_coord > rec_len and (max_coord - rec_len > l_start):
                    max_coord = rec_len + l_start
                
                # If the adjusted regions overlap with each other, then divide the overlap equally between min and max.
                # shift the perspective so that min_coord is at the origin
                min_coord_shifted = 0 # min_coord_shifted = min_coord - min_coord
                max_coord_shifted = max_coord - min_coord
                
                # calculuate overlap size:
                overlap_size = max_coord_shifted - rec_len
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

                if location.strand == -1:
                    upstream_extension = max_coord-l_end
                    downstream_extension = l_start-min_coord 
                else: # location.strand == 1 or location.strand == 0
                    upstream_extension = l_start-min_coord
                    downstream_extension = max_coord-l_end 


                parts = list() # we're touching or wrapping around the origin, so we may have multiple parts to extract
                if min_coord < 0: # min_coord is across the origin from the start of the FeatureLocation, so we need to break it up into two parts
                    parts.append(FeatureLocation(rec_len+min_coord, rec_len, strand=location.strand)) #from the start of the extraction window to the origin
                    if l_end > 0: # not touching the origin, so we need to add a part for the region from the origin to the end of the FeatureLocation
                        parts.append(FeatureLocation(0, l_end, strand=location.strand))
                elif min_coord != l_start: # min_coord is not across the origin from the start of the FeatureLocation
                    parts.append(FeatureLocation(min_coord, l_end, strand=location.strand))
                else: # min_coord == l_start
                    parts.append(FeatureLocation(l_start, l_end, strand=location.strand))

                if max_coord > rec_len: # max_coord is across the origin from the end of the FeatureLocation, so we need to break it up into two parts
                    parts[-1] = FeatureLocation(parts[-1].start, rec_len, strand=location.strand)
                    parts.append(FeatureLocation(0, max_coord-rec_len, strand=location.strand))
                else:
                    parts[-1] = FeatureLocation(parts[-1].start, max_coord, strand=location.strand)

                if location.strand == -1:
                    parts.reverse()
                
                if len(parts) == 1:
                    out_location = parts[0]
                else:
                    out_location = CompoundLocation(parts, operator="join")

                if return_extension_sizes:
                    return out_location, (upstream_extension, downstream_extension)
                else:
                    return out_location
        else: # CompoundLocation with more than one part
            # we treat the CompoundLocation as a single FeatureLocation, so we need to calculate the bounds of the CompoundLocation
            # we consider the Strand of the CompoundLocation to be the Strand of the first part
            # TODO: handle the case where the first part is on the opposite strand of the last part.
            # TODO: handle the case where a middle part is beyond the bounds of the outer parts
            
            part_1_strand = location.parts[0].strand
            envelope = circular_envelope(location, len(r), part_1_strand)
            if len(envelope) == len(r): # envelope covers the entire contig, so just return the starting location
                return location
            elif len(envelope.parts) == 1: # envelope is a single FeatureLocation (i.e. it does not cross the origin, so the padding might cross the origin)
                _, extensions = pad_location(r, envelope.parts[0], pad_upstream, pad_downstream, return_extension_sizes=True)
                
                out_parts = list()
                start_parts = pad_location(r, location.parts[0], extensions[0], 0)
                end_parts = pad_location(r, FeatureLocation(location.parts[-1].start, location.parts[-1].end, part_1_strand), 0, extensions[1])
                out_parts.extend(start_parts.parts)
                out_parts.extend(location.parts[1:-1])
                out_parts.extend(end_parts.parts)
                out_location = CompoundLocation(out_parts, operator="join")
                if return_extension_sizes:
                    return out_location, extensions
                else:
                    return out_location
            else: # envelope crosses the origin, implying that the padding will not cross the origin
                # adjust pad_upstream and pad_downstream so that they won't run into the envelope.
                l_start = stranded_start(envelope)
                l_end = stranded_end(envelope) 
                remaining_span = abs(l_start - l_end) # size of the record minus the size of the envelope
                pad_upstream = min(pad_upstream, remaining_span) # adjust pad_upstream so that it won't run into the envelope
                pad_downstream = min(pad_downstream, remaining_span) # adjust pad_downstream so that it won't run into the envelope
                padding_size = pad_upstream + pad_downstream
                if padding_size > remaining_span:
                    padding_size_diff = padding_size - remaining_span # 
                    # divide the padding_size_diff equally between pad_upstream and pad_downstream
                    pad_upstream -= padding_size_diff // 2
                    pad_downstream -= padding_size_diff // 2
                    # if padding_size_diff is odd, then we need to adjust one of the coordinates by 1
                    if padding_size_diff % 2 == 1:
                        pad_downstream -= 1
                out_parts = list()
                first_part, upstream_extensions = pad_location(r, location.parts[0], pad_upstream, 0, return_extension_sizes=True)
                last_part, downstream_extensions = pad_location(r, FeatureLocation(location.parts[-1].start, location.parts[-1].end, part_1_strand), 0, pad_downstream, return_extension_sizes=True)
                out_parts.extend(first_part.parts) # first part should only have one part
                out_parts.extend(location.parts[1:-1])
                out_parts.extend(last_part.parts) # last part should only have one part
                out_location = CompoundLocation(out_parts, operator="join")
                if return_extension_sizes:
                    return out_location, (upstream_extensions[0], downstream_extensions[1])
                else:
                    return out_location

                
def _merge_intervals(intervals:List[Tuple[int,int]]) -> List[Tuple[int,int]]:
    """
        Merges touching or overlapping (start, end) intervals. Returns them sorted by start.
    """
    merged = list()
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged

def _truncate_cds_qualifiers(feature:SeqFeature, kept_intervals:List[Tuple[int,int]]) -> Dict[str, List]:
    """
        Returns a copy of a CDS feature's qualifiers, updated for a CDS that was cut by slicing.

        codon_start is updated so that the reading frame is preserved, and the translation is trimmed to the residues
        whose codons are completely kept. If the kept part of the CDS is not contiguous (e.g. an interior exon was cut out),
        or the translation can't be placed on the CDS, the translation is removed, so that it will be re-translated.

        Args:
            feature: the original (uncut) CDS feature
            kept_intervals: the kept parts of the CDS, as merged (start, end) intervals along the CDS, 5' to 3'
    """
    qualifiers = feature.qualifiers.copy()
    offset = codon_offset(feature)
    kept_start, kept_end = kept_intervals[0][0], kept_intervals[-1][1]
    new_codon_start = ((offset - kept_start) % 3) + 1
    if new_codon_start != 1 or "codon_start" in qualifiers:
        qualifiers["codon_start"] = [str(new_codon_start)]

    translation = qualifiers.get("translation", [""])[0]
    if len(translation) != 0:
        five_prime, three_prime = location_partial_ends(feature.location)
        translation = on_contig_translation(str(translation), len(feature.location), offset, five_prime, three_prime)
    if len(kept_intervals) != 1 or translation is None:
        qualifiers.pop("translation", None)
    elif len(translation) != 0:
        first_residue = max(0, -(-(kept_start - offset) // 3)) # first codon that starts at or after kept_start
        last_residue = max(0, (kept_end - offset) // 3) # codons that end at or before kept_end
        qualifiers["translation"] = [translation[first_residue:last_residue]]
    return qualifiers

class FeatureIndex:
    """The features of a record, binned by position, for cutting many regions out of one record.

    near(location) returns the features whose outer bounds touch the location, in their original
    order: a superset of the ones slice_record_from_location keeps, found without scanning them all.
    """
    BIN_SIZE = 4096

    def __init__(self, features: Iterable[SeqFeature]):
        self.features = list(features)
        self.bins = defaultdict(list)
        for i, feature in enumerate(self.features):
            for b in range(int(feature.location.start) // self.BIN_SIZE, int(feature.location.end) // self.BIN_SIZE + 1):
                self.bins[b].append(i)

    def near(self, location) -> List[SeqFeature]:
        found = set()
        for part in location.parts:
            for b in range(int(part.start) // self.BIN_SIZE, int(part.end) // self.BIN_SIZE + 1):
                found.update(self.bins.get(b, ()))
        return [self.features[i] for i in sorted(found)]


def slice_record_from_location(r:SeqRecord, location:Union[FeatureLocation,CompoundLocation], features:Iterable[SeqFeature]=None, pad_upstream:int = 0, pad_downstream:int = 0, truncate_features=True) -> SeqRecord:
    """Extract a sub-record from supplied parent record using the location of the SeqFeature object.
        Output record will have the same id, name, description, dbxrefs, and annotations as the parent record.
        Output record will have the same features as the parent record, truncated to the location (see truncate_features).
        Output record will have the same letter annotations as the parent record, except that the letter annotations will be sliced to match the sliced sequence.

        sub-parts of the location will be extracted from the parent record, and then concatenated together to form the output record.
        If the location is a CompoundLocation, then the sub-parts will be extracted in the order they are listed in the CompoundLocation.
        parts on the reverse strand will be extracted in reverse order, and then reversed again to be in the correct orientation.


        

        Adapted from Biopython SeqRecord::SeqRecord::__getitem__
        
        Unlike __getitem__, this function keeps 'source' features (resizing them to the size of the slice),
            as well as other metadata such as dbxrefs, contig 'source', 'organism', and 'taxonomy'.
        
        It also allows you to specify a subset of features to slice, instead of all features.

        r (SeqRecord): parent record
        location (FeatureLocation): location to extract
        features (Iterable[SeqFeature]): subset of features to slice. If None, all features will be sliced.
        pad_upstream (int): number of bases/residues to pad upstream of the location
        pad_downstream (int): number of bases/residues to pad downstream of the location
        truncate_features (bool): if True, then features that are not completely contained in the location will be
            included in the output record, truncated to the size of the location, with '<' or '>' markers on the cut ends.
            Cut CDSs get an updated codon_start and a translation trimmed to the kept codons. Cut Domainator hit features are dropped.
            If False, then features that are not completely contained in the location will be removed from the output record.
            "source" features will always be included (without new markers), regardless of this setting.

        Returns:
            SeqRecord: sliced record

        TODO: should padding be handled here, or should we expect the FeatureLocation to already be padded?
    """

    # TODO: to deal with padding 
    #       linear contigs:   add padding to the first or last part of the location, depending on the strand, and slice the record from there.
    #       circular contigs: add padding to the first and last part of the location, depending on the strand, truncate the padding if there is overlap, 
    #                           and then slice the record from there. 
    #                           No rotation should be necessary!, for example if you are slicing a kb-range across the origin, you just need to define a compound 
    #                           location with one part on either side of the origin, and then add padding to both parts, and then slice the record from there. No rotation needed.)

    if r.seq is None:
        raise ValueError("If the sequence is None, we cannot slice it.")

    if features is None:
        features = r.features

    if pad_upstream < 0 or pad_downstream < 0:
        raise ValueError("Padding must be a positive integer.")
    if pad_upstream != 0 or pad_downstream != 0:
        location = pad_location(r, location, pad_upstream, pad_downstream)
    # get sequence parts
    seq_parts = list()
    seq_parts_len = 0
    new_features = list()
    coord_shifts = list()
    strands = list()
    # location_parts = merge_location_parts(location.parts) #TODO: something like this might fix edge cases where an annotation is completely enclosed by the location, but for some reason the location has a size-zero split with the overlap region.
    for part in location.parts:
        strand = part.strand
        strands.append(strand)
        seq_parts.append(part.extract(r.seq))

        # If the part is on the reverse strand, then to calculate shift for overlaping parts,
        # we need to shift left by the length of the part, flip, then add to the length of the previous parts.
        # 
        # if the part is on the forward strand, then to calculate shift for overlaping parts,
        # we need to shift left by the length of the part, then add to the length of the previous parts.
        # 
        # in either case, we need to add the length of the previous parts to the shift.

        coord_shifts.append(seq_parts_len)         
        
        seq_parts_len += len(seq_parts[-1])
    
    is_protein = r.annotations.get("molecule_type") == "protein"
    # include features that overlap the location. Features cut by the edge of the location are truncated, and the cut ends get '<' or '>' markers.
    # Domain hit features that are cut are dropped, because their scores and coordinates describe the whole hit.
    loc_bounds = [(int(part.start), int(part.end)) for part in location.parts]
    for f in features:
        # Skip features nowhere near the location (the outer bounds of a feature cover all its parts).
        f_outer_start, f_outer_end = int(f.location.start), int(f.location.end)
        if not any(f_outer_start <= loc_end and f_outer_end >= loc_start for loc_start, loc_end in loc_bounds):
            continue
        is_hit_feature = f.type == DOMAIN_FEATURE_NAME or f.type == DOMAIN_SEARCH_BEST_HIT_NAME
        f_parts = f.location.parts
        overlap = [set() for _ in range(len(f_parts))] # list of sets of indices of location parts that overlap with the feature part
        total_overlap = 0
        for feature_part_index, f_part in enumerate(f_parts):
            part_found = False
            f_start = int(f_part.start)
            f_end = int(f_part.end)
            for loc_part_index, (loc_start, loc_end) in enumerate(loc_bounds):
                if ( # part is completely contained in loc_part
                    loc_start <= f_start
                    and loc_end >= f_end
                ):
                    overlap[feature_part_index].add(loc_part_index)
                    part_found = True
                elif f.type == "source" or truncate_features: # part is partially contained in loc_part
                    if regions_overlap((f_start, f_end), (loc_start, loc_end)):
                        overlap[feature_part_index].add(loc_part_index)
                        part_found = True
            if part_found:
                total_overlap += 1
        if total_overlap == len(f_parts) or f.type == "source" or truncate_features:
            # Collect the kept pieces of the feature: one per (feature part, location part) overlap.
            # Each piece also gets its interval in coordinates along the feature (5' to 3'), to find the ends that are really cut,
            # as opposed to pieces that continue in another location part (e.g. a feature split across the origin of a rotated contig).
            pieces = list()
            transcript_offset = 0
            for feature_part_index, f_part in enumerate(f_parts):
                f_start = int(f_part.start)
                f_end = int(f_part.end)
                for loc_part_index in sorted(overlap[feature_part_index]):
                    loc_part = location.parts[loc_part_index]
                    kept_start = max(f_start, int(loc_part.start))
                    kept_end = min(f_end, int(loc_part.end))
                    if f_part.strand == -1:
                        transcript_interval = (transcript_offset + f_end - kept_end, transcript_offset + f_end - kept_start)
                    else:
                        transcript_interval = (transcript_offset + kept_start - f_start, transcript_offset + kept_end - f_start)
                    pieces.append((transcript_interval, f_part, loc_part_index, kept_start, kept_end))
                transcript_offset += f_end - f_start
            transcript_intervals = _merge_intervals([piece[0] for piece in pieces])
            is_cut = sum(end - start for start, end in transcript_intervals) < transcript_offset
            if is_cut and is_hit_feature:
                continue
            cut_5p_ends = {start for start, end in transcript_intervals if start > 0}
            cut_3p_ends = {end for start, end in transcript_intervals if end < transcript_offset}

            new_parts = list()
            for (t_start, t_end), f_part, loc_part_index, kept_start, kept_end in sorted(pieces, key=lambda piece: piece[0]): # in order along the feature
                loc_part = location.parts[loc_part_index]
                loc_start = int(loc_part.start)
                loc_end = int(loc_part.end)
                if f_part.strand == -1:
                    start_is_cut, end_is_cut = t_end in cut_3p_ends, t_start in cut_5p_ends
                else:
                    start_is_cut, end_is_cut = t_start in cut_5p_ends, t_end in cut_3p_ends
                # Ends of the original feature keep their position type (e.g. an existing '<' or '>').
                # Ends that are really cut get a new '<' or '>', except on source features.
                if kept_start == int(f_part.start):
                    new_start = f_part.start + (-loc_start)
                elif start_is_cut and f.type != "source":
                    new_start = BeforePosition(kept_start - loc_start)
                else:
                    new_start = ExactPosition(kept_start - loc_start)
                if kept_end == int(f_part.end):
                    new_end = f_part.end + (-loc_start)
                elif end_is_cut and f.type != "source":
                    new_end = AfterPosition(kept_end - loc_start)
                else:
                    new_end = ExactPosition(kept_end - loc_start)
                new_loc = FeatureLocation(new_start, new_end, strand=f_part.strand)
                if strands[loc_part_index] == -1:
                    new_loc = new_loc._flip(loc_end-loc_start)
                new_loc = new_loc._shift(coord_shifts[loc_part_index])
                new_parts.append(new_loc)
            # Merge parts that are adjacent or overlap
            if getattr(f.location, "operator", None) == "join" or getattr(f.location, "operator", None) is None: # only merge parts if the original location was a join or not a CompoundLocation (for example, an origin-spanning source) TODO: should we merge parts for any other operators?
                new_parts = merge_parts(new_parts) # TODO: There are some edge cases where this will merge things that maybe users won't want to merge, but it's better than not merging at all. Can fix later if anyone complains.

            if len(new_parts) == 0: # this can happen for source features that don't intersect with the location
                pass
            else:
                if len(new_parts) == 1:
                    new_location = new_parts[0]
                elif len(new_parts) > 1:
                    new_operator = getattr(f.location, "operator", "join")
                    new_location = CompoundLocation(new_parts, new_operator) 
                qualifiers = f.qualifiers
                if is_cut and f.type == "CDS" and not is_protein:
                    qualifiers = _truncate_cds_qualifiers(f, transcript_intervals)
                new_feature = f.__class__(
                    new_location,
                    type=f.type,
                    id=f.id,
                    qualifiers=qualifiers,
                )
                new_features.append(new_feature)
    # TODO: handle other partially cut annotations


    answer = r.__class__(
        functools.reduce(lambda x, y: x + y, seq_parts),
        id=r.id,
        name=r.name,
        description=r.description,
        dbxrefs=r.dbxrefs.copy(),
        features=new_features,
    )

    if "molecule_type" in r.annotations:
        # This will still apply, and we need it for GenBank/EMBL etc output
        answer.annotations["molecule_type"] = r.annotations["molecule_type"]

    if "source" in r.annotations:
        answer.annotations["source"] = r.annotations["source"]

    if "organism" in r.annotations:
        answer.annotations["organism"] = r.annotations["organism"]

    if "taxonomy" in r.annotations:
        answer.annotations["taxonomy"] = r.annotations["taxonomy"]

    # # if "references" in r.annotations: #TODO: maybe keep this also?
    # #     answer.annotations["references"] = r.annotations["references"]

    # TODO: handle letter_annotations (genbank files don't have these, so probably not relevant for Domainator)

    return answer
# Positions relative to a feature in its own orientation, e.g. for the --anchor option of rotate_contigs, select_by_cds and domain_search.
WHERE_CHOICES = ("start", "end", "center")


def feature_intervals(feature):
    """(start, end) half-open intervals for a Biopython SeqFeature (``.location.parts``) or a
    lean_record LeanFeature (``.parts`` tuples). Works on both so taxonomy resolution applies
    to materialized SeqRecords and lean fallback records alike."""
    loc = getattr(feature, "location", None)
    if loc is not None:
        return [(int(p.start), int(p.end)) for p in loc.parts]
    return [(int(p[0]), int(p[1])) for p in feature.parts]  # LeanFeature parts


def location_covers(outer_intervals, inner_intervals) -> bool:
    """True iff every `inner` interval is contained in the union of `outer` intervals.

    Intervals are (start, end) half-open, 0-based (the Biopython/lean convention). Origin-
    spanning locations are already split into multiple parts upstream, so per-interval
    containment is sufficient; `outer` is merged to handle adjacent/overlapping parts.
    """
    outer = _merge_intervals([(int(s), int(e)) for s, e in outer_intervals])
    return all(
        any(os <= int(s) and int(e) <= oe for os, oe in outer)
        for s, e in inner_intervals
    )

def circular_dist(start, end, size):
    if end < start:
        end = end + size
    return end - start
