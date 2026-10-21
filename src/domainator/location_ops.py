"""Arithmetic on Biopython feature locations (SimpleLocation and CompoundLocation).

These were methods of Domainator's Biopython fork. As functions they work on ordinary
Biopython locations, using only their public attributes.
"""
from typing import Iterable, List, Optional, Union

from Bio.SeqFeature import CompoundLocation, SimpleLocation

AnyLocation = Union[SimpleLocation, CompoundLocation]


def _parts(location: AnyLocation) -> List[SimpleLocation]:
    return location.parts if isinstance(location, CompoundLocation) else [location]


def copy_location(location: AnyLocation) -> AnyLocation:
    """A copy of ``location``, keeping its (possibly fuzzy) positions."""
    if isinstance(location, CompoundLocation):
        parts = [SimpleLocation(p.start, p.end, p.strand, p.ref, p.ref_db) for p in location.parts]
        return CompoundLocation(parts, location.operator)
    return SimpleLocation(location.start, location.end, location.strand, location.ref, location.ref_db)


def stranded_start(location: AnyLocation):
    """The position where the location starts in its own orientation: the first part's end
    on the minus strand, otherwise its start. Possibly a fuzzy position."""
    first = _parts(location)[0]
    return first.end if first.strand == -1 else first.start


def stranded_end(location: AnyLocation):
    """The position where the location ends in its own orientation: the last part's start
    on the minus strand, otherwise its end. Possibly a fuzzy position."""
    last = _parts(location)[-1]
    return last.start if last.strand == -1 else last.end


def stranded_start_human_readable(location: AnyLocation) -> int:
    """stranded_start as a 1-based coordinate."""
    first = _parts(location)[0]
    return int(first.end) if first.strand == -1 else int(first.start) + 1


def stranded_end_human_readable(location: AnyLocation) -> int:
    """stranded_end as a 1-based coordinate."""
    last = _parts(location)[-1]
    return int(last.start) + 1 if last.strand == -1 else int(last.end)


def overlay(location: AnyLocation, offset: int, length: int) -> AnyLocation:
    """The part of ``location`` that starts ``offset`` positions into it (in its own
    orientation, across the parts of a compound location) and is ``length`` long."""
    if offset < 0:
        raise ValueError("offset must be >= 0")
    if length < 0:
        raise ValueError("length must be >= 0")
    if offset + length > len(location):
        raise ValueError("offset + length must be <= len(location)")
    if not isinstance(location, CompoundLocation):
        if location.strand == -1:
            return SimpleLocation(location.end - offset - length, location.end - offset, location.strand)
        return SimpleLocation(location.start + offset, location.start + offset + length, location.strand)
    new_parts = []
    for part in location.parts:
        part_len = len(part)
        if offset >= part_len:
            offset -= part_len
            continue
        section_len = min(part_len - offset, length)
        new_parts.append(overlay(part, offset, section_len))
        length -= section_len
        offset = 0
        if length == 0:
            break
    if len(new_parts) == 1:
        return new_parts[0]
    return CompoundLocation(new_parts, location.operator)


def touches(location: AnyLocation, other: AnyLocation) -> bool:
    """True if the locations overlap or touch end to end."""
    for part in _parts(location):
        for other_part in _parts(other):
            s_start, s_end = int(part.start), int(part.end)
            o_start, o_end = int(other_part.start), int(other_part.end)
            if o_start <= s_start <= o_end or o_start <= s_end <= o_end or s_start <= o_start <= s_end or s_start <= o_end <= s_end:
                return True
    return False


def overlap_size(location: AnyLocation, other: AnyLocation) -> float:
    """The number of positions the locations share, summed over their parts."""
    size = 0.0
    for part in _parts(location):
        for other_part in _parts(other):
            size += max(0.0, min(int(part.end), int(other_part.end)) - max(int(part.start), int(other_part.start)))
    return float(size)


def overlaps(location: AnyLocation, other: AnyLocation) -> bool:
    """True if the locations share at least one position (just touching doesn't count)."""
    return overlap_size(location, other) > 0


def contains(location: AnyLocation, other: AnyLocation) -> bool:
    """True if every part of ``other`` lies inside one part of ``location``."""
    return all(any(int(part.start) <= int(o.start) and int(part.end) >= int(o.end) for part in _parts(location))
               for o in _parts(other))


def merge_parts(parts: Iterable[SimpleLocation]) -> List[SimpleLocation]:
    """Merge consecutive parts on the same strand that touch or overlap.

    The outer positions are kept, so '<' and '>' markers on the outer ends survive the merge.
    Parts are not sorted, which would scramble the connectivity of a join.
    """
    parts = list(parts)
    if not parts:
        return parts
    merged = [parts[0]]
    for part in parts[1:]:
        if part.strand == merged[-1].strand and touches(part, merged[-1]):
            merged[-1] = SimpleLocation(min(merged[-1].start, part.start, key=int), max(merged[-1].end, part.end, key=int), strand=part.strand)
        else:
            merged.append(part)
    return merged


def circular_envelope(location: AnyLocation, rec_len: int, strand: Optional[int] = None) -> AnyLocation:
    """The smallest gapless location that contains ``location`` and every junction between its parts.

    rec_len: the length of the record the location is on (to tell whether the location wraps
        around the origin).
    strand: the strand of the envelope; by default the strand of the first part.

    If the location does not span the origin, or covers the entire record, the envelope is one
    SimpleLocation; if it spans the origin, it is a CompoundLocation of two parts.

    Assumes a circular record; for a linear one, SimpleLocation(f.start, f.end, f.strand) is
    equivalent. Built for prokaryotic CDSs spanning the origin, so it might not handle edge
    cases that are rare or irrelevant for CDSs the way you expect: unstranded parts are treated
    as forward-strand, and parts that overlap their neighbors are sometimes treated as though an
    origin cross is needed to reach them.
    """
    if not isinstance(location, CompoundLocation):
        return SimpleLocation(location.start, location.end, location.strand if strand is None else strand)
    parts = location.parts
    if len(parts) == 0:
        return SimpleLocation(0, 0, strand=0)
    if len(parts) == 1:
        part = parts[0]
        return SimpleLocation(part.start, part.end, strand=part.strand if strand is None else strand)
    if strand is None:
        strand = parts[0].strand
    loc_parts = [copy_location(parts[0])]
    for next_part in parts[1:]:
        part = loc_parts[-1]
        this_strand, this_start, this_end = part.strand, int(part.start), int(part.end)
        next_strand, next_start, next_end = next_part.strand, int(next_part.start), int(next_part.end)
        if touches(part, next_part):
            loc_parts[-1] = SimpleLocation(min(this_start, next_start), max(this_end, next_end), strand=next_strand)
        elif this_strand == -1:
            # going backward: reaching the next part without crossing the origin means it isn't past this start
            reachable = next_end <= this_start if next_strand == -1 else next_start <= this_start
            if reachable:
                loc_parts[-1] = SimpleLocation(min(this_start, next_start), max(this_end, next_end), strand=next_strand)
            else:
                if this_start != 0:
                    loc_parts[-1] = SimpleLocation(0, this_end, this_strand)
                loc_parts.append(SimpleLocation(next_start, rec_len, next_strand))
        else:
            reachable = next_end >= this_end if next_strand == -1 else next_start >= this_end
            if reachable:
                loc_parts[-1] = SimpleLocation(min(this_start, next_start), max(this_end, next_end), strand=next_strand)
            else:
                if this_end != rec_len:
                    loc_parts[-1] = SimpleLocation(this_start, rec_len, this_strand)
                loc_parts.append(SimpleLocation(0, next_end, next_strand))
    for part in loc_parts:
        part.strand = strand
    loc_parts.sort(key=lambda x: int(x.start))
    loc_parts = merge_parts(loc_parts)
    if len(loc_parts) == 1:
        return SimpleLocation(loc_parts[0].start, loc_parts[0].end, strand=strand)
    # spans the origin: swap the (sorted) parts unless on the minus strand, so the envelope crosses it
    if strand != -1:
        loc_parts = [loc_parts[1], loc_parts[0]]
    return CompoundLocation(loc_parts)


def inverted(location: AnyLocation, length: int) -> Optional[AnyLocation]:
    """The parts of a sequence of ``length`` that ``location`` doesn't cover, on the strand of
    its first part. None if the location covers the whole sequence."""
    if length <= 0:
        return None
    if isinstance(location, CompoundLocation) and len(location.parts) == 0:
        return SimpleLocation(0, length, strand=0)
    loc_parts = [copy_location(x) for x in _parts(location)]
    strand = loc_parts[0].strand
    for part in loc_parts:
        part.strand = strand
    loc_parts = merge_parts(loc_parts)
    loc_parts.sort(key=lambda x: int(x.start))
    loc_parts = merge_parts(loc_parts)
    covered = sum(len(p) for p in loc_parts)
    if covered == length:
        return None
    if covered == 0:
        return SimpleLocation(0, length, strand=strand)
    out_parts = []
    next_start = 0
    for part in loc_parts:
        if part.start > next_start:
            out_parts.append(SimpleLocation(next_start, part.start, strand=strand))
        next_start = part.end
        if next_start >= length:
            break
    if next_start < length:
        out_parts.append(SimpleLocation(next_start, length, strand=strand))
    if len(out_parts) == 0:
        return None
    if len(out_parts) == 1:
        return out_parts[0]
    operator = location.operator if isinstance(location, CompoundLocation) else "join"
    return CompoundLocation(out_parts, operator)
