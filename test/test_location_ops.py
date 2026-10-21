import os
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation, SimpleLocation
from domainator.location_ops import overlaps, overlay, stranded_end, stranded_end_human_readable, stranded_start, stranded_start_human_readable
import tempfile
from glob import glob
import pytest


def test_overlay_1():
    """
    Test CompoundLocation.overlay on a split CompoundLocation
    """
    loc1 = FeatureLocation(1, 10)
    loc2 = FeatureLocation(20, 30)
    loc3 = FeatureLocation(40, 50)
    loc4 = FeatureLocation(60, 70)
    loc5 = FeatureLocation(80, 90)
    loc6 = FeatureLocation(100, 110)

    cloc1 = CompoundLocation([loc1, loc2, loc3, loc4, loc5, loc6])
    cloc3 = overlay(cloc1, 4, 40)
    assert len(cloc3.parts) == 5
    assert cloc3.parts[0].start == 5
    assert cloc3.parts[0].end == 10
    assert cloc3.parts[1].start == 20
    assert cloc3.parts[1].end == 30
    assert cloc3.parts[2].start == 40
    assert cloc3.parts[2].end == 50
    assert cloc3.parts[3].start == 60
    assert cloc3.parts[3].end == 70
    assert cloc3.parts[4].start == 80
    assert cloc3.parts[4].end == 85


def test_overlaps():
    """
    Test SeqFeature.overlaps method
    """
    loc1 = FeatureLocation(1, 10)
    loc2 = FeatureLocation(10, 40)
    loc3 = FeatureLocation(40, 50)
    loc4 = FeatureLocation(50, 70)
    loc5 = FeatureLocation(70, 90)
    loc6 = FeatureLocation(90, 110)

    cloc1 = CompoundLocation([loc1, loc2, loc3, loc4, loc5, loc6])
    cloc2 = CompoundLocation([loc1, loc2, loc3])
    cloc3 = CompoundLocation([loc4, loc5, loc6])
    cloc4 = CompoundLocation([loc2, loc3, loc4])
    cloc5 = CompoundLocation([loc3, loc4, loc5])
    cloc6 = CompoundLocation([loc2, loc3, loc5])
    cloc7 = CompoundLocation([loc1, loc6])

    assert overlaps(cloc1, cloc2) == True
    assert overlaps(cloc1, cloc3) == True
    assert overlaps(cloc2, cloc3) == False
    assert overlaps(cloc2, cloc4) == True
    assert overlaps(cloc3, cloc4) == True
    assert overlaps(cloc3, cloc5) == True
    assert overlaps(cloc4, cloc5) == True
    assert overlaps(cloc4, cloc6) == True
    assert overlaps(cloc5, cloc6) == True
    assert overlaps(cloc1, cloc7) == True
    assert overlaps(cloc2, cloc7) == True
    assert overlaps(cloc3, cloc7) == True
    assert overlaps(cloc4, cloc7) == False
    assert overlaps(cloc5, cloc7) == False
    assert overlaps(cloc6, cloc7) == False


def _assert_stranded_human_readable_consistency(loc):
    """
    The human readable properties should name the same coordinates as their non-human-readable
    counterparts, adding one only to the coordinate that is a 0-based start.
    """
    if loc.parts[0].strand == -1:
        assert stranded_start_human_readable(loc) == int(stranded_start(loc))
    else:
        assert stranded_start_human_readable(loc) == int(stranded_start(loc)) + 1

    if loc.parts[-1].strand == -1:
        assert stranded_end_human_readable(loc) == int(stranded_end(loc)) + 1
    else:
        assert stranded_end_human_readable(loc) == int(stranded_end(loc))


def test_stranded_coordinates_origin_spanning_forward():
    """
    A forward strand feature crossing the origin of a 2552 bp circular contig: join(2438..2552,1..707)
    The end coordinate should come from the last part, not the first.
    """
    loc = CompoundLocation([SimpleLocation(2437, 2552, 1), SimpleLocation(0, 707, 1)])

    assert stranded_start(loc) == 2437
    assert stranded_end(loc) == 707
    assert stranded_start_human_readable(loc) == 2438
    assert stranded_end_human_readable(loc) == 707
    _assert_stranded_human_readable_consistency(loc)


def test_stranded_coordinates_origin_spanning_reverse():
    """
    A reverse strand feature crossing the origin of the same contig: complement(join(1..707,2438..2552))
    Parts are in transcription order, so the last part holds the end coordinate.
    """
    loc = CompoundLocation([SimpleLocation(0, 707, -1), SimpleLocation(2437, 2552, -1)])

    assert stranded_start(loc) == 707
    assert stranded_end(loc) == 2437
    assert stranded_start_human_readable(loc) == 707
    assert stranded_end_human_readable(loc) == 2438
    _assert_stranded_human_readable_consistency(loc)


def test_stranded_coordinates_non_wrapping_join():
    """
    An intron style join that does not cross the origin, on both strands.
    """
    fwd = CompoundLocation([SimpleLocation(100, 200, 1), SimpleLocation(300, 400, 1)])
    assert stranded_start(fwd) == 100
    assert stranded_end(fwd) == 400
    assert stranded_start_human_readable(fwd) == 101
    assert stranded_end_human_readable(fwd) == 400
    _assert_stranded_human_readable_consistency(fwd)

    rev = CompoundLocation([SimpleLocation(300, 400, -1), SimpleLocation(100, 200, -1)])
    assert stranded_start(rev) == 400
    assert stranded_end(rev) == 100
    assert stranded_start_human_readable(rev) == 400
    assert stranded_end_human_readable(rev) == 101
    _assert_stranded_human_readable_consistency(rev)


def test_stranded_coordinates_single_part():
    """
    Single part locations should agree with the SimpleLocation properties.
    """
    for loc in (SimpleLocation(100, 200, 1), SimpleLocation(100, 200, -1)):
        _assert_stranded_human_readable_consistency(loc)

@pytest.mark.parametrize("location", [FeatureLocation(1000, 1300, 1), FeatureLocation(1000, 1300, -1), FeatureLocation(100, 400, 1)])
def test_simple_location_overlay_bounds(location):
    """overlay bounds are relative to the feature, not to its position on the contig."""
    assert overlay(location, 0, 300) == location
    if location.strand == -1:
        assert overlay(location, 30, 60) == FeatureLocation(1210, 1270, -1)
    else:
        assert overlay(location, 30, 60) == FeatureLocation(int(location.start) + 30, int(location.start) + 90, 1)
    for offset, length in [(150, 300), (0, 301), (-1, 10), (0, -1)]:
        with pytest.raises(ValueError):
            overlay(location, offset, length)
