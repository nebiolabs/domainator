import os
from domainator.Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation, SimpleLocation
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
    cloc3 = cloc1.overlay(4, 40)
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

    assert cloc1.overlaps(cloc2) == True
    assert cloc1.overlaps(cloc3) == True
    assert cloc2.overlaps(cloc3) == False
    assert cloc2.overlaps(cloc4) == True
    assert cloc3.overlaps(cloc4) == True
    assert cloc3.overlaps(cloc5) == True
    assert cloc4.overlaps(cloc5) == True
    assert cloc4.overlaps(cloc6) == True
    assert cloc5.overlaps(cloc6) == True
    assert cloc1.overlaps(cloc7) == True
    assert cloc2.overlaps(cloc7) == True
    assert cloc3.overlaps(cloc7) == True
    assert cloc4.overlaps(cloc7) == False
    assert cloc5.overlaps(cloc7) == False
    assert cloc6.overlaps(cloc7) == False


def _assert_stranded_human_readable_consistency(loc):
    """
    The human readable properties should name the same coordinates as their non-human-readable
    counterparts, adding one only to the coordinate that is a 0-based start.
    """
    if loc.parts[0].strand == -1:
        assert loc.stranded_start_human_readable == int(loc.stranded_start)
    else:
        assert loc.stranded_start_human_readable == int(loc.stranded_start) + 1

    if loc.parts[-1].strand == -1:
        assert loc.stranded_end_human_readable == int(loc.stranded_end) + 1
    else:
        assert loc.stranded_end_human_readable == int(loc.stranded_end)


def test_stranded_coordinates_origin_spanning_forward():
    """
    A forward strand feature crossing the origin of a 2552 bp circular contig: join(2438..2552,1..707)
    The end coordinate should come from the last part, not the first.
    """
    loc = CompoundLocation([SimpleLocation(2437, 2552, 1), SimpleLocation(0, 707, 1)])

    assert loc.stranded_start == 2437
    assert loc.stranded_end == 707
    assert loc.stranded_start_human_readable == 2438
    assert loc.stranded_end_human_readable == 707
    _assert_stranded_human_readable_consistency(loc)


def test_stranded_coordinates_origin_spanning_reverse():
    """
    A reverse strand feature crossing the origin of the same contig: complement(join(1..707,2438..2552))
    Parts are in transcription order, so the last part holds the end coordinate.
    """
    loc = CompoundLocation([SimpleLocation(0, 707, -1), SimpleLocation(2437, 2552, -1)])

    assert loc.stranded_start == 707
    assert loc.stranded_end == 2437
    assert loc.stranded_start_human_readable == 707
    assert loc.stranded_end_human_readable == 2438
    _assert_stranded_human_readable_consistency(loc)


def test_stranded_coordinates_non_wrapping_join():
    """
    An intron style join that does not cross the origin, on both strands.
    """
    fwd = CompoundLocation([SimpleLocation(100, 200, 1), SimpleLocation(300, 400, 1)])
    assert fwd.stranded_start == 100
    assert fwd.stranded_end == 400
    assert fwd.stranded_start_human_readable == 101
    assert fwd.stranded_end_human_readable == 400
    _assert_stranded_human_readable_consistency(fwd)

    rev = CompoundLocation([SimpleLocation(300, 400, -1), SimpleLocation(100, 200, -1)])
    assert rev.stranded_start == 400
    assert rev.stranded_end == 100
    assert rev.stranded_start_human_readable == 400
    assert rev.stranded_end_human_readable == 101
    _assert_stranded_human_readable_consistency(rev)


def test_stranded_coordinates_single_part():
    """
    Single part locations should agree with the SimpleLocation properties.
    """
    for loc in (SimpleLocation(100, 200, 1), SimpleLocation(100, 200, -1)):
        _assert_stranded_human_readable_consistency(loc)