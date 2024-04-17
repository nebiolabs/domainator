import os
from domainator.Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
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