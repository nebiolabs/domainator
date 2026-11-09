"""The select_by_system command, and the neighborhood code it shares with select_by_cds."""
import os
import tempfile
import warnings

import pytest

from domainator import SYSTEM_FEATURE_NAME, SYSTEM_SEARCH_BEST_HIT_NAME
from domainator import select_by_system as select_by_system_module
from domainator.errors import InputError
from domainator.find_systems import find_systems
from domainator.model import get_hit_best_score
from domainator.select_by_cds import pad_circular_cds, pad_circular_run, up_down_pair
from domainator.select_by_system import select_by_system
from domainator.seqio import parse_seqfiles, write_genbank

from system_records import contig, layout


def systems(record, expr, window_cds, **kwargs):
    return list(find_systems([record], domain_expr=expr, window_cds=window_cds, **kwargs))


def cds_ids(record):
    """The cds_ids of the record's CDSs, in region order."""
    cdss = sorted((f for f in record.features if f.type == "CDS"), key=lambda f: int(f.location.start))
    return [f.qualifiers["cds_id"][0] for f in cdss]


def search_hits(record):
    return [f for f in record.features if f.type == SYSTEM_SEARCH_BEST_HIT_NAME]


def test_pad_circular_run_is_pad_circular_cds_for_one_cds():
    for size in range(1, 7):
        for center in range(size):
            for left in range(8):
                for right in range(8):
                    assert pad_circular_run(center, center, size, left, right) == pad_circular_cds(center, size, left, right)


def test_pad_circular_run():
    assert pad_circular_run(2, 4, 10, 1, 2) == (1, 6)
    assert pad_circular_run(8, 1, 10, 1, 1) == (7, 2)  # wraps the origin
    assert pad_circular_run(8, 1, 10, 9, 9) == (5, 4)  # the 6 CDSs outside are split between the sides
    assert pad_circular_run(8, 1, 10, 9, 0) == (2, 1)
    with pytest.raises(ValueError):
        pad_circular_run(10, 1, 10, 1, 1)


def test_up_down_pair():
    assert up_down_pair(None, None, None, "a", "a_up", "a_down") is None
    assert up_down_pair(3, None, None, "a", "a_up", "a_down") == (3, 3)
    assert up_down_pair(None, None, 2, "a", "a_up", "a_down") == (0, 2)
    with pytest.raises(ValueError, match="Cannot specify a with a_down or a_up"):
        up_down_pair(3, 1, None, "a", "a_up", "a_down")


def test_region_is_the_system():
    record, = systems(contig(layout({"X": 4, "Y": 6})), "X & Y", 5)
    region, = select_by_system([record])
    assert region.id == "contig_3001:5900"
    assert cds_ids(region) == ["cds3", "cds4", "cds5"]
    hit, = search_hits(region)
    assert hit.qualifiers["name"] == ["domain_expr"]
    assert [f.type for f in region.features].count(SYSTEM_FEATURE_NAME) == 1
    assert get_hit_best_score(region) == 100.0


def test_cds_range_counts_from_the_outermost_members():
    record, = systems(contig(layout({"X": 4, "Y": 6})), "X & Y", 5)
    region, = select_by_system([record], cds_range=(2, 1))
    assert cds_ids(region) == ["cds1", "cds2", "cds3", "cds4", "cds5", "cds6"]


def test_minus_strand_system_pads_in_system_orientation():
    record, = systems(contig(layout({"X": 4, "Y": 6}, strand="-")), "X & Y", 5)
    region, = select_by_system([record], cds_range=(2, 0))
    assert region.id == "contig_7900:3001rc"
    assert sorted(cds_ids(region)) == ["cds3", "cds4", "cds5", "cds6", "cds7"]
    kept, = select_by_system([record], cds_range=(2, 0), normalize_direction=False)
    assert kept.id == "contig_3001:7900"


def test_kb_range():
    record, = systems(contig(layout({"X": 4, "Y": 6})), "X & Y", 5)
    region, = select_by_system([record], kb_range=(0.5, 1.0))
    assert region.id == "contig_2501:6900"


def test_whole_contig():
    record, = systems(contig(layout({"X": 4, "Y": 6})), "X & Y", 5)
    region, = select_by_system([record], whole_contig=True)
    assert len(region) == len(record)


def test_min_context_skips_systems_at_a_contig_end():
    record, = systems(contig(layout({"X": 1, "Y": 2}, length=6)), "X & Y", 5)
    assert list(select_by_system([record], min_context_cds=(1, 1))) == []
    assert len(list(select_by_system([record], min_context_cds=(0, 4)))) == 1
    assert list(select_by_system([record], min_context_cds=(0, 5))) == []
    assert len(list(select_by_system([record], min_context_kb=(0, 4.1)))) == 1  # 1900..6000 is 4.1 kb
    assert list(select_by_system([record], min_context_kb=(0, 4.2))) == []
    # Both must hold.
    assert list(select_by_system([record], min_context_cds=(0, 1), min_context_kb=(0.1, 0))) == []


def test_min_context_upstream_is_in_system_orientation():
    record, = systems(contig(layout({"X": 1, "Y": 2}, length=6, strand="-")), "X & Y", 5)
    # On the minus strand, the system's upstream side is toward the contig's right end.
    assert len(list(select_by_system([record], min_context_cds=(4, 0)))) == 1
    assert list(select_by_system([record], min_context_cds=(0, 1))) == []


def test_min_context_does_not_count_nucleotide_pseudo_cdss():
    record, = systems(contig(["X", "Y", "nc:R", ""]), "X & Y", 5)
    assert len(list(select_by_system([record], min_context_cds=(0, 1)))) == 1
    assert list(select_by_system([record], min_context_cds=(0, 2))) == []


def test_circular_system_across_the_origin():
    record, = systems(contig(layout({"X": 12, "Y": 1}), topology="circular"), "X & Y", 3)
    feature, = [f for f in record.features if f.type == SYSTEM_FEATURE_NAME]
    assert len(feature.location.parts) == 2
    region, = select_by_system([record], cds_range=(1, 1))
    assert region.id == "contig_10001:1900"
    assert cds_ids(region) == ["cds10", "cds11", "cds0", "cds1"]
    # 10 CDSs lie outside the system.
    assert len(list(select_by_system([record], min_context_cds=(5, 5)))) == 1
    assert list(select_by_system([record], min_context_cds=(5, 6))) == []


def test_model_neighborhood():
    record, = systems(contig(layout({"X": 4, "Y": 6})), "X & Y", 5)
    feature, = [f for f in record.features if f.type == SYSTEM_FEATURE_NAME]
    feature.qualifiers["neighborhood_cds"] = ["1"]
    region, = select_by_system([record], model_neighborhood=True)
    assert cds_ids(region) == ["cds2", "cds3", "cds4", "cds5", "cds6"]
    feature.qualifiers["neighborhood_kb"] = ["0.5"]
    region, = select_by_system([record], model_neighborhood=True)
    assert region.id == "contig_2501:6400"
    with pytest.raises(InputError):
        list(select_by_system([record], model_neighborhood=True, cds_range=(1, 1)))


def _two_systems():
    record = contig(layout({"X": 2, "Y": 3, "P": 8, "Q": 9}))
    record, = systems(record, "X & Y", 3, system_name="XY")
    record, = list(find_systems([record], domain_expr="P & Q", window_cds=3, system_name="PQ", clear=False))
    return record


@pytest.mark.parametrize("kwargs, names", [
    ({}, ["XY", "PQ"]),
    ({"systems": {"PQ"}}, ["PQ"]),
    ({"systems": {"command_line:XY"}}, ["XY"]),
    ({"systems": {"PQ"}, "invert": True}, ["XY"]),
    ({"system_categories": {"defense"}}, []),
    ({"databases": {"other"}}, []),
    ({"databases": {"other"}, "invert": True, "systems": {"PQ"}}, []),
    ({"min_score": 101}, []),
    ({"strand": "r"}, []),
    ({"partial_systems": "only"}, []),
])
def test_selectors(kwargs, names):
    regions = list(select_by_system([_two_systems()], **kwargs))
    assert [search_hits(region)[0].qualifiers["name"][0] for region in regions] == names


def test_each_region_has_one_search_hit():
    record, = systems(contig(layout({"X": 2, "Y": 3})), "X & Y", 3, system_name="XY")
    record, = list(find_systems([record], domain_expr="X & Y", window_cds=4, system_name="XY_wide", clear=False))
    regions = list(select_by_system([record]))
    assert len(regions) == 2
    for region in regions:
        assert len(search_hits(region)) == 1
    # Selecting again from System_Search features keeps the one search hit.
    again = list(select_by_system(regions[:1], system_type="system_search", cds_range=(1, 1)))
    assert len(again) == 1 and len(search_hits(again[0])) == 1


def test_partial_systems():
    record, = systems(contig(layout({"X": 1}, length=4)), "X", 3)
    feature, = [f for f in record.features if f.type == SYSTEM_FEATURE_NAME]
    feature.qualifiers["system_partial"] = ["<"]
    assert list(select_by_system([record], partial_systems="exclude")) == []
    assert len(list(select_by_system([record], partial_systems="only"))) == 1


def test_nucleotide_members_are_found_after_slicing():
    record, = systems(contig(["", "X", "nc:R", "", ""]), "X & R", 3)
    region, = select_by_system([record], cds_range=(1, 1))
    assert cds_ids(region) == ["cds0", "cds1", "cds3"]
    # The nucleotide member's id holds its old coordinates; it is found by name instead.
    again, = select_by_system([region], system_type="system_search", cds_range=(0, 0))
    assert cds_ids(again) == ["cds1"]
    assert len(again) == 1200  # 1000..2200: X's CDS through the nucleotide hit


def test_missing_members_warn():
    record, = systems(contig(layout({"X": 2, "Y": 3})), "X & Y", 3)
    record.features = [f for f in record.features if f.qualifiers.get("cds_id") not in (["cds1"], ["cds2"])]
    with pytest.warns(UserWarning, match="none of its member CDSs"):
        assert list(select_by_system([record])) == []


def test_max_region_overlap():
    record, = systems(contig(layout({"X": 2, "Y": 3})), "X & Y", 3, system_name="XY")
    record, = list(find_systems([record], domain_expr="X & Y", window_cds=4, system_name="XY_wide", clear=False))
    assert len(list(select_by_system([record], cds_range=(1, 1), max_region_overlap=0.5))) == 1


def test_main(tmp_path):
    record, = systems(contig(layout({"X": 4, "Y": 6})), "X & Y", 5)
    write_genbank([record], str(tmp_path / "in.gb"))
    out = str(tmp_path / "out.gb")
    select_by_system_module.main(["-i", str(tmp_path / "in.gb"), "-o", out, "--cds_range", "1",
                                  "--min_context_cds_down", "1", "--pad"])
    region, = parse_seqfiles([out])
    assert cds_ids(region) == ["cds2", "cds3", "cds4", "cds5", "cds6"]
    assert len(search_hits(region)) == 1


@pytest.mark.parametrize("args, message", [
    (["--min_context_cds", "1", "--min_context_cds_up", "1"], "Cannot specify --min_context_cds"),
    (["--cds_range", "1", "--model_neighborhood"], "--model_neighborhood can't be combined"),
    (["--cds_range", "1", "--kb_range", "1"], "Cannot specify kb_range and cds_range"),
])
def test_main_rejects_conflicting_options(tmp_path, capsys, args, message):
    with pytest.raises(SystemExit):
        select_by_system_module.main(["-i", "in.gb", "-o", str(tmp_path / "out.gb"), *args])
    assert message in capsys.readouterr().err
