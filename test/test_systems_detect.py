"""System detection (domainator.systems.detect) on synthetic contigs, without HMMER."""
import os
import tempfile
import warnings

import pytest
from Bio.SeqFeature import BeforePosition, CompoundLocation, FeatureLocation

from domainator import SYSTEM_FEATURE_NAME
from domainator.seqio import parse_seqfiles, write_genbank
from domainator.systems.detect import detect_systems
from domainator.systems.features import component_values, system_feature
from domainator.systems.native import collection_from_data, window_model

from system_records import contig, hit_feature, layout, member_ids


def models(*systems, **collection):
    """A collection from model entries (dicts, as in a native file), with a 1-CDS chain window by default."""
    data = {"domainator_systems": 1, "database": collection.pop("database", "test"),
            "defaults": {"window": {"mode": "chain", "max_cds": 1}}, "systems": list(systems), **collection}
    return collection_from_data(data)


def found(record, collection, **options):
    return [member_ids(hit) for hit in detect_systems(record, [collection], **options)]


# The worked example of definition_formats.md: A, B, C and D at CDS positions 1, 4, 7 and 10.
@pytest.mark.parametrize("mode,max_cds,is_found", [
    ("chain", 2, True),
    ("chain", 1, False),
    ("anchor", 2, False),
    ("anchor", 5, True),
    ("anchor", 4, False),
    ("span", 10, True),
    ("span", 9, False),
])
def test_worked_example(mode, max_cds, is_found):
    record = contig(layout({"A": 1, "B": 4, "C": 7, "D": 10}))
    hits = found(record, window_model("A & B & C & D", window_cds=max_cds, window_mode=mode))
    assert hits == ([["cds0", "cds3", "cds6", "cds9"]] if is_found else [])


@pytest.mark.parametrize("mode", ["chain", "anchor", "span"])
def test_kb_limits_at_the_boundary(mode):
    record = contig(layout({"A": 1, "B": 3}))  # A ends at 900, B starts at 2000: a 1100 bp gap, a 2900 bp span
    limit_kb = 2.9 if mode == "span" else 1.1
    assert found(record, window_model("A & B", window_kb=limit_kb, window_mode=mode)) == [["cds0", "cds2"]]
    assert found(record, window_model("A & B", window_kb=limit_kb - 0.001, window_mode=mode)) == []


def test_span_negation_reaches_beyond_the_system():
    record = contig(layout({"X": 1, "Y": 2, "W": 9}))
    assert found(record, window_model("X & Y & ~W", window_cds=10)) == []  # W is within reach
    assert found(record, window_model("X & Y & ~W", window_cds=5)) == [["cds0", "cds1"]]  # W is out of reach


def test_span_system_is_the_tightest_run():
    record = contig(["A", "", "B", "", "", "", "", "", "A"])
    assert found(record, window_model("A & B", window_cds=10)) == [["cds0", "cds2"]]


def test_span_picks_the_best_of_overlapping_runs():
    record = contig(["A:60", "B", "A:90"])
    assert found(record, window_model("A & B", window_cds=2)) == [["cds1", "cds2"]]


def test_tandem_copies_are_separate_systems():
    record = contig(["A", "B", "", "", "", "A", "B"])
    for mode, size in (("chain", 1), ("span", 2), ("anchor", 1)):
        assert found(record, window_model("A & B", window_cds=size, window_mode=mode)) == [["cds0", "cds1"], ["cds5", "cds6"]], mode


def test_nucleotide_hits_dont_count_as_cds():
    record = contig(["A", "nc:R", "B"])
    assert found(record, window_model("A & B", window_cds=0, window_mode="chain")) == [["cds0", "cds2"]]
    hit, = detect_systems(record, [window_model("A & R & B", window_cds=0, window_mode="chain")])
    assert member_ids(hit) == ["cds0", "nuc_1001_1_1200", "cds2"]


def test_attach_to_cds():
    record = contig(["", "B"])
    record.features.append(hit_feature(FeatureLocation(100, 400, 1), "R", ".", program="nhmmer"))
    for attach, expected in ((False, ["nuc_101_1_400", "cds1"]), (True, ["cds0", "cds1"])):
        collection = models({"name": "S", "window": {"mode": "chain", "max_cds": 0},
                             "components": [{"name": "R", "domains": [{"name": "R", "attach_to_cds": attach}]},
                                            {"name": "B", "domains": ["B"]}]})
        hits = detect_systems(record, [collection])
        assert [member_ids(hit) for hit in hits] == [expected]
        assert [hit.location for hit in hits] == [FeatureLocation(0 if attach else 100, 1900, 1)]


def test_no_cds_record():
    record = contig(["nc:R", "nc:S"])
    record.features = [f for f in record.features if f.type != "CDS"]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        assert found(record, window_model("R & S", window_cds=3)) == []
        assert found(record, window_model("R & S", window_kb=1)) == [["nuc_1_1_200", "nuc_301_1_500"]]
    assert sum("only in CDSs" in str(w.message) for w in caught) == 1


def test_forbidden_inside_a_cluster():
    collection = models({"name": "S", "components": [{"name": "A", "domains": ["A"]}, {"name": "B", "domains": ["B"]},
                                                      {"name": "T", "role": "forbidden", "domains": ["T"]}]})
    assert found(contig(["A", "B"]), collection) == [["cds0", "cds1"]]
    assert found(contig(["A", "T", "B"]), collection) == []
    assert found(contig(["A", "", "T"]), collection) == []  # chains in, and vetoes
    assert found(contig(["A", "B", "", "", "T"]), collection) == [["cds0", "cds1"]]


def test_quorum_and_neutral_trimming():
    collection = models({"name": "S", "min_core": 1, "min_total": 2,
                         "components": [{"name": "A", "domains": ["A"]}, {"name": "B", "role": "accessory", "domains": ["B"]},
                                        {"name": "C", "role": "accessory", "domains": ["C"]},
                                        {"name": "N", "role": "neutral", "domains": ["N"]}]})
    hit, = detect_systems(contig(["N", "A", "N", "C", "N"]), [collection])
    assert member_ids(hit) == ["cds1", "cds2", "cds3"]
    assert [v.role for v in component_values(system_feature(hit))] == ["core", "neutral", "accessory"]
    assert (hit.core_count, hit.total_count, hit.completeness, hit.member_count) == (1, 2, 2 / 3, 3)
    assert found(contig(["B", "C"]), collection) == []  # needs a core component
    assert found(contig(["A", "N"]), collection) == []  # neutral components don't count


def test_exchangeable_references():
    collection = models({"name": "Gabija", "components": [{"name": "GajA", "domains": ["GajA"]},
                                                           {"name": "GajB", "domains": ["GajB_1", {"name": "GajB_2", "exchangeable": True}]}]})
    assert found(contig(["GajA", "GajB_2"]), collection) == [["cds0", "cds1"]]
    assert found(contig(["GajA", "GajB_1"]), collection) == [["cds0", "cds1"]]


def test_one_cds_matching_two_components_takes_the_primary_one():
    collection = models({"name": "S", "min_core": 2,
                         "components": [{"name": "A", "domains": ["H"]}, {"name": "B", "domains": ["H", "B"]}]})
    assert found(contig(["H", ""]), collection) == []  # one CDS can't be both core components
    hit, = detect_systems(contig(["H:60", "B:40"]), [collection])
    assert [v.component for v in component_values(system_feature(hit))] == ["A", "B"]


def test_fusion_lists_each_component():
    hit, = detect_systems(contig(["X+Y", "", "Z"]), [window_model("X & Y & Z", window_cds=3)])
    values = component_values(system_feature(hit))
    assert [(v.cds_id, v.component) for v in values] == [("cds0", "X"), ("cds0", "Y"), ("cds2", "Z")]
    assert hit.member_count == 2
    assert hit.score == 150


def test_thresholds():
    record = contig(["A:30", "B"])
    record.features[1].qualifiers["rend"] = ["30"]  # A's hit covers 30% of its reference

    def rule(**thresholds):
        return models({"name": "S", "components": [{"name": "A", "domains": [{"name": "A", **thresholds}]},
                                                   {"name": "B", "domains": ["B"]}]})
    assert found(record, rule(min_score=30)) == [["cds0", "cds1"]]
    assert found(record, rule(min_score=31)) == []
    assert found(record, rule(min_reference_coverage=0.3)) == [["cds0", "cds1"]]
    assert found(record, rule(min_reference_coverage=0.31)) == []
    assert found(record, rule(evalue=1e-11)) == []
    assert found(record, rule(database="other")) == []
    assert found(record, rule(min_target_length=300)) == [["cds0", "cds1"]]  # 900 bp: 300 aa
    assert found(record, rule(max_target_length=299)) == []
    assert found(record, rule(min_target_coverage=1.0)) == [["cds0", "cds1"]]
    assert found(record, window_model("A & B", window_cds=2), evalue=1e-11) == []
    assert found(record, window_model("A & B", window_cds=2), databases={"other"}) == []


@pytest.mark.parametrize("threshold, qualifier", [("min_tmscore", "tmscore"), ("min_lddt", "lddt"), ("min_prob", "prob")])
def test_structural_thresholds(threshold, qualifier):
    """Read from foldseek hits' qualifiers; a hit without one (an HMM hit) fails them."""
    record = contig(["A", "B"])
    record.features[1].qualifiers[qualifier] = ["0.650"]

    def rule(value):
        return models({"name": "S", "components": [{"name": "A", "domains": [{"name": "A", threshold: value}]},
                                                   {"name": "B", "domains": ["B"]}]})
    assert found(record, rule(0.65)) == [["cds0", "cds1"]]
    assert found(record, rule(0.66)) == []
    del record.features[1].qualifiers[qualifier]
    assert found(record, rule(0.1)) == []


def test_partial_cdss():
    record = contig(["A", "B"])
    record.features[0].location = FeatureLocation(BeforePosition(0), 900, 1)
    model = window_model("A & B", window_cds=2)
    assert found(record, model) == [["cds0", "cds1"]]
    assert found(record, model, partial="exclude") == []


def test_competing_references():
    record = contig(["G1:80+G2:70", "B"])
    collection = models({"name": "S", "condition": "G2 & B"}, competing=[["G1", "G2"]])
    assert found(record, collection) == []
    assert found(record, models({"name": "S", "condition": "G2 & B"})) == [["cds0", "cds1"]]
    tie = contig(["G1:80+G2:80", "B"])
    assert found(tie, collection) == [["cds0", "cds1"]]


def test_superiors():
    collection = models({"name": "Big", "condition": "A & B & C"}, {"name": "Small", "condition": "A & B", "superiors": ["Big"]})
    hits = detect_systems(contig(["A", "B", "C", "", "", "A", "B"]), [collection])
    assert [(hit.model.name, member_ids(hit)) for hit in hits] == [("Big", ["cds0", "cds1", "cds2"]), ("Small", ["cds5", "cds6"])]


def test_greedy_overlap():
    collection = models({"name": "One", "condition": "A & B"}, {"name": "Two", "condition": "B & C"})
    record = contig(["A:90", "B", "C"])
    assert [hit.model.name for hit in detect_systems(record, [collection])] == ["One", "Two"]
    assert [hit.model.name for hit in detect_systems(record, [collection], max_system_overlap=0)] == ["One"]
    assert [hit.model.name for hit in detect_systems(record, [collection], max_systems_per_contig=1)] == ["One"]


def test_multi_model_members_are_shared():
    collection = models({"name": "One", "components": [{"name": "A", "domains": ["A"]}, {"name": "B", "domains": ["B"], "multi_model": True}]},
                        {"name": "Two", "components": [{"name": "B", "domains": ["B"], "multi_model": True}, {"name": "C", "domains": ["C"]}]})
    assert [hit.model.name for hit in detect_systems(contig(["A:90", "B", "C"]), [collection], max_system_overlap=0.5)] == ["One", "Two"]


def test_anchor_must_contribute():
    record = contig(["A", "", "C"])
    hit, = detect_systems(record, [window_model("(A & B) | C", window_cds=5, window_mode="anchor")])
    assert member_ids(hit) == ["cds2"]


@pytest.mark.parametrize("slots,strand", [
    (["-A", "-B", "-C"], -1),
    (["A", "-B", "C"], 1),
    (["-A", "B", "-C"], -1),
    (["-A:90", "B"], -1),  # a tie: the best-scoring core member decides
])
def test_orientation(slots, strand):
    model = window_model(" & ".join(s.lstrip("-").split(":")[0] for s in slots), window_cds=3)
    hit, = detect_systems(contig(slots), [model])
    assert hit.location.strand == strand


def test_orientation_ignores_nucleotide_members_unless_asked():
    record = contig(["A", "-nc:R", "B"])
    hit, = detect_systems(record, [window_model("A & R & B", window_cds=2)])
    assert hit.location.strand == 1
    collection = models({"name": "S", "condition": "A & R & B", "orientation_from": "R"})
    hit, = detect_systems(record, [collection])
    assert hit.location.strand == -1


def test_order_follows_system_orientation():
    plus = contig(["HicB", "RelE"])
    minus = contig(["-RelE", "-HicB"])
    model = window_model("HicB < RelE", window_cds=2)
    assert found(plus, model) == [["cds0", "cds1"]]
    assert found(minus, model) == [["cds0", "cds1"]]
    assert found(contig(["RelE", "HicB"]), model) == []
    assert found(plus, window_model("RelE & HicB & ~HicB < RelE", window_cds=2)) == []


def test_same_strand():
    record = contig(["A", "-B", "A", "B"])
    assert found(record, window_model("A & B", window_cds=2)) == [["cds0", "cds1"], ["cds2", "cds3"]]
    assert found(record, window_model("A & B", window_cds=2, same_strand=True)) == [["cds2", "cds3"]]


@pytest.mark.parametrize("mode", ["chain", "anchor", "span"])
def test_circular_wrap(mode):
    record = contig(["B", "", "", "", "", "A"], topology="circular")
    hits = detect_systems(record, [window_model("A & B", window_cds=1 if mode != "span" else 2, window_mode=mode)])
    assert [member_ids(hit) for hit in hits] == [["cds5", "cds0"]]
    location = hits[0].location
    assert isinstance(location, CompoundLocation)
    assert [(int(p.start), int(p.end)) for p in location.parts] == [(5000, len(record)), (0, 900)]
    assert found(contig(["B", "", "", "", "", "A"]), window_model("A & B", window_cds=1, window_mode=mode)) == []


def test_circular_minus_strand_wrap_round_trip():
    record = contig(["-B", "", "", "", "", "-A"], topology="circular")
    hit, = detect_systems(record, [window_model("A & B", window_cds=1, window_mode="chain")])
    record.features.append(system_feature(hit))
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "out.gb")
        write_genbank([record], path)
        with open(path) as handle:
            assert "Domain_System   complement(join(5001..6000,1..900))" in handle.read()
        again, = parse_seqfiles([path])
    feature, = [f for f in again.features if f.type == SYSTEM_FEATURE_NAME]
    assert feature.location == hit.location


def test_system_partial():
    collection = models({"name": "S", "min_core": 2, "min_total": 2, "window": {"mode": "chain", "max_cds": 2},
                         "components": [{"name": "A", "domains": ["A"]}, {"name": "B", "domains": ["B"]},
                                        {"name": "C", "role": "accessory", "domains": ["C"]}]})
    assert [hit.partial for hit in detect_systems(contig(["", "", "A", "B"]), [collection])] == ["<>"]
    assert [hit.partial for hit in detect_systems(contig(["", "", "", "A", "B", "", "", ""]), [collection])] == [""]
    assert [hit.partial for hit in detect_systems(contig(["", "", "", "A", "B", ""]), [collection])] == [">"]
    assert [hit.partial for hit in detect_systems(contig(["", "", "", "-A", "-B", ""]), [collection])] == ["<"]
    assert [hit.partial for hit in detect_systems(contig(["A", "B", "C"]), [collection])] == [""]  # nothing missing
    assert [hit.partial for hit in detect_systems(contig(["", "", "A", "B"], topology="circular"), [collection])] == [""]


def test_max_span_caps_split_a_chain():
    collection = models({"name": "S", "condition": "A & B", "window": {"mode": "chain", "max_cds": 1, "max_span_cds": 2}})
    assert found(contig(["A", "B", "A", "B"]), collection) == [["cds0", "cds1"], ["cds2", "cds3"]]


def test_per_component_gap_override():
    collection = models({"name": "S", "window": {"mode": "chain", "max_cds": 3},
                         "components": [{"name": "A", "domains": ["A"], "max_gap_cds": 1}, {"name": "B", "domains": ["B"]}]})
    assert found(contig(["A", "", "B"]), collection) == [["cds0", "cds2"]]
    assert found(contig(["A", "", "", "B"]), collection) == []


def test_protein_records_have_no_systems():
    record = contig(["A", "B"])
    record.annotations["molecule_type"] = "protein"
    assert found(record, window_model("A & B", window_cds=2)) == []


def test_feature_qualifiers():
    record = contig(["A", "B"])
    record.features[1].qualifiers["name"] = ["A|x"]
    collection = models({"name": "S", "category": "defense", "neighborhood_kb": 5,
                         "components": [{"name": "A", "domains": ["A|x"]}, {"name": "B", "domains": ["B"]}]})
    hit, = detect_systems(record, [collection])
    feature = system_feature(hit)
    assert feature.type == SYSTEM_FEATURE_NAME
    assert {key: values for key, values in feature.qualifiers.items() if key != "component"} == {
        "name": ["S"], "database": ["test"], "program": ["domainator"], "category": ["defense"], "description": ["."],
        "accession": ["."], "window": ["chain:cds=1"], "score": ["100.0"], "evalue": ["1.0e-10"], "completeness": ["1.000"],
        "core_count": ["2"], "total_count": ["2"], "member_count": ["2"], "neighborhood_kb": ["5"]}
    assert feature.qualifiers["component"] == ["core|A|cds0|db|A|x", "core|B|cds1|db|B"]
    assert component_values(feature)[0].reference == "A|x"
