"""Every tool that extracts regions of contigs names them the same way (see utils.region_id).

The name of each extracted region is decoded back into coordinates on its parent contig, and the sequence at those
coordinates is compared to the extracted sequence. Regions extracted by different tools that have the same sequence
must also have the same name.
"""
import os
import re
import tempfile
from collections import defaultdict

import pytest

from domainator import domain_search, extract_domains, extract_unannotated, rotate_contigs, select_by_cds, select_by_coord, trim_contigs
from domainator.Bio.Seq import Seq
from domainator.Bio.SeqFeature import CompoundLocation, FeatureLocation
from domainator.Bio.SeqRecord import SeqRecord
from domainator.seqio import parse_seqfiles
from domainator.model import region_id
from domainator.locations import slice_record_from_location

LINEAR = "pDONR201_multi_genemark_domainator.gb" # CAT, CcdA, CcdB on the reverse strand, APH on the forward strand
CIRCULAR = "pDONR201_domainator_circular.gb" # circular, 4470 bp, same domains as LINEAR
SPR = "bacillus_phage_SPR_with_annotations.gb" # circular, 2552 bp, all features on the reverse strand, Thymidylat_synt spans the origin
SPR_REVERSED = "bacillus_phage_SPR_with_annotations_reversed.gb" # the same contig, reverse complemented
SPR_UNANNOTATED = "bacillus_phage_SPR.gb"
PROTEINS = "simple_genpept.gb"

_REGION_SUFFIX = re.compile(r"_(\d+):(\d+)(rc)?")

# (case id, tool, input file, extra arguments)
CASES = [
    ("select_by_cds_cds_range", select_by_cds, LINEAR, ["--domains", "CAT", "--cds_range", "1"]),
    ("select_by_cds_cds_range_keep_direction", select_by_cds, LINEAR, ["--domains", "CAT", "--cds_range", "1", "--keep_direction"]),
    ("select_by_cds_kb_range", select_by_cds, LINEAR, ["--domains", "CAT", "APH", "--kb_range", "1"]),
    ("select_by_cds_kb_range_keep_direction", select_by_cds, LINEAR, ["--domains", "CAT", "APH", "--kb_range", "1", "--keep_direction"]),
    ("select_by_cds_whole_contig", select_by_cds, LINEAR, ["--domains", "CAT", "APH", "--whole_contig", "--contigs", "pDONR201_1"]),
    ("select_by_cds_whole_contig_keep_direction", select_by_cds, LINEAR, ["--domains", "CAT", "--whole_contig", "--keep_direction", "--contigs", "pDONR201_1"]),
    ("select_by_cds_circular_kb_range", select_by_cds, CIRCULAR, ["--domains", "CAT", "--kb_range", "3"]),
    ("select_by_cds_circular_kb_range_keep_direction", select_by_cds, CIRCULAR, ["--domains", "CAT", "--kb_range", "3", "--keep_direction"]),
    ("select_by_cds_origin_cds", select_by_cds, SPR, ["--domains", "Thymidylat_synt", "--kb_range", "0.1"]),
    ("select_by_cds_origin_cds_keep_direction", select_by_cds, SPR, ["--domains", "Thymidylat_synt", "--kb_range", "0.1", "--keep_direction"]),
    ("select_by_cds_circular_whole_contig", select_by_cds, SPR, ["--domains", "DHFR_1", "--whole_contig"]),
    ("select_by_cds_circular_whole_contig_keep_direction", select_by_cds, SPR, ["--domains", "DHFR_1", "--whole_contig", "--keep_direction"]),
    ("domain_search_whole_contig", domain_search, SPR_UNANNOTATED, ["-r", "SPR.hmm", "--whole_contig", "-Z", "1000", "--cpu", "1"]),
    ("domain_search_whole_contig_keep_direction", domain_search, SPR_UNANNOTATED, ["-r", "SPR.hmm", "--whole_contig", "-Z", "1000", "--cpu", "1", "--keep_direction"]),
    ("domain_search_cds_range", domain_search, SPR_UNANNOTATED, ["-r", "SPR.hmm", "--cds_range", "1", "-Z", "1000", "--cpu", "1"]),
    ("domain_search_cds_range_keep_direction", domain_search, SPR_UNANNOTATED, ["-r", "SPR.hmm", "--cds_range", "1", "-Z", "1000", "--cpu", "1", "--keep_direction"]),
    ("extract_domains", extract_domains, LINEAR, ["--domains", "CAT", "APH"]),
    ("extract_domains_keep_direction", extract_domains, LINEAR, ["--domains", "CAT", "APH", "--keep_direction"]),
    ("extract_domains_pad", extract_domains, LINEAR, ["--domains", "CAT", "APH", "--extend_up", "100", "--extend_down", "50"]),
    ("extract_domains_pad_keep_direction", extract_domains, LINEAR, ["--domains", "CAT", "APH", "--extend_up", "100", "--extend_down", "50", "--keep_direction"]),
    ("extract_domains_combine", extract_domains, LINEAR, ["--domains", "CcdA", "CcdB", "--combine"]),
    ("extract_domains_combine_keep_direction", extract_domains, LINEAR, ["--domains", "CcdA", "CcdB", "--combine", "--keep_direction"]),
    ("extract_domains_origin", extract_domains, SPR, ["--domains", "Thymidylat_synt", "--extend_up", "10"]),
    ("extract_domains_origin_keep_direction", extract_domains, SPR, ["--domains", "Thymidylat_synt", "--extend_up", "10", "--keep_direction"]),
    ("extract_domains_origin_forward", extract_domains, SPR_REVERSED, ["--domains", "Thymidylat_synt"]),
    ("extract_unannotated", extract_unannotated, PROTEINS, []),
    ("extract_unannotated_largest", extract_unannotated, PROTEINS, ["--largest"]),
    ("trim_contigs", trim_contigs, LINEAR, ["--cds_range", "2"]),
    ("rotate_contigs", rotate_contigs, SPR, ["--domains", "DHFR_1", "--anchor", "center"]),
    ("rotate_contigs_keep_direction", rotate_contigs, SPR, ["--domains", "DHFR_1", "--anchor", "center", "--keep_direction"]),
    ("rotate_contigs_forward", rotate_contigs, SPR_REVERSED, ["--domains", "Thymidylat_synt"]),
]

# (case id, input file, rows of contig, start, end, strand)
COORD_CASES = [
    ("select_by_coord_linear", LINEAR, [("pDONR201_1", 1639, 2244, "+"), ("pDONR201_1", 1639, 2244, "-"), ("pDONR201_1", 1, 4470, "+"), ("pDONR201_1", 1, 4470, "-")]),
    ("select_by_coord_circular", SPR, [("OM236515mod", 1846, 115, "+"), ("OM236515mod", 1846, 115, "-"), ("OM236515mod", 1341, 1340, "+"), ("OM236515mod", 1341, 1340, "-")]),
]

# splice output skips the sequence between domains, so only the first and last bases are named
SPLICE_CASES = [
    ("extract_domains_splice", extract_domains, LINEAR, ["--domains", "CcdA", "CcdB", "--splice"]),
    ("extract_domains_splice_keep_direction", extract_domains, LINEAR, ["--domains", "CcdA", "CcdB", "--splice", "--keep_direction"]),
]


def _write_table(path, rows):
    with open(path, "w") as handle:
        print("contig\tstart\tend\tstrand", file=handle)
        for row in rows:
            print("\t".join(str(value) for value in row), file=handle)


def _run(datadir, output_dir, case_id, tool, input_file, args):
    out = os.path.join(output_dir, case_id + ".gb")
    args = [str(datadir / arg) if arg.endswith(".hmm") else arg for arg in args]
    tool.main(["-i", str(datadir / input_file), "-o", out] + args)
    return list(parse_seqfiles([out]))


def _run_coords(datadir, output_dir, case_id, input_file, rows):
    table = os.path.join(output_dir, case_id + ".tsv")
    _write_table(table, rows)
    return _run(datadir, output_dir, case_id, select_by_coord, input_file, ["--coords", table])


def _parents(datadir, input_file):
    return {rec.name: rec for rec in parse_seqfiles([datadir / input_file])}


def _parent_of(rec, parents):
    """The output record keeps the name of its parent (written as the GenBank ACCESSION)."""
    assert rec.name in parents, f"{rec.id} does not keep the name of its parent contig, has name {rec.name}"
    return parents[rec.name]


def _decode(record_id, parent):
    """(first, last, reverse) from a region id, with 1-based coordinates of the first and last bases of the region on the parent."""
    assert record_id.startswith(parent.id), f"{record_id} does not start with its parent id {parent.id}"
    suffix = record_id[len(parent.id):]
    if suffix == "":
        return 1, len(parent), False
    if suffix == "_rc":
        return len(parent), 1, True
    match = _REGION_SUFFIX.fullmatch(suffix)
    assert match is not None, f"can't decode region id {record_id}"
    return int(match.group(1)), int(match.group(2)), match.group(3) is not None


def _forward_region(seq, low, high):
    """bases low to high (1-based, inclusive), wrapping around the origin if low > high."""
    if low <= high:
        return seq[low - 1:high]
    return seq[low - 1:] + seq[:high]


def _expected_seq(record_id, parent):
    first, last, reverse = _decode(record_id, parent)
    seq = str(parent.seq)
    if reverse:
        return str(Seq(_forward_region(seq, last, first)).reverse_complement())
    return _forward_region(seq, first, last)


def _check_name_decodes_to_seq(records):
    assert len(records) > 0
    for rec, parent in records:
        assert str(rec.seq) == _expected_seq(rec.id, parent), f"the name {rec.id} doesn't match its sequence"


def _all_cases(datadir, output_dir):
    """(case id, list of (record, parent)) for every case except splice cases."""
    for case_id, tool, input_file, args in CASES:
        parents = _parents(datadir, input_file)
        yield case_id, [(rec, _parent_of(rec, parents)) for rec in _run(datadir, output_dir, case_id, tool, input_file, args)]
    for case_id, input_file, rows in COORD_CASES:
        parents = _parents(datadir, input_file)
        yield case_id, [(rec, _parent_of(rec, parents)) for rec in _run_coords(datadir, output_dir, case_id, input_file, rows)]


@pytest.mark.parametrize("case_id,tool,input_file,args", CASES, ids=[case[0] for case in CASES])
def test_name_decodes_to_seq(shared_datadir, case_id, tool, input_file, args):
    parents = _parents(shared_datadir, input_file)
    with tempfile.TemporaryDirectory() as output_dir:
        records = _run(shared_datadir, output_dir, case_id, tool, input_file, args)
    _check_name_decodes_to_seq([(rec, _parent_of(rec, parents)) for rec in records])


@pytest.mark.parametrize("case_id,input_file,rows", COORD_CASES, ids=[case[0] for case in COORD_CASES])
def test_select_by_coord_name_decodes_to_seq(shared_datadir, case_id, input_file, rows):
    parents = _parents(shared_datadir, input_file)
    with tempfile.TemporaryDirectory() as output_dir:
        records = _run_coords(shared_datadir, output_dir, case_id, input_file, rows)
    assert len(records) == len(rows)
    _check_name_decodes_to_seq([(rec, _parent_of(rec, parents)) for rec in records])


@pytest.mark.parametrize("case_id,tool,input_file,args", SPLICE_CASES, ids=[case[0] for case in SPLICE_CASES])
def test_splice_name_marks_first_and_last_bases(shared_datadir, case_id, tool, input_file, args):
    parents = _parents(shared_datadir, input_file)
    with tempfile.TemporaryDirectory() as output_dir:
        records = _run(shared_datadir, output_dir, case_id, tool, input_file, args)
    assert len(records) > 0
    for rec in records:
        parent = _parent_of(rec, parents)
        first, last, reverse = _decode(rec.id, parent)
        first_base, last_base = parent.seq[first - 1], parent.seq[last - 1]
        if reverse:
            first_base, last_base = str(Seq(first_base).complement()), str(Seq(last_base).complement())
        assert (rec.seq[0], rec.seq[-1]) == (first_base, last_base), rec.id


def test_same_region_same_name_across_tools(shared_datadir):
    """Regions extracted by different tools with the same sequence from the same parent have the same name."""
    names = defaultdict(set)
    with tempfile.TemporaryDirectory() as output_dir:
        for case_id, records in _all_cases(shared_datadir, output_dir):
            for rec, parent in records:
                names[(parent.id, str(rec.seq))].add((rec.id, case_id))
    for (parent_id, _), named in names.items():
        assert len({rec_id for rec_id, _ in named}) == 1, f"different names for the same region of {parent_id}: {sorted(named)}"


@pytest.mark.parametrize("case_id,tool,input_file,args", [case for case in CASES if case[1] is not extract_unannotated], ids=[case[0] for case in CASES if case[1] is not extract_unannotated])
def test_select_by_coord_reproduces_region(shared_datadir, case_id, tool, input_file, args):
    """Feeding a region's name back into select_by_coord gives the same name and sequence."""
    parents = _parents(shared_datadir, input_file)
    with tempfile.TemporaryDirectory() as output_dir:
        records = _run(shared_datadir, output_dir, case_id, tool, input_file, args)
        rows = list()
        for rec in records:
            parent = _parent_of(rec, parents)
            first, last, reverse = _decode(rec.id, parent)
            rows.append((parent.id, last, first, "-") if reverse else (parent.id, first, last, "+"))
        reextracted = _run_coords(shared_datadir, output_dir, case_id + "_coords", input_file, rows)
    assert [(rec.id, str(rec.seq)) for rec in reextracted] == [(rec.id, str(rec.seq)) for rec in records]


def test_rotate_contigs_matches_select_by_cds_whole_contig(shared_datadir):
    """Rotating on a reverse strand CDS and taking the whole contig around it give the same name for the same rotation."""
    parent = next(parse_seqfiles([shared_datadir / SPR]))
    with tempfile.TemporaryDirectory() as output_dir:
        for keep_direction in ([], ["--keep_direction"]):
            selected = _run(shared_datadir, output_dir, "select", select_by_cds, SPR, ["--domains", "DHFR_1", "--whole_contig"] + keep_direction)
            assert len(selected) == 1
            first, _, reverse = _decode(selected[0].id, parent)
            # rotate so that the same base is first
            cut = first if reverse else first - 1
            rotated = rotate_contigs.rotate_record(parent, cut, reverse)
            assert (rotated.id, str(rotated.seq)) == (selected[0].id, str(selected[0].seq))


def _record(length, molecule_type="DNA"):
    seq = ("ACGTTGCAAC" * (length // 10 + 1))[:length]
    rec = SeqRecord(Seq(seq), id="c", name="c")
    rec.annotations["molecule_type"] = molecule_type
    return rec


@pytest.mark.parametrize("location,flip,expected", [
    (FeatureLocation(10, 20, 1), False, "c_11:20"),
    (FeatureLocation(10, 20, None), False, "c_11:20"),
    (FeatureLocation(10, 20, -1), False, "c_20:11rc"),
    (FeatureLocation(10, 20, -1), True, "c_11:20"),
    (FeatureLocation(10, 20, 1), True, "c_20:11rc"),
    (CompoundLocation([FeatureLocation(90, 100, 1), FeatureLocation(0, 10, 1)]), False, "c_91:10"),
    (CompoundLocation([FeatureLocation(0, 10, -1), FeatureLocation(90, 100, -1)]), False, "c_10:91rc"),
    (CompoundLocation([FeatureLocation(0, 10, -1), FeatureLocation(90, 100, -1)]), True, "c_91:10"),
    (FeatureLocation(0, 100, 1), False, "c"),
    (FeatureLocation(0, 100, -1), False, "c_rc"),
    (FeatureLocation(0, 100, -1), True, "c"),
    (CompoundLocation([FeatureLocation(10, 100, 1), FeatureLocation(0, 10, 1)]), False, "c_11:10"),
    (CompoundLocation([FeatureLocation(0, 10, -1), FeatureLocation(10, 100, -1)]), False, "c_10:11rc"),
])
def test_region_id(location, flip, expected):
    rec = _record(100)
    assert region_id(rec, location, flip=flip) == expected
    sliced = slice_record_from_location(rec, location)
    if flip:
        sliced = sliced.reverse_complement()
    assert str(sliced.seq) == _expected_seq(expected, rec)


def test_region_id_protein():
    rec = _record(100, molecule_type="protein")
    assert region_id(rec, FeatureLocation(10, 20)) == "c_11:20"
    assert region_id(rec, FeatureLocation(10, 20, 1)) == "c_11:20"
    assert region_id(rec, FeatureLocation(0, 100)) == "c"
