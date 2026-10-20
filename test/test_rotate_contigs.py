import os
import subprocess
import sys
import tempfile
import warnings

import pytest

from domainator import rotate_contigs, DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME
from domainator.rotate_contigs import origin_cut, rotation_location, rotate_record
from domainator.domain_search import main as domain_search_main
from domainator.Bio import SeqIO
from domainator.Bio.Seq import Seq
from domainator.Bio.SeqRecord import SeqRecord
from domainator.Bio.SeqFeature import FeatureLocation, CompoundLocation, SeqFeature
from domainator.utils import parse_seqfiles, region_id

SPR = "bacillus_phage_SPR_with_annotations.gb" # circular, 2552 bp, all features on the reverse strand
SPR_REVERSED = "bacillus_phage_SPR_with_annotations_reversed.gb" # the same contig, reverse complemented
SPR_LENGTH = 2552
CIRCULAR_SEQUENCE = "GCTAACCGTTAGCGATCGTACGATCGATGCTAGTCCGATTAACCGGTTAGGCTTACCGATGG"


def _rotate(args):
    with tempfile.TemporaryDirectory() as output_dir:
        out = os.path.join(output_dir, "out.gb")
        rotate_contigs.main(args + ["-o", out])
        return list(parse_seqfiles([out]))


def _domains(rec, name=None, feature_type=DOMAIN_FEATURE_NAME):
    return [f for f in rec.features if f.type == feature_type and (name is None or f.qualifiers["name"][0] == name)]


def _intervals(feature):
    return sorted((int(part.start), int(part.end)) for part in feature.location.parts)


def _rotated_seq(seq, cut, reverse):
    if reverse:
        return str(Seq(seq[cut:] + seq[:cut]).reverse_complement())
    return seq[cut:] + seq[:cut]


def _circular_record(seq, features=()):
    rec = SeqRecord(Seq(seq), id="circ", name="circ", description="circular test", features=list(features))
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = "circular"
    return rec


def _write_circular_genbank(path, seq, rec_id="circular_test"):
    rec = SeqRecord(Seq(seq), id=rec_id, name=rec_id, description="circular test")
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = "circular"
    with open(path, "w") as handle:
        SeqIO.write([rec], handle, "genbank")


@pytest.mark.parametrize("location,where,expected", [
    (FeatureLocation(10, 20, 1), "start", (10, 1)),
    (FeatureLocation(10, 20, 1), "end", (20, 1)),
    (FeatureLocation(10, 20, 1), "center", (15, 1)),
    (FeatureLocation(10, 21, 1), "center", (15, 1)),
    (FeatureLocation(10, 20, -1), "start", (20, -1)),
    (FeatureLocation(10, 20, -1), "end", (10, -1)),
    (FeatureLocation(10, 20, -1), "center", (15, -1)),
    (FeatureLocation(10, 21, -1), "center", (16, -1)),
    (FeatureLocation(90, 100, 1), "end", (0, 1)),
    (FeatureLocation(0, 10, -1), "end", (0, -1)),
    (CompoundLocation([FeatureLocation(95, 100, 1), FeatureLocation(0, 5, 1)]), "start", (95, 1)),
    (CompoundLocation([FeatureLocation(95, 100, 1), FeatureLocation(0, 5, 1)]), "end", (5, 1)),
    (CompoundLocation([FeatureLocation(95, 100, 1), FeatureLocation(0, 5, 1)]), "center", (0, 1)),
    (CompoundLocation([FeatureLocation(0, 5, -1), FeatureLocation(95, 100, -1)]), "start", (5, -1)),
    (CompoundLocation([FeatureLocation(0, 5, -1), FeatureLocation(95, 100, -1)]), "end", (95, -1)),
    (CompoundLocation([FeatureLocation(0, 6, -1), FeatureLocation(95, 100, -1)]), "center", (1, -1)),
])
def test_origin_cut(location, where, expected):
    assert origin_cut(location, 100, where) == expected


@pytest.mark.parametrize("cut", [0, 1, 17, 61])
@pytest.mark.parametrize("reverse", [False, True])
def test_rotation_location_sequence(cut, reverse):
    rec = _circular_record(CIRCULAR_SEQUENCE)
    rotated = rotate_record(rec, cut, reverse)
    assert str(rotated.seq) == _rotated_seq(CIRCULAR_SEQUENCE, cut, reverse)
    assert rotated.annotations["topology"] == "circular"
    assert len(rotation_location(cut, len(rec), reverse)) == len(rec)


@pytest.mark.parametrize("where", ["start", "end", "center"])
@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("keep_direction", [False, True])
def test_rotate_on_feature(where, strand, keep_direction):
    # The domain lands at the start, end, or across the origin, in the domain's own orientation.
    length = len(CIRCULAR_SEQUENCE)
    domain = SeqFeature(FeatureLocation(20, 30, strand), type=DOMAIN_FEATURE_NAME, qualifiers={"name": ["dom"], "evalue": ["1e-10"], "score": ["50"]})
    rec = _circular_record(CIRCULAR_SEQUENCE, [domain])
    rotated = next(rotate_contigs.rotate_contigs([rec], domains=["dom"], where=where, keep_direction=keep_direction))
    feature = _domains(rotated, "dom")[0]
    reverse = strand == -1 and not keep_direction
    domain_seq = str(domain.location.extract(rec.seq))
    assert str(feature.location.extract(rotated.seq)) == domain_seq
    out_strand = -1 if strand == -1 and keep_direction else 1
    assert feature.location.parts[0].strand == out_strand
    domain_first = {"start": True, "end": False}.get(where) if out_strand == 1 else {"start": False, "end": True}.get(where)
    if domain_first is True:
        assert _intervals(feature) == [(0, 10)]
    elif domain_first is False:
        assert _intervals(feature) == [(length - 10, length)]
    else:
        assert _intervals(feature) == [(0, 5), (length - 5, length)]
    assert rotated.id.endswith("rc") == reverse


def test_rotated_name():
    rec = _circular_record("A" * 100)
    rec.id = "c"
    assert region_id(rec, rotation_location(0, 100, False)) == "c"
    assert region_id(rec, rotation_location(10, 100, False)) == "c_11:10"
    assert region_id(rec, rotation_location(0, 100, True)) == "c_rc"
    assert region_id(rec, rotation_location(10, 100, True)) == "c_10:11rc"


@pytest.mark.parametrize("where,expected_interval,cut", [
    ("start", [(0, 504)], 1843),
    ("end", [(SPR_LENGTH - 504, SPR_LENGTH)], 1339),
    ("center", [(0, 252), (SPR_LENGTH - 252, SPR_LENGTH)], 1591),
])
def test_rotate_existing_domain(shared_datadir, where, expected_interval, cut):
    original = next(parse_seqfiles([shared_datadir / SPR]))
    recs = _rotate(["-i", str(shared_datadir / SPR), "--domains", "DHFR_1", "--where", where])
    assert len(recs) == 1
    rec = recs[0]
    assert rec.id == f"OM236515mod_{cut}:{cut + 1}rc"
    assert rec.annotations["topology"] == "circular"
    assert str(rec.seq) == _rotated_seq(str(original.seq), cut, True)
    dhfr = _domains(rec, "DHFR_1")[0]
    assert _intervals(dhfr) == expected_interval
    assert all(part.strand == 1 for part in dhfr.location.parts)
    # the domain sequence is preserved
    original_dhfr = _domains(original, "DHFR_1")[0]
    assert str(dhfr.location.extract(rec.seq)) == str(original_dhfr.location.extract(original.seq))
    # No features are cut by the rotation
    assert len(rec.features) == len(original.features)
    for feature in rec.features:
        for part in feature.location.parts:
            assert str(part.start).isdigit() and str(part.end).isdigit()
    # the origin-spanning CDS and domain are intact
    original_ts = _domains(original, "Thymidylat_synt")[0]
    ts = _domains(rec, "Thymidylat_synt")[0]
    assert str(ts.location.extract(rec.seq)) == str(original_ts.location.extract(original.seq))
    original_cds = sorted(str(f.location.extract(original.seq)) for f in original.features if f.type == "CDS")
    rotated_cds = sorted(str(f.location.extract(rec.seq)) for f in rec.features if f.type == "CDS")
    assert rotated_cds == original_cds


def test_rotate_existing_domain_keep_direction(shared_datadir):
    recs = _rotate(["-i", str(shared_datadir / SPR), "--domains", "DHFR_1", "--keep_direction"])
    rec = recs[0]
    assert rec.id == "OM236515mod_1844:1843"
    dhfr = _domains(rec, "DHFR_1")[0]
    # the domain is on the reverse strand, so its start is at the end of the contig
    assert _intervals(dhfr) == [(SPR_LENGTH - 504, SPR_LENGTH)]
    assert dhfr.location.strand == -1


def test_rotate_forward_strand_domain_spanning_origin(shared_datadir):
    # Thymidylat_synt is join(2438..2552,1..707) in the reversed fixture
    recs = _rotate(["-i", str(shared_datadir / SPR_REVERSED), "--domains", "Thymidylat_synt"])
    rec = recs[0]
    assert rec.id == "OM236515mod_(reversed)_2438:2437"
    assert _intervals(_domains(rec, "Thymidylat_synt")[0]) == [(0, 822)]


def test_rotate_keep_name(shared_datadir):
    recs = _rotate(["-i", str(shared_datadir / SPR), "--domains", "DHFR_1", "--keep_name"])
    assert recs[0].id == "OM236515mod"


def test_rotate_idempotent(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        first = os.path.join(output_dir, "first.gb")
        second = os.path.join(output_dir, "second.gb")
        rotate_contigs.main(["-i", str(shared_datadir / SPR), "--domains", "DHFR_1", "-o", first])
        rotate_contigs.main(["-i", first, "--domains", "DHFR_1", "--keep_name", "-o", second])
        with open(first) as fh1, open(second) as fh2:
            assert fh1.read() == fh2.read()


def test_rotate_already_at_origin_not_renamed(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        first = os.path.join(output_dir, "first.gb")
        rotate_contigs.main(["-i", str(shared_datadir / SPR), "--domains", "DHFR_1", "--keep_name", "-o", first])
        recs = _rotate(["-i", first, "--domains", "DHFR_1"])
        assert recs[0].id == "OM236515mod"


def test_rotate_highest_score_wins(shared_datadir):
    # Thymidylat_synt has the higher score
    recs = _rotate(["-i", str(shared_datadir / SPR), "--domains", "DHFR_1", "Thymidylat_synt"])
    assert recs[0].id == "OM236515mod_115:116rc"
    assert _intervals(_domains(recs[0], "Thymidylat_synt")[0]) == [(0, 822)]


def test_rotate_evalue_filters_existing(shared_datadir):
    # Thymidylat_synt evalue is 4.1e-93, DHFR_1 is 5.6e-50
    recs = _rotate(["-i", str(shared_datadir / SPR), "--domains", "DHFR_1", "Thymidylat_synt", "-e", "1e-60"])
    assert recs[0].id == "OM236515mod_115:116rc"
    with pytest.warns(UserWarning, match="No qualifying domain"):
        recs = _rotate(["-i", str(shared_datadir / SPR), "--domains", "DHFR_1", "-e", "1e-60"])
    assert len(recs) == 0


def test_rotate_databases(shared_datadir):
    recs = _rotate(["-i", str(shared_datadir / SPR), "--domains", "DHFR_1", "--databases", "SPR"])
    assert len(recs) == 1
    with pytest.warns(UserWarning, match="No qualifying domain"):
        recs = _rotate(["-i", str(shared_datadir / SPR), "--domains", "DHFR_1", "--databases", "Pfam"])
    assert len(recs) == 0


def test_rotate_domains_file(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        domains_file = os.path.join(output_dir, "domains.txt")
        with open(domains_file, "w") as fh:
            fh.write("DHFR_1\n")
        recs = _rotate(["-i", str(shared_datadir / SPR), "--domains_file", domains_file])
    assert recs[0].id == "OM236515mod_1843:1844rc"


def test_rotate_search_hits(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        hits = os.path.join(output_dir, "hits.gb")
        domain_search_main(["-i", str(shared_datadir / "bacillus_phage_SPR.gb"), "-r", str(shared_datadir / "SPR.hmm"), "--whole_contig", "-Z", "1000", "--cpu", "1", "-o", hits])
        hits_recs = list(parse_seqfiles([hits]))
        recs = _rotate(["-i", hits, "--search_hits"])
    assert len(recs) == len(hits_recs) == 2
    for rec in recs:
        search_hits = _domains(rec, feature_type=DOMAIN_SEARCH_BEST_HIT_NAME)
        assert len(search_hits) == 1
        assert int(search_hits[0].location.start) == 0


def test_rotate_search_new_annotation_dropped(shared_datadir):
    recs = _rotate(["-i", str(shared_datadir / "bacillus_phage_SPR.gb"), "-r", str(shared_datadir / "SPR.hmm"), "--cpu", "1"])
    assert recs[0].id == "OM236515mod_115:116rc"
    assert len(_domains(recs[0])) == 0


def test_rotate_search_keep_annotation(shared_datadir):
    recs = _rotate(["-i", str(shared_datadir / "bacillus_phage_SPR.gb"), "-r", str(shared_datadir / "SPR.hmm"), "--cpu", "1", "--keep_annotation"])
    domains = _domains(recs[0])
    # only the chosen hit is kept
    assert len(domains) == 1
    assert domains[0].qualifiers["name"] == ["Thymidylat_synt"]
    assert _intervals(domains[0]) == [(0, 822)]


def test_rotate_search_and_existing_domains(shared_datadir):
    # the new Thymidylat_synt hit outscores the pre-existing DHFR_1 annotation, and the pre-existing annotations are kept.
    recs = _rotate(["-i", str(shared_datadir / SPR), "-r", str(shared_datadir / "SPR.hmm"), "--domains", "DHFR_1", "--cpu", "1"])
    assert recs[0].id == "OM236515mod_115:116rc"
    assert len(_domains(recs[0])) == 2


def test_rotate_fasta_gene_call(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        fasta = os.path.join(output_dir, "spr.fasta")
        rec = next(parse_seqfiles([shared_datadir / "bacillus_phage_SPR.gb"]))
        with open(fasta, "w") as fh:
            fh.write(f">{rec.id}\n{rec.seq}\n")
        with pytest.warns(UserWarning, match="not circular"):
            recs = _rotate(["-i", fasta, "-r", str(shared_datadir / "SPR.hmm"), "--gene_call", "unannotated", "--cpu", "1"])
        assert len(recs) == 0
        recs = _rotate(["-i", fasta, "-r", str(shared_datadir / "SPR.hmm"), "--gene_call", "unannotated", "--circularize", "--cpu", "1", "--keep_annotation"])
    assert len(recs) == 1
    assert recs[0].annotations["topology"] == "circular"
    assert len([f for f in recs[0].features if f.type == "CDS"]) > 0
    assert _intervals(_domains(recs[0], "Thymidylat_synt")[0]) == [(0, 822)]


def test_rotate_linear_genbank_skipped(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        linear = os.path.join(output_dir, "linear.gb")
        rec = next(parse_seqfiles([shared_datadir / SPR]))
        rec.annotations["topology"] = "linear"
        with open(linear, "w") as fh:
            SeqIO.write([rec], fh, "genbank")
        with pytest.warns(UserWarning, match="not circular"):
            recs = _rotate(["-i", linear, "--domains", "DHFR_1"])
        assert len(recs) == 0
        recs = _rotate(["-i", linear, "--domains", "DHFR_1", "--circularize"])
    assert len(recs) == 1
    assert recs[0].annotations["topology"] == "circular"


def test_rotate_protein_skipped(shared_datadir):
    with pytest.warns(UserWarning, match="protein"):
        recs = _rotate(["-i", str(shared_datadir / "pdonr_peptides.fasta"), "--domains", "CAT", "--circularize"])
    assert len(recs) == 0


def test_rotate_nucleotide_query_spanning_origin():
    query_sequence = CIRCULAR_SEQUENCE[-18:] + CIRCULAR_SEQUENCE[:18]
    with tempfile.TemporaryDirectory() as output_dir:
        gb_path = os.path.join(output_dir, "circular.gb")
        query_path = os.path.join(output_dir, "origin_query.fna")
        _write_circular_genbank(gb_path, CIRCULAR_SEQUENCE)
        with open(query_path, "w") as handle:
            handle.write(f">origin_query\n{query_sequence}\n")
        recs = _rotate(["-i", gb_path, "-r", query_path, "--cpu", "1", "-e", "0.1", "--keep_annotation"])
    assert len(recs) == 1
    rec = recs[0]
    assert rec.id == "circular_test_45:44"
    assert str(rec.seq) == _rotated_seq(CIRCULAR_SEQUENCE, len(CIRCULAR_SEQUENCE) - 18, False)
    domain = _domains(rec)[0]
    assert domain.qualifiers["program"] == ["nhmmer"]
    assert _intervals(domain) == [(0, 36)]


def test_rotate_nucleotide_fasta_circularize(shared_datadir):
    recs = _rotate(["-i", str(shared_datadir / "simple_dna_target.fna"), "-r", str(shared_datadir / "simple_dna_queries.fna"), "--circularize", "--cpu", "1", "-e", "0.1", "--keep_annotation"])
    assert len(recs) == 1
    domain = _domains(recs[0])[0]
    assert domain.qualifiers["program"] == ["nhmmer"]
    assert int(domain.location.start) == 0


def test_rotate_infernal(shared_datadir):
    recs = _rotate(["-i", str(shared_datadir / "pANT_R100.fa"), "-r", str(shared_datadir / "RF00042.cm"), "--circularize", "--cpu", "1", "--keep_annotation"])
    assert len(recs) == 1
    domain = _domains(recs[0])[0]
    assert domain.qualifiers["program"] == ["infernal"]
    assert int(domain.location.start) == 0


def test_rotate_fasta_out(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = os.path.join(output_dir, "out.fasta")
        rotate_contigs.main(["-i", str(shared_datadir / SPR), "--domains", "DHFR_1", "--fasta_out", "-o", out])
        recs = list(SeqIO.parse(out, "fasta"))
    assert len(recs) == 1
    assert recs[0].id == "OM236515mod_1843:1844rc"
    assert len(recs[0]) == SPR_LENGTH


def test_rotate_stdin_stdout(shared_datadir, monkeypatch, capsys):
    with open(shared_datadir / SPR) as fh:
        monkeypatch.setattr("sys.stdin", fh)
        rotate_contigs.main(["--domains", "DHFR_1"])
    out = capsys.readouterr().out
    assert out.startswith("LOCUS       OM236515mod_1843:1844rc")
    assert "circular" in out.splitlines()[0]


def test_rotate_config_file(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = os.path.join(output_dir, "out.gb")
        config_file = os.path.join(output_dir, "config.yml")
        with open(config_file, "w") as fh:
            subprocess.run([sys.executable, rotate_contigs.__file__, "-i", str(shared_datadir / SPR), "--domains", "DHFR_1", "--where", "end", "-o", out, "--print_config"], stdout=fh, check=True)
        rotate_contigs.main(["--config", config_file])
        recs = list(parse_seqfiles([out]))
    assert recs[0].id == "OM236515mod_1339:1340rc"


@pytest.mark.parametrize("args,match", [
    ([], "Specify references, domains, or search_hits"),
    (["--domains", "DHFR_1", "--gene_call", "unannotated"], "gene_call only applies"),
    (["--domains", "DHFR_1", "-r", "SPR.hmm", "--gene_call", "all"], "gene_call 'all'"),
    (["--domains", "DHFR_1", "-e", "11"], "evalue must be"),
])
def test_rotate_validation(shared_datadir, args, match):
    args = [str(shared_datadir / arg) if arg.endswith(".hmm") else arg for arg in args]
    with tempfile.TemporaryDirectory() as output_dir:
        out = os.path.join(output_dir, "out.gb")
        with pytest.raises(ValueError, match=match):
            rotate_contigs.main(["-i", str(shared_datadir / SPR), "-o", out] + args)
        assert not os.path.exists(out)


def test_rotate_z0_stdin(shared_datadir, monkeypatch):
    with open(shared_datadir / SPR) as fh:
        monkeypatch.setattr("sys.stdin", fh)
        with pytest.raises(ValueError, match="stdin"):
            rotate_contigs.main(["-r", str(shared_datadir / "SPR.hmm"), "-Z", "0"])


def test_entrypoint():
    assert hasattr(rotate_contigs, "_entrypoint")
