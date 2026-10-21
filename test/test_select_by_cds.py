import tempfile
from domainator.Bio import SeqIO
from domainator import DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME
from domainator import select_by_cds
from domainator.domainate import main as domainate_main
from domainator.domain_search import main as domain_search_main
from domainator.seqio import parse_seqfiles
from domainator.model import DomainatorCDS
import pytest
from io import StringIO
import sys
import subprocess
from domainator.Bio.Seq import Seq
from domainator.Bio.SeqRecord import SeqRecord
from helpers import compare_seqfiles


def test_select_by_cds_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        select_by_cds.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--domains", "CAT", ])
        recs = list(parse_seqfiles([out]))
        assert len(recs) == 4
        assert [rec.id for rec in recs] == [
            "pDONR201_1_2265:1606rc",
            "pDONR201_2_2265:1606rc",
            "pDONR201_3_2265:1606rc",
            "pDONR201_4_2265:1606rc",
        ]
        assert all(len(DomainatorCDS.list_from_contig(rec)) == 1 for rec in recs)

@pytest.mark.parametrize("cds_range,expected_cds_count",
[
(1,3),
(2,5),
(3,6),
(4,6),
(100000,6),
])
def test_select_by_cds_cds_range_keep_direction_1(cds_range, expected_cds_count, capsys, shared_datadir):
    
    select_by_cds.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "--domains", "CAT", "--cds_range", str(cds_range), "--keep_direction", "-o", "-"])
    captured = capsys.readouterr()
    # with open(f"test_out/select_by_cds{cds_range}_{expected_cds_count}.gb","w") as outfile:
    #     print(captured.out, file=outfile)
    recs = list(parse_seqfiles([StringIO(captured.out)], filetype_override="genbank"))
    for rec in recs:
        cdss = DomainatorCDS.list_from_contig(rec)
        assert len(DomainatorCDS.list_from_contig(rec)) == expected_cds_count
        found = False
        for cds in cdss:
            if cds.num == '2265_-1_1606':
                found = True
                assert cds.feature.location.strand < 0
        assert found

@pytest.mark.parametrize("cds_range,expected_cds_count",
[
(1,3),
(2,5),
(3,6),
(4,6),
(100000,6),
])
def test_select_by_cds_cds_range_1(cds_range, expected_cds_count, capsys, shared_datadir):
    
    select_by_cds.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "--domains", "CAT", "--cds_range", str(cds_range), "-o", "-"])
    captured = capsys.readouterr()
    # with open(f"test_out/select_by_cds{cds_range}_{expected_cds_count}.gb","w") as outfile:
    #     print(captured.out, file=outfile)
    recs = list(parse_seqfiles([StringIO(captured.out)], filetype_override="genbank"))
    for rec in recs:
        cdss = DomainatorCDS.list_from_contig(rec)
        assert len(DomainatorCDS.list_from_contig(rec)) == expected_cds_count
        found = False
        for cds in cdss:
            if cds.num == '2265_-1_1606':
                found = True
                assert cds.feature.location.strand >= 0
        assert found

# expected_cds_count includes CDSs cut by the edge of the region, which are kept with '<' or '>' markers.
@pytest.mark.parametrize("kb_range,expected_cds_count,expected_size",
[
(0.5,3,1660),
(1,4,2660),
(2,5,4265),
(4,6,4470),
])
def test_select_by_cds_kb_range_keep_direction_1(kb_range, expected_cds_count, expected_size, capsys, shared_datadir):
    
    select_by_cds.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "--domains", "CAT", "--kb_range", str(kb_range), "--keep_direction", "-o", "-"])
    captured = capsys.readouterr()
    # with open(f"test_out/select_by_cds{kb_range}_{expected_cds_count}.gb","w") as outfile:
    #     print(captured.out, file=outfile)
    recs = list(parse_seqfiles([StringIO(captured.out)], filetype_override="genbank"))
    for rec in recs:
        assert len(DomainatorCDS.list_from_contig(rec)) == expected_cds_count
        assert len(rec) == expected_size


def test_select_by_cds_3(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        select_by_cds.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--domains", "CAT", "--contigs", "pDONR201_1", "--cds_up", "1", "--keep_direction"])
        # assert 0
        # compare_seqfiles(out, shared_datadir / "extract_peptides_test_1_out.gb")
        # assert compare_files(out, shared_datadir / "extract_peptides_test_1_out.gb")

def test_select_by_cds_4(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        select_by_cds.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--domains", "CAT", "--contigs", "pDONR201_1", "--kb_down", "2", "--keep_direction"])
        # assert 0
        # compare_seqfiles(out, shared_datadir / "extract_peptides_test_1_out.gb")
        # assert compare_files(out, shared_datadir / "extract_peptides_test_1_out.gb")

def test_select_by_cds_5(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        select_by_cds.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--domains", "CAT_MAGIC", "--contigs", "pDONR201_1", "--cds_range", "1", "--keep_direction"])


def test_select_by_cds_align_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        select_by_cds.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_clipped_domainator.gb"), '--whole_contig', "-o", out, "--domains", "CAT", "--pad"])
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 4
        for record in new_file:
            assert len(record) == len(new_file[0])

def test_select_by_cds_align_2(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        select_by_cds.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_clipped_domainator.gb"), '--cds_range', '1', "-o", out, "--domains", "CAT", "--pad"])
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 4
        for record in new_file:
            assert len(record) == len(new_file[0])

def _cds_anchor(location, where):
    """The coordinate of the start, end, or center of a CDS, in the CDS's own orientation."""
    start, end = int(location.start), int(location.end)
    if location.strand == -1:
        return {"start": end, "end": start, "center": end - (end - start) // 2}[where]
    return {"start": start, "end": end, "center": start + (end - start) // 2}[where]

@pytest.mark.parametrize("where", ["start", "end", "center"])
@pytest.mark.parametrize("keep_direction", [False, True])
def test_select_by_cds_pad_where(where, keep_direction, shared_datadir):
    # CAT and CcdB CDSs have different lengths, so start, end, and center alignments give different paddings.
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/extraction.gb"
        args = ["-i", str(shared_datadir / "pDONR201_multi_genemark_clipped_domainator.gb"), "-o", out, "--domains", "CAT", "CcdB", "--pad", "--anchor", where]
        if keep_direction:
            args.append("--keep_direction")
        select_by_cds.main(args)
        new_file = list(parse_seqfiles([out]))
        assert len(new_file) == 8
        anchors = set()
        strands = set()
        for record in new_file:
            assert len(record) == len(new_file[0])
            cdss = DomainatorCDS.list_from_contig(record)
            assert len(cdss) == 1
            location = cdss[0].feature.location
            strands.add(location.strand)
            anchors.add(_cds_anchor(location, where))
        assert len(anchors) == 1
        assert strands == ({1, -1} if keep_direction else {1})

def test_select_by_cds_pad_where_start_is_default(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out_default = output_dir + "/default.gb"
        out_start = output_dir + "/start.gb"
        args = ["-i", str(shared_datadir / "pDONR201_multi_genemark_clipped_domainator.gb"), "--cds_range", "1", "--domains", "CAT", "CcdB", "--pad"]
        select_by_cds.main(args + ["-o", out_default])
        select_by_cds.main(args + ["-o", out_start, "--anchor", "start"])
        compare_seqfiles(out_default, out_start)

def _padding_record(id, dist_to_start, focus_length, dist_to_end, focus_strand):
    record = SeqRecord(Seq("A" * (dist_to_start + focus_length + dist_to_end)), id=id)
    record.set_dist_to_start(dist_to_start)
    record.set_dist_to_end(dist_to_end)
    record.set_focus_strand(focus_strand)
    return record

@pytest.mark.parametrize("where,expected", [
    # (left pad, right pad) for each of the three records
    ("start", [(7, 0), (7, 5), (0, 4)]),
    ("end", [(2, 5), (12, 0), (0, 4)]),
    ("center", [(5, 0), (9, 1), (0, 2)]),
])
def test_pad_records_where(where, expected):
    # a: forward, focus 5..12 (length 7, center 8, rounded down).
    # b: reverse, focus 2..5 (length 3, start 5, end 2, center 4).
    # c: forward, focus 12..14 (length 2, center 13).
    records = [
        _padding_record("a", 5, 7, 3, 1),
        _padding_record("b", 2, 3, 5, -1),
        _padding_record("c", 12, 2, 4, 1),
    ]
    lengths = [len(x) for x in records]
    select_by_cds.pad_records(records, where=where)
    assert len({len(x) for x in records}) == 1
    for record, length, (left, right) in zip(records, lengths, expected):
        assert str(record.seq) == "N" * left + "A" * length + "N" * right

def test_pad_records_bad_where():
    with pytest.raises(ValueError):
        select_by_cds.pad_records([_padding_record("a", 1, 1, 1, 1)], where="middle")

def test_select_by_cds_domain_expr_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        select_by_cds.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_clipped_domainator.gb"), '--cds_range', '1', "-o", out, "--domain_expr", "CAT", "--pad"])
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 4
        for record in new_file:
            assert len(record) == len(new_file[0])
            
def test_select_by_cds_domain_expr_2(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        select_by_cds.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_clipped_domainator.gb"), '--cds_range', '1', "-o", out, "--domain_expr", "CAT&~CAT", "--pad"])
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 0
        for record in new_file:
            assert len(record) == len(new_file[0])


def test_select_by_cds_domain_expr_config_file(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        config_file = output_dir + "/config.yml"
        with open(config_file, "w") as fh:
            subprocess.run([sys.executable, "-m", "domainator", "select_by_cds", "-i", str(shared_datadir / "pDONR201_multi_genemark_clipped_domainator.gb"), '--cds_range', '1', "-o", out, "--domain_expr", "CAT", "--pad", '--print_config'], stdout=fh, stderr=fh, check=True)
        select_by_cds.main(["--config", config_file])
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 4
        for record in new_file:
            assert len(record) == len(new_file[0])

def test_select_by_cds_invert_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        select_by_cds.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_clipped_domainator.gb"), '--cds_range', '1', "-o", out, "--domain_expr", "CAT", "--pad", "--invert"])
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 17
        for record in new_file:
            if "domainator_Pfam-A" in record.annotations:
                print(record.annotations["domainator_Pfam-A"])
                assert 0
                assert "CAT" not in record.annotations["domainator_Pfam-A"]

@pytest.mark.parametrize("strand", [1, -1])
def test_circular_focus_distances_match_pad_location(shared_datadir, strand):
    contig = next(parse_seqfiles([shared_datadir / "bacillus_phage_SPR_with_annotations.gb"]))
    for feature in contig.features:
        if feature.type != "CDS":
            continue
        envelope = feature.location.circular_envelope(len(contig), strand)
        for pad in [(0, 0), (100, 50), (len(contig), len(contig))]:
            slice_location, distances = select_by_cds.pad_location(contig, envelope, pad[0], pad[1], return_extension_sizes=True)
            assert select_by_cds.circular_focus_distances(slice_location, feature.location, len(contig), strand) == distances


def test_pad_records_circular_becomes_linear():
    from domainator.Bio.Seq import Seq
    from domainator.Bio.SeqRecord import SeqRecord
    records = list()
    for seq, dist_to_start in (("AAAACCCC", 4), ("CCCC", 0)):
        record = SeqRecord(Seq(seq), id="x")
        record.annotations.update(molecule_type="DNA", topology="circular")
        record.set_dist_to_start(dist_to_start)
        record.set_dist_to_end(0)
        record.set_focus_strand(1)
        records.append(record)
    select_by_cds.pad_records(records)
    # only the record that was padded is no longer the circular molecule
    assert [(str(r.seq), r.annotations["topology"]) for r in records] == [("AAAACCCC", "circular"), ("NNNNCCCC", "linear")]


def test_select_by_cds_circular_whole_contig_stays_circular(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/extraction.gb"
        select_by_cds.main(["-i", str(shared_datadir / "bacillus_phage_SPR_with_annotations.gb"), "-o", out, "--domains", "DHFR_1", "--whole_contig"])
        recs = list(parse_seqfiles([out]))
    assert len(recs) == 1
    assert recs[0].annotations["topology"] == "circular"


def test_select_by_cds_circular_cds_range_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        select_by_cds.main(["-i", str(shared_datadir / "bacillus_phage_SPR_with_annotations.gb"), '--cds_range', '1', "-o", out, "--domain_expr", "Thymidylat_synt"])
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        rec = new_file[0]
        assert str(rec.seq) == "ATGTTTATTATAGCGTATCTAACTTTATTTCTAGCAGCTTATCTAATTGCATTAAACATTAATGAGGTCAGATTGATTATTCGAGGAGAAAGTGATACATATAGAAAAGTAAAGAGTGTGCTCGATAATTCAACTATTGAAAACGTGAAACGAAATAAGAATTTGATTTACCTGTTTACTTTGCTCAAGGGAATATCTTTCATTGTCCCTCTGGCTTATATTGGATTAGTTATGCACGATAACATCCTAATGCTTGCGTGGACAGCAGTTTCACTTATTTATGTTGTCCTTAGAATGTTCAAAGTTTTAGATGCATTAGAAGGTGAAAGAATCAAGCAAAACACATATTTTTACTTGTTATTTGTTTGTGGGAATTTTCTTTTTGTTGTATTTTATTTGGCAGATATATTCGTATAAAATATTTAATTAATGCTTTGCTTTTTCAAATACTATTTCCAAAAGGACTGAGAATCATGACTCAATTCGATAAACAATACAGTTCAATCATAAATGACATTATAAATAATGGAATCTCAGACGAAGAATTTCAAGTAAGAACAAAGTGGGACTCAGATGGAACACCGGCTCATACACTAAGCGTGATTAGTAAGCAAATGAGATTCGACAACTCAGAGGTTCCGATTTTAACGACAAAAAAAGTTGCCTGGAAAACAGCCATTAAAGAGTTGCTCTGGATTTGGCAGCTGAAATCTAATGATGTTAATGATTTAAACAAGATGGGCGTACATATTTGGGATCAGTGGAAACAAGAAGACGGCACCATCGGACATGCATATGGATTTCAGCTGGGGAAGAAAAACAGAAGTCTAAATGGAGAAAAAGTGGATCAGGTAGACTATCTTCTTCATCAATTGAAGAACAACCCGTCTTCACGCAGACACATTACAATGCTGTGGAATCCTGATGATTTAGACTCAATGGCCTTAACGCCATGTGTATACGAAACTCAATGGTATGTTAAGCAAGGTAAGCTCCACCTCGAGGTAAGAGCACGGAGCAATGATATGGCGTTGGGGAATCCATTCAATGTATTCCAGTACAATGTGTTGCAGCGCATGATTGCTCAAGTGACTGGTTATAAGCTTGGTGAATATATCTTTAACATTGGGGATTGCCATGTGTACACACGTCATATAGACAATTTGAAAATTCAAATGGAAAGAGAACAGTTTGAAGCACCTGAACTATGGATCAATCCTGAAGTGAAAGATTTTTATCATTTTACAATCGATGATTTCAAATTAATCAACTATAAACATGGGGACAAGCTATTATTTGAGGTAGCGGTTTAATGCTATCTCTTATTGCTTGCTGTGATAAAACTCTGGCCATTGGACATCAAAACAAATTACTGTATCATGTGCCTGCTGACATGAAACATTTCAAAGAAAAAACTGAGGGGAAAATATGTATTCAAGGAAGATCAACATACGAATCAATTATCGGTATGACAGGTAAGCCTCTACAGAATAGAAGGAATATTATACTTACTAGGGATCAGAACTTTAAGCCAGATTATTCATCCTTTGTGTATCATTCAATTGAGGAGGTCTTAAAGCTCATTCAAGGTCAAGTTAACACTGATGAGGAAGTGATGGTGATAGGGGGAAGTATGATTTACAAAGCATTCTTGCCCTACGCTGATAAAGTGTATTTGACAATTGTTGATTCAGAGTCAAATGAAGCAGATTCATATTTCCCAATGTTAGATGATCATTGGAAAGTGACTAATAAACAACATAATGAAGCCGATGAAAAGAACAAATACAATTATTCCTTCCTAACTTTTGAAAATAAATATAGACAAAAATAA"
        assert len(rec.features) == 9
        assert rec.features[0].type == "source"

def test_select_by_cds_circular_cds_range_2(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        select_by_cds.main(["-i", str(shared_datadir / "bacillus_phage_SPR_with_annotations.gb"), '--cds_range', '100', "-o", out, "--domain_expr", "Thymidylat_synt"])
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        rec = new_file[0]
        assert str(rec.seq) == "GAGATGACTACAAAAATCAGGAAAACTGGTCAGAAACGGATTGTAAAAGGGATTGATCAGTTGCCATATACAATTAAAGAAAACATTGCTGATCAAGTGAATGAATTAAAGGGGGAGCTTTTTGGATGAAAGCAAATATTAATGGTATTGAATTTGAGGGTACGCCCGAAGAAATAAACGAATTGATTAATTTACATGGATATAAGAATATGCTGCAAGATGATTTATTAATGAGGAACTTTGAACAGCTTAGCCGAGTGAATAGATCAAGACCGACTAATTTATTTAAAGTAGGGGATTATGATTTCCATGACTCGCCTAAATGCTTAATTATCACTTGATTAGAGTAGTAATAAATCAAAAGACAAATATAAAATAAGGAGATGTTTATTATAGCGTATCTAACTTTATTTCTAGCAGCTTATCTAATTGCATTAAACATTAATGAGGTCAGATTGATTATTCGAGGAGAAAGTGATACATATAGAAAAGTAAAGAGTGTGCTCGATAATTCAACTATTGAAAACGTGAAACGAAATAAGAATTTGATTTACCTGTTTACTTTGCTCAAGGGAATATCTTTCATTGTCCCTCTGGCTTATATTGGATTAGTTATGCACGATAACATCCTAATGCTTGCGTGGACAGCAGTTTCACTTATTTATGTTGTCCTTAGAATGTTCAAAGTTTTAGATGCATTAGAAGGTGAAAGAATCAAGCAAAACACATATTTTTACTTGTTATTTGTTTGTGGGAATTTTCTTTTTGTTGTATTTTATTTGGCAGATATATTCGTATAAAATATTTAATTAATGCTTTGCTTTTTCAAATACTATTTCCAAAAGGACTGAGAATCATGACTCAATTCGATAAACAATACAGTTCAATCATAAATGACATTATAAATAATGGAATCTCAGACGAAGAATTTCAAGTAAGAACAAAGTGGGACTCAGATGGAACACCGGCTCATACACTAAGCGTGATTAGTAAGCAAATGAGATTCGACAACTCAGAGGTTCCGATTTTAACGACAAAAAAAGTTGCCTGGAAAACAGCCATTAAAGAGTTGCTCTGGATTTGGCAGCTGAAATCTAATGATGTTAATGATTTAAACAAGATGGGCGTACATATTTGGGATCAGTGGAAACAAGAAGACGGCACCATCGGACATGCATATGGATTTCAGCTGGGGAAGAAAAACAGAAGTCTAAATGGAGAAAAAGTGGATCAGGTAGACTATCTTCTTCATCAATTGAAGAACAACCCGTCTTCACGCAGACACATTACAATGCTGTGGAATCCTGATGATTTAGACTCAATGGCCTTAACGCCATGTGTATACGAAACTCAATGGTATGTTAAGCAAGGTAAGCTCCACCTCGAGGTAAGAGCACGGAGCAATGATATGGCGTTGGGGAATCCATTCAATGTATTCCAGTACAATGTGTTGCAGCGCATGATTGCTCAAGTGACTGGTTATAAGCTTGGTGAATATATCTTTAACATTGGGGATTGCCATGTGTACACACGTCATATAGACAATTTGAAAATTCAAATGGAAAGAGAACAGTTTGAAGCACCTGAACTATGGATCAATCCTGAAGTGAAAGATTTTTATCATTTTACAATCGATGATTTCAAATTAATCAACTATAAACATGGGGACAAGCTATTATTTGAGGTAGCGGTTTAATGCTATCTCTTATTGCTTGCTGTGATAAAACTCTGGCCATTGGACATCAAAACAAATTACTGTATCATGTGCCTGCTGACATGAAACATTTCAAAGAAAAAACTGAGGGGAAAATATGTATTCAAGGAAGATCAACATACGAATCAATTATCGGTATGACAGGTAAGCCTCTACAGAATAGAAGGAATATTATACTTACTAGGGATCAGAACTTTAAGCCAGATTATTCATCCTTTGTGTATCATTCAATTGAGGAGGTCTTAAAGCTCATTCAAGGTCAAGTTAACACTGATGAGGAAGTGATGGTGATAGGGGGAAGTATGATTTACAAAGCATTCTTGCCCTACGCTGATAAAGTGTATTTGACAATTGTTGATTCAGAGTCAAATGAAGCAGATTCATATTTCCCAATGTTAGATGATCATTGGAAAGTGACTAATAAACAACATAATGAAGCCGATGAAAAGAACAAATACAATTATTCCTTCCTAACTTTTGAAAATAAATATAGACAAAAATAAAAAATATGTATAATATAATAACAAGAGGTTGATGAAATAGCAGTTTTTTATTTAAAATACATTGAACATATGAGACAAGTCGCTGAGTACGAGAAAGACCTGAGATACAAAAGTGCAGCAACTAATTTCTTAAAGATGATCGACCATAATTAAAATACGCGGGTGATTAATTGTTAAAGGATAAAAATAAATTATTAAAGAGTATTGAAAAGATCAACAAACTTGAAGAAGGGTTGTCACTATTTGAAGAAGGTGACGAAGAATATTTAAGTGTATTAGTGAAAATTCAGGGGCTATATGATGAAATCTCAGATACTGCTTTAGAGTGTTTTAAA"
        # 13 features: the CDS and gene spanning the origin of the circular input are split across the ends of the rotated record.
        assert len(rec.features) == 13
        assert rec.features[0].type == "source"
        origin_cds = [f for f in rec.features if f.type == "CDS" and len(f.location.parts) == 2]
        assert len(origin_cds) == 1
        assert str(origin_cds[0].location) == "join{[2387:2552](+), [0:129](+)}" # not cut, so no '<' or '>'

def test_select_by_cds_circular_kb_range_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        select_by_cds.main(["-i", str(shared_datadir / "bacillus_phage_SPR_with_annotations.gb"), '--kb_range', '0.5', "-o", out, "--domain_expr", "Thymidylat_synt"])
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        assert new_file[0].name == "OM236515mod_630:1343rc"
        assert str(new_file[0].seq) == "ATCAAAAGACAAATATAAAATAAGGAGATGTTTATTATAGCGTATCTAACTTTATTTCTAGCAGCTTATCTAATTGCATTAAACATTAATGAGGTCAGATTGATTATTCGAGGAGAAAGTGATACATATAGAAAAGTAAAGAGTGTGCTCGATAATTCAACTATTGAAAACGTGAAACGAAATAAGAATTTGATTTACCTGTTTACTTTGCTCAAGGGAATATCTTTCATTGTCCCTCTGGCTTATATTGGATTAGTTATGCACGATAACATCCTAATGCTTGCGTGGACAGCAGTTTCACTTATTTATGTTGTCCTTAGAATGTTCAAAGTTTTAGATGCATTAGAAGGTGAAAGAATCAAGCAAAACACATATTTTTACTTGTTATTTGTTTGTGGGAATTTTCTTTTTGTTGTATTTTATTTGGCAGATATATTCGTATAAAATATTTAATTAATGCTTTGCTTTTTCAAATACTATTTCCAAAAGGACTGAGAATCATGACTCAATTCGATAAACAATACAGTTCAATCATAAATGACATTATAAATAATGGAATCTCAGACGAAGAATTTCAAGTAAGAACAAAGTGGGACTCAGATGGAACACCGGCTCATACACTAAGCGTGATTAGTAAGCAAATGAGATTCGACAACTCAGAGGTTCCGATTTTAACGACAAAAAAAGTTGCCTGGAAAACAGCCATTAAAGAGTTGCTCTGGATTTGGCAGCTGAAATCTAATGATGTTAATGATTTAAACAAGATGGGCGTACATATTTGGGATCAGTGGAAACAAGAAGACGGCACCATCGGACATGCATATGGATTTCAGCTGGGGAAGAAAAACAGAAGTCTAAATGGAGAAAAAGTGGATCAGGTAGACTATCTTCTTCATCAATTGAAGAACAACCCGTCTTCACGCAGACACATTACAATGCTGTGGAATCCTGATGATTTAGACTCAATGGCCTTAACGCCATGTGTATACGAAACTCAATGGTATGTTAAGCAAGGTAAGCTCCACCTCGAGGTAAGAGCACGGAGCAATGATATGGCGTTGGGGAATCCATTCAATGTATTCCAGTACAATGTGTTGCAGCGCATGATTGCTCAAGTGACTGGTTATAAGCTTGGTGAATATATCTTTAACATTGGGGATTGCCATGTGTACACACGTCATATAGACAATTTGAAAATTCAAATGGAAAGAGAACAGTTTGAAGCACCTGAACTATGGATCAATCCTGAAGTGAAAGATTTTTATCATTTTACAATCGATGATTTCAAATTAATCAACTATAAACATGGGGACAAGCTATTATTTGAGGTAGCGGTTTAATGCTATCTCTTATTGCTTGCTGTGATAAAACTCTGGCCATTGGACATCAAAACAAATTACTGTATCATGTGCCTGCTGACATGAAACATTTCAAAGAAAAAACTGAGGGGAAAATATGTATTCAAGGAAGATCAACATACGAATCAATTATCGGTATGACAGGTAAGCCTCTACAGAATAGAAGGAATATTATACTTACTAGGGATCAGAACTTTAAGCCAGATTATTCATCCTTTGTGTATCATTCAATTGAGGAGGTCTTAAAGCTCATTCAAGGTCAAGTTAACACTGATGAGGAAGTGATGGTGATAGGGGGAAGTATGATTTACAAAGCATTCTTGCCCTACGCTGATAAAGTGTATTTGACAATTGTTGATTCAGAGTCAAATGAAGCAGATTCATATTTCCCAATGTTAGATGATCATTGGAAAGTGACTAATAAACAACATAATGAAGCCGATGAAAAGAACAAATACAATTATTCCTTCCTAACTTTTGAA"
        cdss = DomainatorCDS.list_from_contig(new_file[0])
        assert len(cdss) == 3 # including one cut by the end of the region
        assert str(cdss[-1].feature.location) == "[1339:>1840](+)"

def test_select_by_cds_circular_kb_range_2(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        select_by_cds.main(["-i", str(shared_datadir / "bacillus_phage_SPR_with_annotations.gb"), '--kb_range', '10', "-o", out, "--domain_expr", "Thymidylat_synt"])
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        assert new_file[0].name == "OM236515mod_986:987rc"
        assert str(new_file[0].seq) == "GAGATGACTACAAAAATCAGGAAAACTGGTCAGAAACGGATTGTAAAAGGGATTGATCAGTTGCCATATACAATTAAAGAAAACATTGCTGATCAAGTGAATGAATTAAAGGGGGAGCTTTTTGGATGAAAGCAAATATTAATGGTATTGAATTTGAGGGTACGCCCGAAGAAATAAACGAATTGATTAATTTACATGGATATAAGAATATGCTGCAAGATGATTTATTAATGAGGAACTTTGAACAGCTTAGCCGAGTGAATAGATCAAGACCGACTAATTTATTTAAAGTAGGGGATTATGATTTCCATGACTCGCCTAAATGCTTAATTATCACTTGATTAGAGTAGTAATAAATCAAAAGACAAATATAAAATAAGGAGATGTTTATTATAGCGTATCTAACTTTATTTCTAGCAGCTTATCTAATTGCATTAAACATTAATGAGGTCAGATTGATTATTCGAGGAGAAAGTGATACATATAGAAAAGTAAAGAGTGTGCTCGATAATTCAACTATTGAAAACGTGAAACGAAATAAGAATTTGATTTACCTGTTTACTTTGCTCAAGGGAATATCTTTCATTGTCCCTCTGGCTTATATTGGATTAGTTATGCACGATAACATCCTAATGCTTGCGTGGACAGCAGTTTCACTTATTTATGTTGTCCTTAGAATGTTCAAAGTTTTAGATGCATTAGAAGGTGAAAGAATCAAGCAAAACACATATTTTTACTTGTTATTTGTTTGTGGGAATTTTCTTTTTGTTGTATTTTATTTGGCAGATATATTCGTATAAAATATTTAATTAATGCTTTGCTTTTTCAAATACTATTTCCAAAAGGACTGAGAATCATGACTCAATTCGATAAACAATACAGTTCAATCATAAATGACATTATAAATAATGGAATCTCAGACGAAGAATTTCAAGTAAGAACAAAGTGGGACTCAGATGGAACACCGGCTCATACACTAAGCGTGATTAGTAAGCAAATGAGATTCGACAACTCAGAGGTTCCGATTTTAACGACAAAAAAAGTTGCCTGGAAAACAGCCATTAAAGAGTTGCTCTGGATTTGGCAGCTGAAATCTAATGATGTTAATGATTTAAACAAGATGGGCGTACATATTTGGGATCAGTGGAAACAAGAAGACGGCACCATCGGACATGCATATGGATTTCAGCTGGGGAAGAAAAACAGAAGTCTAAATGGAGAAAAAGTGGATCAGGTAGACTATCTTCTTCATCAATTGAAGAACAACCCGTCTTCACGCAGACACATTACAATGCTGTGGAATCCTGATGATTTAGACTCAATGGCCTTAACGCCATGTGTATACGAAACTCAATGGTATGTTAAGCAAGGTAAGCTCCACCTCGAGGTAAGAGCACGGAGCAATGATATGGCGTTGGGGAATCCATTCAATGTATTCCAGTACAATGTGTTGCAGCGCATGATTGCTCAAGTGACTGGTTATAAGCTTGGTGAATATATCTTTAACATTGGGGATTGCCATGTGTACACACGTCATATAGACAATTTGAAAATTCAAATGGAAAGAGAACAGTTTGAAGCACCTGAACTATGGATCAATCCTGAAGTGAAAGATTTTTATCATTTTACAATCGATGATTTCAAATTAATCAACTATAAACATGGGGACAAGCTATTATTTGAGGTAGCGGTTTAATGCTATCTCTTATTGCTTGCTGTGATAAAACTCTGGCCATTGGACATCAAAACAAATTACTGTATCATGTGCCTGCTGACATGAAACATTTCAAAGAAAAAACTGAGGGGAAAATATGTATTCAAGGAAGATCAACATACGAATCAATTATCGGTATGACAGGTAAGCCTCTACAGAATAGAAGGAATATTATACTTACTAGGGATCAGAACTTTAAGCCAGATTATTCATCCTTTGTGTATCATTCAATTGAGGAGGTCTTAAAGCTCATTCAAGGTCAAGTTAACACTGATGAGGAAGTGATGGTGATAGGGGGAAGTATGATTTACAAAGCATTCTTGCCCTACGCTGATAAAGTGTATTTGACAATTGTTGATTCAGAGTCAAATGAAGCAGATTCATATTTCCCAATGTTAGATGATCATTGGAAAGTGACTAATAAACAACATAATGAAGCCGATGAAAAGAACAAATACAATTATTCCTTCCTAACTTTTGAAAATAAATATAGACAAAAATAAAAAATATGTATAATATAATAACAAGAGGTTGATGAAATAGCAGTTTTTTATTTAAAATACATTGAACATATGAGACAAGTCGCTGAGTACGAGAAAGACCTGAGATACAAAAGTGCAGCAACTAATTTCTTAAAGATGATCGACCATAATTAAAATACGCGGGTGATTAATTGTTAAAGGATAAAAATAAATTATTAAAGAGTATTGAAAAGATCAACAAACTTGAAGAAGGGTTGTCACTATTTGAAGAAGGTGACGAAGAATATTTAAGTGTATTAGTGAAAATTCAGGGGCTATATGATGAAATCTCAGATACTGCTTTAGAGTGTTTTAAA"
        cdss = DomainatorCDS.list_from_contig(new_file[0])
        assert len(cdss) == 5 # including the CDS split across the ends of the rotated record

# pDONR_201_domain_search.gb

def test_domain_search_hit_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        select_by_cds.main(["-i", str(shared_datadir / "pDONR_201_domain_search.gb"), "-o", out, "--annotation_type", "domain_search"])
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        cdss = DomainatorCDS.list_from_contig(new_file[0])
        assert len(cdss) == 1
        assert cdss[0].name == "pDONR201_3"


def test_select_by_cds_include_nucleic_acids_flag(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        annotated = output_dir + "/annotated.gb"
        out_without = output_dir + "/without.gb"
        out_with = output_dir + "/with.gb"

        domainate_main([
            "--input", str(shared_datadir / "simple_dna_target.fna"),
            "--alphabet", "nucleotide",
            "-r", str(shared_datadir / "simple_dna_queries.fna"),
            "--genbank_out", annotated,
            "--evalue", "0.1",
        ])

        select_by_cds.main(["-i", annotated, "-o", out_without, "--domains", "dna_query_1"])
        assert list(SeqIO.parse(out_without, "genbank")) == []

        select_by_cds.main(["-i", annotated, "-o", out_with, "--domains", "dna_query_1", "--include_nucleic_acids"])
        records = list(SeqIO.parse(out_with, "genbank"))

        assert len(records) == 1
        domainator_features = [feature for feature in records[0].features if feature.type == DOMAIN_FEATURE_NAME]
        assert len(domainator_features) == 1
        assert domainator_features[0].qualifiers["name"] == ["dna_query_1"]


def test_select_by_cds_search_hits_selects_nucleic_acid_annotations(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        annotated = output_dir + "/search.gb"
        out = output_dir + "/extraction.gb"

        domain_search_main([
            "--input", str(shared_datadir / "simple_dna_target.fna"),
            "--alphabet", "nucleotide",
            "-r", str(shared_datadir / "simple_dna_queries.fna"),
            "-o", annotated,
            "--whole_contig",
            "--cpu", "1",
            "--evalue", "0.1",
        ])

        select_by_cds.main(["-i", annotated, "-o", out, "--annotation_type", "domain_search"])
        records = list(SeqIO.parse(out, "genbank"))

        assert len(records) == 1
        best_hit_features = [feature for feature in records[0].features if feature.type == DOMAIN_SEARCH_BEST_HIT_NAME]
        assert len(best_hit_features) == 1
        assert best_hit_features[0].qualifiers["name"] == ["dna_query_1"]


def test_select_by_cds_no_overlap_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        select_by_cds.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--domains", "CAT", "CcdB", "CcdA", "--cds_range", "1", "--max_region_overlap", "0.4"])
        
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 8

def test_select_by_cds_strand_f_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        select_by_cds.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--domains", "CAT", "CcdB", "--strand", "f"])
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 0

def test_select_by_cds_strand_r_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        select_by_cds.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--domains", "CAT", "CcdB", "--strand", "r"])
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 8

def test_select_by_cds_strand_f_2(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        select_by_cds.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--domains", "APH", "TCAD9", "--strand", "f"])
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 4

def test_select_by_cds_strand_r_2(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        select_by_cds.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--domains", "APH", "TCAD9", "--strand", "r"])
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 0
        
def test_select_by_contig_database_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_cds.main( ["-i", str(shared_datadir /  "pDONR201_multi_genemark_domainator_multi_hmm_2.gb"), "-o", out, "--databases", "pdonr_hmms_1", "--domains", "APH", "--whole_contig"] )
        seqs = list(SeqIO.parse(out, "genbank"))
        assert len(seqs) == 2
        assert seqs[0].id == "pDONR201_3"
        assert seqs[1].id == "pDONR201_4"

def test_select_by_cds_unannoated_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        select_by_cds.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--unannotated"])
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 8


# --- --partial ---

import helpers as _helpers


def test_select_by_cds_partial(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        locations = dict()
        for partial in ("include", "exclude", "only"):
            out = output_dir + f"/{partial}.gb"
            select_by_cds.main(["-i", str(shared_datadir / _helpers.PARTIAL_FIXTURE), "-o", out, "--partial", partial])
            locations[partial] = [str(f.location) for r in parse_seqfiles([out]) for f in r.features if f.type == "CDS"]
        assert locations["only"] == ["[<0:1886](+)"]
        assert locations["exclude"] == ["[0:501](+)"]
        assert sorted(locations["include"]) == sorted(locations["only"] + locations["exclude"])
        # --invert doesn't apply to --partial: inverting "every CDS" selects nothing
        out = output_dir + "/invert.gb"
        select_by_cds.main(["-i", str(shared_datadir / _helpers.PARTIAL_FIXTURE), "-o", out, "--partial", "only", "--cds", _helpers.COMPLETE_CDS_LOCUS_TAG, "--invert"])
        assert [str(f.location) for r in parse_seqfiles([out]) for f in r.features if f.type == "CDS"] == ["[<0:1886](+)"]


def test_select_by_cds_writes_compressed_output(shared_datadir, tmp_path):
    import gzip
    plain, compressed = str(tmp_path / "out.gb"), str(tmp_path / "out.gb.gz")
    for out in (plain, compressed):
        select_by_cds.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--domains", "CAT"])
    with gzip.open(compressed, "rt") as handle:
        assert handle.read() == open(plain).read()
