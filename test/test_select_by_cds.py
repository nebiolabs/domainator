import tempfile
from domainator.Bio import SeqIO
from domainator import select_by_cds
from domainator.utils import parse_seqfiles, DomainatorCDS
import pytest
from io import StringIO
import sys
import subprocess


#TODO: replace writing to tempdir with capturing stdout for most tests

def test_select_by_cds_1(shared_datadir):
    #TODO: add some assertions
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        select_by_cds.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--domains", "CAT", ])
        # assert 0
        # compare_seqfiles(out, shared_datadir / "extract_peptides_test_1_out.gb")
        # assert compare_files(out, shared_datadir / "extract_peptides_test_1_out.gb")

@pytest.mark.parametrize("cds_range,expected_cds_count",
[
(1,3),
(2,5),
(3,6),
(4,6),
(100000,6),
])
def test_select_by_cds_cds_range_keep_direction_1(cds_range, expected_cds_count, capsys, shared_datadir):
    
    select_by_cds.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "--domains", "CAT", "--cds_range", str(cds_range), "--keep_direction"])
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
    
    select_by_cds.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "--domains", "CAT", "--cds_range", str(cds_range)])
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

@pytest.mark.parametrize("kb_range,expected_cds_count,expected_size",
[
(0.5,2,1660),
(1,3,2660),
(2,5,4265),
(4,6,4470),
])
def test_select_by_cds_kb_range_keep_direction_1(kb_range, expected_cds_count, expected_size, capsys, shared_datadir):
    
    select_by_cds.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "--domains", "CAT", "--kb_range", str(kb_range), "--keep_direction"])
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
            subprocess.run([sys.executable, select_by_cds.__file__, "-i", str(shared_datadir / "pDONR201_multi_genemark_clipped_domainator.gb"), '--cds_range', '1', "-o", out, "--domain_expr", "CAT", "--pad", '--print_config'], stdout=fh, stderr=fh, check=True)
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
        assert len(rec.features) == 11
        assert rec.features[0].type == "source"

def test_select_by_cds_circular_kb_range_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        select_by_cds.main(["-i", str(shared_datadir / "bacillus_phage_SPR_with_annotations.gb"), '--kb_range', '0.5', "-o", out, "--domain_expr", "Thymidylat_synt"])
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        assert new_file[0].name == "OM236515mod_631:1342rc"
        assert str(new_file[0].seq) == "ATCAAAAGACAAATATAAAATAAGGAGATGTTTATTATAGCGTATCTAACTTTATTTCTAGCAGCTTATCTAATTGCATTAAACATTAATGAGGTCAGATTGATTATTCGAGGAGAAAGTGATACATATAGAAAAGTAAAGAGTGTGCTCGATAATTCAACTATTGAAAACGTGAAACGAAATAAGAATTTGATTTACCTGTTTACTTTGCTCAAGGGAATATCTTTCATTGTCCCTCTGGCTTATATTGGATTAGTTATGCACGATAACATCCTAATGCTTGCGTGGACAGCAGTTTCACTTATTTATGTTGTCCTTAGAATGTTCAAAGTTTTAGATGCATTAGAAGGTGAAAGAATCAAGCAAAACACATATTTTTACTTGTTATTTGTTTGTGGGAATTTTCTTTTTGTTGTATTTTATTTGGCAGATATATTCGTATAAAATATTTAATTAATGCTTTGCTTTTTCAAATACTATTTCCAAAAGGACTGAGAATCATGACTCAATTCGATAAACAATACAGTTCAATCATAAATGACATTATAAATAATGGAATCTCAGACGAAGAATTTCAAGTAAGAACAAAGTGGGACTCAGATGGAACACCGGCTCATACACTAAGCGTGATTAGTAAGCAAATGAGATTCGACAACTCAGAGGTTCCGATTTTAACGACAAAAAAAGTTGCCTGGAAAACAGCCATTAAAGAGTTGCTCTGGATTTGGCAGCTGAAATCTAATGATGTTAATGATTTAAACAAGATGGGCGTACATATTTGGGATCAGTGGAAACAAGAAGACGGCACCATCGGACATGCATATGGATTTCAGCTGGGGAAGAAAAACAGAAGTCTAAATGGAGAAAAAGTGGATCAGGTAGACTATCTTCTTCATCAATTGAAGAACAACCCGTCTTCACGCAGACACATTACAATGCTGTGGAATCCTGATGATTTAGACTCAATGGCCTTAACGCCATGTGTATACGAAACTCAATGGTATGTTAAGCAAGGTAAGCTCCACCTCGAGGTAAGAGCACGGAGCAATGATATGGCGTTGGGGAATCCATTCAATGTATTCCAGTACAATGTGTTGCAGCGCATGATTGCTCAAGTGACTGGTTATAAGCTTGGTGAATATATCTTTAACATTGGGGATTGCCATGTGTACACACGTCATATAGACAATTTGAAAATTCAAATGGAAAGAGAACAGTTTGAAGCACCTGAACTATGGATCAATCCTGAAGTGAAAGATTTTTATCATTTTACAATCGATGATTTCAAATTAATCAACTATAAACATGGGGACAAGCTATTATTTGAGGTAGCGGTTTAATGCTATCTCTTATTGCTTGCTGTGATAAAACTCTGGCCATTGGACATCAAAACAAATTACTGTATCATGTGCCTGCTGACATGAAACATTTCAAAGAAAAAACTGAGGGGAAAATATGTATTCAAGGAAGATCAACATACGAATCAATTATCGGTATGACAGGTAAGCCTCTACAGAATAGAAGGAATATTATACTTACTAGGGATCAGAACTTTAAGCCAGATTATTCATCCTTTGTGTATCATTCAATTGAGGAGGTCTTAAAGCTCATTCAAGGTCAAGTTAACACTGATGAGGAAGTGATGGTGATAGGGGGAAGTATGATTTACAAAGCATTCTTGCCCTACGCTGATAAAGTGTATTTGACAATTGTTGATTCAGAGTCAAATGAAGCAGATTCATATTTCCCAATGTTAGATGATCATTGGAAAGTGACTAATAAACAACATAATGAAGCCGATGAAAAGAACAAATACAATTATTCCTTCCTAACTTTTGAA"
        cdss = DomainatorCDS.list_from_contig(new_file[0])
        assert len(cdss) == 2

def test_select_by_cds_circular_kb_range_2(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        select_by_cds.main(["-i", str(shared_datadir / "bacillus_phage_SPR_with_annotations.gb"), '--kb_range', '10', "-o", out, "--domain_expr", "Thymidylat_synt"])
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        assert new_file[0].name == "OM236515mod_987:986rc"
        assert str(new_file[0].seq) == "GAGATGACTACAAAAATCAGGAAAACTGGTCAGAAACGGATTGTAAAAGGGATTGATCAGTTGCCATATACAATTAAAGAAAACATTGCTGATCAAGTGAATGAATTAAAGGGGGAGCTTTTTGGATGAAAGCAAATATTAATGGTATTGAATTTGAGGGTACGCCCGAAGAAATAAACGAATTGATTAATTTACATGGATATAAGAATATGCTGCAAGATGATTTATTAATGAGGAACTTTGAACAGCTTAGCCGAGTGAATAGATCAAGACCGACTAATTTATTTAAAGTAGGGGATTATGATTTCCATGACTCGCCTAAATGCTTAATTATCACTTGATTAGAGTAGTAATAAATCAAAAGACAAATATAAAATAAGGAGATGTTTATTATAGCGTATCTAACTTTATTTCTAGCAGCTTATCTAATTGCATTAAACATTAATGAGGTCAGATTGATTATTCGAGGAGAAAGTGATACATATAGAAAAGTAAAGAGTGTGCTCGATAATTCAACTATTGAAAACGTGAAACGAAATAAGAATTTGATTTACCTGTTTACTTTGCTCAAGGGAATATCTTTCATTGTCCCTCTGGCTTATATTGGATTAGTTATGCACGATAACATCCTAATGCTTGCGTGGACAGCAGTTTCACTTATTTATGTTGTCCTTAGAATGTTCAAAGTTTTAGATGCATTAGAAGGTGAAAGAATCAAGCAAAACACATATTTTTACTTGTTATTTGTTTGTGGGAATTTTCTTTTTGTTGTATTTTATTTGGCAGATATATTCGTATAAAATATTTAATTAATGCTTTGCTTTTTCAAATACTATTTCCAAAAGGACTGAGAATCATGACTCAATTCGATAAACAATACAGTTCAATCATAAATGACATTATAAATAATGGAATCTCAGACGAAGAATTTCAAGTAAGAACAAAGTGGGACTCAGATGGAACACCGGCTCATACACTAAGCGTGATTAGTAAGCAAATGAGATTCGACAACTCAGAGGTTCCGATTTTAACGACAAAAAAAGTTGCCTGGAAAACAGCCATTAAAGAGTTGCTCTGGATTTGGCAGCTGAAATCTAATGATGTTAATGATTTAAACAAGATGGGCGTACATATTTGGGATCAGTGGAAACAAGAAGACGGCACCATCGGACATGCATATGGATTTCAGCTGGGGAAGAAAAACAGAAGTCTAAATGGAGAAAAAGTGGATCAGGTAGACTATCTTCTTCATCAATTGAAGAACAACCCGTCTTCACGCAGACACATTACAATGCTGTGGAATCCTGATGATTTAGACTCAATGGCCTTAACGCCATGTGTATACGAAACTCAATGGTATGTTAAGCAAGGTAAGCTCCACCTCGAGGTAAGAGCACGGAGCAATGATATGGCGTTGGGGAATCCATTCAATGTATTCCAGTACAATGTGTTGCAGCGCATGATTGCTCAAGTGACTGGTTATAAGCTTGGTGAATATATCTTTAACATTGGGGATTGCCATGTGTACACACGTCATATAGACAATTTGAAAATTCAAATGGAAAGAGAACAGTTTGAAGCACCTGAACTATGGATCAATCCTGAAGTGAAAGATTTTTATCATTTTACAATCGATGATTTCAAATTAATCAACTATAAACATGGGGACAAGCTATTATTTGAGGTAGCGGTTTAATGCTATCTCTTATTGCTTGCTGTGATAAAACTCTGGCCATTGGACATCAAAACAAATTACTGTATCATGTGCCTGCTGACATGAAACATTTCAAAGAAAAAACTGAGGGGAAAATATGTATTCAAGGAAGATCAACATACGAATCAATTATCGGTATGACAGGTAAGCCTCTACAGAATAGAAGGAATATTATACTTACTAGGGATCAGAACTTTAAGCCAGATTATTCATCCTTTGTGTATCATTCAATTGAGGAGGTCTTAAAGCTCATTCAAGGTCAAGTTAACACTGATGAGGAAGTGATGGTGATAGGGGGAAGTATGATTTACAAAGCATTCTTGCCCTACGCTGATAAAGTGTATTTGACAATTGTTGATTCAGAGTCAAATGAAGCAGATTCATATTTCCCAATGTTAGATGATCATTGGAAAGTGACTAATAAACAACATAATGAAGCCGATGAAAAGAACAAATACAATTATTCCTTCCTAACTTTTGAAAATAAATATAGACAAAAATAAAAAATATGTATAATATAATAACAAGAGGTTGATGAAATAGCAGTTTTTTATTTAAAATACATTGAACATATGAGACAAGTCGCTGAGTACGAGAAAGACCTGAGATACAAAAGTGCAGCAACTAATTTCTTAAAGATGATCGACCATAATTAAAATACGCGGGTGATTAATTGTTAAAGGATAAAAATAAATTATTAAAGAGTATTGAAAAGATCAACAAACTTGAAGAAGGGTTGTCACTATTTGAAGAAGGTGACGAAGAATATTTAAGTGTATTAGTGAAAATTCAGGGGCTATATGATGAAATCTCAGATACTGCTTTAGAGTGTTTTAAA"
        cdss = DomainatorCDS.list_from_contig(new_file[0])
        assert len(cdss) == 4

# pDONR_201_domain_search.gb

def test_domain_search_hit_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        select_by_cds.main(["-i", str(shared_datadir / "pDONR_201_domain_search.gb"), "-o", out, "--search_hits"])
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        cdss = DomainatorCDS.list_from_contig(new_file[0])
        assert len(cdss) == 1
        assert cdss[0].name == "pDONR201_3"


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