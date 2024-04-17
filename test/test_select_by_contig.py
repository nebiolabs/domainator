import tempfile

from domainator import select_by_contig
from domainator import utils
import subprocess
import sys

#TODO: add more assertions

def test_select_by_contig_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main(["-i", str(shared_datadir / "simple_genpept.gb"), "-o", out, "--domains", "CAT", ])
        # assert 0
        # compare_seqfiles(out, shared_datadir / "extract_peptides_test_1_out.gb")
        # assert compare_files(out, shared_datadir / "extract_peptides_test_1_out.gb")

def test_select_by_contig_2(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main(["-i", str(shared_datadir / "simple_genpept.gb"), "-o", out, "--first", "1"])
        seqs = list(utils.parse_seqfiles([out]))
        assert len(seqs) == 1

def test_select_by_contig_3(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main(["-i", str(shared_datadir / "simple_genpept.gb"), "-o", out, "--first", "3"])
        seqs = list(utils.parse_seqfiles([out]))
        assert len(seqs) == 3


def test_select_by_contig_4(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main(["-i", str(shared_datadir / "simple_genpept.gb"), "-o", out])
        seqs = list(utils.parse_seqfiles([out]))
        assert len(seqs) == 5

def test_select_by_contig_5(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main(["-i", str(shared_datadir / "simple_genpept.gb"), "--contigs_file", str(shared_datadir / "simple_genpept_contigs.txt"), "-o", out])
        seqs = list(utils.parse_seqfiles([out]))
        assert len(seqs) == 2

def test_select_by_contig_6(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main(["-i", str(shared_datadir / "simple_genpept.gb"), "--length_lb", "210", "--length_ub", "230","-o", out])
        seqs = list(utils.parse_seqfiles([out]))
        assert len(seqs) == 1

def test_select_by_conting_domain_expression_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main(["-i", str(shared_datadir / "simple_genpept.gb"), "--domain_expr", "(CcdB) | (CAT & Condensation)", "-o", out])
        seqs = list(utils.parse_seqfiles([out]))
        assert len(seqs) == 2

def test_select_by_contig_definition_regex_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main( ["-i", str(shared_datadir / "FeSOD_20.gb"), "--definition_regex", "A0A1.*", "-o", out] )
        seqs = list(utils.parse_seqfiles([out]))
        assert len(seqs) == 3


def test_select_by_conting_domain_expression_config_file(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        config_file = output_dir + "/config.yml"
        with open(config_file, "w") as fh:
            subprocess.run([sys.executable, select_by_contig.__file__, "-i", str(shared_datadir / "simple_genpept.gb"), "--domain_expr", "(CcdB) | (CAT & Condensation)", "-o", out, '--print_config'], stdout=fh, stderr=fh, check=True)
        select_by_contig.main(["--config", config_file])
        seqs = list(utils.parse_seqfiles([out]))
        assert len(seqs) == 2

#TODO: test evalue

#TODO: test invert

#TODO: test phylogeny filtering
        
def test_select_by_contig_unannotated_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main(["-i", str(shared_datadir / "simple_genpept.gb"), "--unannotated", "-o", out])
        seqs = list(utils.parse_seqfiles([out]))
        assert len(seqs) == 1
        assert seqs[0].id == "pDONR201_1"

def test_select_by_contig_database_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main( ["-i", str(shared_datadir /  "pDONR201_multi_genemark_domainator_multi_hmm_2.gb"), "-o", out, "--databases", "pdonr_hmms_1", "--domains", "APH"] )
        seqs = list(utils.parse_seqfiles([out]))
        assert len(seqs) == 2
        assert seqs[0].id == "pDONR201_3"
        assert seqs[1].id == "pDONR201_4"

def test_select_by_contig_fasta_io_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction_pept.fasta"
        select_by_contig.main(["-i", str(shared_datadir / "pdonr_peptides.fasta"), "-o", out, "--contigs", "pDONR201_3", "pDONR201_5", "--fasta_out"])

        seqs = list(utils.parse_seqfiles([out]))
        assert len(seqs) == 2
        assert seqs[0].id == "pDONR201_3"
        assert seqs[1].id == "pDONR201_5"


def test_select_by_contig_name_regex_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main(["-i", str(shared_datadir / "pdonr_peptides.fasta"), "-o", out, "--contigs_regex", "pDONR201_(3|5)"])

        seqs = list(utils.parse_seqfiles([out]))
        assert len(seqs) == 2
        assert seqs[0].id == "pDONR201_3"
        assert seqs[1].id == "pDONR201_5"
# def test_select_by_cds_2(shared_datadir):
    
#     with tempfile.TemporaryDirectory() as output_dir:
#         output_dir = "test_out"
#         out = output_dir + "/extraction.gb"
#         select_by_cds.main(["", "-g", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--domains", "CAT", "--contigs", "pDONR201_1", "--cds_range", "1", "--keep_direction"])
#         assert 0
#         # compare_seqfiles(out, shared_datadir / "extract_peptides_test_1_out.gb")
#         # assert compare_files(out, shared_datadir / "extract_peptides_test_1_out.gb")

def test_select_by_contig_search_score_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main(["-i", str(shared_datadir / "CuSOD_enum_report_test.gb"), "-o", out, "--search_score_lb", "400", "--search_score_ub", "500"])

        seqs = list(utils.parse_seqfiles([out]))
        assert len(seqs) == 1
        assert seqs[0].id == "sp|O31851|YOJM_BACSU"

def test_select_by_contig_search_name_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main(["-i", str(shared_datadir / "CuSOD_enum_report_test.gb"), "-o", out, "--domain_type", "search", "--domains", "sp|O31851|YOJM_BACSU"])

        seqs = list(utils.parse_seqfiles([out]))
        assert len(seqs) == 1
        assert seqs[0].id == "sp|O31851|YOJM_BACSU"

def test_select_by_contig_search_name_2(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main(["-i", str(shared_datadir / "CuSOD_enum_report_test.gb"), "-o", out, "--domain_type", "domain", "--domains", "sp|O31851|YOJM_BACSU"])

        seqs = list(utils.parse_seqfiles([out]))
        assert len(seqs) == 2
