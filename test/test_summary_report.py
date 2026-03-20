from domainator import summary_report
from domainator.domainate import main as domainate_main
import tempfile

def test_contig_stats_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out_html = output_dir + "/contig_stats_test.html"
        out_txt = output_dir + "/contig_stats_test.txt"
        summary_report.main(["-i", str(shared_datadir / "FeSOD_20_pfam.gb"), "-o", out_txt, "--html", out_html, "--domains", "Sod_Fe_C", "Sod_Fe_N" ])
        for fh in (out_html, out_txt):
            f_txt = open(fh).read()
            assert "Domain Stats" in f_txt
            assert "Sod_Fe_C" in f_txt
            assert "avg score" in f_txt
            assert "100.0" in f_txt
            assert "101" in f_txt



def test_contig_stats_empty_input(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out_html = output_dir + "/contig_stats_test.html"
        out_txt = output_dir + "/contig_stats_test.txt"
        summary_report.main(["-i", str(shared_datadir / "empty.gb"), "-o", out_txt, "--html", out_html, ])
        for fh in (out_html, out_txt):
            f_txt = open(fh).read()
            assert "LOCUS" not in f_txt
            assert "avg score" in f_txt
            assert "Domain Stats" in f_txt
        # assert 0
        # compare_seqfiles(out, shared_datadir / "extract_peptides_test_1_out.gb")
        # assert compare_files(out, shared_datadir / "extract_peptides_test_1_out.gb")

def test_contig_stats_taxonomy_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out_html = output_dir + "/contig_stats_test.html"
        out_txt = output_dir + "/contig_stats_test.txt"
        summary_report.main(["-i", str(shared_datadir / "swissprot_CuSOD_subset.fasta"), "-o", out_txt, "--html", out_html, "--taxonomy", "--ncbi_taxonomy_path", str(shared_datadir / "taxdmp")])
        # TODO: test output


def test_summary_report_database_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out_html = output_dir + "/contig_stats_test.html"
        out_txt = output_dir + "/contig_stats_test.txt"
        summary_report.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator_multi_hmm_2.gb"), "-o", out_txt, "--html", out_html, "--databases", "pdonr_hmms_1"])
        for fh in (out_html, out_txt):
            f_txt = open(fh).read()
            assert "pdonr_hmms_1" in f_txt
            assert "pdonr_hmms_2" not in f_txt


def test_summary_report_nucleotide_annotations(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        annotated = output_dir + "/annotated.gb"
        out_html = output_dir + "/summary.html"
        out_txt = output_dir + "/summary.txt"

        domainate_main([
            "--input", str(shared_datadir / "simple_dna_target.fna"),
            "--fasta_type", "nucleotide",
            "-r", str(shared_datadir / "simple_dna_queries.fna"),
            "--output", annotated,
            "--evalue", "0.1",
        ])

        summary_report.main([
            "-i", annotated,
            "-o", out_txt,
            "--html", out_html,
        ])

        for path in (out_txt, out_html):
            text = open(path).read()
            assert "Domain Stats" in text
            assert "dna_query_1" in text
