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
        txt = open(out_txt).read()
        html = open(out_html).read()

        assert "Contig Stats" in txt
        assert "Taxonomy" in txt
        assert "Bacillus" in txt
        assert "Escherichia" in txt
        assert "taxid" in txt

        assert "Summary Report" in html
        assert "Taxonomy" in html
        assert "Bacillus" in html


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


# --- --partial ---

import pytest
import helpers as _helpers


@pytest.mark.parametrize("partial,expected_cdss,expected_domains", [("include", 2, {"HOV79_30120", "HOV79_30125"}), ("exclude", 1, {"HOV79_30120"}), ("only", 1, {"HOV79_30125"})])
def test_summary_report_partial(shared_datadir, partial, expected_cdss, expected_domains):
    with tempfile.TemporaryDirectory() as output_dir:
        annotated = _helpers.domainate_partial_fixture(shared_datadir, output_dir)
        out = output_dir + "/summary.txt"
        table = output_dir + "/domains.tsv"
        summary_report.main(["-i", annotated, "-o", out, "--domains_table", table, "--partial", partial])
        with open(out) as handle:
            text = handle.read()
        assert f"CDSs: {expected_cdss}\n" in text
        assert f"partial CDSs: {1 if partial != 'exclude' else 0}\n" in text
        assert "contigs: 1\n" in text # the contig is counted either way
        with open(table) as handle:
            domains = {line.split("\t")[0] for line in handle.read().splitlines()[1:]}
        assert domains == expected_domains


def test_summary_report_counts_protein_fragments(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        fasta = output_dir + "/fragments.fasta"
        _helpers.write_uniprot_fragment_fasta(shared_datadir, fasta)
        out = output_dir + "/summary.json"
        summary_report.main(["-i", fasta, "--json", out])
        import json
        with open(out) as handle:
            stats = json.load(handle)["contig_stats"]
        assert stats["fragment_proteins"] == 1
        assert stats["partial_cdss"] == 0
