from domainator import summary_report
from domainator.domainate import main as domainate_main
from helpers import report_data
from types import SimpleNamespace
import tempfile

def test_contig_stats_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out_html = output_dir + "/contig_stats_test.html"
        out_txt = output_dir + "/contig_stats_test.txt"
        summary_report.main(["-i", str(shared_datadir / "FeSOD_20_pfam.gb"), "-o", out_txt, "--html_out", out_html, "--domains", "Sod_Fe_C", "Sod_Fe_N" ])
        f_txt = open(out_txt).read()
        assert "Domain Frequency" in f_txt
        assert "Sod_Fe_C" in f_txt
        assert "avg score" in f_txt
        assert "100.0" in f_txt
        assert "101" in f_txt

        envelope = report_data(open(out_html).read())
        assert envelope["report"] == "summary_report"
        data = envelope["data"]
        assert data["summary"]["contig_stats"]["contigs"] == 20
        assert "Sod_Fe_C" in {row["domain"] for row in data["summary"]["domain_frequency"]}
        assert data["summary"]["domain_cooccurrence"]["Sod_Fe_C"]["Sod_Fe_C"] == 100.0
        lengths = data["length_counts"]
        assert sum(lengths["counts"]) == 20
        assert lengths["lengths"] == sorted(set(lengths["lengths"]))
        assert data["taxonomy_tree"] is None



def test_contig_stats_empty_input(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out_html = output_dir + "/contig_stats_test.html"
        out_txt = output_dir + "/contig_stats_test.txt"
        summary_report.main(["-i", str(shared_datadir / "empty.gb"), "-o", out_txt, "--html_out", out_html, ])
        f_txt = open(out_txt).read()
        assert "LOCUS" not in f_txt
        assert "avg score" in f_txt
        assert "Domain Frequency" in f_txt
        data = report_data(open(out_html).read())["data"]
        assert data["summary"]["contig_stats"]["contigs"] == 0
        assert data["summary"]["domain_frequency"] == []
        assert data["length_counts"] is None
        # assert 0
        # compare_seqfiles(out, shared_datadir / "extract_peptides_test_1_out.gb")
        # assert compare_files(out, shared_datadir / "extract_peptides_test_1_out.gb")

def test_contig_stats_taxonomy_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out_html = output_dir + "/contig_stats_test.html"
        out_txt = output_dir + "/contig_stats_test.txt"
        summary_report.main(["-i", str(shared_datadir / "swissprot_CuSOD_subset.fasta"), "-o", out_txt, "--html_out", out_html, "--taxonomy", "--ncbi_taxonomy_path", str(shared_datadir / "taxdmp")])
        txt = open(out_txt).read()
        html = open(out_html).read()

        assert "Contig Stats" in txt
        assert "Taxonomy" in txt
        assert "Bacillus" in txt
        assert "Escherichia" in txt
        assert "taxid" in txt

        assert "<title>Summary Report</title>" in html
        data = report_data(html)["data"]
        assert data["summary"]["taxonomy"]

        def names(node):
            yield node["name"]
            for child in node["children"]:
                yield from names(child)

        def total(node):
            return node["count"] + sum(total(child) for child in node["children"])

        tree = data["taxonomy_tree"]
        assert tree["name"] == "All"
        assert total(tree) == data["summary"]["contig_stats"]["contigs"]
        assert {"Bacillus subtilis", "Escherichia coli"} <= set(names(tree))


def test_summary_report_database_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out_html = output_dir + "/contig_stats_test.html"
        out_txt = output_dir + "/contig_stats_test.txt"
        summary_report.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator_multi_hmm_2.gb"), "-o", out_txt, "--html_out", out_html, "--databases", "pdonr_hmms_1"])
        f_txt = open(out_txt).read()
        assert "pdonr_hmms_1" in f_txt
        assert "pdonr_hmms_2" not in f_txt
        frequency = report_data(open(out_html).read())["data"]["summary"]["domain_frequency"]
        assert {row["database"] for row in frequency} == {"pdonr_hmms_1"}


def test_summary_report_nucleotide_annotations(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        annotated = output_dir + "/annotated.gb"
        out_html = output_dir + "/summary.html"
        out_txt = output_dir + "/summary.txt"

        domainate_main([
            "--input", str(shared_datadir / "simple_dna_target.fna"),
            "--alphabet", "nucleotide",
            "-r", str(shared_datadir / "simple_dna_queries.fna"),
            "-o", annotated,
            "--evalue", "0.1",
        ])

        summary_report.main([
            "-i", annotated,
            "-o", out_txt,
            "--html_out", out_html,
        ])

        text = open(out_txt).read()
        assert "Domain Frequency" in text
        assert "dna_query_1" in text
        frequency = report_data(open(out_html).read())["data"]["summary"]["domain_frequency"]
        assert "dna_query_1" in {row["domain"] for row in frequency}



# --- the HTML page's data ---

def test_text_table_aligns_columns():
    from domainator.text_plot import text_table
    table = text_table(["name", "count", "note"], [["a", "8", "x"], ["longer", "12", ""]], [False, True, False])
    assert table == "name    count  note\na           8  x\nlonger     12\n"
    assert text_table(["name", "n"], [], [False, True]) == "name  n\n"


def test_text_report_layout(shared_datadir, tmp_path):
    out = tmp_path / "summary.txt"
    summary_report.main(["-i", str(shared_datadir / "MT_nbs.gb"), "-o", str(out), "--domains", "N6_N4_Mtase"])
    lines = out.read_text().splitlines()
    assert lines[:7] == ["Contig Stats", "", "contigs: 20", "CDSs: 180", "partial CDSs: 0", "protein fragments: 0",
                         "CDSs per 10 kb: 10.61"]
    table = lines.index("Domain Frequency") + 1
    assert lines[table].split() == ["domain", "database", "count", "description", "avg", "score"]
    # Text columns are left-aligned, numbers right-aligned under their headers.
    assert lines[table + 1].startswith("N6_N4_Mtase  ")
    assert len(lines[table + 1]) == len(lines[table])
    assert all(line == line.rstrip() for line in lines if not line.startswith(" "))


def test_length_counts():
    assert summary_report.length_counts([]) is None
    assert summary_report.length_counts([7]) == {"lengths": [7], "counts": [1]}
    assert summary_report.length_counts([140, 95, 113, 95]) == {"lengths": [95, 113, 140], "counts": [2, 1, 1]}


def _taxon(lineage, names, ranks, count):
    return {"count": count, "data": SimpleNamespace(lineage=lineage, names=names, ranks=ranks)}


def test_taxonomy_tree_skips_unranked_and_single_child_taxa():
    taxa = {
        # root > cellular organisms (no rank) > Bacteria (domain) > Bacillota (phylum) > Bacilli (class) > Bacillus (genus)
        1386: _taxon([1, 131567, 2, 1239, 91061, 1386], ["root", "cellular organisms", "Bacteria", "Bacillota", "Bacilli", "Bacillus"],
                     ["no rank", "no rank", "domain", "phylum", "class", "genus"], 3),
        561: _taxon([1, 131567, 2, 1224, 1236, 561], ["root", "cellular organisms", "Bacteria", "Pseudomonadota", "Gammaproteobacteria", "Escherichia"],
                    ["no rank", "no rank", "domain", "phylum", "class", "genus"], 5),
        # An unranked taxon of its own is kept.
        12908: _taxon([1, 12908], ["root", "unclassified sequences"], ["no rank", "no rank"], 1),
    }
    tree = summary_report.taxonomy_tree(taxa)
    assert (tree["taxid"], tree["name"], tree["count"]) == (1, "All", 0)
    bacteria, unclassified = tree["children"]
    # Bacteria (a domain) is kept; the single-child phyla and classes are not.
    assert (bacteria["name"], bacteria["rank"], bacteria["count"]) == ("Bacteria", "domain", 0)
    assert [(c["name"], c["count"], c["children"]) for c in bacteria["children"]] == [("Escherichia", 5, []), ("Bacillus", 3, [])]
    assert (unclassified["taxid"], unclassified["rank"], unclassified["count"]) == (12908, "no rank", 1)


# --- --partial ---

import pytest
import helpers as _helpers


@pytest.mark.parametrize("partial,expected_cdss,expected_domains", [("include", 2, {"HOV79_30120", "HOV79_30125"}), ("exclude", 1, {"HOV79_30120"}), ("only", 1, {"HOV79_30125"})])
def test_summary_report_partial(shared_datadir, partial, expected_cdss, expected_domains):
    with tempfile.TemporaryDirectory() as output_dir:
        annotated = _helpers.domainate_partial_fixture(shared_datadir, output_dir)
        out = output_dir + "/summary.txt"
        table = output_dir + "/domains.tsv"
        summary_report.main(["-i", annotated, "-o", out, "--domains_tsv_out", table, "--partial", partial])
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
        summary_report.main(["-i", fasta, "--json_out", out])
        import json
        with open(out) as handle:
            stats = json.load(handle)["contig_stats"]
        assert stats["fragment_proteins"] == 1
        assert stats["partial_cdss"] == 0


def test_co_occurence_count_is_unaffected_by_focus_domains(shared_datadir, tmp_path):
    # The focus-domain percentages shown in the report used to overwrite those domains'
    # rows of the count matrix.
    import numpy as np
    from domainator.data_matrix import DataMatrix

    gb = str(shared_datadir / "FeSOD_20_pfam.gb")
    plain, focused = tmp_path / "plain.tsv", tmp_path / "focused.tsv"
    plain_fraction, focused_fraction = tmp_path / "plain_fraction.tsv", tmp_path / "focused_fraction.tsv"
    summary_report.main(["-i", gb, "-o", str(tmp_path / "plain.txt"), "--co_occurrence_count_out", str(plain),
                         "--co_occurrence_fraction_out", str(plain_fraction)])
    summary_report.main(["-i", gb, "-o", str(tmp_path / "focused.txt"), "--co_occurrence_count_out", str(focused),
                         "--co_occurrence_fraction_out", str(focused_fraction), "--domains", "Sod_Fe_C", "Sod_Fe_N"])

    counts = DataMatrix.from_file(plain)
    focused_counts = DataMatrix.from_file(focused)
    assert counts.rows == focused_counts.rows
    np.testing.assert_array_equal(focused_counts.toarray(), counts.toarray())
    assert counts.symmetric
    np.testing.assert_array_equal(DataMatrix.from_file(focused_fraction).toarray(), DataMatrix.from_file(plain_fraction).toarray())


def test_summary_report_output_dash_writes_stdout(shared_datadir, tmp_path, capsys, monkeypatch):
    monkeypatch.chdir(tmp_path)
    summary_report.main(["-i", str(shared_datadir / "FeSOD_20_pfam.gb"), "-o", "-"])
    assert "Sod_Fe" in capsys.readouterr().out
    assert not (tmp_path / "-").exists()


def test_tied_domain_counts_keep_first_seen_order(shared_datadir, tmp_path):
    """Equal counts are listed in the order the domains were first seen. pandas' default
    sort is unstable, which reordered them once there were more than 16 domains."""
    import json
    from domainator import DOMAIN_FEATURE_NAME
    from domainator.seqio import parse_seqfiles

    gb = shared_datadir / "MT_nbs.gb"
    first_seen = dict()  # name -> count, in first-seen order
    for record in parse_seqfiles([str(gb)]):
        for feature in record.features:
            if feature.type == DOMAIN_FEATURE_NAME:
                name = feature.qualifiers["name"][0]
                first_seen[name] = first_seen.get(name, 0) + 1
    assert len(first_seen) > 16 and len(set(first_seen.values())) < len(first_seen)

    out = tmp_path / "summary.json"
    summary_report.main(["-i", str(gb), "--json_out", str(out)])
    rows = json.loads(out.read_text())["domain_frequency"]
    assert [row["domain"] for row in rows] == sorted(first_seen, key=lambda name: -first_seen[name])


def test_summary_report_taxid_missing_from_the_database(shared_datadir, tmp_path):
    """MT_nbs.gb has taxids the test taxonomy lacks; they are counted, as unknown."""
    import json
    out = tmp_path / "summary.json"
    summary_report.main(["-i", str(shared_datadir / "MT_nbs.gb"), "--json_out", str(out),
                         "--html_out", str(tmp_path / "summary.html"), "--taxonomy",
                         "--ncbi_taxonomy_path", str(shared_datadir / "taxdmp")])
    data = json.loads(out.read_text())
    assert data["contig_stats"]["contigs"] == 20
    assert {"taxid": 59919, "name": "", "rank": "", "count": 1} in data["taxonomy"]
