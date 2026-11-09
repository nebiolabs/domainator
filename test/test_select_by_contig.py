import tempfile

import pytest

from domainator import select_by_contig
from domainator.seqio import parse_seqfiles
import subprocess
import sys

def test_select_by_contig_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main(["-i", str(shared_datadir / "simple_genpept.gb"), "-o", out, "--domains", "CAT", ])
        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 1
        assert seqs[0].id == "pDONR201_4"

def test_select_by_contig_2(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main(["-i", str(shared_datadir / "simple_genpept.gb"), "-o", out, "--max_records", "1"])
        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 1

def test_select_by_contig_3(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main(["-i", str(shared_datadir / "simple_genpept.gb"), "-o", out, "--max_records", "3"])
        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 3


def test_select_by_contig_4(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main(["-i", str(shared_datadir / "simple_genpept.gb"), "-o", out])
        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 5

def test_select_by_contig_5(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main(["-i", str(shared_datadir / "simple_genpept.gb"), "--contigs_file", str(shared_datadir / "simple_genpept_contigs.txt"), "-o", out])
        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 2

def test_select_by_contig_6(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main(["-i", str(shared_datadir / "simple_genpept.gb"), "--min_length", "210", "--max_length", "230","-o", out])
        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 1

def test_select_by_contig_cds_count_lb_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "--min_cds", "6", "-o", out])
        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 4

def test_select_by_contig_cds_count_ub_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "--max_cds", "5", "-o", out])
        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 0

def test_select_by_conting_domain_expression_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main(["-i", str(shared_datadir / "simple_genpept.gb"), "--domain_expr", "(CcdB) | (CAT & Condensation)", "-o", out])
        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 2

def test_select_by_contig_definition_regex_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main( ["-i", str(shared_datadir / "FeSOD_20.gb"), "--definition_regex", "A0A1.*", "-o", out] )
        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 3


def test_select_by_conting_domain_expression_config_file(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        config_file = output_dir + "/config.yml"
        with open(config_file, "w") as fh:
            subprocess.run([sys.executable, "-m", "domainator", "select_by_contig", "-i", str(shared_datadir / "simple_genpept.gb"), "--domain_expr", "(CcdB) | (CAT & Condensation)", "-o", out, '--print_config'], stdout=fh, check=True)
        select_by_contig.main(["--config", config_file])
        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 2

def test_select_by_contig_evalue_filter_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main([
            "-i", str(shared_datadir / "simple_genpept.gb"),
            "-o", out,
            "--domains", "CAT",
            "--evalue", "1e-200",
        ])

        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 0


def test_select_by_contig_invert_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main([
            "-i", str(shared_datadir / "simple_genpept.gb"),
            "-o", out,
            "--domains", "CAT",
            "--invert",
        ])

        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 4
        assert [seq.id for seq in seqs] == ["pDONR201_1", "pDONR201_2", "pDONR201_3", "pDONR201_5"]


def test_select_by_contig_phylogeny_filter_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main([
            "-i", str(shared_datadir / "swissprot_CuSOD_subset.fasta"),
            "-o", out,
            "--include_taxids", "1224",
            "--ncbi_taxonomy_path", str(shared_datadir / "taxdmp"),
        ])

        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 1
        assert seqs[0].id == "sp|P0AGD1|SODC_ECOLI"


def test_select_by_contig_taxonomy_expr(shared_datadir):
    # "1224" (Proteobacteria) is equivalent to --include_taxids 1224.
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main([
            "-i", str(shared_datadir / "swissprot_CuSOD_subset.fasta"),
            "-o", out,
            "--taxonomy_expr", "1224",
            "--ncbi_taxonomy_path", str(shared_datadir / "taxdmp"),
        ])

        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 1
        assert seqs[0].id == "sp|P0AGD1|SODC_ECOLI"

def test_select_by_contig_unannotated_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main(["-i", str(shared_datadir / "simple_genpept.gb"), "--unannotated", "-o", out])
        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 1
        assert seqs[0].id == "pDONR201_1"

def test_select_by_contig_database_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main( ["-i", str(shared_datadir /  "pDONR201_multi_genemark_domainator_multi_hmm_2.gb"), "-o", out, "--databases", "pdonr_hmms_1", "--domains", "APH"] )
        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 2
        assert seqs[0].id == "pDONR201_3"
        assert seqs[1].id == "pDONR201_4"

def test_select_by_contig_fasta_io_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction_pept.fasta"
        select_by_contig.main(["-i", str(shared_datadir / "pdonr_peptides.fasta"), "--fasta_out", out, "--contigs", "pDONR201_3", "pDONR201_5"])

        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 2
        assert seqs[0].id == "pDONR201_3"
        assert seqs[1].id == "pDONR201_5"


def test_select_by_contig_name_regex_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main(["-i", str(shared_datadir / "pdonr_peptides.fasta"), "-o", out, "--contigs_regex", "pDONR201_(3|5)"])

        seqs = list(parse_seqfiles([out]))
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
        select_by_contig.main(["-i", str(shared_datadir / "CuSOD_enum_report_test.gb"), "-o", out, "--min_search_score", "400", "--max_search_score", "500"])

        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 1
        assert seqs[0].id == "sp|O31851|YOJM_BACSU"

def test_select_by_contig_search_name_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main(["-i", str(shared_datadir / "CuSOD_enum_report_test.gb"), "-o", out, "--annotation_type", "domain_search", "--domains", "sp|O31851|YOJM_BACSU"])

        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 1
        assert seqs[0].id == "sp|O31851|YOJM_BACSU"

def test_select_by_contig_search_name_2(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction_pept.gb"
        select_by_contig.main(["-i", str(shared_datadir / "CuSOD_enum_report_test.gb"), "-o", out, "--annotation_type", "domainator", "--domains", "sp|O31851|YOJM_BACSU"])

        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 2


# --- --partial ---

import helpers as _helpers


def test_select_by_contig_partial(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # pDONR201_multi_genemark.gb has GeneMark /partial="10" and "01" codes, read as partial CDSs
        inputs = [str(shared_datadir / _helpers.PARTIAL_FIXTURE), str(shared_datadir / "pDONR201.gb"), str(shared_datadir / "pDONR201_multi_genemark.gb")]
        ids = dict()
        for partial in ("include", "exclude", "only"):
            out = output_dir + f"/{partial}.gb"
            select_by_contig.main(["-i"] + inputs + ["-o", out, "--partial", partial])
            ids[partial] = [r.id for r in parse_seqfiles([out])]
        assert ids["only"] == ["JABFVH010000506_extraction", "pDONR201_1", "pDONR201_2", "pDONR201_3", "pDONR201_4"]
        assert ids["exclude"] == ["pDONR201"]
        assert sorted(ids["include"]) == sorted(ids["only"] + ids["exclude"])
        # --invert doesn't apply to --partial
        out = output_dir + "/invert.gb"
        select_by_contig.main(["-i"] + inputs + ["-o", out, "--partial", "only", "--contigs", "pDONR201_1", "--invert"])
        assert [r.id for r in parse_seqfiles([out])] == ["JABFVH010000506_extraction", "pDONR201_2", "pDONR201_3", "pDONR201_4"]


def test_select_by_contig_partial_uniprot_fasta(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        fasta = output_dir + "/fragments.fasta"
        fragment_id = _helpers.write_uniprot_fragment_fasta(shared_datadir, fasta)
        out = output_dir + "/only.gb"
        select_by_contig.main(["-i", fasta, "-o", out, "--partial", "only"])
        assert [r.id for r in parse_seqfiles([out])] == [fragment_id]


def _records_with_systems():
    from domainator.find_systems import find_systems
    from system_records import contig, layout
    gabija = contig(layout({"GajA": 2, "GajB": 3}), record_id="gabija")
    both = contig(layout({"GajA": 2, "GajB": 3, "RelE": 8, "HicB": 9}), record_id="both")
    neither = contig(layout({"GajA": 2}), record_id="neither")
    records = list(find_systems([gabija, both, neither], domain_expr="GajA & GajB", window_cds=3, system_name="Gabija"))
    records = list(find_systems(records, domain_expr="RelE & HicB", window_cds=3, system_name="TA", clear=False))
    for record in records:
        for feature in record.features:
            if feature.type == "Domain_System" and feature.qualifiers["name"] == ["TA"]:
                feature.qualifiers["category"] = ["toxin_antitoxin"]
                feature.qualifiers["score"] = ["10.0"]
    return records


@pytest.mark.parametrize("kwargs, selected", [
    ({"systems": {"TA"}}, ["both"]),
    ({"systems": {"command_line:Gabija"}}, ["gabija", "both"]),
    ({"system_expr": "Gabija & ~TA"}, ["gabija"]),
    ({"system_categories": {"toxin_antitoxin"}}, ["both"]),
    ({"min_system_score": 50}, ["gabija", "both"]),
    ({"min_system_score": 50, "systems": {"TA"}}, []),
    ({"systems": {"TA"}, "invert": True}, ["gabija", "neither"]),
    ({"system_type": "system_search", "systems": {"TA"}}, []),
])
def test_select_by_systems(kwargs, selected):
    assert [record.id for record in select_by_contig.select_by_contig(_records_with_systems(), **kwargs)] == selected


def test_domain_expr_within_a_window():
    """With a window, the expression must hold within it, and can use <."""
    from system_records import contig
    near = contig(["A", "B", "", "", "", "", ""], record_id="near")
    far = contig(["A", "", "", "", "", "", "B"], record_id="far")
    reversed_order = contig(["B", "A"], record_id="reversed")

    def selected(**kwargs):
        return [r.id for r in select_by_contig.select_by_contig([near, far, reversed_order], **kwargs)]

    assert selected(domain_expr="A & B") == ["near", "far", "reversed"]
    assert selected(domain_expr="A & B", window_cds=2) == ["near", "reversed"]
    assert selected(domain_expr="A < B", window_cds=2) == ["near"]
    assert selected(domain_expr="A & B", window_kb=10) == ["near", "far", "reversed"]
    assert selected(domain_expr="A & B", window_cds=2, invert=True) == ["far"]
    with pytest.raises(select_by_contig.InputError, match="apply to domain_expr"):
        selected(window_cds=2)
