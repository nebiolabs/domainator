from pathlib import Path
import tempfile
from glob import glob
from helpers import compare_files
import json
import pytest

from domainator import hmmer_report

def test_hmmer_report_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/hmmer_report.tsv"
        hmmer_report.main(["-i", str(shared_datadir / "pdonr_hmms.hmm"), "-o", out, '--source', '--acc', '--desc', '--length', '--consensus', '--append', 'one', 'int', '1', '--append', 'two', 'float', '2.0', '--append', 'three', 'str', 'three'])
        assert Path(out).is_file()
        
        with open(out) as f:
            lines = f.readlines()
            assert len(lines) == 8
            assert lines[0].strip().split("\t") == ["name","source","acc","desc","length","consensus","one","two","three"]
            assert lines[1].strip() == "2-oxoacid_dh\tpdonr_hmms\tPF00198.25\t2-oxoacid dehydrogenases acyltransferase (catalytic domain)\t233\teqeeervplsgirkaiakrlteskqeiphftlsdevdvtallalrkelkedeakeekakltlldflikavalAlkefPelnasvdeeekeivlkkhvniGvAvatprGLlvPviknadkkslleiakelkelaeraregklkpedleggtftisNlGmlGvtsftPiinppqvaIlgvgrikerpvvkegelvarkvmplslsaDHRvidGaeaarFlntlkkllenpeelll\t1\t2.0\tthree"
            assert lines[-1].strip() == "TCAD9\tpdonr_hmms\tPF19974.1\tTernary complex associated domain 9\t437\tdqvevvrvLtgGrSGaqVlevtvfvkeknqalrhVlKigsaseiakEweAyqrliqpllnalfatIiavsesvlengdqvldelgavvYshagqfagepgeklrsLedlfqealrgpeaadravallerlletllnllYagateeplqtlreelnsrLGpdlvvevkevdseqlvvypdDllqakmssysaseynskvagilvsvelsrlevkvrgprlsavdddvrvevllsggalseleeqgdefleGsvvatranlrlrllkeledelvleetllevdglqlahPfaalrsaLtealearvtssvHGDLNprNiLlaeedrvyLIDfartreggpllsDlAwLevnLlrtvladrldlqellrLqrlLalasrllelealaealagesealakafrllaaiRrfarkqyplerrelwwreylaaLllaahrtLk\t1\t2.0\tthree"


# --- DNA / RNA alphabets --------------------------------------------------

def test_hmmer_report_alphabet_column(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out_path = output_dir + "/out.tsv"
        hmmer_report.main(["-i", str(shared_datadir / "CcdB.hmm"), str(shared_datadir / "dna_profiles_1.hmm"),
                           str(shared_datadir / "rna_profiles.hmm"),
                           "--source", "--alphabet", "--length", "-o", out_path])
        rows = [line.split("\t") for line in Path(out_path).read_text().splitlines()]
        assert rows[0] == ["name", "source", "alphabet", "length"]
        by_name = {row[0]: row for row in rows[1:]}
        assert by_name["CcdB"][2] == "amino"
        assert by_name["dna_prof_1"][2] == "DNA"
        assert by_name["rna_prof_1"][2] == "RNA"


def test_hmmer_report_alphabet_column_json(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        json_path = output_dir + "/out.ndjson"
        hmmer_report.main(["-i", str(shared_datadir / "dna_profiles.hmm"), "--alphabet", "--json", json_path])
        records = [json.loads(line) for line in Path(json_path).read_text().splitlines() if line.strip()]
        assert [rec["name"] for rec in records] == ["dna_prof_1", "dna_prof_2", "dna_prof_3"]
        assert {rec["alphabet"] for rec in records} == {"DNA"}


def test_hmmer_report_dna_consensus_and_length(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out_path = output_dir + "/out.tsv"
        hmmer_report.main(["-i", str(shared_datadir / "dna_profiles_1.hmm"), "--length", "--consensus", "-o", out_path])
        row = Path(out_path).read_text().splitlines()[1].split("\t")
        assert row[0] == "dna_prof_1"
        assert int(row[1]) == 72
        assert set(row[2]) <= set("acgt")
        assert len(row[2]) == 72


def test_hmmer_report_reports_mixed_alphabet_file(shared_datadir, tmp_path):
    mixed = tmp_path / "mixed.hmm"
    mixed.write_bytes(
        (shared_datadir / "dna_profiles_1.hmm").read_bytes()
        + (shared_datadir / "CcdB.hmm").read_bytes()
    )
    with tempfile.TemporaryDirectory() as output_dir:
        with pytest.raises(ValueError, match="more than one alphabet"):
            hmmer_report.main(["-i", str(mixed), "--alphabet", "-o", output_dir + "/out.tsv"])
