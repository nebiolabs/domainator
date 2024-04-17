import tempfile
from domainator import hmmer_search
from pathlib import Path

#TODO: better tests!

def test_hmmer_search_1(shared_datadir):

    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out_path = output_dir + f"/out_scores.hmm"
        hmmer_search.main(["-i", str(shared_datadir / "pdonr_hmms_1.hmm"),  "-r", str(shared_datadir / "pdonr_hmms.hmm"), "-o", out_path, "--score_cutoff", "13"])
        file_contents = Path(out_path).read_text()
        assert "NAME  CAT" in file_contents
        assert "NAME  2-oxoacid_dh" in file_contents
        assert "NAME  APH" in file_contents
        assert "NAME  CcdA" not in file_contents
        assert "NAME  CcdB" not in file_contents
        assert "NAME  Condensation" not in file_contents
        assert "NAME  TCAD9" not in file_contents
        
def test_hmmer_search_2(shared_datadir):

    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out_path = output_dir + f"/out_scores.hmm"
        hmmer_search.main(["-i", str(shared_datadir / "pdonr_hmms_1.hmm"),  "-r", str(shared_datadir / "pdonr_hmms.hmm"), "-o", out_path, "--score_cutoff", "13", "--max_hits", "2"])
        file_contents = Path(out_path).read_text()
        assert "NAME  CAT" in file_contents
        assert "NAME  2-oxoacid_dh" in file_contents
        assert "NAME  APH" not in file_contents
        assert "NAME  CcdA" not in file_contents
        assert "NAME  CcdB" not in file_contents
        assert "NAME  Condensation" not in file_contents
        assert "NAME  TCAD9" not in file_contents
        