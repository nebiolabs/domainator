import tempfile
from domainator import hmmer_compare
from helpers import compare_files
from pathlib import Path
import pytest

#TODO: better tests!

def test_hmmer_compare_1(shared_datadir):

    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out_path = output_dir + f"/out_scores.tsv"
        hmmer_compare.main(["-i", str(shared_datadir / "pdonr_hmms.hmm"),  "-r", str(shared_datadir / "pdonr_hmms.hmm"), "-o", out_path, "--alignment", "--score_cutoff", "13", "--cpu", "10"])
        compare_files(out_path, shared_datadir / "pDONR_201_hmm_scores.tsv")


def test_hmmer_compare_max_output_gb_blocks_large_tsv_output(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out_path = output_dir + "/out_scores.tsv"
        with pytest.raises(SystemExit, match="--max_output_gb"):
            hmmer_compare.main([
                "-i", str(shared_datadir / "pdonr_hmms.hmm"),
                "-r", str(shared_datadir / "pdonr_hmms.hmm"),
                "-o", out_path,
                "--alignments",
                "--score_cutoff", "13",
                "--cpu", "2",
                "--max_output_gb", "0.000001",
            ])
        assert not Path(out_path).exists()
        
