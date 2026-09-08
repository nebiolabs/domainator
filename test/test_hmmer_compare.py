import tempfile
from domainator import hmmer_compare
from helpers import compare_files
from pathlib import Path
import pytest

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
        


# --- DNA / RNA alphabets --------------------------------------------------
#
# test_hmmer_compare_1 above byte-compares against a golden file of protein
# alignments, so it is also the guard that utils.alphabet_score_scale() is exactly
# 1.0 for amino acids and leaves protein output untouched.

def test_hmmer_compare_dna(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out_path = output_dir + "/out_scores.tsv"
        hmmer_compare.main(["-i", str(shared_datadir / "dna_profiles.hmm"), "-r", str(shared_datadir / "dna_profiles.hmm"),
                            "-o", out_path, "--alignments", "--score_cutoff", "13", "--cpu", "2"])
        compare_files(out_path, shared_datadir / "dna_profiles_hmm_scores.tsv")


def test_hmmer_compare_dna_midline_reaches_strong_symbols(shared_datadir):
    # A DNA column score tops out near log(4)=1.39, so with the unscaled protein
    # cutoffs '|' (score > 1.5) was unreachable and every nucleotide midline collapsed
    # into '.' and '+'. This is the regression test for the alphabet-aware thresholds.
    with tempfile.TemporaryDirectory() as output_dir:
        out_path = output_dir + "/out_scores.tsv"
        hmmer_compare.main(["-i", str(shared_datadir / "dna_profiles_1.hmm"), "-r", str(shared_datadir / "dna_profiles_1.hmm"),
                            "-o", out_path, "--alignments", "--cpu", "2"])
        assert "|" in Path(out_path).read_text()


def test_hmmer_compare_rna(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out_path = output_dir + "/out_scores.tsv"
        hmmer_compare.main(["-i", str(shared_datadir / "rna_profiles.hmm"), "-r", str(shared_datadir / "rna_profiles.hmm"),
                            "-o", out_path, "--score_cutoff", "13", "--cpu", "2"])
        rows = [line.split("\t") for line in Path(out_path).read_text().splitlines()[1:] if line]
        assert rows == [["rna_prof_1", "rna_prof_1", "52.72"]]


def test_hmmer_compare_rejects_mismatched_alphabets(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out_path = output_dir + "/out_scores.tsv"
        with pytest.raises(ValueError, match="must use the same alphabet"):
            hmmer_compare.main(["-i", str(shared_datadir / "dna_profiles.hmm"), "-r", str(shared_datadir / "pdonr_hmms.hmm"),
                                "-o", out_path, "--cpu", "2"])
        # reported from the parent process, before the pool starts, leaving no output behind
        assert not Path(out_path).exists()
