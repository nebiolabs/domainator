import tempfile
from domainator import hmmer_search
from pathlib import Path
import pytest
from pyhmmer import easel

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


def test_hmmer_search_max_output_gb_blocks_large_hmm_output(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out_path = output_dir + "/out_scores.hmm"
        with pytest.raises(SystemExit, match="--max_output_gb"):
            hmmer_search.main([
                "-i", str(shared_datadir / "pdonr_hmms_1.hmm"),
                "-r", str(shared_datadir / "pdonr_hmms.hmm"),
                "-o", out_path,
                "--score_cutoff", "13",
                "--max_output_gb", "0.000001",
            ])
        assert not Path(out_path).exists()
        

# --- DNA / RNA alphabets --------------------------------------------------

def test_hmmer_search_dna(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out_path = output_dir + "/out.hmm"
        hmmer_search.main(["-i", str(shared_datadir / "dna_profiles.hmm"), "-r", str(shared_datadir / "dna_profiles_1.hmm"),
                           "-o", out_path, "--score_cutoff", "13", "--cpu", "2"])
        file_contents = Path(out_path).read_text()
        assert "ALPH  DNA" in file_contents
        assert "NAME  dna_prof_1" in file_contents # self hit
        assert "NAME  dna_prof_2" in file_contents # diverged relative
        assert "NAME  dna_prof_3" not in file_contents # unrelated


def test_hmmer_search_rna(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out_path = output_dir + "/out.hmm"
        hmmer_search.main(["-i", str(shared_datadir / "rna_profiles.hmm"), "-r", str(shared_datadir / "rna_profiles.hmm"),
                           "-o", out_path, "--score_cutoff", "13", "--cpu", "2"])
        file_contents = Path(out_path).read_text()
        assert "ALPH  RNA" in file_contents
        assert "NAME  rna_prof_1" in file_contents


def test_hmmer_search_default_score_cutoff_is_alphabet_aware(shared_datadir, capsys):
    with tempfile.TemporaryDirectory() as output_dir:
        hmmer_search.main(["-i", str(shared_datadir / "dna_profiles_1.hmm"), "-r", str(shared_datadir / "dna_profiles_1.hmm"),
                           "-o", output_dir + "/dna.hmm", "--cpu", "2"])
        dna_message = capsys.readouterr().err
        hmmer_search.main(["-i", str(shared_datadir / "CcdB.hmm"), "-r", str(shared_datadir / "CcdB.hmm"),
                           "-o", output_dir + "/aa.hmm", "--cpu", "2"])
        amino_message = capsys.readouterr().err
    assert "default for DNA profiles" in dna_message
    assert "default for amino profiles" in amino_message
    assert f"{hmmer_search.DEFAULT_SCORE_CUTOFF_AMINO:g}" in amino_message
    # the nucleotide default sits lower on the nucleotide score scale
    assert hmmer_search.default_score_cutoff(easel.Alphabet.dna()) < hmmer_search.DEFAULT_SCORE_CUTOFF_AMINO
    assert hmmer_search.default_score_cutoff(easel.Alphabet.amino()) == hmmer_search.DEFAULT_SCORE_CUTOFF_AMINO


def test_hmmer_search_explicit_score_cutoff_is_not_rescaled(shared_datadir):
    # dna_prof_1 vs dna_prof_3 scores just under 5, so a literal 4.0 keeps the unrelated
    # profile while the amino-scale 10.0 drops it. Only the default is ever rescaled.
    with tempfile.TemporaryDirectory() as output_dir:
        loose_path = output_dir + "/loose.hmm"
        hmmer_search.main(["-i", str(shared_datadir / "dna_profiles.hmm"), "-r", str(shared_datadir / "dna_profiles_1.hmm"),
                           "-o", loose_path, "--score_cutoff", "4.0", "--cpu", "2"])
        strict_path = output_dir + "/strict.hmm"
        hmmer_search.main(["-i", str(shared_datadir / "dna_profiles.hmm"), "-r", str(shared_datadir / "dna_profiles_1.hmm"),
                           "-o", strict_path, "--score_cutoff", "10.0", "--cpu", "2"])
        assert "NAME  dna_prof_3" in Path(loose_path).read_text()
        assert "NAME  dna_prof_3" not in Path(strict_path).read_text()


def test_hmmer_search_rejects_mismatched_input_and_reference_alphabets(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out_path = output_dir + "/out.hmm"
        with pytest.raises(ValueError, match="must use the same alphabet"):
            hmmer_search.main(["-i", str(shared_datadir / "dna_profiles_1.hmm"), "-r", str(shared_datadir / "CcdB.hmm"),
                               "-o", out_path, "--cpu", "2"])
        # the check runs before the worker pool, so no partial output survives
        assert not Path(out_path).exists()


def test_hmmer_search_rejects_mismatched_reference_files(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        with pytest.raises(ValueError, match="Mismatched alphabets among the reference hmm files"):
            hmmer_search.main(["-i", str(shared_datadir / "dna_profiles_1.hmm"),
                               "-r", str(shared_datadir / "dna_profiles.hmm"), str(shared_datadir / "CcdB.hmm"),
                               "-o", output_dir + "/out.hmm", "--cpu", "2"])


def test_hmmer_search_rejects_dna_against_rna(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        with pytest.raises(ValueError, match="are DNA, but the reference .* are RNA"):
            hmmer_search.main(["-i", str(shared_datadir / "dna_profiles_1.hmm"), "-r", str(shared_datadir / "rna_profiles.hmm"),
                               "-o", output_dir + "/out.hmm", "--cpu", "2"])
