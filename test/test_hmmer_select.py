from domainator.hmmer_select import main, hmmer_select
from domainator.utils import pyhmmer_decode
import tempfile
import pyhmmer
import os
import pytest

def test_hmmer_select_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/out.hmm"
        main(['--input', str(shared_datadir/"pdonr_hmms.hmm"), "--output", out, "--field", "all", "--regex", "dehyd.*"])
        output_hmms = list(pyhmmer.plan7.HMMFile(out))
        assert len(output_hmms) == 1

def test_hmmer_select_2(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/out.hmm"
        main(['--input', str(shared_datadir/"pdonr_hmms.hmm"), "--output", out, "--field", "name", "--exact", "TCAD9"])
        output_hmms = list(pyhmmer.plan7.HMMFile(out))
        assert len(output_hmms) == 1

def test_hmmer_select_3(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/out.hmm"
        main(['--input', str(shared_datadir/"pdonr_hmms.hmm"), "--output", out, "--field", "name", "--exact", "TCAD"])
        
        # check that the file size of out is 0
        assert os.path.getsize(out) == 0
def test_hmmer_select_4(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/out.hmm"
        main(['--input', str(shared_datadir/"pdonr_hmms.hmm"), "--output", out, "--field", "acc", "--contains", "PF19974"])
        
        output_hmms = list(pyhmmer.plan7.HMMFile(out))
        assert len(output_hmms) == 1


def test_hmmer_select_case_sensitivity_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/out.hmm"
        main(['--input', str(shared_datadir/"pdonr_hmms.hmm"), "--output", out, "--field", "acc", "--contains", "pf19974"])
        
        output_hmms = list(pyhmmer.plan7.HMMFile(out))
        assert len(output_hmms) == 1

def test_hmmer_select_case_sensitivity_2(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/out.hmm"
        main(['--input', str(shared_datadir/"pdonr_hmms.hmm"), "--output", out, "--field", "acc", "--regex", "pf19974"])
        
        output_hmms = list(pyhmmer.plan7.HMMFile(out))
        assert len(output_hmms) == 1

def test_hmmer_select_case_sensitivity_3(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/out.hmm"
        main(['--input', str(shared_datadir/"pdonr_hmms.hmm"), "--output", out, "--field", "acc", "--contains", "pf19974", "--case_sensitive"])
        
        assert os.path.getsize(out) == 0

def test_hmmer_select_case_sensitivity_4(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/out.hmm"
        main(['--input', str(shared_datadir/"pdonr_hmms.hmm"), "--output", out, "--field", "acc", "--regex", "pf19974", "--case_sensitive"])
        
        assert os.path.getsize(out) == 0

# --- DNA / RNA alphabets --------------------------------------------------

def _profile_names(path):
    # Reading the whole file back is itself the assertion that the output is valid:
    # a mixed-alphabet .hmm raises AlphabetMismatch past its first profile.
    with pyhmmer.plan7.HMMFile(path) as hmm_file:
        return [pyhmmer_decode(hmm.name) for hmm in hmm_file]


def test_hmmer_select_alphabet_splits_mixed_file_set(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out_path = os.path.join(output_dir, "dna_only.hmm")
        main(["-i", str(shared_datadir / "pdonr_hmms.hmm"), str(shared_datadir / "dna_profiles.hmm"),
              "--alphabet", "dna", "-o", out_path])
        assert _profile_names(out_path) == ["dna_prof_1", "dna_prof_2", "dna_prof_3"]

        amino_path = os.path.join(output_dir, "amino_only.hmm")
        main(["-i", str(shared_datadir / "pdonr_hmms.hmm"), str(shared_datadir / "dna_profiles.hmm"),
              "--alphabet", "amino", "-o", amino_path])
        assert "dna_prof_1" not in _profile_names(amino_path)
        assert len(_profile_names(amino_path)) == 7


def test_hmmer_select_alphabet_narrows_text_criteria(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out_path = os.path.join(output_dir, "out.hmm")
        main(["-i", str(shared_datadir / "pdonr_hmms.hmm"), str(shared_datadir / "dna_profiles.hmm"),
              "--alphabet", "dna", "--regex", "prof_2", "-o", out_path])
        assert _profile_names(out_path) == ["dna_prof_2"]


def test_hmmer_select_alphabet_rna(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out_path = os.path.join(output_dir, "out.hmm")
        main(["-i", str(shared_datadir / "dna_profiles.hmm"), str(shared_datadir / "rna_profiles.hmm"),
              "--alphabet", "rna", "-o", out_path])
        assert _profile_names(out_path) == ["rna_prof_1"]


def test_hmmer_select_rejects_mixed_alphabet_inputs(shared_datadir):
    # Without this check the tool writes an .hmm that cannot be read past its first
    # profile, because HMMER locks a file to the alphabet of that profile.
    with tempfile.TemporaryDirectory() as output_dir:
        out_path = os.path.join(output_dir, "out.hmm")
        with pytest.raises(ValueError, match="different alphabets"):
            main(["-i", str(shared_datadir / "dna_profiles.hmm"), str(shared_datadir / "pdonr_hmms.hmm"),
                  "--contains", "prof", "-o", out_path])


def test_hmmer_select_no_criteria_selects_nothing(shared_datadir):
    # Unchanged historical behavior: --alphabet is what makes a criteria-free run useful.
    with tempfile.TemporaryDirectory() as output_dir:
        out_path = os.path.join(output_dir, "out.hmm")
        main(["-i", str(shared_datadir / "dna_profiles.hmm"), "-o", out_path])
        assert os.path.getsize(out_path) == 0
