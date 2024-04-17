from domainator.hmmer_select import main, hmmer_select
import tempfile
import pyhmmer
import os

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