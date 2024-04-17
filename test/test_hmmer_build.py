import pytest
import tempfile
import os
from io import BytesIO, StringIO
from pyhmmer import easel
from pyhmmer.plan7 import HMM
from domainator.hmmer_build import hmmer_build, main


@pytest.fixture
def msa_content():
    return b">seq1\nACGTACGT\n>seq2\nACGTACGA\n"

def create_msa_file(content):
    with tempfile.NamedTemporaryFile(delete=False, mode="wb") as f:
        f.write(content)
        return f.name

@pytest.fixture
def msa_file(msa_content):
    file_path = create_msa_file(msa_content)
    yield file_path
    os.unlink(file_path)

def test_hmmer_build_with_required_params():
    msa_content = b">seq1\nACGTACGT\n>seq2\nACGTACGA\n"
    msa_file = create_msa_file(msa_content)

    with open(msa_file, "rb") as file:
        hmm = hmmer_build(file, name="test_profile")
        assert isinstance(hmm, HMM)
        assert hmm.name.decode() == "test_profile"

    os.unlink(msa_file)

def test_hmmer_build_with_optional_params():
    msa_content = b">seq1\nACGTACGT\n>seq2\nACGTACGA\n"
    msa_file = create_msa_file(msa_content)

    with open(msa_file, "rb") as file:
        hmm = hmmer_build(file, name="test_profile", acc="P12345", desc="Test profile description", alphabet=easel.Alphabet.dna())
        assert isinstance(hmm, HMM)
        assert hmm.name.decode() == "test_profile"
        assert hmm.accession.decode() == "P12345"
        assert hmm.description.decode() == "Test profile description"

    os.unlink(msa_file)

def test_hmmer_build_with_binaryio():
    msa_content = b">seq1\nACGTACGT\n>seq2\nACGTACGA\n"
    msa_file = BytesIO(msa_content)

    hmm = hmmer_build(msa_file, name="test_profile", acc="P12345", desc="Test profile description", alphabet=easel.Alphabet.dna())
    assert isinstance(hmm, HMM)
    assert hmm.name.decode() == "test_profile"
    assert hmm.accession.decode() == "P12345"
    assert hmm.description.decode() == "Test profile description"


def test_main_with_required_params(msa_file, capsys):
    main(["--name", "test_profile", "--input", msa_file])
    captured = capsys.readouterr()
    assert "HMMER3/f" in captured.out
    assert "NAME  test_profile" in captured.out

def test_main_with_optional_params(msa_file, capsys):
        main(["--name", "test_profile", "--acc", "P12345", "--desc", "Test profile description", "--input", msa_file, "--alphabet", "dna"])
        caputured = capsys.readouterr()
        output = caputured.out
        assert "HMMER3/f" in output
        assert "NAME  test_profile" in output
        assert "ACC   P12345" in output
        assert "DESC  Test profile description" in output

