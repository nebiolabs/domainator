import os
import shutil
import tempfile
from io import BytesIO

import pytest
import pyinfernal.cm

from domainator.infernal_build import infernal_build, main


CMBUILD_AVAILABLE = shutil.which("cmbuild") is not None

STOCKHOLM_CONTENT = b"# STOCKHOLM 1.0\nseq1         GCAAAC\nseq2         GCUAAC\n#=GC SS_cons <<..>>\n//\n"
FASTA_CONTENT = b">seq1\nGCAAAC\n>seq2\nGCUAAC\n"


@pytest.fixture
def stockholm_file():
    with tempfile.NamedTemporaryFile(delete=False, mode="wb") as handle:
        handle.write(STOCKHOLM_CONTENT)
        path = handle.name
    yield path
    os.unlink(path)


@pytest.fixture
def fasta_file():
    with tempfile.NamedTemporaryFile(delete=False, mode="wb") as handle:
        handle.write(FASTA_CONTENT)
        path = handle.name
    yield path
    os.unlink(path)


@pytest.mark.skipif(not CMBUILD_AVAILABLE, reason="Infernal cmbuild not installed")
def test_infernal_build_with_required_params(stockholm_file):
    cm_text = infernal_build(stockholm_file, name="test_cm")
    assert "INFERNAL1/a" in cm_text
    assert "NAME     test_cm" in cm_text

    with tempfile.NamedTemporaryFile(delete=False, mode="w", encoding="utf-8") as handle:
        handle.write(cm_text)
        output_path = handle.name

    try:
        with pyinfernal.cm.CMFile(output_path) as cm_file:
            cm = cm_file.read()
        assert cm.name == "test_cm"
        assert cm.accession is None
        assert cm.description is None
    finally:
        os.unlink(output_path)


@pytest.mark.skipif(not CMBUILD_AVAILABLE, reason="Infernal cmbuild not installed")
def test_infernal_build_with_optional_params(stockholm_file):
    cm_text = infernal_build(stockholm_file, name="test_cm", acc="RFTEST1", desc="Test CM description")
    assert "NAME     test_cm" in cm_text
    assert "ACC      RFTEST1" in cm_text
    assert "DESC     Test CM description" in cm_text

    with tempfile.NamedTemporaryFile(delete=False, mode="w", encoding="utf-8") as handle:
        handle.write(cm_text)
        output_path = handle.name

    try:
        with pyinfernal.cm.CMFile(output_path) as cm_file:
            cm = cm_file.read()
        assert cm.name == "test_cm"
        assert cm.accession == "RFTEST1"
        assert cm.description == "Test CM description"
    finally:
        os.unlink(output_path)


@pytest.mark.skipif(not CMBUILD_AVAILABLE, reason="Infernal cmbuild not installed")
def test_infernal_build_with_binaryio():
    cm_text = infernal_build(BytesIO(STOCKHOLM_CONTENT), name="test_cm", acc="RFTEST1", desc="Binary IO description")
    assert "NAME     test_cm" in cm_text
    assert "ACC      RFTEST1" in cm_text
    assert "DESC     Binary IO description" in cm_text


@pytest.mark.skipif(not CMBUILD_AVAILABLE, reason="Infernal cmbuild not installed")
def test_infernal_build_with_fasta_input(fasta_file):
    cm_text = infernal_build(fasta_file, name="test_cm", acc="RFTEST1", desc="FASTA input description", noss=True)
    assert "NAME     test_cm" in cm_text
    assert "ACC      RFTEST1" in cm_text
    assert "DESC     FASTA input description" in cm_text

    with tempfile.NamedTemporaryFile(delete=False, mode="w", encoding="utf-8") as handle:
        handle.write(cm_text)
        output_path = handle.name

    try:
        with pyinfernal.cm.CMFile(output_path) as cm_file:
            cm = cm_file.read()
        assert cm.name == "test_cm"
        assert cm.accession == "RFTEST1"
        assert cm.description == "FASTA input description"
    finally:
        os.unlink(output_path)


@pytest.mark.skipif(not CMBUILD_AVAILABLE, reason="Infernal cmbuild not installed")
def test_main_with_required_params(stockholm_file, capsys):
    main(["--name", "test_cm", "--input", stockholm_file])
    captured = capsys.readouterr()
    assert "INFERNAL1/a" in captured.out
    assert "NAME     test_cm" in captured.out


@pytest.mark.skipif(not CMBUILD_AVAILABLE, reason="Infernal cmbuild not installed")
def test_main_with_optional_params(stockholm_file, capsys):
    main(["--name", "test_cm", "--acc", "RFTEST1", "--desc", "Test CM description", "--input", stockholm_file])
    captured = capsys.readouterr()
    output = captured.out
    assert "INFERNAL1/a" in output
    assert "NAME     test_cm" in output
    assert "ACC      RFTEST1" in output
    assert "DESC     Test CM description" in output


@pytest.mark.skipif(not CMBUILD_AVAILABLE, reason="Infernal cmbuild not installed")
def test_main_with_fasta_input(fasta_file, capsys):
    main(["--name", "test_cm", "--acc", "RFTEST1", "--desc", "FASTA input description", "--input", fasta_file, "--noss"])
    captured = capsys.readouterr()
    output = captured.out
    assert "NAME     test_cm" in output
    assert "ACC      RFTEST1" in output
    assert "DESC     FASTA input description" in output


@pytest.mark.skipif(not CMBUILD_AVAILABLE, reason="Infernal cmbuild not installed")
def test_main_writes_output_file(stockholm_file):
    with tempfile.NamedTemporaryFile(delete=False, mode="w", encoding="utf-8") as handle:
        output_path = handle.name

    try:
        main(["--name", "test_cm", "--acc", "RFTEST1", "--desc", "Test CM description", "--input", stockholm_file, "--output", output_path])
        with open(output_path, "r", encoding="utf-8") as output_handle:
            output = output_handle.read()
        assert "NAME     test_cm" in output
        assert "ACC      RFTEST1" in output
        assert "DESC     Test CM description" in output
    finally:
        os.unlink(output_path)