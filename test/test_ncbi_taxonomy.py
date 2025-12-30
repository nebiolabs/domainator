from domainator.Taxonomy import NCBITaxonomy
import pytest
import shutil


def test_NCBItaxonomy(shared_datadir):
    tx = NCBITaxonomy(shared_datadir / "taxdmp", overwrite=False)
	
    assert tx.lineage(562) == [562, 561, 543, 91347, 1236, 1224, 2, 131567, 1]
    assert tx.lineage(1423) == [1423, 653685, 1386, 186817, 1385, 91061, 1239, 1783272, 2, 131567, 1]
    assert tx.rank(562) == "species"
    assert tx.rank(1423) == "species"
    assert tx.name(562) == "Escherichia coli"
    assert tx.name(1423) == "Bacillus subtilis"

def test_NCBItaxonomy_2(shared_datadir):
    tx = NCBITaxonomy(shared_datadir / "taxdmp", overwrite=False)
    with pytest.warns(UserWarning):
        assert tx.lineage(1985417) == []


def test_NCBItaxonomy_extracted_only(tmp_path, shared_datadir):
    """If taxdump.tar.gz is missing, we should still load from extracted .dmp files."""
    src = shared_datadir / "taxdmp"
    for name in ("nodes.dmp", "names.dmp", "merged.dmp", "delnodes.dmp"):
        shutil.copyfile(src / name, tmp_path / name)

    tx = NCBITaxonomy(tmp_path, overwrite=False)
    assert tx.name(562) == "Escherichia coli"
    assert tx.rank(562) == "species"


def test_NCBItaxonomy_fallback_when_tar_corrupt(tmp_path, shared_datadir):
    """If taxdump.tar.gz exists but is unreadable, warn and fall back to extracted .dmp files."""
    src = shared_datadir / "taxdmp"
    for name in ("nodes.dmp", "names.dmp", "merged.dmp", "delnodes.dmp"):
        shutil.copyfile(src / name, tmp_path / name)

    # Create an invalid/empty tarball.
    (tmp_path / "taxdump.tar.gz").write_bytes(b"")

    with pytest.warns(UserWarning):
        tx = NCBITaxonomy(tmp_path, overwrite=False)
    assert tx.name(1423) == "Bacillus subtilis"
