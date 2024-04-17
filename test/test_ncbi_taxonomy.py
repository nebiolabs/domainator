from domainator.Taxonomy import NCBITaxonomy
import pytest


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
