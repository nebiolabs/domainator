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


def _lineage_match(tx, taxid, include, exclude):
    """Reference implementation: filter via a full lineage walk."""
    lineage = set(tx.lineage(taxid))
    if include and not lineage.intersection(include):
        return False
    if exclude and lineage.intersection(exclude):
        return False
    return True


def test_compile_taxonomy_filter_matches_lineage(shared_datadir):
    tx = NCBITaxonomy(shared_datadir / "taxdmp", overwrite=False)

    # No filter -> None.
    assert tx.compile_taxonomy_filter(None, None) is None

    # include=Bacteria(2), exclude=Proteobacteria(1224): keeps B. subtilis (1423),
    # drops E. coli (562, which descends from 1224).
    cases = [
        ({2}, None),
        (None, {1224}),
        ({2}, {1224}),
        ({1239}, None),  # Firmicutes
    ]
    for include, exclude in cases:
        allowed = tx.compile_taxonomy_filter(include, exclude)
        for taxid in tx._nodes_tab:
            assert (taxid in allowed) == _lineage_match(tx, taxid, include, exclude), (
                f"mismatch for taxid {taxid} with include={include} exclude={exclude}"
            )

    allowed = tx.compile_taxonomy_filter({2}, {1224})
    assert 1423 in allowed   # Bacillus subtilis kept
    assert 562 not in allowed  # Escherichia coli excluded


def _expr_match(tx, taxid, marker_taxids, predicate):
    """Reference: evaluate `predicate` against which markers are in the lineage."""
    lineage = set(tx.lineage(taxid))
    present = {m: (m in lineage) for m in marker_taxids}
    return predicate(present)


def test_compile_taxonomy_expression_semantics(shared_datadir):
    tx = NCBITaxonomy(shared_datadir / "taxdmp", overwrite=False)

    # Empty / None -> no filter.
    assert tx.compile_taxonomy_expression(None) is None
    assert tx.compile_taxonomy_expression("   ") is None

    # "2 & ~1224": within Bacteria but not Proteobacteria. Equivalent to the
    # include={2}/exclude={1224} filter.
    expr_allowed = tx.compile_taxonomy_expression("2 & ~1224")
    filter_allowed = tx.compile_taxonomy_filter({2}, {1224})
    assert expr_allowed == filter_allowed
    assert 1423 in expr_allowed   # Bacillus subtilis (Bacteria, not Proteobacteria)
    assert 562 not in expr_allowed  # Escherichia coli (Proteobacteria)

    # Parity against a reference evaluator for several expressions.
    checks = [
        ("2 | 1224", {2, 1224}, lambda p: p[2] or p[1224]),
        ("~1224", {1224}, lambda p: not p[1224]),
        ("1239 & ~1423", {1239, 1423}, lambda p: p[1239] and not p[1423]),
        ("(2 & ~1224) | 2157", {2, 1224, 2157}, lambda p: (p[2] and not p[1224]) or p[2157]),
    ]
    for expr, markers, predicate in checks:
        allowed = tx.compile_taxonomy_expression(expr)
        for taxid in tx._nodes_tab:
            assert (taxid in allowed) == _expr_match(tx, taxid, markers, predicate), (
                f"mismatch for taxid {taxid} with expr {expr!r}"
            )


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
