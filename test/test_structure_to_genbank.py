"""Tests for structure_to_genbank.py."""
import shutil

import pytest

from domainator import structure_to_genbank, domainate, utils


pytestmark = pytest.mark.skipif(shutil.which("foldseek") is None, reason="foldseek is not installed")


EXPECTED_NAMES = ["1CRN_A", "1IGD_A", "1NDD_A", "1UBQ_A",
                  "1ZNI_A", "1ZNI_B", "1ZNI_C", "1ZNI_D", "some.long.name_A"]


def _read(path):
    return list(utils.parse_seqfiles((str(path),), default_molecule_type="protein"))


def test_one_record_per_chain(shared_datadir, tmp_path):
    out = tmp_path / "chains.gb"
    structure_to_genbank.main(["-i", str(shared_datadir / "structures" / "inputs"),
                               "-o", str(out), "--cpu", "2"])
    records = _read(out)
    assert [record.id for record in records] == EXPECTED_NAMES
    lengths = {record.id: len(record) for record in records}
    # 1ZNI is a multi-chain structure: each chain becomes its own record
    assert lengths["1ZNI_A"] == 21 and lengths["1ZNI_B"] == 30
    assert lengths["1UBQ_A"] == 76
    for record in records:
        assert record.annotations["molecule_type"] == "protein"
        assert str(record.seq)


def test_long_and_dotted_names_round_trip(shared_datadir, tmp_path):
    """A dotted basename must survive naming, the GenBank LOCUS line, and reparsing."""
    out = tmp_path / "chains.gb"
    structure_to_genbank.main(["-i", str(shared_datadir / "structures" / "inputs" / "some.long.name.pdb.gz"),
                               "-o", str(out), "--cpu", "2"])
    records = _read(out)
    assert [record.id for record in records] == ["some.long.name_A"]
    assert "LOCUS       some.long.name_A" in out.read_text()


def test_provenance_in_description(shared_datadir, tmp_path):
    out = tmp_path / "chains.gb"
    structure_to_genbank.main(["-i", str(shared_datadir / "structures" / "inputs"),
                               "-o", str(out), "--cpu", "2"])
    for record in _read(out):
        assert "structure=" in record.description


def test_store_3di(shared_datadir, tmp_path):
    out = tmp_path / "chains.gb"
    structure_to_genbank.main(["-i", str(shared_datadir / "structures" / "inputs"),
                               "-o", str(out), "--store_3di", "--cpu", "2"])
    records = _read(out)
    for record in records:
        threedi = [f.qualifiers["threedi"][0] for f in record.features if "threedi" in f.qualifiers]
        assert len(threedi) == 1, f"{record.id} should carry exactly one /threedi qualifier"
        # foldseek's 3Di string is per-residue, so it must round trip at full length
        assert len(threedi[0]) == len(record)


def test_no_3di_by_default(shared_datadir, tmp_path):
    out = tmp_path / "chains.gb"
    structure_to_genbank.main(["-i", str(shared_datadir / "structures" / "inputs"),
                               "-o", str(out), "--cpu", "2"])
    assert all(not record.features for record in _read(out))


def test_keep_db_is_reusable_as_input(shared_datadir, tmp_path):
    """The persisted database must carry coordinates and reproduce the same records."""
    from domainator import structure_lib

    first = tmp_path / "first.gb"
    prefix = tmp_path / "kept"
    structure_to_genbank.main(["-i", str(shared_datadir / "structures" / "inputs"),
                               "-o", str(first), "--keep_db", str(prefix), "--cpu", "2"])
    assert structure_lib.is_aligner_db(str(prefix))
    assert structure_lib.StructureDB(prefix=str(prefix)).has_coordinates()

    second = tmp_path / "second.gb"
    structure_to_genbank.main(["-i", str(prefix), "-o", str(second), "--cpu", "2"])
    first_records = _read(first)
    second_records = _read(second)
    assert [r.id for r in first_records] == [r.id for r in second_records]
    assert [str(r.seq) for r in first_records] == [str(r.seq) for r in second_records]


def test_keep_db_rejected_for_database_input(shared_datadir, tmp_path):
    with pytest.raises(RuntimeError, match="already a database"):
        structure_to_genbank.main(["-i", str(shared_datadir / "foldseek" / "FeSOD"),
                                   "-o", str(tmp_path / "x.gb"),
                                   "--keep_db", str(tmp_path / "kept")])


def test_database_input(shared_datadir, tmp_path):
    """A prebuilt sequence-derived database converts too, even with no coordinates."""
    out = tmp_path / "fesod.gb"
    structure_to_genbank.main(["-i", str(shared_datadir / "foldseek" / "FeSOD"),
                               "-o", str(out), "--cpu", "2"])
    records = _read(out)
    assert len(records) == 20
    assert all(record.annotations["molecule_type"] == "protein" for record in records)


def test_downstream_domainate_annotates_converted_structures(shared_datadir, tmp_path):
    """The point of the converter: structures become ordinary domainator inputs.

    Converts the reference structure to a protein FASTA, then annotates the converted
    input chains with it through domainate.py's phmmer path -- no structural alignment
    involved -- and checks the expected homologs are found.
    """
    chains = tmp_path / "chains.gb"
    structure_to_genbank.main(["-i", str(shared_datadir / "structures" / "inputs"),
                               "-o", str(chains), "--cpu", "2"])
    reference_gb = tmp_path / "ref.gb"
    structure_to_genbank.main(["-i", str(shared_datadir / "structures" / "refs"),
                               "-o", str(reference_gb), "--cpu", "2"])

    reference_fasta = tmp_path / "ref.fasta"
    with open(reference_fasta, "w") as handle:
        for record in _read(reference_gb):
            handle.write(f">{record.id}\n{str(record.seq).upper()}\n")

    annotated = tmp_path / "annotated.gb"
    domainate.main(["-i", str(chains), "-r", str(reference_fasta),
                    "-o", str(annotated), "-e", "0.001", "--cpu", "2"])

    hits = {record.id for record in _read(annotated)
            if any(f.type == "Domainator" for f in record.features)}
    # ubiquitin and its homolog NEDD8, and nothing else
    assert hits == {"1UBQ_A", "1NDD_A"}
