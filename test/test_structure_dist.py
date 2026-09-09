"""Tests for structure_dist.py."""
import shutil

import numpy as np
import pytest
import scipy.sparse
import scipy.sparse.csgraph

from domainator import structure_dist, structure_to_genbank, transform_matrix, utils
from domainator.data_matrix import DataMatrix


needs_foldseek = pytest.mark.skipif(shutil.which("foldseek") is None,
                                    reason="foldseek is not installed")
needs_reseek = pytest.mark.skipif(shutil.which("reseek") is None,
                                  reason="reseek is not installed")
needs_both = pytest.mark.skipif(shutil.which("foldseek") is None or shutil.which("reseek") is None,
                                reason="foldseek and reseek are not both installed")


# Every chain foldseek finds in test/data/structures/inputs, in name order.
ALL_NAMES = ["1CRN_A", "1IGD_A", "1NDD_A", "1UBQ_A",
             "1ZNI_A", "1ZNI_B", "1ZNI_C", "1ZNI_D", "some.long.name_A"]
# reseek drops chains shorter than 32 residues, so the four insulin chains are absent.
RESEEK_NAMES = {"1CRN_A", "1IGD_A", "1NDD_A", "1UBQ_A", "some.long.name_A"}

CHAIN_LENGTHS = {"1UBQ_A": 76, "1NDD_A": 74, "1IGD_A": 61, "1CRN_A": 46,
                 "1ZNI_A": 21, "1ZNI_B": 30, "1ZNI_C": 21, "1ZNI_D": 30,
                 "some.long.name_A": 61}


@pytest.fixture
def inputs(shared_datadir):
    return str(shared_datadir / "structures" / "inputs")


@pytest.fixture
def references(shared_datadir):
    return str(shared_datadir / "structures" / "refs")


def _cell(matrix, row, col):
    array = matrix.toarray()
    return array[matrix.row_to_idx[row], matrix.column_to_idx[col]]


# --------------------------------------------------------------------------------------
# Validation and pure logic: no external binary required.
# --------------------------------------------------------------------------------------

def test_mode_sets_do_not_drift_from_transform_matrix():
    # If transform_matrix ever grows a mode, structure_dist should offer it too rather
    # than silently falling behind.
    assert structure_dist.BITSCORE_MODES == frozenset(transform_matrix.MODES)
    assert structure_dist.ALL_MODES >= structure_dist.BITSCORE_MODES
    assert structure_dist.MST_KNN_STREAMABLE_MODES <= structure_dist.ALL_MODES
    assert structure_dist.DIST_MODES <= structure_dist.ALL_MODES


@pytest.mark.parametrize("mode", sorted(structure_dist.DIST_MODES))
def test_sparse_rejected_for_every_distance_mode(inputs, tmp_path, mode):
    with pytest.raises(ValueError, match="Sparse distance matrices not implemented"):
        structure_dist.main(["-i", inputs, "--mode", mode,
                             "--sparse", str(tmp_path / "out.hdf5")])


def test_symmetrize_requires_a_square_comparison(inputs, references, tmp_path):
    with pytest.raises(ValueError, match="square, symmetric comparison"):
        structure_dist.main(["-i", inputs, "-r", references, "--symmetrize", "max",
                             "--dense_text", str(tmp_path / "out.tsv")])


def test_symmetrize_rejected_with_streaming_pruning(inputs, tmp_path):
    with pytest.raises(ValueError, match="always max-symmetric"):
        structure_dist.main(["-i", inputs, "--symmetrize", "min", "--knn", "2",
                             "--sparse", str(tmp_path / "out.hdf5")])


def test_streaming_pruning_rejects_a_non_streamable_mode(inputs, tmp_path):
    with pytest.raises(ValueError, match="only supported with --mode"):
        structure_dist.main(["-i", inputs, "--mode", "norm_score", "--mst_knn", "2",
                             "--sparse", str(tmp_path / "out.hdf5")])


def test_streaming_pruning_requires_a_square_comparison(inputs, references, tmp_path):
    with pytest.raises(ValueError, match="square, symmetric comparison"):
        structure_dist.main(["-i", inputs, "-r", references, "--knn", "2",
                             "--sparse", str(tmp_path / "out.hdf5")])


def test_metrics_is_rejected_because_mode_derives_it(inputs, tmp_path):
    with pytest.raises(ValueError, match="derives the structural metric from --mode"):
        structure_dist.main(["-i", inputs, "--metrics", "tmscore",
                             "--dense_text", str(tmp_path / "out.tsv")])


def test_no_output_specified(inputs):
    with pytest.raises(ValueError, match="No output specified"):
        structure_dist.main(["-i", inputs])


def test_dense_output_needs_an_hdf5_extension(inputs, tmp_path):
    with pytest.raises(ValueError, match="hdf5 related extension"):
        structure_dist.main(["-i", inputs, "--dense", str(tmp_path / "out.txt")])


def test_structural_mode_bounds_the_lower_bound(inputs, tmp_path):
    with pytest.raises(ValueError, match="--lb must be between 0 and 1"):
        structure_dist.main(["-i", inputs, "--mode", "tmscore", "--lb", "5",
                             "--dense_text", str(tmp_path / "out.tsv")])


def test_symmetrize_max_min_and_mean_on_a_one_directional_pair():
    # A pair the aligner reported in only one direction: 0 in the matrix means "no
    # alignment reported", so min deletes it and mean halves it. Only max keeps it whole.
    matrix = scipy.sparse.csr_array(np.array([[0.0, 8.0, 0.0],
                                              [0.0, 0.0, 0.0],
                                              [0.0, 0.0, 4.0]]))
    assert structure_dist.symmetrize_matrix(matrix, "max").toarray()[0, 1] == 8.0
    assert structure_dist.symmetrize_matrix(matrix, "max").toarray()[1, 0] == 8.0
    assert structure_dist.symmetrize_matrix(matrix, "min").toarray()[0, 1] == 0.0
    assert structure_dist.symmetrize_matrix(matrix, "mean").toarray()[0, 1] == 4.0
    # a symmetric cell is untouched by all three
    for how in ("max", "min", "mean"):
        assert structure_dist.symmetrize_matrix(matrix, how).toarray()[2, 2] == 4.0
    # "none" is the identity
    assert structure_dist.symmetrize_matrix(matrix, "none") is matrix


def test_symmetrize_works_on_a_dense_matrix():
    dense = np.array([[0.0, 8.0], [2.0, 0.0]])
    assert structure_dist.symmetrize_matrix(dense, "max").tolist() == [[0.0, 8.0], [8.0, 0.0]]
    assert structure_dist.symmetrize_matrix(dense, "min").tolist() == [[0.0, 2.0], [2.0, 0.0]]
    assert structure_dist.symmetrize_matrix(dense, "mean").tolist() == [[0.0, 5.0], [5.0, 0.0]]


def test_structural_distance_treats_an_absent_pair_as_maximally_distant():
    matrix = scipy.sparse.csr_array(np.array([[0.9, 0.0], [0.0, 0.9]]))
    out = structure_dist._apply_mode(matrix, "tmscore_dist", None, None,
                                     self_comparison=True)
    # absent pair -> maximum distance, self -> zero even though the stored value was 0.9
    assert out[0, 1] == pytest.approx(1.0)
    assert out[0, 0] == pytest.approx(0.0)


def test_metric_reader_clamps_quietly_within_tolerance_and_loudly_outside():
    # foldseek's alntmscore overshoots 1.0 slightly on a perfect self-alignment, which is
    # routine; a wildly out-of-range value is not.
    class _Hit:
        def __init__(self, tmscore):
            self.tmscore = tmscore

    reader = structure_dist._MetricValueReader("tmscore")
    assert reader.value(_Hit(1.05)) == pytest.approx(1.0)
    assert reader.far_out_of_range == 0
    assert reader.value(_Hit(7.0)) == pytest.approx(1.0)
    assert reader.far_out_of_range == 1
    assert reader.value(_Hit(None)) is None
    assert reader.missing == 1


# --------------------------------------------------------------------------------------
# foldseek
# --------------------------------------------------------------------------------------

@needs_foldseek
def test_self_comparison_matrix(inputs, tmp_path):
    out = tmp_path / "d.hdf5"
    structure_dist.main(["-i", inputs, "--dense", str(out), "--cpu", "2"])
    matrix = DataMatrix.from_file(str(out))

    assert sorted(matrix.rows) == ALL_NAMES
    assert matrix.rows == matrix.columns
    # build_ssn and build_ssn_viewer both require this
    assert matrix.symmetric_labels is True
    assert matrix.data_type == "score"
    assert matrix.shape == (9, 9)

    array = matrix.toarray()
    # every chain aligns to itself
    assert (np.diag(array) > 0).all()
    # ubiquitin and NEDD8 are homologs; crambin is an unrelated fold
    assert _cell(matrix, "1UBQ_A", "1NDD_A") > 0
    assert _cell(matrix, "1UBQ_A", "1CRN_A") == 0


@needs_foldseek
def test_row_lengths_are_the_chain_lengths(inputs, tmp_path):
    out = tmp_path / "d.hdf5"
    structure_dist.main(["-i", inputs, "--dense", str(out), "--cpu", "2"])
    matrix = DataMatrix.from_file(str(out))
    lengths = dict(zip(matrix.rows, matrix.row_lengths))
    assert lengths == CHAIN_LENGTHS


@needs_foldseek
def test_labels_match_structure_to_genbank_record_ids(inputs, tmp_path):
    # The load-bearing interoperability property: a matrix and a genbank of the same
    # structures must use the same names, or nothing downstream lines up.
    matrix_path = tmp_path / "d.hdf5"
    genbank_path = tmp_path / "chains.gb"
    structure_dist.main(["-i", inputs, "--dense", str(matrix_path), "--cpu", "2"])
    structure_to_genbank.main(["-i", inputs, "-o", str(genbank_path), "--cpu", "2"])

    matrix = DataMatrix.from_file(str(matrix_path))
    record_ids = [record.id for record in
                  utils.parse_seqfiles((str(genbank_path),), default_molecule_type="protein")]
    assert sorted(matrix.rows) == sorted(record_ids)


@needs_foldseek
def test_the_raw_matrix_is_asymmetric_and_symmetrize_fixes_it(inputs, tmp_path):
    raw_path = tmp_path / "raw.hdf5"
    sym_path = tmp_path / "sym.hdf5"
    structure_dist.main(["-i", inputs, "--dense", str(raw_path), "--cpu", "2"])
    structure_dist.main(["-i", inputs, "--symmetrize", "max", "--dense", str(sym_path),
                         "--cpu", "2"])
    raw = DataMatrix.from_file(str(raw_path)).toarray()
    sym = DataMatrix.from_file(str(sym_path)).toarray()

    # foldseek scores each direction separately, so the raw matrix really is asymmetric
    assert not np.allclose(raw, raw.T)
    assert np.allclose(sym, sym.T)
    # max never loses a value
    assert (sym >= raw).all()


@needs_foldseek
def test_symmetrize_min_keeps_a_subset_of_max(inputs, tmp_path):
    paths = dict()
    for how in ("min", "max"):
        path = tmp_path / f"{how}.hdf5"
        structure_dist.main(["-i", inputs, "--symmetrize", how, "--dense", str(path),
                             "--cpu", "2"])
        paths[how] = DataMatrix.from_file(str(path)).toarray()
    assert np.allclose(paths["min"], paths["min"].T)
    # a pair found in only one direction is dropped by min but kept by max
    assert ((paths["min"] != 0) <= (paths["max"] != 0)).all()


@needs_foldseek
def test_bool_mode(inputs, tmp_path):
    score_path = tmp_path / "score.hdf5"
    bool_path = tmp_path / "bool.hdf5"
    structure_dist.main(["-i", inputs, "--dense", str(score_path), "--cpu", "2"])
    structure_dist.main(["-i", inputs, "--mode", "bool", "--dense", str(bool_path),
                         "--cpu", "2"])
    score = DataMatrix.from_file(str(score_path)).toarray()
    booled = DataMatrix.from_file(str(bool_path))
    assert booled.data_type == "bool"
    assert set(np.unique(booled.toarray())) <= {0, 1}
    assert ((booled.toarray() != 0) == (score != 0)).all()


@needs_foldseek
def test_fident_mode_is_a_fraction(inputs, tmp_path):
    out = tmp_path / "f.hdf5"
    structure_dist.main(["-i", inputs, "--mode", "fident", "--dense", str(out), "--cpu", "2"])
    matrix = DataMatrix.from_file(str(out))
    assert matrix.data_type == "fident"
    array = matrix.toarray()
    assert array.min() >= 0 and array.max() <= 1
    assert _cell(matrix, "1UBQ_A", "1UBQ_A") == pytest.approx(1.0)


@needs_foldseek
def test_tmscore_and_lddt_modes(inputs, tmp_path):
    tm_path = tmp_path / "tm.hdf5"
    lddt_path = tmp_path / "lddt.hdf5"
    structure_dist.main(["-i", inputs, "--mode", "tmscore", "--alignment_type", "1",
                         "-e", "10", "--dense", str(tm_path), "--cpu", "2"])
    structure_dist.main(["-i", inputs, "--mode", "lddt", "-e", "10",
                         "--dense", str(lddt_path), "--cpu", "2"])
    tm = DataMatrix.from_file(str(tm_path))
    lddt = DataMatrix.from_file(str(lddt_path))

    assert tm.data_type == "tmscore"
    # foldseek's alntmscore overshoots 1.0 on a self-alignment; we clamp it
    assert tm.toarray().max() <= 1.0
    assert _cell(tm, "1UBQ_A", "1UBQ_A") == pytest.approx(1.0)
    # ubiquitin and NEDD8 are the same fold
    assert _cell(tm, "1UBQ_A", "1NDD_A") > 0.9
    assert _cell(lddt, "1UBQ_A", "1UBQ_A") == pytest.approx(1.0, abs=0.01)


@needs_foldseek
def test_tmscore_dist_is_one_minus_tmscore_with_a_zero_diagonal(inputs, tmp_path):
    tm_path = tmp_path / "tm.hdf5"
    dist_path = tmp_path / "dist.hdf5"
    common = ["-i", inputs, "--alignment_type", "1", "-e", "10", "--cpu", "2"]
    structure_dist.main(common + ["--mode", "tmscore", "--dense", str(tm_path)])
    structure_dist.main(common + ["--mode", "tmscore_dist", "--dense", str(dist_path)])
    tm = DataMatrix.from_file(str(tm_path)).toarray()
    dist = DataMatrix.from_file(str(dist_path))

    assert dist.data_type == "tmscore_dist"
    array = dist.toarray()
    assert np.allclose(np.diag(array), 0.0)
    off_diagonal = ~np.eye(len(ALL_NAMES), dtype=bool)
    assert np.allclose(array[off_diagonal], (1.0 - tm)[off_diagonal])
    # a pair that never aligned is maximally distant, not zero
    assert _cell(dist, "1UBQ_A", "1CRN_A") == pytest.approx(1.0)


@needs_foldseek
def test_k_keeps_one_entry_per_row(inputs, tmp_path):
    out = tmp_path / "k.hdf5"
    structure_dist.main(["-i", inputs, "-k", "1", "--dense", str(out), "--cpu", "2"])
    array = DataMatrix.from_file(str(out)).toarray()
    assert ((array != 0).sum(axis=1) == 1).all()


@needs_foldseek
def test_lower_bound_drops_weak_pairs(inputs, tmp_path):
    unbounded = tmp_path / "a.hdf5"
    bounded = tmp_path / "b.hdf5"
    structure_dist.main(["-i", inputs, "--dense", str(unbounded), "--cpu", "2"])
    structure_dist.main(["-i", inputs, "--lb", "400", "--dense", str(bounded), "--cpu", "2"])
    loose = DataMatrix.from_file(str(unbounded)).toarray()
    tight = DataMatrix.from_file(str(bounded)).toarray()
    assert (tight[tight != 0] > 400).all()
    assert (tight != 0).sum() < (loose != 0).sum()


@needs_foldseek
def test_rectangular_comparison(inputs, references, tmp_path):
    out = tmp_path / "r.hdf5"
    structure_dist.main(["-i", inputs, "-r", references, "--dense", str(out), "--cpu", "2"])
    matrix = DataMatrix.from_file(str(out))
    assert sorted(matrix.rows) == ALL_NAMES
    assert matrix.columns == ["UBQref_A"]
    assert matrix.shape == (9, 1)
    assert matrix.symmetric_labels is False
    assert _cell(matrix, "1UBQ_A", "UBQref_A") > 0
    assert _cell(matrix, "1CRN_A", "UBQref_A") == 0


@needs_foldseek
def test_reference_identical_to_input_is_treated_as_a_self_comparison(inputs, tmp_path):
    out = tmp_path / "s.hdf5"
    structure_dist.main(["-i", inputs, "-r", inputs, "--dense", str(out), "--cpu", "2"])
    matrix = DataMatrix.from_file(str(out))
    assert matrix.symmetric_labels is True


@needs_foldseek
def test_mst_knn_streaming_matches_the_batch_transform(inputs, tmp_path):
    full = tmp_path / "full.hdf5"
    batch = tmp_path / "batch.hdf5"
    stream = tmp_path / "stream.hdf5"
    structure_dist.main(["-i", inputs, "--sparse", str(full), "--cpu", "2"])
    transform_matrix.main(["-i", str(full), "--mst_knn", "2", "--sparse", str(batch)])
    structure_dist.main(["-i", inputs, "--mst_knn", "2", "--sparse", str(stream),
                         "--cpu", "2"])

    batch_matrix = DataMatrix.from_file(str(batch))
    stream_matrix = DataMatrix.from_file(str(stream))
    assert batch_matrix.rows == stream_matrix.rows

    # Compare by connected-components partition, which is invariant to which of several
    # equally-valid spanning forests each path happened to pick.
    batch_labels = scipy.sparse.csgraph.connected_components(
        batch_matrix.toarray(), directed=False)[1]
    stream_labels = scipy.sparse.csgraph.connected_components(
        stream_matrix.toarray(), directed=False)[1]
    assert (batch_labels == stream_labels).all()


@needs_foldseek
def test_knn_output_is_max_symmetric(inputs, tmp_path):
    out = tmp_path / "knn.hdf5"
    structure_dist.main(["-i", inputs, "--knn", "2", "--sparse", str(out), "--cpu", "2"])
    array = DataMatrix.from_file(str(out)).toarray()
    assert np.allclose(array, array.T)


@needs_foldseek
def test_mst_knn_with_a_structural_mode(inputs, tmp_path):
    out = tmp_path / "tm_knn.hdf5"
    structure_dist.main(["-i", inputs, "--mode", "tmscore", "--alignment_type", "1",
                         "-e", "10", "--mst_knn", "2", "--sparse", str(out), "--cpu", "2"])
    matrix = DataMatrix.from_file(str(out))
    assert matrix.data_type == "tmscore"
    array = matrix.toarray()
    assert np.allclose(array, array.T)
    assert array.max() <= 1.0


@needs_foldseek
def test_default_run_does_not_warn_about_saturation(inputs, tmp_path, recwarn):
    # Regression guard: foldseek's own default caps hits at 1000 per query, which would
    # silently truncate a large all-vs-all. structure_dist derives its own cap instead.
    structure_dist.main(["-i", inputs, "--dense", str(tmp_path / "d.hdf5"), "--cpu", "2"])
    assert not [w for w in recwarn if "max_seqs" in str(w.message)]


@needs_foldseek
def test_max_seqs_saturation_warns(inputs, tmp_path):
    with pytest.warns(RuntimeWarning, match="was reached by"):
        structure_dist.main(["-i", inputs, "--max_seqs", "1", "-e", "10",
                             "--dense", str(tmp_path / "d.hdf5"), "--cpu", "2"])


@needs_foldseek
def test_max_output_gb(inputs, tmp_path):
    out = tmp_path / "d.hdf5"
    with pytest.raises(SystemExit, match="--max_output_gb"):
        structure_dist.main(["-i", inputs, "--dense", str(out), "--cpu", "2",
                             "--max_output_gb", "0.000000001"])
    assert not out.exists()


@needs_foldseek
def test_max_output_gb_zero_disables_the_guardrail(inputs, tmp_path):
    out = tmp_path / "d.hdf5"
    structure_dist.main(["-i", inputs, "--dense", str(out), "--cpu", "2",
                         "--max_output_gb", "0"])
    assert out.exists()


@needs_foldseek
def test_prebuilt_database_input_matches_structure_file_input(inputs, tmp_path):
    from_files = tmp_path / "files.hdf5"
    from_db = tmp_path / "db.hdf5"
    db_prefix = tmp_path / "mydb"
    structure_dist.main(["-i", inputs, "--dense", str(from_files),
                         "--keep_db", str(db_prefix), "--cpu", "2"])
    structure_dist.main(["-i", str(db_prefix), "--dense", str(from_db), "--cpu", "2"])

    files_matrix = DataMatrix.from_file(str(from_files))
    db_matrix = DataMatrix.from_file(str(from_db))
    assert files_matrix.rows == db_matrix.rows
    assert np.allclose(files_matrix.toarray(), db_matrix.toarray())


@needs_foldseek
def test_hits_tsv(inputs, tmp_path):
    hits = tmp_path / "hits.tsv"
    structure_dist.main(["-i", inputs, "--dense", str(tmp_path / "d.hdf5"),
                         "--hits_tsv", str(hits), "--cpu", "2"])
    assert hits.is_file()
    queries = {line.split("\t")[0] for line in hits.read_text().splitlines() if line.strip()}
    assert queries <= set(ALL_NAMES)


@needs_foldseek
def test_matrix_feeds_the_downstream_tools(inputs, tmp_path):
    from domainator import build_ssn, build_tree

    score_path = tmp_path / "score.hdf5"
    dist_path = tmp_path / "dist.hdf5"
    structure_dist.main(["-i", inputs, "--sparse", str(score_path), "--cpu", "2"])
    structure_dist.main(["-i", inputs, "--mode", "tmscore_dist", "--alignment_type", "1",
                         "-e", "10", "--dense", str(dist_path), "--cpu", "2"])

    xgmml = tmp_path / "ssn.xgmml"
    build_ssn.main(["-i", str(score_path), "--xgmml", str(xgmml)])
    assert xgmml.is_file() and "<edge" in xgmml.read_text()

    # A data_type of "tmscore_dist" must not be rejected anywhere downstream.
    newick = tmp_path / "tree.nwk"
    build_tree.main(["-i", str(dist_path), "--newick", str(newick)])
    tree = newick.read_text()
    assert tree.strip() and "1UBQ_A" in tree


@needs_foldseek
def test_transform_matrix_accepts_a_structural_matrix_for_pruning_only(inputs, tmp_path):
    tm_path = tmp_path / "tm.hdf5"
    structure_dist.main(["-i", inputs, "--mode", "tmscore", "--alignment_type", "1",
                         "-e", "10", "--sparse", str(tm_path), "--cpu", "2"])
    # pass-through pruning is fine
    transform_matrix.main(["-i", str(tm_path), "--knn", "2",
                           "--sparse", str(tmp_path / "pruned.hdf5")])
    # but re-deriving a score transform from it is not
    with pytest.raises(ValueError, match="only supported from a 'score' matrix"):
        transform_matrix.main(["-i", str(tm_path), "--mode", "efi_score",
                               "--dense", str(tmp_path / "bad.hdf5")])


# --------------------------------------------------------------------------------------
# reseek
# --------------------------------------------------------------------------------------

@needs_reseek
def test_reseek_self_comparison(inputs, tmp_path):
    out = tmp_path / "r.hdf5"
    with pytest.warns(RuntimeWarning, match="shorter than 32 residues"):
        structure_dist.main(["-i", inputs, "--algorithm", "reseek",
                             "--dense", str(out), "--cpu", "2"])
    matrix = DataMatrix.from_file(str(out))
    assert set(matrix.rows) == RESEEK_NAMES
    assert (np.diag(matrix.toarray()) > 0).all()
    assert _cell(matrix, "1UBQ_A", "1NDD_A") > 0


@needs_reseek
def test_reseek_supports_fident(inputs, tmp_path):
    out = tmp_path / "f.hdf5"
    with pytest.warns(RuntimeWarning, match="shorter than 32 residues"):
        structure_dist.main(["-i", inputs, "--algorithm", "reseek", "--mode", "fident",
                             "--dense", str(out), "--cpu", "2"])
    matrix = DataMatrix.from_file(str(out))
    # reseek reports pctid, which structure_lib converts to a fraction
    assert matrix.toarray().max() <= 1.0
    assert _cell(matrix, "1UBQ_A", "1UBQ_A") == pytest.approx(1.0)


@needs_reseek
def test_reseek_rejects_tmscore_and_efi_score(inputs, tmp_path):
    with pytest.raises(RuntimeError, match="tmscore"):
        structure_dist.main(["-i", inputs, "--algorithm", "reseek", "--mode", "tmscore",
                             "--dense", str(tmp_path / "a.hdf5")])
    with pytest.raises(RuntimeError, match="bit score"):
        structure_dist.main(["-i", inputs, "--algorithm", "reseek", "--mode", "efi_score",
                             "--dense", str(tmp_path / "b.hdf5")])


@needs_both
def test_backends_agree_on_which_pairs_are_related(inputs, tmp_path):
    foldseek_path = tmp_path / "fs.hdf5"
    reseek_path = tmp_path / "rs.hdf5"
    structure_dist.main(["-i", inputs, "--mode", "bool", "--dense", str(foldseek_path),
                         "--cpu", "2"])
    with pytest.warns(RuntimeWarning, match="shorter than 32 residues"):
        structure_dist.main(["-i", inputs, "--algorithm", "reseek", "--mode", "bool",
                             "--dense", str(reseek_path), "--cpu", "2"])

    foldseek = DataMatrix.from_file(str(foldseek_path))
    reseek = DataMatrix.from_file(str(reseek_path))
    # compare only the chains reseek keeps, and only the relationships that matter
    for a, b in [("1UBQ_A", "1NDD_A"), ("1IGD_A", "some.long.name_A"),
                 ("1UBQ_A", "1CRN_A")]:
        assert bool(_cell(foldseek, a, b)) == bool(_cell(reseek, a, b))
