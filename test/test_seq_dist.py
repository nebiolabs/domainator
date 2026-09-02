import warnings
warnings.filterwarnings("ignore", module='numpy')
from domainator import seq_dist
import tempfile
import pandas as pd
import numpy as np
import scipy
import scipy.sparse
import pytest
from domainator.data_matrix import DataMatrix
from helpers import compare_iterables


def test_diamond_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        results = seq_dist.diamond( str(shared_datadir / "FeSOD_20.fasta"), str(shared_datadir / "FeSOD_20.fasta"), max_target_seqs=1, threads=8, tmpdir=output_dir, mode="us", max_hsps=1)
        assert len(list(results)) == 20
        for score,q,s in results:
            assert q == s
            assert score > 100

def test_seq_dist_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"
        seq_dist.main(["-i", str(shared_datadir / "FeSOD_20.fasta"), "-r", str(shared_datadir / "FeSOD_20.fasta"), "--dense", out, "-k", "1", "--mode", "score"])
        dense_matrix = DataMatrix.from_file(out)
        out_mat = dense_matrix.toarray()
        assert (out_mat > 0).sum() == 20
        assert (out_mat.diagonal() == [410., 429., 405., 425., 411., 410., 400., 434., 426., 414., 398., 424., 430., 413., 402., 417., 451., 419., 429., 422.]).all()
        assert dense_matrix.data_type == "score"
        

def test_seq_dist_2(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"
        seq_dist.main(["-i", str(shared_datadir / "FeSOD_20.fasta"), "-r", str(shared_datadir / "FeSOD_20.fasta"), "--dense", out, "-k", "3","--mode","bool" ])

        dense_matrix = DataMatrix.from_file(out)
        out_mat = dense_matrix.toarray()
        assert (out_mat == 1).sum() == 60
        assert ((out_mat).sum(axis=1) == 3).all()
        assert dense_matrix.data_type == "bool"

def test_seq_dist_mst_knn_matches_batch(shared_datadir):
    from domainator import transform_matrix
    from scipy.sparse.csgraph import connected_components

    fasta = str(shared_datadir / "FeSOD_20.fasta")
    with tempfile.TemporaryDirectory() as output_dir:
        full = output_dir + "/full.hdf5"
        batch = output_dir + "/batch.hdf5"
        stream = output_dir + "/stream.hdf5"

        # batch reference: full matrix, then post-hoc mst_knn
        seq_dist.main(["-i", fasta, "-r", fasta, "--sparse", full, "--mode", "score"])
        transform_matrix.main(["-i", full, "--sparse", batch, "--mst_knn", "3"])

        # streaming path
        seq_dist.main(["-i", fasta, "-r", fasta, "--sparse", stream, "--mode", "score", "--mst_knn", "3"])

        batch_dm = DataMatrix.from_file(batch)
        stream_dm = DataMatrix.from_file(stream)
        b = batch_dm.toarray()
        s = stream_dm.toarray()
        assert batch_dm.rows == stream_dm.rows
        # diagonal (self-scores) preserved by the streaming path
        np.testing.assert_array_equal(b.diagonal(), s.diagonal())
        # tie-invariant equivalence: identical connected-components partition
        _, bl = connected_components(b > 0, directed=False)
        _, sl = connected_components(s > 0, directed=False)
        np.testing.assert_array_equal(bl[:, None] == bl[None, :], sl[:, None] == sl[None, :])


def test_seq_dist_mst_knn_rejects_non_streamable_mode(shared_datadir):
    fasta = str(shared_datadir / "FeSOD_20.fasta")
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"
        with pytest.raises(ValueError):
            seq_dist.main(["-i", fasta, "-r", fasta, "--sparse", out, "--mode", "norm_score", "--mst_knn", "2"])


def test_seq_dist_mst_knn_requires_symmetric_reference(shared_datadir):
    fasta = str(shared_datadir / "FeSOD_20.fasta")
    ref = str(shared_datadir / "pdonr_peptides.fasta")
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"
        with pytest.raises(ValueError):
            seq_dist.main(["-i", fasta, "-r", ref, "--sparse", out, "--mode", "score", "--mst_knn", "2"])


def test_seq_dist_mst_knn_accepts_k_one(shared_datadir):
    # k == 1 reduces to the MST plus each node's single nearest neighbor.
    fasta = str(shared_datadir / "FeSOD_20.fasta")
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"
        seq_dist.main(["-i", fasta, "-r", fasta, "--sparse", out, "--mode", "score", "--mst_knn", "1"])


def test_seq_dist_mst_knn_k_zero_is_mst_only(shared_datadir):
    # k == 0 keeps only the maximum spanning forest: at most n - 1 off-diagonal edges.
    fasta = str(shared_datadir / "FeSOD_20.fasta")
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"
        seq_dist.main(["-i", fasta, "-r", fasta, "--sparse", out, "--mode", "score", "--mst_knn", "0"])
        matrix = DataMatrix.from_file(out)
        array = matrix.toarray()
        n = len(matrix.rows)
        off_diagonal = np.count_nonzero(array) - np.count_nonzero(array.diagonal())
        assert 0 < off_diagonal <= 2 * (n - 1)


def test_seq_dist_mst_knn_rejects_negative_k(shared_datadir):
    fasta = str(shared_datadir / "FeSOD_20.fasta")
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"
        with pytest.raises(SystemExit):
            seq_dist.main(["-i", fasta, "-r", fasta, "--sparse", out, "--mode", "score", "--mst_knn", "-1"])


def test_seq_hmmsearch_1(shared_datadir):
    expected = np.array([
     [0.0,0.0,0.0,0.0,0.0,0.0,0.0],
     [0.0,0.0,0.0,0.0,103.15,0.0,0.0],
     [0.0,0.0,0.0,50.15,0.0,0.0,0.0],
     [16.27,0.0,329.9,0.0,0.0,12.7,0.0],
     [0.0,90.6,0.0,0.0,0.0,0.0,10.96],
    ])
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "tmp_out"
        out = output_dir + "/out.hdf5"
        seq_dist.main(["-i", str(shared_datadir / "pdonr_peptides.fasta"), "-r", str(shared_datadir / "pdonr_hmms.hmm"), "--dense", out, "--algorithm", "hmmer","--mode", "score" ])

        out_table = DataMatrix.from_file(out)
        compare_iterables(out_table.rows, ["pDONR201_1", "pDONR201_2", "pDONR201_3","pDONR201_4","pDONR201_5"])
        compare_iterables(out_table.columns,  ["2-oxoacid_dh", "APH","CAT","CcdA","CcdB","Condensation","TCAD9"])

        compare_iterables(out_table.column_lengths , [233, 239, 204, 71, 100, 457, 437])
        compare_iterables(out_table.row_lengths , [35, 102, 42, 220, 254])
        assert out_table.data_type == "score"
        out_mat = out_table.toarray()
        assert (out_mat == expected).all()

def test_seq_hmmsearch_sparse_1(shared_datadir):
    expected = np.array([103.15, 50.15, 16.27, 329.9, 12.7, 90.6, 10.96])
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/out.hdf5"
        seq_dist.main(["-i", str(shared_datadir / "pdonr_peptides.fasta"), "-r", str(shared_datadir / "pdonr_hmms.hmm"), "--sparse", out, "--algorithm", "hmmer","--mode", "score", ])

        out_table = DataMatrix.from_file(out)
        compare_iterables(out_table.rows, ["pDONR201_1", "pDONR201_2", "pDONR201_3","pDONR201_4","pDONR201_5"])
        compare_iterables(out_table.columns, ["2-oxoacid_dh", "APH","CAT","CcdA","CcdB","Condensation","TCAD9"])
        assert out_table.data_type == "score"
        compare_iterables(out_table.column_lengths , [233, 239, 204, 71, 100, 457, 437])
        compare_iterables(out_table.row_lengths , [35, 102, 42, 220, 254])
        out_data = list([x[2] for x in out_table.iter_data()])

        assert (out_data == expected).all()


def test_seq_dist_max_output_gb_blocks_dense_output(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"

        with pytest.raises(SystemExit, match="--max_output_gb"):
            seq_dist.main([
                "-i", str(shared_datadir / "FeSOD_20.fasta"),
                "-r", str(shared_datadir / "FeSOD_20.fasta"),
                "--dense", out,
                "-k", "1",
                "--mode", "score",
                "--max_output_gb", "0.000001",
            ])


def test_seq_dist_max_output_gb_zero_disables_guardrail(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"

        seq_dist.main([
            "-i", str(shared_datadir / "FeSOD_20.fasta"),
            "-r", str(shared_datadir / "FeSOD_20.fasta"),
            "--dense", out,
            "-k", "1",
            "--mode", "score",
            "--max_output_gb", "0",
        ])

        assert DataMatrix.from_file(out).shape == (20, 20)


def test_seq_dist_progress_flag(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"

        seq_dist.main([
            "-i", str(shared_datadir / "FeSOD_20.fasta"),
            "-r", str(shared_datadir / "FeSOD_20.fasta"),
            "--dense", out,
            "-k", "1",
            "--mode", "score",
            "--progress",
        ])

        dense_matrix = DataMatrix.from_file(out)
        assert dense_matrix.shape == (20, 20)
        assert (dense_matrix.toarray() > 0).sum() == 20


def test_seq_hmmsearch_progress_flag(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"

        seq_dist.main([
            "-i", str(shared_datadir / "pdonr_peptides.fasta"),
            "-r", str(shared_datadir / "pdonr_hmms.hmm"),
            "--dense", out,
            "--algorithm", "hmmer",
            "--mode", "score",
            "--progress",
        ])

        assert DataMatrix.from_file(out).shape == (5, 7)


def test_seq_hmmsearch_bool_lb_filters_raw_scores(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"

        seq_dist.main([
            "-i", str(shared_datadir / "pdonr_peptides.fasta"),
            "-r", str(shared_datadir / "pdonr_hmms.hmm"),
            "--sparse", out,
            "--algorithm", "hmmer",
            "--mode", "bool",
            "--lb", "100",
        ])

        out_table = DataMatrix.from_file(out)
        out_data = list(x[2] for x in out_table.iter_data())

        assert out_table.data.nnz == 2
        assert out_data == [1.0, 1.0]


class _FakeProgressBar:
    def __init__(self):
        self.updates = []
        self.closed = False

    def update(self, amount):
        self.updates.append(amount)

    def close(self):
        self.closed = True


def test_diamond_query_progress_accounts_for_skipped_queries():
    progress_bar = _FakeProgressBar()
    query_positions = {"q1": 0, "q2": 1, "q3": 2}

    current_query_idx = None
    current_query_idx = seq_dist._advance_query_progress(progress_bar, query_positions, current_query_idx, "q2")
    current_query_idx = seq_dist._advance_query_progress(progress_bar, query_positions, current_query_idx, "q3")
    seq_dist._finish_query_progress(progress_bar, current_query_idx, 3)

    assert progress_bar.updates == [2, 1]
    assert progress_bar.closed is True


def test_diamond_query_progress_advances_on_first_seen_query():
    progress_bar = _FakeProgressBar()
    query_positions = {"q1": 0, "q2": 1, "q3": 2}

    current_query_idx = None
    current_query_idx = seq_dist._advance_query_progress(progress_bar, query_positions, current_query_idx, "q1")
    seq_dist._finish_query_progress(progress_bar, current_query_idx, 3)

    assert progress_bar.updates == [1, 2]
    assert progress_bar.closed is True


def test_sparse_max_result_builder_keeps_best_duplicate_score():
    builder = seq_dist._SparseMaxResultBuilder({"q1": 0}, {"s1": 0, "s2": 1})

    builder.add_result("q1", "s1", 5.0)
    builder.add_result("q1", "s1", 9.0)
    builder.add_result("q1", "s2", 3.0)

    matrix = builder.build()

    assert np.allclose(matrix.toarray(), np.array([[9.0, 3.0]]))


def test_grouped_query_sparse_builder_rejects_noncontiguous_repeated_query():
    builder = seq_dist._GroupedQuerySparseMaxResultBuilder(
        {"q1": 0, "q2": 1},
        {"s1": 0},
        stream_name="test stream",
        expected_query_order=["q1", "q2"],
    )

    builder.add_result("q1", "s1", 1.0)
    builder.add_result("q2", "s1", 2.0)

    with pytest.raises(RuntimeError, match="grouped by query"):
        builder.add_result("q1", "s1", 3.0)


def test_row_norm_score_mode():
    mode="row_norm_score"
    input = np.array([
     [-1,2,1,3,0],   
     [1,2,3,4,5],
     [5,4,3,2,1],
     [1,3,3,3,2],
     [1,1,1,3,2],
    ])
    expected = np.array([
     [-1/3.,2/3.,1/3.,1,0],   
     [1/5.,2/5.,3/5.,4/5.,1],
     [1,4/5.,3/5.,2/5.,1/5.],
     [1/3.,1,1,1,2/3.],
     [1/3.,1/3.,1/3.,1,2/3.],
    ])

    assert (seq_dist.MODES[mode](input) == expected).all()

def test_row_norm_score_mode_sparse():
    mode="row_norm_score"
    input = np.array([
     [-1,2,1,3,0],   
     [1,2,3,4,5],
     [5,4,3,2,1],
     [1,3,3,3,2],
     [1,1,1,3,2],
    ])
    input = scipy.sparse.dok_array(input, dtype=np.float64)

    expected = np.array([
     [-1/3.,2/3.,1/3.,1,0],   
     [1/5.,2/5.,3/5.,4/5.,1],
     [1,4/5.,3/5.,2/5.,1/5.],
     [1/3.,1,1,1,2/3.],
     [1/3.,1/3.,1/3.,1,2/3.],
    ])

    #expected = scipy.sparse.dok_array(expected, dtype=np.float64)
    assert (seq_dist.MODES[mode](input).toarray() == expected).all()

def test_norm_score_mode():
    mode="norm_score"
    input = np.array([
     [-1,2,1,3,0],   
     [1,2,3,4,5],
     [5,4,3,2,1],
     [1,3,3,3,2],
     [1,1,1,3,2],
    ])
    expected = np.array([
     [1 + 1/3., 1/3.,   2/3.,  0.,    1.],   
     [4/5.,     1/2.,   0.,    0.,    0.],
     [0.,       0.,     0.,    1/2., 4/5.],
     [2/3.,     0.,     0.,    0.,   1/3.],
     [2/3.,     2/3.,   2/3.,  0.,   1/3.],
    ])
    expected = -1 * (expected - 1)
    
    assert np.allclose(seq_dist.MODES[mode](input), expected)

def test_norm_score_mode_sparse():
    mode="norm_score"
    input = np.array([
     [-1,2,1,3,0],   
     [1,2,3,4,5],
     [5,4,3,2,1],
     [1,3,3,3,2],
     [1,1,1,3,2],
    ])
    input = scipy.sparse.dok_array(input, dtype=np.float64)
    expected = np.array([
     [1 + 1/3., 1/3.,   2/3.,  0.,    1.],   
     [4/5.,     1/2.,   0.,    0.,    0.],
     [0.,       0.,     0.,    1/2., 4/5.],
     [2/3.,     0.,     0.,    0.,   1/3.],
     [2/3.,     2/3.,   2/3.,  0.,   1/3.],
    ])
    expected = -1 * (expected - 1)

    assert np.allclose(seq_dist.MODES[mode](input).toarray(), expected)

def test_score_dist_mode():
    mode="score_dist"
    input = np.array([
     [-1,2,1,3,0],   
     [1,2,3,4,5],
     [5,4,3,2,1],
     [1,3,3,3,2],
     [1,1,1,3,2],
    ])
    expected = np.array([
     [1 + 1/3., 1/3.,   2/3.,  0.,    1.],   
     [4/5.,     1/2.,   0.,    0.,    0.],
     [0.,       0.,     0.,    1/2., 4/5.],
     [2/3.,     0.,     0.,    0.,   1/3.],
     [2/3.,     2/3.,   2/3.,  0.,   1/3.],
    ])
    assert np.allclose(seq_dist.MODES[mode](input), expected)

def test_efi_score_mode():
    mode="efi_score"
    input = 50 * np.array([
     [-1,2,1,3,0],   
     [1,2,3,4,5],
     [5,4,3,2,1],
     [1,3,3,3,2],
     [1,1,1,3,2],
    ])
    input_sparse = scipy.sparse.dok_array(input, dtype=np.float64)
    row_lens = np.array([1, 2, 3, 4, 5])
    col_lens = np.array([7, 8, 9, 10, 11])
    expected = np.array([[-0., 29.19990958, 14.09725727, 44.15449935, -0., ],
        [13.90537175, 28.89887958, 43.89922684, 58.90496914, 73.91507624],
        [73.93527962, 58.82578789, 43.72313559, 28.62587831, 13.53298584],
        [13.60434175, 43.64934937, 43.59819685, 43.55243936, 28.45954689],
        [13.50743174, 13.44943979, 13.39828727, 43.45552935, 28.36263688]]
    )
    
    dense_out = seq_dist.MODES[mode](input, row_lens, col_lens)
    sparse_out = seq_dist.MODES[mode](input_sparse, row_lens, col_lens)
    
    assert np.allclose(sparse_out.toarray(), dense_out)
    assert np.allclose(dense_out, expected)
    #assert np.allclose(seq_dist.MODES[mode](input, row_lens, col_lens), expected)

# def test_score_dist_mode_sparse():
#     mode="score_dist"
#     input = np.array([
#      [-1,2,1,3,0],   
#      [1,2,3,4,5],
#      [5,4,3,2,1],
#      [1,3,3,3,2],
#      [1,1,1,3,2],
#     ])
#     input = scipy.sparse.dok_array(input, dtype=np.float64)
#     expected = np.array([
#      [1 + 1/3., 1/3.,   2/3.,  0.,    1.],   
#      [4/5.,     1/2.,   0.,    0.,    0.],
#      [0.,       0.,     0.,    1/2., 4/5.],
#      [2/3.,     0.,     0.,    0.,   1/3.],
#      [2/3.,     2/3.,   2/3.,  0.,   1/3.],
#     ])
#     # expected = scipy.sparse.dok_array(expected, dtype=np.float64)
#     print(seq_dist.MODES[mode](input).toarray())
#     print(seq_dist.MODES[mode](np.array(input.toarray())) )
#     assert np.allclose(seq_dist.MODES[mode](input).toarray(), expected)

def test_hmm_compare_1(shared_datadir):
    expected = np.array([
        [146.4,5.05,12.46,6.41,5.82,9.62,5.55],
        [5.05,119.66,4.03,5.21,4.48,5.67,27.21],
        [12.46,4.03,143.71,7.34,6.45,13.24,6.15],
        [6.41,5.21,7.34,60.93,5.21,7.31,5.8],
        [5.82,4.48,6.45,5.21,59.86,6.14,5.68],
        [9.62,5.67,13.24,7.31,6.14,267.28,7.1],
        [5.55,27.21,6.15,5.8,5.68,7.1,244.69],
    ])
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "tmp_out"
        out = output_dir + "/out.tsv"
        seq_dist.main(["-i", str(shared_datadir / "pdonr_hmms.hmm"), "-r", str(shared_datadir / "pdonr_hmms.hmm"), "--dense_text", out, "--algorithm", "hmmer_compare","--mode", "score", ])

        out_table = DataMatrix.from_file(out)
        compare_iterables(out_table.rows, ["2-oxoacid_dh", "APH","CAT","CcdA","CcdB","Condensation","TCAD9"])
        compare_iterables(out_table.columns, ["2-oxoacid_dh", "APH","CAT","CcdA","CcdB","Condensation","TCAD9"])
        
        out_mat = out_table.toarray().round(2)
        assert (out_mat == expected).all()




# MODES = {"rank": (lambda x: scipy.stats.rankdata(x, method='min', axis=1)), 
#         "score": (lambda x: x), 
#         "bool": (lambda x: x), 
#         "norm_score": (lambda x: x/np.amax(x, axis=1))}
# def test_rank_mode():
#     mode="rank"
#     input = np.array([
#      [-1,2,1,3,0],   
#      [1,2,3,4,5],
#      [5,4,3,2,1],
#      [1,3,3,3,2],
#      [1,1,1,3,2],
#     ])
#     expected = np.array([
#      [1,4,3,5,2],   
#      [1,2,3,4,5],
#      [5,4,3,2,1],
#      [1,3,3,3,2],
#      [1,1,1,5,4],
#     ])

#     assert (seq_dist.MODES[mode](input) == expected).all()


def test_seq_dist_knn_matches_batch_transform(shared_datadir):
    """Streaming --knn equals a full comparison followed by transform_matrix --knn."""
    from domainator import transform_matrix

    fasta = str(shared_datadir / "FeSOD_20.fasta")
    with tempfile.TemporaryDirectory() as output_dir:
        full = output_dir + "/full.hdf5"
        batch = output_dir + "/batch.hdf5"
        stream = output_dir + "/stream.hdf5"

        seq_dist.main(["-i", fasta, "-r", fasta, "--sparse", full, "--mode", "score"])
        transform_matrix.main(["-i", full, "--sparse", batch, "--knn", "3"])
        seq_dist.main(["-i", fasta, "-r", fasta, "--sparse", stream, "--mode", "score", "--knn", "3"])

        batch_dm = DataMatrix.from_file(batch)
        stream_dm = DataMatrix.from_file(stream)
        assert batch_dm.rows == stream_dm.rows
        np.testing.assert_array_equal(batch_dm.toarray(), stream_dm.toarray())


def test_seq_dist_knn_keeps_fewer_edges_than_mst_knn(shared_datadir):
    """--knn is a subset of --mst_knn: no spanning tree edges are added."""
    fasta = str(shared_datadir / "FeSOD_20.fasta")
    with tempfile.TemporaryDirectory() as output_dir:
        knn_out = output_dir + "/knn.hdf5"
        mst_knn_out = output_dir + "/mst_knn.hdf5"
        seq_dist.main(["-i", fasta, "-r", fasta, "--sparse", knn_out, "--mode", "score", "--knn", "1"])
        seq_dist.main(["-i", fasta, "-r", fasta, "--sparse", mst_knn_out, "--mode", "score", "--mst_knn", "1"])

        def off_diagonal_pairs(path):
            array = DataMatrix.from_file(path).toarray()
            return {(min(i, j), max(i, j)) for i, j in zip(*np.nonzero(array)) if i != j}

        knn_pairs = off_diagonal_pairs(knn_out)
        mst_knn_pairs = off_diagonal_pairs(mst_knn_out)
        assert len(knn_pairs) > 0
        assert knn_pairs <= mst_knn_pairs


def test_seq_dist_knn_rejects_invalid_use(shared_datadir):
    fasta = str(shared_datadir / "FeSOD_20.fasta")
    ref = str(shared_datadir / "pdonr_peptides.fasta")
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"
        with pytest.raises(SystemExit):  # k must be >= 1
            seq_dist.main(["-i", fasta, "-r", fasta, "--sparse", out, "--mode", "score", "--knn", "0"])
        with pytest.raises(SystemExit):  # --knn and --mst_knn are mutually exclusive
            seq_dist.main(["-i", fasta, "-r", fasta, "--sparse", out, "--mode", "score",
                           "--knn", "2", "--mst_knn", "2"])
        with pytest.raises(ValueError):  # same streaming requirements as --mst_knn
            seq_dist.main(["-i", fasta, "-r", fasta, "--sparse", out, "--mode", "norm_score", "--knn", "2"])
        with pytest.raises(ValueError):
            seq_dist.main(["-i", fasta, "-r", ref, "--sparse", out, "--mode", "score", "--knn", "2"])


def test_parse_diamond_params_accepts_flags_and_values():
    parsed = seq_dist.parse_diamond_params('"-c":4, "-b":0.4, "--target-indexed":null')
    assert list(parsed.items()) == [("-c", 4), ("-b", 0.4), ("--target-indexed", None)]

    argv = []
    seq_dist.add_params_to_args_list(argv, parsed)
    assert argv == ["-c", "4", "-b", "0.4", "--target-indexed"]

    assert seq_dist.parse_diamond_params(None) == {}
    assert seq_dist.parse_diamond_params("") == {}


def test_parse_diamond_params_rejects_bad_input():
    for reserved in ('"--outfmt":6', '"-p":4', '"--max-target-seqs":10', '"--ultra-sensitive":null'):
        with pytest.raises(ValueError, match="may not set"):
            seq_dist.parse_diamond_params(reserved)
    with pytest.raises(ValueError, match="must be diamond option flags"):
        seq_dist.parse_diamond_params('"c":4')
    with pytest.raises(ValueError, match="Could not parse"):
        seq_dist.parse_diamond_params("not json at all")


def test_diamond_params_reach_the_command_line(shared_datadir):
    """--params must be appended, and must replace (not duplicate) seq_dist's own defaults."""
    import subprocess

    fasta = str(shared_datadir / "FeSOD_20.fasta")
    captured = []
    real_popen = subprocess.Popen

    class _Stop(Exception):
        pass

    def fake_popen(cmd, *args, **kwargs):
        if isinstance(cmd, list) and "blastp" in cmd:
            captured.append(cmd)
            raise _Stop()
        return real_popen(cmd, *args, **kwargs)

    with tempfile.TemporaryDirectory() as tmpdir:
        subprocess.Popen = fake_popen
        try:
            for diamond_params in (None, {"--algo": 1}, {"--tmpdir": "/scratch/x"}, {"-b": 0.4, "-c": 4}):
                try:
                    list(seq_dist.diamond(fasta, fasta, None, 1, tmpdir, "us", diamond_params=diamond_params))
                except _Stop:
                    pass
        finally:
            subprocess.Popen = real_popen

        default_cmd, algo_cmd, tmpdir_cmd, extra_cmd = captured

        # seq_dist's defaults, and diamond's scratch kept inside the managed tmpdir
        assert default_cmd[default_cmd.index("--algo") + 1] == "0"
        assert default_cmd[default_cmd.index("--tmpdir") + 1] == tmpdir

        # overrides replace the default rather than appending a second copy
        assert algo_cmd.count("--algo") == 1
        assert algo_cmd[algo_cmd.index("--algo") + 1] == "1"
        assert tmpdir_cmd.count("--tmpdir") == 1
        assert tmpdir_cmd[tmpdir_cmd.index("--tmpdir") + 1] == "/scratch/x"

        # unrelated options are passed through
        assert extra_cmd[extra_cmd.index("-b") + 1] == "0.4"
        assert extra_cmd[extra_cmd.index("-c") + 1] == "4"


def test_seq_dist_params_does_not_change_results(shared_datadir):
    """Block/chunk tuning is a resource knob, so it must not alter the matrix."""
    fasta = str(shared_datadir / "FeSOD_20.fasta")
    with tempfile.TemporaryDirectory() as output_dir:
        base = output_dir + "/base.hdf5"
        tuned = output_dir + "/tuned.hdf5"
        seq_dist.main(["-i", fasta, "-r", fasta, "--sparse", base, "--mode", "score", "--mst_knn", "2"])
        seq_dist.main(["-i", fasta, "-r", fasta, "--sparse", tuned, "--mode", "score", "--mst_knn", "2",
                       "--params", '"-c":4, "-b":0.4'])
        base_data = DataMatrix.from_file(base).data.toarray()
        tuned_data = DataMatrix.from_file(tuned).data.toarray()
        assert np.allclose(base_data, tuned_data)


def test_seq_dist_params_rejects_non_diamond_algorithm(shared_datadir):
    fasta = str(shared_datadir / "FeSOD_20.fasta")
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"
        with pytest.raises(ValueError, match="only applies to the diamond"):
            seq_dist.main(["-i", fasta, "-r", fasta, "--sparse", out, "--mode", "score",
                           "--algorithm", "hmmer", "--params", '"-c":4'])


def test_diamond_failure_message_diagnoses_disk_exhaustion(tmp_path):
    """The GT2 crash: diamond fills its scratch during the seed search."""
    stderr = ("Searching alignments... No space left on device\n"
              + "No space left on device\n" * 27
              + "Error: Write error in HitBuffer\n")
    message = seq_dist.build_diamond_failure_message(
        "blastp", 1, stderr, "GT2.fasta", tmpdir=str(tmp_path),
        command=["diamond", "blastp", "-q", "GT2.fasta"])
    assert "ran out of disk space" in message
    assert str(tmp_path) in message
    # the mitigations a user can actually act on
    assert "--algorithm diamond_s" in message
    assert '"-b"' in message
    assert "TMPDIR" in message
    # it must not mislead the user into thinking streaming can help
    assert "streaming cannot reduce it" in message
    # 28 identical lines collapse to one
    assert "No space left on device (x27)" in message
    assert message.count("No space left on device") == 2


def test_diamond_failure_message_diagnoses_memory_and_options():
    killed = seq_dist.build_diamond_failure_message("blastp", -9, "some output\n", "in.fasta")
    assert "SIGKILL" in killed and "out-of-memory" in killed
    assert "--cpu" in killed

    alloc = seq_dist.build_diamond_failure_message(
        "blastp", 1, "Error: std::bad_alloc\n", "in.fasta")
    assert "could not allocate memory" in alloc

    bad_option = seq_dist.build_diamond_failure_message(
        "blastp", 1, "Error: Invalid option: no-such-flag\n", "in.fasta",
        diamond_params={"--no-such-flag": 1})
    assert "rejected one of its command line options" in bad_option
    assert "--no-such-flag" in bad_option

    bad_input = seq_dist.build_diamond_failure_message(
        "makedb", 1, "Error: Error detecting input file format. Input file seems to be empty.\n",
        "empty.fasta")
    assert "could not read an input or database file" in bad_input


def test_diamond_failure_message_falls_back_to_stderr():
    """An unrecognized failure must still surface diamond's own output."""
    message = seq_dist.build_diamond_failure_message(
        "blastp", 3, "Error: something entirely new\n", "in.fasta")
    assert "something entirely new" in message
    assert "exit code 3" in message
    assert "Suggested mitigations" not in message  # no guess when the cause is unknown


def test_seq_dist_reports_rejected_diamond_option(shared_datadir, capsys):
    """End-to-end: a bad --params flag reaches diamond and is reported usefully."""
    fasta = str(shared_datadir / "FeSOD_20.fasta")
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"
        with pytest.raises(SystemExit):
            seq_dist.main(["-i", fasta, "-r", fasta, "--sparse", out, "--mode", "score",
                           "--params", '"--no-such-flag":1'])
    captured = capsys.readouterr()
    assert "rejected one of its command line options" in captured.err
    assert "--no-such-flag" in captured.err
