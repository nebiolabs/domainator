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


def test_seq_dist_mst_knn_rejects_k_below_one(shared_datadir):
    fasta = str(shared_datadir / "FeSOD_20.fasta")
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"
        with pytest.raises(SystemExit):
            seq_dist.main(["-i", fasta, "-r", fasta, "--sparse", out, "--mode", "score", "--mst_knn", "0"])


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
