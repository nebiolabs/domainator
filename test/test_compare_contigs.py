from domainator import compare_contigs, utils
import tempfile
import pandas as pd
from domainator.data_matrix import DataMatrix
import pytest

#import filecmp

def test_compare_contigs_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out_gb = output_dir + "/sorted.gb"
        out_dense = output_dir + "/dense.hdf5"
        out_sparse = output_dir + "/sparse.hdf5"
        compare_contigs.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "--ji", "0.5", "--ai", "0.1", "--dense", out_dense, "--sparse", out_sparse, "-o", out_gb])
        
        sparse_matrix = DataMatrix.from_file(out_sparse)
        assert sparse_matrix.shape[0] == 4
        assert sparse_matrix.shape[1] == 4
        assert all([x == 0.6 for x in sparse_matrix.itervalues()])
        dense_matrix = DataMatrix.from_file(out_dense)
        assert len(dense_matrix) == 4
        assert dense_matrix.shape[1] == 4
        assert all([x == 0.6 for x in dense_matrix.itervalues()])

def test_compare_contigs_2(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out_gb = output_dir + "/sorted.gb"
        out_dense = output_dir + "/dense.hdf5"
        out_sparse = output_dir + "/sparse.hdf5"
        compare_contigs.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "--ji", "0.5", "--ai", "0.1", "--dense", out_dense, "--sparse", out_sparse, "-o", out_gb, "--name_by_order"])
        sparse_matrix = DataMatrix.from_file(out_sparse)
        assert len(sparse_matrix) == 4
        assert sparse_matrix.shape[1] == 4
        assert all([x == 0.6 for x in sparse_matrix.itervalues()])
        dense_matrix = DataMatrix.from_file(out_dense)
        assert len(dense_matrix) == 4
        assert dense_matrix.shape[1] == 4
        assert all([x == 0.6 for x in dense_matrix.itervalues()])
        
        genbanks = list(utils.parse_seqfiles([out_gb]))
        for i in range(len(genbanks)): #TODO: this is a bad test because it doesn't test changes in order, because in this case the output is the same order as the input.
            prefix = str(i).zfill(1)
            assert genbanks[i].id.startswith(prefix+"_")

def test_compare_no_annotations_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out_gb = output_dir + "/sorted.gb"
        out_dense = output_dir + "/dense.hdf5"
        out_sparse = output_dir + "/sparse.hdf5"
        compare_contigs.main(["-i", str(shared_datadir / "pDONR201_multi_genemark.gb"), "--ji", "0.5", "--ai", "0.1", "--dense", out_dense, "--sparse", out_sparse, "-o", out_gb, "--name_by_order"])
        sparse_matrix = DataMatrix.from_file(out_sparse)
        assert len(sparse_matrix) == 4
        assert sparse_matrix.shape[1] == 4
        print(list(sparse_matrix.itervalues()))
        assert pytest.approx(sum(sparse_matrix.itervalues()),0.001) == (0.6 * 16) # every value should be 0.6
        # assert all([x == 0.6 for x in sparse_matrix.itervalues()])
        dense_matrix = DataMatrix.from_file(out_dense)
        assert len(dense_matrix) == 4
        assert dense_matrix.shape[1] == 4
        assert  pytest.approx(sum(dense_matrix.itervalues()),0.001) == (0.6 * 16) # every value should be 0.6
        
        # assert all([x == 0.6 for x in dense_matrix.itervalues()])
        
        genbanks = list(utils.parse_seqfiles([out_gb]))
        for i in range(len(genbanks)): #TODO: this is a bad test because it doesn't test changes in order, because in this case the output is the same order as the input.
            prefix = str(i).zfill(1)
            assert genbanks[i].id.startswith(prefix+"_")

def test_compare_contigs_database_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out_gb = output_dir + "/sorted.gb"
        out_dense = output_dir + "/dense.hdf5"
        out_dense_text = output_dir + "/dense.txt"
        out_sparse = output_dir + "/sparse.hdf5"

        compare_contigs.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator_multi_hmm_2.gb"), "--dense_text", out_dense_text, "--ji", "0.5", "--ai", "0.5", "--dense", out_dense, "--sparse", out_sparse, "-o", out_gb, "--name_by_order", "--databases", "pdonr_hmms_1"])

        sparse_matrix = DataMatrix.from_file(out_sparse)
        assert len(sparse_matrix) == 4
        assert sparse_matrix.shape[1] == 4
        assert  pytest.approx(sum(sparse_matrix.itervalues()),0.001) == 8.0
        assert sparse_matrix.data[0,0] == sparse_matrix.data[0,1]
        assert sparse_matrix.data[1,0] == sparse_matrix.data[1,1]
        assert sparse_matrix.data[2,2] == sparse_matrix.data[2,3]
        assert sparse_matrix.data[3,2] == sparse_matrix.data[3,3]

        dense_matrix = DataMatrix.from_file(out_dense)
        assert len(dense_matrix) == 4
        assert dense_matrix.shape[1] == 4
        assert  pytest.approx(sum(dense_matrix.itervalues()),0.001) == 8.0
        assert dense_matrix.data[0,0] == dense_matrix.data[0,1]
        assert dense_matrix.data[1,0] == dense_matrix.data[1,1]
        assert dense_matrix.data[2,2] == dense_matrix.data[2,3]
        assert dense_matrix.data[3,2] == dense_matrix.data[3,3]
        
        gb_out_text = open(out_gb).read()
        assert "pdonr_hmms_1" in gb_out_text
        assert "pdonr_hmms_2" in gb_out_text
        
        genbanks = list(utils.parse_seqfiles([out_gb]))
        for i in range(len(genbanks)): #TODO: this is a bad test because it doesn't test changes in order, because in this case the output is the same order as the input.
            prefix = str(i).zfill(1)
            assert genbanks[i].id.startswith(prefix+"_")
