from domainator import compare_contigs, utils
import tempfile
import pandas as pd
from domainator.data_matrix import DataMatrix
import pytest
from domainator import DOMAIN_FEATURE_NAME
from domainator.Bio.Seq import Seq
from domainator.Bio.SeqFeature import SeqFeature, FeatureLocation
from domainator.Bio.SeqRecord import SeqRecord
from domainator.utils import write_genbank

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
        
        _assert_name_by_order_prefixes(list(utils.parse_seqfiles([out_gb])))

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
        
        _assert_name_by_order_prefixes(list(utils.parse_seqfiles([out_gb])))

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
        
        _assert_name_by_order_prefixes(list(utils.parse_seqfiles([out_gb])))


def test_compare_contigs_max_output_gb_blocks_dense_output(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out_dense = output_dir + "/dense.hdf5"

        with pytest.raises(SystemExit, match="--max_output_gb"):
            compare_contigs.main([
                "-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"),
                "--ji", "0.5",
                "--ai", "0.1",
                "--dense", out_dense,
                "--max_output_gb", "0.000001",
            ])


def _make_compare_contig_record(record_id, domain_names):
    record = SeqRecord(Seq("ATGCATGCATGC"), id=record_id, name=record_id, description=record_id)
    record.annotations["molecule_type"] = "DNA"
    for index, domain_name in enumerate(domain_names):
        record.features.append(
            SeqFeature(
                FeatureLocation(index * 2, index * 2 + 2, strand=1),
                type=DOMAIN_FEATURE_NAME,
                qualifiers={
                    "database": ["test_db"],
                    "name": [domain_name],
                    "cds_id": [str(index)],
                    "description": [f"{domain_name} description"],
                    "evalue": ["1.0e-5"],
                    "score": ["10.0"],
                },
            )
        )
    return record


def _assert_name_by_order_prefixes(records):
    digits = len(str(len(records) - 1)) if len(records) > 1 else 1
    for index, record in enumerate(records):
        prefix, stripped_id = record.id.split("_", 1)
        assert prefix == str(index).zfill(digits)
        assert stripped_id


def test_compare_contigs_name_by_order_reorders_clustered_records():
    records = [
        _make_compare_contig_record("different", ["C"]),
        _make_compare_contig_record("similar_1", ["A", "B"]),
        _make_compare_contig_record("similar_2", ["A", "B"]),
    ]

    with tempfile.TemporaryDirectory() as output_dir:
        input_path = output_dir + "/input.gb"
        out_gb = output_dir + "/sorted.gb"

        with open(input_path, "w") as handle:
            write_genbank(records, handle)

        compare_contigs.main([
            "-i", input_path,
            "--ji", "1.0",
            "--ai", "0.0",
            "-o", out_gb,
            "--name_by_order",
        ])

        ordered_records = list(utils.parse_seqfiles([out_gb]))
        stripped_ids = [record.id.split("_", 1)[1] for record in ordered_records]

        _assert_name_by_order_prefixes(ordered_records)
        assert stripped_ids != [record.id for record in records]
        assert set(stripped_ids[:2]) == {"similar_1", "similar_2"}
        assert stripped_ids[2] == "different"


def test_compare_contigs_lb_thresholds_small_scores():
    records = [
        _make_compare_contig_record("rec1", ["A", "B"]),
        _make_compare_contig_record("rec2", ["A", "C"]),
    ]

    with tempfile.TemporaryDirectory() as output_dir:
        input_path = output_dir + "/input.gb"
        out_dense = output_dir + "/dense.hdf5"
        out_sparse = output_dir + "/sparse.hdf5"

        with open(input_path, "w") as handle:
            write_genbank(records, handle)

        compare_contigs.main([
            "-i", input_path,
            "--ji", "1.0",
            "--ai", "0.0",
            "--lb", "0.5",
            "--dense", out_dense,
            "--sparse", out_sparse,
        ])

        sparse_matrix = DataMatrix.from_file(out_sparse)
        dense_matrix = DataMatrix.from_file(out_dense)

        assert sparse_matrix.data.nnz == 2
        assert sparse_matrix.data[0, 0] == 1.0
        assert sparse_matrix.data[1, 1] == 1.0
        assert sparse_matrix.data[0, 1] == 0.0
        assert sparse_matrix.data[1, 0] == 0.0

        assert dense_matrix.data[0, 0] == 1.0
        assert dense_matrix.data[1, 1] == 1.0
        assert dense_matrix.data[0, 1] == 0.0
        assert dense_matrix.data[1, 0] == 0.0
