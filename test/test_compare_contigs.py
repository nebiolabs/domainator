from domainator import compare_contigs, utils
import tempfile
import pandas as pd
import numpy as np
import scipy.sparse
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
        assert sparse_matrix.data.nnz == 4
        assert pytest.approx(sum(sparse_matrix.itervalues()),0.001) == (0.6 * 4)
        dense_matrix = DataMatrix.from_file(out_dense)
        assert len(dense_matrix) == 4
        assert dense_matrix.shape[1] == 4
        assert pytest.approx(sum(dense_matrix.itervalues()),0.001) == (0.6 * 4)
        assert (dense_matrix.toarray() == np.eye(4) * 0.6).all()
        
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


def test_compare_contigs_progress_flag(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out_dense = output_dir + "/dense.hdf5"

        compare_contigs.main([
            "-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"),
            "--ji", "0.5",
            "--ai", "0.1",
            "--dense", out_dense,
            "--progress",
        ])

        dense_matrix = DataMatrix.from_file(out_dense)
        assert dense_matrix.shape == (4, 4)
        assert all([x == 0.6 for x in dense_matrix.itervalues()])


def test_compare_contigs_cpu_flag_matches_serial(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        serial_out = output_dir + "/serial.hdf5"
        parallel_out = output_dir + "/parallel.hdf5"

        args = [
            "-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"),
            "--ji", "0.5",
            "--ai", "0.1",
        ]

        compare_contigs.main(args + ["--dense", serial_out, "--cpu", "1"])
        compare_contigs.main(args + ["--dense", parallel_out, "--cpu", "2"])

        serial_matrix = DataMatrix.from_file(serial_out)
        parallel_matrix = DataMatrix.from_file(parallel_out)

        assert np.array_equal(serial_matrix.toarray(), parallel_matrix.toarray())


def test_compare_contigs_cpu_flag_matches_serial_for_single_metric(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        serial_out = output_dir + "/serial_single_metric.hdf5"
        parallel_out = output_dir + "/parallel_single_metric.hdf5"

        args = [
            "-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"),
            "--ji", "0.5",
            "--ai", "0.0",
        ]

        compare_contigs.main(args + ["--dense", serial_out, "--cpu", "1"])
        compare_contigs.main(args + ["--dense", parallel_out, "--cpu", "2"])

        serial_matrix = DataMatrix.from_file(serial_out)
        parallel_matrix = DataMatrix.from_file(parallel_out)

        assert np.array_equal(serial_matrix.toarray(), parallel_matrix.toarray())


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


def test_compare_contigs_singleton_adjacency_scores_zero():
    records = [
        _make_compare_contig_record("rec1", ["A"]),
        _make_compare_contig_record("rec2", ["A"]),
    ]

    with tempfile.TemporaryDirectory() as output_dir:
        input_path = output_dir + "/input.gb"
        out_dense = output_dir + "/dense.hdf5"

        with open(input_path, "w") as handle:
            write_genbank(records, handle)

        compare_contigs.main([
            "-i", input_path,
            "--ji", "0.0",
            "--ai", "1.0",
            "--dense", out_dense,
        ])

        dense_matrix = DataMatrix.from_file(out_dense)

        assert dense_matrix.data[0, 0] == 1.0
        assert dense_matrix.data[1, 1] == 1.0
        assert dense_matrix.data[0, 1] == 0.0
        assert dense_matrix.data[1, 0] == 0.0


def test_compare_contigs_lb_preserves_pairs_that_only_clear_threshold_after_combining_metrics():
    records = [
        _make_compare_contig_record("rec1", ["A", "B", "C"]),
        _make_compare_contig_record("rec2", ["A", "B", "D"]),
    ]

    with tempfile.TemporaryDirectory() as output_dir:
        input_path = output_dir + "/input.gb"
        out_dense = output_dir + "/dense.hdf5"

        with open(input_path, "w") as handle:
            write_genbank(records, handle)

        compare_contigs.main([
            "-i", input_path,
            "--ji", "0.5",
            "--ai", "0.5",
            "--lb", "0.3",
            "--dense", out_dense,
        ])

        dense_matrix = DataMatrix.from_file(out_dense)

        assert dense_matrix.data[0, 1] == pytest.approx((0.5 * 0.5) + (0.5 * (1.0 / 3.0)))
        assert dense_matrix.data[1, 0] == pytest.approx((0.5 * 0.5) + (0.5 * (1.0 / 3.0)))


class _FixedMatrixMetric(compare_contigs.ContigMetric):
    def __init__(self, weight, matrix):
        self.weight = weight
        self.matrix = scipy.sparse.csr_array(matrix, dtype=np.float64)

    def compute(self, recs, cpu=1, progress=False):
        assert len(recs) == self.matrix.shape[0]
        return self.matrix


class _CaptureScoresReport(compare_contigs.DistanceReport):
    def __init__(self):
        self.scores = None

    def write(self, recs, scores, options):
        self.scores = scipy.sparse.csr_array(scores)


def test_compare_contigs_mst_knn_matches_batch_transform():
    from domainator import transform_matrix
    from scipy.sparse.csgraph import connected_components

    # two clusters: {rec1..rec4} share domains, {rec5,rec6} share a disjoint set.
    records = [
        _make_compare_contig_record("rec1", ["A", "B", "C"]),
        _make_compare_contig_record("rec2", ["A", "B", "D"]),
        _make_compare_contig_record("rec3", ["A", "C", "E"]),
        _make_compare_contig_record("rec4", ["B", "C", "F"]),
        _make_compare_contig_record("rec5", ["X", "Y", "Z"]),
        _make_compare_contig_record("rec6", ["X", "Y", "W"]),
    ]

    with tempfile.TemporaryDirectory() as output_dir:
        input_path = output_dir + "/input.gb"
        with open(input_path, "w") as handle:
            write_genbank(records, handle)

        base_args = ["-i", input_path, "--ji", "0.5", "--ai", "0.5"]

        # batch reference: full matrix, then post-hoc mst_knn
        full_sparse = output_dir + "/full.hdf5"
        batch_sparse = output_dir + "/batch.hdf5"
        compare_contigs.main(base_args + ["--sparse", full_sparse])
        transform_matrix.main(["-i", full_sparse, "--sparse", batch_sparse, "--mst_knn", "2"])

        # streaming path
        stream_sparse = output_dir + "/stream.hdf5"
        compare_contigs.main(base_args + ["--sparse", stream_sparse, "--mst_knn", "2"])

        batch_dm = DataMatrix.from_file(batch_sparse)
        stream_dm = DataMatrix.from_file(stream_sparse)
        batch = batch_dm.data.toarray()
        stream = stream_dm.data.toarray()
        assert batch_dm.rows == stream_dm.rows

        # Equal-weight similarity edges make the maximum spanning forest non-unique (and a
        # tie-chosen MST edge may or may not coincide with a kNN edge, changing the union
        # size). The tie-invariant property is the connected-components partition, which is
        # exactly what mst_knn is used for in SSN clustering. Byte-exact equivalence on
        # distinct-weight inputs is covered in test_transform_matrix.
        _, batch_labels = connected_components(batch > 0, directed=False)
        _, stream_labels = connected_components(stream > 0, directed=False)
        # same partition: two nodes share a component in batch iff they do in stream
        batch_same = batch_labels[:, None] == batch_labels[None, :]
        stream_same = stream_labels[:, None] == stream_labels[None, :]
        np.testing.assert_array_equal(batch_same, stream_same)
        # sanity: the two domain clusters stay separate
        assert connected_components(stream > 0, directed=False)[0] == 2


def test_compare_contigs_k_applies_to_combined_scores():
    records = [
        _make_compare_contig_record("rec1", []),
        _make_compare_contig_record("rec2", []),
        _make_compare_contig_record("rec3", []),
        _make_compare_contig_record("rec4", []),
    ]

    metric_1 = _FixedMatrixMetric(
        0.5,
        np.array([
            [1.0, 0.9, 0.4, 0.6],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]),
    )
    metric_2 = _FixedMatrixMetric(
        0.5,
        np.array([
            [1.0, 0.1, 0.7, 0.6],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]),
    )

    capture = _CaptureScoresReport()

    with tempfile.TemporaryDirectory() as output_dir:
        input_path = output_dir + "/input.gb"

        with open(input_path, "w") as handle:
            write_genbank(records, handle)

        compare_contigs.compare_contigs(
            [input_path],
            [metric_1, metric_2],
            [capture],
            k=2,
            contigs=[record.id for record in records],
            name_by_order=False,
        )

    out = capture.scores.toarray()

    assert np.count_nonzero(out[0]) == 2
    assert out[0, 0] == pytest.approx(1.0)
    assert out[0, 1] == 0.0
    assert out[0, 2] == 0.0
    assert out[0, 3] == pytest.approx(0.6)


def test_streaming_csr_chunk_builder_rejects_out_of_order_chunks():
    builder = compare_contigs._StreamingCSRChunkBuilder((3, 3))

    with pytest.raises(RuntimeError, match="out of order"):
        builder.append_chunk(1, 2, scipy.sparse.csr_array([[1.0, 0.0, 0.0]]))


def test_iter_query_chunks_caps_max_rows_per_chunk():
    chunks = compare_contigs._iter_query_chunks(100, cpu=1, chunks_per_worker=1, max_rows_per_chunk=15)

    assert max(end - start for start, end in chunks) == 15
    assert chunks[0] == (0, 15)
    assert chunks[-1] == (90, 100)


def test_compare_contigs_combined_sparse_metric_path_matches_fallback():
    records = [
        _make_compare_contig_record("rec1", ["A", "B", "C"]),
        _make_compare_contig_record("rec2", ["A", "B", "D"]),
        _make_compare_contig_record("rec3", ["A", "E"]),
    ]

    with tempfile.TemporaryDirectory() as output_dir:
        input_path = output_dir + "/input.gb"
        with open(input_path, "w") as handle:
            write_genbank(records, handle)

        recs = list(utils.parse_seqfiles([input_path]))
        prepared_metrics = [
            compare_contigs.JaccardIndex(0.5).prepare_sparse_metric(recs),
            compare_contigs.AdjacencyIndex(0.5).prepare_sparse_metric(recs),
        ]

        combined = compare_contigs._compute_combined_sparse_similarity_matrix(
            prepared_metrics,
            k=2,
            lb=0.2,
            cpu=1,
        )

        fallback = scipy.sparse.csr_array((len(recs), len(recs)), dtype=np.float64)
        for metric in [compare_contigs.JaccardIndex(0.5), compare_contigs.AdjacencyIndex(0.5)]:
            fallback = fallback + (metric.weight * metric.compute(recs, cpu=1, progress=False))
        fallback = compare_contigs._keep_top_k_per_row(fallback, 2)
        fallback = compare_contigs._prune_scores_inplace(fallback, 0.2)

        assert np.allclose(combined.toarray(), fallback.toarray())
