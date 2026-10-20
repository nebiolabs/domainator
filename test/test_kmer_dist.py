import tempfile

import numpy as np
import pytest

from domainator import kmer_dist
from domainator.data_matrix import DataMatrix


@pytest.mark.parametrize(
    "metric,expected",
    [
        (
            "jaccard",
            np.array(
                [
                    [1.0, 0.5, 0.75],
                    [0.5, 1.0, 2.0 / 3.0],
                    [0.75, 2.0 / 3.0, 1.0],
                ]
            ),
        ),
        (
            "max_containment",
            np.array(
                [
                    [1.0, 1.0, 1.0],
                    [1.0, 1.0, 1.0],
                    [1.0, 1.0, 1.0],
                ]
            ),
        ),
        (
            "min_containment",
            np.array(
                [
                    [1.0, 0.5, 0.75],
                    [0.5, 1.0, 2.0 / 3.0],
                    [0.75, 2.0 / 3.0, 1.0],
                ]
            ),
        ),
    ],
)
def test_kmer_dist_protein_metrics(shared_datadir, metric, expected):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"
        kmer_dist.main([
            "-i", str(shared_datadir / "kmer_dist_proteins.faa"),
            "--dense", out,
            "--metric", metric,
            "--ksize", "2",
        ])

        matrix = DataMatrix.from_file(out)
        assert matrix.data_type == "score"
        assert np.allclose(matrix.toarray(), expected)


def test_kmer_dist_protein_reference_matrix(shared_datadir):
    expected = np.array(
        [
            [0.75, 0.0],
            [2.0 / 3.0, 0.0],
            [1.0, 0.0],
        ]
    )

    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"
        kmer_dist.main([
            "-i", str(shared_datadir / "kmer_dist_proteins.faa"),
            "-r", str(shared_datadir / "kmer_dist_reference.faa"),
            "--dense", out,
            "--metric", "jaccard",
            "--ksize", "2",
        ])

        matrix = DataMatrix.from_file(out)
        assert matrix.shape == (3, 2)
        assert np.allclose(matrix.toarray(), expected)


def test_kmer_dist_sparse_top_k(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"
        kmer_dist.main([
            "-i", str(shared_datadir / "kmer_dist_proteins.faa"),
            "--sparse", out,
            "--metric", "jaccard",
            "--ksize", "2",
            "--max_hits_per_query", "1",
        ])

        matrix = DataMatrix.from_file(out)
        dense = matrix.toarray()
        assert np.array_equal(dense, np.eye(3))


def test_kmer_dist_bool_mode_with_lb(shared_datadir):
    expected = np.array(
        [
            [1, 0, 1],
            [0, 1, 0],
            [1, 0, 1],
        ]
    )

    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"
        kmer_dist.main([
            "-i", str(shared_datadir / "kmer_dist_proteins.faa"),
            "--dense", out,
            "--metric", "jaccard",
            "--mode", "bool",
            "--ksize", "2",
            "--lb", "0.7",
        ])

        matrix = DataMatrix.from_file(out)
        assert matrix.data_type == "bool"
        assert np.array_equal(matrix.toarray(), expected)


def test_kmer_dist_nucleotide_kmers_are_strand_invariant(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"
        kmer_dist.main([
            "-i", str(shared_datadir / "kmer_dist_dna.fna"),
            "--dense", out,
            "--metric", "jaccard",
            "--fasta_type", "nucleotide",
            "--nucleotide_kmers",
            "--ksize", "3",
        ])

        matrix = DataMatrix.from_file(out)
        assert np.array_equal(matrix.toarray(), np.ones((2, 2)))


def test_kmer_dist_nucleotide_fasta_requires_gene_call_or_nucleotide_kmers(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"
        with pytest.raises(ValueError, match="Nucleotide FASTA input requires either --gene_call"):
            kmer_dist.main([
                "-i", str(shared_datadir / "pDONR201.fasta"),
                "--dense", out,
                "--fasta_type", "nucleotide",
            ])


def test_kmer_dist_gene_call_all_for_nucleotide_fasta(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"
        kmer_dist.main([
            "-i", str(shared_datadir / "pDONR201.fasta"),
            "--dense", out,
            "--fasta_type", "nucleotide",
            "--gene_call", "all",
        ])

        matrix = DataMatrix.from_file(out)
        assert matrix.shape == (1, 1)
        assert matrix.toarray()[0, 0] == 1.0


def test_kmer_dist_translated_cds_genbank_runs(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"
        kmer_dist.main([
            "-i", str(shared_datadir / "pDONR201_domainator_circular.gb"),
            "--dense", out,
        ])

        matrix = DataMatrix.from_file(out)
        assert matrix.shape == (1, 1)
        assert matrix.toarray()[0, 0] == 1.0


def test_kmer_dist_ignores_noncanonical_protein_kmers(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        input_fasta = output_dir + "/noncanonical.faa"
        out = output_dir + "/out.hdf5"

        with open(input_fasta, "w", encoding="utf-8") as handle:
            handle.write(">seq1\nACXDE\n>seq2\nACYDE\n")

        kmer_dist.main([
            "-i", input_fasta,
            "--dense", out,
            "--metric", "jaccard",
            "--ksize", "2",
        ])

        matrix = DataMatrix.from_file(out)
        expected = np.array(
            [
                [1.0, 0.5],
                [0.5, 1.0],
            ]
        )
        assert np.allclose(matrix.toarray(), expected)


def test_kmer_dist_max_output_gb_blocks_dense_output(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"

        with pytest.raises(SystemExit, match="--max_output_gb"):
            kmer_dist.main([
                "-i", str(shared_datadir / "kmer_dist_proteins.faa"),
                "--dense", out,
                "--metric", "jaccard",
                "--ksize", "2",
                "--max_output_gb", "0.000001",
            ])


def test_kmer_dist_progress_flag(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"
        kmer_dist.main([
            "-i", str(shared_datadir / "kmer_dist_proteins.faa"),
            "--dense", out,
            "--metric", "jaccard",
            "--ksize", "2",
            "--progress",
        ])

        matrix = DataMatrix.from_file(out)
        assert matrix.shape == (3, 3)

def _off_diagonal_pairs(array):
    return {(i, j) for i, j in zip(*np.nonzero(array)) if i < j}


@pytest.mark.parametrize("sparsify_args", [["--knn", "2"], ["--mst"], ["--mst", "--knn", "2"]])
@pytest.mark.parametrize("cap", [None, "3"])
def test_kmer_dist_sparsify_matches_batch_transform(shared_datadir, sparsify_args, cap):
    from domainator import transform_matrix
    from scipy.sparse.csgraph import connected_components

    fasta = str(shared_datadir / "FeSOD_20.fasta")
    cap_args = [] if cap is None else ["--max_hits_per_query", cap]
    with tempfile.TemporaryDirectory() as output_dir:
        full = output_dir + "/full.hdf5"
        batch = output_dir + "/batch.hdf5"
        stream = output_dir + "/stream.hdf5"
        base = ["-i", fasta, "--ksize", "3", "--cpu", "1"] + cap_args
        kmer_dist.main(base + ["--sparse", full])
        transform_matrix.main(["-i", full, "--sparse", batch] + sparsify_args)
        kmer_dist.main(base + ["--sparse", stream] + sparsify_args)

        batch_array = DataMatrix.from_file(batch).toarray()
        stream_array = DataMatrix.from_file(stream).toarray()
        np.testing.assert_array_equal(stream_array, stream_array.T)
        np.testing.assert_array_equal(batch_array.diagonal(), stream_array.diagonal())
        if "--mst" in sparsify_args:
            # Tied k-mer similarities make the spanning forest non-unique; the partition is not.
            _, batch_labels = connected_components(batch_array > 0, directed=False)
            _, stream_labels = connected_components(stream_array > 0, directed=False)
            np.testing.assert_array_equal(batch_labels[:, None] == batch_labels[None, :],
                                          stream_labels[:, None] == stream_labels[None, :])
        else:
            np.testing.assert_array_equal(stream_array, batch_array)


def test_kmer_dist_sparsify_bool_mode_keeps_selected_edges(shared_datadir):
    fasta = str(shared_datadir / "FeSOD_20.fasta")
    with tempfile.TemporaryDirectory() as output_dir:
        score_out = output_dir + "/score.hdf5"
        bool_out = output_dir + "/bool.hdf5"
        base = ["-i", fasta, "--ksize", "3", "--cpu", "1", "--knn", "2"]
        kmer_dist.main(base + ["--sparse", score_out])
        kmer_dist.main(base + ["--sparse", bool_out, "--mode", "bool"])
        score_array = DataMatrix.from_file(score_out).toarray()
        bool_matrix = DataMatrix.from_file(bool_out)
        assert bool_matrix.data_type == "bool"
        np.testing.assert_array_equal(bool_matrix.toarray() != 0, score_array != 0)


def test_kmer_dist_sparsify_rejects_invalid_use(shared_datadir):
    proteins = str(shared_datadir / "kmer_dist_proteins.faa")
    reference = str(shared_datadir / "kmer_dist_reference.faa")
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.hdf5"
        with pytest.raises(ValueError, match="square"):
            kmer_dist.main(["-i", proteins, "-r", reference, "--sparse", out, "--mst"])
        with pytest.raises(ValueError, match="only supported with --mode"):
            kmer_dist.main(["-i", proteins, "--dense", out, "--mode", "score_dist", "--knn", "1"])
        with pytest.raises(SystemExit):
            kmer_dist.main(["-i", proteins, "--sparse", out, "--knn", "0"])
