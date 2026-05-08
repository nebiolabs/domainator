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
            "-k", "1",
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