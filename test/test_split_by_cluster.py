from pathlib import Path

import numpy as np
import pyhmmer
import pytest

from domainator import build_ssn, split_by_cluster
from domainator.Bio import SeqIO
from domainator.data_matrix import DenseDataMatrix
from domainator.utils import parse_seqfiles, pyhmmer_decode


def _write_test_matrix(matrix_path, labels, output_type="dense_text"):
    data = np.zeros((len(labels), len(labels)), dtype=float)
    if len(labels) >= 2:
        data[0, 1] = data[1, 0] = 10.0
    if len(labels) >= 4:
        data[2, 3] = data[3, 2] = 9.0

    matrix = DenseDataMatrix(data, labels, labels, data_type="score")
    matrix.write(matrix_path, output_type)
    return matrix


def _expected_cluster_members(labels, matrix, lb, min_cluster_size):
    cluster_labels = build_ssn.cluster_labels_from_graph(matrix, lb)
    cluster_members = {}
    for label, cluster_id in zip(labels, cluster_labels):
        cluster_members.setdefault(int(cluster_id), []).append(label)

    regular_clusters = {
        cluster_id: members
        for cluster_id, members in cluster_members.items()
        if len(members) >= min_cluster_size
    }
    small_clusters = [
        label
        for cluster_id, members in cluster_members.items()
        if len(members) < min_cluster_size
        for label in members
    ]
    return regular_clusters, small_clusters


def _read_seq_ids(path, file_type):
    return [record.id for record in parse_seqfiles([path], filetype_override=file_type)]


def _read_hmm_names(path):
    with pyhmmer.plan7.HMMFile(path) as hmm_file:
        return [pyhmmer_decode(model.name) for model in hmm_file]


def _read_cm_names(path):
    return [name for name, _ in split_by_cluster._iter_cm_blocks(path)]


@pytest.mark.parametrize(
    "matrix_output_type,matrix_name",
    [("dense_text", "matrix.tsv"), ("dense", "matrix.dense.hdf5"), ("sparse", "matrix.sparse.hdf5")],
)
def test_split_by_cluster_genbank_matrix_formats(shared_datadir, tmp_path, matrix_output_type, matrix_name):
    input_path = shared_datadir / "FeSOD_20.gb"
    labels = [record.id for record in parse_seqfiles([input_path], filetype_override="genbank")]
    matrix_path = tmp_path / matrix_name
    matrix = _write_test_matrix(matrix_path, labels, output_type=matrix_output_type)
    expected_clusters, expected_small = _expected_cluster_members(labels, matrix, lb=0, min_cluster_size=2)

    split_by_cluster.main([
        "-i", str(input_path),
        "--matrix", str(matrix_path),
        "--outdir", str(tmp_path / "out"),
        "--lb", "0",
        "--min_cluster_size", "2",
        "--write_small_clusters",
    ])

    outdir = tmp_path / "out"
    for cluster_id, members in expected_clusters.items():
        assert _read_seq_ids(outdir / f"{cluster_id}.gb", "genbank") == members
    assert _read_seq_ids(outdir / "small_clusters.gb", "genbank") == expected_small


def test_split_by_cluster_fasta(shared_datadir, tmp_path):
    input_path = shared_datadir / "FeSOD_20.fasta"
    labels = [record.id for record in parse_seqfiles([input_path], filetype_override="fasta")]
    matrix_path = tmp_path / "matrix.tsv"
    matrix = _write_test_matrix(matrix_path, labels, output_type="dense_text")
    expected_clusters, expected_small = _expected_cluster_members(labels, matrix, lb=0, min_cluster_size=2)

    split_by_cluster.main([
        "-i", str(input_path),
        "--matrix", str(matrix_path),
        "--outdir", str(tmp_path / "out"),
        "--lb", "0",
        "--min_cluster_size", "2",
        "--write_small_clusters",
    ])

    outdir = tmp_path / "out"
    for cluster_id, members in expected_clusters.items():
        assert _read_seq_ids(outdir / f"{cluster_id}.fasta", "fasta") == members
    assert _read_seq_ids(outdir / "small_clusters.fasta", "fasta") == expected_small


def test_split_by_cluster_hmm(shared_datadir, tmp_path):
    input_path = shared_datadir / "pdonr_hmms.hmm"
    labels = _read_hmm_names(input_path)
    matrix_path = tmp_path / "matrix.tsv"
    matrix = _write_test_matrix(matrix_path, labels, output_type="dense_text")
    expected_clusters, expected_small = _expected_cluster_members(labels, matrix, lb=0, min_cluster_size=2)

    split_by_cluster.main([
        "-i", str(input_path),
        "--matrix", str(matrix_path),
        "--outdir", str(tmp_path / "out"),
        "--lb", "0",
        "--min_cluster_size", "2",
        "--write_small_clusters",
    ])

    outdir = tmp_path / "out"
    for cluster_id, members in expected_clusters.items():
        assert _read_hmm_names(outdir / f"{cluster_id}.hmm") == members
    assert _read_hmm_names(outdir / "small_clusters.hmm") == expected_small


def test_split_by_cluster_cm_small_clusters(shared_datadir, tmp_path):
    source_cm = shared_datadir / "RF00042.cm"
    input_path = tmp_path / "single.cm"
    first_name, first_block = next(split_by_cluster._iter_cm_blocks(source_cm))
    input_path.write_text(first_block, encoding="utf-8")
    labels = [first_name]
    matrix_path = tmp_path / "matrix.tsv"
    _write_test_matrix(matrix_path, labels, output_type="dense_text")

    split_by_cluster.main([
        "-i", str(input_path),
        "--matrix", str(matrix_path),
        "--outdir", str(tmp_path / "out"),
        "--lb", "0",
        "--min_cluster_size", "2",
        "--write_small_clusters",
    ])

    outdir = tmp_path / "out"
    assert not (outdir / "1.cm").exists()
    assert _read_cm_names(outdir / "small_clusters.cm") == labels


def test_split_by_cluster_requires_exact_label_match(shared_datadir, tmp_path):
    input_path = shared_datadir / "FeSOD_20.gb"
    labels = [record.id for record in parse_seqfiles([input_path], filetype_override="genbank")][:3]
    matrix_path = tmp_path / "matrix.tsv"
    _write_test_matrix(matrix_path, labels, output_type="dense_text")

    with pytest.raises(ValueError, match="do not match exactly"):
        split_by_cluster.main([
            "-i", str(input_path),
            "--matrix", str(matrix_path),
            "--outdir", str(tmp_path / "out"),
        ])


def test_split_by_cluster_skips_small_clusters_by_default(shared_datadir, tmp_path):
    input_path = shared_datadir / "FeSOD_20.gb"
    labels = [record.id for record in parse_seqfiles([input_path], filetype_override="genbank")]
    matrix_path = tmp_path / "matrix.tsv"
    matrix = _write_test_matrix(matrix_path, labels, output_type="dense_text")
    expected_clusters, _expected_small = _expected_cluster_members(labels, matrix, lb=0, min_cluster_size=2)

    split_by_cluster.main([
        "-i", str(input_path),
        "--matrix", str(matrix_path),
        "--outdir", str(tmp_path / "out"),
        "--lb", "0",
        "--min_cluster_size", "2",
    ])

    outdir = tmp_path / "out"
    for cluster_id, members in expected_clusters.items():
        assert _read_seq_ids(outdir / f"{cluster_id}.gb", "genbank") == members
    assert not (outdir / "small_clusters.gb").exists()


def test_split_by_cluster_pad_on_search_hits_centers_genbank(shared_datadir, tmp_path):
    input_path = shared_datadir / "pDONR_201_domain_search.gb"
    labels = [record.id for record in parse_seqfiles([input_path], filetype_override="genbank")]
    matrix_path = tmp_path / "matrix.tsv"
    _write_test_matrix(matrix_path, labels, output_type="dense_text")

    split_by_cluster.main([
        "-i", str(input_path),
        "--matrix", str(matrix_path),
        "--outdir", str(tmp_path / "out"),
        "--pad_on_search_hits",
    ])

    records = list(SeqIO.parse(tmp_path / "out" / "1.gb", "genbank"))
    assert len(records) == 1
    assert records[0].id == "pDONR201_1"
    assert len(records[0]) == 4470
    best_hit_features = [feature for feature in records[0].features if feature.type == "Domain_Search"]
    assert len(best_hit_features) == 1
    assert best_hit_features[0].qualifiers["name"] == ["CcdA"]


def test_split_by_cluster_pad_on_search_hits_requires_domain_search(shared_datadir, tmp_path):
    input_path = shared_datadir / "FeSOD_20.gb"
    labels = [record.id for record in parse_seqfiles([input_path], filetype_override="genbank")]
    matrix_path = tmp_path / "matrix.tsv"
    _write_test_matrix(matrix_path, labels, output_type="dense_text")

    with pytest.raises(ValueError, match="does not contain any Domain_Search"):
        split_by_cluster.main([
            "-i", str(input_path),
            "--matrix", str(matrix_path),
            "--outdir", str(tmp_path / "out"),
            "--pad_on_search_hits",
        ])


def test_split_by_cluster_pad_on_search_hits_requires_genbank(shared_datadir, tmp_path):
    input_path = shared_datadir / "FeSOD_20.fasta"
    labels = [record.id for record in parse_seqfiles([input_path], filetype_override="fasta")]
    matrix_path = tmp_path / "matrix.tsv"
    _write_test_matrix(matrix_path, labels, output_type="dense_text")

    with pytest.raises(ValueError, match="only supported for GenBank"):
        split_by_cluster.main([
            "-i", str(input_path),
            "--matrix", str(matrix_path),
            "--outdir", str(tmp_path / "out"),
            "--pad_on_search_hits",
        ])


def test_split_by_cluster_refuses_to_overwrite(shared_datadir, tmp_path):
    input_path = shared_datadir / "FeSOD_20.gb"
    labels = [record.id for record in parse_seqfiles([input_path], filetype_override="genbank")]
    matrix_path = tmp_path / "matrix.tsv"
    _write_test_matrix(matrix_path, labels, output_type="dense_text")
    outdir = tmp_path / "out"
    outdir.mkdir()
    (outdir / "1.gb").write_text("already here", encoding="utf-8")

    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        split_by_cluster.main([
            "-i", str(input_path),
            "--matrix", str(matrix_path),
            "--outdir", str(outdir),
            "--lb", "0",
            "--min_cluster_size", "2",
        ])


def test_split_by_cluster_reuses_build_ssn_mst_cut_logic(tmp_path):
    labels = ["A", "B", "C", "D", "E"]
    matrix = DenseDataMatrix(
        np.array(
            [
                [0, 9, 0, 0, 0],
                [9, 0, 5, 0, 0],
                [0, 5, 0, 0, 0],
                [0, 0, 0, 0, 7],
                [0, 0, 0, 7, 0],
            ],
            dtype=float,
        ),
        labels,
        labels,
    )
    matrix_path = tmp_path / "matrix.tsv"
    matrix.write(matrix_path, "dense_text")

    loaded_matrix, label_to_cluster, cluster_sizes = split_by_cluster._load_cluster_membership(matrix_path, lb=5)
    expected = build_ssn.cluster_labels_from_graph(loaded_matrix, 5)

    assert [label_to_cluster[label] for label in labels] == [int(x) for x in expected]
    assert sorted(cluster_sizes.values()) == [1, 2, 2]