from domainator import build_ssn
from domainator.utils import OTHER_COLOR
from domainator.data_matrix import DenseDataMatrix, SparseDataMatrix, MaxTree, mst_knn_edge_index_dict
import pytest
import tempfile
import pandas as pd
import numpy as np
from pathlib import Path
from helpers import compare_files
import re
import scipy.sparse
from scipy.sparse.csgraph import connected_components

@pytest.mark.parametrize("input_file,expected_output",
[
    ["FeSOD_dist.tsv","ssn_FeSOD.xgmml"],
    ["FeSOD_dist.sparse.hdf5","ssn_FeSOD.sparse.xgmml"],
    ["FeSOD_dist.dense.hdf5","ssn_FeSOD.xgmml"]
])
def test_build_ssn(input_file, expected_output, shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        metadata = str(shared_datadir / "FeSOD_metadata.tsv")
        out_clusters = output_dir + f"/{input_file}_out_clusters.tsv"
        out_cytoscape = output_dir + f"/{input_file}_out.xgmml"
        build_ssn.main(["-i", str(shared_datadir / input_file),"--xgmml", out_cytoscape, "--lb", "175", "--color_by", "SSN_cluster", "--cluster_tsv", out_clusters, "--no_cluster_header", "--metadata", metadata])
        assert Path(out_cytoscape).is_file()
        assert Path(out_clusters).is_file()
        compare_files(out_clusters,shared_datadir/'ssn_FeSOD_clusters.tsv')
        compare_files(out_cytoscape, shared_datadir/expected_output, skip_lines=2)

@pytest.mark.parametrize("input_file,expected_output",
[
    ["FeSOD_dist.tsv","ssn_FeSOD.xgmml"],
])
def test_build_ssn_2(input_file, expected_output, shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        metadata = str(shared_datadir / "FeSOD_metadata.tsv")
        out_clusters = output_dir + f"/{input_file}_out_clusters.tsv"
        out_cytoscape = output_dir + f"/{input_file}_out.xgmml"
        build_ssn.main(["-i", str(shared_datadir / input_file),"--xgmml", out_cytoscape, "--lb", "175", "--color_by", "SSN_cluster", "--cluster_tsv", out_clusters, "--metadata", metadata])
        assert Path(out_cytoscape).is_file()
        assert Path(out_clusters).is_file()
        compare_files(out_clusters,shared_datadir/'ssn_FeSOD_clusters_header.tsv')
        compare_files(out_cytoscape, shared_datadir/expected_output, skip_lines=2)


def test_build_ssn_3(shared_datadir):
    input_file = "FeSOD_dist.tsv"
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        metadata = str(shared_datadir / "FeSOD_metadata.tsv")
        out_clusters = output_dir + f"/{input_file}_out_clusters.tsv"
        out_cytoscape = output_dir + f"/{input_file}_out.xgmml"
        build_ssn.main(["-i", str(shared_datadir / input_file),"--xgmml", out_cytoscape, "--lb", "175", "--color_by", "SSN_cluster",
                        "--cluster_tsv", out_clusters, "--metadata", metadata, "--color_table_out", output_dir + "/color_table.tsv"])
        assert Path(out_cytoscape).is_file()
        assert Path(out_clusters).is_file()
        assert Path(output_dir + "/color_table.tsv").is_file()
        compare_files(out_clusters,shared_datadir/'ssn_FeSOD_clusters_header.tsv')
        compare_files(out_cytoscape, shared_datadir/"ssn_FeSOD.xgmml", skip_lines=2)

        color_table_dict = {}
        with open(output_dir + "/color_table.tsv", "r") as f:
            for line in f:
                domain, color = line.strip().split("\t")
                color_table_dict[domain] = color
        assert len(color_table_dict) == 3
        assert set(color_table_dict.keys()) == {"1","2","3"}
        assert all([re.match(r"#[0-9a-fA-F]{6}",x) for x in color_table_dict.values()])
        assert len(set(color_table_dict.values())) == 3

def test_build_ssn_max_color_groups(shared_datadir):
    """Only the largest clusters get their own color; the tail shares OTHER_COLOR."""
    input_file = "FeSOD_dist.tsv"
    with tempfile.TemporaryDirectory() as output_dir:
        metadata = str(shared_datadir / "FeSOD_metadata.tsv")
        out_color_table = output_dir + "/color_table.tsv"
        build_ssn.main(["-i", str(shared_datadir / input_file), "--xgmml", output_dir + "/out.xgmml",
                        "--lb", "175", "--color_by", "SSN_cluster", "--cluster_tsv", output_dir + "/clusters.tsv",
                        "--metadata", metadata, "--color_table_out", out_color_table,
                        "--max_color_groups", "2"])

        color_table = {}
        with open(out_color_table, "r") as f:
            for line in f:
                cluster, color = line.strip().split("\t")
                color_table[cluster] = color

        # all three clusters are still listed, so the table still round-trips
        assert list(color_table.keys()) == ["1", "2", "3"]
        assert len({color_table["1"], color_table["2"]}) == 2
        assert OTHER_COLOR not in {color_table["1"], color_table["2"]}
        assert color_table["3"] == OTHER_COLOR


def test_build_ssn_4(shared_datadir):
    input_file = "FeSOD_dist.tsv"
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        metadata = str(shared_datadir / "FeSOD_metadata.tsv")
        out_clusters = output_dir + f"/{input_file}_out_clusters.tsv"
        out_cytoscape = output_dir + f"/{input_file}_out.xgmml"
        build_ssn.main(["-i", str(shared_datadir / input_file),"--xgmml", out_cytoscape, "--lb", "175", "--color_by", "SSN_cluster",
                        "--cluster_tsv", out_clusters, "--metadata", metadata, "--color_table_out", output_dir + "/color_table.tsv", "--color_table", str(shared_datadir / "color_table_123.tsv")])
        assert Path(out_cytoscape).is_file()
        assert Path(out_clusters).is_file()
        assert Path(output_dir + "/color_table.tsv").is_file()
        compare_files(out_clusters,shared_datadir/'ssn_FeSOD_clusters_header.tsv')
        compare_files(output_dir + "/color_table.tsv", shared_datadir/"color_table_123.tsv")


def test_build_ssn_mst(shared_datadir):
    """Test that MST mode produces correct clusters with fewer edges"""
    input_file = "FeSOD_dist.tsv"
    with tempfile.TemporaryDirectory() as output_dir:
        metadata = str(shared_datadir / "FeSOD_metadata.tsv")
        
        # Run without MST
        out_clusters_full = output_dir + f"/{input_file}_out_clusters_full.tsv"
        out_cytoscape_full = output_dir + f"/{input_file}_out_full.xgmml"
        build_ssn.main(["-i", str(shared_datadir / input_file), "--xgmml", out_cytoscape_full, "--lb", "175", 
                        "--color_by", "SSN_cluster", "--cluster_tsv", out_clusters_full, "--metadata", metadata])
        
        # Run with MST
        out_clusters_mst = output_dir + f"/{input_file}_out_clusters_mst.tsv"
        out_cytoscape_mst = output_dir + f"/{input_file}_out_mst.xgmml"
        build_ssn.main(["-i", str(shared_datadir / input_file), "--xgmml", out_cytoscape_mst, "--lb", "175", 
                        "--color_by", "SSN_cluster", "--cluster_tsv", out_clusters_mst, "--metadata", metadata, "--mst"])
        
        # Both files should be created
        assert Path(out_cytoscape_full).is_file()
        assert Path(out_cytoscape_mst).is_file()
        assert Path(out_clusters_full).is_file()
        assert Path(out_clusters_mst).is_file()
        
        # Clusters should be identical
        compare_files(out_clusters_mst, out_clusters_full)
        
        # MST should have fewer edges than the full graph
        # Count edges in each file by counting lines with "<edge"
        with open(out_cytoscape_full, 'r') as f:
            full_edges = sum(1 for line in f if '<edge' in line)
        
        with open(out_cytoscape_mst, 'r') as f:
            mst_edges = sum(1 for line in f if '<edge' in line)
        
        assert mst_edges < full_edges, f"MST should have fewer edges ({mst_edges}) than full graph ({full_edges})"
        assert mst_edges > 0, "MST should have at least one edge"


def test_build_ssn_mst_knn(shared_datadir):
    input_file = "FeSOD_dist.tsv"
    with tempfile.TemporaryDirectory() as output_dir:
        metadata = str(shared_datadir / "FeSOD_metadata.tsv")

        out_clusters_full = output_dir + f"/{input_file}_out_clusters_full.tsv"
        out_cytoscape_full = output_dir + f"/{input_file}_out_full.xgmml"
        build_ssn.main([
            "-i", str(shared_datadir / input_file),
            "--xgmml", out_cytoscape_full,
            "--lb", "175",
            "--color_by", "SSN_cluster",
            "--cluster_tsv", out_clusters_full,
            "--metadata", metadata,
        ])

        out_clusters_mst = output_dir + f"/{input_file}_out_clusters_mst.tsv"
        out_cytoscape_mst = output_dir + f"/{input_file}_out_mst.xgmml"
        build_ssn.main([
            "-i", str(shared_datadir / input_file),
            "--xgmml", out_cytoscape_mst,
            "--lb", "175",
            "--color_by", "SSN_cluster",
            "--cluster_tsv", out_clusters_mst,
            "--metadata", metadata,
            "--mst",
        ])

        out_clusters_mst_knn = output_dir + f"/{input_file}_out_clusters_mst_knn.tsv"
        out_cytoscape_mst_knn = output_dir + f"/{input_file}_out_mst_knn.xgmml"
        build_ssn.main([
            "-i", str(shared_datadir / input_file),
            "--xgmml", out_cytoscape_mst_knn,
            "--lb", "175",
            "--color_by", "SSN_cluster",
            "--cluster_tsv", out_clusters_mst_knn,
            "--metadata", metadata,
            "--mst_knn", "2",
        ])

        compare_files(out_clusters_mst_knn, out_clusters_full)

        with open(out_cytoscape_full, 'r') as f:
            full_edges = sum(1 for line in f if '<edge' in line)
        with open(out_cytoscape_mst, 'r') as f:
            mst_edges = sum(1 for line in f if '<edge' in line)
        with open(out_cytoscape_mst_knn, 'r') as f:
            mst_knn_edges = sum(1 for line in f if '<edge' in line)

        assert mst_edges < mst_knn_edges < full_edges


def test_build_ssn_mst_knn_preserves_ssn_cluster_colors(shared_datadir):
    input_file = "FeSOD_dist.tsv"
    expected_clusters = {"1", "2", "3"}

    with tempfile.TemporaryDirectory() as output_dir:
        metadata = str(shared_datadir / "FeSOD_metadata.tsv")

        out_clusters_baseline = output_dir + f"/{input_file}_out_clusters_baseline.tsv"
        out_cytoscape_baseline = output_dir + f"/{input_file}_out_baseline.xgmml"
        out_color_table_baseline = output_dir + "/color_table_baseline.tsv"
        build_ssn.main([
            "-i", str(shared_datadir / input_file),
            "--xgmml", out_cytoscape_baseline,
            "--lb", "175",
            "--color_by", "SSN_cluster",
            "--color_table_out", out_color_table_baseline,
            "--cluster_tsv", out_clusters_baseline,
            "--metadata", metadata,
        ])

        out_clusters_mst_knn = output_dir + f"/{input_file}_out_clusters_mst_knn_colors.tsv"
        out_cytoscape_mst_knn = output_dir + f"/{input_file}_out_mst_knn_colors.xgmml"
        out_color_table_mst_knn = output_dir + "/color_table_mst_knn.tsv"
        build_ssn.main([
            "-i", str(shared_datadir / input_file),
            "--xgmml", out_cytoscape_mst_knn,
            "--lb", "175",
            "--color_by", "SSN_cluster",
            "--color_table_out", out_color_table_mst_knn,
            "--cluster_tsv", out_clusters_mst_knn,
            "--metadata", metadata,
            "--mst_knn", "2",
        ])

        compare_files(out_clusters_baseline, out_clusters_mst_knn)

        baseline_color_table = {}
        with open(out_color_table_baseline, "r") as handle:
            for line in handle:
                cluster, color = line.strip().split("\t")
                baseline_color_table[cluster] = color

        mst_knn_color_table = {}
        with open(out_color_table_mst_knn, "r") as handle:
            for line in handle:
                cluster, color = line.strip().split("\t")
                mst_knn_color_table[cluster] = color

        assert baseline_color_table == mst_knn_color_table
        assert set(baseline_color_table.keys()) == expected_clusters
        assert all(re.match(r"#[0-9a-fA-F]{6}", color) for color in baseline_color_table.values())
        assert len(set(baseline_color_table.values())) == len(expected_clusters)


def _read_edges(xgmml_path):
    """Return [(mst_att_value_or_None, score), ...] for every edge in an xgmml file."""
    text = Path(xgmml_path).read_text()
    edges = []
    for edge in re.findall(r"<edge .*?</edge>", text, re.S):
        mst = re.search(r'name="MST" value="(\w+)"', edge)
        score = re.search(r'name="SSN_SCORE" value="([\d.]+)"', edge)
        edges.append((mst.group(1) if mst else None, score.group(1)))
    return edges


@pytest.mark.parametrize("extra_args,expect_non_mst", [
    ([], True),                  # full graph: MST backbone plus the rest
    (["--mst"], False),          # only MST edges are emitted, so all are marked true
    (["--mst_knn", "2"], True),  # MST plus kNN edges around it
])
def test_build_ssn_mark_mst(extra_args, expect_non_mst, shared_datadir):
    """--mark_mst tags every edge, and finds the same MST on all three edge paths."""
    with tempfile.TemporaryDirectory() as output_dir:
        out_cytoscape = output_dir + "/out.xgmml"
        build_ssn.main(["-i", str(shared_datadir / "FeSOD_dist.tsv"), "--xgmml", out_cytoscape,
                        "--lb", "175", "--metadata", str(shared_datadir / "FeSOD_metadata.tsv"),
                        "--mark_mst"] + extra_args)

        edges = _read_edges(out_cytoscape)
        assert len(edges) > 0
        assert all(mst in ("true", "false") for mst, _ in edges)
        # 20 nodes in 3 clusters, so the maximum spanning forest has 20 - 3 = 17 edges
        assert sum(1 for mst, _ in edges if mst == "true") == 17
        assert any(mst == "false" for mst, _ in edges) == expect_non_mst


def test_build_ssn_mark_mst_is_additive(shared_datadir):
    """--mark_mst only adds the attribute: clusters and edges are otherwise unchanged."""
    with tempfile.TemporaryDirectory() as output_dir:
        common = ["-i", str(shared_datadir / "FeSOD_dist.tsv"), "--lb", "175",
                  "--metadata", str(shared_datadir / "FeSOD_metadata.tsv"),
                  "--color_by", "SSN_cluster"]
        build_ssn.main(common + ["--xgmml", output_dir + "/off.xgmml",
                                 "--cluster_tsv", output_dir + "/off.tsv"])
        build_ssn.main(common + ["--xgmml", output_dir + "/on.xgmml",
                                 "--cluster_tsv", output_dir + "/on.tsv", "--mark_mst"])

        compare_files(output_dir + "/off.tsv", output_dir + "/on.tsv")
        off, on = _read_edges(output_dir + "/off.xgmml"), _read_edges(output_dir + "/on.xgmml")
        assert all(mst is None for mst, _ in off)
        assert [score for _, score in off] == [score for _, score in on]


def test_build_ssn_max_output_gb_blocks_xgmml(shared_datadir):
    input_file = "FeSOD_dist.tsv"
    with tempfile.TemporaryDirectory() as output_dir:
        output_path = Path(output_dir)
        metadata = str(shared_datadir / "FeSOD_metadata.tsv")
        out_clusters = output_dir + f"/{input_file}_out_clusters.tsv"
        out_cytoscape = output_dir + f"/{input_file}_out.xgmml"

        with pytest.raises(SystemExit, match="--max_output_gb"):
            build_ssn.main([
                "-i", str(shared_datadir / input_file),
                "--xgmml", out_cytoscape,
                "--lb", "175",
                "--color_by", "SSN_cluster",
                "--cluster_tsv", out_clusters,
                "--metadata", metadata,
                "--max_output_gb", "0.000001",
            ])

            assert not Path(out_cytoscape).exists()
            assert not Path(out_clusters).exists()
            assert not any(path.suffix == ".tmp" for path in output_path.iterdir())


@pytest.mark.parametrize("subset_mode", ["subset", "subset_file"])
def test_build_ssn_subset(shared_datadir, subset_mode):
    input_file = "FeSOD_dist.tsv"
    subset_labels = [
        "FeSOD_A0A1F4ZT98|unreviewed|Superoxide",
        "FeSOD_A0A067LT26|unreviewed|Superoxide",
        "FeSOD_A0A538G8K1|unreviewed|Superoxide",
        "FeSOD_B8LFE6|unreviewed|Superoxide",
    ]

    with tempfile.TemporaryDirectory() as output_dir:
        metadata = str(shared_datadir / "FeSOD_metadata.tsv")
        out_clusters = output_dir + f"/{input_file}_subset_clusters.tsv"
        out_cytoscape = output_dir + f"/{input_file}_subset.xgmml"

        subset_args = ["--subset", *subset_labels]
        if subset_mode == "subset_file":
            subset_file = Path(output_dir) / "subset.txt"
            subset_file.write_text("\n".join(subset_labels) + "\n")
            subset_args = ["--subset_file", str(subset_file)]

        build_ssn.main([
            "-i", str(shared_datadir / input_file),
            "--xgmml", out_cytoscape,
            "--lb", "175",
            "--color_by", "SSN_cluster",
            "--cluster_tsv", out_clusters,
            "--metadata", metadata,
            *subset_args,
        ])

        clusters = pd.read_csv(out_clusters, sep="\t")
        assert clusters["contig"].tolist() == subset_labels
        assert clusters["cluster"].tolist() == [1, 2, 1, 2]

        with open(out_cytoscape, "r") as handle:
            xgmml = handle.read()

        assert sum(1 for line in xgmml.splitlines() if "<node " in line) == 4
        assert sum(1 for line in xgmml.splitlines() if "<edge " in line) == 2
        assert "FeSOD_A0A2E1RF15|unreviewed|Superoxide" not in xgmml


def test_build_ssn_mst_knn_subset_after_filtering():
    labels = ["A", "B", "C", "D", "E"]
    matrix = DenseDataMatrix(np.array([
        [0, 100, 99, 98, 0],
        [100, 0, 97, 1, 96],
        [99, 97, 0, 95, 94],
        [98, 1, 95, 0, 200],
        [0, 96, 94, 200, 0],
    ]), labels, labels)

    subset_labels = {"A", "B", "C", "D"}
    subset_matrix = build_ssn.subset_matrix_by_labels(matrix, subset_labels)
    edge_dict = mst_knn_edge_index_dict(subset_matrix, 2, lower_bound=0)
    subset_edges = {
        tuple(sorted((subset_matrix.rows[source_idx], subset_matrix.rows[target_idx])))
        for source_idx, target_idx in edge_dict.keys()
    }

    assert subset_edges == {
        ("A", "B"),
        ("A", "C"),
        ("A", "D"),
        ("B", "C"),
        ("C", "D"),
    }


def test_cluster_labels_from_tree_matches_connected_components():
    labels = ["A", "B", "C", "D", "E"]
    array = np.array([
        [0, 9, 0, 0, 0],
        [9, 0, 5, 0, 0],
        [0, 5, 0, 0, 0],
        [0, 0, 0, 0, 7],
        [0, 0, 0, 7, 0],
    ], dtype=float)

    for matrix in (
        DenseDataMatrix(array, labels, labels),
        SparseDataMatrix(scipy.sparse.csr_array(array), labels, labels),
    ):
        tree = MaxTree(matrix)
        expected = connected_components(matrix.data > 5, directed=False, return_labels=True)[1]
        expected = build_ssn.rename_labels_by_frequency(expected) + 1
        observed = build_ssn.cluster_labels_from_tree(tree, 5)
        assert np.array_equal(observed, expected)


def test_iter_default_ssn_edges_uses_indices_and_lexical_endpoint_order():
    labels = ["zeta", "alpha", "beta"]
    matrix = DenseDataMatrix(np.array([
        [0, 10, 0],
        [10, 0, 7],
        [0, 7, 0],
    ], dtype=float), labels, labels)

    edges = list(build_ssn.iter_default_ssn_edges(matrix, lb=0))

    assert edges == [
        (1, 0, 10.0),
        (1, 2, 7.0),
    ]


def test_build_ssn_mst_knn_accepts_k_one(shared_datadir):
    # k == 1 reduces to the MST plus each node's single nearest neighbor.
    with tempfile.TemporaryDirectory() as output_dir:
        build_ssn.main([
            "-i", str(shared_datadir / "FeSOD_dist.tsv"),
            "--xgmml", str(Path(output_dir) / "out.xgmml"),
            "--mst_knn", "1",
        ])


def test_build_ssn_mst_knn_rejects_negative_k(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        with pytest.raises(SystemExit):
            build_ssn.main([
                "-i", str(shared_datadir / "FeSOD_dist.tsv"),
                "--xgmml", str(Path(output_dir) / "out.xgmml"),
                "--mst_knn", "-1",
            ])


def test_build_ssn_mst_knn_zero_matches_mst(shared_datadir):
    """--mst_knn 0 emits the same edges as --mst."""
    input_file = "FeSOD_dist.tsv"
    with tempfile.TemporaryDirectory() as output_dir:
        out_mst = output_dir + "/out_mst.xgmml"
        out_mst_knn_0 = output_dir + "/out_mst_knn_0.xgmml"

        for out_path, sparsify_args in ((out_mst, ["--mst"]), (out_mst_knn_0, ["--mst_knn", "0"])):
            build_ssn.main([
                "-i", str(shared_datadir / input_file),
                "--xgmml", out_path,
                "--lb", "175",
                *sparsify_args,
            ])

        def edge_lines(path):
            with open(path) as f:
                return sorted(line.strip() for line in f if '<edge' in line)

        mst_lines = edge_lines(out_mst)
        assert len(mst_lines) > 0
        assert edge_lines(out_mst_knn_0) == mst_lines


def _bridged_triangles_matrix():
    """Two triangles joined by a weak bridge that is nobody's nearest neighbor."""
    labels = list("abcdef")
    weights = {("a", "b"): 10.0, ("a", "c"): 9.0, ("b", "c"): 8.0,
               ("d", "e"): 7.0, ("d", "f"): 6.0, ("e", "f"): 5.0,
               ("c", "d"): 1.0}
    data = np.zeros((6, 6))
    for (left, right), value in weights.items():
        i, j = labels.index(left), labels.index(right)
        data[i, j] = data[j, i] = value
    return DenseDataMatrix(data, labels, labels, data_type="score"), labels


def _xgmml_edge_names(path, labels):
    """Return the sorted (source, target) label pairs of an xgmml network.

    Node ids are written as n<row index> in matrix row order by both the index-style and
    name-style edge writers, so the row labels translate them back.
    """
    pairs = re.findall(r'<edge[^>]*source="n(\d+)"[^>]*target="n(\d+)"', Path(path).read_text())
    return sorted(tuple(sorted((labels[int(source)], labels[int(target)]))) for source, target in pairs)


def test_build_ssn_knn_and_mst_combinations(shared_datadir):
    """--mst --knn K == --mst_knn K, and --knn alone drops the connecting bridge."""
    matrix, labels = _bridged_triangles_matrix()

    with tempfile.TemporaryDirectory() as output_dir:
        input_file = output_dir + "/input.hdf5"
        matrix.write(input_file, "dense")

        outputs = {}
        clusters = {}
        for tag, args in (("full", []), ("mst", ["--mst"]), ("knn1", ["--knn", "1"]),
                          ("knn2", ["--knn", "2"]), ("mst_knn1", ["--mst_knn", "1"]),
                          ("mst+knn1", ["--mst", "--knn", "1"]), ("mst_knn0", ["--mst_knn", "0"])):
            xgmml = output_dir + f"/{tag}.xgmml"
            cluster_tsv = output_dir + f"/{tag}.tsv"
            build_ssn.main(["-i", input_file, "--xgmml", xgmml, "--cluster_tsv", cluster_tsv] + args)
            outputs[tag] = _xgmml_edge_names(xgmml, labels)
            clusters[tag] = len({line.split("\t")[1] for line in
                                 Path(cluster_tsv).read_text().splitlines()[1:]})

        # the documented equivalences
        assert outputs["mst+knn1"] == outputs["mst_knn1"]
        assert outputs["mst_knn0"] == outputs["mst"]

        # a pure kNN graph keeps neighbors but loses the bridge, so it splits the cluster
        assert ("c", "d") in outputs["full"]
        assert ("c", "d") in outputs["mst_knn1"]
        assert ("c", "d") not in outputs["knn1"]
        assert ("c", "d") not in outputs["knn2"]
        assert clusters["full"] == 1
        assert clusters["mst"] == 1
        assert clusters["mst_knn1"] == 1
        assert clusters["knn1"] == 2
        assert clusters["knn2"] == 2

        # kNN edges are always a subset of the corresponding union
        assert set(outputs["knn1"]) < set(outputs["mst_knn1"])
        assert outputs["knn2"] == [("a", "b"), ("a", "c"), ("b", "c"), ("d", "e"), ("d", "f"), ("e", "f")]


def test_build_ssn_knn_marks_mst_edges():
    """--mark_mst still labels which of the kept kNN edges are MST edges."""
    matrix, labels = _bridged_triangles_matrix()
    with tempfile.TemporaryDirectory() as output_dir:
        input_file = output_dir + "/input.hdf5"
        matrix.write(input_file, "dense")
        xgmml = output_dir + "/knn_marked.xgmml"
        build_ssn.main(["-i", input_file, "--xgmml", xgmml, "--knn", "2", "--mark_mst"])

        text = Path(xgmml).read_text()
        assert 'name="MST"' in text
        # a, b, c triangle: a-b and a-c are MST edges, b-c is not
        assert text.count('name="MST" value="true"') == 4
        assert text.count('name="MST" value="false"') == 2


def test_build_ssn_knn_rejects_invalid_combinations():
    matrix, labels = _bridged_triangles_matrix()
    with tempfile.TemporaryDirectory() as output_dir:
        input_file = output_dir + "/input.hdf5"
        matrix.write(input_file, "dense")
        xgmml = output_dir + "/out.xgmml"

        # k must select at least one neighbor when there is no spanning tree to fall back on
        with pytest.raises(SystemExit):
            build_ssn.main(["-i", input_file, "--xgmml", xgmml, "--knn", "0"])

        # --mst_knn is shorthand for --mst --knn, so it cannot be combined with either
        with pytest.raises(SystemExit):
            build_ssn.main(["-i", input_file, "--xgmml", xgmml, "--knn", "2", "--mst_knn", "2"])
        with pytest.raises(ValueError):
            build_ssn.main(["-i", input_file, "--xgmml", xgmml, "--mst", "--mst_knn", "2"])
