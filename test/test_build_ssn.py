from domainator import build_ssn
from domainator.data_matrix import DenseDataMatrix, mst_knn_edge_index_dict
import pytest
import tempfile
import pandas as pd
import numpy as np
from pathlib import Path
from helpers import compare_files
import re

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


def test_build_ssn_mst_knn_requires_k_gt_one(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        with pytest.raises(SystemExit):
            build_ssn.main([
                "-i", str(shared_datadir / "FeSOD_dist.tsv"),
                "--xgmml", str(Path(output_dir) / "out.xgmml"),
                "--mst_knn", "1",
            ])
