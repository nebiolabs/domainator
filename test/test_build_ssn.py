from domainator import build_ssn
import pytest
import tempfile
import pandas as pd
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
