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
        output_dir = "test_out"
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
