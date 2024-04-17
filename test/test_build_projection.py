from domainator import build_projection
import tempfile
import pandas as pd
from pathlib import Path
import pytest

#TODO: test from coords file


@pytest.mark.parametrize("input_file",
[
    "scorefull.tsv",
    "scorefull.dense.hdf5"
])
def test_build_projection_1(input_file, shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "tmp_out"
        out_png = output_dir + "/out.png"
        out_coords = output_dir + "/out.tsv"
        out_cytoscape = output_dir + "/out.xgmml"
        build_projection.main(["-i", str(shared_datadir / input_file), "--algorithm", "pca", "--png_out", out_png, "--coords_out", out_coords, "--metadata",  str(shared_datadir / "metadata_FeSOD_20.tsv"), "--color_by", "group","--xgmml", out_cytoscape ])
        out_table = pd.read_csv(out_coords, sep="\t", index_col=0)
        out_mat = out_table.to_numpy()
        assert out_mat.shape == (20,2)
        assert Path(out_png).is_file()
        assert Path(out_cytoscape).is_file()
        #assert (out_mat.diagonal() == [410., 429., 405., 425., 411., 410., 400., 434., 426., 414., 398., 424., 430., 413., 402., 417., 451., 419., 429., 422.]).all()

@pytest.mark.parametrize("input_file",
[
    "scorefull.tsv",
    "scorefull.dense.hdf5"
])
def test_build_projection_umap(input_file, shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "tmp_out"
        out_coords = output_dir + "/out.tsv"
        build_projection.main(["-i", str(shared_datadir / input_file), "--algorithm", "umap", "--coords_out", out_coords])
        out_table = pd.read_csv(out_coords, sep="\t", index_col=0)
        out_mat = out_table.to_numpy()
        assert out_mat.shape == (20,2)

@pytest.mark.parametrize("input_file",
[
    "scorefull.tsv",
    "scorefull.dense.hdf5"
])
def test_build_projection_tsne(input_file, shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "tmp_out"
        out_coords = output_dir + "/out.tsv"
        build_projection.main(["-i", str(shared_datadir / input_file), "--algorithm", "tsne", "--coords_out", out_coords])
        out_table = pd.read_csv(out_coords, sep="\t", index_col=0)
        out_mat = out_table.to_numpy()
        assert out_mat.shape == (20,2)

@pytest.mark.parametrize("input_file",
[
    "scorefull.tsv",
    "scorefull.dense.hdf5"
])
def test_build_projection_3d(input_file, shared_datadir): 
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "tmp_out"
        out_coords = output_dir + "/out.tsv"
        out_cytoscape = output_dir + "/out.xgmml"
        build_projection.main(["-i", str(shared_datadir / input_file), "--3d", "--algorithm", "pca", "--coords_out", out_coords, "--metadata",  str(shared_datadir / "metadata_FeSOD_20.tsv"), "--color_by", "group","--xgmml", out_cytoscape ])
        out_table = pd.read_csv(out_coords, sep="\t", index_col=0)
        out_mat = out_table.to_numpy()
        assert out_mat.shape == (20,3)
        assert Path(out_cytoscape).is_file()

def test_build_projection_sparse_pca(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "tmp_out"
        out_sparse = output_dir + "/out_sparse.tsv"
        build_projection.main(["-i", str(shared_datadir / "bin3.sparse.hdf5"), "--algorithm", "pca", "--coords_out", out_sparse])
        out_sparse_table = pd.read_csv(out_sparse, sep="\t", index_col=0)
        out_mat = out_sparse_table.to_numpy()
        assert out_mat.shape == (20,2)

def test_build_projection_sparse_umap(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "tmp_out"
        out_sparse = output_dir + "/out_sparse.tsv"
        build_projection.main(["-i", str(shared_datadir / "bin3.sparse.hdf5"), "--algorithm", "umap", "--coords_out", out_sparse,  ])
        out_sparse_table = pd.read_csv(out_sparse, sep="\t", index_col=0)
        out_mat = out_sparse_table.to_numpy()
        assert out_mat.shape == (20,2)

def test_build_projection_sparse_tsne(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "tmp_out"
        out_sparse = output_dir + "/out_sparse.tsv"
        build_projection.main(["-i", str(shared_datadir / "bin3.sparse.hdf5"), "--algorithm", "tsne", "--coords_out", out_sparse])
        out_sparse_table = pd.read_csv(out_sparse, sep="\t", index_col=0)
        out_mat = out_sparse_table.to_numpy()
        assert out_mat.shape == (20,2)

def test_build_projection_color_table_out(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "tmp_out"
        out_sparse = output_dir + "/out_sparse.tsv"
        out_color_table = output_dir + "/color_table.tsv"
        build_projection.main(["-i", str(shared_datadir / "bin3.sparse.hdf5"), "--algorithm", "pca", "--coords_out", out_sparse, "--color_by", "active",
                               "--color_table_out", out_color_table, "--metadata", str(shared_datadir / "FeSOD_metadata.tsv")])
        out_sparse_table = pd.read_csv(out_sparse, sep="\t", index_col=0)
        out_mat = out_sparse_table.to_numpy()
        assert out_mat.shape == (20,2)
        assert Path(out_color_table).is_file()
        color_table_file = pd.read_csv(out_color_table, sep="\t", header=None)
        assert color_table_file.shape == (2,2)