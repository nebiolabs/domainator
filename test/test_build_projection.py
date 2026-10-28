from domainator import build_projection
from domainator.colors import DISTINCT_COLORS, OTHER_COLOR
import io
import tempfile
import pandas as pd
from pathlib import Path
import pytest
from helpers import report_data


def _page_data(html):
    envelope = report_data(html)
    assert envelope["report"] == "build_projection"
    return envelope["data"]


@pytest.mark.parametrize("input_file",
[
    "scorefull.tsv",
    "scorefull.dense.hdf5"
])
def test_build_projection_1(input_file, shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "tmp_out"
        out_html = output_dir + "/out.html"
        out_coords = output_dir + "/out.tsv"
        out_cytoscape = output_dir + "/out.xgmml"
        out_color_table = output_dir + "/colors.tsv"
        build_projection.main(["-i", str(shared_datadir / input_file), "--algorithm", "pca", "--html_out", out_html, "--coords_out", out_coords, "--metadata",  str(shared_datadir / "metadata_FeSOD_20.tsv"), "--color_by", "group","--xgmml_out", out_cytoscape, "--color_table_out", out_color_table])
        out_table = pd.read_csv(out_coords, sep="\t", index_col=0)
        out_mat = out_table.to_numpy()
        assert out_mat.shape == (20,2)
        assert Path(out_cytoscape).is_file()
        #assert (out_mat.diagonal() == [410., 429., 405., 425., 411., 410., 400., 434., 426., 414., 398., 424., 430., 413., 402., 417., 451., 419., 429., 422.]).all()

        # The page has the same points, coordinates and colors as the other outputs.
        data = _page_data(Path(out_html).read_text())
        assert data["algorithm"] == "pca"
        assert data["axes"] == ["projection_X", "projection_Y"]
        assert data["ids"] == list(out_table.index)
        assert data["coords"] == [out_table["projection_X"].tolist(), out_table["projection_Y"].tolist()]
        assert data["color_by"] == "group"
        color_table = dict(line.split("\t") for line in Path(out_color_table).read_text().splitlines())
        assert [group["value"] for group in data["groups"]] == ["g1", "g2", "g3"]
        assert {group["value"]: group["color"] for group in data["groups"]} == color_table
        metadata = pd.read_csv(shared_datadir / "metadata_FeSOD_20.tsv", sep="\t", index_col=0)["group"]
        assert [data["groups"][g]["value"] for g in data["group"]] == [metadata[name] for name in data["ids"]]
        assert [group["count"] for group in data["groups"]] == [data["group"].count(g) for g in range(3)]
        assert data["metadata"] == {"columns": ["group"], "values": [[metadata[name] for name in data["ids"]]]}


def test_png_out_is_gone(shared_datadir, tmp_path):
    """--png_out was replaced by --html_out, whose page saves PNG and SVG pictures."""
    with pytest.raises(SystemExit):
        build_projection.main(["-i", str(shared_datadir / "scorefull.tsv"), "--png_out", str(tmp_path / "out.png")])


def test_build_projection_html_3d(shared_datadir):
    page = io.StringIO()
    build_projection.build_projection(build_projection.DataMatrix.from_file(shared_datadir / "scorefull.tsv"),
                                      three_d=True, html_out=page)
    data = _page_data(page.getvalue())
    assert data["axes"] == ["projection_X", "projection_Y", "projection_Z"]
    assert len(data["coords"]) == 3 and all(len(axis) == 20 for axis in data["coords"])
    assert data["color_by"] is None and data["groups"] is None and data["group"] is None
    assert data["metadata"] == {"columns": [], "values": []}


def test_build_projection_html_colors(shared_datadir, tmp_path):
    """Values a --color_table does not list, and points with no value, are gray; table order
    comes first and the points with no value last."""
    names = pd.read_csv(shared_datadir / "metadata_FeSOD_20.tsv", sep="\t", index_col=0).index
    metadata = tmp_path / "metadata.tsv"
    # Every fifth point has no value; the others are a, b or c.
    metadata.write_text("name\tkind\tsize\n" + "".join(
        f"{name}\t{'' if i % 5 == 0 else 'abc'[i % 3]}\t{i}\n" for i, name in enumerate(names)))
    color_table = tmp_path / "colors.tsv"
    color_table.write_text("c\t#123456\na\t#ABCDEF\n")
    out_html = tmp_path / "out.html"
    build_projection.main(["-i", str(shared_datadir / "scorefull.tsv"), "--html_out", str(out_html),
                           "--metadata", str(metadata), "--color_by", "kind", "--color_table", str(color_table)])
    data = _page_data(out_html.read_text())
    assert [(group["value"], group["color"]) for group in data["groups"]] == [
        ("c", "#123456"), ("a", "#ABCDEF"), ("b", OTHER_COLOR), (None, OTHER_COLOR)]
    assert data["groups"][-1]["count"] == 4
    kinds = dict(zip(data["ids"], data["metadata"]["values"][0]))
    assert [data["groups"][g]["value"] for g in data["group"]] == [kinds[name] for name in data["ids"]]
    assert data["metadata"]["columns"] == ["kind", "size"]
    assert all(isinstance(value, int) for value in data["metadata"]["values"][1])


@pytest.mark.parametrize("input_file",
[
    "scorefull.tsv",
    "scorefull.dense.hdf5"
])
def test_build_projection_umap(input_file, shared_datadir):
    pytest.importorskip("umap", reason="needs the projection extra (umap-learn)")
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
        build_projection.main(["-i", str(shared_datadir / input_file), "--three_d", "--algorithm", "pca", "--coords_out", out_coords, "--metadata",  str(shared_datadir / "metadata_FeSOD_20.tsv"), "--color_by", "group","--xgmml_out", out_cytoscape ])
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
    pytest.importorskip("umap", reason="needs the projection extra (umap-learn)")
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


def test_build_projection_palette_offset(shared_datadir, tmp_path):
    """--palette_offset skips that many palette colors before coloring the values."""
    out_color_table = tmp_path / "color_table.tsv"
    build_projection.main(["-i", str(shared_datadir / "bin3.sparse.hdf5"), "--algorithm", "pca",
                           "--coords_out", str(tmp_path / "coords.tsv"), "--color_by", "active",
                           "--color_table_out", str(out_color_table), "--palette_offset", "5",
                           "--metadata", str(shared_datadir / "FeSOD_metadata.tsv")])
    color_table = dict(line.split("\t") for line in out_color_table.read_text().splitlines())
    assert color_table == {"0": DISTINCT_COLORS[5], "1": DISTINCT_COLORS[6]}


def test_build_projection_from_coords_file(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        first_coords = output_dir + "/first.tsv"
        second_coords = output_dir + "/second.tsv"

        build_projection.main([
            "-i", str(shared_datadir / "scorefull.tsv"),
            "--algorithm", "pca",
            "--coords_out", first_coords,
        ])
        build_projection.main([
            "-i", first_coords,
            "--algorithm", "coords",
            "--coords_out", second_coords,
        ])

        first_table = pd.read_csv(first_coords, sep="\t", index_col=0)
        second_table = pd.read_csv(second_coords, sep="\t", index_col=0)
        pd.testing.assert_frame_equal(first_table, second_table)

def test_build_projection_coords_rejects_params(shared_datadir, tmp_path):
    first_coords = str(tmp_path / "first.tsv")
    build_projection.main(["-i", str(shared_datadir / "scorefull.tsv"), "--algorithm", "pca", "--coords_out", first_coords])
    with pytest.raises(SystemExit):
        build_projection.main(["-i", first_coords, "--algorithm", "coords", "--coords_out", str(tmp_path / "second.tsv"),
                               "--tsne_args=perplexity=5"])
