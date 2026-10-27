import pandas as pd
import pytest
from domainator import plot_contigs
from domainator.colors import DISTINCT_COLORS, get_palette
import tempfile

def test_plot_contigs_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        input = str(shared_datadir / "MT_nbs.gb")
        output = output_dir + "/contigs.html"
        plot_contigs.main(["-i", input, "--html_out", output])
        output_text = open(output).read()
        assert "<title>Domainator Contigs Plot</title>" in output_text
        assert "BX548174_369054:361090rc" in output_text
        assert not ('"type": "source"' in output_text)

def test_plot_contigs_2(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        input = str(shared_datadir / "MT_nbs.gb")
        output = output_dir + "/contigs.html"
        plot_contigs.main(["-i", input, "--html_out", output, "--height", "1000", "--width", "800"])
        output_text = open(output).read()
        assert '<svg width="800" height="1000' in output_text
        assert "BX548174_369054:361090rc" in output_text
        assert not ('"type": "source"' in output_text)

def test_plot_contigs_3(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        input = str(shared_datadir / "MT_nbs.gb")
        output = output_dir + "/contigs.html"
        plot_contigs.main(["-i", input, "--html_out", output, "--height", "1000", "--width", "800", "--types", "source"])
        output_text = open(output).read()
        assert '<svg width="800" height="1000' in output_text
        assert "BX548174_369054:361090rc" in output_text
        assert '"type": "source"' in output_text

def test_plot_contigs_4(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        #plot_contigs.py -i test/data/MT_nbs.gb --metadata MT_nbs.enum_report.tsv --color_by taxid_species --html out.html
        input = str(shared_datadir / "MT_nbs.gb")
        output = output_dir + "/contigs.html"
        metadata = str(shared_datadir / "MT_nbs.enum_report.tsv")
        plot_contigs.main(["-i", input, "--html_out", output, "--metadata", metadata, "--color_by", "taxid_species"])
        output_text = open(output).read()
        assert "BX548174_369054:361090rc" in output_text
        assert not ('"type": "source"' in output_text)
        assert '<b>taxid_species</b>: <span style=\\"color:#' in output_text

def test_plot_contigs_protein_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        input = str(shared_datadir / "FeSOD_20_pfam.gb")
        output = output_dir + "/contigs.html"
        plot_contigs.main(["-i", input, "--html_out", output])
        output_text = open(output).read()
        assert "<title>Domainator Contigs Plot</title>" in output_text
        assert "FeSOD_A0A1F4ZT98|unreviewed|Superoxide" in output_text
        assert not ('"type": "source"' in output_text)

def test_plot_contigs_color_table(shared_datadir, tmp_path):
    color_table = tmp_path / "colors.tsv"
    color_table.write_text("1219\tAB12CD\n")
    output = str(tmp_path / "contigs.html")
    plot_contigs.main(["-i", str(shared_datadir / "MT_nbs.gb"), "--html_out", output,
                       "--metadata", str(shared_datadir / "MT_nbs.enum_report.tsv"), "--color_by", "taxid_species",
                       "--color_table", str(color_table)])
    output_text = open(output).read()
    assert '"color": "#AB12CD"' in output_text  # taxid 1219 comes from the table
    assert '"color": "#000000"' in output_text  # taxids missing from the table keep the default


def test_plot_contigs_palette_offset(shared_datadir, tmp_path):
    """--palette_offset skips that many palette colors before coloring the values."""
    metadata = shared_datadir / "MT_nbs.enum_report.tsv"
    output = tmp_path / "contigs.html"
    plot_contigs.main(["-i", str(shared_datadir / "MT_nbs.gb"), "--html_out", str(output),
                       "--metadata", str(metadata), "--color_by", "taxid_species", "--palette_offset", "3"])
    output_text = output.read_text()
    species = pd.read_csv(metadata, sep="\t", index_col=0)["taxid_species"].unique()
    assert len(species) == 7
    for color in get_palette(species, offset=3).values():
        assert f'"color": "{color}"' in output_text
    # The first three colors are the ones skipped.
    for color in DISTINCT_COLORS[:3]:
        assert f'"color": "{color}"' not in output_text


def test_plot_contigs_color_table_requires_color_by(shared_datadir, tmp_path, capsys):
    color_table = tmp_path / "colors.tsv"
    color_table.write_text("1219\tAB12CD\n")
    with pytest.raises(SystemExit):
        plot_contigs.main(["-i", str(shared_datadir / "MT_nbs.gb"), "--html_out", str(tmp_path / "out.html"),
                           "--color_table", str(color_table)])
    assert "--color_table requires --color_by" in capsys.readouterr().err


def test_feature_to_dict_accepts_one_hex_default_color():
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from domainator.plot_contigs import feature_to_dict
    feature = SeqFeature(FeatureLocation(0, 30, 1), type="misc_feature")
    assert feature_to_dict(feature, default_color="CC00CC")["color"] == "#CC00CC"
