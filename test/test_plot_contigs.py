from domainator import plot_contigs
import tempfile

def test_plot_contigs_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        input = str(shared_datadir / "MT_nbs.gb")
        output = output_dir + "/contigs.html"
        plot_contigs.main(["-i", input, "--html", output])
        output_text = open(output).read()
        assert "<title>Domainator Contigs Plot</title>" in output_text
        assert "BX548174_369054:361090rc" in output_text
        assert not ('"type": "source"' in output_text)

def test_plot_contigs_2(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        input = str(shared_datadir / "MT_nbs.gb")
        output = output_dir + "/contigs.html"
        plot_contigs.main(["-i", input, "--html", output, "--height", "1000", "--width", "800"])
        output_text = open(output).read()
        assert '<svg width="800" height="1000' in output_text
        assert "BX548174_369054:361090rc" in output_text
        assert not ('"type": "source"' in output_text)

def test_plot_contigs_3(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        input = str(shared_datadir / "MT_nbs.gb")
        output = output_dir + "/contigs.html"
        plot_contigs.main(["-i", input, "--html", output, "--height", "1000", "--width", "800", "--types", "source"])
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
        plot_contigs.main(["-i", input, "--html", output, "--metadata", metadata, "--color_by", "taxid_species"])
        output_text = open(output).read()
        assert "BX548174_369054:361090rc" in output_text
        assert not ('"type": "source"' in output_text)
        assert '<b>taxid_species</b>: <span style=\\"color:#' in output_text

def test_plot_contigs_protein_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        input = str(shared_datadir / "FeSOD_20_pfam.gb")
        output = output_dir + "/contigs.html"
        plot_contigs.main(["-i", input, "--html", output])
        output_text = open(output).read()
        assert "<title>Domainator Contigs Plot</title>" in output_text
        assert "FeSOD_A0A1F4ZT98|unreviewed|Superoxide" in output_text
        assert not ('"type": "source"' in output_text)