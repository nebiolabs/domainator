from collections import OrderedDict
from domainator import color_table_to_legend

def test_color_table_legend(shared_datadir, tmp_path):
    svg_file = tmp_path / "test.svg"
    title = "Test Legend"

    color_table_to_legend.main(["-i", str(shared_datadir / "color_specification.tsv"), "--svg", str(svg_file), "--title", title])
    """
    CcdB	#ff0000
    APH	#00ff00
    CAT	#0000ff
    Condensation	#ff00ff
    2-oxoacid_dh	#ffffff

    """

    with open(svg_file, "r") as f:
        # read entire file into string
        text = f.read()

    assert "Test Legend" in text
    assert "#FF0000" in text
    assert "#00FF00" in text
    assert "#0000FF" in text
    assert "#FF00FF" in text
    assert "#FFFFFF" in text
    assert "CcdB" in text
    assert "APH" in text
    assert "CAT" in text
    assert "Condensation" in text
    assert "2-oxoacid_dh" in text

