import pytest
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature

from domainator import plot_contigs
from domainator.colors import DISTINCT_COLORS, OTHER_COLOR, get_palette
from domainator.errors import InputError
from helpers import report_data


def _page_data(path):
    envelope = report_data(open(path).read())
    assert envelope["report"] == "plot_contigs"
    assert envelope["encoding"] == "gzip+base64"
    return envelope["data"]


def _contig(data, name):
    return next(contig for contig in data["contigs"] if contig["name"] == name)


def _details(contig):
    return dict(contig["details"])


def test_plot_contigs_1(shared_datadir, tmp_path):
    output = str(tmp_path / "contigs.html")
    plot_contigs.main(["-i", str(shared_datadir / "MT_nbs.gb"), "--html_out", output])
    assert "<title>Contig Plot</title>" in open(output).read()
    data = _page_data(output)
    assert data["types"] == plot_contigs.DEFAULT_TYPES
    assert data["color_by"] is None
    assert data["layout"]["width"] == 1200 and data["layout"]["height"] == 2400
    assert data["layout"]["details_width"] == 1200
    contig = _contig(data, "BX548174_369054:361090rc")
    assert contig["length"] == 7963
    assert contig["color"] == plot_contigs.DEFAULT_CONTIG_COLOR
    assert {feature["type"] for contig in data["contigs"] for feature in contig["features"]} <= set(plot_contigs.DEFAULT_TYPES)
    details = _details(contig)
    assert details["Name"] == "BX548174_369054:361090rc"
    assert details["Length"] == "7963"
    assert details["taxonomy"].startswith("Bacteria; Cyanobacteriota")
    assert details["BioProject"] == "PRJNA213"
    assert "keywords" not in details  # Biopython's [""] has no value
    first = contig["features"][0]
    assert first["name"] == "PMM0387"
    assert first["parts"] == [[0, 441, "-"]]
    assert ["gene", ["aroD"]] in first["qualifiers"]
    assert first["color"] == "#AE9D31"  # its Color qualifier


def test_plot_contigs_2(shared_datadir, tmp_path):
    output = str(tmp_path / "contigs.html")
    plot_contigs.main(["-i", str(shared_datadir / "MT_nbs.gb"), "--html_out", output, "--height", "1000", "--width", "800"])
    data = _page_data(output)
    assert data["layout"]["width"] == 800 and data["layout"]["height"] == 1000
    assert data["layout"]["details_width"] == 800
    assert any(contig["name"] == "BX548174_369054:361090rc" for contig in data["contigs"])


def test_plot_contigs_3(shared_datadir, tmp_path):
    output = str(tmp_path / "contigs.html")
    plot_contigs.main(["-i", str(shared_datadir / "MT_nbs.gb"), "--html_out", output, "--height", "1000", "--width", "800", "--types", "source"])
    data = _page_data(output)
    assert data["types"] == ["source"]
    for contig in data["contigs"]:
        assert [feature["type"] for feature in contig["features"]] == ["source"]
    # A source feature is named by its organism.
    assert _contig(data, "BX548174_369054:361090rc")["features"][0]["name"] == "Prochlorococcus marinus subsp. pastoris str. CCMP1986"


def test_plot_contigs_4(shared_datadir, tmp_path):
    output = str(tmp_path / "contigs.html")
    plot_contigs.main(["-i", str(shared_datadir / "MT_nbs.gb"), "--html_out", output, "--metadata", str(shared_datadir / "MT_nbs.enum_report.tsv"), "--color_by", "taxid_species"])
    data = _page_data(output)
    assert data["color_by"] == "taxid_species"
    species = [_details(contig)["taxid_species"] for contig in data["contigs"]]
    palette = get_palette(list(dict.fromkeys(species)))
    assert [contig["color"] for contig in data["contigs"]] == [palette[value] for value in species]
    # The metadata follows the record's own details, as text.
    contig = _contig(data, "BX548174_369054:361090rc")
    assert contig["details"][-3:] == [["cds_count", "9"], ["domain_count", "6"], ["taxid_species", "1219"]]


def test_plot_contigs_protein_1(shared_datadir, tmp_path):
    output = str(tmp_path / "contigs.html")
    plot_contigs.main(["-i", str(shared_datadir / "FeSOD_20_pfam.gb"), "--html_out", output])
    data = _page_data(output)
    assert data["contigs"][0]["name"].startswith("FeSOD_A0A1F4ZT98|unreviewed|Superoxide")
    features = data["contigs"][0]["features"]
    assert [feature["type"] for feature in features] == ["Domainator", "Domainator"]
    assert {feature["name"] for feature in features} == {"Sod_Fe_N", "Sod_Fe_C"}
    assert all(feature["color"] == "#00FFFF" for feature in features)


def test_plot_contigs_color_table(shared_datadir, tmp_path):
    color_table = tmp_path / "colors.tsv"
    color_table.write_text("1219\tAB12CD\n")
    output = str(tmp_path / "contigs.html")
    plot_contigs.main(["-i", str(shared_datadir / "MT_nbs.gb"), "--html_out", output, "--metadata", str(shared_datadir / "MT_nbs.enum_report.tsv"), "--color_by", "taxid_species", "--color_table", str(color_table)])
    data = _page_data(output)
    for contig in data["contigs"]:
        # taxid 1219 comes from the table; taxids missing from it keep the default.
        expected = "#AB12CD" if _details(contig)["taxid_species"] == "1219" else plot_contigs.DEFAULT_CONTIG_COLOR
        assert contig["color"] == expected
    assert {contig["color"] for contig in data["contigs"]} == {"#AB12CD", plot_contigs.DEFAULT_CONTIG_COLOR}


def test_plot_contigs_palette_offset(shared_datadir, tmp_path):
    """--palette_offset skips that many palette colors before coloring the values."""
    output = str(tmp_path / "contigs.html")
    plot_contigs.main(["-i", str(shared_datadir / "MT_nbs.gb"), "--html_out", output, "--metadata", str(shared_datadir / "MT_nbs.enum_report.tsv"), "--color_by", "taxid_species", "--palette_offset", "3"])
    data = _page_data(output)
    species = list(dict.fromkeys(_details(contig)["taxid_species"] for contig in data["contigs"]))
    assert len(species) == 7
    colors = {contig["color"] for contig in data["contigs"]}
    assert colors == set(get_palette(species, offset=3).values())
    # The first three colors are the ones skipped.
    assert colors.isdisjoint(DISTINCT_COLORS[:3])


def test_plot_contigs_color_table_requires_color_by(shared_datadir, tmp_path, capsys):
    color_table = tmp_path / "colors.tsv"
    color_table.write_text("1219\tAB12CD\n")
    with pytest.raises(SystemExit):
        plot_contigs.main(["-i", str(shared_datadir / "MT_nbs.gb"), "--html_out", str(tmp_path / "out.html"), "--color_table", str(color_table)])
    assert "--color_table requires --color_by" in capsys.readouterr().err


def test_plot_contigs_color_by_requires_metadata(shared_datadir, tmp_path, capsys):
    with pytest.raises(SystemExit):
        plot_contigs.main(["-i", str(shared_datadir / "MT_nbs.gb"), "--html_out", str(tmp_path / "out.html"), "--color_by", "taxid_species"])
    assert "--color_by requires --metadata" in capsys.readouterr().err


def test_plot_contigs_color_by_must_be_a_metadata_column(shared_datadir, tmp_path):
    from domainator.seqio import parse_seqfiles
    records = list(parse_seqfiles([str(shared_datadir / "MT_nbs.gb")]))
    with pytest.raises(InputError, match="not a column"):
        plot_contigs.plot_contigs(records, html_out=str(tmp_path / "out.html"), metadata=[shared_datadir / "MT_nbs.enum_report.tsv"], color_by="species")


def test_contigs_without_a_value_are_drawn_in_the_other_color(shared_datadir, tmp_path):
    from domainator.seqio import parse_seqfiles
    records = list(parse_seqfiles([str(shared_datadir / "MT_nbs.gb")]))
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(f"contig\tgroup\n{records[0].id}\tA\n{records[1].id}\t\n")
    output = tmp_path / "out.html"
    plot_contigs.plot_contigs(records, html_out=str(output), metadata=[metadata], color_by="group")
    data = _page_data(output)
    assert data["contigs"][0]["color"] == DISTINCT_COLORS[0]
    assert {contig["color"] for contig in data["contigs"][1:]} == {OTHER_COLOR}
    assert "group" not in _details(data["contigs"][1])


def test_read_metadata_keeps_text_and_the_first_value(tmp_path):
    first = tmp_path / "first.tsv"
    first.write_text("id\tcount\tgroup\na\t007\t\na\t8\tlater row\nb\t1\tB\n")
    second = tmp_path / "second.tsv"
    second.write_text("name\tgroup\textra\na\tA\tx\nb\tnot this\t\n")
    columns, values = plot_contigs.read_metadata([first, second])
    assert columns == ["count", "group", "extra"]
    assert values == {"a": {"count": "007", "group": "A", "extra": "x"}, "b": {"count": "1", "group": "B"}}


def test_feature_to_dict_accepts_one_hex_default_color():
    feature = SeqFeature(FeatureLocation(0, 30, 1), type="misc_feature")
    assert plot_contigs.feature_to_dict(feature, default_color="CC00CC")["color"] == "#CC00CC"


def test_feature_to_dict_keeps_each_part():
    location = CompoundLocation([FeatureLocation(50, 80, -1), FeatureLocation(10, 30, -1)])
    feature = SeqFeature(location, type="CDS", qualifiers={"locus_tag": ["X1"], "note": ["a", "b"]})
    assert plot_contigs.feature_to_dict(feature) == {
        "name": "X1", "type": "CDS", "color": "#CCCCCC",
        "parts": [[50, 80, "-"], [10, 30, "-"]],
        "qualifiers": [["locus_tag", ["X1"]], ["note", ["a", "b"]]],
    }
