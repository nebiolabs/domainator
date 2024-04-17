from domainator import utils, DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME
from pathlib import Path
import tempfile

from domainator import color_genbank

def test_color_genbank_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/tmp_colored_genbank.gb"
        
        color_genbank.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--color_table", str(shared_datadir / "color_specification.tsv"), "--color_both"])
        assert Path(out).is_file()
        cds_colors = {"1264_-1_959": "#FF0000", "2916_1_3677": "#00FF00", "2265_-1_1606": "#0000FF"}
        domain_colors = {"CcdB": "#FF0000", "APH": "#00FF00", "CAT": "#0000FF"}
        seen = 0
        for rec in utils.parse_seqfiles([out]):
            for feature in rec.features:
                if feature.type == "CDS":
                    if feature.qualifiers['cds_id'][0] in cds_colors:
                        seen += 1
                        assert feature.qualifiers['Color'][0] == cds_colors[feature.qualifiers['cds_id'][0]]
                elif feature.type == DOMAIN_FEATURE_NAME:
                    if feature.qualifiers['name'][0] in domain_colors:
                        assert feature.qualifiers['Color'][0] == domain_colors[feature.qualifiers['name'][0]]
        assert seen == 12


def test_color_genbank_2(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/tmp_colored_genbank.gb"
        
        color_genbank.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--color_both"])
        assert Path(out).is_file()
        cds_colors = {"1264_-1_959": "#3BA3EC", "2916_1_3677": "#B287F4", "2265_-1_1606": "#F77468"}
        domain_colors = {"CcdB": "#3BA3EC", "APH": "#B287F4", "CAT": "#F77468"}
        seen = 0
        for rec in utils.parse_seqfiles([out]):
            for feature in rec.features:
                if feature.type == "CDS":
                    if feature.qualifiers['cds_id'][0] in cds_colors:
                        seen += 1
                        assert feature.qualifiers['Color'][0] == cds_colors[feature.qualifiers['cds_id'][0]]
                elif feature.type == DOMAIN_FEATURE_NAME:
                    if feature.qualifiers['name'][0] in domain_colors:
                        assert feature.qualifiers['Color'][0] == domain_colors[feature.qualifiers['name'][0]]
        assert seen == 12

def test_color_genbank_3(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/tmp_colored_genbank.gb"
        
        color_genbank.main(["-i", str(shared_datadir / "pDONR_201_domainator_domain_reorder.gb"), "-o", out, "--color_both", "--color_table", str(shared_datadir / "color_specification.tsv")])
        assert Path(out).is_file()
        cds_colors = {"1605_-1_2265": "#0000FF"}
        domain_colors = {"CcdB": "#FF0000", "APH": "#00FF00", "CAT": "#0000FF", "Condensation":	"#FF00FF", "2-oxoacid_dh":	"#FFFFFF"}
        seen = 0
        for rec in utils.parse_seqfiles([out]):
            for feature in rec.features:
                if feature.type == "CDS":
                    if feature.qualifiers['cds_id'][0] in cds_colors:
                        seen += 1
                        assert feature.qualifiers['Color'][0] == cds_colors[feature.qualifiers['cds_id'][0]]
                elif feature.type == DOMAIN_FEATURE_NAME:
                    if feature.qualifiers['name'][0] in domain_colors:
                        assert feature.qualifiers['Color'][0] == domain_colors[feature.qualifiers['name'][0]]
        assert seen == 1

def test_color_genbank_pseudo_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/tmp_colored_genbank.gb"
        
        color_genbank.main(["-i", str(shared_datadir / "pDONR201_pseudo.gb"), "-o", out, "--color_table", str(shared_datadir / "color_specification.tsv"), "--color_both"])
        assert Path(out).is_file()
        # that's the only assertion we need because we're just testing that it doesn't crash.

def test_color_genbank_domain_search_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/tmp_colored_genbank.gb"
        
        color_genbank.main(["-i", str(shared_datadir / "color_domain_search_test.gb"), "-o", out, "--color_both", "--color_table", str(shared_datadir / "color_specification.tsv"), "--search_hit_color", "#555500"])
        assert Path(out).is_file()
        cds_colors = {"1605_-1_2265": "#0000FF", "2265_-1_1606": "#555500", "1264_-1_959":"#555500"}
        domain_colors = {"CcdB": "#FF0000", "APH": "#00FF00", "CAT": "#0000FF", "Condensation":	"#FF00FF", "2-oxoacid_dh":	"#FFFFFF"}
        for rec in utils.parse_seqfiles([out]):
            for feature in rec.features:
                if feature.type == "CDS":
                    if feature.qualifiers['cds_id'][0] in cds_colors:
                        assert feature.qualifiers['Color'][0] == cds_colors[feature.qualifiers['cds_id'][0]]
                elif feature.type == DOMAIN_FEATURE_NAME:
                    if feature.qualifiers['name'][0] in domain_colors:
                        assert feature.qualifiers['Color'][0] == domain_colors[feature.qualifiers['name'][0]]
                elif feature.type == DOMAIN_SEARCH_BEST_HIT_NAME:
                    assert feature.qualifiers['Color'][0] == "#555500"

def test_color_genbank_clear_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/tmp_colored_genbank.gb"
        
        color_genbank.main(["-i", str(shared_datadir / "MT_nbs.gb"), "-o", out, "--clear", "--color_domains"])
        assert Path(out).is_file()
        for rec in utils.parse_seqfiles([out]):
            for feature in rec.features:
                if feature.type == "CDS":
                        assert not ('Color' in feature.qualifiers)
                elif feature.type == DOMAIN_SEARCH_BEST_HIT_NAME:
                    assert not ('Color' in feature.qualifiers)


def test_color_genbank_color_table_out_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/tmp_colored_genbank.gb"
        color_table_out = output_dir + "/color_table.tsv"
        color_genbank.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--color_table", str(shared_datadir / "color_specification.tsv"), "--color_both", "--color_table_out", color_table_out])
        assert Path(out).is_file()
        color_table_dict = {}
        with open(color_table_out, "r") as f:
            for line in f:
                domain, color = line.strip().split("\t")
                color_table_dict[domain] = color
        assert color_table_dict["CcdB"] == "#FF0000"
        assert color_table_dict["APH"] == "#00FF00"
        assert color_table_dict["CAT"] == "#0000FF"

def test_color_genbank_color_table_out_2(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/tmp_colored_genbank.gb"
        color_table_out = output_dir + "/color_table.tsv"
        color_genbank.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--color_both", "--color_table_out", color_table_out])
        assert Path(out).is_file()
        color_table_dict = {}
        
        with open(color_table_out, "r") as f:
            for line in f:
                domain, color = line.strip().split("\t")
                color_table_dict[domain] = color
        assert len(color_table_dict) == 7
        assert color_table_dict["CcdB"] == "#3BA3EC"
        assert color_table_dict["APH"] == "#B287F4"
        assert color_table_dict["CAT"] == "#F77468"
        assert len(set(color_table_dict.values())) == 7

