from domainator import DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME, SYSTEM_FEATURE_NAME
from domainator.seqio import parse_seqfiles
import os
from pathlib import Path
import tempfile

from domainator import color_genbank, find_systems

def test_color_genbank_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/tmp_colored_genbank.gb"
        
        color_genbank.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--color_table", str(shared_datadir / "color_specification.tsv"), "--target", "domains_and_cdss"])
        assert Path(out).is_file()
        cds_colors = {"1264_-1_959": "#FF0000", "2916_1_3677": "#00FF00", "2265_-1_1606": "#0000FF"}
        domain_colors = {"CcdB": "#FF0000", "APH": "#00FF00", "CAT": "#0000FF"}
        seen = 0
        for rec in parse_seqfiles([out]):
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
        
        color_genbank.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--target", "domains_and_cdss"])
        assert Path(out).is_file()
        cds_colors = {"1264_-1_959": "#3BA3EC", "2916_1_3677": "#B287F4", "2265_-1_1606": "#F77468"}
        domain_colors = {"CcdB": "#3BA3EC", "APH": "#B287F4", "CAT": "#F77468"}
        seen = 0
        for rec in parse_seqfiles([out]):
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
        
        color_genbank.main(["-i", str(shared_datadir / "pDONR_201_domainator_domain_reorder.gb"), "-o", out, "--target", "domains_and_cdss", "--color_table", str(shared_datadir / "color_specification.tsv")])
        assert Path(out).is_file()
        cds_colors = {"1605_-1_2265": "#0000FF"}
        domain_colors = {"CcdB": "#FF0000", "APH": "#00FF00", "CAT": "#0000FF", "Condensation":	"#FF00FF", "2-oxoacid_dh":	"#FFFFFF"}
        seen = 0
        for rec in parse_seqfiles([out]):
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
        
        color_genbank.main(["-i", str(shared_datadir / "pDONR201_pseudo.gb"), "-o", out, "--color_table", str(shared_datadir / "color_specification.tsv"), "--target", "domains_and_cdss"])
        assert Path(out).is_file()
        # that's the only assertion we need because we're just testing that it doesn't crash.

def test_color_genbank_domain_search_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/tmp_colored_genbank.gb"
        
        color_genbank.main(["-i", str(shared_datadir / "color_domain_search_test.gb"), "-o", out, "--target", "domains_and_cdss", "--color_table", str(shared_datadir / "color_specification.tsv"), "--search_hit_color", "#555500"])
        assert Path(out).is_file()
        cds_colors = {"1605_-1_2265": "#0000FF", "2265_-1_1606": "#555500", "1264_-1_959":"#555500"}
        domain_colors = {"CcdB": "#FF0000", "APH": "#00FF00", "CAT": "#0000FF", "Condensation":	"#FF00FF", "2-oxoacid_dh":	"#FFFFFF"}
        for rec in parse_seqfiles([out]):
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
        
        color_genbank.main(["-i", str(shared_datadir / "MT_nbs.gb"), "-o", out, "--clear", "--target", "domains"])
        assert Path(out).is_file()
        for rec in parse_seqfiles([out]):
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
        color_genbank.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--color_table", str(shared_datadir / "color_specification.tsv"), "--target", "domains_and_cdss", "--color_table_out", color_table_out])
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
        color_genbank.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--target", "domains_and_cdss", "--color_table_out", color_table_out])
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

def test_color_genbank_databases_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/tmp_colored_genbank.gb"
        
        color_genbank.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator_multi_hmm_2.gb"), "-o", out, "--color_table", str(shared_datadir / "color_specification.tsv"), "--target", "domains", "--databases", "pdonr_hmms_1"])
        assert Path(out).is_file()
        cds_colors = {"1264_-1_959": "#FF0000", "2916_1_3677": "#00FF00", "2265_-1_1606": "#0000FF"}
        domain_colors = {"CcdB": "#FF0000", "APH": "#00FF00", "CAT": "#0000FF"}
        seen = 0
        for rec in parse_seqfiles([out]):
            for feature in rec.features:
                if feature.type == DOMAIN_FEATURE_NAME:
                    if feature.qualifiers['database'][0] == "pdonr_hmms_1":
                        if feature.qualifiers['name'][0] in domain_colors:
                            seen += 1
                            assert feature.qualifiers['Color'][0] == domain_colors[feature.qualifiers['name'][0]]
                    else:
                        assert not ('Color' in feature.qualifiers)
        assert seen == 6

def test_color_genbank_mode_domains_and_cdss(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/tmp_colored_genbank.gb"

        color_genbank.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--target", "domains_and_cdss"])
        assert Path(out).is_file()

        saw_domain_color = False
        saw_cds_color = False
        for rec in parse_seqfiles([out]):
            for feature in rec.features:
                if feature.type == DOMAIN_FEATURE_NAME and 'Color' in feature.qualifiers:
                    saw_domain_color = True
                if feature.type == "CDS" and 'Color' in feature.qualifiers:
                    saw_cds_color = True

        assert saw_domain_color
        assert saw_cds_color

def test_color_genbank_domain_database_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/tmp_colored_genbank.gb"
        color_table_out = output_dir + "/color_table.tsv"

        color_genbank.main([
            "-i", str(shared_datadir / "pDONR201_multi_genemark_domainator_multi_hmm_2.gb"),
            "-o", out,
            "--target", "domains",
            "--domain_database",
            "--color_table_out", color_table_out,
        ])
        assert Path(out).is_file()

        database_colors = {}
        for rec in parse_seqfiles([out]):
            for feature in rec.features:
                if feature.type != DOMAIN_FEATURE_NAME:
                    continue
                db = feature.qualifiers['database'][0]
                color = feature.qualifiers['Color'][0]
                if db in database_colors:
                    assert database_colors[db] == color
                else:
                    database_colors[db] = color

        assert set(database_colors) == {"pdonr_hmms_1", "pdonr_hmms_2"}
        assert len(set(database_colors.values())) == 2

        color_table_dict = {}
        with open(color_table_out, "r") as f:
            for line in f:
                annotation, color = line.strip().split("\t")
                color_table_dict[annotation] = color

        assert color_table_dict == database_colors

def _systems_record(shared_datadir):
    records = find_systems.find_systems(parse_seqfiles([shared_datadir / "systems" / "contigs.gb"]),
                           system_models=[shared_datadir / "systems" / "defense.yaml"])
    return next(r for r in records if r.id == "gabija_ta")


def _colors(record, feature_type):
    return {f.qualifiers.get("cds_id", f.qualifiers.get("name"))[0]: f.qualifiers.get("Color", [None])[0]
            for f in record.features if f.type == feature_type}


def test_color_genbank_systems(shared_datadir):
    record = _systems_record(shared_datadir)
    table = {"Gabija": "#FF0000", "TA": "#0000FF", "defense": "#00FF00", "toxin_antitoxin": "#FFFF00"}
    colored, = color_genbank.color_genbank([record], color_table=table, targets=color_genbank.COLOR_MODES["systems"])
    assert _colors(colored, SYSTEM_FEATURE_NAME) == {"Gabija": "#FF0000", "TA": "#0000FF"}
    assert set(_colors(colored, "CDS").values()) == {None}
    assert set(_colors(colored, DOMAIN_FEATURE_NAME).values()) == {None}

    colored, = color_genbank.color_genbank([_systems_record(shared_datadir)], color_table=table,
                                           targets=color_genbank.COLOR_MODES["systems_and_cdss"], system_color_by="category")
    assert _colors(colored, SYSTEM_FEATURE_NAME) == {"Gabija": "#00FF00", "TA": "#FFFF00"}
    assert _colors(colored, "CDS") == {"cds0": "#00FF00", "cds1": None, "cds2": "#00FF00", "cds3": None, "cds4": None,
                                       "cds5": None, "cds6": "#FFFF00", "cds7": "#FFFF00"}


def test_color_genbank_systems_cli(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        systems = os.path.join(output_dir, "systems.gb")
        out = os.path.join(output_dir, "colored.gb")
        table_out = os.path.join(output_dir, "colors.tsv")
        find_systems.main(["-i", str(shared_datadir / "systems" / "contigs.gb"), "--system_models",
                           str(shared_datadir / "systems" / "defense.yaml"), "-o", systems])
        color_genbank.main(["-i", systems, "-o", out, "--target", "systems", "--color_table_out", table_out])
        with open(table_out) as handle:
            assert sorted(line.split("\t")[0] for line in handle) == ["Gabija", "TA"]  # domains get no colors
