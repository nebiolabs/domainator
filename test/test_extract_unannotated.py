import tempfile
import importlib
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from domainator import extract_unannotated
from domainator.seqio import parse_seqfiles, write_genbank
from domainator.model import DomainatorCDS
import pytest
from io import StringIO
import sys
import subprocess


@pytest.mark.parametrize("tool,old,new,value", [
    ("extract_domains", "keep_name", "keep_parent_id", True),
    ("extract_peptides", "keep_name", "keep_parent_id", True),
    ("extract_unannotated", "keep_name", "keep_parent_id", True),
    ("extract_unannotated", "min_length", "min_region_length", 10),
    ("extract_unannotated", "min_score", "min_annotation_score", 50),
])
def test_extraction_parameter_renames(tool, old, new, value, tmp_path, capsys):
    module = importlib.import_module(f"domainator.{tool}")
    base = ["-i", "input.gb", "-o", str(tmp_path / "output.gb")]
    option = [f"--{new}"] if value is True else [f"--{new}", str(value)]
    assert getattr(module.build_parser().parse_args(base + option), new) == value
    assert list(getattr(module, tool)([], **{new: value})) == []

    config = tmp_path / "config.json"
    config.write_text(json.dumps({new: value}))
    params = module.build_parser().parse_args(base + ["--config", str(config)])
    assert getattr(params, new) == value

    old_option = [f"--{old}"] if value is True else [f"--{old}", str(value)]
    with pytest.raises(SystemExit) as error:
        module.build_parser().parse_args(base + old_option)
    assert error.value.code == 2
    assert old in capsys.readouterr().err

    config.write_text(json.dumps({old: value}))
    with pytest.raises(SystemExit) as error:
        module.build_parser().parse_args(base + ["--config", str(config)])
    assert error.value.code == 2
    assert old in capsys.readouterr().err
    with pytest.raises(TypeError, match=old):
        list(getattr(module, tool)([], **{old: value}))


@pytest.mark.parametrize("score,length,expected", [
    (50, 10, ["M" * 10, "G" * 10]),
    (51, 11, ["M" * 10 + "A" * 10 + "G" * 10]),
    (50, 11, []),
])
def test_unannotated_renamed_cutoffs_cli_config_library(score, length, expected, tmp_path):
    record = SeqRecord(Seq("M" * 10 + "A" * 10 + "G" * 10), id="parent",
                       annotations={"molecule_type": "protein"})
    record.features = [SeqFeature(FeatureLocation(10, 20), type="Domainator",
                                  qualifiers={"score": ["50"], "evalue": ["0.001"],
                                              "cds_id": ["0_1_30"], "name": ["mask"],
                                              "database": ["test"], "description": ["mask"]})]
    options = dict(min_annotation_score=score, min_region_length=length, keep_parent_id=True)
    library = list(extract_unannotated.extract_unannotated(
        [record], min_annotation_score=score, min_region_length=length, keep_parent_id=True))
    expected_records = [("parent", sequence) for sequence in expected]
    assert [(r.id, str(r.seq)) for r in library] == expected_records

    source = tmp_path / "input.gb"
    with source.open("w") as handle:
        write_genbank([record], handle)
    config = tmp_path / "config.json"
    config.write_text(json.dumps(options))
    for mode in ("cli", "config"):
        output = tmp_path / f"{mode}.gb"
        args = ["-i", str(source), "-o", str(output)]
        if mode == "cli":
            args += ["--min_annotation_score", str(score), "--min_region_length",
                     str(length), "--keep_parent_id"]
        else:
            args += ["--config", str(config)]
        extract_unannotated.main(args)
        assert [(r.id, str(r.seq)) for r in parse_seqfiles([output])] == expected_records


def test_extract_unannotated_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        extract_unannotated.main(["-i", str(shared_datadir / "simple_genpept.gb"), "-o", out])
        # assert 0
        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 8
        [str(seq) for seq in seqs] == ["FPALSPDSVDNRITASQEEFVETQKGHPSGWPSA*","M","I*","*","M","DEWQGGA*","MDADLYGYKWARDNVGQSGATIYRLYGKPDAPELFLKHGKGSVANDV","DMNKLQFHLMLDEFF*"]

def test_extract_unannotated_largest_keep_parent_id_2(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        extract_unannotated.main(["-i", str(shared_datadir / "simple_genpept.gb"), "-o", out, "--largest", "--keep_parent_id"])
        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 5
        [str(seq) for seq in seqs] == ["FPALSPDSVDNRITASQEEFVETQKGHPSGWPSA*","I*","*","DEWQGGA*","MDADLYGYKWARDNVGQSGATIYRLYGKPDAPELFLKHGKGSVANDV"]

def test_extract_unannotated_min_region_length_2(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        extract_unannotated.main(["-i", str(shared_datadir / "simple_genpept.gb"), "-o", out, "--min_region_length", "10"])
        seqs = list(parse_seqfiles([out]))
        assert len(seqs) == 3
        [str(seq) for seq in seqs] == ["FPALSPDSVDNRITASQEEFVETQKGHPSGWPSA*", "MDADLYGYKWARDNVGQSGATIYRLYGKPDAPELFLKHGKGSVANDV", "DMNKLQFHLMLDEFF*"]