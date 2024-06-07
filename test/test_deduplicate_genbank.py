from domainator import utils
from pathlib import Path
import tempfile
from domainator import deduplicate_genbank

def test_deduplicate_genbank_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/deduplicate_out.gb"
        
        deduplicate_genbank.main(["-i", str(shared_datadir / "extract_peptides_test_2.gb"), "--id", "0.99", "-o", out, "--params", "\"-s\":0.9", "--prefix_count"])
        assert Path(out).is_file()
        recs  = list(utils.parse_seqfiles([out]))
        for rec in recs:
            assert rec.id[:2] == "4-"
        assert len(recs) == 6


def test_deduplicate_genbank_cdhit_low_id_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/deduplicate_out.gb"
        
        deduplicate_genbank.main(["-i", str(shared_datadir / "FeSOD_20.gb"), "--id", "0.4", "--cpu", "1", "-o", out, "--prefix_count"])
        assert Path(out).is_file()
        recs  = list(utils.parse_seqfiles([out]))
        assert recs[0].id[:3] == "10-"
        assert len(recs) == 4

def test_deduplicate_genbank_cdhit_low_id_2(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/deduplicate_out.gb"
        
        deduplicate_genbank.main(["-i", str(shared_datadir / "FeSOD_20.gb"), "--algorithm", "usearch", "--id", "0.4", "--cpu", "1", "-o", out, "--prefix_count"])
        assert Path(out).is_file()
        recs  = list(utils.parse_seqfiles([out]))
        assert recs[0].id[:2] == "6-"
        assert len(recs) == 5

def test_deduplicate_genbank_usearch(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/deduplicate_out.gb"
        log = output_dir + "/deduplicate_out.log"

        deduplicate_genbank.main( ["-i", str(shared_datadir / "extract_peptides_test_2.gb"), "--algorithm", "usearch", "--id", "0.3", "-o", out, "--cpu", "1", "--suffix_count", "--log", log] )
        assert Path(out).is_file()
        recs  = list(utils.parse_seqfiles([out]))
        for rec in recs:
            assert rec.id[-2:] == "-4"
        assert len(recs) == 6

        log_file = open(log, "r").read()
        assert "Input  size: 24" in log_file
        assert "Output size: 6" in log_file

def test_deduplicate_genbank_hash(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/deduplicate_out.gb"
        log = output_dir + "/deduplicate_out.log"
        
        deduplicate_genbank.main(["-i", str(shared_datadir / "extract_peptides_test_2.gb"), "--id", "1", "--algorithm", "hash", "-o", out, "--prefix_count", "--log", log])
        assert Path(out).is_file()
        recs  = list(utils.parse_seqfiles([out]))
        for rec in recs:
            assert rec.id[:2] == "4-"
        assert len(recs) == 6
        
        log_file = open(log, "r").read()
        assert "Input  size: 24" in log_file
        assert "Output size: 6" in log_file

def test_deduplicate_genbank_cdhit_low_id_nt(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/deduplicate_out.gb"
        
        deduplicate_genbank.main(["-i", str(shared_datadir / "pDONR201_multi_genemark.gb"), "--algorithm", "cd-hit", "--id", "0.8", "--cpu", "1", "-o", out, "--prefix_count"])
        assert Path(out).is_file()
        recs  = list(utils.parse_seqfiles([out]))
        assert recs[0].id[:2] == "4-"
        assert len(recs) == 1

def test_deduplicate_genbank_cdhit_fasta_input(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/deduplicate_out.gb"
        log = output_dir + "/deduplicate_out.log"
        
        deduplicate_genbank.main(["-i", str(shared_datadir / "FeSOD_20.fasta"), "--algorithm", "cd-hit", "--id", "0.8", "--cpu", "1", "-o", out, "--prefix_count", "--log", log])
        assert Path(out).is_file()
        recs  = list(utils.parse_seqfiles([out]))
        assert recs[0].id[:2] == "1-"
        assert len(recs) == 20

        log_file = open(log, "r").read()
        assert "Input  size: 20" in log_file
        assert "Output size: 20" in log_file

def test_deduplicate_genbank_cdhit_fasta_input_2(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/deduplicate_out.fasta"
        log = output_dir + "/deduplicate_out.log"
        
        deduplicate_genbank.main(["-i", str(shared_datadir / "pDONR201_multi.fasta"), "--algorithm", "cd-hit", "--id", "0.8", "--cpu", "1", "-o", out, "--prefix_count", "--fasta_type", "nucleotide", "--fasta_out", "--log", log])
        assert Path(out).is_file()
        recs  = list(utils.parse_seqfiles([out]))
        assert recs[0].id[:2] == "4-"
        assert len(recs) == 1

        log_file = open(log, "r").read()
        assert "Input  size: 4" in log_file
        assert "Output size: 1" in log_file
