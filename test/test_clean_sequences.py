import tempfile

from domainator import clean_sequences
from domainator import utils
from domainator.Bio.SeqRecord import SeqRecord
from domainator.Bio.Seq import Seq


def test_clean_name_for_newick():
    """Test that problematic characters are replaced with underscores"""
    assert clean_sequences.clean_name_for_newick("simple") == "simple"
    assert clean_sequences.clean_name_for_newick("with space") == "with_space"
    assert clean_sequences.clean_name_for_newick("has;semicolon") == "has_semicolon"
    assert clean_sequences.clean_name_for_newick("has:colon") == "has_colon"
    assert clean_sequences.clean_name_for_newick("has,comma") == "has_comma"
    assert clean_sequences.clean_name_for_newick("has(parens)") == "has_parens_"
    assert clean_sequences.clean_name_for_newick("has'quote") == "has_quote"
    assert clean_sequences.clean_name_for_newick('has"double') == "has_double"


def test_deduplicate_name():
    """Test that duplicate names get unique suffixes"""
    seen = {}
    assert clean_sequences.deduplicate_name("seq1", seen) == "seq1"
    assert clean_sequences.deduplicate_name("seq2", seen) == "seq2"
    assert clean_sequences.deduplicate_name("seq1", seen) == "seq1_1"
    assert clean_sequences.deduplicate_name("seq1", seen) == "seq1_2"
    assert clean_sequences.deduplicate_name("seq1_1", seen) == "seq1_1_1"


def test_strip_non_canonical_protein():
    """Test stripping non-canonical amino acids from ends"""
    canonical = clean_sequences.CANONICAL_PROTEIN
    
    # Simple case with * at end
    assert clean_sequences.strip_non_canonical("MKTLV*", canonical) == "MKTLV"
    
    # * at both ends
    assert clean_sequences.strip_non_canonical("*MKTLV*", canonical) == "MKTLV"
    
    # Non-canonical at start (X is non-canonical in standard set)
    assert clean_sequences.strip_non_canonical("XMKTLV", canonical) == "MKTLV"
    
    # Multiple non-canonical at ends
    assert clean_sequences.strip_non_canonical("XXMKTLVXX", canonical) == "MKTLV"
    
    # No stripping needed
    assert clean_sequences.strip_non_canonical("MKTLV", canonical) == "MKTLV"
    
    # All non-canonical
    assert clean_sequences.strip_non_canonical("***", canonical) == ""


def test_strip_non_canonical_nucleotide():
    """Test stripping non-canonical nucleotides from ends"""
    canonical = clean_sequences.CANONICAL_NUCLEOTIDE
    
    # N is ambiguous, not canonical
    assert clean_sequences.strip_non_canonical("NATGCN", canonical) == "ATGC"
    
    # Simple case
    assert clean_sequences.strip_non_canonical("ATGC", canonical) == "ATGC"


def test_has_non_canonical_protein():
    """Test detection of non-canonical amino acids"""
    canonical = clean_sequences.CANONICAL_PROTEIN
    
    assert clean_sequences.has_non_canonical("MKTLV", canonical) == False
    assert clean_sequences.has_non_canonical("MKTLV*", canonical) == True
    assert clean_sequences.has_non_canonical("MXKTLV", canonical) == True
    assert clean_sequences.has_non_canonical("MBKTLV", canonical) == True  # B is ambiguous


def test_has_non_canonical_nucleotide():
    """Test detection of non-canonical nucleotides"""
    canonical = clean_sequences.CANONICAL_NUCLEOTIDE
    
    assert clean_sequences.has_non_canonical("ATGC", canonical) == False
    assert clean_sequences.has_non_canonical("ATGCN", canonical) == True  # N is ambiguous


def test_clean_sequences_clean_name(shared_datadir):
    """Test cleaning sequence names"""
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/clean_out.gb"
        clean_sequences.main(["-i", str(shared_datadir / "simple_genpept.gb"), "-o", out, "--clean_name"])
        seqs = list(utils.parse_seqfiles([out]))
        # Check that no names have problematic characters
        for seq in seqs:
            assert " " not in seq.id
            assert ";" not in seq.id
            assert ":" not in seq.id


def test_clean_sequences_deduplicate(shared_datadir):
    """Test deduplicating sequence names"""
    # Create input with duplicate names
    records = [
        SeqRecord(Seq("MKTLV"), id="seq1", name="seq1", description="test"),
        SeqRecord(Seq("AKSLV"), id="seq1", name="seq1", description="test2"),
        SeqRecord(Seq("GKPLV"), id="seq1", name="seq1", description="test3"),
    ]
    
    result = list(clean_sequences.clean_sequences(records, deduplicate_names=True))
    
    ids = [r.id for r in result]
    assert len(ids) == len(set(ids))  # All unique
    assert "seq1" in ids
    assert "seq1_1" in ids
    assert "seq1_2" in ids


def test_clean_sequences_strip():
    """Test stripping non-canonical residues"""
    records = [
        SeqRecord(Seq("MKTLV*"), id="seq1", name="seq1", description="test"),
        SeqRecord(Seq("*AKSLV*"), id="seq2", name="seq2", description="test2"),
    ]
    
    result = list(clean_sequences.clean_sequences(records, strip=True))
    
    assert str(result[0].seq) == "MKTLV"
    assert str(result[1].seq) == "AKSLV"


def test_clean_sequences_filter():
    """Test filtering sequences with non-canonical residues"""
    records = [
        SeqRecord(Seq("MKTLV"), id="seq1", name="seq1", description="test"),
        SeqRecord(Seq("MKTLVX"), id="seq2", name="seq2", description="test2"),
        SeqRecord(Seq("AKSLV"), id="seq3", name="seq3", description="test3"),
    ]
    
    result = list(clean_sequences.clean_sequences(records, filter_by_aa=True))
    
    assert len(result) == 2
    ids = [r.id for r in result]
    assert "seq1" in ids
    assert "seq3" in ids
    assert "seq2" not in ids


def test_clean_sequences_strip_then_filter():
    """Test that strip is applied before filter"""
    records = [
        SeqRecord(Seq("MKTLV*"), id="seq1", name="seq1", description="test"),  # Should pass after strip
        SeqRecord(Seq("MXTLV"), id="seq2", name="seq2", description="test2"),  # Should fail filter (X in middle)
    ]
    
    result = list(clean_sequences.clean_sequences(records, strip=True, filter_by_aa=True))
    
    assert len(result) == 1
    assert result[0].id == "seq1"
    assert str(result[0].seq) == "MKTLV"


def test_clean_sequences_empty_after_strip():
    """Test that empty sequences after strip are filtered out"""
    records = [
        SeqRecord(Seq("***"), id="seq1", name="seq1", description="test"),
        SeqRecord(Seq("MKTLV"), id="seq2", name="seq2", description="test2"),
    ]
    
    result = list(clean_sequences.clean_sequences(records, strip=True))
    
    assert len(result) == 1
    assert result[0].id == "seq2"


def test_clean_sequences_nucleotide():
    """Test with nucleotide sequences"""
    records = [
        SeqRecord(Seq("NATGCN"), id="seq1", name="seq1", description="test"),
        SeqRecord(Seq("ATGC"), id="seq2", name="seq2", description="test2"),
    ]
    
    result = list(clean_sequences.clean_sequences(records, fasta_type="nucleotide", strip=True))
    
    assert str(result[0].seq) == "ATGC"
    assert str(result[1].seq) == "ATGC"


def test_clean_sequences_fasta_output(shared_datadir):
    """Test fasta output"""
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/clean_out.fasta"
        clean_sequences.main(["-i", str(shared_datadir / "simple_genpept.gb"), "-o", out, "--fasta_out"])
        # Just verify the file exists and can be read
        seqs = list(utils.parse_seqfiles([out]))
        assert len(seqs) > 0


def test_clean_sequences_combined_options():
    """Test multiple options combined"""
    records = [
        SeqRecord(Seq("MKTLV*"), id="seq (1)", name="seq (1)", description="test"),
        SeqRecord(Seq("AKSLV*"), id="seq (1)", name="seq (1)", description="test2"),
    ]
    
    result = list(clean_sequences.clean_sequences(
        records, 
        clean_name=True, 
        deduplicate_names=True, 
        strip=True
    ))
    
    assert len(result) == 2
    # Names should be cleaned
    assert "(" not in result[0].id
    assert "(" not in result[1].id
    # Names should be deduplicated
    assert result[0].id != result[1].id
    # Sequences should be stripped
    assert str(result[0].seq) == "MKTLV"
    assert str(result[1].seq) == "AKSLV"


def test_clean_sequences_stats():
    """Test that statistics are tracked correctly"""
    records = [
        SeqRecord(Seq("MKTLV*"), id="seq (1)", name="seq (1)", description="test"),
        SeqRecord(Seq("AKSLV*"), id="seq (1)", name="seq (1)", description="test2"),
        SeqRecord(Seq("GXPLV"), id="seq3", name="seq3", description="test3"),  # Will be filtered
        SeqRecord(Seq("***"), id="seq4", name="seq4", description="test4"),  # Will be empty after strip
    ]
    
    stats = clean_sequences.CleaningStats()
    result = list(clean_sequences.clean_sequences(
        records, 
        clean_name=True, 
        deduplicate_names=True, 
        strip=True,
        filter_by_aa=True,
        stats=stats
    ))
    
    assert stats.input_count == 4
    assert stats.output_count == 2
    assert stats.filtered_count == 1  # seq3 filtered for having X
    assert stats.empty_after_strip_count == 1  # seq4 empty after strip
    assert stats.names_changed_count == 2  # Both output seqs had names changed
    assert stats.sequences_changed_count == 2  # Both output seqs had sequences changed (strip)


def test_clean_sequences_log_file(shared_datadir):
    """Test that log file is written correctly"""
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/clean_out.gb"
        log = output_dir + "/clean_log.txt"
        clean_sequences.main(["-i", str(shared_datadir / "simple_genpept.gb"), "-o", out, "--strip", "--log", log])
        
        with open(log, "r") as f:
            log_content = f.read()
        
        assert "Input sequences:" in log_content
        assert "Output sequences:" in log_content
        assert "Names changed:" in log_content
        assert "Sequences changed:" in log_content
