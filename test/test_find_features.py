"""Tests for find_features.py"""

import os
import pytest
import tempfile
from domainator.Bio import SeqIO
from domainator import DOMAIN_FEATURE_NAME


# Check if tmbed is available
try:
    import tmbed
    TMBED_AVAILABLE = True
except ImportError:
    TMBED_AVAILABLE = False


def test_find_features_import():
    """Test that find_features module can be imported."""
    from domainator import find_features
    assert hasattr(find_features, 'main')
    assert hasattr(find_features, '_entrypoint')
    assert hasattr(find_features, 'find_features')


def test_score_to_pseudo_evalue():
    """Test pseudo-evalue calculation."""
    from domainator.find_features import score_to_pseudo_evalue
    import math
    
    # Perfect score should give very low evalue
    assert score_to_pseudo_evalue(100.0) < 1e-40
    
    # Zero score should give evalue of 1
    assert abs(score_to_pseudo_evalue(0.0) - 1.0) < 1e-10
    
    # Higher scores should give lower evalues
    assert score_to_pseudo_evalue(90.0) < score_to_pseudo_evalue(50.0)
    assert score_to_pseudo_evalue(50.0) < score_to_pseudo_evalue(10.0)
    
    # Check specific values: e^(-90) ≈ 1.2e-39
    assert 1e-40 < score_to_pseudo_evalue(90.0) < 1e-38


def test_feature_result_namedtuple():
    """Test FeatureResult namedtuple."""
    from domainator.find_features import FeatureResult
    
    result = FeatureResult(
        name="transmembrane_alpha_helix",
        desc="transmembrane alpha helix",
        score=95.0,
        evalue=1e-20,
        start=10,
        end=30,
        program="TMbed",
        identity=95.0,
        probabilities={"P(H)": [0.95, 0.96, 0.94]}
    )
    
    assert result.name == "transmembrane_alpha_helix"
    assert result.score == 95.0
    assert result.start == 10
    assert result.end == 30
    assert result.program == "TMbed"


def test_tmbed_labels():
    """Test TMbed label mapping."""
    from domainator.find_features import TMBED_LABELS
    
    # Check that all expected labels are present
    assert 0 in TMBED_LABELS  # B
    assert 1 in TMBED_LABELS  # b -> B
    assert 2 in TMBED_LABELS  # H
    assert 3 in TMBED_LABELS  # h -> H
    assert 4 in TMBED_LABELS  # S
    assert 5 in TMBED_LABELS  # i
    assert 6 in TMBED_LABELS  # o
    
    # Check label structure
    for label_id, (char, name, desc) in TMBED_LABELS.items():
        assert isinstance(char, str)
        assert isinstance(name, str)
        assert isinstance(desc, str)


@pytest.mark.skipif(not TMBED_AVAILABLE, reason="TMbed not installed")
def test_find_features_protein_fasta(shared_datadir):
    """Test find_features with a protein fasta file."""
    from domainator.find_features import main
    
    # Use a simple protein fasta
    fasta = shared_datadir / "FeSOD_20.fasta"
    
    with tempfile.TemporaryDirectory() as output_dir:
        out = os.path.join(output_dir, "out.gb")
        args = [
            '--input', str(fasta),
            '-o', str(out),
            '--algorithms', 'TMbed',
            '--gpu_device', 'None',  # Use CPU for testing
        ]
        
        main(args)
        
        assert os.path.isfile(out)
        
        # Parse output and check for annotations
        records = list(SeqIO.parse(out, "genbank"))
        assert len(records) > 0


@pytest.mark.skipif(not TMBED_AVAILABLE, reason="TMbed not installed")
def test_find_features_genbank(shared_datadir):
    """Test find_features with a genbank file."""
    from domainator.find_features import main
    
    # Use a genbank file with CDS annotations
    gb = shared_datadir / "FeSOD_20.gb"
    
    with tempfile.TemporaryDirectory() as output_dir:
        out = os.path.join(output_dir, "out.gb")
        args = [
            '--input', str(gb),
            '-o', str(out),
            '--algorithms', 'TMbed',
            '--gpu_device', 'None',  # Use CPU for testing
        ]
        
        main(args)
        
        assert os.path.isfile(out)
        
        # Parse output
        records = list(SeqIO.parse(out, "genbank"))
        assert len(records) > 0


@pytest.mark.skipif(not TMBED_AVAILABLE, reason="TMbed not installed")
def test_find_features_algorithms_validation():
    """Test that invalid algorithms raise errors."""
    from domainator.find_features import main
    import tempfile
    
    with tempfile.NamedTemporaryFile(suffix='.fasta', mode='w') as f:
        f.write(">test\nMKFLILLFNILCLFPVLAADNHGVGPQGASGVWD\n")
        f.flush()
        
        with tempfile.TemporaryDirectory() as output_dir:
            out = os.path.join(output_dir, "out.gb")
            args = [
                '--input', f.name,
                '-o', str(out),
                '--algorithms', 'InvalidAlgorithm',
            ]
            
            with pytest.raises(ValueError, match="Unknown algorithm"):
                main(args)


# ============================================================================
# PROSITE Motif Search Tests
# ============================================================================

def test_prosite_to_regex():
    """Test PROSITE pattern to regex conversion."""
    from domainator.find_features import prosite_to_regex
    
    # Basic amino acid
    assert prosite_to_regex("A") == "A"
    assert prosite_to_regex("ALA") == "ALA"
    
    # Any amino acid (x)
    assert prosite_to_regex("x") == "."
    assert prosite_to_regex("X") == "."
    assert prosite_to_regex("A-x-G") == "A.G"
    
    # Character class [ABC]
    assert prosite_to_regex("[ALT]") == "[ALT]"
    assert prosite_to_regex("A-[LM]-G") == "A[LM]G"
    
    # Negated class {ABC}
    assert prosite_to_regex("{AM}") == "[^AM]"
    assert prosite_to_regex("A-{P}-G") == "A[^P]G"
    
    # Repetition (n)
    assert prosite_to_regex("x(3)") == ".{3}"
    assert prosite_to_regex("A(2)") == "A{2}"
    
    # Range repetition (n,m)
    assert prosite_to_regex("x(2,4)") == ".{2,4}"
    
    # Anchors
    assert prosite_to_regex("<M") == "^M"
    assert prosite_to_regex("K>") == "K$"
    
    # Trailing period (should be ignored)
    assert prosite_to_regex("A-G.") == "AG"
    
    # Complex pattern from PROSITE docs
    pattern = "[DE](2)HS{P}x(2)Px(2,4)C"
    regex = prosite_to_regex(pattern)
    assert regex == "[DE]{2}HS[^P].{2}P.{2,4}C"


def test_sanitize_motif_name():
    """Test motif name sanitization for boolean evaluator compatibility."""
    from domainator.find_features import sanitize_motif_name
    
    # Parentheses should be replaced
    assert "(" not in sanitize_motif_name("x(2)")
    assert ")" not in sanitize_motif_name("x(2)")
    
    # Dashes should be replaced
    assert "-" not in sanitize_motif_name("A-G-C")
    
    # Trailing period should be removed
    assert sanitize_motif_name("A.") == "A"
    
    # Check specific transformations
    name = sanitize_motif_name("[DE](2)HS{P}")
    assert "(" not in name
    assert ")" not in name
    # Brackets are allowed (not special in boolean evaluator)
    assert "[DE]" in name


def test_search_motif_simple():
    """Test simple motif search."""
    from domainator.find_features import search_motif
    
    # Simple exact match
    results = search_motif("ACDEFGHIKLMNPQRSTVWY", "DEF")
    assert len(results) == 1
    assert results[0].start == 2
    assert results[0].end == 5
    assert results[0].score == 100.0
    assert results[0].evalue == 0.0
    assert results[0].program == "MOTIF"
    
    # No match
    results = search_motif("AAAAAAA", "DEF")
    assert len(results) == 0


def test_search_motif_character_class():
    """Test motif search with character classes."""
    from domainator.find_features import search_motif
    
    # [DE] matches D or E
    results = search_motif("AADAAEAA", "[DE]")
    assert len(results) == 2
    
    # {P} matches anything except P
    results = search_motif("APG", "A{P}G")
    assert len(results) == 0  # P is excluded
    
    results = search_motif("AXG", "A{P}G")
    assert len(results) == 1  # X is allowed


def test_search_motif_repetition():
    """Test motif search with repetition."""
    from domainator.find_features import search_motif
    
    # x(3) matches any 3 residues
    # In "AABCDEF", "Ax(3)E" matches at position 1: A-BCD-E
    results = search_motif("AABCDEF", "Ax(3)E")
    assert len(results) == 1
    assert results[0].start == 1  # Matches "ABCDE" starting at index 1
    assert results[0].end == 6
    
    # (2,4) range
    results = search_motif("AXXB", "Ax(2,4)B")
    assert len(results) == 1


def test_search_motif_case_insensitive():
    """Test that motif search is case insensitive."""
    from domainator.find_features import search_motif
    
    results_upper = search_motif("ACDEFG", "CDE")
    results_lower = search_motif("acdefg", "CDE")
    results_mixed = search_motif("AcDeFg", "cde")
    
    assert len(results_upper) == 1
    assert len(results_lower) == 1
    assert len(results_mixed) == 1


def test_search_all_motifs():
    """Test searching for multiple motifs."""
    from domainator.find_features import search_all_motifs
    
    sequence = "ACDEFGHIKLMNPQRSTVWY"
    motifs = ["CDE", "MNP"]
    
    results = search_all_motifs(sequence, motifs)
    assert len(results) == 2
    
    names = [r.name for r in results]
    assert "CDE" in names
    assert "MNP" in names


def test_motif_integration():
    """Test motif search integration with find_features main."""
    from domainator.find_features import main
    
    with tempfile.NamedTemporaryFile(suffix='.fasta', mode='w', delete=False) as f:
        # Create a test sequence with a known motif
        # Pattern [DE](2)HS matches "DDHS" or "DEHS" or "EDHS" or "EEHS"
        f.write(">test\nMKFLDDHSAXXPXXXXCILLFNILCLFPVLA\n")
        f.flush()
        
        try:
            with tempfile.TemporaryDirectory() as output_dir:
                out = os.path.join(output_dir, "out.gb")
                args = [
                    '--input', f.name,
                    '-o', str(out),
                    '--algorithms', 'none',  # Don't run ML algorithms
                    '--motif', '[DE](2)HS',
                ]
                
                main(args)
                
                # Verify output was created
                assert os.path.isfile(out)
                
                # Parse and check for motif annotations
                records = list(SeqIO.parse(out, "genbank"))
                assert len(records) == 1
                
                # Find Domainator features
                motif_features = [f for f in records[0].features 
                                  if f.type == DOMAIN_FEATURE_NAME 
                                  and f.qualifiers.get('program', [''])[0] == 'MOTIF']
                
                assert len(motif_features) == 1
                feat = motif_features[0]
                assert feat.qualifiers['database'][0] == 'MOTIF'
                assert feat.qualifiers['evalue'][0] == '0.0e+00'
                assert feat.qualifiers['score'][0] == '100.0'
                assert feat.qualifiers['identity'][0] == '100.0'
        finally:
            os.unlink(f.name)


def test_motif_annotation_format():
    """Test that motif annotations have correct format."""
    from domainator.find_features import search_motif, FeatureResult
    
    pattern = "[FY]-[LIV]-G-[DE]-E-A-Q-x-[RKQ](2)-G"
    results = search_motif("MFLGDEAQARRGXXX", pattern)
    
    if len(results) > 0:
        result = results[0]
        assert result.program == "MOTIF"
        assert result.desc == pattern  # Description is the raw pattern
        assert result.score == 100.0
        assert result.evalue == 0.0
        assert result.identity == 100.0
