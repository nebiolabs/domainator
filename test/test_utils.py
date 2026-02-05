from domainator import utils
import tempfile
from domainator.Bio import SeqIO, SeqRecord, Seq
import io
import pytest
from array import array
import re

def test_split_string_list():
    data=["abcde   asdfasdf   asdf", " abcdef::GACAF", ""]
    out = utils.split_string_list(data)
    assert out[0] == ["abcde   asdfasdf   asdf"]
    assert out[1] == ["abcdef", "GACAF"]
    assert out[2] == [""]


def test_write_genbank_1(shared_datadir):
    rec = SeqRecord.SeqRecord(Seq.Seq("GACT"),id="BIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAME",
    name="BIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAME",
    description="BIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAME",
    )
    buf = io.StringIO()
    utils.write_genbank([rec], buf)

def test_write_genbank_space_name(shared_datadir):
    rec = SeqRecord.SeqRecord(Seq.Seq("GACT"),id="BIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAME with_space",
    name="BIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAME with_space",
    description="BIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAMEBIGNAME with_space",
    )
    buf = io.StringIO()
    utils.write_genbank([rec], buf)


def test_list_and_file_to_dict_keys(shared_datadir):
    keys = utils.list_and_file_to_dict_keys(None, str(shared_datadir / "CcdB.hmm"))
    print(keys)
    assert 'CcdB' in keys


# regions are tuples of start and stop coordinates
# returns true if a fraction of region2 >= min_overlap_fraction overlaps with region1
# coordinates within regions must be sorted low to high
@pytest.mark.parametrize("region1,region2,min_overlap_fraction,expected",
[((327,503),(325,507),0.6,True),
((325,507),(327,503),0.6,True),
])
def test_regions_overlap(region1, region2, min_overlap_fraction, expected):

    assert utils.regions_overlap(region1, region2, min_overlap_fraction) == expected

@pytest.mark.parametrize("files,offset,read_count,expected_record_ct,rec0_name",
[(["pDONR201.gb"],0,1,1,"pDONR201"),
(["pDONR201.gb"],0,10,1,"pDONR201"),
(["pDONR201.gb"],0,0,0,""),
(["pDONR201_multi_genemark.gb","pDONR201.gb"],16382,10,2, "pDONR201_3"), #seeks past the end of pDONR201.gb
# (["simple_genpept_equals_second_line.gb"],0,float("inf"),5,"pDONR201_1") #uncomment to test multiline qualifier name handling
])
def test_parse_seqfiles(files,offset,read_count,expected_record_ct,rec0_name,shared_datadir):
    files = [str(shared_datadir / x) for x in files]
    recs = list( utils.parse_seqfiles(files,None,None,offset,read_count) )
    assert len(recs) == expected_record_ct
    if len(recs) > 0:
        assert recs[0].id == rec0_name

@pytest.mark.parametrize("file,offsets,num_proteins",
[("pDONR201.gb",[0],[3]),
("pDONR201_multigenemark_partition.gb",[0,8191,16382,24573],[6,6,6,6]),
("pdonr_peptides.fasta",[0,50,169,226,465],[1,1,1,1,1]),
("pDONR201_empty.gb",[0],[0]),
("simple_genpept.gb",[0,296,987,1620,3904],[1,1,1,1,1]),
])
def test_get_offsets(file,offsets,num_proteins,shared_datadir):
    new_offsets, new_num_proteins = utils.get_offsets(str(shared_datadir / file))
    assert len(new_offsets) == len(new_num_proteins)
    assert new_offsets == array('Q', offsets)
    assert new_num_proteins == array('Q',num_proteins)


def test_get_palette_1():
    palette = utils.get_palette(["A","B","C"])
    assert set(palette.keys()) == {"A","B","C"}
    assert len(palette) == 3
    assert len(set(palette.values())) == 3
    assert all([re.match(r"#[0-9a-fA-F]{6}",x) for x in palette.values()])


class TestBooleanEvaluatorSanitizeIdentifier:
    """Tests for BooleanEvaluator.sanitize_identifier"""
    
    def test_special_characters_replaced(self):
        """Test that BooleanEvaluator special characters are replaced."""
        sanitize = utils.BooleanEvaluator.sanitize_identifier
        
        # Parentheses (grouping)
        assert "(" not in sanitize("x(2)")
        assert ")" not in sanitize("x(2)")
        
        # Operators
        assert "&" not in sanitize("A&B")
        assert "|" not in sanitize("A|B")
        assert "~" not in sanitize("~A")
        
        # Space
        assert " " not in sanitize("A B C")
    
    def test_prosite_patterns(self):
        """Test sanitization of PROSITE-style patterns."""
        sanitize = utils.BooleanEvaluator.sanitize_identifier
        
        # Dashes should be replaced
        assert "-" not in sanitize("A-G-C")
        
        # Anchors
        assert "<" not in sanitize("<M")
        assert ">" not in sanitize("K>")
        
        # Trailing period removed
        result = sanitize("A.")
        assert not result.endswith(".")
    
    def test_complex_pattern(self):
        """Test sanitization of complex patterns."""
        sanitize = utils.BooleanEvaluator.sanitize_identifier
        
        pattern = "[DE](2)HS{P}"
        result = sanitize(pattern)
        
        # Should not contain boolean operator characters
        assert "(" not in result
        assert ")" not in result
        
        # Brackets are allowed (not special in BooleanEvaluator)
        assert "[DE]" in result
    
    def test_round_trip_uniqueness(self):
        """Test that different patterns produce different sanitized names."""
        sanitize = utils.BooleanEvaluator.sanitize_identifier
        
        patterns = ["A(2)", "A(3)", "[DE]", "[EF]", "A-B", "A-C"]
        sanitized = [sanitize(p) for p in patterns]
        
        # All sanitized names should be unique
        assert len(set(sanitized)) == len(patterns)
