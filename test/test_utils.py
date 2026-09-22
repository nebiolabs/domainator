from domainator import utils, DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME
import tempfile
from domainator.Bio import SeqIO, SeqRecord, Seq
from domainator.Bio.SeqFeature import FeatureLocation, CompoundLocation, SeqFeature
import io
import pytest
from array import array
import re
import math
import pandas as pd
import helpers
from pyhmmer import easel


def test_location_covers():
    assert utils.location_covers([(0, 100)], [(10, 20)])
    assert utils.location_covers([(0, 100)], [(0, 100)])
    assert not utils.location_covers([(0, 50)], [(10, 60)])
    # gappy outer (order with a gap 50..60): a region in a part is covered, one spanning the gap isn't
    assert utils.location_covers([(0, 50), (60, 100)], [(10, 20)])
    assert not utils.location_covers([(0, 50), (60, 100)], [(40, 70)])
    # origin-spanning region expressed as two parts, both inside outer
    assert utils.location_covers([(0, 100)], [(90, 100), (0, 10)])
    # adjacent outer parts merge
    assert utils.location_covers([(0, 50), (50, 100)], [(40, 60)])


def _source_feature(taxid, parts, operator=None):
    if operator:
        loc = CompoundLocation([FeatureLocation(s, e) for s, e in parts], operator=operator)
    else:
        s, e = parts[0]
        loc = FeatureLocation(s, e)
    return SeqFeature(location=loc, type="source", qualifiers={"db_xref": [f"taxon:{taxid}"]})


def test_location_taxid_longest_covering():
    # Gappy host (order) covering everything except a prophage gap (400..600); prophage fills it.
    host = _source_feature(100, [(0, 400), (600, 1000)], operator="order")
    phage = _source_feature(200, [(400, 600)])
    sources = [host, phage]
    assert utils.location_taxid([(100, 200)], sources) == 100   # in a host part
    assert utils.location_taxid([(450, 550)], sources) == 200   # in the prophage gap
    assert utils.location_taxid([(0, 1000)], sources) is None   # no single source covers the whole contig
    # Nested: a full-genome host source covers everything, so it wins as the longest covering.
    host_full = _source_feature(100, [(0, 1000)])
    assert utils.location_taxid([(450, 550)], [host_full, phage]) == 100
    assert utils.location_taxid([(0, 1000)], [host_full, phage]) == 100  # whole-contig source exists
    # Covered by a source with no taxon -> unidentified (32644), not None.
    notax = SeqFeature(location=FeatureLocation(0, 1000), type="source", qualifiers={})
    assert utils.location_taxid([(10, 20)], [notax]) == 32644
    # Uncovered region -> None.
    assert utils.location_taxid([(5000, 5001)], sources) is None


def test_circular_wrapped_contig_hit_taxid():
    # A nucleotide hit spanning the origin of a circular contig maps to a join; resolving it
    # against the whole-contig source must keep it (not drop it as "uncovered"). This guards
    # the nucleotide post-filter against using raw (start, end) for wrapped hits.
    from types import SimpleNamespace
    from domainator.domainate import build_contig_hit_location
    contig = SeqRecord.SeqRecord(Seq.Seq("A" * 100), id="c", description="c")
    contig.annotations = {"molecule_type": "DNA", "topology": "circular"}
    src = SeqFeature(FeatureLocation(0, 100), type="source", qualifiers={"db_xref": ["taxon:42"]})
    hit = SimpleNamespace(start=90, end=110, strand=1)  # wraps past the origin
    loc = build_contig_hit_location(contig, hit)
    region = [(int(p.start), int(p.end)) for p in loc.parts]
    assert region == [(90, 100), (0, 10)]
    assert utils.location_taxid(region, [src]) == 42

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


def test_domainator_cds_creates_one_pseudo_cds_per_nucleic_acid_annotation():
    record = SeqRecord.SeqRecord(Seq.Seq("ATGCGTACGTAA"), id="dna_target")
    record.annotations["molecule_type"] = "DNA"
    shared_location = FeatureLocation(2, 8, strand=1)
    record.features = [
        SeqFeature(FeatureLocation(0, len(record.seq)), type="source", qualifiers={}),
        SeqFeature(
            shared_location,
            type=DOMAIN_FEATURE_NAME,
            qualifiers={
                "name": ["dna_query_1"],
                "description": ["query 1"],
                "database": ["nuc_db"],
                "cds_id": ["."],
                "evalue": ["1e-20"],
                "score": ["50"],
            },
        ),
        SeqFeature(
            shared_location,
            type=DOMAIN_FEATURE_NAME,
            qualifiers={
                "name": ["dna_query_2"],
                "description": ["query 2"],
                "database": ["nuc_db"],
                "cds_id": ["."],
                "evalue": ["1e-10"],
                "score": ["40"],
            },
        ),
        SeqFeature(
            shared_location,
            type=DOMAIN_SEARCH_BEST_HIT_NAME,
            qualifiers={
                "name": ["dna_query_1"],
                "description": ["query 1"],
                "cds_id": ["."],
                "evalue": ["1e-30"],
                "score": ["60"],
                "rstart": ["1"],
                "rend": ["6"],
                "rlen": ["6"],
            },
        ),
    ]

    cdss = utils.DomainatorCDS.list_from_contig(record, include_nucleic_acid_annotations=True)

    assert len(cdss) == 3
    assert all(cds.is_nucleic_acid for cds in cdss)
    assert [cds.num for cds in cdss] == ["nuc_0", "nuc_1", "nuc_2"]
    assert [cds.name for cds in cdss] == ["dna_query_1", "dna_query_2", "dna_query_1"]
    assert [len(cds.domain_features) for cds in cdss] == [1, 1, 0]
    assert cdss[0].domain_features[0].qualifiers["name"] == ["dna_query_1"]
    assert cdss[1].domain_features[0].qualifiers["name"] == ["dna_query_2"]
    assert cdss[2].domain_search_feature is not None
    assert cdss[2].domain_search_feature.qualifiers["name"] == ["dna_query_1"]


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


def test_get_palette_cycles_over_full_palette():
    """More groups than colors cycles, but adjacent numbers never collide."""
    palette = utils.get_palette(list(range(1, 201)))
    assert len(set(palette.values())) == len(utils.DISTINCT_COLORS)
    assert all(palette[i] != palette[i + 1] for i in range(1, 200))
    # the only collisions are a full palette apart
    assert palette[1] == palette[1 + len(utils.DISTINCT_COLORS)]


def test_get_palette_assignment_is_sorted():
    """Colors follow the values, not the order they happen to appear in."""
    assert utils.get_palette([3, 1, 2]) == utils.get_palette([1, 2, 3])
    assert utils.get_palette([3, 1, 2])[1] == utils.DISTINCT_COLORS[0]
    # numeric-looking strings sort numerically, not lexically
    assert utils.get_palette(["10", "2", "1"])["2"] == utils.DISTINCT_COLORS[1]


def test_get_palette_max_groups():
    values = pd.Series(["a"] * 10 + ["b"] * 5 + ["c"] * 2 + ["d"])
    palette = utils.get_palette(values, max_groups=2)
    assert len(palette) == 4
    assert len({palette["a"], palette["b"]}) == 2
    assert palette["c"] == palette["d"] == utils.OTHER_COLOR
    # max_groups at or above the group count changes nothing
    assert utils.get_palette(values, max_groups=4) == utils.get_palette(values)


def test_get_palette_missing_values():
    """NaN gets the neutral color under key None, without consuming a palette slot."""
    palette = utils.get_palette(pd.Series(["A", "B", None]))
    assert palette == {"A": utils.DISTINCT_COLORS[0], "B": utils.DISTINCT_COLORS[1],
                       None: utils.OTHER_COLOR}


def test_sort_palette_values():
    """Mixed-type columns sort without raising: numbers, then strings, then missing."""
    assert utils.sort_palette_values(["10", "2", "zed", "1", None, "apple"]) == \
        ["1", "2", "10", "apple", "zed", None]


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


# --- hmm alphabet helpers -------------------------------------------------

def test_alphabet_score_scale_amino_is_exactly_one():
    # Every committed protein fixture (e.g. pDONR_201_hmm_scores.tsv, which
    # test_hmmer_compare_1 byte-compares) depends on the amino scale being an
    # exact no-op, so assert exactly rather than approximately.
    assert utils.alphabet_score_scale(easel.Alphabet.amino()) == 1.0


def test_alphabet_score_scale_nucleotide():
    expected = math.log(4) / math.log(20)
    assert utils.alphabet_score_scale(easel.Alphabet.dna()) == pytest.approx(expected)
    assert utils.alphabet_score_scale(easel.Alphabet.rna()) == pytest.approx(expected)


def test_alphabet_name():
    assert utils.alphabet_name(easel.Alphabet.amino()) == "amino"
    assert utils.alphabet_name(easel.Alphabet.dna()) == "DNA"
    assert utils.alphabet_name(easel.Alphabet.rna()) == "RNA"


def test_is_nucleic_acid_alphabet():
    assert not utils.is_nucleic_acid_alphabet(easel.Alphabet.amino())
    assert utils.is_nucleic_acid_alphabet(easel.Alphabet.dna())
    assert utils.is_nucleic_acid_alphabet(easel.Alphabet.rna())


def test_get_alphabet():
    assert utils.get_alphabet(None) is None
    assert utils.get_alphabet("DNA") == easel.Alphabet.dna()
    assert utils.get_alphabet("rna") == easel.Alphabet.rna()
    assert utils.get_alphabet("amino") == easel.Alphabet.amino()
    with pytest.raises(ValueError, match="Unknown alphabet"):
        utils.get_alphabet("protein")


def test_pyhmmer_alphabet_is_unhashable():
    # The alphabet helpers must compare alphabets pairwise rather than putting them
    # in a set or dict, or using them as an lru_cache key.
    with pytest.raises(TypeError):
        hash(easel.Alphabet.dna())


def test_peek_hmm_alphabet(shared_datadir):
    assert utils.peek_hmm_alphabet(shared_datadir / "CcdB.hmm").is_amino()
    assert utils.peek_hmm_alphabet(shared_datadir / "dna_profiles.hmm").is_dna()
    assert utils.peek_hmm_alphabet(shared_datadir / "rna_profiles.hmm").is_rna()


def test_peek_hmm_alphabet_rejects_streams(shared_datadir):
    with open(shared_datadir / "CcdB.hmm", "rb") as handle:
        with pytest.raises(TypeError, match="file paths"):
            utils.peek_hmm_alphabet(handle)


def test_common_hmm_alphabet(shared_datadir):
    alphabet = utils.common_hmm_alphabet([shared_datadir / "CcdB.hmm", shared_datadir / "pdonr_hmms.hmm"])
    assert alphabet.is_amino()


def test_common_hmm_alphabet_raises_on_disagreement(shared_datadir):
    with pytest.raises(ValueError, match="Mismatched alphabets among the reference hmm files"):
        utils.common_hmm_alphabet(
            [shared_datadir / "CcdB.hmm", shared_datadir / "dna_profiles.hmm"],
            role="reference",
        )


def test_iter_hmms_with_alphabet_yields_every_profile(shared_datadir):
    # The first profile is read eagerly to learn the alphabet, then chained back on,
    # so guard against it being dropped or duplicated.
    alphabet, profiles = utils.iter_hmms_with_alphabet(shared_datadir / "pdonr_hmms.hmm")
    assert alphabet.is_amino()
    names = [utils.pyhmmer_decode(hmm.name) for hmm in profiles]
    assert len(names) == len(set(names)) == 7


def test_iter_hmms_with_alphabet_empty_file(tmp_path):
    # pyhmmer refuses to open a completely empty hmm file; the helper does not mask that.
    empty = tmp_path / "empty.hmm"
    empty.write_bytes(b"")
    with pytest.raises(EOFError):
        utils.iter_hmms_with_alphabet(empty)


def test_iter_hmms_reports_mixed_alphabet_file(shared_datadir, tmp_path):
    # pyhmmer locks an HMMFile to the alphabet of its first profile and then raises a
    # bare "Expected DNA alphabet"; the wrapped error must name the file.
    mixed = tmp_path / "mixed.hmm"
    mixed.write_bytes(
        (shared_datadir / "dna_profiles_1.hmm").read_bytes()
        + (shared_datadir / "CcdB.hmm").read_bytes()
    )
    with pytest.raises(ValueError, match="more than one alphabet"):
        list(utils.iter_hmms(mixed))

    _alphabet, profiles = utils.iter_hmms_with_alphabet(mixed)
    with pytest.raises(ValueError, match="mixed.hmm"):
        list(profiles)


# ---------------------------------------------------------------------------
# Compressed hmm files, pressed-sidecar handling, and database naming.
# ---------------------------------------------------------------------------

def _names(hmm_file):
    with hmm_file as handle:
        return [utils.pyhmmer_decode(h.name) for h in handle]


def test_open_hmm_file_reads_gzip_and_bgzf(shared_datadir, tmp_path):
    src = shared_datadir / "FeSOD_pfam.hmm"
    plain = _names(utils.open_hmm_file(src))
    assert plain  # sanity
    assert _names(utils.open_hmm_file(helpers.gzip_file(src, tmp_path / "r.hmm.gz"))) == plain
    # BGZF used to raise "format not recognized by HMMER": easel rejects the
    # FEXTRA header, so it has to be read through gzip.open rather than by path.
    assert _names(utils.open_hmm_file(helpers.bgzip_file(src, tmp_path / "r.hmm.bgz"))) == plain


def test_open_hmm_file_detects_compression_by_content_not_extension(shared_datadir, tmp_path):
    # easel dispatches on the extension and would fail on both of these.
    src = shared_datadir / "FeSOD_pfam.hmm"
    plain = _names(utils.open_hmm_file(src))
    misnamed_gz = helpers.gzip_file(src, tmp_path / "compressed.hmm")   # gzip content, no .gz
    assert _names(utils.open_hmm_file(misnamed_gz)) == plain
    uncompressed_but_named_gz = tmp_path / "plain.hmm.gz"
    uncompressed_but_named_gz.write_bytes(src.read_bytes())
    assert _names(utils.open_hmm_file(uncompressed_but_named_gz)) == plain


def test_open_hmm_file_rewind_and_close_compressed(shared_datadir, tmp_path):
    gz = helpers.gzip_file(shared_datadir / "FeSOD_pfam.hmm", tmp_path / "r.hmm.gz")
    handle = utils.open_hmm_file(gz)
    first = [utils.pyhmmer_decode(h.name) for h in handle]
    handle.rewind()
    assert [utils.pyhmmer_decode(h.name) for h in handle] == first
    handle.close()
    handle.close()  # must be idempotent: pyhmmer 0.12 double-frees on a second close


def test_open_hmm_file_closes_underlying_handle(shared_datadir, tmp_path):
    gz = helpers.gzip_file(shared_datadir / "FeSOD_pfam.hmm", tmp_path / "r.hmm.gz")
    with utils.open_hmm_file(gz) as handle:
        inner = handle._handle
        assert not inner.closed
    assert inner.closed


def test_open_hmm_file_passes_handles_through(shared_datadir):
    with open(shared_datadir / "FeSOD_pfam.hmm", "rb") as fh:
        assert _names(utils.open_hmm_file(fh))


def test_pressed_sidecars_are_ignored(shared_datadir, tmp_path):
    """An .hmm edited after hmmpress must not be read from its stale sidecars."""
    import pyhmmer
    target = tmp_path / "ref.hmm"
    target.write_bytes((shared_datadir / "FeSOD_pfam.hmm").read_bytes())
    pyhmmer.hmmer.hmmpress(list(utils.open_hmm_file(target)), str(target))
    assert (tmp_path / "ref.hmm.h3m").exists()

    # Replace the contents entirely; the sidecars now describe the old profiles.
    target.write_bytes((shared_datadir / "pdonr_hmms.hmm").read_bytes())
    expected = [utils.pyhmmer_decode(h.name) for h in utils.iter_hmms(shared_datadir / "pdonr_hmms.hmm")]
    assert _names(utils.open_hmm_file(target)) == expected
    assert list(utils.iter_hmm_names(target)) == expected


@pytest.mark.parametrize("path,expected", [
    ("a.hmm", "a"),
    ("a.hmm.gz", "a"),
    ("a.hmm.bgz", "a"),
    ("/p/some.long.name.hmm", "some.long.name"),
    ("refs.v1", "refs"),
    ("db.gz", "db"),
    ("plain", "plain"),
])
def test_db_name_from_path(path, expected):
    assert utils.db_name_from_path(path) == expected


def test_read_hmms_db_name_ignores_compression_suffix(shared_datadir, tmp_path):
    gz = helpers.gzip_file(shared_datadir / "FeSOD_pfam.hmm", tmp_path / "FeSOD_pfam.hmm.gz")
    assert set(utils.read_hmms([gz])) == {"FeSOD_pfam"}


def test_open_writable_hmm_file_rejects_bgzf(tmp_path):
    with pytest.raises(ValueError, match="BGZF"):
        utils.open_writable_hmm_file(tmp_path / "out.hmm.bgz")


def test_open_writable_hmm_file_writes_plain_gzip(shared_datadir, tmp_path):
    out = tmp_path / "out.hmm.gz"
    with utils.open_writable_hmm_file(out) as handle:
        for hmm in utils.iter_hmms(shared_datadir / "FeSOD_pfam.hmm"):
            hmm.write(handle)
    # Plain gzip, deliberately not BGZF: easel cannot open BGZF by path, so a
    # BGZF .hmm would be unreadable by hmmsearch and by Domainator itself.
    assert utils.detect_compression(out) == "gzip"
    import pyhmmer
    with pyhmmer.plan7.HMMFile(str(out)) as handle:  # the HMMER interop proof
        assert [utils.pyhmmer_decode(h.name) for h in handle] == \
            [utils.pyhmmer_decode(h.name) for h in utils.iter_hmms(shared_datadir / "FeSOD_pfam.hmm")]


@pytest.mark.parametrize("fixture", [
    "CcdB.hmm", "FeSOD_pfam.hmm", "Peptidase_M28.hmm", "SPR.hmm",
    "pdonr_hmms.hmm", "dna_profiles.hmm", "rna_profiles.hmm",
])
def test_iter_hmm_names_matches_parsed_names(shared_datadir, fixture):
    path = shared_datadir / fixture
    scanned = list(utils.iter_hmm_names(path))
    parsed = [utils.pyhmmer_decode(h.name) for h in utils.iter_hmms(path)]
    assert scanned == parsed


@pytest.mark.parametrize("fixture", ["CcdB.hmm", "Peptidase_M28.hmm"])
def test_committed_fixtures_with_sidecars_read_as_text(shared_datadir, fixture):
    """test/data ships .h3* sidecars for these two; they must not be read.

    They were committed in the repo's first commit and no code references them.
    Before hmm reads were pinned to db=False, pyhmmer preferred the sidecars, so
    these fixtures exercised the binary path everywhere they were used. Keeping
    them means this invariant is checked against files that really do have
    sidecars on disk, which is the situation users hit after running hmmpress.
    """
    import pyhmmer
    path = shared_datadir / fixture
    assert (shared_datadir / (fixture + ".h3m")).exists(), "fixture lost its sidecars"
    with pyhmmer.plan7.HMMFile(path) as stock:
        assert stock.is_pressed(), "fixture no longer pressed; this test is meaningless"
    with utils.open_hmm_file(path) as handle:
        assert not handle.is_pressed()
    assert [utils.pyhmmer_decode(h.name) for h in utils.iter_hmms(path)] == \
        [utils.pyhmmer_decode(h.name) for h in pyhmmer.plan7.HMMFile(path, db=False)]
