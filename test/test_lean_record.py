"""Tests for the lean record model and its conversion boundary (fast_parser_plan.md).

Two-tier model: the native (Rust) ``_gbfast`` parser (the "lean" backend) and the
Biopython parser. These tests validate that the native lean path produces records
equivalent to Biopython, and that the lean<->SeqRecord conversions are lossless.
"""
import glob
import warnings

import pytest

from domainator import utils
from domainator import lean_record
from domainator.lean_record import (
    LeanContig,
    LeanFeature,
    seqrecord_to_lean,
    lean_to_seqrecord,
)
from helpers import compare_seqrecords

native_only = pytest.mark.skipif(
    lean_record._gbfast is None, reason="native _gbfast extension not built"
)


def _feat_key(f):
    return (
        f.type,
        [(int(p.start), int(p.end), p.strand) for p in f.location.parts],
        {k: list(v) for k, v in f.qualifiers.items()},
    )


def _feat_key_lean(f):
    return (f.type, f.parts, f.operator, f.between, {k: list(v) for k, v in f.qualifiers.items()})


def _collect_partial(producer):
    """Collect the prefix a native producer yields before any early-stop fallback."""
    out = []
    try:
        for x in producer:
            out.append(x)
    except lean_record.LeanParseError:
        pass
    return out


def _native_lean_seqrecords(path):
    """Materialize the native search path to SeqRecords (post-swap), prefix only."""
    out = []
    for sc in _collect_partial(lean_record.iter_lean_search_native(str(path), 0, -1, "protein")):
        out.append(lean_to_seqrecord(lean_record.search_contig_to_lean(sc, set(), "protein")))
    return out


# Annotations the lean model intentionally carries (it drops comment/references/etc.).
_LEAN_ANNOTATION_KEYS = {
    "molecule_type", "topology", "data_file_division", "accessions",
    "sequence_version", "keywords", "source", "organism", "taxonomy",
}


@native_only
def test_native_lean_matches_biopython(shared_datadir):
    """The native lean path must match Biopython feature-for-feature and on the
    lean-carried annotations, for the records it parses (others fall back)."""
    warnings.simplefilter("ignore")
    checked = 0
    for fn in sorted(glob.glob(str(shared_datadir / "*.gb"))):
        try:
            bp = list(utils.parse_seqfiles([fn], default_molecule_type="protein", genbank_parser="biopython"))
        except Exception:
            continue
        native = _native_lean_seqrecords(fn)
        for a, b in zip(native, bp):  # zip stops at the native prefix
            assert a.id == b.id
            assert str(a.seq) == str(b.seq)
            assert [_feat_key(f) for f in a.features] == [_feat_key(f) for f in b.features]
            for key in _LEAN_ANNOTATION_KEYS:
                if key in b.annotations:
                    assert a.annotations.get(key) == b.annotations[key], (fn, key)
            checked += 1
    assert checked > 0


def test_seqrecord_lean_roundtrip(shared_datadir):
    """SeqRecord -> lean -> SeqRecord (the fallback path) must be lossless."""
    warnings.simplefilter("ignore")
    for fn in sorted(glob.glob(str(shared_datadir / "*.gb"))):
        try:
            recs = list(utils.parse_seqfiles([fn], genbank_parser="biopython"))
        except Exception:
            continue
        for r in recs:
            r.annotations.setdefault("molecule_type", "protein")
            back = lean_to_seqrecord(seqrecord_to_lean(r))
            compare_seqrecords(r, back)


def test_lean_feature_len_and_span():
    f = LeanFeature("CDS", ((10, 40, -1, False, False),), None, False, {"locus_tag": ["x"]})
    assert len(f) == 30
    assert f.start == 10 and f.end == 40 and f.strand == -1
    cf = LeanFeature("source", ((0, 100, 1, False, False), (200, 250, 1, False, False)), "join", False, {})
    assert len(cf) == 150
    assert cf.start == 0 and cf.end == 250


def test_swap_name_id_on_lean():
    c = LeanContig("ACC.1", "LOCUSNAME", "d", "ACGT", {"molecule_type": "DNA"}, [])
    utils.swap_name_id(c)
    assert c.id == "LOCUSNAME" and c.name == "ACC.1"


# GenBank location strings (adapted from Biopython Tests/test_SeqIO_features.py).
# The native lean parser must produce the same location parts as Biopython.
_LOCATION_CASES = [
    "6..10",
    "complement(6..10)",
    "5^6",
    "<6..10",
    "6..>10",
    "<6..>10",
    "join(6..10,13..15)",
    "complement(join(6..10,13..15))",
    "order(6..10,13..15)",
    "join(6..10,complement(13..15))",  # mixed-strand (trans-spliced)
]


@native_only
@pytest.mark.parametrize("loc", _LOCATION_CASES)
def test_native_lean_location_matches_biopython(loc, tmp_path):
    import io
    template = (
        "LOCUS       T                         60 bp    DNA     linear   UNK 01-JAN-2024\n"
        "FEATURES             Location/Qualifiers\n"
        "     misc_feature    {loc}\n"
        '                     /note="x"\n'
        "ORIGIN\n"
        "        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt\n"
        "//\n"
    ).format(loc=loc)
    gb_path = tmp_path / "loc.gb"
    gb_path.write_text(template)

    bp = list(utils.parse_seqfiles([io.StringIO(template)], filetype_override="genbank", genbank_parser="biopython"))[0]
    native = _native_lean_seqrecords(gb_path)
    assert native, "native parser unexpectedly fell back on this location case"
    bpf = [f for f in bp.features if f.type == "misc_feature"][0]
    lf = [f for f in native[0].features if f.type == "misc_feature"][0]
    assert [(int(p.start), int(p.end), p.strand) for p in bpf.location.parts] == \
           [(int(p.start), int(p.end), p.strand) for p in lf.location.parts]


@native_only
def test_lean_search_native_matches_parse_lean(shared_datadir):
    """The peptides-only search path must materialize to the SAME record the native
    per-feature parser builds (both are the Rust record builder), and its
    cds_peptides indices must point at CDS features in the materialized record."""
    for fn in ("206.gb", "pDONR201.gb", "Staph_phages.gb", "MT_nbs.gb"):
        path = str(shared_datadir / fn)
        ref = _collect_partial(lean_record.iter_lean_genbank_native(path, 0, -1, "protein"))
        search = _collect_partial(lean_record.iter_lean_search_native(path, 0, -1, "protein"))
        n = min(len(ref), len(search))
        if n == 0:
            continue
        for i in range(n):
            mat = lean_record.search_contig_to_lean(search[i], set(), "protein")
            # parse_lean yields pre-swap; the search path is post-swap.
            r = ref[i]
            r.id, r.name = r.name, r.id
            assert mat.id == r.id
            assert mat.seq == r.seq
            assert [_feat_key_lean(f) for f in mat.features] == [_feat_key_lean(f) for f in r.features]
            for fidx, pep in search[i].cds_peptides(set()):
                assert mat.features[fidx].type == "CDS"


@native_only
def test_native_translation_matches_biopython(shared_datadir):
    """Rust CDS translation (used for search peptides when /translation is absent
    or empty) must match Biopython's feature.translate, so hits are identical."""
    # pDONR201's ccdB CDS stores an empty /translation, so the native path
    # translates it in Rust; compare every CDS peptide to Biopython.
    for fn in ("pDONR201.gb",):
        path = str(shared_datadir / fn)
        search = _collect_partial(lean_record.iter_lean_search_native(path, 0, -1, "nucleotide"))
        bp = list(utils.parse_seqfiles([path], genbank_parser="biopython"))
        for sc, rec in zip(search, bp):
            rust_peps = dict(sc.cds_peptides(set()))
            for j, f in enumerate(rec.features):
                if f.type == "CDS" and "pseudo" not in f.qualifiers and "pseudogene" not in f.qualifiers:
                    tr = f.qualifiers.get("translation")
                    if tr and tr[0]:
                        bp_pep = "".join(c for c in tr[0] if c.isalnum())
                    else:
                        bp_pep = "".join(c for c in str(f.translate(rec.seq, cds=False)) if c.isalnum())
                    assert rust_peps.get(j) == bp_pep


@native_only
def test_native_fasta_matches_biopython(shared_datadir):
    """The native FASTA parser must materialize to records equal to Biopython's
    (id == name == first token, full-title description, case-preserved seq,
    annotations == {molecule_type})."""
    for fn, mt in [("FeSOD_20.fasta", "protein"), ("simple_dna_target.fna", "nucleotide")]:
        path = str(shared_datadir / fn)
        bp = list(utils.parse_seqfiles([path], default_molecule_type=mt, genbank_parser="biopython"))
        native = list(lean_record.iter_lean_fasta_native(path, 0, -1, mt))
        assert len(native) == len(bp) > 0
        for sc, b in zip(native, bp):
            sr = lean_to_seqrecord(lean_record.search_fasta_to_lean(sc, set(), mt))
            assert sr.id == b.id and sr.name == b.name
            assert sr.description == b.description
            assert str(sr.seq) == str(b.seq)
            assert dict(sr.annotations) == dict(b.annotations)
            assert sc.cds_peptides(set()) == []  # FASTA has no CDS features


@native_only
def test_native_bgzf_matches_uncompressed(shared_datadir, tmp_path):
    """Parsing a BGZF copy via the native parser (virtual offsets) must yield the
    same records as parsing the uncompressed original."""
    from domainator.Bio import bgzf
    src = shared_datadir / "206.gb"
    bgz = tmp_path / "206.gb.bgz"
    with open(src, "rb") as fh, bgzf.BgzfWriter(str(bgz)) as w:
        w.write(fh.read())

    plain = _collect_partial(lean_record.iter_lean_search_native(str(src), 0, -1, "protein"))
    comp = _collect_partial(lean_record.iter_lean_search_native(str(bgz), 0, -1, "protein", compressed=True))
    assert len(plain) == len(comp) > 0
    for a, b in zip(plain, comp):
        ma = lean_record.search_contig_to_lean(a, set(), "protein")
        mb = lean_record.search_contig_to_lean(b, set(), "protein")
        assert ma.id == mb.id
        assert ma.seq == mb.seq
        assert [_feat_key_lean(f) for f in ma.features] == [_feat_key_lean(f) for f in mb.features]


@native_only
def test_native_per_cds_taxid_matches_python(shared_datadir):
    """The Rust per-CDS taxid (longest covering source) must equal the Python
    location_taxid on the materialized record, on a gappy multi-source genome
    (Staph_phages.gb: a host `order(...)` source + several prophage sources)."""
    from domainator import utils
    sc = list(lean_record.iter_lean_search_native(str(shared_datadir / "Staph_phages.gb"), 0, -1, "nucleotide"))[0]
    rust = {idx: (32644 if t is None else t) for idx, _, t in sc.cds_peptides_with_taxid(set())}
    rec = lean_record.lean_to_seqrecord(lean_record.materialize_lean_search(sc, None))
    sources = utils.get_sources(rec)
    checked = 0
    for j, f in enumerate(rec.features):
        if f.type == "CDS" and "pseudo" not in f.qualifiers and "pseudogene" not in f.qualifiers and j in rust:
            py = utils.location_taxid(utils.feature_intervals(f), sources)
            assert (32644 if py is None else py) == rust[j]
            checked += 1
    assert checked > 100
    # Per-CDS resolution actually distinguishes host from prophages (the bug fix):
    taxids = set(rust.values())
    assert 198466 in taxids and 198538 in taxids and len(taxids) > 2
    # No single source covers the whole contig -> gappy, needs per-region nucleotide filtering.
    assert sc.whole_contig_taxid() is None


@native_only
def test_native_whole_contig_taxid_single_source(shared_datadir):
    """A record whose single source covers the whole contig reports that taxid for the
    bypass test (206.gb has one source covering the contig)."""
    sc = list(lean_record.iter_lean_search_native(str(shared_datadir / "206.gb"), 0, -1, "nucleotide"))[0]
    whole = sc.whole_contig_taxid()
    assert whole is not None and whole == (sc.taxid())  # single source -> per-record == whole-contig


@native_only
def test_native_fasta_taxid(shared_datadir):
    """LeanFastaContig.taxid parses ` OX=` from the description, matching get_taxid."""
    from domainator import utils
    recs = list(lean_record.iter_lean_fasta_native(str(shared_datadir / "swissprot_CuSOD_subset.fasta"), 0, -1, "protein"))
    assert recs
    assert any(r.taxid() is not None for r in recs)
    for r in recs:
        got = r.taxid()
        assert (32644 if got is None else got) == utils.get_taxid(r)


@native_only
def test_get_prot_list_genbank_per_cds_filter(shared_datadir):
    """get_prot_list with allowed_taxids yields only CDSs whose longest-covering-source
    taxid is allowed -- per CDS, so a prophage filter selects only prophage CDSs."""
    from domainator.domainate import get_prot_list
    sc = list(lean_record.iter_lean_search_native(str(shared_datadir / "Staph_phages.gb"), 0, -1, "nucleotide"))[0]
    taxid_of = {idx: (32644 if t is None else t) for idx, _, t in sc.cds_peptides_with_taxid(set())}

    phage_names = [name for name, _ in get_prot_list(sc, 0, allowed_taxids=frozenset({198538}))]
    assert phage_names
    for name in phage_names:
        assert taxid_of[int(name.split(",")[1])] == 198538

    host_names = [name for name, _ in get_prot_list(sc, 0, allowed_taxids=frozenset({198466}))]
    assert len(host_names) > len(phage_names)  # host has far more CDSs than one prophage
    # no filter -> every searchable CDS
    assert len(list(get_prot_list(sc, 0, allowed_taxids=None))) == len(taxid_of)


@native_only
def test_nucleotide_contig_searchable_bypass(shared_datadir):
    """Pre-search nucleotide decision: single-source contig -> searched iff its taxid is
    allowed (bypass); gappy multi-source contig -> searched iff some source (or 32644) is
    allowed, else skipped pre-search."""
    from domainator.domainate import _nucleotide_contig_searchable
    staph = list(lean_record.iter_lean_search_native(str(shared_datadir / "Staph_phages.gb"), 0, -1, "nucleotide"))[0]
    assert staph.whole_contig_taxid() is None  # gappy
    assert _nucleotide_contig_searchable(staph, frozenset({198538}))      # a prophage source allowed -> search
    assert not _nucleotide_contig_searchable(staph, frozenset({999999}))  # no source / 32644 allowed -> skip
    assert _nucleotide_contig_searchable(staph, frozenset({32644}))       # unidentified allowed -> search (uncovered regions possible)
    # single whole-contig source (206.gb) -> bypass: searchable iff that taxid is allowed
    s206 = list(lean_record.iter_lean_search_native(str(shared_datadir / "206.gb"), 0, -1, "nucleotide"))[0]
    whole = s206.whole_contig_taxid()
    assert whole is not None
    assert _nucleotide_contig_searchable(s206, frozenset({whole}))
    assert not _nucleotide_contig_searchable(s206, frozenset({whole + 1}))
