"""Tests for the reseek backend of the structure tools.

reseek differs from foldseek in ways the interface has to absorb -- no alignment types, no
structural metrics, no target sequence in the search output, a p-value instead of an
E-value, and a minimum chain length -- so these tests pin each of those, plus the
cross-backend agreement that makes the abstraction worth having.
"""
import contextlib
import shutil

import pytest

from domainator import structure_lib, structure_search, structure_domainate, structure_to_genbank, utils


def _has(binary):
    return shutil.which(binary) is not None


needs_reseek = pytest.mark.skipif(not _has("reseek"), reason="reseek is not installed")
needs_both = pytest.mark.skipif(not (_has("reseek") and _has("foldseek")),
                                reason="reseek and foldseek are both required")

# 1UBQ is a copy of the reference; 1NDD (NEDD8) is a genuine remote homolog of it.
EXPECTED_HITS = {"1UBQ_A", "1NDD_A"}
# reseek drops chains shorter than 32 residues, which removes all four insulin chains.
RESEEK_NAMES = {"1CRN_A", "1IGD_A", "1NDD_A", "1UBQ_A", "some.long.name_A"}


def _maybe_warns(algorithm):
    """reseek warns about the short chains it drops; foldseek keeps them."""
    if algorithm == "reseek":
        return pytest.warns(RuntimeWarning)
    return contextlib.nullcontext()


def _read(path):
    return list(utils.parse_seqfiles((str(path),), default_molecule_type="protein"))


def _features(record, feature_type):
    return [f for f in record.features if f.type == feature_type]


@pytest.fixture
def inputs(shared_datadir):
    return str(shared_datadir / "structures" / "inputs")


@pytest.fixture
def references(shared_datadir):
    return str(shared_datadir / "structures" / "refs")


# ---------------------------------------------------------------- capabilities (no binary)


def test_reseek_is_registered():
    assert "reseek" in structure_lib.ALIGNERS


def test_reseek_capability_declarations():
    """These declarations are what let the tools reject unsupported requests up front."""
    reseek = structure_lib.ReseekAligner
    assert reseek.supports_alignment_types == frozenset()
    assert reseek.supports_metrics == frozenset()
    assert reseek.supports_3di is False
    assert reseek.provides_target_sequence is False
    # reseek has no native cap, so --max_seqs is honored by keeping the best hits itself
    assert reseek.supports_max_seqs is True
    assert reseek.default_max_seqs is None
    # foldseek, by contrast, supports all of it
    assert structure_lib.FoldseekAligner.provides_target_sequence is True
    assert structure_lib.FoldseekAligner.supports_3di is True
    assert structure_lib.FoldseekAligner.default_alignment_type == 2
    assert structure_lib.FoldseekAligner.supports_max_seqs is True
    assert structure_lib.FoldseekAligner.default_max_seqs == 1000


def test_effective_max_seqs_defaults():
    """No cap for reseek unless asked for; foldseek keeps its own default of 1000."""
    reseek = structure_lib.ReseekAligner.__new__(structure_lib.ReseekAligner)
    assert reseek.effective_max_seqs(None) is None
    assert reseek.effective_max_seqs(50) == 50

    foldseek = structure_lib.FoldseekAligner.__new__(structure_lib.FoldseekAligner)
    assert foldseek.effective_max_seqs(None) == 1000
    assert foldseek.effective_max_seqs(50) == 50


def test_a_backend_without_max_seqs_support_refuses_it():
    """The capability is still enforced for any future backend that cannot honor a cap."""
    class _NoCap(structure_lib.StructureAligner):
        name = "nocap"
        supports_max_seqs = False

        def _resolve_bin(self, bin_path):
            return "nocap"

        def build_db(self, inputs, out_prefix):
            raise NotImplementedError

        def search(self, *args, **kwargs):
            raise NotImplementedError

        def iter_sequences(self, db):
            raise NotImplementedError

    aligner = _NoCap.__new__(_NoCap)
    assert aligner.effective_max_seqs(None) is None
    with pytest.raises(RuntimeError, match="no per-reference hit cap"):
        aligner.effective_max_seqs(50)


def test_apply_max_seqs_keeps_the_most_significant_hits(tmp_path):
    """The load-bearing property: reseek's output is not ordered by significance, so a cap
    taken from file order would discard the best hits."""
    aligner = structure_lib.ReseekAligner.__new__(structure_lib.ReseekAligner)
    columns = structure_lib.ReseekAligner.COLUMNS

    def row(target, pvalue, query="Q"):
        values = {c: "1" for c in columns}
        values.update(query=query, target=target, pvalue=str(pvalue))
        return "\t".join(values[c] for c in columns) + "\n"

    raw = tmp_path / "raw.tsv"
    # deliberately worst-first, mimicking reseek's unordered output
    raw.write_text(row("worst", 1e-2) + row("best", 1e-30)
                   + row("mid", 1e-10) + row("worse", 1e-3))

    kept = [line.split("\t")[1] for line in aligner._apply_max_seqs(str(raw), 2)]
    assert sorted(kept) == ["best", "mid"]


def test_apply_max_seqs_returns_best_first(tmp_path):
    """No file is written, so the returned order is what downstream sees."""
    aligner = structure_lib.ReseekAligner.__new__(structure_lib.ReseekAligner)
    columns = structure_lib.ReseekAligner.COLUMNS

    def row(target, pvalue):
        values = {c: "1" for c in columns}
        values.update(query="Q", target=target, pvalue=str(pvalue))
        return "\t".join(values[c] for c in columns) + "\n"

    raw = tmp_path / "raw.tsv"
    raw.write_text(row("mid", 1e-10) + row("worst", 1e-2) + row("best", 1e-30))
    kept = [line.split("\t")[1] for line in aligner._apply_max_seqs(str(raw), 3)]
    assert kept == ["best", "mid", "worst"]


def test_apply_max_seqs_is_per_query(tmp_path):
    aligner = structure_lib.ReseekAligner.__new__(structure_lib.ReseekAligner)
    columns = structure_lib.ReseekAligner.COLUMNS

    def row(query, target, pvalue):
        values = {c: "1" for c in columns}
        values.update(query=query, target=target, pvalue=str(pvalue))
        return "\t".join(values[c] for c in columns) + "\n"

    raw = tmp_path / "raw.tsv"
    raw.write_text(row("A", "a1", 1e-9) + row("A", "a2", 1e-8)
                   + row("B", "b1", 1e-7) + row("B", "b2", 1e-6))
    lines = aligner._apply_max_seqs(str(raw), 1)
    kept = {line.split("\t")[0]: line.split("\t")[1] for line in lines}
    # one hit per query, and the better one of each pair
    assert kept == {"A": "a1", "B": "b1"}


def test_reseek_db_header_parsing(tmp_path):
    """The chain count is read from the .bcb header, so E-values need no database pass."""
    good = tmp_path / "db.bcb"
    good.write_bytes(structure_lib.RESEEK_BCB_MAGIC.to_bytes(4, "little")
                     + (1234).to_bytes(8, "little") + b"rest of file")
    assert structure_lib.is_reseek_db(good)
    assert structure_lib.reseek_db_size(good) == 1234
    assert structure_lib.resolve_structure_input(str(good)) == ("db", [str(good)])

    bca = tmp_path / "old.bca"
    bca.write_bytes(structure_lib.RESEEK_BCA_MAGIC.to_bytes(4, "little")
                    + (7).to_bytes(8, "little") + b"x")
    assert structure_lib.reseek_db_size(bca) == 7


def test_reseek_db_detection_rejects_other_files(tmp_path, shared_datadir):
    not_a_db = tmp_path / "thing.bcb"
    not_a_db.write_bytes(b"HEADER    not a reseek database at all")
    assert not structure_lib.is_reseek_db(not_a_db)
    assert not structure_lib.is_reseek_db(shared_datadir / "structures" / "inputs" / "1UBQ.pdb.gz")
    with pytest.raises(RuntimeError, match="not a readable reseek database"):
        structure_lib.reseek_db_size(not_a_db)


def test_neg_log10_score_derivation():
    """reseek reports no bitscore, so score is -log10(p): higher must mean better."""
    assert structure_lib._neg_log10(1e-12) == pytest.approx(12.0)
    assert structure_lib._neg_log10(0.1) == pytest.approx(1.0)
    # an underflowed p-value must stay finite and still sort as the best hit
    assert structure_lib._neg_log10(0.0) > structure_lib._neg_log10(1e-300)


# ---------------------------------------------------------------- needs reseek


@pytest.fixture
def aligner():
    return structure_lib.ReseekAligner(cpu=2)


@needs_reseek
def test_evalue_threshold_is_converted_to_a_pvalue(aligner, tmp_path):
    """domainator's -e is an E-value against the searched database; reseek wants a
    p-value. E = p * N, so the threshold has to be divided by the database size."""
    db = tmp_path / "db.bcb"
    db.write_bytes(structure_lib.RESEEK_BCB_MAGIC.to_bytes(4, "little")
                   + (1000).to_bytes(8, "little"))
    target = structure_lib.StructureDB(prefix=str(db))
    pvalue, size = aligner._search_pvalue(0.001, target)
    assert size == 1000
    assert pvalue == pytest.approx(1e-6)
    # a p-value is a probability, so the conversion must never exceed 1
    assert aligner._search_pvalue(1e9, target)[0] == 1.0


@needs_reseek
def test_build_db_warns_about_short_chains(aligner, inputs, tmp_path):
    """The insulin chains are 21-30 residues, below reseek's minimum, and would otherwise
    vanish from the output with no explanation."""
    _, resolved = structure_lib.resolve_structure_inputs([inputs])
    with pytest.warns(RuntimeWarning, match=r"reseek skipped 4 chain\(s\) shorter than 32"):
        db = aligner.build_db(resolved, str(tmp_path / "db"))
    assert structure_lib.is_reseek_db(db.prefix)
    assert db.prefix.endswith(".bcb")
    assert aligner.database_has_coordinates(db)


@needs_reseek
def test_iter_sequences(aligner, inputs, tmp_path):
    _, resolved = structure_lib.resolve_structure_inputs([inputs])
    with pytest.warns(RuntimeWarning):
        db = aligner.build_db(resolved, str(tmp_path / "db"))
    entries = list(aligner.iter_sequences(db))
    assert {name for name, _, _ in entries} == RESEEK_NAMES
    lengths = {name: len(sequence) for name, sequence, _ in entries}
    assert lengths["1UBQ_A"] == 76 and lengths["1NDD_A"] == 74
    # reseek records no source-file map, so provenance is unavailable rather than wrong
    assert all(source is None for _, _, source in entries)


@needs_reseek
def test_iter_3di_unsupported(aligner, inputs, tmp_path):
    _, resolved = structure_lib.resolve_structure_inputs([inputs])
    with pytest.warns(RuntimeWarning):
        db = aligner.build_db(resolved, str(tmp_path / "db"))
    with pytest.raises(NotImplementedError, match="3Di"):
        list(aligner.iter_3di(db))


@needs_reseek
def test_search_hits(aligner, inputs, references, tmp_path):
    _, resolved_inputs = structure_lib.resolve_structure_inputs([inputs])
    _, resolved_refs = structure_lib.resolve_structure_inputs([references])
    with pytest.warns(RuntimeWarning):
        indb = aligner.build_db(resolved_inputs, str(tmp_path / "indb"))
    refdb = aligner.build_db(resolved_refs, str(tmp_path / "refdb"))
    work = tmp_path / "work"
    work.mkdir()

    hits = list(aligner.search(refdb, indb, evalue=0.001, max_seqs=None, want_tseq=False,
                               alignment_type=None, metrics=[], work_dir=str(work)))
    by_target = {hit.target: hit for hit in hits}
    assert set(by_target) == EXPECTED_HITS

    self_hit = by_target["1UBQ_A"]
    assert (self_hit.tstart, self_hit.tend, self_hit.tlen) == (1, 76, 76)
    assert (self_hit.qstart, self_hit.qend, self_hit.qlen) == (1, 76, 76)
    # pctid is already a percentage, so fident must be the fraction
    assert self_hit.fident == pytest.approx(1.0)
    assert self_hit.tseq == "", "reseek reports no target sequence"
    assert self_hit.bits > 0
    assert self_hit.evalue < 1e-6


@needs_reseek
def test_search_rejects_metrics(aligner, inputs, references, tmp_path):
    _, resolved = structure_lib.resolve_structure_inputs([references])
    db = aligner.build_db(resolved, str(tmp_path / "db"))
    work = tmp_path / "work"
    work.mkdir()
    with pytest.raises(RuntimeError, match="cannot compute"):
        list(aligner.search(db, db, evalue=1.0, max_seqs=None, want_tseq=False,
                            alignment_type=None, metrics=["lddt"], work_dir=str(work)))


@needs_reseek
def test_check_capabilities_rejects_alignment_type_and_metrics(aligner, tmp_path):
    db = tmp_path / "db.bcb"
    db.write_bytes(structure_lib.RESEEK_BCB_MAGIC.to_bytes(4, "little")
                   + (1).to_bytes(8, "little"))
    target = structure_lib.StructureDB(prefix=str(db))
    with pytest.raises(RuntimeError, match="no alignment-type option at all"):
        aligner.check_capabilities(1, [], [target])
    with pytest.raises(RuntimeError, match="cannot compute"):
        aligner.check_capabilities(None, ["tmscore"], [target])
    # the default (unset) alignment type must be accepted
    aligner.check_capabilities(None, [], [target])


# ---------------------------------------------------------------- through the CLIs


@needs_reseek
def test_structure_search_with_reseek(inputs, references, tmp_path):
    out = tmp_path / "hits.gb"
    with pytest.warns(RuntimeWarning):
        structure_search.main(["-i", inputs, "-r", references, "--algorithm", "reseek",
                               "-o", str(out), "--cpu", "2"])
    records = _read(out)
    assert {record.id for record in records} == EXPECTED_HITS
    for record in records:
        feature = _features(record, "Domain_Search")[0]
        assert feature.qualifiers["program"] == ["reseek"]
        # records are built from database sequences, since reseek reports none
        assert len(record) == {"1UBQ_A": 76, "1NDD_A": 74}[record.id]
        # the reference is 76 aa regardless of which input record was hit
        assert int(feature.qualifiers["rlen"][0]) == 76
        assert 0 <= int(feature.location.start) < int(feature.location.end) <= len(record)
        for metric in ("tmscore", "lddt", "rmsd", "prob"):
            assert metric not in feature.qualifiers


@needs_reseek
def test_structure_search_max_hits_with_reseek(inputs, references, tmp_path):
    out = tmp_path / "one.gb"
    with pytest.warns(RuntimeWarning):
        structure_search.main(["-i", inputs, "-r", references, "--algorithm", "reseek",
                               "-o", str(out), "--max_hits", "1", "--cpu", "2"])
    records = _read(out)
    assert len(records) == 1
    assert records[0].id == "1UBQ_A"


@needs_reseek
def test_structure_domainate_with_reseek(inputs, references, tmp_path):
    out = tmp_path / "annotated.gb"
    with pytest.warns(RuntimeWarning):
        structure_domainate.main(["-i", inputs, "-r", references, "--algorithm", "reseek",
                                  "-o", str(out), "--cpu", "2"])
    records = _read(out)
    # the short insulin chains are absent because reseek dropped them
    assert {record.id for record in records} == RESEEK_NAMES
    annotated = {record.id for record in records if _features(record, "Domainator")}
    assert annotated == EXPECTED_HITS


@needs_reseek
def test_structure_to_genbank_with_reseek(inputs, tmp_path):
    out = tmp_path / "chains.gb"
    with pytest.warns(RuntimeWarning):
        structure_to_genbank.main(["-i", inputs, "--algorithm", "reseek",
                                   "-o", str(out), "--cpu", "2"])
    assert {record.id for record in _read(out)} == RESEEK_NAMES


@needs_reseek
def test_store_3di_rejected_for_reseek(inputs, tmp_path):
    with pytest.raises(RuntimeError, match="--store_3di is not available"):
        structure_to_genbank.main(["-i", inputs, "--algorithm", "reseek",
                                   "--store_3di", "-o", str(tmp_path / "x.gb")])


@needs_reseek
def test_reseek_stats_level_changes_significance(inputs, references, tmp_path):
    """--reseek_stats moves p-values by orders of magnitude, so it changes what passes."""
    counts = {}
    for level in ("family", "superfamily"):
        out = tmp_path / f"{level}.gb"
        with pytest.warns(RuntimeWarning):
            structure_search.main(["-i", inputs, "-r", references, "--algorithm", "reseek",
                                   "--reseek_stats", level, "-e", "10",
                                   "-o", str(out), "--cpu", "2"])
        counts[level] = len(_read(out))
    # the family-level model is less stringent, so it admits at least as many hits
    assert counts["family"] >= counts["superfamily"]


@needs_reseek
def test_max_seqs_end_to_end_keeps_best_hits(inputs, references, tmp_path):
    """At the family level all four fixtures hit, and reseek reports them out of order, so
    a --max_seqs of 2 has to return the two most significant of them."""
    common = ["-i", inputs, "-r", references, "--algorithm", "reseek",
              "--reseek_stats", "family", "-e", "100000", "--cpu", "2"]

    uncapped = tmp_path / "all.gb"
    with pytest.warns(RuntimeWarning):
        structure_search.main(common + ["-o", str(uncapped)])
    ranked = sorted(
        ((float(_features(r, "Domain_Search")[0].qualifiers["evalue"][0]), r.id)
         for r in _read(uncapped)))
    assert len(ranked) == 4
    best_two = {name for _, name in ranked[:2]}

    capped = tmp_path / "capped.gb"
    with pytest.warns(RuntimeWarning):
        structure_search.main(common + ["--max_seqs", "2", "-o", str(capped)])
    assert {r.id for r in _read(capped)} == best_two


@needs_reseek
def test_max_seqs_saturation_warns_for_reseek(inputs, references, tmp_path):
    with pytest.warns(RuntimeWarning, match=r"--max_seqs \(1\) was reached"):
        structure_search.main(["-i", inputs, "-r", references, "--algorithm", "reseek",
                               "--reseek_stats", "family", "-e", "100000",
                               "--max_seqs", "1", "-o", str(tmp_path / "x.gb"), "--cpu", "2"])


@needs_reseek
def test_hits_tsv_for_reseek(inputs, references, tmp_path):
    """Each backend names its own hit table, so --hits_tsv has to know which to look for."""
    out = tmp_path / "hits.gb"
    tsv = tmp_path / "raw.tsv"
    with pytest.warns(RuntimeWarning):
        structure_search.main(["-i", inputs, "-r", references, "--algorithm", "reseek",
                               "--hits_tsv", str(tsv), "-o", str(out), "--cpu", "2"])
    lines = [line for line in tsv.read_text().split("\n") if line]
    assert len(lines) == 2
    assert {line.split("\t")[1] for line in lines} == EXPECTED_HITS


@needs_reseek
def test_prebuilt_reseek_database_as_input(inputs, references, tmp_path):
    prefix = tmp_path / "kept"
    first = tmp_path / "first.gb"
    with pytest.warns(RuntimeWarning):
        structure_search.main(["-i", inputs, "-r", references, "--algorithm", "reseek",
                               "--keep_db", str(prefix), "-o", str(first), "--cpu", "2"])
    second = tmp_path / "second.gb"
    structure_search.main(["-i", f"{prefix}.bcb", "-r", references, "--algorithm", "reseek",
                           "-o", str(second), "--cpu", "2"])
    assert ([(r.id, str(r.seq)) for r in _read(first)]
            == [(r.id, str(r.seq)) for r in _read(second)])


@needs_reseek
def test_database_format_mismatch_is_explained(shared_datadir, references, tmp_path):
    with pytest.raises(RuntimeError, match="not a reseek database"):
        structure_search.main(["-i", str(shared_datadir / "foldseek" / "FeSOD"),
                               "-r", references, "--algorithm", "reseek",
                               "-o", str(tmp_path / "x.gb")])


@needs_reseek
def test_missing_reseek_binary(inputs, references, tmp_path, monkeypatch):
    monkeypatch.setattr(structure_lib.shutil, "which", lambda name: None)
    with pytest.raises(RuntimeError, match=r"--reseek_path"):
        structure_search.main(["-i", inputs, "-r", references, "--algorithm", "reseek",
                               "-o", str(tmp_path / "x.gb")])


# ---------------------------------------------------------------- cross-backend agreement


@needs_both
def test_backends_agree_on_hits_and_coordinates(inputs, references, tmp_path):
    """The two backends are independent implementations, so agreeing on which structures
    are homologous and where the alignment sits is the real check on the abstraction."""
    results = {}
    for algorithm in ("foldseek", "reseek"):
        out = tmp_path / f"{algorithm}.gb"
        with _maybe_warns(algorithm):
            structure_search.main(["-i", inputs, "-r", references,
                                   "--algorithm", algorithm, "-o", str(out), "--cpu", "2"])
        results[algorithm] = {
            record.id: (
                len(record),
                int(_features(record, "Domain_Search")[0].location.start),
                int(_features(record, "Domain_Search")[0].location.end),
                int(_features(record, "Domain_Search")[0].qualifiers["rlen"][0]),
                round(float(_features(record, "Domain_Search")[0].qualifiers["identity"][0])),
            )
            for record in _read(out)
        }

    assert set(results["foldseek"]) == EXPECTED_HITS
    assert set(results["reseek"]) == EXPECTED_HITS
    assert results["foldseek"] == results["reseek"]
