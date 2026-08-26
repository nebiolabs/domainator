"""Tests for domainator.structure_lib.

The pure-Python half (input resolution, naming, binary lookup) is deliberately not guarded
by a foldseek skipif, so that logic stays covered on machines without foldseek installed.
Anything that shells out is guarded.
"""
import gzip
import os
import shutil
import tempfile

import pytest

from domainator import structure_lib


def _has_foldseek():
    return shutil.which("foldseek") is not None


needs_foldseek = pytest.mark.skipif(not _has_foldseek(), reason="foldseek is not installed")


# ---------------------------------------------------------------- pure python


@pytest.mark.parametrize("path,expected", [
    ("a.pdb", True),
    ("a.cif", True),
    ("a.mmcif", True),
    ("a.pdb.gz", True),
    ("a.cif.gz", True),
    ("some.long.name.pdb.gz", True),
    ("/x/y/a.ent", True),
    ("a.gb", False),
    ("a.fasta", False),
    ("a.hmm", False),
    ("a.gb.gz", False),
    ("noextension", False),
])
def test_is_structure_path(path, expected):
    assert structure_lib.is_structure_path(path) is expected


def test_resolve_structure_input_db_prefix(shared_datadir):
    """A database prefix is not a file, so it has to be recognized before the file checks."""
    prefix = str(shared_datadir / "foldseek" / "FeSOD")
    assert structure_lib.is_aligner_db(prefix)
    assert structure_lib.resolve_structure_input(prefix) == ("db", [prefix])


def test_resolve_structure_input_directory(shared_datadir):
    kind, resolved = structure_lib.resolve_structure_input(str(shared_datadir / "structures" / "inputs"))
    assert kind == "structures"
    assert len(resolved) == 6
    assert all(structure_lib.is_structure_path(p) for p in resolved)


def test_resolve_structure_input_single_file(shared_datadir):
    path = str(shared_datadir / "structures" / "refs" / "UBQref.pdb.gz")
    assert structure_lib.resolve_structure_input(path) == ("structures", [path])


def test_resolve_structure_input_glob(shared_datadir):
    """Globs must be expanded internally: --config values are never shell-expanded."""
    kind, resolved = structure_lib.resolve_structure_input(
        str(shared_datadir / "structures" / "inputs" / "1*.pdb.gz"))
    assert kind == "structures"
    assert len(resolved) == 5  # everything but some.long.name.pdb.gz


def test_resolve_structure_input_rejects_non_structure(shared_datadir):
    with pytest.raises(structure_lib.StructureInputError, match="not a recognized structure file"):
        structure_lib.resolve_structure_input(str(shared_datadir / "FeSOD_20.fasta"))


def test_resolve_structure_input_missing():
    with pytest.raises(structure_lib.StructureInputError, match="Could not interpret"):
        structure_lib.resolve_structure_input("/nonexistent/path/thing")


def test_resolve_structure_input_empty_glob(shared_datadir):
    with pytest.raises(structure_lib.StructureInputError, match="matched no structure files"):
        structure_lib.resolve_structure_input(str(shared_datadir / "structures" / "*.xyz"))


def test_resolve_structure_input_empty_directory(tmp_path):
    (tmp_path / "empty").mkdir()
    with pytest.raises(structure_lib.StructureInputError, match="No structure files found"):
        structure_lib.resolve_structure_input(str(tmp_path / "empty"))


def test_resolve_structure_inputs_rejects_mixing(shared_datadir):
    with pytest.raises(structure_lib.StructureInputError, match="Cannot mix"):
        structure_lib.resolve_structure_inputs([
            str(shared_datadir / "foldseek" / "FeSOD"),
            str(shared_datadir / "structures" / "refs" / "UBQref.pdb.gz"),
        ])


def test_resolve_structure_inputs_deduplicates(shared_datadir):
    """Overlapping globs and directories make duplicate paths easy to produce."""
    path = str(shared_datadir / "structures" / "refs" / "UBQref.pdb.gz")
    kind, resolved = structure_lib.resolve_structure_inputs([path, path])
    assert (kind, resolved) == ("structures", [path])


def test_database_label_keeps_dotted_prefixes():
    """Path.stem would collapse refs.v1 and refs.v2 onto the same label."""
    assert structure_lib.database_label("/a/b/refs.v1") == "refs.v1"
    assert structure_lib.database_label("/a/b/refs.v2") == "refs.v2"
    assert structure_lib.database_label("/a/b/refs.v1") != structure_lib.database_label("/a/b/refs.v2")


def test_database_label_strips_structure_extensions():
    assert structure_lib.database_label("/a/b/UBQref.pdb.gz") == "UBQref"
    assert structure_lib.database_label("/a/b/UBQref.pdb") == "UBQref"


@pytest.mark.parametrize("name,expected", [
    ("1abc.pdb_A", "1abc_A"),
    ("x.cif_B", "x_B"),
    ("1abc.pdb", "1abc"),
    ("plain_A", "plain_A"),
    ("some.long.name_A", "some.long.name_A"),
])
def test_strip_structure_extension(name, expected):
    assert structure_lib.strip_structure_extension(name) == expected


@pytest.mark.parametrize("path,expected", [
    ("x.pdb", "x"),
    ("x.pdb.gz", "x"),
    ("some.long.name.pdb", "some.long"),
    ("some.long.name.pdb.gz", "some.long"),
    ("/a/b/1UBQ.pdb.gz", "1UBQ"),
])
def test_foldseek_source_key(path, expected):
    """Must match how foldseek truncates names in <prefix>.source, or provenance is lost."""
    assert structure_lib.FoldseekAligner._source_key(path) == expected


def test_resolve_binary_missing(monkeypatch):
    monkeypatch.setattr(structure_lib.shutil, "which", lambda name: None)
    with pytest.raises(RuntimeError, match=r"--foldseek_path"):
        structure_lib.resolve_binary("foldseek", None)


def test_resolve_binary_rejects_non_executable(tmp_path):
    bogus = tmp_path / "notabinary"
    bogus.write_text("hello")
    with pytest.raises(RuntimeError, match="not an executable file"):
        structure_lib.resolve_binary("foldseek", str(bogus))


def test_build_protein_record():
    record = structure_lib.build_protein_record("1UBQ_A", "MQIFV", "/x/1UBQ.pdb")
    assert record.id == "1UBQ_A" and record.name == "1UBQ_A"
    assert str(record.seq) == "MQIFV"
    assert record.annotations["molecule_type"] == "protein"
    assert "structure=/x/1UBQ.pdb" in record.description
    assert not record.features


def test_build_protein_record_with_3di():
    record = structure_lib.build_protein_record("x", "MQIFV", None, "DWEWE")
    assert [f.type for f in record.features] == ["misc_feature"]
    assert record.features[0].qualifiers["threedi"] == ["DWEWE"]


def test_build_protein_record_rejects_mismatched_3di():
    """A wrong-length 3Di string must warn rather than be silently stored."""
    with pytest.warns(UserWarning, match="does not match sequence length"):
        record = structure_lib.build_protein_record("x", "MQIFV", None, "DW")
    assert not record.features


def test_warn_on_saturation():
    with pytest.warns(RuntimeWarning, match=r"--max_seqs \(10\) was reached"):
        structure_lib.warn_on_saturation({"ref1": 10, "ref2": 3}, 10)


def test_no_warning_when_not_saturated():
    import warnings as warnings_module
    with warnings_module.catch_warnings():
        warnings_module.simplefilter("error")
        structure_lib.warn_on_saturation({"ref1": 9}, 10)


def test_structure_hits_to_search_results_uses_target_coordinates():
    """References are the query and inputs the target, so start/end come from t* and
    rstart/rend/rlen from q*. This is inverted relative to domainate's older foldseek path,
    so it is asserted directly on a synthetic hit."""
    hit = structure_lib.StructureHit(
        query="REF_A", target="IN_A",
        qstart=5, qend=40, qlen=100,
        tstart=11, tend=60, tlen=200,
        evalue=1e-9, bits=123.0, fident=0.75, alnlen=50,
        qheader="REF_A some description",
    )
    grouped = structure_lib.structure_hits_to_search_results([hit], "mydb", evalue=1.0)
    assert list(grouped) == ["IN_A"]
    result = grouped["IN_A"][0]
    assert (result.start, result.end) == (10, 60)      # on the input, 0-based half open
    assert (result.rstart, result.rend, result.rlen) == (5, 40, 100)  # on the reference
    assert result.name == "REF_A"
    assert result.desc == "some description"
    assert result.identity == pytest.approx(75.0)      # fident is a fraction
    assert result.database == "mydb"
    assert result.program == "foldseek"


def test_structure_hits_to_search_results_applies_both_evalue_bounds():
    def hit(evalue):
        return structure_lib.StructureHit(
            query="R", target="T", qstart=1, qend=2, qlen=2, tstart=1, tend=2, tlen=2,
            evalue=evalue, bits=1.0, fident=1.0, alnlen=2)

    hits = [hit(1e-20), hit(1e-5), hit(10.0)]
    grouped = structure_lib.structure_hits_to_search_results(
        hits, "db", evalue=1.0, min_evalue=1e-10)
    assert len(grouped["T"]) == 1
    assert grouped["T"][0].evalue == 1e-5


def test_structure_hits_to_search_results_handles_descending_target_range():
    """FeatureLocation rejects start > end, so a reversed range must be normalized."""
    hit = structure_lib.StructureHit(
        query="R", target="T", qstart=1, qend=10, qlen=10, tstart=50, tend=20, tlen=100,
        evalue=1e-9, bits=1.0, fident=1.0, alnlen=10)
    result = structure_lib.structure_hits_to_search_results([hit], "db", evalue=1.0)["T"][0]
    assert result.start < result.end
    assert (result.start, result.end) == (19, 50)


def test_group_hits_by_target_streams_groups():
    def hit(target, evalue):
        return structure_lib.StructureHit(
            query="R", target=target, qstart=1, qend=2, qlen=2, tstart=1, tend=2, tlen=2,
            evalue=evalue, bits=1.0, fident=1.0, alnlen=2)

    # target-sorted, as StructureAligner.search(sort_by_target=True) guarantees
    hits = [hit("A", 1e-9), hit("A", 1e-8), hit("B", 1e-7)]
    groups = list(structure_lib.group_hits_by_target(hits, "db", evalue=1.0))
    assert [(target, len(results)) for target, results, _ in groups] == [("A", 2), ("B", 1)]


def test_group_hits_by_target_skips_filtered_targets():
    def hit(target, evalue):
        return structure_lib.StructureHit(
            query="R", target=target, qstart=1, qend=2, qlen=2, tstart=1, tend=2, tlen=2,
            evalue=evalue, bits=1.0, fident=1.0, alnlen=2)

    hits = [hit("A", 10.0), hit("B", 1e-7)]  # A is above the threshold
    groups = list(structure_lib.group_hits_by_target(hits, "db", evalue=1.0))
    assert [target for target, _, _ in groups] == ["B"]


# ---------------------------------------------------------------- needs foldseek


@pytest.fixture
def aligner():
    return structure_lib.FoldseekAligner(cpu=2)


@needs_foldseek
def test_build_db_and_iter_sequences(aligner, shared_datadir, tmp_path):
    _, inputs = structure_lib.resolve_structure_inputs([str(shared_datadir / "structures" / "inputs")])
    db = aligner.build_db(inputs, str(tmp_path / "indb"))
    assert db.owned
    assert db.has_coordinates(), "createdb from real structures must produce a _ca sub-database"

    entries = list(aligner.iter_sequences(db))
    names = [name for name, _, _ in entries]
    # one record per chain: 1ZNI contributes four
    assert names == ["1CRN_A", "1IGD_A", "1NDD_A", "1UBQ_A",
                     "1ZNI_A", "1ZNI_B", "1ZNI_C", "1ZNI_D", "some.long.name_A"]
    lengths = {name: len(sequence) for name, sequence, _ in entries}
    assert lengths["1UBQ_A"] == 76 and lengths["1NDD_A"] == 74 and lengths["1CRN_A"] == 46
    for _, sequence, _ in entries:
        assert sequence and set(sequence) <= set("ACDEFGHIKLMNPQRSTVWYXBZUO")
    # provenance resolves for every entry, including the dotted filename
    assert all(source is not None for _, _, source in entries)
    assert os.path.basename(dict((n, s) for n, _, s in entries)["some.long.name_A"]) == "some.long.name.pdb.gz"


@needs_foldseek
def test_iter_3di_matches_sequence_lengths(aligner, shared_datadir, tmp_path):
    _, inputs = structure_lib.resolve_structure_inputs([str(shared_datadir / "structures" / "inputs")])
    db = aligner.build_db(inputs, str(tmp_path / "indb"))
    sequences = {name: sequence for name, sequence, _ in aligner.iter_sequences(db)}
    threedi = dict(aligner.iter_3di(db))
    assert set(threedi) == set(sequences)
    for name, states in threedi.items():
        assert len(states) == len(sequences[name])


@needs_foldseek
def test_prebuilt_3di_database_has_no_coordinates(shared_datadir):
    """The committed FeSOD database was built from sequences, not structures."""
    db = structure_lib.StructureDB(prefix=str(shared_datadir / "foldseek" / "FeSOD"))
    assert not db.has_coordinates()


@needs_foldseek
def test_check_capabilities_rejects_coordinate_work_without_coordinates(aligner, shared_datadir):
    db = structure_lib.StructureDB(prefix=str(shared_datadir / "foldseek" / "FeSOD"))
    with pytest.raises(RuntimeError, match="requires C-alpha coordinates"):
        aligner.check_capabilities(1, [], [db])
    with pytest.raises(RuntimeError, match="C-alpha coordinates"):
        aligner.check_capabilities(2, ["tmscore"], [db])
    # metrics that need no coordinates are fine
    aligner.check_capabilities(2, ["prob"], [db])


@needs_foldseek
def test_check_capabilities_rejects_unknown_requests(aligner, shared_datadir):
    db = structure_lib.StructureDB(prefix=str(shared_datadir / "foldseek" / "FeSOD"))
    with pytest.raises(RuntimeError, match="does not support alignment type 7"):
        aligner.check_capabilities(7, [], [db])


@needs_foldseek
def test_search_reports_self_hit_over_full_length(aligner, shared_datadir, tmp_path):
    _, inputs = structure_lib.resolve_structure_inputs([str(shared_datadir / "structures" / "inputs")])
    _, refs = structure_lib.resolve_structure_inputs([str(shared_datadir / "structures" / "refs")])
    indb = aligner.build_db(inputs, str(tmp_path / "indb"))
    refdb = aligner.build_db(refs, str(tmp_path / "refdb"))
    work = tmp_path / "work"
    work.mkdir()

    hits = list(aligner.search(refdb, indb, evalue=0.001, max_seqs=1000, want_tseq=True,
                               alignment_type=2, metrics=[], work_dir=str(work)))
    by_target = {hit.target: hit for hit in hits}
    assert set(by_target) == {"1UBQ_A", "1NDD_A"}

    self_hit = by_target["1UBQ_A"]
    assert self_hit.fident == pytest.approx(1.0)
    assert (self_hit.tstart, self_hit.tend, self_hit.tlen) == (1, 76, 76)
    assert len(self_hit.tseq) == 76, "tseq must carry the full target sequence"


@needs_foldseek
def test_search_sorted_by_target(aligner, shared_datadir, tmp_path):
    _, inputs = structure_lib.resolve_structure_inputs([str(shared_datadir / "structures" / "inputs")])
    _, refs = structure_lib.resolve_structure_inputs([str(shared_datadir / "structures" / "refs")])
    indb = aligner.build_db(inputs, str(tmp_path / "indb"))
    refdb = aligner.build_db(refs, str(tmp_path / "refdb"))
    work = tmp_path / "work"
    work.mkdir()
    hits = list(aligner.search(refdb, indb, evalue=10, max_seqs=1000, want_tseq=False,
                               alignment_type=2, metrics=[], work_dir=str(work),
                               sort_by_target=True))
    targets = [hit.target for hit in hits]
    assert targets == sorted(targets)


@needs_foldseek
def test_search_with_metrics(aligner, shared_datadir, tmp_path):
    """lddt and rmsd need the alignment backtrace, which search only writes with -a."""
    _, inputs = structure_lib.resolve_structure_inputs([str(shared_datadir / "structures" / "inputs")])
    _, refs = structure_lib.resolve_structure_inputs([str(shared_datadir / "structures" / "refs")])
    indb = aligner.build_db(inputs, str(tmp_path / "indb"))
    refdb = aligner.build_db(refs, str(tmp_path / "refdb"))
    work = tmp_path / "work"
    work.mkdir()
    hits = list(aligner.search(refdb, indb, evalue=10, max_seqs=1000, want_tseq=False,
                               alignment_type=1,
                               metrics=["tmscore", "lddt", "rmsd", "prob"],
                               work_dir=str(work)))
    self_hit = next(hit for hit in hits if hit.target == "1UBQ_A")
    # foldseek's approximate TM-score can exceed 1.0, so this is not bounded above by 1
    assert self_hit.tmscore > 0.95
    assert self_hit.lddt == pytest.approx(1.0, abs=0.01)
    assert self_hit.rmsd < 0.1
    assert self_hit.prob > 0.9


@needs_foldseek
def test_prepared_databases_rejects_keep_db_for_database_input(aligner, shared_datadir, tmp_path):
    with pytest.raises(RuntimeError, match="already a database"):
        with structure_lib.prepared_databases(
            aligner,
            [str(shared_datadir / "foldseek" / "FeSOD")],
            [str(shared_datadir / "structures" / "refs")],
            keep_db=str(tmp_path / "kept"),
        ):
            pass


@needs_foldseek
def test_prepared_databases_labels(aligner, shared_datadir, tmp_path):
    with structure_lib.prepared_databases(
        aligner,
        [str(shared_datadir / "structures" / "inputs")],
        [str(shared_datadir / "structures" / "refs" / "UBQref.pdb.gz")],
        tmp_dir=str(tmp_path),
    ) as prepared:
        assert prepared.input_label == "inputs"
        assert prepared.reference_label == "UBQref"


def test_group_hits_by_target_carries_the_target_sequence():
    """The sequence travels with its group so a bounded heap bounds memory as well."""
    def hit(target, tseq):
        return structure_lib.StructureHit(
            query="R", target=target, qstart=1, qend=2, qlen=2, tstart=1, tend=2, tlen=2,
            evalue=1e-9, bits=1.0, fident=1.0, alnlen=2, tseq=tseq)

    hits = [hit("A", "MQIF"), hit("A", ""), hit("B", "GGGG")]
    groups = list(structure_lib.group_hits_by_target(hits, "db", evalue=1.0))
    assert [(target, sequence) for target, _, sequence in groups] == [("A", "MQIF"), ("B", "GGGG")]


@needs_foldseek
def test_iter_flat_db_streams_without_slurping(aligner, shared_datadir, tmp_path):
    """Entries are seeked to individually, so memory does not scale with database size."""
    _, inputs = structure_lib.resolve_structure_inputs([str(shared_datadir / "structures" / "inputs")])
    db = aligner.build_db(inputs, str(tmp_path / "indb"))
    entries = dict(structure_lib.FoldseekAligner._iter_flat_db(db.prefix))
    assert len(entries) == 9
    assert all(sequence and "\x00" not in sequence and "\n" not in sequence
               for sequence in entries.values())
