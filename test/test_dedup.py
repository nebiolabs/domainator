import json
import gc
import random
import weakref
from contextlib import nullcontext

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation

from domainator import domain_search, select_by_cds
from domainator.cli.manifest import describe_tool
from domainator.dedup import DEDUP_CHOICES, RecordDeduplicator, TopHits, check_dedup
from domainator.errors import InputError
from domainator.model import get_hit_best_score, set_dist_to_start, set_dist_to_end, set_focus_strand, set_hit_best_score
from domainator.seqio import parse_seqfiles, write_genbank


def record(name, sequence, molecule_type="DNA"):
    return SeqRecord(Seq(sequence), id=name, annotations={"molecule_type": molecule_type})


@pytest.mark.parametrize("dedup,expected", [
    ("none", [0, 1, 2, 3, 4, 5]),
    ("id", [0, 2, 3, 4, 5]),
    ("sequence", [0, 3, 4, 5]),
    ("sequence_rc", [0, 4]),
])
def test_policies_with_consistent_ids(dedup, expected):
    records = [
        record("a", "AACG"),
        record("a", "AACG"),
        record("b", "AACG"),
        record("c", "CGTT"),
        record("d", "TTGC"),
        record("e", "GCAA"),
    ]
    deduplicator = RecordDeduplicator(dedup)
    assert [i for i, rec in enumerate(records) if deduplicator.keep(rec)] == expected


@pytest.mark.parametrize("dedup", ["unknown", "", None, [], ["sequence"], ["id", "sequence_rc"]])
def test_invalid_policies(dedup):
    with pytest.raises(InputError):
        check_dedup(dedup)
    with pytest.raises(InputError):
        RecordDeduplicator(dedup)
    with pytest.raises(InputError):
        TopHits(1, dedup)


@pytest.mark.parametrize("dedup", ["sequence", "sequence_rc"])
def test_sequence_normalization(dedup):
    deduplicator = RecordDeduplicator(dedup)
    assert deduplicator.keep(record("a", "AACG"))
    assert not deduplicator.keep(record("b", "aacg*"))


def test_protein_rc_rejected():
    with pytest.raises(InputError, match="only supported for nucleotide output"):
        RecordDeduplicator("sequence_rc").keep(record("a", "MA", "protein"))


def test_dedup_before_padding():
    records = [record("a", "AACG"), record("b", "AACG")]
    for rec, start in zip(records, (0, 2)):
        set_dist_to_start(rec, start)
        set_dist_to_end(rec, 0)
        set_focus_strand(rec, 1)
    result = list(select_by_cds.deduplicate_and_pad(records, dedup="sequence", pad=True))
    assert [(rec.id, str(rec.seq)) for rec in result] == [("a", "AACG")]


@pytest.mark.parametrize("module,default", [(domain_search, "none"), (select_by_cds, "id")])
def test_parser_policies_and_config(module, default, tmp_path):
    args = ["-i", "input.gb", "-o", "output.gb"]
    if module is domain_search:
        args += ["-r", "query.hmm"]
    parser = module.build_parser()
    assert parser.parse_args(args).dedup == default
    assert parser.parse_args(args + ["--dedup", "sequence_rc"]).dedup == "sequence_rc"
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"dedup": "sequence_rc"}))
    assert parser.parse_args(args + ["--config", str(config)]).dedup == "sequence_rc"
    assert parser.parse_args(args + ["--config", str(config), "--dedup", "none"]).dedup == "none"


@pytest.mark.parametrize("module", [domain_search, select_by_cds])
def test_multiple_policies_rejected(module, tmp_path):
    args = ["-i", "input.gb", "-o", "output.gb"]
    if module is domain_search:
        args += ["-r", "query.hmm"]
    parser = module.build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(args + ["--dedup", "id", "sequence_rc"])
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"dedup": ["sequence_rc"]}))
    with pytest.raises(SystemExit):
        parser.parse_args(args + ["--config", str(config)])


@pytest.mark.parametrize("command,default", [("domain_search", "none"), ("select_by_cds", "id")])
def test_manifest_describes_scalar_policy(command, default):
    parameter = next(p for p in describe_tool(command)["parameters"] if p["name"] == "dedup")
    assert parameter["kind"] == "value"
    assert parameter["type"] == "string"
    assert parameter["choices"] == list(DEDUP_CHOICES)
    assert parameter["default"] == default
    assert not parameter.get("multiple", False)


@pytest.mark.parametrize("module", [domain_search, select_by_cds])
@pytest.mark.parametrize("flag", ["--deduplicate", "--no_deduplicate",
                                 "--deduplicate_by_sequence", "--both_strands"])
def test_old_flags_rejected(module, flag):
    args = ["-i", "input.gb", "-o", "output.gb"]
    if module is domain_search:
        args += ["-r", "query.hmm"]
    with pytest.raises(SystemExit):
        module.build_parser().parse_args(args + [flag])


@pytest.mark.parametrize("dedup,expected", [
    (None, [0, 2, 3, 4, 5]),
    ("none", [0, 1, 2, 3, 4, 5]),
    ("id", [0, 2, 3, 4, 5]),
    ("sequence", [0, 3, 4, 5]),
    ("sequence_rc", [0, 4]),
])
def test_select_by_cds_cli_policies(tmp_path, dedup, expected):
    records = [record("a", "AACG"), record("a", "AACG"),
               record("b", "AACG"), record("c", "CGTT"),
               record("d", "TTGC"), record("e", "GCAA")]
    for rec in records:
        rec.features.append(SeqFeature(SimpleLocation(0, len(rec), strand=1), type="CDS",
                                       qualifiers={"locus_tag": ["cds"], "translation": ["M"]}))
    input_path = tmp_path / "input.gb"
    output_path = tmp_path / "output.gb"
    write_genbank(records, input_path)
    args = ["-i", str(input_path), "-o", str(output_path), "--whole_contig"]
    if dedup is not None:
        args += ["--dedup", dedup]
    select_by_cds.main(args)
    result = list(parse_seqfiles([output_path]))
    assert [str(rec.seq) for rec in result] == [str(records[i].seq) for i in expected]


def full_sort_reference(records, max_hits, dedup):
    seen = set()
    result = []
    for _, rec in sorted(enumerate(records), key=lambda item: (-get_hit_best_score(item[1]), item[0])):
        sequence = str(rec.seq).upper().strip("*")
        if dedup == "sequence_rc":
            sequence = min(sequence, str(Seq(sequence).reverse_complement()))
        key = rec.id if dedup == "id" else sequence
        if dedup != "none" and key in seen:
            continue
        seen.add(key)
        result.append(rec)
        if len(result) == max_hits:
            break
    return result


@pytest.mark.parametrize("dedup", DEDUP_CHOICES)
@pytest.mark.parametrize("max_hits", [1, 7, 100])
def test_streaming_top_hits_match_full_sort(dedup, max_hits):
    rng = random.Random(42)
    sequences = {}
    for i in range(20):
        sequence = "".join(rng.choices("ACGT", k=12))
        sequences[f"{i}_forward"] = sequence
        sequences[f"{i}_alias"] = sequence
        sequences[f"{i}_reverse"] = str(Seq(sequence).reverse_complement())
    names = list(sequences)
    records = []
    winners = TopHits(max_hits, dedup)
    for i in range(800):
        name = rng.choice(names)
        rec = record(name, sequences[name])
        set_hit_best_score(rec, rng.randrange(20))
        records.append(rec)
        winners.add(rec)
        assert len(winners) <= max_hits
        assert len(winners._heap) <= 2 * max_hits
        if i < 20 or i % 101 == 0:
            expected = full_sort_reference(records, max_hits, dedup)
            assert [id(hit) for hit in winners.records()] == [id(hit) for hit in expected]
    expected = full_sort_reference(records, max_hits, dedup)
    assert [id(hit) for hit in winners.records()] == [id(hit) for hit in expected]


def test_top_hits_release_losers_and_bound_heap_metadata():
    winners = TopHits(3, "sequence")
    references = []
    for i in range(2000):
        group = i % 23 if i < 1000 else 0
        rec = record(str(group), format(group, "08b").translate(str.maketrans("01", "AC")))
        set_hit_best_score(rec, i)
        references.append(weakref.ref(rec))
        winners.add(rec)
        assert len(winners) <= 3
        assert len(winners._heap) <= 6
    del rec
    gc.collect()
    assert sum(ref() is not None for ref in references) == len(winners) == 3


@pytest.mark.parametrize("max_hits", [0, -1])
def test_top_hits_reject_invalid_limit(max_hits):
    with pytest.raises(InputError, match="at least 1"):
        TopHits(max_hits)


@pytest.mark.parametrize("dedup", DEDUP_CHOICES)
@pytest.mark.parametrize("batch_size", [1, 17, 1000])
@pytest.mark.parametrize("max_hits", [1, 5, 100])
def test_partitioned_top_hits_match_full_sort(monkeypatch, dedup, batch_size, max_hits):
    rng = random.Random(123)
    sequences = {}
    for i in range(20):
        sequence = "".join(rng.choices("ACGT", k=12))
        sequences[f"{i}_forward"] = sequence
        sequences[f"{i}_alias"] = sequence
        sequences[f"{i}_reverse"] = str(Seq(sequence).reverse_complement())
    names = list(sequences)
    records = []
    for score in rng.sample(range(10000), 300):
        name = rng.choice(names)
        rec = record(name, sequences[name])
        set_hit_best_score(rec, score)
        records.append(rec)

    monkeypatch.setattr(domain_search, "make_pool", lambda **kwargs: nullcontext())
    monkeypatch.setattr(domain_search.domainate, "read_references", lambda *args, **kwargs: {})
    monkeypatch.setattr(domain_search.domainate, "domainate", lambda records, **kwargs: records)
    monkeypatch.setattr(select_by_cds, "select_by_cds", lambda records, **kwargs: records)

    def run_tasks(pool, worker, tasks, **kwargs):
        assert worker.max_hits == max_hits
        batches = list(tasks)
        rng.shuffle(batches)
        for batch in batches:
            hits = worker(batch)
            assert len(hits) <= max_hits
            yield hits

    monkeypatch.setattr(domain_search, "imap_unordered_bounded", run_tasks)
    result = list(domain_search.domain_search(
        records, references=[], z=1000, cpu=2, batch_size=batch_size,
        max_hits=max_hits, dedup=dedup,
    ))
    expected = full_sort_reference(records, max_hits, dedup)
    assert [id(hit) for hit in result] == [id(hit) for hit in expected]
