import tempfile

from domainator import DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from domainator.domain_search import main as domain_search_main
from domainator.sort_contigs import main as sort_contigs_main
from domainator.seqio import write_genbank

from helpers import compare_seqfiles
import pytest


def test_sort_contigs_matches_domain_search_max_hits(shared_datadir):
    gb = shared_datadir / "pDONR201.gb"
    hmms = shared_datadir / "pdonr_peptides.fasta"

    with tempfile.TemporaryDirectory() as output_dir:
        unsorted_out = output_dir + "/domain_search_unsorted.gb"
        sorted_out = output_dir + "/domain_search_sorted.gb"
        resorted_out = output_dir + "/sort_contigs.gb"

        domain_search_main([
            "--input", str(gb),
            "-r", str(hmms),
            "--evalue", "0.1",
            "-o", str(unsorted_out),
            "--max_overlap", "1",
            "-Z", "1000",
            "--cpu", "1",
        ])

        domain_search_main([
            "--input", str(gb),
            "-r", str(hmms),
            "--evalue", "0.1",
            "-o", str(sorted_out),
            "--max_overlap", "1",
            "-Z", "1000",
            "--cpu", "1",
            "--max_hits", "2",
        ])

        sort_contigs_main([
            "--input", str(unsorted_out),
            "-o", str(resorted_out),
            "--max_hits", "2",
        ])

        compare_seqfiles(resorted_out, sorted_out, skip_qualifiers={"accession"})


def test_sort_contigs_rejects_domain_search_feature_without_score():
    record = SeqRecord(Seq("ATGC"), id="bad_score", name="bad_score", description="bad score")
    record.annotations["molecule_type"] = "DNA"
    record.features.append(
        SeqFeature(
            FeatureLocation(0, 4, strand=1),
            type=DOMAIN_SEARCH_BEST_HIT_NAME,
            qualifiers={
                "program": ["hmmsearch"],
                "database": ["test"],
                "evalue": ["1.0e-5"],
                "name": ["test_hit"],
            },
        )
    )

    with tempfile.TemporaryDirectory() as output_dir:
        input_path = output_dir + "/bad_score.gb"
        output_path = output_dir + "/sorted.gb"
        with open(input_path, "w") as handle:
            write_genbank((record,), handle)

        with pytest.raises(ValueError, match="without a valid score qualifier"):
            sort_contigs_main([
                "--input", str(input_path),
                "-o", str(output_path),
            ])


def make_scored_record(record_id, domainator_hits=None, search_hits=None):
    record = SeqRecord(Seq("ATGCATGC"), id=record_id, name=record_id, description=record_id)
    record.annotations["molecule_type"] = "DNA"

    if domainator_hits is not None:
        for index, (name, database, score) in enumerate(domainator_hits):
            record.features.append(
                SeqFeature(
                    FeatureLocation(index, index + 2, strand=1),
                    type=DOMAIN_FEATURE_NAME,
                    qualifiers={
                        "program": ["hmmsearch"],
                        "database": [database],
                        "description": [f"{name} description"],
                        "accession": [name],
                        "evalue": ["1.0e-5"],
                        "score": [str(score)],
                        "name": [name],
                        "cds_id": [str(index)],
                    },
                )
            )

    if search_hits is not None:
        for index, score in enumerate(search_hits):
            record.features.append(
                SeqFeature(
                    FeatureLocation(index + 2, index + 4, strand=1),
                    type=DOMAIN_SEARCH_BEST_HIT_NAME,
                    qualifiers={
                        "program": ["hmmsearch"],
                        "database": ["search_db"],
                        "description": [f"search_{index} description"],
                        "accession": [f"search_{index}"],
                        "evalue": ["1.0e-5"],
                        "score": [str(score)],
                        "name": [f"search_{index}"],
                        "cds_id": [str(index)],
                    },
                )
            )

    return record


def test_sort_contigs_domains_sorts_by_best_matching_domainator_score():
    records = [
        make_scored_record("rec1", domainator_hits=[("target", "db1", 40), ("other", "db1", 90)], search_hits=[5]),
        make_scored_record("rec2", domainator_hits=[("target", "db2", 70)], search_hits=[1]),
    ]

    with tempfile.TemporaryDirectory() as output_dir:
        input_path = output_dir + "/domainator.gb"
        output_path = output_dir + "/sorted.gb"
        with open(input_path, "w") as handle:
            write_genbank(records, handle)

        sort_contigs_main([
            "--input", str(input_path),
            "-o", str(output_path),
            "--domains", "target",
        ])

        out_records = list(SeqIO.parse(output_path, "genbank"))
        assert [record.id for record in out_records] == ["rec2", "rec1"]


def test_sort_contigs_domains_databases_filters_domainator_scores():
    records = [
        make_scored_record("rec1", domainator_hits=[("target", "db1", 80), ("target", "db2", 20)]),
        make_scored_record("rec2", domainator_hits=[("target", "db2", 70)]),
    ]

    with tempfile.TemporaryDirectory() as output_dir:
        input_path = output_dir + "/domainator.gb"
        output_path = output_dir + "/sorted.gb"
        with open(input_path, "w") as handle:
            write_genbank(records, handle)

        sort_contigs_main([
            "--input", str(input_path),
            "-o", str(output_path),
            "--domains", "target",
            "--databases", "db2",
        ])

        out_records = list(SeqIO.parse(output_path, "genbank"))
        assert [record.id for record in out_records] == ["rec2", "rec1"]


def test_sort_contigs_domains_raises_when_no_matching_score():
    record = make_scored_record("rec1", domainator_hits=[("other", "db1", 80)], search_hits=[100])

    with tempfile.TemporaryDirectory() as output_dir:
        input_path = output_dir + "/domainator.gb"
        output_path = output_dir + "/sorted.gb"
        with open(input_path, "w") as handle:
            write_genbank((record,), handle)

        with pytest.raises(ValueError, match="matching the requested domains"):
            sort_contigs_main([
                "--input", str(input_path),
                "-o", str(output_path),
                "--domains", "target",
            ])

def test_sort_contigs_by_system_score():
    import pytest
    from domainator import SYSTEM_SEARCH_BEST_HIT_NAME
    from domainator.find_systems import find_systems
    from domainator.sort_contigs import sort_contigs
    from system_records import contig
    records = [contig(["GajA:20", "GajB:20"], record_id="weak"), contig(["GajA:90", "GajB:90"], record_id="strong"),
               contig(["GajA:50", "GajB:50"], record_id="middle")]
    records = list(find_systems(records, domain_expr="GajA & GajB", window_cds=2, system_name="Gabija"))
    assert [r.id for r in sort_contigs(records, system_type="domain_system")] == ["strong", "middle", "weak"]
    assert [r.id for r in sort_contigs(records, system_type="all", max_hits=1)] == ["strong"]
    with pytest.raises(ValueError, match="does not contain any system_search systems"):
        list(sort_contigs(records, system_type="system_search"))
    for record in records:
        for feature in record.features:
            if feature.type == "Domain_System":
                feature.type = SYSTEM_SEARCH_BEST_HIT_NAME
    assert [r.id for r in sort_contigs(records, system_type="system_search")] == ["strong", "middle", "weak"]
