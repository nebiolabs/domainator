import tempfile

from domainator import DOMAIN_SEARCH_BEST_HIT_NAME
from domainator.Bio.Seq import Seq
from domainator.Bio.SeqFeature import SeqFeature, FeatureLocation
from domainator.Bio.SeqRecord import SeqRecord
from domainator.domain_search import main as domain_search_main
from domainator.sort_contigs import main as sort_contigs_main
from domainator.utils import write_genbank

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