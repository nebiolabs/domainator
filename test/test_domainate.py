import os
import domainator.domainate as domainate_module
from domainator.domainate import main, filter_by_overlap, SearchResult
import tempfile
from glob import glob
from domainator.Bio.Seq import Seq
from domainator.Bio.SeqRecord import SeqRecord
from domainator.Bio import SeqIO
from domainator import DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME
from domainator.utils import DomainatorCDS, count_peptides_in_record
from domainator.domainate import read_references
import pyhmmer
import pytest
from helpers import compare_seqfiles, compare_seqrecords, gzip_file, bgzip_file

@pytest.mark.parametrize("Z,expected_string",
[(0, "CcdB (CcdB protein, 1.3e-33, 103.1)"),
(10000000, "CcdB (CcdB protein, 4.4e-27, 103.1)"),
(1, "CcdB (CcdB protein, 4.4e-34, 103.1)"),

])
def test_domainator_one_file(Z, expected_string, shared_datadir):
    gb = shared_datadir / "pDONR201.gb"
    hmms = shared_datadir / "CcdB.hmm"

    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir="test_out"
        out = output_dir + f"/out{Z}.gb"
        args = ['--input', str(gb) , "-r", str(hmms), "--evalue", str(0.1), "-o", str(out), "--max_domains", str(1), "--max_overlap", str(1)]
        args += ["-Z", str(Z)]
        main(args)
        assert os.path.isfile(out)
        
        #assert len(glob(output_dir+"/*.gb")) == 1
        new_file = list(SeqIO.parse(out, "genbank"))
        domainator_features = [x for x in new_file[0].features if x.type == DOMAIN_FEATURE_NAME]
        assert len(domainator_features) == 1
        assert domainator_features[0].qualifiers["database"] == ["CcdB"]
        
        CDS_features = {x.qualifiers["cds_id"][0]:x for x in new_file[0].features if x.type == "CDS"}
        assert len(CDS_features) == 3
        assert CDS_features["1264_-1_959"].qualifiers["domainator_CcdB"][0] == expected_string

def test_domainator_multi(shared_datadir):
    gbs = [shared_datadir / "pDONR201.gb", shared_datadir / "Polymorphism_feature.gb",
           shared_datadir / "Staph_phages.gb"]
    hmms = shared_datadir / "CcdB.hmm"

    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/out1.gb"
        args = ['--input'] + [str(x) for x in gbs] + [ "-r", str(hmms), "--evalue", str(0.1), "-o", str(out), "--max_domains", str(1), "--max_overlap", str(1)]

        main(args)
        assert len(glob(output_dir+"/*.gb")) == 1

def test_domainator_no_cdss(shared_datadir):
    gb = shared_datadir / "pDONR201_empty.gb"
    hmms = shared_datadir / "CcdB.hmm"


    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/out1.gb"
        args = ['--input', str(gb),  "-r", str(hmms), "--evalue", str(0.1), "-o", str(out), "--max_domains", str(1), "--max_overlap", str(1)]

        main(args)
        assert len(glob(output_dir+"/*.gb")) == 1

def test_domainator_pseudo(shared_datadir):
    gb = shared_datadir / "pDONR201_pseudo.gb"
    hmms = shared_datadir / "CcdB.hmm"


    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/out1.gb"
        args = ['--input', str(gb),  "-r", str(hmms), "--evalue", str(0.1), "-o", str(out), "--hits_only", "--max_domains", str(1), "--max_overlap", str(1)]

        main(args)
        assert len(glob(output_dir+"/*.gb")) == 1
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 0


def test_domainator_multi_hmm(shared_datadir):
    query_seqs = shared_datadir / "pDONR201_multi_genemark.gb"
    hmms = shared_datadir / "pdonr_hmms.hmm"
    ref_file = shared_datadir / "pDONR201_multi_genemark_domainator.gb"

    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/out1.gb"
        args = ['--input'] + [str(query_seqs)] + [ "-r", str(hmms), "-o", str(out), "-e", "10", "-Z", "0"]

        main(args)

        # The reference keeps GeneMark's /partial="10" and "01" codes, so it also tests reading them. Domainator converts them
        # to '<' / '>' markers when reading, so compare against the reference as Domainator reads it.
        # The reference's cds_id for the /partial="10" CDS was generated from its old location (2..106); new runs generate it from <1..106.
        ref_records = list(SeqIO.parse(str(ref_file), "genbank"))
        for record in ref_records:
            _utils.normalize_partial_codes(record)
            for feature in record.features:
                if feature.qualifiers.get("cds_id") == ["2_1_106"]:
                    feature.qualifiers["cds_id"] = ["1_1_106"]
        normalized_ref = output_dir + "/ref.gb"
        SeqIO.write(ref_records, normalized_ref, "genbank")
        compare_seqfiles(out, normalized_ref, skip_qualifiers={"identity", "accession"})

def test_domainator_multi_hmm_2(shared_datadir):
    query_seqs = shared_datadir / "pDONR201_multi_genemark.gb"
    hmms = [str(shared_datadir / "pdonr_hmms_1.hmm"), str(shared_datadir / "pdonr_hmms_2.hmm")]
    ref_file = shared_datadir / "pDONR201_multi_genemark_domainator_multi_hmm.gb"

    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/out1.gb"
        args = ['--input'] + [str(query_seqs)] + [ "-r"] + hmms + ["-o", str(out), "-e", "10", "-Z", "0"]

        main(args)

        compare_seqfiles(out, ref_file, skip_qualifiers={"identity", "accession"})


def test_domainator_multi_hmm_3(shared_datadir):
    query_seqs = shared_datadir / "pDONR201_multi_genemark.gb"
    hmms = [str(shared_datadir / "pdonr_hmms_1.hmm"), str(shared_datadir / "pdonr_hmms_2.hmm")]
    ref_file = shared_datadir / "pDONR201_multi_genemark_domainator_multi_hmm.gb"

    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/out3.gb"
        args = ['--input'] + [str(query_seqs)] + [ "-r"] + hmms + ["-o", str(out), "-e", "10", "-Z", "0", "--max_overlap", "0.6"]

        main(args)

        new_file = list(SeqIO.parse(out, "genbank"))
        domainator_features = [x for x in new_file[0].features if x.type == DOMAIN_FEATURE_NAME]
        assert len(domainator_features) == 4

def test_domainator_multi_hmm_4(shared_datadir):
    query_seqs = shared_datadir / "pDONR201_multi_genemark.gb"
    hmms = [str(shared_datadir / "pdonr_hmms_1.hmm"), str(shared_datadir / "pdonr_hmms_2.hmm")]

    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/out4.gb"
        args = ['--input'] + [str(query_seqs)] + [ "-r"] + hmms + ["-o", str(out), "-e", "10", "-Z", "0", "--max_overlap", "0.6", "--overlap_by_db"]
        main(args)

        new_file = list(SeqIO.parse(out, "genbank"))
        domainator_features = [x for x in new_file[0].features if x.type == DOMAIN_FEATURE_NAME]
        assert len(domainator_features) == 7
 

def test_domain_overlap():
    hits = [SearchResult('0', '1', 'PF01', 2, 3, 1, 10,"db", 80.0, 1, 10, 10), #SearchResult(name,desc,evalue,score,start,end,database)
            SearchResult('0', '1', 'PF01', 2, 3, 9, 19,"db", 80.0, 1, 10, 10),
            SearchResult('0', '1', 'PF01', 2, 3, 17, 24,"db", 80.0, 1, 10, 10),
            SearchResult('0', '1', 'PF01', 2, 3, 22, 25,"db", 80.0, 1, 10, 10),
            ]

    non_overlap_hits = [SearchResult('0', '1', 'PF01', 2, 3, 1, 10,"db", 80.0, 1, 10, 10),
                        SearchResult('0', '1', 'PF01', 2, 3, 9, 19,"db", 80.0, 1, 10, 10),
                        SearchResult('0', '1', 'PF01', 2, 3, 22, 25,"db", 80.0, 1, 10, 10)]
    actual = filter_by_overlap(hits, .2)
    assert actual == non_overlap_hits


def test_annotate_peptides(shared_datadir):

    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/out_fa"
        main(['-i', str(shared_datadir/"FeSOD_20.gb"), '--references', str(shared_datadir/"FeSOD_pfam.hmm"), "--output", out + ".gb"])
        assert os.path.isfile(out + ".gb")

#(genbanks, hmm, evalue, output, max_overlap, cpu, max_domains, format="genbank", gff=False):
def test_annotate_peptides_fasta(shared_datadir):

    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/out_fa"

        main(['--input', str(shared_datadir/"FeSOD_20.fasta"), '--references', str(shared_datadir/"FeSOD_pfam.hmm"), "--output", out + ".gb"])
        assert os.path.isfile(out + ".gb")


def test_domainator_phmmer(shared_datadir):
    gbs = [shared_datadir / "pDONR201.gb"]
    references = shared_datadir / "pdonr_peptides.fasta"

    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/out1.gb"
        args = ['--input'] + [str(x) for x in gbs] + [ "-r", str(references), "--evalue", str(0.1), "-o", str(out), "--max_domains", str(1), "--max_overlap", str(1)]

        main(args)
        new_file = list(SeqIO.parse(out, "genbank"))
        domainator_features = [x for x in new_file[0].features if x.type == DOMAIN_FEATURE_NAME]
        assert len(domainator_features) == 3


def test_domainator_evalue_cutoff_filters_hits(shared_datadir):
    gb = shared_datadir / "pDONR201.gb"
    hmms = shared_datadir / "CcdB.hmm"

    with tempfile.TemporaryDirectory() as output_dir:
        permissive_out = output_dir + "/permissive.gb"
        strict_out = output_dir + "/strict.gb"

        main([
            "--input", str(gb),
            "-r", str(hmms),
            "--evalue", "0.1",
            "-o", permissive_out,
            "--max_domains", "1",
            "--max_overlap", "1",
            "-Z", "0",
        ])
        main([
            "--input", str(gb),
            "-r", str(hmms),
            "--evalue", "1e-40",
            "-o", strict_out,
            "--max_domains", "1",
            "--max_overlap", "1",
            "-Z", "0",
        ])

        permissive_record = list(SeqIO.parse(permissive_out, "genbank"))[0]
        strict_record = list(SeqIO.parse(strict_out, "genbank"))[0]

        assert len([f for f in permissive_record.features if f.type == DOMAIN_FEATURE_NAME]) == 1
        assert len([f for f in strict_record.features if f.type == DOMAIN_FEATURE_NAME]) == 0


def test_domainator_overlap_prefers_higher_scoring_smaller_hit():
    hits = [
        SearchResult("small_high", "desc", "", 1e-30, 100.0, 10, 40, "db", 0.0, 1, 10, 30),
        SearchResult("large_low", "desc", "", 1e-20, 90.0, 5, 60, "db", 0.0, 1, 10, 55),
    ]

    filtered_hits = filter_by_overlap(hits, 0.5)

    assert [hit.name for hit in filtered_hits] == ["small_high"]


def test_read_references_nucleotide_fasta_routes_to_nhmmer(shared_datadir):
    references = [str(shared_datadir / "simple_dna_queries.fna")]

    parsed = read_references(references, foldseek=None)

    assert "nhmmer" in parsed
    assert "phmmer" not in parsed
    assert "simple_dna_queries" in parsed["nhmmer"]


def test_read_references_foldseek_paths_are_registered():
    parsed = read_references(None, foldseek=["/tmp/example_db"])

    assert parsed["foldseek"] == {"example_db": "/tmp/example_db"}


def test_domainator_rejects_protein_input_with_nucleotide_references(shared_datadir):
    protein_input = shared_datadir / "simple_genpept.gb"
    nucleotide_reference = shared_datadir / "simple_dna_queries.fna"

    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.gb"
        with pytest.raises(ValueError, match="nucleotide reference"):
            main([
                "--input", str(protein_input),
                "-r", str(nucleotide_reference),
                "--output", str(out),
            ])


def test_domainator_nucleotide_fasta_query_uses_nhmmer(shared_datadir):
    nucleotide_input = shared_datadir / "simple_dna_target.fna"
    nucleotide_reference = shared_datadir / "simple_dna_queries.fna"

    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.gb"
        main([
            "--input", str(nucleotide_input),
            "--fasta_type", "nucleotide",
            "-r", str(nucleotide_reference),
            "--output", str(out),
            "--evalue", "0.1",
        ])

        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1

        domainator_features = [x for x in new_file[0].features if x.type == DOMAIN_FEATURE_NAME]
        assert len(domainator_features) == 1
        assert domainator_features[0].qualifiers["program"] == ["nhmmer"]
        assert domainator_features[0].qualifiers["cds_id"] == ["."]



def test_domainator_nucleotide_hmm_query_uses_nhmmer(shared_datadir):
    nucleotide_input = shared_datadir / "simple_dna_target.fna"

    builder = pyhmmer.plan7.Builder(pyhmmer.easel.Alphabet.dna())
    background = pyhmmer.plan7.Background(pyhmmer.easel.Alphabet.dna())
    sequence = pyhmmer.easel.TextSequence(
        name=b"dna_hmm_query",
        sequence="TTGACCGATGCTAGTCGATCGTAGCTAGGCTAACCGTTAGCGATCGTACGATCGATGCTAGT",
    ).digitize(pyhmmer.easel.Alphabet.dna())
    hmm, _, _ = builder.build(sequence, background)

    with tempfile.TemporaryDirectory() as output_dir:
        hmm_path = os.path.join(output_dir, "dna_query.hmm")
        with open(hmm_path, "wb") as handle:
            hmm.write(handle)

        out = output_dir + "/out.gb"
        main([
            "--input", str(nucleotide_input),
            "--fasta_type", "nucleotide",
            "-r", str(hmm_path),
            "--output", str(out),
            "--evalue", "0.1",
        ])

        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        domainator_features = [x for x in new_file[0].features if x.type == DOMAIN_FEATURE_NAME]
        assert len(domainator_features) == 1
        assert domainator_features[0].qualifiers["program"] == ["nhmmer"]


def test_domainator_infernal_cm_query(shared_datadir):
    nucleotide_input = shared_datadir / "pANT_R100.fa"
    cm_reference = shared_datadir / "RF00042.cm"

    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.gb"
        main([
            "--input", str(nucleotide_input),
            "--fasta_type", "nucleotide",
            "-r", str(cm_reference),
            "--output", str(out),
            "--hits_only",
            "--evalue", "0.1",
        ])

        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) > 0
        domainator_features = [x for x in new_file[0].features if x.type == DOMAIN_FEATURE_NAME]
        assert len(domainator_features) > 0
        assert domainator_features[0].qualifiers["program"] == ["infernal"]


def _assert_mixed_contig_and_cds_hits_are_both_annotated(shared_datadir, monkeypatch, cds_hit_name, cds_hit_desc, cds_hit_database, cds_hit_program):
    contig = next(SeqIO.parse(shared_datadir / "pDONR201.gb", "genbank"))
    domainate_module.clean_rec(contig)
    cds_index = next(index for index, feature in enumerate(contig.features) if feature.type == "CDS")
    cds_id = contig.features[cds_index].qualifiers["cds_id"][0]

    def fake_run_search(*args, **kwargs):
        return {
            0: {
                -1: [
                    SearchResult(
                        "RF00042",
                        "cm hit",
                        "RF00042",
                        1e-20,
                        42.0,
                        10,
                        40,
                        "rfam",
                        0.0,
                        1,
                        30,
                        120,
                        "infernal",
                    )
                ],
                cds_index: [
                    SearchResult(
                        cds_hit_name,
                        cds_hit_desc,
                        "ACC",
                        1e-30,
                        55.0,
                        2,
                        20,
                        cds_hit_database,
                        80.0,
                        1,
                        18,
                        90,
                        cds_hit_program,
                    )
                ],
            }
        }

    monkeypatch.setattr(domainate_module, "run_search", fake_run_search)

    domainate_module.domainator_inner(
        [contig],
        [],
        [],
        [],
        [],
        {},
        0.1,
        1,
        None,
        False,
        False,
        10,
        1,
        False,
    )

    domainator_features = [feature for feature in contig.features if feature.type == DOMAIN_FEATURE_NAME]
    infernal_features = [
        feature
        for feature in domainator_features
        if feature.qualifiers["program"] == ["infernal"] and feature.qualifiers["cds_id"] == ["."]
    ]
    hmmsearch_features = [
        feature
        for feature in domainator_features
        if feature.qualifiers["program"] == [cds_hit_program] and feature.qualifiers["cds_id"] == [cds_id]
    ]

    assert len(infernal_features) == 1
    assert len(hmmsearch_features) == 1


def test_domainator_mixed_cm_and_protein_hmm_hits_are_both_annotated(shared_datadir, monkeypatch):
    _assert_mixed_contig_and_cds_hits_are_both_annotated(
        shared_datadir,
        monkeypatch,
        "CcdB",
        "protein hmm hit",
        "hmms",
        "hmmsearch",
    )


def test_domainator_mixed_cm_and_protein_fasta_hits_are_both_annotated(shared_datadir, monkeypatch):
    _assert_mixed_contig_and_cds_hits_are_both_annotated(
        shared_datadir,
        monkeypatch,
        "pDONR201_2",
        "protein fasta hit",
        "peptides",
        "phmmer",
    )


def test_domainator_nucleotide_multi_hit_annotations(shared_datadir):
    """Multiple non-overlapping nucleotide hits should each get a DOMAIN_FEATURE_NAME feature."""
    nucleotide_input = shared_datadir / "multi_hit_dna_target.fna"
    nucleotide_reference = shared_datadir / "multi_hit_dna_queries.fna"

    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + "/out.gb"
        main([
            "--input", str(nucleotide_input),
            "--fasta_type", "nucleotide",
            "-r", str(nucleotide_reference),
            "--output", str(out),
            "--evalue", "0.1",
        ])

        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        domainator_features = [x for x in new_file[0].features if x.type == DOMAIN_FEATURE_NAME]
        assert len(domainator_features) >= 2
        query_names = {f.qualifiers["name"][0] for f in domainator_features}
        assert "dna_query_A" in query_names
        assert "dna_query_B" in query_names


def test_domainator_nucleotide_query_spanning_circular_origin(shared_datadir):
    circular_sequence = "GCTAACCGTTAGCGATCGTACGATCGATGCTAGTCCGATTAACCGGTTAGGCTTACCGATGG"
    query_sequence = circular_sequence[-18:] + circular_sequence[:18]

    with tempfile.TemporaryDirectory() as output_dir:
        gb_path = os.path.join(output_dir, "circular.gb")
        query_path = os.path.join(output_dir, "origin_query.fna")
        out = os.path.join(output_dir, "out.gb")

        record = SeqRecord(Seq(circular_sequence), id="circular_test", name="circular_test", description="circular test")
        record.annotations["molecule_type"] = "DNA"
        record.annotations["topology"] = "circular"
        with open(gb_path, "w") as handle:
            SeqIO.write([record], handle, "genbank")
        with open(query_path, "w") as handle:
            handle.write(">origin_query\n")
            handle.write(query_sequence + "\n")

        main([
            "--input", gb_path,
            "-r", query_path,
            "--output", out,
            "--hits_only",
            "--evalue", "0.1",
        ])

        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        domainator_features = [x for x in new_file[0].features if x.type == DOMAIN_FEATURE_NAME]
        assert len(domainator_features) == 1
        assert len(domainator_features[0].location.parts) == 2


def test_domainator_min_evalue(shared_datadir):
    gbs = [shared_datadir / "pDONR201.gb"]
    references = shared_datadir / "pdonr_peptides.fasta"

    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/out1.gb"
        args = ['--input'] + [str(x) for x in gbs] + [ "-r", str(references), "--evalue", str(0.1), "-o", str(out), "--max_domains", str(1), "--max_overlap", str(1), "--min_evalue", "2e-160", "-Z", "0"]

        main(args)
        new_file = list(SeqIO.parse(out, "genbank"))
        domainator_features = [x for x in new_file[0].features if x.type == DOMAIN_FEATURE_NAME]
        assert len(domainator_features) == 1
        assert domainator_features[0].qualifiers["name"][0] == "pDONR201_2"



@pytest.mark.parametrize("files,offset,read_count,expected_record_ct,rec0_name",
[(["pDONR201.gb"],0,1,1,"pDONR201"),
(["pDONR201.gb"],0,10,1,"pDONR201"),
(["pDONR201.gb"],0,0,0,""),
(["pDONR201_multi_genemark.gb","pDONR201.gb"],16382,10,2, "pDONR201_3"), #seeks past the end of pDONR201.gb
])
def test_domainator_seek(files,offset,read_count,expected_record_ct,rec0_name,shared_datadir):
    hmms = shared_datadir / "CcdB.hmm"


    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/out1.gb"
        args = ['--input'] + [str(shared_datadir / x) for x in files] + [ "-r", str(hmms), "--evalue", str(0.1), "-o", str(out), "--max_domains", str(1), "--max_overlap", str(1)]
        args += ['--offset', str(offset), "--recs_to_read", str(read_count)]
        main(args)
        recs = list(SeqIO.parse(out, "genbank"))
        assert len(recs) == expected_record_ct
        if len(recs) > 0:
            assert recs[0].name == rec0_name

def test_domainator_origin_spanning_1(shared_datadir):
    gbs = [shared_datadir / "bacillus_phage_SPR.gb"]
    references = shared_datadir / "SPR.hmm"

    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/out1.gb"
        args = ['--input'] + [str(x) for x in gbs] + [ "-r", str(references), "-Z", "1000", "--evalue", str(0.1), "-o", str(out), "--max_domains", str(1), "--max_overlap", str(1)]

        main(args)
        compare_seqfiles(out, shared_datadir / "bacillus_phage_SPR_with_annotations.gb", skip_qualifiers={"identity", "accession"})


def test_domainator_intron_1(shared_datadir):
    gbs = [shared_datadir / "saccharomyces_extraction.gb"]
    references = shared_datadir / "saccharomyces_defense_finder.hmm"


    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/out1.gb"
        args = ['--input'] + [str(x) for x in gbs] + [ "-r", str(references), "-Z", "1000", "--evalue", str(0.1), "-o", str(out), "--max_domains", str(1), "--max_overlap", str(1)]

        main(args)
        recs = list(SeqIO.parse(out, "genbank"))
        assert len(recs) == 1
        assert len(recs[0].features) == 23

def test_domainator_intron_2(shared_datadir):
    gbs = [shared_datadir / "saccharomyces_extraction_circular.gb"]
    references = shared_datadir / "saccharomyces_defense_finder.hmm"

    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/out1.gb"
        args = ['--input'] + [str(x) for x in gbs] + [ "-r", str(references), "-Z", "1000", "--evalue", str(0.1), "-o", str(out), "--max_domains", str(1), "--max_overlap", str(1)]

        main(args)
        recs = list(SeqIO.parse(out, "genbank"))
        assert len(recs) == 1
        assert len(recs[0].features) == 23

def test_domainator_gene_annotate_1(shared_datadir):
    hmms = str(shared_datadir / "pdonr_hmms.hmm")

    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out1 = output_dir + "/out1.gb"
        query_seqs = shared_datadir / "pDONR201_no_CDSs.gb"
        args = ['--input'] + [str(query_seqs)] + [ "-r"] + [hmms] + ["-o", str(out1)] + ["--gene_call", "all", "-Z", "1000"]
        main(args) # This should run gene calling and domain prediction on the unannotated file

        out2 = output_dir + "/out2.gb"
        query_seqs = shared_datadir / "pDONR201_domainator_circular.gb"
        args = ['--input'] + [str(query_seqs)] + [ "-r"] +[hmms] + ["-o", str(out2)] + ["--gene_call", "all", "-Z", "1000"]
        main(args) # This should re-run gene calling and domain prediction on the unannotated file

        compare_seqfiles(out1, out2, skip_attrs={"id", "description", "name"})

def test_domainator_gene_annotate_2(shared_datadir):
    hmms = str(shared_datadir / "pdonr_hmms.hmm")

    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out1 = output_dir + "/out1.gb"
        query_seqs = shared_datadir / "pDONR201_no_CDSs.gb"
        args = ['--input'] + [str(query_seqs)] + [ "-r"] + [hmms] + ["-o", str(out1)] + ["--gene_call", "unannotated", "-Z", "1000"]
        main(args) # This should run gene calling and domain prediction on the unannotated file

        out2 = output_dir + "/out2.gb"
        query_seqs = shared_datadir / "pDONR201_partly_CDSs.gb"
        args = ['--input'] + [str(query_seqs)] + [ "-r"] +[hmms] + ["-o", str(out2)] + ["--gene_call", "unannotated", "-Z", "1000"]
        main(args) # This should re-run gene calling and domain prediction on the unannotated file

        out1 = list(SeqIO.parse(out1, "genbank"))
        out2 = list(SeqIO.parse(out2, "genbank"))
        compare_seqrecords(out1[0], out2[0], skip_attrs={"id", "description", "name"})
        out1_CDSs = DomainatorCDS.list_from_contig(out1[0])
        out2_CDSs = DomainatorCDS.list_from_contig(out2[1])
        assert len(out1_CDSs) == len(out2_CDSs) == 3
        assert 'gene_id' in out1_CDSs[0].feature.qualifiers
        assert 'gene_id' not in out2_CDSs[0].feature.qualifiers

def test_domainator_gene_annotate_fasta_1(shared_datadir):
    hmms = str(shared_datadir / "pdonr_hmms.hmm")

    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out1 = output_dir + "/out1.gb"
        query_seqs = shared_datadir / "pDONR201.fasta"
        args = ['--input'] + [str(query_seqs)] + [ "-r"] + [hmms] + ["-o", str(out1)] + ["--gene_call", "all", "-Z", "1000", "--fasta_type", "nucleotide"]
        main(args) # This should run gene calling and domain prediction on the unannotated file

        f_txt = open(out1).read()
        assert "/gene_id=\"pDONR201_3\"" in f_txt
        assert "/score=\"90.2\"" in f_txt

def test_domainate_taxonomy_1(shared_datadir):
    input = shared_datadir / "swissprot_CuSOD_subset.fasta"
    hmms = shared_datadir / "swissprot_CuSOD_subset.fasta"
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + f"/out.gb"
        args = ['--input', str(input) , "-r", str(hmms), "--evalue", str(0.1), "-o", str(out), "--max_overlap", str(1), "-Z", "1000", "--ncbi_taxonomy_path", str(shared_datadir / "taxdmp"), "--include_taxids", "2", "--exclude_taxids", "1224"]
        main(args)
        assert os.path.isfile(out)
        
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        assert count_peptides_in_record(new_file[0]) == 1
        assert new_file[0].id == "sp|O31851|YOJM_BACSU"

def test_domainate_taxonomy_expr(shared_datadir):
    # "2 & ~1224" is equivalent to --include_taxids 2 --exclude_taxids 1224.
    input = shared_datadir / "swissprot_CuSOD_subset.fasta"
    hmms = shared_datadir / "swissprot_CuSOD_subset.fasta"
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + f"/out.gb"
        args = ['--input', str(input), "-r", str(hmms), "--evalue", str(0.1), "-o", str(out), "--max_overlap", str(1), "-Z", "1000", "--ncbi_taxonomy_path", str(shared_datadir / "taxdmp"), "--taxonomy_expr", "2 & ~1224"]
        main(args)
        assert os.path.isfile(out)

        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        assert count_peptides_in_record(new_file[0]) == 1
        assert new_file[0].id == "sp|O31851|YOJM_BACSU"

def test_domainate_mixed_query_1(shared_datadir):
    input = shared_datadir / "pDONR201.gb"
    hmms = shared_datadir / "pdonr_hmms.hmm"
    fasta = shared_datadir / "pdonr_peptides.fasta"
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + f"/out.gb"
        args = ['--input', str(input) , "-r", str(hmms), str(fasta), "--evalue", str(0.1), "-o", str(out), "--max_overlap", str(1), "-Z", "1000"]
        main(args)
        assert os.path.isfile(out)
        
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        assert count_peptides_in_record(new_file[0]) == 3
        cdss = DomainatorCDS.list_from_contig(new_file[0])
        assert len(cdss) == 3
        assert len(cdss[0].domain_features) == 2
        assert len(cdss[1].domain_features) == 3
        assert len(cdss[2].domain_features) == 3


def _database_qualifiers(path):
    dbs = set()
    for rec in SeqIO.parse(path, "genbank"):
        for feat in rec.features:
            if feat.type == DOMAIN_FEATURE_NAME:
                dbs.update(feat.qualifiers.get("database", []))
    return dbs


@pytest.mark.parametrize("compressor,suffix", [
    (gzip_file, ".hmm.gz"),
    (bgzip_file, ".hmm.bgz"),   # raised "format not recognized by HMMER" before
])
def test_annotate_with_compressed_references(shared_datadir, tmp_path, compressor, suffix):
    """A compressed reference must give byte-for-byte the same result as a plain one.

    It previously did not: the database qualifier was derived with a single-suffix
    stem, so a .gz reference was labeled "FeSOD_pfam.hmm" instead of "FeSOD_pfam".
    """
    plain_ref = shared_datadir / "FeSOD_pfam.hmm"
    compressed_ref = compressor(plain_ref, tmp_path / ("FeSOD_pfam" + suffix))

    plain_out = str(tmp_path / "plain.gb")
    compressed_out = str(tmp_path / "compressed.gb")
    main(['-i', str(shared_datadir / "FeSOD_20.gb"), '--references', str(plain_ref), '-o', plain_out])
    main(['-i', str(shared_datadir / "FeSOD_20.gb"), '--references', str(compressed_ref), '-o', compressed_out])

    assert _database_qualifiers(plain_out) == {"FeSOD_pfam"}
    assert _database_qualifiers(compressed_out) == {"FeSOD_pfam"}
    compare_seqfiles(plain_out, compressed_out)


def test_annotate_ignores_stale_pressed_sidecars(shared_datadir, tmp_path):
    """hmmpress sidecars left next to an edited .hmm must not be used."""
    from domainator import utils
    ref = tmp_path / "ref.hmm"
    ref.write_bytes((shared_datadir / "FeSOD_pfam.hmm").read_bytes())
    pyhmmer.hmmer.hmmpress(list(utils.open_hmm_file(ref)), str(ref))
    # Replace the profiles; the sidecars still describe the FeSOD ones.
    ref.write_bytes((shared_datadir / "pdonr_hmms.hmm").read_bytes())

    out = str(tmp_path / "out.gb")
    main(['-i', str(shared_datadir / "FeSOD_20.gb"), '--references', str(ref), '-o', out])
    # The stale sidecars would have annotated Sod_Fe_C / Sod_Fe_N.
    names = {q for rec in SeqIO.parse(out, "genbank") for f in rec.features
             if f.type == DOMAIN_FEATURE_NAME for q in f.qualifiers.get("name", [])}
    assert not (names & {"Sod_Fe_C", "Sod_Fe_N"})


# --- partial CDSs: '<' / '>' markers, codon_start, and the --partial filter ---

from domainator import utils as _utils

_PARTIAL_CDS_LOCUS_TAG = "HOV79_30125" # complement(803..>2688), /codon_start=3, translation of the on-contig part


def _write_partial_cds_hmm(shared_datadir, path):
    """Builds a single-sequence HMM from the translation of the partial CDS in JABFVH010000506_extraction.gb."""
    record = next(SeqIO.parse(str(shared_datadir / "JABFVH010000506_extraction.gb"), "genbank"))
    cds = [f for f in record.features if f.type == "CDS" and f.qualifiers.get("locus_tag") == [_PARTIAL_CDS_LOCUS_TAG]][0]
    alphabet = pyhmmer.easel.Alphabet.amino()
    sequence = pyhmmer.easel.TextSequence(name=b"partial_cds", sequence=cds.qualifiers["translation"][0]).digitize(alphabet)
    hmm, _, _ = pyhmmer.plan7.Builder(alphabet).build(sequence, pyhmmer.plan7.Background(alphabet))
    with open(path, "wb") as handle:
        hmm.write(handle)


def _domainator_features(record):
    return [f for f in record.features if f.type == DOMAIN_FEATURE_NAME]


@pytest.mark.parametrize("parser", ["biopython", "lean"])
def test_domainate_partial_cds_codon_start(shared_datadir, monkeypatch, parser):
    monkeypatch.setenv("DOMAINATOR_GB_PARSER", parser)
    with tempfile.TemporaryDirectory() as output_dir:
        hmm = output_dir + "/partial.hmm"
        _write_partial_cds_hmm(shared_datadir, hmm)
        out = output_dir + "/out.gb"
        main(["-i", str(shared_datadir / "JABFVH010000506_extraction.gb"), "-r", hmm, "-o", out])
        record = next(SeqIO.parse(out, "genbank"))
        cds = [f for f in record.features if f.type == "CDS" and f.qualifiers.get("locus_tag") == [_PARTIAL_CDS_LOCUS_TAG]][0]
        assert str(cds.location) == "[802:>2688](-)" # markers survive
        source = [f for f in record.features if f.type == "source"][0]
        assert str(source.location) == "[<0:2688](+)"
        hits = _domainator_features(record)
        assert len(hits) == 1
        # residue 0 starts at codon_start (3), 2 nucleotides into the CDS, which runs from the high end on the minus strand
        assert hits[0].location.end == 2686
        assert (hits[0].location.end - hits[0].location.start) % 3 == 0


@pytest.mark.parametrize("partial,expected_hits", [("include", 1), ("exclude", 0), ("only", 1)])
@pytest.mark.parametrize("parser", ["biopython", "lean"])
def test_domainate_partial_filter(shared_datadir, monkeypatch, parser, partial, expected_hits):
    monkeypatch.setenv("DOMAINATOR_GB_PARSER", parser)
    with tempfile.TemporaryDirectory() as output_dir:
        hmm = output_dir + "/partial.hmm"
        _write_partial_cds_hmm(shared_datadir, hmm)
        out = output_dir + "/out.gb"
        main(["-i", str(shared_datadir / "JABFVH010000506_extraction.gb"), "-r", hmm, "-o", out, "--partial", partial])
        records = list(SeqIO.parse(out, "genbank"))
        assert len(records) == 1 # written with or without hits
        assert len(_domainator_features(records[0])) == expected_hits


def test_domainate_partial_filter_excludes_complete_with_only(shared_datadir):
    # the complete CDS in the fixture is not searched with --partial only
    with tempfile.TemporaryDirectory() as output_dir:
        record = next(_utils.parse_seqfiles([str(shared_datadir / "JABFVH010000506_extraction.gb")]))
        cds = [f for f in record.features if f.type == "CDS" and f.qualifiers.get("locus_tag") != [_PARTIAL_CDS_LOCUS_TAG]][0]
        alphabet = pyhmmer.easel.Alphabet.amino()
        sequence = pyhmmer.easel.TextSequence(name=b"complete_cds", sequence=cds.qualifiers["translation"][0]).digitize(alphabet)
        hmm, _, _ = pyhmmer.plan7.Builder(alphabet).build(sequence, pyhmmer.plan7.Background(alphabet))
        hmm_path = output_dir + "/complete.hmm"
        with open(hmm_path, "wb") as handle:
            hmm.write(handle)
        for partial, expected_hits in (("include", 1), ("exclude", 1), ("only", 0)):
            out = output_dir + f"/out_{partial}.gb"
            main(["-i", str(shared_datadir / "JABFVH010000506_extraction.gb"), "-r", hmm_path, "-o", out, "--partial", partial])
            assert len(_domainator_features(next(_utils.parse_seqfiles([out])))) == expected_hits


@pytest.mark.parametrize("parser", ["biopython", "lean"])
def test_domainate_trims_overlong_partial_translation(shared_datadir, monkeypatch, parser):
    """A 5'-partial CDS whose translation includes residues off the edge of the contig is trimmed to the on-contig residues."""
    monkeypatch.setenv("DOMAINATOR_GB_PARSER", parser)
    with tempfile.TemporaryDirectory() as output_dir:
        record = next(SeqIO.parse(str(shared_datadir / "JABFVH010000506_extraction.gb"), "genbank"))
        cds = [f for f in record.features if f.type == "CDS" and f.qualifiers.get("locus_tag") == [_PARTIAL_CDS_LOCUS_TAG]][0]
        on_contig_translation = cds.qualifiers["translation"][0]
        cds.qualifiers["translation"] = ["MKRFSLAIL" + on_contig_translation]
        overlong = output_dir + "/overlong.gb"
        SeqIO.write([record], overlong, "genbank")
        hmm = output_dir + "/partial.hmm"
        _write_partial_cds_hmm(shared_datadir, hmm)
        out = output_dir + "/out.gb"
        main(["-i", overlong, "-r", hmm, "-o", out])
        out_record = next(SeqIO.parse(out, "genbank"))
        out_cds = [f for f in out_record.features if f.type == "CDS" and f.qualifiers.get("locus_tag") == [_PARTIAL_CDS_LOCUS_TAG]][0]
        assert out_cds.qualifiers["translation"][0] == on_contig_translation
        hits = _domainator_features(out_record)
        assert len(hits) == 1
        assert hits[0].location.end == 2686


def test_domainate_uniprot_fragment_filter(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        records = list(_utils.parse_seqfiles([str(shared_datadir / "swissprot_CuSOD_subset.fasta")]))
        fragment_id = records[0].id
        fasta = output_dir + "/fragments.fasta"
        with open(fasta, "w") as handle:
            for i, record in enumerate(records):
                description = record.description
                if i == 0:
                    description = description.replace(" OS=", " (Fragment) OS=", 1)
                handle.write(f">{description}\n{str(record.seq)}\n")
        alphabet = pyhmmer.easel.Alphabet.amino()
        sequence = pyhmmer.easel.TextSequence(name=b"sod", sequence=str(records[0].seq)).digitize(alphabet)
        hmm, _, _ = pyhmmer.plan7.Builder(alphabet).build(sequence, pyhmmer.plan7.Background(alphabet))
        hmm_path = output_dir + "/sod.hmm"
        with open(hmm_path, "wb") as handle:
            hmm.write(handle)
        hit_ids = dict()
        for partial in ("include", "exclude", "only"):
            out = output_dir + f"/out_{partial}.gb"
            main(["-i", fasta, "-r", hmm_path, "-o", out, "--partial", partial, "--hits_only"])
            hit_ids[partial] = {r.id for r in _utils.parse_seqfiles([out])}
        assert fragment_id in hit_ids["include"]
        assert fragment_id not in hit_ids["exclude"]
        assert hit_ids["only"] == {fragment_id}
        assert hit_ids["include"] == hit_ids["exclude"] | hit_ids["only"]


@pytest.mark.parametrize("input_file", ["pDONR201_multi_genemark.gb", "FeSOD_20.gb", "FeSOD_20.fasta", "pDONR201_multi_genemark_domainator.gb"])
def test_domainate_lean_parser_writes_contigs_without_hits(shared_datadir, monkeypatch, input_file):
    """Without --hits_only, contigs without hits are written too, and the lean parser gives the same output as Biopython."""
    with tempfile.TemporaryDirectory() as output_dir:
        outputs = dict()
        for parser in ("biopython", "lean"):
            monkeypatch.setenv("DOMAINATOR_GB_PARSER", parser)
            outputs[parser] = output_dir + f"/{parser}.gb"
            main(["-i", str(shared_datadir / input_file), "-r", str(shared_datadir / "pdonr_hmms.hmm"), "-o", outputs[parser]])
        monkeypatch.delenv("DOMAINATOR_GB_PARSER")
        # Biopython output first: the lean parser adds topology="linear" when the LOCUS line has none, which is ignored here.
        compare_seqfiles(outputs["biopython"], outputs["lean"])
        assert len(list(SeqIO.parse(outputs["lean"], "genbank"))) == len(list(SeqIO.parse(str(shared_datadir / input_file), "genbank" if input_file.endswith(".gb") else "fasta")))


def test_prodigal_partial_genes(shared_datadir):
    """Gene calling keeps genes that run off the edge of a linear contig, annotated like GenBank partial CDSs."""
    from domainator.domainate import prodigal_CDS_annotate
    record = next(_utils.parse_seqfiles([str(shared_datadir / "JABFVH010000506_extraction.gb")]))
    ncbi = [(str(f.location), f.qualifiers.get("codon_start", ["1"])) for f in record.features if f.type == "CDS"]
    record.features = [f for f in record.features if f.type != "CDS"]
    prodigal_CDS_annotate(record)
    cdss = [f for f in record.features if f.type == "CDS"]
    # same as the NCBI annotation: complement(803..>2688) with /codon_start=3 (prodigal itself stops at the last whole codon, 2686)
    assert [(str(f.location), f.qualifiers.get("codon_start", ["1"])) for f in cdss] == ncbi
    for cds in cdss: # prodigal's translation is of the codons on the contig
        assert cds.qualifiers["translation"][0][1:].rstrip("*") == str(cds.translate(record.seq, cds=False))[1:].rstrip("*")

    # on a circular contig, a gene running off the end continues across the origin, so it isn't partial: it's skipped
    record.features = [f for f in record.features if f.type != "CDS"]
    record.annotations["topology"] = "circular"
    prodigal_CDS_annotate(record)
    assert [str(f.location) for f in record.features if f.type == "CDS"] == ["[279:780](-)"]


@pytest.mark.parametrize("partial,expected_hits", [("include", 1), ("exclude", 0), ("only", 1)])
def test_domainate_gene_call_partial_filter(shared_datadir, partial, expected_hits):
    with tempfile.TemporaryDirectory() as output_dir:
        hmm = output_dir + "/partial.hmm"
        _write_partial_cds_hmm(shared_datadir, hmm)
        out = output_dir + "/out.gb"
        main(["-i", str(shared_datadir / "JABFVH010000506_extraction.gb"), "-r", hmm, "-o", out, "--gene_call", "all", "-Z", "1000", "--partial", partial])
        record = next(SeqIO.parse(out, "genbank"))
        assert any(f.type == "CDS" and str(f.location) == "[802:>2688](-)" for f in record.features)
        hits = _domainator_features(record)
        assert len(hits) == expected_hits
        if expected_hits:
            assert hits[0].location.end == 2686 # placed in the CDS's reading frame (codon_start=3)


def test_prodigal_circular_origin_genes(shared_datadir):
    """On circular contigs, genes crossing the origin are called and annotated as joins across it."""
    from domainator.domainate import prodigal_CDS_annotate
    record = next(_utils.parse_seqfiles([str(shared_datadir / "bacillus_phage_SPR.gb")]))
    assert record.annotations["topology"] == "circular"
    annotated = [str(f.location) for f in record.features if f.type == "CDS"]
    assert "join{[1842:2552](-), [0:130](-)}" in annotated or "join{[0:130](-), [1842:2552](-)}" in annotated
    record.features = [f for f in record.features if f.type != "CDS"]
    prodigal_CDS_annotate(record)
    origin_genes = [f for f in record.features if f.type == "CDS" and len(f.location.parts) == 2]
    assert len(origin_genes) == 1
    # matches the GenBank annotation complement(join(1843..2552,1..130)), with parts in order along the gene
    assert [(int(p.start), int(p.end), p.strand) for p in origin_genes[0].location.parts] == [(0, 130, -1), (1842, 2552, -1)]
    assert origin_genes[0].qualifiers["translation"][0][1:].rstrip("*") == str(origin_genes[0].translate(record.seq, cds=False))[1:].rstrip("*")

    # rotating the contig moves the origin; the gene that crossed it is then called by the main pass, at the same place
    shift = 1000
    rotated = record[shift:] + record[:shift]
    rotated.annotations["topology"] = "circular"
    rotated.features = [f for f in rotated.features if f.type != "CDS"]
    prodigal_CDS_annotate(rotated)
    def gene_ends(rec, offset):
        return sorted(((int(f.location.stranded_start) + offset) % len(rec), (int(f.location.stranded_end) + offset) % len(rec)) for f in rec.features if f.type == "CDS")
    assert gene_ends(rotated, shift) == gene_ends(record, 0)
