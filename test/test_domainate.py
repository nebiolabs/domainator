import os
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
from helpers import compare_seqfiles, compare_seqrecords

#TODO: test evalue cutoffs
#TODO: look into situation where higher scoring domain is smaller than lower-scoring domain, but lower scoring domain is not filtered.
#      see "bad_genbank_overlap.gb" in test data

#TODO: test foldseek, but somehow make it optional, and opt-in

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

        compare_seqfiles(out, ref_file, skip_qualifiers={"identity", "accession"})

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


def test_read_references_nucleotide_fasta_routes_to_nhmmer(shared_datadir):
    references = [str(shared_datadir / "simple_dna_queries.fna")]

    parsed = read_references(references, foldseek=None)

    assert "nhmmer" in parsed
    assert "phmmer" not in parsed
    assert "simple_dna_queries" in parsed["nhmmer"]


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
