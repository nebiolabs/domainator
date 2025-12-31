import os
from domainator.domain_search import main
import tempfile
from glob import glob
from domainator.Bio import SeqIO
from domainator import utils, DOMAIN_FEATURE_NAME
import pytest
from helpers import compare_seqfiles, compare_seqrecords


def test_domain_search_one_file(shared_datadir):
    gb = shared_datadir / "pDONR201.gb"
    hmms = shared_datadir / "CcdB.hmm"
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + f"/out.gb"
        args = ['--input', str(gb) , "-r", str(hmms), "--evalue", str(0.1), "-o", str(out), "--max_overlap", str(1), "-Z", "1000", "--add_annotations"]
        main(args)
        assert os.path.isfile(out)
        
        #assert len(glob(output_dir+"/*.gb")) == 1
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        assert utils.count_peptides_in_record(new_file[0]) == 1
        compare_seqfiles(out, shared_datadir / "ccdb.gb", skip_qualifiers={"accession"})

def test_domain_search_fasta_query(shared_datadir):
    gb = shared_datadir / "pDONR201.gb"
    hmms = shared_datadir / "pdonr_peptides.fasta"
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + f"/out.gb"
        args = ['--input', str(gb) , "-r", str(hmms), "--evalue", str(0.1), "-o", str(out), "--max_overlap", str(1), "-Z", "1000"]
        main(args)
        assert os.path.isfile(out)
        
        #assert len(glob(output_dir+"/*.gb")) == 1
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 3
        assert utils.count_peptides_in_record(new_file[0]) == 1
        assert utils.count_peptides_in_record(new_file[1]) == 1
        assert utils.count_peptides_in_record(new_file[2]) == 1
        #compare_seqfiles(out, shared_datadir / "ccdb.gb")

def test_domain_search_one_file_no_z(shared_datadir):
    gb = shared_datadir / "pDONR201.gb"
    hmms = shared_datadir / "CcdB.hmm"
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir="test_out"
        out = output_dir + f"/out.gb"
        args = ['--input', str(gb) , "-r", str(hmms), "--evalue", str(0.1), "-o", str(out), "--max_overlap", str(1)]
        main(args)
        assert os.path.isfile(out)
        
        #assert len(glob(output_dir+"/*.gb")) == 1
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        assert utils.count_peptides_in_record(new_file[0]) == 1

@pytest.mark.parametrize("max_hits,expected_hits",
[(1,1),
(10000,4),
(2,2),
])
def test_domain_search_max_hits(max_hits, expected_hits, shared_datadir):
    gb = shared_datadir / "pDONR201_multi_genemark.gb"
    hmms = shared_datadir / "CcdB.hmm"
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir="test_out"
        out = output_dir + f"/out.gb"
        args = ['--input', str(gb) , "-r", str(hmms), "--evalue", str(0.1), "-o", str(out), "--max_overlap", str(1), "-Z", "1000", "--max_hits", str(max_hits)]
        main(args)
        assert os.path.isfile(out)
        
        #assert len(glob(output_dir+"/*.gb")) == 1
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == expected_hits
        assert utils.count_peptides_in_record(new_file[0]) == 1

def test_no_annotations_whole_contig_keep_direction(shared_datadir):
    gb = shared_datadir / "pDONR201_multi_genemark.gb"
    hmms = shared_datadir / "CcdB.hmm"
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir="test_out"
        out = output_dir + f"/out.gb"
        args = ['--input', str(gb) , "-r", str(hmms), "--evalue", str(0.1), "-o", str(out), "-Z", "1000", "--whole_contig", "--keep_direction"]
        main(args)
        assert os.path.isfile(out)
        
        compare_seqfiles(out, shared_datadir / "domain_search_test_out1.gb", skip_qualifiers={"accession"})

def test_pad_cds_range(shared_datadir):
    gb = shared_datadir / "pDONR201_multi_genemark.gb"
    hmms = shared_datadir / "pdonr_hmms.hmm"
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir="test_out"
        out = output_dir + f"/out.gb"
        args = ['--input', str(gb) , "-r", str(hmms), "--evalue", str(0.1), "-o", str(out), "-Z", "1000", "--cds_range", str(1), "--pad", "--cpu", "1", "--add_annotations"] 
        main(args)
        assert os.path.isfile(out)
        
        compare_seqfiles(out, shared_datadir / "domain_search_test_out2.gb", skip_qualifiers={"accession"})

def test_pad_cds_range_top_hits(shared_datadir):
    gb = shared_datadir / "pDONR201_multi_genemark.gb"
    hmms = shared_datadir / "pdonr_hmms.hmm"
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir="test_out"
        out = output_dir + f"/out.gb"
        args = ['--input', str(gb) , "-r", str(hmms), "--evalue", str(0.1), "-o", str(out), "-Z", "1000", "--cds_range", str(1), "--pad", "--cpu", "1", "--max_hits", "5", "--add_annotations"]
        main(args)
        assert os.path.isfile(out)
        
        compare_seqfiles(out, shared_datadir / "domain_search_test_out3.gb", skip_qualifiers={"accession"})

def test_deduplicate(shared_datadir):
    gb = shared_datadir / "pDONR201_multi_genemark_domainator.gb"
    hmms = shared_datadir / "pdonr_hmms.hmm"
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir="test_out"
        out = output_dir + f"/out.gb"
        args = ['--input', str(gb) , "-r", str(hmms), "--evalue", str(0.1), "-o", str(out), "-Z", "1000", "--whole_contig", "--cpu", "1", "--deduplicate", "--keep_direction"]
        main(args)
        assert os.path.isfile(out)
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 4
        
        # compare_seqfiles(out, shared_datadir / "pDONR201_multi_genemark_domainator.gb")

def test_domain_search_peptides(shared_datadir):
    input_files = [str(shared_datadir / "simple_genpept.gb"),  str(shared_datadir / "pdonr_peptides.fasta")]
    hmms = shared_datadir / "CcdB.hmm"
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir="test_out"
        out = output_dir + f"/out.gb"
        args = ['--input'] + input_files + ["-r", str(hmms), "--evalue", str(0.1), "-o", str(out), "--max_overlap", str(1), "-Z", "1000"]
        main(args)
        assert os.path.isfile(out)
        
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 2
        assert utils.count_peptides_in_record(new_file[0]) == 1
        assert utils.count_peptides_in_record(new_file[1]) == 1

def test_domain_search_one_file_translate(shared_datadir):
    gb = shared_datadir / "pDONR201.gb"
    hmms = shared_datadir / "CcdB.hmm"
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + f"/out.gb"
        args = ['--input', str(gb) , "-r", str(hmms), "--evalue", str(0.1), "-o", str(out), "--max_overlap", str(1), "-Z", "1000", "--translate"]
        main(args)
        assert os.path.isfile(out)
        
        #assert len(glob(output_dir+"/*.gb")) == 1
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        compare_seqfiles(out, shared_datadir / "domain_search_translate_out.gb", skip_qualifiers={"accession"})

def test_domain_search_gene_annotate_1(shared_datadir):
    hmms = str(shared_datadir / "CcdB.hmm")

    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out1 = output_dir + "/out1.gb"
        query_seqs = shared_datadir / "pDONR201_no_CDSs.gb"
        args = ['--input'] + [str(query_seqs)] + [ "-r"] + [hmms] + ["-o", str(out1)] + ["--gene_call", "all", "-Z", "1000"]
        main(args) # This should run gene calling and domain prediction on the unannotated file

        out2 = output_dir + "/out2.gb"
        query_seqs = shared_datadir / "pDONR201_domainator_circular.gb"
        args = ['--input'] + [str(query_seqs)] + [ "-r"] +[hmms] + ["-o", str(out2)] + ["--gene_call", "all", "-Z", "1000"]
        main(args) # This should re-run gene calling and domain prediction on the unannotated file

        compare_seqfiles(out1, out2, skip_attrs={"id", "description", "name"})

def test_domain_search_gene_annotate_fasta_1(shared_datadir):
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

def test_domain_search_gene_annotate_2(shared_datadir):
    hmms = str(shared_datadir / "CcdB.hmm")

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
        out1_CDSs = utils.DomainatorCDS.list_from_contig(out1[0])
        out2_CDSs = utils.DomainatorCDS.list_from_contig(out2[1])
        assert len(out1_CDSs) == len(out2_CDSs) == 1
        assert 'gene_id' in out1_CDSs[0].feature.qualifiers
        assert 'gene_id' not in out2_CDSs[0].feature.qualifiers

def test_domain_search_min_evalue(shared_datadir):
    gbs = [shared_datadir / "pDONR201_multi_genemark.gb"]
    references = shared_datadir / "pdonr_peptides.fasta"

    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/out1.gb"
        args = ['--input'] + [str(x) for x in gbs] + [ "-r", str(references), "--evalue", str(0.1), "-o", str(out), "--max_overlap", str(1), "--min_evalue", "3.5e-27", "-Z", "0", "--add_annotations"]

        main(args)
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 4
        domainator_features = [x for x in new_file[0].features if x.type == DOMAIN_FEATURE_NAME]
        assert domainator_features[0].qualifiers["name"][0] == "pDONR201_1"

def test_domain_search_taxonomy_1(shared_datadir):
    input = shared_datadir / "swissprot_CuSOD_subset.fasta"
    hmms = shared_datadir / "swissprot_CuSOD_subset.fasta"
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + f"/out.gb"
        args = ['--input', str(input) , "-r", str(hmms), "--evalue", str(0.1), "-o", str(out), "--max_overlap", str(1), "-Z", "1000", "--ncbi_taxonomy_path", str(shared_datadir / "taxdmp"), "--include_taxids", "2", "--exclude_taxids", "1224"]
        main(args)
        assert os.path.isfile(out)
        
        #assert len(glob(output_dir+"/*.gb")) == 1
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        assert utils.count_peptides_in_record(new_file[0]) == 1
        assert new_file[0].id == "sp|O31851|YOJM_BACSU"
        


def test_domain_search_strand_f_1(shared_datadir):
    gb = shared_datadir / "pDONR201.gb"
    hmms = shared_datadir / "pdonr_peptides.fasta"
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + f"/out.gb"
        args = ['--input', str(gb) , "-r", str(hmms), "--evalue", str(0.1), "-o", str(out), "--max_overlap", str(1), "-Z", "1000", "--strand", "f"]
        main(args)
        assert os.path.isfile(out)
        
        #assert len(glob(output_dir+"/*.gb")) == 1
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        assert utils.count_peptides_in_record(new_file[0]) == 1
        #compare_seqfiles(out, shared_datadir / "ccdb.gb")

def test_domain_search_strand_r_1(shared_datadir):
    gb = shared_datadir / "pDONR201.gb"
    hmms = shared_datadir / "pdonr_peptides.fasta"
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + f"/out.gb"
        args = ['--input', str(gb) , "-r", str(hmms), "--evalue", str(0.1), "-o", str(out), "--max_overlap", str(1), "-Z", "1000", "--strand", "r"]
        main(args)
        assert os.path.isfile(out)
        
        #assert len(glob(output_dir+"/*.gb")) == 1
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 2
        assert utils.count_peptides_in_record(new_file[0]) == 1
        #compare_seqfiles(out, shared_datadir / "ccdb.gb")

def test_decoys_nucleotide_1(shared_datadir):
    gb = shared_datadir / "pDONR201_multi_genemark_domainator.gb"
    hmms = shared_datadir / "pdonr_hmms.hmm"
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir="test_out"
        out = output_dir + f"/out.gb"
        args = ['--input', str(gb) , "-r", str(hmms), "--evalue", str(0.1), "-o", str(out), "-Z", "1000", "--cpu", "1", "--keep_direction", "--decoys", "CAT", "APH"]
        main(args)
        assert os.path.isfile(out)
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 8

def test_decoys_translate_1(shared_datadir):
    gb = shared_datadir / "pDONR201_multi_genemark_domainator.gb"
    hmms = shared_datadir / "pdonr_hmms.hmm"
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir="test_out"
        out = output_dir + f"/out.gb"
        args = ['--input', str(gb) , "-r", str(hmms), "--evalue", str(0.1), "-o", str(out), "-Z", "1000", "--cpu", "1", "--keep_direction", "--decoys", "CAT", "APH", "--translate"]
        main(args)
        assert os.path.isfile(out)
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 8

def test_decoys_proteins_1(shared_datadir):
    gb = shared_datadir / "pdonr_peptides.fasta"
    hmms = shared_datadir / "pdonr_hmms.hmm"
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir="test_out"
        out = output_dir + f"/out.gb"
        args = ['--input', str(gb) , "-r", str(hmms), "--evalue", str(0.1), "-o", str(out), "-Z", "1000", "--cpu", "1", "--keep_direction", "--decoys", "CAT", "APH", "--translate"]
        main(args)
        assert os.path.isfile(out)
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 2
# def test_domain_search_annotations_circular_contig(shared_datadir): #TODO: test on circular contigs!
#     input_files = [str(shared_datadir / "bacillus_phage_SPR.gb")]
#     hmms = shared_datadir / "thymidylate_synthase.fasta"
#     with tempfile.TemporaryDirectory() as output_dir:
#         # output_dir="test_out"
#         out = output_dir + f"/out.gb"
#         args = ['--input'] + input_files + ["-r", str(hmms), "--evalue", str(0.1), "-o", str(out), "--max_overlap", str(1), "-Z", "1000"]
#         main(args)
#         assert os.path.isfile(out)
        
#         new_file = list(SeqIO.parse(out, "genbank"))


@pytest.mark.parametrize("params, expected_hit_size",
[(["--kb_up", "0.0", "--kb_down", "3.0"],3306),
(["--kb_up", "3.0", "--kb_down", "0.0"],3306),
(["--kb_up", "0.0", "--kb_down", "100.0"],4470),
(["--kb_up", "100.0", "--kb_down", "0.0"],4470),
(["--cds_up", "1", "--cds_down", "0"],1307),
(["--cds_up", "0", "--cds_down", "1"],2867),
(["--cds_up", "100", "--cds_down", "0"],2719),
(["--cds_up", "0", "--cds_down", "100"],4129),
])
def test_domain_search_range_up_down_1(params, expected_hit_size, shared_datadir):
    gb = shared_datadir / "pDONR201.gb"
    hmms = shared_datadir / "CcdB.hmm"
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir="test_out"
        out = output_dir + f"/out.gb"
        args = ['--input', str(gb) , "-r", str(hmms), "--evalue", str(0.1), "-o", str(out), "--max_overlap", str(1), "-Z", "1000", "--add_annotations" ] + params
        main(args)
        assert os.path.isfile(out)
        
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        assert len(new_file[0]) == expected_hit_size

#TODO: test keeping annotations like source, organism, etc.


def test_domain_search_single_pass(shared_datadir):
    """Test the single-pass mode which uses raw text partitioning."""
    gb = shared_datadir / "pDONR201.gb"
    hmms = shared_datadir / "CcdB.hmm"
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + f"/out.gb"
        args = ['--input', str(gb) , "-r", str(hmms), "--evalue", str(0.1), "-o", str(out), 
                "--max_overlap", str(1), "-Z", "1000", "--add_annotations", "--single_pass"]
        main(args)
        assert os.path.isfile(out)
        
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        assert utils.count_peptides_in_record(new_file[0]) == 1
        compare_seqfiles(out, shared_datadir / "ccdb.gb", skip_qualifiers={"accession"})


def test_domain_search_single_pass_no_z(shared_datadir):
    """Test single-pass mode with Z=0 falls back to offset-based mode with warning."""
    gb = shared_datadir / "pDONR201.gb"
    hmms = shared_datadir / "CcdB.hmm"
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + f"/out.gb"
        args = ['--input', str(gb) , "-r", str(hmms), "--evalue", str(0.1), "-o", str(out), 
                "--max_overlap", str(1), "-Z", "0", "--single_pass"]
        # Should issue a warning about single_pass with Z=0
        import warnings
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            main(args)
            # Check that a warning was issued
            assert any("single_pass" in str(warning.message).lower() for warning in w)
        
        assert os.path.isfile(out)
        
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        assert utils.count_peptides_in_record(new_file[0]) == 1


def test_domain_search_single_pass_multiple_files(shared_datadir):
    """Test single-pass mode with multiple input files."""
    gb1 = shared_datadir / "pDONR201.gb"
    gb2 = shared_datadir / "pDONR201_multi_genemark.gb"
    hmms = shared_datadir / "CcdB.hmm"
    with tempfile.TemporaryDirectory() as output_dir:
        out = output_dir + f"/out.gb"
        args = ['--input', str(gb1), str(gb2), "-r", str(hmms), "--evalue", str(0.1), 
                "-o", str(out), "--max_overlap", str(1), "-Z", "1000", "--single_pass"]
        main(args)
        assert os.path.isfile(out)
        
        new_file = list(SeqIO.parse(out, "genbank"))
        # Should have 1 hit from pDONR201 and 4 hits from pDONR201_multi_genemark
        assert len(new_file) == 5