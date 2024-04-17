from domainator import domainator_db_download
import tempfile
from domainator import utils
from pathlib import Path



def test_uniprot_download_fasta(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        outfile = Path(output_dir) / "uniprot_sprot.fasta"
        domainator_db_download.main(["--db", "swissprot", "--num_recs", "2", "--output", str(outfile)])
        assert outfile.exists()
        assert outfile.stat().st_size > 0
        recs = list(utils.parse_seqfiles([str(outfile)]))
        assert len(recs) == 2

def test_uniprot_download_genbank(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        outfile = Path(output_dir) / "uniprot_sprot.gb"
        domainator_db_download.main(["--db", "swissprot_gb", "--num_recs", "2", "--output", str(outfile)])
        assert outfile.exists()
        assert outfile.stat().st_size > 0
        recs = list(utils.parse_seqfiles([str(outfile)]))
        assert len(recs) == 2

#TODO: add tests for genbank downloads



def test_genbank_download_genbank_1(shared_datadir):
    small_genbanks=['https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/766/775/GCA_008766775.1_ASM876677v1','https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/031/580/395/GCA_031580395.1_ASM3158039v1']
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        outfile = Path(output_dir) / "gb.gb"
        domainator_db_download.process_genbank_accessions([{'ftp_path':small_genbank} for small_genbank in small_genbanks], outfile, gene_call=None, num_recs=1, cpus=3)
        assert outfile.exists()
        # read output file
        recs = list(utils.parse_seqfiles([str(outfile)]))
        assert len(recs) == 1

def test_genbank_download_genbank_2(shared_datadir):
    small_genbanks=['https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/766/775/GCA_008766775.1_ASM876677v1','https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/031/580/395/GCA_031580395.1_ASM3158039v1']
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        outfile = Path(output_dir) / "gb.gb"
        domainator_db_download.process_genbank_accessions([{'ftp_path':small_genbank} for small_genbank in small_genbanks], outfile, gene_call=None, num_recs=None, cpus=3)
        assert outfile.exists()
        # read output file
        recs = list(utils.parse_seqfiles([str(outfile)]))
        assert len(recs) == 2

def test_genbank_download_genbank_3(shared_datadir):
    small_genbanks=['https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/766/775/GCA_008766775.1_ASM876677v1', 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/987/885/GCA_002987885.1_ASM298788v1']
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        outfile = Path(output_dir) / "gb.gb"
        domainator_db_download.process_genbank_accessions([{'ftp_path':small_genbank} for small_genbank in small_genbanks], outfile, gene_call="all", num_recs=None, cpus=2)
        assert outfile.exists()
        # read output file
        recs = list(utils.parse_seqfiles([str(outfile)]))
        assert len(recs) == 2
        outfile_text = outfile.read_text()
        assert "CDS" in outfile_text
        assert '/gene_id="AM260465_1"' in outfile_text

# def test_genbank_download_genbank_skipped_record_log(shared_datadir):
#     small_genbanks=['https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/766/775/GCA_008766775.1_ASM87667v1', 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/987/885/GCA_002987885.1_ASM298788v1']
#     with tempfile.TemporaryDirectory() as output_dir:
#         #output_dir = "test_out"
#         outfile = Path(output_dir) / "gb.gb"
#         skipped_record_log = Path(output_dir) / "skipped_record_log.txt"
#         domainator_db_download.process_genbank_accessions([{'ftp_path':small_genbank} for small_genbank in small_genbanks], outfile, gene_call="all", num_recs=None, cpus=2, skipped_record_log=skipped_record_log)
#         assert outfile.exists()
#         # read output file
#         recs = list(utils.parse_seqfiles([str(outfile)]))
#         assert len(recs) == 1
#         outfile_text = outfile.read_text()
#         assert "CDS" in outfile_text
#         assert '/gene_id="AM260465_1"' in outfile_text
#         assert skipped_record_log.exists()
#         skipped_record_log_text = skipped_record_log.read_text()
#         assert "GCA_008766775.1_ASM87667v1" in skipped_record_log_text
#         assert "GCA_002987885.1_ASM298788v1" not in skipped_record_log_text