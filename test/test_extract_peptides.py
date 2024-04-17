import tempfile
import io
from helpers import compare_seqfiles, compare_seqrecords
from domainator.Bio import SeqIO
from domainator import extract_peptides
from domainator.utils import parse_seqfiles


#TODO: test long filenames
#TODO: test domain locations on forward and reverse-complement translations

def compare_files(f1,f2):
    with open(f1,"r") as newfile:
        with open(f2, "r") as oldfile:
            return newfile.read() == oldfile.read()

def test_extract_peptides_1(shared_datadir):
    """
        Test read and write genbank files.
    """
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/peptides.gb"
        extract_peptides.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out])
        compare_seqfiles(out, shared_datadir / "extract_peptides_test_1_out.gb" )

def test_extract_peptides_2(shared_datadir):
    """
        Test select only one cds_name and one contig
    """
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/peptides.gb"
        extract_peptides.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--contigs", "pDONR201_1", "pDONR201_3", "--cdss", "pDONR201_2"])
        
        #compare_seqfiles(out, shared_datadir / "extract_peptides_test_1_out.gb")
        recs = list(SeqIO.parse(out, "genbank"))
        assert len(recs) == 2
        compare_seqrecords(recs[0], recs[1], skip_attrs = {"description"}, skip_qualifiers = {"source_contig"})
        assert str(recs[0].seq).upper() == "mqfkvytykresryrlfvdvqsdiidtpgrrmviplasarllsdkvsrelypvvhigdeswrmmttdmasvpvsvigeevadlshrendiknainlmfwgi*".upper()

def test_extract_peptides_3(shared_datadir):
    """
        Test select only one contig, and one domain
    """
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/peptides.gb"
        extract_peptides.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--contigs", "pDONR201_1", "pDONR201_3", "--domains", "CAT"])
        
        #compare_seqfiles(out, shared_datadir / "extract_peptides_test_1_out.gb")
        recs = list(SeqIO.parse(out, "genbank"))
        assert len(recs) == 2
        compare_seqrecords(recs[0], recs[1], skip_attrs = {"description"}, skip_qualifiers = {"source_contig"})
        assert str(recs[0].seq).upper() == "MEKKITGYTTVDISQWHRKEHFEAFQSVAQCTYNQTVQLDITAFLKTVKKNKHKFYPAFIHILARLMNAHPEFRMAMKDGELVIWDSVHPCYTVFHEQTETFSSLWSEYHDDFRQFLHIYSQDVACYGENLAYFPKGFIENMFFVSANPWVSFTSFDLNVANMDNFFAPVFTMGKYYTQGDKVLMPLAIQVHHAVCDGFHVGRMLNELQQYCDEWQGGA*"

def test_extract_peptides_origin_spanning_1(shared_datadir):
    """
        Test that origin-spanning peptides are correctly extracted
    """
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/peptides.gb"
        extract_peptides.main(["-i", str(shared_datadir / "bacillus_phage_SPR_with_annotations.gb"), "-o", out])
        recs = list(SeqIO.parse(out, "genbank"))
        assert len(recs) == 5
        assert str(recs[0].seq).upper() == "MTQFDKQYSSIINDIINNGISDEEFQVRTKWDSDGTPAHTLSVISKQMRFDNSEVPILTTKKVAWKTAIKELLWIWQLKSNDVNDLNKMGVHIWDQWKQEDGTIGHAYGFQLGKKNRSLNGEKVDQVDYLLHQLKNNPSSRRHITMLWNPDDLDSMALTPCVYETQWYVKQGKLHLEVRARSNDMALGNPFNVFQYNVLQRMIAQVTGYKLGEYIFNIGDCHVYTRHIDNLKIQMEREQFEAPELWINPEVKDFYHFTIDDFKLINYKHGDKLLFEVAV"
        assert str(recs[4].seq).upper() == "MLSLIACCDKTLAIGHQNKLLYHVPADMKHFKEKTEGKICIQGRSTYESIIGMTGKPLQNRRNIILTRDQNFKPDYSSFVYHSIEEVLKLIQGQVNTDEEVMVIGGSMIYKAFLPYADKVYLTIVDSESNEADSYFPMLDDHWKVTNKQHNEADEKNKYNYSFLTFENKYRQK"
        #TODO: test that it keeps all the domain and source annotations

def test_domain_search_hit_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        extract_peptides.main(["-i", str(shared_datadir / "pDONR_201_domain_search.gb"), "-o", out, "--search_hits"])
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 1
        new_file[0].name = "pDONR201_1"

def test_extract_peptides_strand_r_1(shared_datadir):
    """
        Test read and write genbank files.
    """
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/peptides.gb"
        extract_peptides.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--strand", "r"])
        recs = list(SeqIO.parse(out, "genbank"))
        assert len(recs) == 12
        rec_ids = {rec.id for rec in recs}
        assert rec_ids == {"pDONR201_2", "pDONR201_3", "pDONR201_4"}

def test_extract_peptides_strand_f_1(shared_datadir):
    """
        Test read and write genbank files.
    """
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/peptides.gb"
        extract_peptides.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--strand", "f"])
        recs = list(SeqIO.parse(out, "genbank"))
        assert len(recs) == 12
        rec_ids = {rec.id for rec in recs}
        assert rec_ids == {"pDONR201_1", "pDONR201_5", "pDONR201_6"}

def test_domain_search_unannotated_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        extract_peptides.main(["-i", str(shared_datadir / "pDONR_201_domain_search.gb"), "-o", out, "--unannotated"])
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 3
        rec_ids = {rec.id for rec in new_file}
        assert rec_ids == {"pDONR201_1", "pDONR201_3", "pDONR201_6"}

def test_domain_search_unannotated_invert_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        extract_peptides.main(["-i", str(shared_datadir / "pDONR_201_domain_search.gb"), "-o", out, "--unannotated", "--invert"])
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 3
        rec_ids = {rec.id for rec in new_file}
        assert rec_ids == {"pDONR201_2", "pDONR201_4", "pDONR201_5"}

def test_extract_peptides_end_of_contig(shared_datadir):
    """
        Test read and write genbank files.
    """
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/peptides.gb"
        extract_peptides.main(["-i", str(shared_datadir / "JABFVH010000506_extraction.gb"), "-o", out, "--strand", "r"])
        recs = list(SeqIO.parse(out, "genbank"))
        assert len(recs) == 2
        rec_ids = {rec.id for rec in recs}
        rec_seqs = {str(rec.seq) for rec in recs}
        assert rec_ids == {"HOV79_30120", "HOV79_30125"}
        assert rec_seqs == {"MILNAIIHREWKYAETSAEMAVLVTQVMEDLRHEGVVTIGGLDYPTAGETARLAFSENRHSGEYDWPAKELVVAVNTETGYGALTWRSVGTVAGGDLVHETWVSENPEPPEADPRVVADPHVPTFHDRRNALPLMRIREALEEYCRVGTGDRPRNIAWTRPVNQRV", 
                            "AEAAAAAAADALAAQQSVAAAAQSAAVAAQEAARASAAAARTAQLNAQAQQDAALAAQYAAEAANEAAAAVAAADEAERDAAAARQAAQQAEDAAAAARAAADQAERDAAAAEEAAQHALADAQQAQDAAALAQENADAQARAALGTSSPTGEAGVQALPRVNAEVVSRTVIQCPPLTDSDYCEFTVTYRITGSIDYVLVTCPDFNDLYCPGEQITDHLKSQPVDMTHEQLVQLTREDINDLLKRLATSLISDYVDCAKGLGILDGKAGETPDNWGVSCAWVAADLVLPAVAAVIARSIKALRIAMRTGDGIVEAYTALKATEISAGTLAKIGDDVYRTVLNFCLGHSFAADTPVLMADHSFKPIGEVVPGDAVLAFDPASGQSVAGTVTRQFVNKDTRLTDVTVRGDNGDLAVLNTTPTHPFWVEGTVSGWVAAQDLDPGDPLRTADGQRATVVSERTFDGAQTMYDLTVEDAHTYYVGVGDAGVLVHNATCPTWVYNALSALRTKALTSGRLFAPNGTELYTEIRSGTDASTSAIDTYLKTVPGWPSNADGFWVATHAETKYAWWMRNNGVKDADIVINNVDGPCAGQYSCPMAVQAILPEGSKVRIWWPGRTTPMEVIGTGVLP"}
        

def test_domain_search_hit_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        extract_peptides.main(["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator_multi_hmm_2.gb"), "-o", out, "--databases", "pdonr_hmms_1", "--domains", "APH", "CcdB"])
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 4
        assert new_file[0].id == "pDONR201_2"
        assert new_file[1].id == "pDONR201_2"
        assert new_file[2].id == "pDONR201_5"
        assert new_file[3].id == "pDONR201_5"

def test_extract_peptides_keep_name_1(shared_datadir):
    """
        Test read and write genbank files.
    """
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/peptides.gb"
        extract_peptides.main(["-i", str(shared_datadir / "JABFVH010000506_extraction.gb"), "-o", out, "--keep_name"])
        recs = list(parse_seqfiles([out]))
        assert len(recs) == 2
        rec_ids = [rec.id for rec in recs]
        rec_seqs = {str(rec.seq) for rec in recs}
        assert rec_ids == ["JABFVH010000506_extraction", "JABFVH010000506_extraction"]
        assert rec_seqs == {"MILNAIIHREWKYAETSAEMAVLVTQVMEDLRHEGVVTIGGLDYPTAGETARLAFSENRHSGEYDWPAKELVVAVNTETGYGALTWRSVGTVAGGDLVHETWVSENPEPPEADPRVVADPHVPTFHDRRNALPLMRIREALEEYCRVGTGDRPRNIAWTRPVNQRV", 
                            "AEAAAAAAADALAAQQSVAAAAQSAAVAAQEAARASAAAARTAQLNAQAQQDAALAAQYAAEAANEAAAAVAAADEAERDAAAARQAAQQAEDAAAAARAAADQAERDAAAAEEAAQHALADAQQAQDAAALAQENADAQARAALGTSSPTGEAGVQALPRVNAEVVSRTVIQCPPLTDSDYCEFTVTYRITGSIDYVLVTCPDFNDLYCPGEQITDHLKSQPVDMTHEQLVQLTREDINDLLKRLATSLISDYVDCAKGLGILDGKAGETPDNWGVSCAWVAADLVLPAVAAVIARSIKALRIAMRTGDGIVEAYTALKATEISAGTLAKIGDDVYRTVLNFCLGHSFAADTPVLMADHSFKPIGEVVPGDAVLAFDPASGQSVAGTVTRQFVNKDTRLTDVTVRGDNGDLAVLNTTPTHPFWVEGTVSGWVAAQDLDPGDPLRTADGQRATVVSERTFDGAQTMYDLTVEDAHTYYVGVGDAGVLVHNATCPTWVYNALSALRTKALTSGRLFAPNGTELYTEIRSGTDASTSAIDTYLKTVPGWPSNADGFWVATHAETKYAWWMRNNGVKDADIVINNVDGPCAGQYSCPMAVQAILPEGSKVRIWWPGRTTPMEVIGTGVLP"}

def test_overhanging_annotations_1(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/extraction.gb"
        extract_peptides.main(["-i", str(shared_datadir / "pDONR_201_domain_search_long_annotations.gb"), "-o", out])
        new_file = list(SeqIO.parse(out, "genbank"))
        assert len(new_file) == 2
        rec_ids = {rec.id for rec in new_file}
        recs = {rec.id:rec for rec in new_file}
        assert rec_ids == {"pDONR201_2", "pDONR201_3"}
        assert len(recs["pDONR201_2"].features) == 2
        assert len(recs["pDONR201_3"].features) == 2
        assert recs["pDONR201_2"].features[0].location.start == 0
        assert recs["pDONR201_2"].features[0].location.end == 102
        assert recs["pDONR201_3"].features[0].location.start == 0
        assert recs["pDONR201_3"].features[0].location.end == 42

# TODO: add a test for source compound locations.