import tempfile
import io
from helpers import compare_seqfiles, compare_seqrecords
from domainator.Bio import SeqIO
from domainator import extract_peptides
from domainator import DOMAIN_FEATURE_NAME
from domainator.Bio.Seq import Seq
from domainator.Bio.SeqRecord import SeqRecord
from domainator.Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from domainator.utils import parse_seqfiles

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
        assert recs[0].annotations["source"] == "Bacillus phage SPR"
        assert recs[0].annotations["organism"] == "Bacillus phage SPR"
        assert [feature.type for feature in recs[0].features] == ["source", "CDS", DOMAIN_FEATURE_NAME]
        assert recs[0].features[0].qualifiers["organism"] == ["Bacillus phage SPR"]
        assert recs[0].features[2].qualifiers["name"] == ["Thymidylat_synt"]

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
        extract_peptides.main(["-i", str(shared_datadir / "JABFVH010000506_extraction.gb"), "-o", out, "--strand", "r", "--name_field", "locus_tag"])
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


def _make_extract_peptide_record():
    record = SeqRecord(
        Seq("ATGGCCGCCGCCTTAGGCATGGCC"),
        id="synthetic_contig",
        name="synthetic_contig",
        description="synthetic contig",
    )
    record.annotations["molecule_type"] = "DNA"
    record.annotations["source"] = "synthetic source"
    record.annotations["organism"] = "synthetic organism"
    record.features.append(
        SeqFeature(
            FeatureLocation(0, 12, strand=1),
            type="CDS",
            qualifiers={
                "cds_id": ["cds_f"],
                "protein_id": ["forward_protein"],
                "translation": ["MAAA"],
            },
        )
    )
    record.features.append(
        SeqFeature(
            FeatureLocation(3, 9, strand=1),
            type=DOMAIN_FEATURE_NAME,
            qualifiers={
                "cds_id": ["cds_f"],
                "name": ["forward_domain"],
                "database": ["db"],
                "evalue": ["1e-20"],
                "score": ["50"],
            },
        )
    )
    record.features.append(
        SeqFeature(
            FeatureLocation(12, 24, strand=-1),
            type="CDS",
            qualifiers={
                "cds_id": ["cds_r"],
                "protein_id": ["reverse_protein"],
                "translation": ["GHLP"],
            },
        )
    )
    record.features.append(
        SeqFeature(
            FeatureLocation(15, 21, strand=-1),
            type=DOMAIN_FEATURE_NAME,
            qualifiers={
                "cds_id": ["cds_r"],
                "name": ["reverse_domain"],
                "database": ["db"],
                "evalue": ["1e-20"],
                "score": ["50"],
            },
        )
    )
    return record


def test_extract_peptides_long_protein_id_is_preserved():
    long_name = "protein_" + ("x" * 80)
    record = SeqRecord(Seq("ATGGCCGCCGCC"), id="contig", name="contig", description="contig")
    record.annotations["molecule_type"] = "DNA"
    record.features.append(
        SeqFeature(
            FeatureLocation(0, 12, strand=1),
            type="CDS",
            qualifiers={
                "cds_id": ["cds1"],
                "protein_id": [long_name],
                "translation": ["MAAA"],
            },
        )
    )

    peptides = list(extract_peptides.extract_peptides([record], 1e9, None, None, extract_all=True))

    assert len(peptides) == 1
    assert peptides[0].id == long_name
    assert len(peptides[0].id) == len(long_name)


def test_extract_peptides_translates_domain_locations_on_both_strands():
    peptides = list(extract_peptides.extract_peptides([_make_extract_peptide_record()], 1e9, None, None, extract_all=True))
    peptides_by_id = {record.id: record for record in peptides}

    assert peptides_by_id["forward_protein"].features[0].location == FeatureLocation(1, 3, strand=1)
    assert peptides_by_id["reverse_protein"].features[0].location == FeatureLocation(1, 3, strand=1)
    assert str(peptides_by_id["forward_protein"].seq) == "MAAA"
    assert str(peptides_by_id["reverse_protein"].seq) == "GHLP"


def test_extract_peptides_source_compound_locations_are_translated():
    record = SeqRecord(Seq("ATGGCCGCCGCC"), id="contig", name="contig", description="contig")
    record.annotations["molecule_type"] = "DNA"
    record.annotations["source"] = "synthetic source"
    record.features.append(
        SeqFeature(
            CompoundLocation(
                [FeatureLocation(0, 3, strand=1), FeatureLocation(9, 12, strand=1)],
                operator="join",
            ),
            type="source",
            qualifiers={"organism": ["synthetic organism"]},
        )
    )
    record.features.append(
        SeqFeature(
            FeatureLocation(0, 12, strand=1),
            type="CDS",
            qualifiers={
                "cds_id": ["cds1"],
                "protein_id": ["prot1"],
                "translation": ["MAAA"],
            },
        )
    )

    peptide = list(extract_peptides.extract_peptides([record], 1e9, None, None, extract_all=True))[0]
    source_feature = peptide.features[0]

    assert isinstance(source_feature.location, CompoundLocation)
    assert [(int(part.start), int(part.end)) for part in source_feature.location.parts] == [(0, 1), (3, 4)]