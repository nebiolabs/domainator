from pathlib import Path
import tempfile
from glob import glob
from helpers import compare_files
import pytest

from domainator import hmmer_report

def test_hmmer_report_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/hmmer_report.tsv"
        hmmer_report.main(["-i", str(shared_datadir / "pdonr_hmms.hmm"), "-o", out, '--source', '--acc', '--desc', '--length', '--consensus', '--append', 'one', 'int', '1', '--append', 'two', 'float', '2.0', '--append', 'three', 'str', 'three'])
        assert Path(out).is_file()
        
        with open(out) as f:
            lines = f.readlines()
            assert len(lines) == 8
            assert lines[0].strip().split("\t") == ["name","source","acc","desc","length","consensus","one","two","three"]
            assert lines[1].strip() == "2-oxoacid_dh\tpdonr_hmms\tPF00198.25\t2-oxoacid dehydrogenases acyltransferase (catalytic domain)\t233\teqeeervplsgirkaiakrlteskqeiphftlsdevdvtallalrkelkedeakeekakltlldflikavalAlkefPelnasvdeeekeivlkkhvniGvAvatprGLlvPviknadkkslleiakelkelaeraregklkpedleggtftisNlGmlGvtsftPiinppqvaIlgvgrikerpvvkegelvarkvmplslsaDHRvidGaeaarFlntlkkllenpeelll\t1\t2.0\tthree"
            assert lines[-1].strip() == "TCAD9\tpdonr_hmms\tPF19974.1\tTernary complex associated domain 9\t437\tdqvevvrvLtgGrSGaqVlevtvfvkeknqalrhVlKigsaseiakEweAyqrliqpllnalfatIiavsesvlengdqvldelgavvYshagqfagepgeklrsLedlfqealrgpeaadravallerlletllnllYagateeplqtlreelnsrLGpdlvvevkevdseqlvvypdDllqakmssysaseynskvagilvsvelsrlevkvrgprlsavdddvrvevllsggalseleeqgdefleGsvvatranlrlrllkeledelvleetllevdglqlahPfaalrsaLtealearvtssvHGDLNprNiLlaeedrvyLIDfartreggpllsDlAwLevnLlrtvladrldlqellrLqrlLalasrllelealaealagesealakafrllaaiRrfarkqyplerrelwwreylaaLllaahrtLk\t1\t2.0\tthree"
