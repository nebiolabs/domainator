from domainator import matrix_report
import tempfile
import pytest

@pytest.mark.parametrize("input_file",
[
    "scorefull.tsv",
    "scorefull.dense.hdf5"
])
def test_matrix_report_1(shared_datadir, input_file):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out_html = output_dir + "/matrix_report_test.html"
        out_txt = output_dir + "/matrix_report_test.txt"
        matrix_report.main(["-i", str(shared_datadir / input_file), "-o", out_txt, "--html", out_html])
        for fh in (out_html, out_txt):
            f_txt = open(fh).read()
            assert "Matrix Report" in f_txt
            assert "Min" in f_txt
            assert "152.0" in f_txt
            assert "451.0" in f_txt
            assert "Total values" in f_txt
            assert f_txt.count("400") == 2


# def test_matrix_report_empty_input(shared_datadir):
#     pass
#     #TODO