from domainator.filter_domains import main, find_kept_annotations, domain_overlap, filter_domains
from pathlib import Path
import tempfile
import pytest
from domainator.Bio import SeqIO
from domainator.Bio.Seq import Seq
from domainator import DOMAIN_FEATURE_NAME

from domainator.Bio.SeqFeature import SeqFeature, FeatureLocation
from domainator.Bio.SeqRecord import SeqRecord

# Mocking a SeqFeature
def mock_seq_feature(start, end, score, evalue, db, name, cds_id):
    location = FeatureLocation(start, end)
    qualifiers = {
        'score': [str(score)],
        'evalue': [str(evalue)],
        'database': [db],
        'name': [name],
        'cds_id': [cds_id]
    }
    return SeqFeature(location=location, type=DOMAIN_FEATURE_NAME, qualifiers=qualifiers)

# Mocking a SeqRecord
def mock_seq_record(features):
    return SeqRecord(seq=Seq("A"*1000),features=features, id="test_record")

# Example test: Check if evalue filtering works correctly
def test_evalue_filtering():
    features = [
        mock_seq_feature(100, 200, 50, 0.001, "db1", "domain1", "cds1"),
        mock_seq_feature(150, 250, 60, 0.1, "db1", "domain2", "cds1")  # Higher evalue
    ]
    record = mock_seq_record(features)
    result = find_kept_annotations(0.05, 1.0, record, None, None, None)
    assert len(result["cds1"]) == 1
    assert list(result["cds1"].values())[0].qualifiers['name'][0] == "domain1"

# Example test: Check for domain overlap
def test_domain_overlap_1():
    features = [
        mock_seq_feature(100, 200, 50, 0.001, "db1", "domain1", "cds1"),
        mock_seq_feature(180, 280, 70, 0.001, "db1", "domain2", "cds1")  # Overlapping domain
    ]
    record = mock_seq_record(features)
    result = find_kept_annotations(0.05, 0.1, record, None, None, None)
    assert len(result["cds1"]) == 1  # Only one should be kept due to overlap

# Example test: Check for domain overlap
def test_domain_db_keep():
    features = [
        mock_seq_feature(100, 200, 50, 0.001, "db1", "domain1", "cds1"),
        mock_seq_feature(180, 280, 70, 0.001, "db2", "domain2", "cds2")  # Overlapping domain
    ]
    record = mock_seq_record(features)
    result = find_kept_annotations(1, 1, record, set(["db1"]), None, None)
    print(result)
    assert len(result["cds1"]) == 1  # Only one should be kept due to overlap
    assert len(result) == 1  # Only one should be kept due to overlap

def test_domain_overlap_1():
    features = [
        mock_seq_feature(100, 200, 50, 0.001, "db1", "domain1", "cds1"),
        mock_seq_feature(150, 250, 60, 0.001, "db1", "domain2", "cds1"),
        mock_seq_feature(300, 400, 70, 0.001, "db1", "domain3", "cds1")
    ]
    hits = {id(f): f for f in features}
    allowed_fraction_overlap = 0.2
    result = domain_overlap(hits, allowed_fraction_overlap)
    assert len(result) == 2  # Expecting 2 domains, as one should be filtered out due to overlap

def test_filter_domains_high_evalue_1():
    features = [
        mock_seq_feature(100, 200, 50, 0.001, "db1", "domain1", "cds1"),
        mock_seq_feature(250, 350, 10, 10, "db1", "domain2", "cds1")  # High evalue
    ]
    record = mock_seq_record(features)
    sequence_iterator = [record]
    result = list(filter_domains(sequence_iterator, 0.05, 1.0, None, None, None))
    assert len(result[0].features) == 1  # Only one feature should be kept due to evalue filtering

def test_main_1(tmpdir):
    # Create a temporary input file
    input_file = tmpdir.join("input.gb")
    input_file.write("mock genbank data")

    # Simulate command-line arguments
    test_args = [
        "--input", str(input_file),
        "--evalue", "0.05",
        "--max_overlap", "1.0"
    ]

    # Call main function
    main(test_args)
    # Add assertions to verify the expected behavior


#TODO: test content of output files!
@pytest.mark.parametrize("input_files,overlap,evalue,input_type,database",
[(["pDONR201_multi_genemark_domainator.gb"], 1, 1000, "gb", None),
 (["pDONR201_multi_genemark_domainator.gb", "pDONR201_multi_genemark_domainator.gb"], 1, "1e-10", "gb", "Pfam-A"),

])
def test_filter_domainator_1(input_files,overlap,evalue,input_type,database,shared_datadir):
    input_file_paths = [str(shared_datadir / x) for x in input_files] 
    arguments = [ "--max_overlap", str(overlap)]
    arguments = arguments + ["--evalue", str(evalue)]
    if input_type == "gb":
        arguments = arguments + ["-i"] + input_file_paths
    if database is not None:
        arguments = arguments + ["--databases_filter", database]
    with tempfile.TemporaryDirectory() as output_dir:
        outfile = output_dir + "/out.gb"
        arguments = arguments + ["-o"] + [outfile]
        main(arguments)
        assert Path(outfile).is_file()


def test_filter_domainator_2(shared_datadir):

    with tempfile.TemporaryDirectory() as output_dir:
        outfile = output_dir + "/out.gb"
        arguments = ["-i", str(shared_datadir / "pDONR_201_domainator.gb"),"-o", outfile, "--evalue", str("1e-150")]
        main(arguments)
        assert Path(outfile).is_file()
        recs = list(SeqIO.parse(outfile, "genbank"))
        assert len(recs) == 1
        for feature in recs[0].features:
            if feature.type == "CDS":
                assert "domainator_Pfam-A" not in feature.qualifiers


#TODO: more tests