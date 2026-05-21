import io
import os
import pytest
import shutil
import subprocess
from domainator import utils
from pathlib import Path
import tempfile
from domainator import deduplicate_genbank
from domainator.Bio import Seq


def _probe_diamond_cluster(bin_path):
    with tempfile.TemporaryDirectory() as output_dir:
        input_fasta = Path(output_dir) / "input.fasta"
        input_fasta.write_text(
            ">seq1\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR\n"
            ">seq2\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR\n"
            ">seq3\nMVLSAADKTNVKAAWGKVGGHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYA\n",
            encoding="utf-8",
        )
        cluster_output = Path(output_dir) / "clusters.tsv"

        try:
            subprocess.run(
                [bin_path, "cluster", "-d", str(input_fasta), "-o", str(cluster_output), "--approx-id", "99", "--threads", "1"],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        except (FileNotFoundError, subprocess.CalledProcessError):
            return False

    return True


def _find_working_diamond():
    candidates = list()
    which_diamond = shutil.which("diamond")
    if which_diamond is not None:
        candidates.append(which_diamond)

    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix is not None:
        prefix_path = Path(conda_prefix).resolve()
        package_dirs = [prefix_path / "pkgs", prefix_path.parent / "pkgs"]
        for package_dir in package_dirs:
            if package_dir.is_dir():
                candidates.extend(str(path) for path in sorted(package_dir.glob("diamond-*/bin/diamond"), reverse=True))

    seen = set()
    for candidate in candidates:
        if candidate in seen:
            continue
        seen.add(candidate)
        if _probe_diamond_cluster(candidate):
            return candidate

    return None


def test_deduplicate_genbank_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/deduplicate_out.gb"
        cluster_table_out = output_dir + "/cluster_table.tsv"

        deduplicate_genbank.main(["-i", str(shared_datadir / "extract_peptides_test_2.gb"), "--id", "0.99", "-o", out, "--params", "\"-s\":0.9", "--prefix_count", "--cluster_table", cluster_table_out])
        assert Path(out).is_file()
        recs  = list(utils.parse_seqfiles([out]))
        for rec in recs:
            assert rec.id[:2] == "4-"
        assert len(recs) == 6
#cluster_table_out:
#         representative	contigs
# 4-pDONR201_5_1	pDONR201_5_1 ; pDONR201_5_2 ; pDONR201_5_3 ; pDONR201_5_4
# 4-pDONR201_4_1	pDONR201_4_1 ; pDONR201_4_2 ; pDONR201_4_3 ; pDONR201_4_4
# 4-pDONR201_2_1	pDONR201_2_1 ; pDONR201_2_2 ; pDONR201_2_3 ; pDONR201_2_4
# 4-pDONR201_3_1	pDONR201_3_1 ; pDONR201_3_2 ; pDONR201_3_3 ; pDONR201_3_4
# 4-pDONR201_1_1	pDONR201_1_1 ; pDONR201_1_2 ; pDONR201_1_3 ; pDONR201_1_4
# 4-pDONR201_6_1	pDONR201_6_1 ; pDONR201_6_2 ; pDONR201_6_3 ; pDONR201_6_4
        assert Path(cluster_table_out).is_file()
        cluster_table = open(cluster_table_out, "r").read().split("\n")
        assert "representative" in cluster_table[0]
        assert "contigs" in cluster_table[0]
        assert "4-pDONR201_5_1\tpDONR201_5_1 ; pDONR201_5_2 ; pDONR201_5_3 ; pDONR201_5_4" in cluster_table
        assert "4-pDONR201_4_1\tpDONR201_4_1 ; pDONR201_4_2 ; pDONR201_4_3 ; pDONR201_4_4" in cluster_table
        assert len(cluster_table) == 8




def test_deduplicate_genbank_cdhit_low_id_1(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/deduplicate_out.gb"
        
        deduplicate_genbank.main(["-i", str(shared_datadir / "FeSOD_20.gb"), "--id", "0.4", "--cpu", "1", "-o", out, "--prefix_count"])
        assert Path(out).is_file()
        recs  = list(utils.parse_seqfiles([out]))
        assert recs[0].id[:3] == "10-"
        assert len(recs) == 4

def test_deduplicate_genbank_cdhit_low_id_2(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/deduplicate_out.gb"
        
        deduplicate_genbank.main(["-i", str(shared_datadir / "FeSOD_20.gb"), "--algorithm", "usearch", "--id", "0.4", "--cpu", "1", "-o", out, "--prefix_count"])
        assert Path(out).is_file()
        recs  = list(utils.parse_seqfiles([out]))
        assert recs[0].id[:2] == "6-"
        assert len(recs) == 5

def test_deduplicate_genbank_usearch(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/deduplicate_out.gb"
        log = output_dir + "/deduplicate_out.log"
        cluster_table_out = output_dir + "/cluster_table.tsv"

        deduplicate_genbank.main( ["-i", str(shared_datadir / "extract_peptides_test_2.gb"), "--algorithm", "usearch", "--id", "0.3", "-o", out, "--cpu", "1", "--suffix_count", "--log", log, "--cluster_table", cluster_table_out] )
        assert Path(out).is_file()
        recs  = list(utils.parse_seqfiles([out]))
        for rec in recs:
            assert rec.id[-2:] == "-4"
        assert len(recs) == 6

        log_file = open(log, "r").read()
        assert "Input  size: 24" in log_file
        assert "Output size: 6" in log_file

        #cluster_table_out:
#         representative	contigs
# pDONR201_5_1-4	pDONR201_5_1 ; pDONR201_5_2 ; pDONR201_5_3 ; pDONR201_5_4
# pDONR201_4_1-4	pDONR201_4_1 ; pDONR201_4_2 ; pDONR201_4_3 ; pDONR201_4_4
# pDONR201_2_1-4	pDONR201_2_1 ; pDONR201_2_2 ; pDONR201_2_3 ; pDONR201_2_4
# pDONR201_3_1-4	pDONR201_3_1 ; pDONR201_3_2 ; pDONR201_3_3 ; pDONR201_3_4
# pDONR201_1_1-4	pDONR201_1_1 ; pDONR201_1_2 ; pDONR201_1_3 ; pDONR201_1_4
# pDONR201_6_1-4	pDONR201_6_1 ; pDONR201_6_2 ; pDONR201_6_3 ; pDONR201_6_4
        assert Path(cluster_table_out).is_file()
        cluster_table = open(cluster_table_out, "r").read().split("\n")
        assert "representative" in cluster_table[0]
        assert "contigs" in cluster_table[0]
        assert "pDONR201_5_1-4\tpDONR201_5_1 ; pDONR201_5_2 ; pDONR201_5_3 ; pDONR201_5_4" in cluster_table
        assert "pDONR201_4_1-4\tpDONR201_4_1 ; pDONR201_4_2 ; pDONR201_4_3 ; pDONR201_4_4" in cluster_table
        assert len(cluster_table) == 8

def test_deduplicate_genbank_hash(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/deduplicate_out.gb"
        log = output_dir + "/deduplicate_out.log"
        cluster_table_out = output_dir + "/cluster_table.tsv"
        
        deduplicate_genbank.main(["-i", str(shared_datadir / "extract_peptides_test_2.gb"), "--id", "1", "--algorithm", "hash", "-o", out, "--prefix_count", "--log", log, "--cluster_table", cluster_table_out])
        assert Path(out).is_file()
        recs  = list(utils.parse_seqfiles([out]))
        for rec in recs:
            assert rec.id[:2] == "4-"
        assert len(recs) == 6
        
        log_file = open(log, "r").read()
        assert "Input  size: 24" in log_file
        assert "Output size: 6" in log_file

        #cluster_table_out:
#         representative	contigs
# 4-pDONR201_5_1	pDONR201_5_1 ; pDONR201_5_2 ; pDONR201_5_3 ; pDONR201_5_4
# 4-pDONR201_4_1	pDONR201_4_1 ; pDONR201_4_2 ; pDONR201_4_3 ; pDONR201_4_4
# 4-pDONR201_2_1	pDONR201_2_1 ; pDONR201_2_2 ; pDONR201_2_3 ; pDONR201_2_4
# 4-pDONR201_3_1	pDONR201_3_1 ; pDONR201_3_2 ; pDONR201_3_3 ; pDONR201_3_4
# 4-pDONR201_1_1	pDONR201_1_1 ; pDONR201_1_2 ; pDONR201_1_3 ; pDONR201_1_4
# 4-pDONR201_6_1	pDONR201_6_1 ; pDONR201_6_2 ; pDONR201_6_3 ; pDONR201_6_4
        assert Path(cluster_table_out).is_file()
        cluster_table = open(cluster_table_out, "r").read().split("\n")
        assert "representative" in cluster_table[0]
        assert "contigs" in cluster_table[0]
        assert "4-pDONR201_5_1\tpDONR201_5_1 ; pDONR201_5_2 ; pDONR201_5_3 ; pDONR201_5_4" in cluster_table
        assert "4-pDONR201_4_1\tpDONR201_4_1 ; pDONR201_4_2 ; pDONR201_4_3 ; pDONR201_4_4" in cluster_table
        assert len(cluster_table) == 8


def test_deduplicate_genbank_hash_both_strands(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        input_fasta = output_dir + "/input.fasta"
        out = output_dir + "/deduplicate_out.fasta"
        log = output_dir + "/deduplicate_out.log"

        seq = "ACGTAAGT"
        rev_comp = str(Seq.Seq(seq).reverse_complement())

        with open(input_fasta, "w") as handle:
            handle.write(f">seq1\n{seq}\n")
            handle.write(f">seq2\n{rev_comp}\n")

        deduplicate_genbank.main([
            "-i", input_fasta,
            "--id", "1",
            "--algorithm", "hash",
            "--fasta_type", "nucleotide",
            "--both_strands",
            "--prefix_count",
            "--fasta_out",
            "-o", out,
            "--log", log,
        ])

        assert Path(out).is_file()
        recs = list(utils.parse_seqfiles([out]))
        assert len(recs) == 1
        assert recs[0].id[:2] == "2-"


def test_deduplicate_genbank_usearch_both_strands_sets_flag():
    with tempfile.TemporaryDirectory() as output_dir:
        input_fasta = output_dir + "/input.fasta"

        with open(input_fasta, "w") as handle:
            handle.write(">seq1\nACGTAAGT\n")

        captured = {}

        def fake_run_clustering(algorithm, input_fasta_path, id, params, bin_path=None, add_count=None, cpus=1, log_handle=None):
            captured["algorithm"] = algorithm
            captured["params"] = dict(params)
            return {"0"}, {"0": ["0"]}

        original_run_clustering = deduplicate_genbank.run_clustering
        deduplicate_genbank.run_clustering = fake_run_clustering
        try:
            dedup = deduplicate_genbank.DeduplicateGenbank()
            recs = list(dedup.deduplicate_genbank(
                [input_fasta],
                "usearch",
                {},
                1,
                None,
                None,
                1,
                "nucleotide",
                both_strands=True,
            ))
        finally:
            deduplicate_genbank.run_clustering = original_run_clustering

        assert captured["algorithm"] == "usearch"
        assert captured["params"].get("-strand") == "both"
        assert len(recs) == 1


def test_diamond_builds_cluster_command():
    search = deduplicate_genbank.Diamond(
        "/tmp/input.fasta",
        0.99,
        {"--member-cover": 90},
        bin_path="/opt/diamond",
        cpus=3,
    )

    assert search.search_args[:2] == ["/opt/diamond", "cluster"]
    assert search.search_args[search.search_args.index("-d") + 1] == "/tmp/input.fasta"
    assert search.search_args[search.search_args.index("-o") + 1] == search.cluster_assignment_path
    assert search.search_args[search.search_args.index("--approx-id") + 1] == "99.0"
    assert search.search_args[search.search_args.index("--threads") + 1] == "3"
    assert search.search_args[-2:] == ["--member-cover", "90"]


def test_diamond_run_falls_back_to_longest_sequence(monkeypatch):
    with tempfile.TemporaryDirectory() as output_dir:
        input_fasta = Path(output_dir) / "input.fasta"
        input_fasta.write_text(
            ">seq1\nMFTLPPLPYPTNALEPYLDTQTLEIHFGKHHATYLKNLNDLLPEKSDADLIPVLQHLDDLPQDIRVKVRNNAGGVYNHNLYWQCMSPKSKSPSPRLLSSIESGFGTLDAFKEKFSQAALTHFGSGWAWLVKGTKGLEIVTTPNQDSPVSTGLTPILGLDVWEHAYYLKYQNRRVEYIQAWWNVVNWDYVSSLLADR\n"
            ">seq2\nMFTLPPLPYPTNALEPYLDTQTLEIHFGKHHATYLKNLNDLLPEKSDADLIPVLQHLDDLPQDIRVKVRNNAGGVYNHNLYWQCMSPKSKSPSPRLLSSIESGFGTLDAFKEKFSQAALTHFGSGWAWLVKGTKGLEIVTTPNQDSPVSTGLTPILGLDVWEHAYYLKYQNRRVEYIQAWWNVVNWDYVSSLLADRAAAAA\n"
            ">seq3\nMFTLPPLPYPTNALEPYLDTQTLEIHFGKHHATYLKNLNDLLPEKSDADLIPVLQHLDDLPQDIRVKVRNNAGGVYNHNLYWQCMSPKSKSPSPRLLSSIESGFGTLDAFKEKFSQAALTHFGSGWAWLVKGTKGLEIVTTPNQDSPVSTGLTPILGLDVWEHAYYLKYQNRRVEYIQAWWNVVNWDYVSSLL\n",
            encoding="utf-8",
        )

        log_handle = io.StringIO()
        search = deduplicate_genbank.Diamond(str(input_fasta), 0.0, {}, cpus=1, log_handle=log_handle)
        calls = list()

        def fake_run(args, stdout=None, stderr=None, encoding=None):
            calls.append(list(args))
            return subprocess.CompletedProcess(args, 1, "", "Error: Record count not set\n")

        monkeypatch.setattr(deduplicate_genbank.subprocess, "run", fake_run)

        search.run()

        assert len(calls) == 1
        assert calls[0][calls[0].index("-d") + 1] == str(input_fasta)
        assert search.get_cluster_members() == {"seq2": ["seq1", "seq2", "seq3"]}
        assert "Warning: diamond cluster failed with 'Record count not set'" in log_handle.getvalue()


def test_diamond_parses_cluster_table():
    with tempfile.TemporaryDirectory() as output_dir:
        input_fasta = Path(output_dir) / "input.fasta"
        input_fasta.write_text(">0\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR\n", encoding="utf-8")

        search = deduplicate_genbank.Diamond(str(input_fasta), 0.99, {}, cpus=1)
        Path(search.cluster_assignment_path).write_text(
            "representative\tmember\n1\t0\n2\t2\n",
            encoding="utf-8",
        )

        assert search.get_cluster_members() == {"1": ["1", "0"], "2": ["2"]}
        assert search.get_reduced_names() == {"1", "2"}


def test_deduplicate_genbank_diamond_rejects_nucleotide_input():
    with tempfile.TemporaryDirectory() as output_dir:
        input_fasta = Path(output_dir) / "input.fasta"
        input_fasta.write_text(">seq1\nACGTAAGT\n", encoding="utf-8")

        dedup = deduplicate_genbank.DeduplicateGenbank()
        with pytest.raises(ValueError, match="diamond deduplication is only supported for protein input"):
            list(dedup.deduplicate_genbank(
                [str(input_fasta)],
                "diamond",
                {},
                0.99,
                None,
                None,
                1,
                "nucleotide",
            ))


def test_deduplicate_genbank_diamond():
    working_diamond = _find_working_diamond()
    if working_diamond is None:
        pytest.skip("No working diamond cluster binary available")

    with tempfile.TemporaryDirectory() as output_dir:
        input_fasta = Path(output_dir) / "input.fasta"
        out = Path(output_dir) / "deduplicate_out.fasta"
        log = Path(output_dir) / "deduplicate_out.log"
        cluster_table_out = Path(output_dir) / "cluster_table.tsv"

        input_fasta.write_text(
            ">seq1\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR\n"
            ">seq2\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR\n"
            ">seq3\nMVLSAADKTNVKAAWGKVGGHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYA\n",
            encoding="utf-8",
        )

        deduplicate_genbank.main([
            "-i", str(input_fasta),
            "--algorithm", "diamond",
            "--bin_path", working_diamond,
            "--id", "0.99",
            "--cpu", "1",
            "--prefix_count",
            "--fasta_out",
            "-o", str(out),
            "--log", str(log),
            "--cluster_table", str(cluster_table_out),
        ])

        assert out.is_file()
        recs = list(utils.parse_seqfiles([str(out)]))
        assert len(recs) == 2
        assert any(rec.id.startswith("2-") for rec in recs)
        assert any(rec.id.startswith("1-") for rec in recs)

        cluster_table = cluster_table_out.read_text(encoding="utf-8").splitlines()
        assert cluster_table[0] == "representative\tcontigs"
        cluster_members = [set(line.split("\t")[1].split(" ; ")) for line in cluster_table[1:] if line]
        assert {"seq1", "seq2"} in cluster_members
        assert {"seq3"} in cluster_members

        log_file = log.read_text(encoding="utf-8")
        assert "Input  size: 3" in log_file
        assert "Output size: 2" in log_file


def test_deduplicate_genbank_diamond_single_centroid():
    working_diamond = _find_working_diamond()
    if working_diamond is None:
        pytest.skip("No working diamond cluster binary available")

    with tempfile.TemporaryDirectory() as output_dir:
        input_fasta = Path(output_dir) / "input.fasta"
        out = Path(output_dir) / "deduplicate_out.fasta"
        log = Path(output_dir) / "deduplicate_out.log"
        cluster_table_out = Path(output_dir) / "cluster_table.tsv"

        input_fasta.write_text(
            ">seq1\nMFTLPPLPYPTNALEPYLDTQTLEIHFGKHHATYLKNLNDLLPEKSDADLIPVLQHLDDLPQDIRVKVRNNAGGVYNHNLYWQCMSPKSKSPSPRLLSSIESGFGTLDAFKEKFSQAALTHFGSGWAWLVKGTKGLEIVTTPNQDSPVSTGLTPILGLDVWEHAYYLKYQNRRVEYIQAWWNVVNWDYVSSLLADR\n"
            ">seq2\nMFTLPPLPYPTNALEPYLDTQTLEIHFGKHHATYLKNLNDLLPEKSDADLIPVLQHLDDLPQDIRVKVRNNAGGVYNHNLYWQCMSPKSKSPSPRLLSSIESGFGTLDAFKEKFSQAALTHFGSGWAWLVKGTKGLEIVTTPNQDSPVSTGLTPILGLDVWEHAYYLKYQNRRVEYIQAWWNVVNWDYVSSLLADR\n",
            encoding="utf-8",
        )

        deduplicate_genbank.main([
            "-i", str(input_fasta),
            "--algorithm", "diamond",
            "--bin_path", working_diamond,
            "--id", "0",
            "--cpu", "1",
            "--prefix_count",
            "--fasta_out",
            "-o", str(out),
            "--log", str(log),
            "--cluster_table", str(cluster_table_out),
        ])

        assert out.is_file()
        recs = list(utils.parse_seqfiles([str(out)]))
        assert len(recs) == 1
        assert recs[0].id == "2-seq1"

        cluster_table = cluster_table_out.read_text(encoding="utf-8").splitlines()
        assert cluster_table[0] == "representative\tcontigs"
        assert cluster_table[1].endswith("seq1 ; seq2") or cluster_table[1].endswith("seq2 ; seq1")

        log_file = log.read_text(encoding="utf-8")
        assert "Warning: diamond cluster failed with 'Record count not set'" in log_file
        assert "Input  size: 2" in log_file
        assert "Output size: 1" in log_file


def test_deduplicate_genbank_cdhit_est_both_strands_sets_flag():
    with tempfile.TemporaryDirectory() as output_dir:
        input_fasta = output_dir + "/input.fasta"

        with open(input_fasta, "w") as handle:
            handle.write(">seq1\nACGTAAGT\n")

        captured = {}

        def fake_run_clustering(algorithm, input_fasta_path, id, params, bin_path=None, add_count=None, cpus=1, log_handle=None):
            captured["algorithm"] = algorithm
            captured["params"] = dict(params)
            return {"0"}, {"0": ["0"]}

        original_run_clustering = deduplicate_genbank.run_clustering
        deduplicate_genbank.run_clustering = fake_run_clustering
        try:
            dedup = deduplicate_genbank.DeduplicateGenbank()
            recs = list(dedup.deduplicate_genbank(
                [input_fasta],
                "cd-hit",
                {},
                0.8,
                None,
                None,
                1,
                "nucleotide",
                both_strands=True,
            ))
        finally:
            deduplicate_genbank.run_clustering = original_run_clustering

        assert captured["algorithm"] == "cd-hit-est"
        assert captured["params"].get("-r") == 1
        assert len(recs) == 1


def test_deduplicate_genbank_cdhit_low_id_nt(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/deduplicate_out.gb"
        
        deduplicate_genbank.main(["-i", str(shared_datadir / "pDONR201_multi_genemark.gb"), "--algorithm", "cd-hit", "--id", "0.8", "--cpu", "1", "-o", out, "--prefix_count"])
        assert Path(out).is_file()
        recs  = list(utils.parse_seqfiles([out]))
        assert recs[0].id[:2] == "4-"
        assert len(recs) == 1

def test_deduplicate_genbank_cdhit_fasta_input(shared_datadir):
    
    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "test_out"
        out = output_dir + "/deduplicate_out.gb"
        log = output_dir + "/deduplicate_out.log"
        
        deduplicate_genbank.main(["-i", str(shared_datadir / "FeSOD_20.fasta"), "--algorithm", "cd-hit", "--id", "0.8", "--cpu", "1", "-o", out, "--prefix_count", "--log", log])
        assert Path(out).is_file()
        recs  = list(utils.parse_seqfiles([out]))
        assert recs[0].id[:2] == "1-"
        assert len(recs) == 20

        log_file = open(log, "r").read()
        assert "Input  size: 20" in log_file
        assert "Output size: 20" in log_file

def test_deduplicate_genbank_cdhit_fasta_input_2(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        # output_dir = "test_out"
        out = output_dir + "/deduplicate_out.fasta"
        log = output_dir + "/deduplicate_out.log"
        
        deduplicate_genbank.main(["-i", str(shared_datadir / "pDONR201_multi.fasta"), "--algorithm", "cd-hit", "--id", "0.8", "--cpu", "1", "-o", out, "--prefix_count", "--fasta_type", "nucleotide", "--fasta_out", "--log", log])
        assert Path(out).is_file()
        recs  = list(utils.parse_seqfiles([out]))
        assert recs[0].id[:2] == "4-"
        assert len(recs) == 1

        log_file = open(log, "r").read()
        assert "Input  size: 4" in log_file
        assert "Output size: 1" in log_file
