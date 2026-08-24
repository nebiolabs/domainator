"""Tests for the agent-friendly JSON outputs and the SSN navigator."""
import json
import os
import tempfile

import numpy as np
import pandas as pd
import pytest

from domainator import enum_report, hmmer_report, summary_report, matrix_report
from domainator import build_ssn_viewer, ssn_navigator
from domainator.data_matrix import DenseDataMatrix
from domainator.ssn_bundle import (
    clusters_at_threshold,
    component_members,
    load_bundle,
    summarize_cluster_metadata,
)

ANNOTATED_GB = "pDONR201_multi_genemark_domainator.gb"


def _read_ndjson(path):
    with open(path) as handle:
        return [json.loads(line) for line in handle if line.strip()]


def _read_json(path):
    with open(path) as handle:
        return json.load(handle)


def test_enum_report_json(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = os.path.join(output_dir, "enum.ndjson")
        tsv = os.path.join(output_dir, "enum.tsv")
        enum_report.main(["-i", str(shared_datadir / ANNOTATED_GB), "-o", tsv,
                          "--domains", "--length", "--json", out])
        records = _read_ndjson(out)
        assert len(records) == 4
        for record in records:
            assert set(record.keys()) == {"contig", "domains", "length"}
            assert isinstance(record["length"], int)
            assert record["domains"] == "2-oxoacid_dh; APH; CAT; CcdA; CcdB; Condensation; TCAD9"


def test_enum_report_json_only(shared_datadir):
    # --json is a valid sole output (no -o / --html required)
    with tempfile.TemporaryDirectory() as output_dir:
        out = os.path.join(output_dir, "enum.ndjson")
        enum_report.main(["-i", str(shared_datadir / ANNOTATED_GB), "--json", out, "--domain_count"])
        records = _read_ndjson(out)
        assert len(records) == 4
        assert all("domain_count" in r for r in records)


def test_hmmer_report_json(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = os.path.join(output_dir, "hmm.ndjson")
        hmmer_report.main(["-i", str(shared_datadir / "FeSOD_pfam.hmm"), "--acc", "--length", "--json", out])
        records = _read_ndjson(out)
        assert len(records) == 2
        for record in records:
            assert set(record.keys()) == {"name", "acc", "length"}
            assert isinstance(record["length"], int)


def test_summary_report_json(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = os.path.join(output_dir, "summary.json")
        summary_report.main(["-i", str(shared_datadir / "FeSOD_20_pfam.gb"), "--json", out])
        payload = _read_json(out)
        assert payload["contig_stats"]["contigs"] == 20
        assert payload["contig_stats"]["length"]["min"] > 0
        domains = {row["domain"] for row in payload["domain_frequency"]}
        assert "Sod_Fe_N" in domains
        assert payload["domain_cooccurrence"] is None


def test_matrix_report_json(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        out = os.path.join(output_dir, "matrix.json")
        matrix_report.main(["-i", str(shared_datadir / "FeSOD_dist.dense.hdf5"), "--json", out])
        payload = _read_json(out)
        assert payload["nodes"] == 20
        assert payload["non_zero_edges"] > 0
        assert payload["edge_scores"]["max"] >= payload["edge_scores"]["min"]
        assert payload["connected_components"] >= 1
        assert isinstance(payload["split_events"], list)


# ---- SSN navigator / bundle ----

def _build_bundle(output_dir):
    data = np.array([
        [0, 10, 6, 0],
        [10, 0, 7, 0],
        [6, 7, 0, 4],
        [0, 0, 4, 0],
    ], dtype=float)
    names = ["A", "B", "C", "D"]
    matrix = DenseDataMatrix(data, names, names)
    matrix_file = os.path.join(output_dir, "m.hdf5")
    matrix.write(matrix_file, output_type="dense")
    meta = os.path.join(output_dir, "meta.tsv")
    pd.DataFrame(
        {"category": ["alpha", "alpha", "beta", "gamma"], "count": [1, 2, 3, 4]},
        index=names,
    ).to_csv(meta, sep="\t")
    bundle_path = os.path.join(output_dir, "net.ssnv")
    build_ssn_viewer.main(["-i", matrix_file, "-o", bundle_path, "--metadata", meta])
    return bundle_path


def test_ssn_bundle_clusters_and_metadata(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        bundle = load_bundle(_build_bundle(output_dir))
        hierarchy = bundle["graph"]["hierarchy"]

        # At the coarsest cut (threshold None) everything is one component.
        coarse = clusters_at_threshold(hierarchy, None)
        assert len(coarse) == 1
        assert sorted(component_members(hierarchy, coarse[0])) == [0, 1, 2, 3]

        # At threshold 5, the 4-4 edge (weight 4) is cut: {A,B,C} and {D}.
        active = clusters_at_threshold(hierarchy, 5.0)
        sizes = sorted(hierarchy["nodes"][cid]["size"] for cid in active)
        assert sizes == [1, 3]

        big = max(active, key=lambda cid: hierarchy["nodes"][cid]["size"])
        members = component_members(hierarchy, big)
        summary = summarize_cluster_metadata(members, bundle["metadata"], top_n=10)
        by_name = {col["name"]: col for col in summary}
        assert by_name["category"]["type"] == "str"
        top = {entry["value"]: entry["count"] for entry in by_name["category"]["top"]}
        assert top["alpha"] == 2
        assert by_name["count"]["type"] == "int"
        assert by_name["count"]["min"] == 1.0
        assert by_name["count"]["max"] == 3.0


def test_ssn_navigator_modes(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        bundle_path = _build_bundle(output_dir)

        def run(args):
            out = os.path.join(output_dir, "out.json")
            ssn_navigator.main(["-i", bundle_path, "-o", out] + args)
            return _read_json(out)

        overview = run(["--mode", "overview"])
        assert overview["nodes"] == 4
        assert {c["name"] for c in overview["metadata_columns"]} == {"category", "count"}

        thresholds = run(["--mode", "thresholds"])
        assert len(thresholds["thresholds"]) >= 1

        clusters = run(["--mode", "clusters", "--threshold", "5"])
        assert clusters["clusters"] == 2
        assert clusters["singletons"] == 1

        # Locate the cluster containing node A and confirm its members.
        cluster = run(["--mode", "cluster", "--threshold", "5", "--node", "A", "--members"])
        assert cluster["size"] == 3
        assert set(cluster["members"]) == {"A", "B", "C"}

        node = run(["--mode", "node", "--node", "A"])
        assert node["metadata"]["category"] == "alpha"
        assert len(node["memberships"]) == len(thresholds["thresholds"])


def test_ssn_navigator_rejects_bad_bundle(shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        bad = os.path.join(output_dir, "bad.ssnv")
        with open(bad, "w") as handle:
            json.dump({"format": "not_a_bundle"}, handle)
        with pytest.raises(ValueError):
            load_bundle(bad)
