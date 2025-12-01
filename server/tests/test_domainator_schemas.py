from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict

import pytest

from domainator_server.app import create_app
from domainator_server.config import ServerConfig
from domainator_server.tool_executor import JobStatus, ToolExecutor
from domainator.Bio import SeqIO


DATA_DIR = Path(__file__).resolve().parents[2] / "test" / "data"
REFERENCE_PEPTIDE_PATH = DATA_DIR / "extract_peptides_test_1_out.gb"

if REFERENCE_PEPTIDE_PATH.exists():
    with REFERENCE_PEPTIDE_PATH.open(encoding="utf-8") as handle:
        _REFERENCE_PEPTIDE_RECORDS = list(SeqIO.parse(handle, "genbank"))
    REFERENCE_PEPTIDE_IDS = {record.id for record in _REFERENCE_PEPTIDE_RECORDS}
    REFERENCE_PEPTIDE_COUNT = len(_REFERENCE_PEPTIDE_RECORDS)
    REFERENCE_PEPTIDE_FIRST_SEQ = str(_REFERENCE_PEPTIDE_RECORDS[0].seq) if _REFERENCE_PEPTIDE_RECORDS else ""
else:  # pragma: no cover - defensive guard if fixture missing
    REFERENCE_PEPTIDE_IDS = set()
    REFERENCE_PEPTIDE_COUNT = 0
    REFERENCE_PEPTIDE_FIRST_SEQ = ""


@dataclass(frozen=True)
class SchemaIntegrationCase:
    tool_id: str
    parameters: Callable[[Path], dict]
    expected_outputs: Dict[str, Callable[[Path], None]]


def _text_contains(expected: str) -> Callable[[Path], None]:
    def _checker(path: Path) -> None:
        content = path.read_text(encoding="utf-8")
        assert expected in content

    return _checker


def _starts_with(prefix: str) -> Callable[[Path], None]:
    def _checker(path: Path) -> None:
        content = path.read_text(encoding="utf-8")
        assert content.startswith(prefix)

    return _checker


def _matches_reference_peptides() -> Callable[[Path], None]:
    def _checker(path: Path) -> None:
        with path.open(encoding="utf-8") as handle:
            records = list(SeqIO.parse(handle, "genbank"))
        assert len(records) == REFERENCE_PEPTIDE_COUNT
        assert {record.id for record in records} == REFERENCE_PEPTIDE_IDS
        if records:
            assert str(records[0].seq) == REFERENCE_PEPTIDE_FIRST_SEQ

    return _checker


def _has_feature_qualifier(qualifier: str) -> Callable[[Path], None]:
    def _checker(path: Path) -> None:
        target = qualifier.lower()
        with path.open(encoding="utf-8") as handle:
            records = list(SeqIO.parse(handle, "genbank"))
        assert any(
            target in {key.lower() for key in feature.qualifiers}
            for record in records
            for feature in record.features
        ), f"expected qualifier '{qualifier}' not present in output"

    return _checker


def _check_trimmed_contig(expected_length: int, expected_ids: set[str]) -> Callable[[Path], None]:
    def _checker(path: Path) -> None:
        with path.open(encoding="utf-8") as handle:
            records = list(SeqIO.parse(handle, "genbank"))
        assert len(records) == 1
        record = records[0]
        assert len(record) == expected_length
        if expected_ids:
            assert record.id in expected_ids

    return _checker


SCHEMA_CASES = [
    SchemaIntegrationCase(
        tool_id="matrix_report",
        parameters=lambda data_dir: {
            "input": str(data_dir / "scorefull.tsv"),
            "output": "report.txt",
            "html": "report.html",
        },
        expected_outputs={
            "report.txt": _text_contains("Matrix Report"),
            "report.html": _text_contains("Matrix Report"),
        },
    ),
    SchemaIntegrationCase(
        tool_id="color_table_to_legend",
        parameters=lambda data_dir: {
            "input": [str(data_dir / "color_table_123.tsv")],
            "svg": "legend.svg",
            "title": "Integration Legend",
        },
        expected_outputs={
            "legend.svg": _text_contains("<svg"),
        },
    ),
    SchemaIntegrationCase(
        tool_id="genbank_to_fasta",
        parameters=lambda data_dir: {
            "input": [str(data_dir / "pDONR201.gb")],
            "output": "output.fasta",
        },
        expected_outputs={
            "output.fasta": _starts_with(">"),
        },
    ),
    SchemaIntegrationCase(
        tool_id="color_genbank",
        parameters=lambda data_dir: {
            "input": [str(data_dir / "color_domain_search_test.gb")],
            "output": "colored.gb",
            "color_both": True,
        },
        expected_outputs={
            "colored.gb": _has_feature_qualifier("Color"),
        },
    ),
    SchemaIntegrationCase(
        tool_id="extract_peptides",
        parameters=lambda data_dir: {
            "input": [str(data_dir / "pDONR201_multi_genemark_domainator.gb")],
            "output": "peptides.gb",
        },
        expected_outputs={
            "peptides.gb": _matches_reference_peptides(),
        },
    ),
    SchemaIntegrationCase(
        tool_id="trim_contigs",
        parameters=lambda data_dir: {
            "input": [str(data_dir / "pDONR201_multi_genemark_domainator.gb")],
            "output": "trimmed.gb",
            "contigs": ["pDONR201_1"],
            "cds_both": 2,
        },
        expected_outputs={
            "trimmed.gb": _check_trimmed_contig(1000, {"pDONR201_1", "pDONR201_1_1266:2265"}),
        },
    ),
]


@pytest.fixture()
def executor(tmp_path) -> ToolExecutor:
    config = ServerConfig(data_dir=tmp_path / "server-data")
    app = create_app(config)
    return app.extensions["domainator_server"]["executor"]


def _run_tool(executor: ToolExecutor, case: SchemaIntegrationCase, data_dir: Path, timeout: float = 120.0):
    job_id = executor.execute_tool(case.tool_id, case.parameters(data_dir))
    deadline = time.time() + timeout
    while time.time() < deadline:
        executor.poll_jobs()
        job = executor.get_job(job_id)
        assert job is not None
        if job.status not in {JobStatus.QUEUED, JobStatus.RUNNING}:
            break
        time.sleep(0.1)

    job = executor.get_job(job_id)
    assert job is not None
    assert job.status == JobStatus.COMPLETED, job.error or f"job failed with status {job.status}"
    return job


@pytest.mark.parametrize("case", SCHEMA_CASES, ids=lambda c: c.tool_id)
def test_domainator_schema_integration(executor: ToolExecutor, case: SchemaIntegrationCase) -> None:
    job = _run_tool(executor, case, DATA_DIR)

    output_set = set(job.output_files)
    for relative, assertion in case.expected_outputs.items():
        assert relative in output_set, f"expected output '{relative}' not found in job outputs"
        output_path = job.work_dir / relative
        assert output_path.exists(), f"output '{relative}' missing from workspace"
        assertion(output_path)

    assert len(job.output_artifacts) >= len(case.expected_outputs)
