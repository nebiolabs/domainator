from __future__ import annotations

import json
import stat
import time
from pathlib import Path

import pytest

from domainator_server.config import ServerConfig
from domainator_server.file_manager import FileManager
from domainator_server.tool_executor import JobStatus, ToolExecutor
from domainator_server.tool_registry import ToolRegistry
from werkzeug.utils import secure_filename


SCRIPT_TEMPLATE = """#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--config', required=True)
args = parser.parse_args()
with open(args.config, 'r', encoding='utf-8') as handle:
    data = json.load(handle)
output = data.get('output', 'result.txt')
Path(output).write_text('done')
"""


@pytest.fixture()
def schema_dir(tmp_path: Path) -> Path:
    directory = tmp_path / "schemas"
    directory.mkdir()
    return directory


@pytest.fixture()
def tool_script(tmp_path: Path) -> Path:
    script = tmp_path / "echo_tool.py"
    script.write_text(SCRIPT_TEMPLATE)
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    return script


def _write_schema(schema_dir: Path, tool_id: str, entry_point: Path) -> None:
    payload = {
        "id": tool_id,
        "runner": "python",
        "entry_point": str(entry_point),
        "parameters": [
            {"name": "output", "type": "string"},
        ],
    }
    with (schema_dir / f"{tool_id}.json").open("w", encoding="utf-8") as handle:
        json.dump(payload, handle)


def test_tool_executor_runs_process(tmp_path: Path, schema_dir: Path, tool_script: Path) -> None:
    _write_schema(schema_dir, "echo", tool_script)
    config_dir = tmp_path / "data"
    config = ServerConfig(data_dir=config_dir)

    registry = ToolRegistry([schema_dir])
    file_manager = FileManager(config)
    executor = ToolExecutor(config, registry, file_manager)

    job_id = executor.execute_tool("echo", {"output": "result.txt"})

    for _ in range(20):
        executor.poll_jobs()
        job = executor.get_job(job_id)
        if job and job.status != JobStatus.RUNNING:
            break
        time.sleep(0.1)

    job = executor.get_job(job_id)
    assert job is not None
    assert job.status == JobStatus.COMPLETED
    assert "result.txt" in job.output_files
    assert (job.work_dir / "result.txt").read_text() == "done"
    assert job.output_artifacts
    stored = Path(job.output_artifacts[0]["path"])
    assert stored.exists()


def test_executor_builds_shell_command(tmp_path: Path, schema_dir: Path) -> None:
    payload = {
        "id": "shell_tool",
        "runner": "shell",
        "entry_point": "echo hello",
        "parameters": [],
    }
    with (schema_dir / "shell_tool.json").open("w", encoding="utf-8") as handle:
        json.dump(payload, handle)

    config = ServerConfig(data_dir=tmp_path / "data")
    registry = ToolRegistry([schema_dir])
    file_manager = FileManager(config)
    executor = ToolExecutor(config, registry, file_manager)

    command = executor._build_command("shell_tool", tmp_path / "args.json")
    assert command[:2] == ["/bin/sh", "-c"]


def test_executor_resolves_file_parameters(tmp_path: Path, schema_dir: Path, tool_script: Path) -> None:
    schema_payload = {
        "id": "file_tool",
        "runner": "python",
        "entry_point": str(tool_script),
        "parameters": [
            {"name": "Input file", "parameter": "input", "type": "file", "multiple": True},
            {"name": "Output file", "parameter": "output", "type": "output", "default": "result.txt"},
        ],
    }
    with (schema_dir / "file_tool.json").open("w", encoding="utf-8") as handle:
        json.dump(schema_payload, handle)

    config = ServerConfig(data_dir=tmp_path / "data")
    registry = ToolRegistry([schema_dir])
    file_manager = FileManager(config)
    executor = ToolExecutor(config, registry, file_manager)

    source_file = tmp_path / "source.gb"
    source_file.write_text("CONTENT")
    artifact = file_manager.ingest(source_file)

    job_id = executor.execute_tool("file_tool", {"input": [artifact.file_id], "output": "output.gb"})
    job = executor.get_job(job_id)
    assert job is not None
    assert job.parameters["input"] == [str(artifact.path)]
    assert (job.work_dir / "output.gb").parent == job.work_dir


def test_executor_accepts_display_name_payload(tmp_path: Path, schema_dir: Path, tool_script: Path) -> None:
    schema_payload = {
        "id": "alias_tool",
        "runner": "python",
        "entry_point": str(tool_script),
        "parameters": [
            {"name": "Reference databases", "parameter": "references", "type": "file", "multiple": True},
            {"name": "CPUs", "parameter": "cpu", "type": "integer"},
        ],
    }
    with (schema_dir / "alias_tool.json").open("w", encoding="utf-8") as handle:
        json.dump(schema_payload, handle)

    config = ServerConfig(data_dir=tmp_path / "data")
    registry = ToolRegistry([schema_dir])
    file_manager = FileManager(config)
    executor = ToolExecutor(config, registry, file_manager)

    source_file = tmp_path / "ref.hmm"
    source_file.write_text("HMM")
    artifact = file_manager.ingest(source_file)

    job_id = executor.execute_tool(
        "alias_tool",
        {
            "Reference databases": [artifact.file_id],
            "CPUs": "2",
        },
    )
    job = executor.get_job(job_id)
    assert job is not None
    assert job.parameters["references"] == [str(artifact.path)]
    assert job.parameters["cpu"] == 2


def test_executor_sanitizes_output_filename(tmp_path: Path, schema_dir: Path) -> None:
    script = tmp_path / "writer.py"
    script.write_text(SCRIPT_TEMPLATE)
    script.chmod(script.stat().st_mode | stat.S_IEXEC)

    schema_payload = {
        "id": "sanitized",
        "runner": "python",
        "entry_point": str(script),
        "parameters": [
            {"name": "Output file", "parameter": "output", "type": "output", "default": "result.txt"},
        ],
    }
    with (schema_dir / "sanitized.json").open("w", encoding="utf-8") as handle:
        json.dump(schema_payload, handle)

    config = ServerConfig(data_dir=tmp_path / "data")
    registry = ToolRegistry([schema_dir])
    file_manager = FileManager(config)
    executor = ToolExecutor(config, registry, file_manager)

    requested_name = "../unsafe/output?.gbk"
    job_id = executor.execute_tool("sanitized", {"output": requested_name})

    for _ in range(20):
        executor.poll_jobs()
        job = executor.get_job(job_id)
        if job and job.status != JobStatus.RUNNING:
            break
        time.sleep(0.1)

    job = executor.get_job(job_id)
    assert job is not None
    output_path = Path(job.parameters["output"])
    assert output_path.parent == job.work_dir
    assert secure_filename(output_path.name) == output_path.name
    assert output_path.is_file()
    assert ".." not in output_path.parts


def test_executor_loads_existing_manifests(tmp_path: Path, schema_dir: Path, tool_script: Path) -> None:
    _write_schema(schema_dir, "echo", tool_script)
    config = ServerConfig(data_dir=tmp_path / "data")

    registry = ToolRegistry([schema_dir])
    file_manager = FileManager(config)
    executor = ToolExecutor(config, registry, file_manager)

    job_id = executor.execute_tool("echo", {"output": "result.txt"})

    for _ in range(20):
        executor.poll_jobs()
        job = executor.get_job(job_id)
        if job and job.status != JobStatus.RUNNING:
            break
        time.sleep(0.1)

    job = executor.get_job(job_id)
    assert job is not None
    assert job.status == JobStatus.COMPLETED
    assert job.output_artifacts

    # Recreate executor to force manifest reload
    registry_reloaded = ToolRegistry([schema_dir])
    file_manager_reloaded = FileManager(config)
    executor_reloaded = ToolExecutor(config, registry_reloaded, file_manager_reloaded)

    restored = executor_reloaded.get_job(job_id)
    assert restored is not None
    assert restored.status == JobStatus.COMPLETED
    assert restored.output_files == job.output_files
    assert restored.output_artifacts == job.output_artifacts
    assert restored.config_path is not None
    assert restored.config_path.exists()
