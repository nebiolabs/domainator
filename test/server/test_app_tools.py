from __future__ import annotations

from pathlib import Path

import pytest

from unittest.mock import MagicMock

from domainator.server.app import create_app
from domainator.server.config import ServerConfig
from domainator.server.tool_executor import Job, JobStatus


@pytest.fixture()
def flask_app(tmp_path: Path):
    config = ServerConfig(data_dir=tmp_path / "data")
    app = create_app(config, [])
    app.config.update(TESTING=True)
    return app


@pytest.fixture()
def client(flask_app):
    return flask_app.test_client()


def test_list_tools(client):
    resp = client.get("/api/tools")
    assert resp.status_code == 200
    tools = resp.get_json()
    assert isinstance(tools, list)
    ids = {tool["id"] for tool in tools}
    assert "domainate" in ids

    detail = client.get("/api/tools/domainate")
    assert detail.status_code == 200
    assert detail.get_json()["id"] == "domainate"


def test_cancel_job_endpoint(client, flask_app, tmp_path):
    executor = flask_app.extensions["domainator_server"]["executor"]

    job_id = "test-cancel"
    work_dir = executor.config.paths.outputs_dir / job_id
    work_dir.mkdir(parents=True, exist_ok=True)

    process = MagicMock()
    process.terminate.return_value = None
    process.wait.return_value = None
    process.kill.return_value = None
    process.poll.return_value = None

    job = Job(
        job_id=job_id,
        tool_name="dummy",
        parameters={},
        status=JobStatus.RUNNING,
        process=process,
        work_dir=work_dir,
    )
    executor.jobs[job_id] = job

    resp = client.post(f"/api/jobs/{job_id}/cancel")
    assert resp.status_code == 200
    payload = resp.get_json()
    assert payload["job_id"] == job_id
    assert payload["status"] == JobStatus.CANCELLED.value
