from __future__ import annotations

import json
import stat
import time
from pathlib import Path

import pytest

from domainator_server.app import create_app
from domainator_server.config import ServerConfig

SCRIPT = """#!/usr/bin/env python3
import argparse
import json
from pathlib import Path
parser = argparse.ArgumentParser()
parser.add_argument('--config', required=True)
args = parser.parse_args()
with open(args.config, 'r', encoding='utf-8') as handle:
    data = json.load(handle)
output = data.get('output', 'out.txt')
Path(output).write_text('ok')
"""


@pytest.fixture()
def job_client(tmp_path: Path):
    schema_dir = tmp_path / "schemas"
    schema_dir.mkdir()
    script_path = tmp_path / "tool.py"
    script_path.write_text(SCRIPT)
    script_path.chmod(script_path.stat().st_mode | stat.S_IEXEC)

    schema = {
        "id": "echo_tool",
        "runner": "python",
        "entry_point": str(script_path),
        "parameters": [
            {"name": "output", "type": "output", "default": "out.txt"}
        ]
    }
    with (schema_dir / "echo_tool.json").open("w", encoding="utf-8") as handle:
        json.dump(schema, handle)

    config = ServerConfig(data_dir=tmp_path / "data", schema_dirs=[schema_dir])
    app = create_app(config)
    app.config.update(TESTING=True)
    return app.test_client()


def wait_for_completion(client, job_id: str, timeout: float = 5.0):
    start = time.time()
    while time.time() - start < timeout:
        resp = client.get(f"/api/jobs/{job_id}")
        assert resp.status_code == 200
        data = resp.get_json()
        if data["status"] != "running":
            return data
        time.sleep(0.1)
    raise TimeoutError("job did not complete")


def test_job_updates_endpoint(job_client):
    client = job_client

    initial = client.get("/api/jobs/updates")
    assert initial.status_code == 200
    assert initial.get_json() == []

    resp = client.post("/api/tools/echo_tool/execute", json={"parameters": {"output": "result.txt"}})
    assert resp.status_code == 202
    job_id = resp.get_json()["job_id"]

    updates = client.get("/api/jobs/updates?since=0")
    assert updates.status_code == 200
    jobs = updates.get_json()
    assert any(job["job_id"] == job_id for job in jobs)
    first_updated = max(job["updated_at"] for job in jobs)

    final_state = wait_for_completion(client, job_id)
    second_updates = client.get(f"/api/jobs/updates?since={first_updated}")
    assert second_updates.status_code == 200
    payload = second_updates.get_json()
    assert any(job["job_id"] == job_id for job in payload)
    latest_updated = max(job["updated_at"] for job in payload)

    empty = client.get(f"/api/jobs/updates?since={latest_updated}")
    assert empty.status_code == 200
    assert empty.get_json() == []

    artifacts = final_state.get("output_artifacts", [])
    assert artifacts
    file_id = artifacts[0]["file_id"]

    file_listing = client.get(f"/api/files/{file_id}")
    assert file_listing.status_code == 200
    download = client.get(f"/api/files/{file_id}/download")
    assert download.status_code == 200
    assert download.data == b"ok"