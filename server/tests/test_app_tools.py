from __future__ import annotations

from pathlib import Path

import pytest

from domainator_server.app import create_app
from domainator_server.config import ServerConfig


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
    assert "genome_annotation" in ids

    detail = client.get("/api/tools/genome_annotation")
    assert detail.status_code == 200
    assert detail.get_json()["id"] == "genome_annotation"
