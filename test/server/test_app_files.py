from __future__ import annotations

from io import BytesIO
from pathlib import Path

import pytest

from domainator.server.app import create_app
from domainator.server.config import ServerConfig


@pytest.fixture()
def flask_app(tmp_path: Path):
    config = ServerConfig(data_dir=tmp_path / "data")
    app = create_app(config, [])
    app.config.update(TESTING=True)
    return app


@pytest.fixture()
def client(flask_app):
    return flask_app.test_client()


def test_upload_list_download_delete(client):
    response = client.post(
        "/api/files/upload",
        data={"files": (BytesIO(b"hello"), "example.gb")},
        content_type="multipart/form-data",
    )
    assert response.status_code == 201
    payload = response.get_json()
    assert payload and "files" in payload
    file_id = payload["files"][0]["file_id"]

    list_response = client.get("/api/files")
    assert list_response.status_code == 200
    assert any(entry["file_id"] == file_id for entry in list_response.get_json())

    metadata = client.get(f"/api/files/{file_id}")
    assert metadata.status_code == 200
    assert metadata.get_json()["file_id"] == file_id

    download = client.get(f"/api/files/{file_id}/download")
    assert download.status_code == 200
    assert download.data == b"hello"

    delete = client.delete(f"/api/files/{file_id}")
    assert delete.status_code == 204

    missing = client.get(f"/api/files/{file_id}")
    assert missing.status_code == 404
