from __future__ import annotations

from io import BytesIO
from pathlib import Path

from werkzeug.datastructures import FileStorage

from domainator_server.config import ServerConfig
from domainator_server.file_manager import FileManager


def test_ingest_filestorage(tmp_path: Path) -> None:
    config = ServerConfig(data_dir=tmp_path / "data")
    manager = FileManager(config)

    storage = FileStorage(stream=BytesIO(b"contents"), filename="sample.gb")
    artifact = manager.ingest_filestorage(storage)

    assert artifact.file_type == "genbank"
    assert artifact.path.exists()
    listed = manager.list_uploads()
    assert any(item.file_id == artifact.file_id for item in listed)

    fetched = manager.get_upload(artifact.file_id)
    assert fetched.file_id == artifact.file_id

    manager.delete_upload(artifact.file_id)
    assert not (config.paths.uploads_dir / artifact.file_id).exists()


def test_ingest_preserves_reserved_metadata_filename(tmp_path: Path) -> None:
    config = ServerConfig(data_dir=tmp_path / "data")
    manager = FileManager(config)

    storage = FileStorage(stream=BytesIO(b"payload"), filename="metadata.json")
    artifact = manager.ingest_filestorage(storage)

    assert artifact.path.name == "metadata.json"
    assert artifact.path.parent.name == "files"

    payload = artifact.path.read_bytes()
    assert payload == b"payload"

    metadata_path = config.paths.uploads_dir / artifact.file_id / "metadata.json"
    assert metadata_path.exists()
    assert metadata_path != artifact.path

    loaded = manager.get_upload(artifact.file_id)
    assert loaded.metadata_path == metadata_path
    assert loaded.path == artifact.path


def test_ingest_sanitizes_filenames(tmp_path: Path) -> None:
    config = ServerConfig(data_dir=tmp_path / "data")
    manager = FileManager(config)

    odd_name = "  weird file name ?!.gbk  "
    source = tmp_path / odd_name
    source.write_text("hello")

    artifact = manager.ingest(source)
    safe_name = manager._sanitize_filename(odd_name)

    assert artifact.path.name == safe_name
    assert artifact.original_name == safe_name
    assert " " not in artifact.path.name
    assert "?" not in artifact.path.name

    storage = FileStorage(stream=BytesIO(b"payload"), filename="../tricky?.hmm")
    uploaded = manager.ingest_filestorage(storage)
    safe_upload = manager._sanitize_filename("../tricky?.hmm")
    assert uploaded.path.name == safe_upload
    assert uploaded.original_name == safe_upload
