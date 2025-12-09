from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from shutil import copy2
from typing import List, Optional
from uuid import uuid4

from werkzeug.datastructures import FileStorage
from werkzeug.utils import secure_filename

from .config import ServerConfig

_METADATA_FILENAME = "metadata.json"
_PAYLOAD_SUBDIR = "files"


@dataclass(slots=True)
class Artifact:
    file_id: str
    original_name: str
    path: Path
    metadata_path: Path
    size: int
    file_type: Optional[str]
    uploaded_at: datetime
    checksum: Optional[str] = None
    description: Optional[str] = None

    def to_payload(self) -> dict:
        return {
            "file_id": self.file_id,
            "original_name": self.original_name,
            "path": str(self.path),
            "metadata_path": str(self.metadata_path),
            "size": self.size,
            "type": self.file_type,
            "uploaded_at": self.uploaded_at.isoformat(),
            "checksum": self.checksum,
            "description": self.description,
        }

    @classmethod
    def from_metadata(cls, metadata_path: Path) -> Artifact:
        with metadata_path.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
        file_id = data["file_id"]
        original = data["original_name"]
        path = Path(data["path"])
        uploaded_at = datetime.fromisoformat(data["uploaded_at"])
        return cls(
            file_id=file_id,
            original_name=original,
            path=path,
            metadata_path=metadata_path,
            size=data.get("size", 0),
            file_type=data.get("type"),
            uploaded_at=uploaded_at,
            checksum=data.get("checksum"),
            description=data.get("description"),
        )


class FileManager:
    def __init__(self, config: ServerConfig):
        self.config = config

    def ingest(self, source: Path, file_type: Optional[str] = None, description: Optional[str] = None) -> Artifact:
        if not source.exists():
            raise FileNotFoundError(source)
        file_id = uuid4().hex
        dest_dir = self.config.paths.uploads_dir / file_id
        files_dir = dest_dir / _PAYLOAD_SUBDIR
        files_dir.mkdir(parents=True, exist_ok=True)
        sanitized_name = self._sanitize_filename(source.name)
        dest_path = files_dir / sanitized_name
        copy2(source, dest_path)
        checksum = self._sha256(dest_path)
        resolved_type = file_type or self._guess_file_type(sanitized_name)
        metadata = Artifact(
            file_id=file_id,
            original_name=sanitized_name,
            path=dest_path,
            metadata_path=dest_dir / _METADATA_FILENAME,
            size=dest_path.stat().st_size,
            file_type=resolved_type,
            uploaded_at=datetime.now(timezone.utc),
            checksum=checksum,
            description=description,
        )
        self._write_metadata(metadata)
        return metadata

    def ingest_filestorage(
        self,
        file_storage: FileStorage,
        file_type: Optional[str] = None,
        description: Optional[str] = None,
    ) -> Artifact:
        filename = self._sanitize_filename(file_storage.filename or "")

        file_id = uuid4().hex
        dest_dir = self.config.paths.uploads_dir / file_id
        files_dir = dest_dir / _PAYLOAD_SUBDIR
        files_dir.mkdir(parents=True, exist_ok=True)
        dest_path = files_dir / filename

        file_storage.save(dest_path)
        checksum = self._sha256(dest_path)
        resolved_type = file_type or self._guess_file_type(filename)

        metadata = Artifact(
            file_id=file_id,
            original_name=filename,
            path=dest_path,
            metadata_path=dest_dir / _METADATA_FILENAME,
            size=dest_path.stat().st_size,
            file_type=resolved_type,
            uploaded_at=datetime.now(timezone.utc),
            checksum=checksum,
            description=description,
        )
        self._write_metadata(metadata)
        return metadata

    def list_uploads(self) -> List[Artifact]:
        artifacts: List[Artifact] = []
        for entry in self.config.paths.uploads_dir.iterdir():
            meta_path = entry / _METADATA_FILENAME
            if not meta_path.exists():
                continue
            try:
                artifacts.append(Artifact.from_metadata(meta_path))
            except Exception:
                continue
        return artifacts

    def get_upload(self, file_id: str) -> Artifact:
        meta_path = self.config.paths.uploads_dir / file_id / _METADATA_FILENAME
        if not meta_path.exists():
            raise FileNotFoundError(file_id)
        return Artifact.from_metadata(meta_path)

    def delete_upload(self, file_id: str) -> None:
        target = self.config.paths.uploads_dir / file_id
        if not target.exists():
            return
        for item in target.iterdir():
            if item.is_file():
                item.unlink(missing_ok=True)
            else:
                self._remove_tree(item)
        target.rmdir()

    def rename_upload(self, file_id: str, new_name: str) -> Artifact:
        artifact = self.get_upload(file_id)
        sanitized = self._sanitize_filename(new_name)
        if not sanitized:
            raise ValueError("invalid filename")

        current_path = artifact.path
        dest_path = current_path.parent / sanitized

        if dest_path.exists() and dest_path != current_path:
            raise FileExistsError(dest_path.name)

        if dest_path != current_path:
            current_path.rename(dest_path)

        updated_type = self._guess_file_type(sanitized) or artifact.file_type
        updated = Artifact(
            file_id=artifact.file_id,
            original_name=sanitized,
            path=dest_path,
            metadata_path=artifact.metadata_path,
            size=dest_path.stat().st_size if dest_path.exists() else artifact.size,
            file_type=updated_type,
            uploaded_at=artifact.uploaded_at,
            checksum=artifact.checksum,
            description=artifact.description,
        )

        self._write_metadata(updated)
        return updated

    def _write_metadata(self, artifact: Artifact) -> None:
        payload = artifact.to_payload()
        with artifact.metadata_path.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)

    def _remove_tree(self, path: Path) -> None:
        for entry in path.iterdir():
            if entry.is_dir():
                self._remove_tree(entry)
            else:
                entry.unlink(missing_ok=True)
        path.rmdir()

    def _sha256(self, path: Path) -> str:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return f"sha256:{digest.hexdigest()}"

    def _guess_file_type(self, filename: str) -> Optional[str]:
        suffix = Path(filename).suffix.lower()
        mapping = {
            ".gb": "genbank",
            ".gbk": "genbank",
            ".fasta": "fasta",
            ".fa": "fasta",
            ".hmm": "hmm",
            ".hdf5": "hdf5",
            ".tsv": "tsv",
            ".html": "html",
        }
        return mapping.get(suffix) or (suffix[1:] if suffix else None)

    def _sanitize_filename(self, name: str | None) -> str:
        original = name or ""
        sanitized = secure_filename(original)
        if sanitized:
            return sanitized
        suffix = Path(original).suffix if original else ""
        fallback = secure_filename(f"file{suffix}")
        if fallback:
            return fallback
        return "file"
