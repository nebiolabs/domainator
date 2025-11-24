from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence
import sys

DEFAULT_DATA_DIR = Path.home() / "domainator_server"


@dataclass(slots=True)
class ServerPaths:
    base_dir: Path
    uploads_dir: Path
    outputs_dir: Path
    jobs_dir: Path
    logs_dir: Path


@dataclass(slots=True)
class ServerConfig:
    data_dir: Optional[Path] = None
    debug: bool = False
    python_executable: Optional[str] = None
    schema_dirs: Optional[Sequence[Path]] = None
    _paths: ServerPaths = field(init=False, repr=False)

    def __post_init__(self) -> None:
        base_dir = (self.data_dir or DEFAULT_DATA_DIR).expanduser().resolve()
        uploads_dir = base_dir / "uploads"
        outputs_dir = base_dir / "outputs"
        jobs_dir = base_dir / "jobs"
        logs_dir = base_dir / "logs"
        for directory in (base_dir, uploads_dir, outputs_dir, jobs_dir, logs_dir):
            directory.mkdir(parents=True, exist_ok=True)
        self._paths = ServerPaths(
            base_dir=base_dir,
            uploads_dir=uploads_dir,
            outputs_dir=outputs_dir,
            jobs_dir=jobs_dir,
            logs_dir=logs_dir,
        )
        if self.python_executable is None:
            self.python_executable = self._infer_python()
        if self.schema_dirs is None:
            base_schema = Path(__file__).resolve().parent / "schemas"
            self.schema_dirs = [base_schema / "sample_tools", base_schema / "sample_workflows"]

    @property
    def paths(self) -> ServerPaths:
        return self._paths

    def _infer_python(self) -> str:
        return sys.executable
