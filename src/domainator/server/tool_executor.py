from __future__ import annotations

import json
import math
import re
import subprocess
import time
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Dict, Iterable, List, Optional, TextIO, Union
from uuid import uuid4

from .config import ServerConfig
from .file_manager import FileManager
from werkzeug.utils import secure_filename

from .tool_registry import ToolRegistry


class JobStatus(str, Enum):
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


@dataclass(slots=True)
class Job:
    job_id: str
    tool_name: str
    parameters: Dict[str, object]
    status: JobStatus
    process: Optional[subprocess.Popen] = None
    started_at: Optional[float] = None
    completed_at: Optional[float] = None
    output_files: List[str] = field(default_factory=list)
    output_artifacts: List[Dict[str, object]] = field(default_factory=list)
    error: Optional[str] = None
    log_path: Optional[Path] = None
    work_dir: Optional[Path] = None
    updated_at: Optional[float] = None
    config_path: Optional[Path] = None
    _log_handle: Optional[TextIO] = field(default=None, repr=False)


class ToolExecutor:
    def __init__(self, config: ServerConfig, registry: ToolRegistry, file_manager: FileManager):
        self.config = config
        self.registry = registry
        self.file_manager = file_manager
        self.jobs: Dict[str, Job] = {}
        self._load_existing_jobs()

    def execute_tool(self, tool_name: str, parameters: Dict[str, object]) -> str:
        job_id = self._generate_job_id()
        work_dir = self._prepare_job_dir(job_id)
        log_path = self.config.paths.logs_dir / f"{job_id}.log"
        log_path.parent.mkdir(parents=True, exist_ok=True)

        job = Job(
            job_id=job_id,
            tool_name=tool_name,
            parameters={},
            status=JobStatus.QUEUED,
            log_path=log_path,
            work_dir=work_dir,
        )
        resolved_parameters = self._prepare_parameters(tool_name, parameters, work_dir)
        job.parameters = resolved_parameters
        now = time.time()
        job.updated_at = now
        self.jobs[job_id] = job
        cfg_path = self._write_args_file(resolved_parameters, work_dir)
        job.parameters = resolved_parameters
        job.config_path = cfg_path
        self._write_job_manifest(job)
        command = self._build_command(tool_name, cfg_path)

        log_handle = log_path.open("w", encoding="utf-8")
        process = subprocess.Popen(
            command,
            cwd=work_dir,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            text=True,
        )

        job.process = process
        job.status = JobStatus.RUNNING
        job.started_at = time.time()
        job.updated_at = job.started_at
        job._log_handle = log_handle
        self._write_job_manifest(job)
        return job.job_id

    def get_job(self, job_id: str) -> Optional[Job]:
        return self.jobs.get(job_id)

    def all_jobs(self) -> List[Job]:
        return list(self.jobs.values())

    def updates_since(self, timestamp: float) -> List[Job]:
        return [job for job in self.jobs.values() if job.updated_at and job.updated_at > timestamp]

    def poll_jobs(self) -> None:
        for job in self.jobs.values():
            if job.status != JobStatus.RUNNING or job.process is None:
                continue
            return_code = job.process.poll()
            if return_code is None:
                continue
            job.completed_at = time.time()
            if return_code == 0:
                job.status = JobStatus.COMPLETED
                job.output_files = self._discover_outputs(job)
                if not job.output_artifacts:
                    job.output_artifacts = self._ingest_outputs(job)
            else:
                job.status = JobStatus.FAILED
                job.error = f"Process exited with code {return_code}"
            if job._log_handle:
                job._log_handle.close()
                job._log_handle = None
            job.updated_at = time.time()
            self._write_job_manifest(job)

    def cancel_job(self, job_id: str) -> None:
        job = self.jobs.get(job_id)
        if not job or job.status != JobStatus.RUNNING or job.process is None:
            return
        job.process.terminate()
        try:
            job.process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            job.process.kill()
            job.process.wait()
        job.status = JobStatus.CANCELLED
        job.completed_at = time.time()
        if job._log_handle:
            job._log_handle.close()
            job._log_handle = None
        job.updated_at = time.time()
        self._write_job_manifest(job)

    def _generate_job_id(self) -> str:
        return uuid4().hex

    def _prepare_job_dir(self, job_id: str) -> Path:
        work_dir = self.config.paths.outputs_dir / job_id
        work_dir.mkdir(parents=True, exist_ok=True)
        return work_dir

    def _write_args_file(self, parameters: Dict[str, object], work_dir: Path) -> Path:
        with NamedTemporaryFile("w", delete=False, suffix=".json", dir=work_dir) as handle:
            json.dump(parameters, handle, indent=2)
            return Path(handle.name)

    def _build_command(self, tool_name: str, cfg_path: Path) -> List[str]:
        schema = self._schema_for(tool_name)
        runner = schema.get("runner", "binary")
        entry = schema["entry_point"]
        config_flag = f"--config={cfg_path}"
        if runner == "python":
            return [self.config.python_executable, entry, config_flag]
        if runner == "module":
            return [self.config.python_executable, "-m", entry, config_flag]
        if runner == "shell":
            return ["/bin/sh", "-c", f"{entry} {config_flag}"]
        return [entry, config_flag]

    def _schema_for(self, tool_name: str) -> Dict[str, object]:
        schema = self.registry.get(tool_name)
        if schema is None:
            raise KeyError(f"Unknown tool id: {tool_name}")
        return schema

    def _prepare_parameters(self, tool_name: str, parameters: Dict[str, object], work_dir: Path) -> Dict[str, object]:
        schema = self._schema_for(tool_name)
        definitions: Dict[str, dict] = {}
        for entry in schema.get("parameters", []):
            if not isinstance(entry, dict):
                continue
            display_name = entry.get("name")
            actual_key = entry.get("parameter") or display_name
            if actual_key:
                definitions[actual_key] = entry
            if display_name and display_name != actual_key:
                definitions.setdefault(display_name, entry)

        resolved: Dict[str, object] = {}
        for key, value in parameters.items():
            definition = definitions.get(key)
            if not definition:
                resolved[key] = value
                continue
            param_type = str(definition.get("type", "")).lower()
            multiple = bool(definition.get("multiple"))
            behavior = definition.get("behavior")
            target_key = (
                definition.get("target_parameter")
                or (definition.get("append_const") or {}).get("target")
                or (definition.get("dynamic_action") or {}).get("target")
                or definition.get("parameter")
                or definition.get("name")
                or key
            )

            if behavior == "append_const":
                if self._ensure_boolean(value, False):
                    self._append_to_list(resolved, target_key, [definition.get("append_const", {}).get("value")])
                continue

            if behavior == "dynamic":
                groups = self._coerce_dynamic_values(definition, value)
                self._append_to_list(resolved, target_key, groups)
                continue

            if param_type == "file":
                resolved[target_key] = self._resolve_file_value(value, multiple)
            elif param_type == "output":
                resolved[target_key] = str(self._resolve_output_value(target_key, value, work_dir, definition))
            elif param_type in {"integer", "number", "float"}:
                resolved[target_key] = self._convert_numeric_value(param_type, value, multiple)
            elif param_type == "boolean":
                resolved[target_key] = self._ensure_boolean(value, multiple)
            else:
                resolved[target_key] = self._ensure_sequence(value, multiple)
        return resolved

    def _append_to_list(self, resolved: Dict[str, object], key: str, values: Iterable[object]) -> None:
        items = [item for item in values if item is not None]
        if not items:
            return
        existing = resolved.get(key)
        if existing is None:
            resolved[key] = list(items)
            return
        if not isinstance(existing, list):
            existing = [existing]
        existing.extend(items)
        resolved[key] = existing

    def _ensure_sequence(self, value: object, multiple: bool) -> object:
        if multiple:
            values = self._as_iterable(value)
            return [self._clean_scalar(item) for item in values if item is not None and item != ""]
        return self._clean_scalar(value)

    def _clean_scalar(self, value: object) -> object:
        if isinstance(value, str):
            return value.strip()
        return value

    def _convert_numeric_value(self, param_type: str, value: object, multiple: bool) -> object:
        def _convert_single(item: object) -> Union[int, float]:
            if item is None or item == "":
                raise ValueError("numeric parameter requires a value")
            number = float(item)
            if not math.isfinite(number):
                raise ValueError("numeric parameter must be finite")
            if param_type == "integer":
                if not float(number).is_integer():
                    raise ValueError("integer parameter must be a whole number")
                return int(number)
            return number

        if multiple:
            return [_convert_single(item) for item in self._as_iterable(value)]
        return _convert_single(value)

    def _ensure_boolean(self, value: object, multiple: bool) -> object:
        def _convert(item: object) -> bool:
            if isinstance(item, bool):
                return item
            if isinstance(item, str):
                lowered = item.strip().lower()
                if lowered in {"true", "1", "yes", "on"}:
                    return True
                if lowered in {"false", "0", "no", "off"}:
                    return False
            return bool(item)

        if multiple:
            return [_convert(entry) for entry in self._as_iterable(value)]
        return _convert(value)

    def _coerce_dynamic_values(self, definition: dict, value: object) -> List[List[object]]:
        meta = definition.get("dynamic_action") or {}
        const_value = meta.get("const")
        nargs = meta.get("nargs")
        if const_value is None:
            return []

        groups: List[List[object]] = []
        raw_groups: List[List[object]]

        if value is None or value == "":
            return []
        if isinstance(value, str):
            lines = [line.strip() for line in value.splitlines() if line.strip()]
            raw_groups = [re.split(r"[\s,]+", line) for line in lines]
        elif isinstance(value, (list, tuple)):
            if value and all(isinstance(item, (list, tuple)) for item in value):
                raw_groups = [list(item) for item in value]
            else:
                raw_groups = [list(value)]
        else:
            raw_groups = [[value]]

        for group in raw_groups:
            cleaned = [self._clean_scalar(item) for item in group if item is not None and str(item).strip() != ""]
            if not cleaned:
                continue
            if not self._validate_dynamic_group(cleaned, nargs):
                requirement = self._describe_dynamic_requirement(nargs)
                name = definition.get("name") or definition.get("parameter") or "dynamic parameter"
                raise ValueError(f"{name} expects {requirement}")
            groups.append([const_value, cleaned])
        return groups

    @staticmethod
    def _validate_dynamic_group(group: List[object], nargs: object) -> bool:
        if not group:
            return False
        if nargs in (None, "*"):
            return True
        if nargs == "+":
            return len(group) >= 1
        if nargs == "?":
            return len(group) <= 1
        try:
            count = int(nargs)
        except (TypeError, ValueError):
            return True
        if count <= 0:
            return True
        return len(group) == count

    @staticmethod
    def _describe_dynamic_requirement(nargs: object) -> str:
        if nargs in (None, "*"):
            return "any number of values"
        if nargs == "+":
            return "at least one value"
        if nargs == "?":
            return "zero or one value"
        try:
            count = int(nargs)
        except (TypeError, ValueError):
            return "valid values"
        if count <= 0:
            return "valid values"
        return "exactly one value" if count == 1 else f"exactly {count} values"

    def _resolve_file_value(self, value: object, multiple: bool) -> Union[str, List[str]]:
        def _resolve_single(item: object) -> str:
            if item is None or item == "":
                raise ValueError("file parameter requires a value")
            if isinstance(item, str):
                path_candidate = Path(item)
                if path_candidate.exists():
                    return str(path_candidate.resolve())
                try:
                    artifact = self.file_manager.get_upload(item)
                except FileNotFoundError as exc:
                    raise ValueError(f"unknown file id '{item}'") from exc
                return str(artifact.path)
            raise ValueError("unsupported file reference")

        if multiple:
            return [_resolve_single(item) for item in self._as_iterable(value)]
        return _resolve_single(value)

    def _resolve_output_value(self, name: str, value: object, work_dir: Path, definition: dict) -> Path:
        default_name = definition.get("default") or f"{name}.out"
        raw_name = str(value).strip() if isinstance(value, str) and value.strip() else default_name
        safe_name = secure_filename(raw_name)
        if not safe_name:
            suffix = Path(raw_name).suffix if raw_name else ""
            fallback = f"{name}{suffix}" if suffix else f"{name}.out"
            safe_name = secure_filename(fallback) or fallback

        candidate = Path(safe_name)
        if candidate.is_absolute():
            target = candidate
        else:
            target = (work_dir / candidate.name).resolve()

        work_dir_resolved = work_dir.resolve()
        if work_dir_resolved not in target.parents and target != work_dir_resolved:
            raise ValueError("output path must reside within the job workspace")

        target.parent.mkdir(parents=True, exist_ok=True)
        return target

    def _as_iterable(self, value: object) -> Iterable[object]:
        if value is None:
            return []
        if isinstance(value, (list, tuple)):
            return value
        return [value]

    def _discover_outputs(self, job: Job) -> List[str]:
        if not job.work_dir:
            return []
        outputs: List[str] = []
        for path in job.work_dir.glob("**/*"):
            if path.is_dir():
                continue
            if path.name.endswith(".json") and path.name.startswith("tmp"):
                continue
            outputs.append(str(path.relative_to(job.work_dir)))
        return outputs

    def _ingest_outputs(self, job: Job) -> List[Dict[str, object]]:
        artifacts: List[Dict[str, object]] = []
        if not job.work_dir:
            return artifacts
        for relative in job.output_files:
            abs_path = (job.work_dir / relative).resolve()
            if not abs_path.exists() or not abs_path.is_file():
                continue
            try:
                artifact = self.file_manager.ingest(
                    abs_path,
                    description=f"Output from job {job.job_id} ({job.tool_name})",
                )
            except Exception:
                continue
            file_type = artifact.file_type or (abs_path.suffix.lower()[1:] if abs_path.suffix else None)
            artifacts.append(
                {
                    "file_id": artifact.file_id,
                    "name": artifact.original_name,
                    "path": str(artifact.path),
                    "size": artifact.size,
                    "relative_path": relative,
                    "type": file_type,
                }
            )
        return artifacts

    def _write_job_manifest(self, job: Job) -> None:
        if not job.work_dir:
            raise ValueError("Job workspace not initialised")
        manifest_dir = self.config.paths.jobs_dir / job.job_id
        manifest_dir.mkdir(parents=True, exist_ok=True)
        payload = {
            "job_id": job.job_id,
            "tool": job.tool_name,
            "status": job.status.value,
            "parameters": job.parameters,
            "started_at": job.started_at,
            "completed_at": job.completed_at,
            "output_files": job.output_files,
            "output_artifacts": job.output_artifacts,
            "error": job.error,
            "log_file": str(job.log_path) if job.log_path else None,
            "work_dir": str(job.work_dir) if job.work_dir else None,
            "updated_at": job.updated_at,
            "config_path": str(job.config_path) if job.config_path else None,
        }
        manifest_path = manifest_dir / "job.json"
        with manifest_path.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)

    def _load_existing_jobs(self) -> None:
        jobs_dir = self.config.paths.jobs_dir
        if not jobs_dir.exists():
            return
        for manifest_path in jobs_dir.glob("*/job.json"):
            try:
                with manifest_path.open("r", encoding="utf-8") as handle:
                    data = json.load(handle)
            except Exception:
                continue

            job_id = data.get("job_id")
            if not job_id or job_id in self.jobs:
                continue

            status_value = data.get("status", JobStatus.FAILED.value)
            try:
                status = JobStatus(status_value)
            except ValueError:
                status = JobStatus.FAILED

            job = Job(
                job_id=job_id,
                tool_name=data.get("tool", ""),
                parameters=data.get("parameters") or {},
                status=status,
                started_at=data.get("started_at"),
                completed_at=data.get("completed_at"),
                output_files=data.get("output_files") or [],
                output_artifacts=data.get("output_artifacts") or [],
                error=data.get("error"),
                log_path=Path(data["log_file"]) if data.get("log_file") else None,
                work_dir=Path(data["work_dir"]) if data.get("work_dir") else None,
                updated_at=data.get("updated_at"),
                config_path=Path(data["config_path"]) if data.get("config_path") else None,
            )

            if job.work_dir is None:
                job.work_dir = self.config.paths.outputs_dir / job_id

            self.jobs[job_id] = job

            if job.status in {JobStatus.RUNNING, JobStatus.QUEUED}:
                job.status = JobStatus.FAILED
                if not job.error:
                    job.error = "Job interrupted by server restart"
                job.completed_at = job.completed_at or time.time()
                job.updated_at = time.time()
                self._write_job_manifest(job)