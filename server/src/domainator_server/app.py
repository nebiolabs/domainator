from __future__ import annotations

from pathlib import Path
from typing import Iterable

from flask import Flask, jsonify, render_template, request, send_file

from .config import ServerConfig
from .file_manager import FileManager
from .tool_executor import Job, JobStatus, ToolExecutor
from .tool_registry import ToolRegistry


def create_app(config: ServerConfig, schema_dirs: Iterable[Path] | None = None) -> Flask:
    app = Flask(__name__)
    file_manager = FileManager(config)
    dirs = list(schema_dirs or config.schema_dirs or [])
    registry = ToolRegistry(dirs)
    executor = ToolExecutor(config, registry, file_manager)

    app.config["SERVER_CONFIG"] = config
    app.extensions["domainator_server"] = {
        "file_manager": file_manager,
        "registry": registry,
        "executor": executor,
    }

    @app.before_request
    def _poll_jobs() -> None:
        executor.poll_jobs()

    @app.get("/")
    def index():
        return render_template("index.html")

    @app.get("/api/health")
    def health_check():
        return jsonify({"status": "ok"})

    def error_response(message: str, status: int = 400):
        return jsonify({"error": message}), status

    @app.get("/api/files")
    def list_files():
        artifacts = file_manager.list_uploads()
        return jsonify([artifact.to_payload() for artifact in artifacts])

    @app.post("/api/files/upload")
    def upload_files():
        if not request.files:
            return error_response("no files uploaded", 400)
        uploaded = []
        for storage in request.files.getlist("files"):
            if not storage.filename:
                continue
            try:
                artifact = file_manager.ingest_filestorage(storage)
            except ValueError as exc:
                return error_response(str(exc), 400)
            uploaded.append(artifact.to_payload())
        if not uploaded:
            return error_response("no valid files provided", 400)
        return jsonify({"files": uploaded}), 201

    @app.get("/api/files/<file_id>")
    def get_file(file_id: str):
        try:
            artifact = file_manager.get_upload(file_id)
        except FileNotFoundError:
            return error_response("file not found", 404)
        return jsonify(artifact.to_payload())

    @app.get("/api/files/<file_id>/download")
    def download_file(file_id: str):
        try:
            artifact = file_manager.get_upload(file_id)
        except FileNotFoundError:
            return error_response("file not found", 404)
        if not artifact.path.exists():
            return error_response("file missing on disk", 410)
        return send_file(artifact.path, as_attachment=True, download_name=artifact.original_name)

    @app.delete("/api/files/<file_id>")
    def delete_file(file_id: str):
        try:
            file_manager.get_upload(file_id)
        except FileNotFoundError:
            return error_response("file not found", 404)
        file_manager.delete_upload(file_id)
        return ("", 204)

    @app.get("/api/tools")
    def list_tools():
        schemas = [schema.payload for schema in registry.list() if schema.category != "Workflows"]
        return jsonify(schemas)

    @app.get("/api/tools/<tool_id>")
    def get_tool(tool_id: str):
        schema = registry.get(tool_id)
        if not schema:
            return error_response("not found", 404)
        return jsonify(schema)

    @app.get("/api/workflows")
    def list_workflows():
        workflows = [schema.payload for schema in registry.list() if schema.category == "Workflows"]
        return jsonify(workflows)

    @app.post("/api/tools/<tool_id>/execute")
    def execute_tool(tool_id: str):
        payload = request.get_json(silent=True) or {}
        parameters = payload.get("parameters", {})
        if not isinstance(parameters, dict):
            return jsonify({"error": "parameters must be an object"}), 400
        try:
            job_id = executor.execute_tool(tool_id, parameters)
        except KeyError:
            return jsonify({"error": f"unknown tool '{tool_id}'"}), 404
        except ValueError as exc:
            return error_response(str(exc), 400)
        return jsonify({"job_id": job_id, "status": JobStatus.RUNNING.value}), 202

    @app.get("/api/jobs")
    def list_jobs():
        jobs = [serialize_job(job) for job in executor.all_jobs()]
        return jsonify(jobs)

    @app.get("/api/jobs/<job_id>")
    def get_job(job_id: str):
        job = executor.get_job(job_id)
        if not job:
            return jsonify({"error": "not found"}), 404
        return jsonify(serialize_job(job))

    @app.post("/api/jobs/<job_id>/cancel")
    def cancel_job(job_id: str):
        job = executor.get_job(job_id)
        if not job:
            return error_response("not found", 404)
        executor.cancel_job(job_id)
        return jsonify(serialize_job(job))

    @app.get("/api/jobs/<job_id>/logs")
    def job_logs(job_id: str):
        job = executor.get_job(job_id)
        if not job:
            return error_response("not found", 404)
        if not job.log_path or not Path(job.log_path).exists():
            return error_response("log not available", 404)
        return send_file(job.log_path, mimetype="text/plain")

    @app.get("/api/jobs/<job_id>/outputs/<path:filename>")
    def job_output(job_id: str, filename: str):
        job = executor.get_job(job_id)
        if not job or not job.work_dir:
            return error_response("not found", 404)
        work_dir = Path(job.work_dir).resolve()
        target = (work_dir / filename).resolve()
        if work_dir not in target.parents and target != work_dir:
            return error_response("output not found", 404)
        if not target.exists() or not target.is_file():
            return error_response("output not found", 404)
        return send_file(target, as_attachment=True, download_name=target.name)

    @app.get("/api/jobs/updates")
    def job_updates():
        since_param = request.args.get("since")
        if since_param is None:
            return jsonify([serialize_job(job) for job in executor.all_jobs()])
        try:
            since_value = float(since_param)
        except ValueError:
            return error_response("invalid 'since' parameter", 400)
        jobs = [serialize_job(job) for job in executor.updates_since(since_value)]
        return jsonify(jobs)

    return app


def serialize_job(job: Job) -> dict:
    return {
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
    }
