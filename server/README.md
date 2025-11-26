# Domainator Server

Prototype Flask-based web UI for running Domainator tools.

## Project status

- HTTP API for uploads, tool discovery, job execution, and log retrieval.
- Initial web UI for file management, tool browsin, and job monitoring.
- Subprocess-backed executor that consumes jsonargparse-compatible config files.
- Schema registry that loads bundled sample tools.
- Pytest coverage for executor, registry, file management, and HTTP endpoints (see `pytest server/tests`).

## Setup

```bash
cd server
pip install -e .[dev]  # or use conda env with requirements from pyproject
domainator-server --data-dir /tmp/domainator_server
```

Use `--schema-dir` to point at additional tool schemas.
Open http://127.0.0.1:5000/ to use the browser interface.
