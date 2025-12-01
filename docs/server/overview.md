# Domainator Web UI

The Domainator web interface bundles a lightweight Flask application, REST API,
and single-page UI for running Domainator tools without leaving the browser. The
server ships as part of the main `domainator` package and reuses the same
schemas, executors, and file-management utilities as the CLI workflows.

## Features
- Upload, catalogue, and download data sets through a simple file manager.
- Discover available tools directly from the schema registry and inspect their
  parameters.
- Launch jobs that execute Domainator CLIs in the background and stream logs in
  real time.
- Retrieve generated outputs and configuration files from the job workspace.

## Quick start
1. Install Domainator with the default extras (they now include the web UI
   dependencies):
   ```bash
   pip install domainator
   # or, from a checkout
   pip install -e .
   ```
   When using the provided `conda_env.yml`, the server requirements are already
   included.
2. Run the server:
   ```bash
   domainator-server --data-dir /tmp/domainator_server
   ```
   Replace the path with a writable directory to store uploads, jobs, and logs.
3. Visit <http://127.0.0.1:8080/> to open the browser UI. Use `--host` and
   `--port` to change the bind address.
4. Pass `--schema-dir` one or more times to load additional tool schema
   directories alongside the bundled defaults.

## Tool schemas and customisation
Schema files live under `src/domainator/server/schemas/` inside the project.
Use [`scripts/server/generate_tool_schemas.py`](../../scripts/server/generate_tool_schemas.py)
to capture argument definitions from CLI entry points and emit JSON schemas for
the web UI. For step-by-step instructions on extending the registry see
[Adding Tools to the Domainator Server](adding_tools.md).

## Development notes
- Configuration defaults are derived from `domainator.server.config.ServerConfig`.
  It provisions upload/output/log directories and infers the active Python
  executable when launching jobs.
- Server-specific tests live under `test/server/`. Run the suite with:
  ```bash
  python -m pytest test/server
  ```
- Static assets (HTML, CSS, JS) are located in `src/domainator/server/templates`
  and `src/domainator/server/static`.
