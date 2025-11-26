# Adding Tools to the Domainator Server

This guide walks through the steps required to surface new Domainator tools in the server UI. It covers generating schemas from existing CLIs, crafting schemas manually, and developing tool entry points that work with the server's executors.

## 1. Generate Schemas from Existing CLIs

Most Domainator CLIs already expose a `main(argv)` entry point that the server can introspect. The `scripts/generate_tool_schemas.py` helper captures each CLI's argument parser and emits a JSON schema compatible with the web app.

### Prerequisites

- A Python environment with Domainator (and its optional extras) installed.
- Access to the project `pyproject.toml` that declares the CLI entry points under `[project.scripts]`.

### Command

Run the generator from the `server/` directory:

```
python scripts/generate_tool_schemas.py --pyproject ../pyproject.toml --output-dir src/domainator_server/schemas/generated/
```

Key options:

- `--pyproject`: path to the Domainator `pyproject.toml` containing the CLI script definitions.
- `--output-dir`: destination folder that will receive one `<tool-id>.json` schema per CLI.
- `--include`: optional whitelist of script names (matching `[project.scripts]` keys) when you only need a subset.

Generated schemas can be committed directly or used as a starting point for manual edits.

## 2. Author Schemas Manually

Schemas live under `src/domainator_server/schemas/` and have the following shape:

```
{
  "id": "tool-id",
  "display_name": "Tool Display Name",
  "category": "Domainator",
  "description": "One or more lines describing the tool.",
  "runner": "module",
  "entry_point": "domainator.some_cli",
  "parameters": [ ... ],
  "advanced_parameters": [ ... ]
}
```

### Parameters vs. Advanced Parameters

- `parameters`: primary inputs shown immediately in the UI.
- `advanced_parameters`: optional or niche inputs. These render inside a collapsed "Advanced Parameters" panel.

Each parameter entry supports:

- `name`: human-friendly label.
- `parameter`: the key passed to the CLI when the job runs (defaults to `name` if omitted).
- `type`: `string`, `integer`, `number`, `boolean`, `file`, or `output`.
- `help`: short tip displayed inline in the form.
- `required`: mark required inputs.
- `multiple`: accept a list of values (UI renders a multi-select or collects repeated values).
- `choices`: restrict to specific values.
- `default`: pre-populated default.
- `flags`: optional list of CLI flags (for reference only).
- `file_types`: optional list of MIME-like extensions that filter the file picker when `type` is `file`.

### Tips

- Keep descriptions concise—the first line is used in the tool list; subsequent lines appear in the detail pane with preserved formatting.
- When introducing advanced options, prefer copying the generator's output and moving the relevant entries into `advanced_parameters`.
- Verify the JSON with `python -m json.tool path/to/schema.json` before committing.

## 3. Implement Compatible Runners

The `runner` field determines how the server launches the tool. Supported values are:

| Runner  | Behavior | Entry Point Expectations |
|---------|----------|---------------------------|
| `module` | Invokes `python -m <entry_point> --config=<args.json>` | `<entry_point>` must be an importable module that processes a `--config` file. All Domainator CLIs expose `_entrypoint/main` helpers that accept this flag via `jsonargparse`. |
| `python` | Executes `python <entry_point> --config=<args.json>` | Use when the CLI is a Python script rather than a module. Ensure it accepts the `--config` argument. |
| `shell`  | Runs `/bin/sh -c "<entry_point> --config=<args.json>"` | Provide a shell command string. Useful for wrappers around non-Python binaries. |
| `binary` (default) | Launches `<entry_point> --config=<args.json>` directly | For standalone executables already on `PATH`. |

### Writing Tool Code

1. Accept configuration via JSON:
   - The server writes the collected parameters to a JSON file and passes `--config=/path/to/file` when launching the tool.
   - Domainator CLIs typically use `jsonargparse` or `domainator.utils.configure_cli` to load this config automatically.

2. Emit outputs into the working directory:
   - The executor creates a per-job workspace (`$DATA_DIR/outputs/<job-id>`). Write derived files there.
   - Any produced files are offered to the user. For convenience, call the shared file manager APIs when you want artifacts to appear in the "Stored Files" list immediately.

3. Log to stdout/stderr:
   - All console output is captured in `<job-id>.log`. Avoid using interactive prompts.

4. Exit codes matter:
   - Return `0` on success. Non-zero codes mark the job as failed and surface the log to the user.

### Testing a New Tool

1. Create or update the schema.
2. Restart the development server or call the tool registry refresh endpoint if available.
3. Visit the web UI, locate the tool (verify alphabetical sorting), and inspect the parameters.
4. Run a test job. Confirm:
   - Required and advanced parameters behave correctly.
   - Outputs appear in the job results and optional file manager.
   - Logs contain helpful information on failure.

## 4. Troubleshooting

- **Schema not appearing:** Ensure the `.json` file resides under one of the schema directories configured in `ServerConfig.schema_dirs`.
- **Missing parameter UI control:** Confirm `type` matches one of the supported values and that `parameter`/`name` are present.
- **Job fails instantly:** Check the job log in the UI. Validate the runner/entry point combination and confirm the CLI handles `--config`.
- **Output file not listed:** Verify the tool wrote files to the job workspace and that they are not temporary JSON configs.

With these steps, you can confidently add and maintain tools accessible through the Domainator Server interface.
