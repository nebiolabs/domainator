# Domainator Web Server Design Document

**Version:** 1.0  
**Date:** November 24, 2025  
**Status:** Draft

## Table of Contents
1. [Executive Summary](#executive-summary)
2. [Goals and Non-Goals](#goals-and-non-goals)
3. [Architecture Overview](#architecture-overview)
4. [Technical Stack](#technical-stack)
5. [Installation and Deployment](#installation-and-deployment)
6. [Core Features](#core-features)
7. [User Interface Design](#user-interface-design)
8. [Backend API Design](#backend-api-design)
9. [File Management](#file-management)
10. [Workflow Execution](#workflow-execution)
11. [Workflow System](#workflow-system)
12. [Developer Guide](#developer-guide)
13. [Security Considerations](#security-considerations)
14. [Testing Strategy](#testing-strategy)
15. [Future Enhancements](#future-enhancements)
16. [Open Questions](#open-questions)

---

## Executive Summary

The Domainator Web Server is a lightweight Flask-based local web application that provides a graphical user interface for running Domainator workflows without command-line knowledge. The server will be distributed via conda/bioconda, requiring only conda installation and three simple commands to get started.

**Target user experience:**
```bash
conda create --name domainator_server -c bioconda domainator_server
conda activate domainator_server
domainator_server --port 8080
```

The server will support file upload/download, in-browser visualization of HTML outputs, markdown-based help documentation, execution of individual Domainator tools, and pre-configured multi-step workflows defined by developers.

---

## Goals and Non-Goals

### Goals
1. **Ease of installation**: Single conda package installation
2. **Ease of use**: Web interface accessible to non-programmers
3. **Local execution**: All computation happens on user's machine
4. **File management**: Upload inputs, download outputs, view HTML reports
5. **Tool execution**: Run any Domainator CLI tool via web form
6. **Workflow support**: Pre-defined multi-step workflows created by developers
7. **Developer-friendly**: Easy to add new tools and workflows
8. **Documentation**: Integrated help system with markdown files
9. **No configuration**: Works out-of-the-box with sensible defaults

### Non-Goals
1. **User authentication**: Not implementing multi-user accounts (single-user local server)
2. **Cloud deployment**: Not designing for remote/production deployment
3. **Database backend**: No persistent database (filesystem-based)
4. **Advanced scheduling**: No job queuing beyond simple async execution
5. **Collaboration features**: No sharing, commenting, or version control
6. **Custom visualization**: Will leverage existing Domainator HTML outputs

---

## Architecture Overview

### High-Level Architecture

```
┌─────────────────────────────────────────────────────┐
│                   Web Browser                       │
│  ┌─────────────┐  ┌──────────────┐  ┌────────────┐  │
│  │  File Mgmt  │  │ Tool Executor│  │  Workflow  │  │
│  │     UI      │  │      UI      │  │     UI     │  │
│  └─────────────┘  └──────────────┘  └────────────┘  │
└────────────┬────────────────────────────────────────┘
             │ HTTP/WebSocket
┌────────────▼────────────────────────────────────────┐
│               Flask Web Server                      │
│  ┌──────────────────────────────────────────────┐   │
│  │           REST API Endpoints                 │   │
│  │  /api/files  /api/tools  /api/workflows      │   │
│  └──────────────────────────────────────────────┘   │
│  ┌──────────────────────────────────────────────┐   │
│  │         Business Logic Layer                 │   │
│  │  FileManager  ToolRegistry  WorkflowRegistry │   │
│  └──────────────────────────────────────────────┘   │
└────────────┬────────────────────────────────────────┘
             │
┌────────────▼────────────────────────────────────────┐
│           Domainator Core Library                   │
│  (from domainator.tool import main)                 │
└────────────┬────────────────────────────────────────┘
             │
┌────────────▼────────────────────────────────────────┐
│           Filesystem Storage                        │
│  ~/domainator_server/                               │
│    ├── uploads/                                     │
│    ├── outputs/                                     │
│    └── logs/                                        │
│    ...                                              │
└─────────────────────────────────────────────────────┘
```

### Component Responsibilities

- **Web Browser**: Renders UI, handles user interactions
- **Flask Server**: Routes requests, manages state, orchestrates execution
- **Domainator CLI Tools**: jsonargparse-powered executables run via subprocess with per-job config files
- **Filesystem**: Source of truth for uploads, outputs, logs, and JSON metadata that replaces a traditional database

---

## Technical Stack

### Backend
- **Flask** (2.3+): Web framework
- **Flask-CORS**: Cross-origin support if needed
- **Werkzeug**: File uploads and security
- **jsonschema**: Parameter validation
- **domainator**: Core functionality (required dependency)

### Frontend
- **HTML5/CSS3**: Structure and styling
- **Bootstrap 5**: Responsive UI framework
- **Vanilla JavaScript** or **Alpine.js**: Lightweight interactivity
- **Marked.js**: Markdown rendering for help docs
- **Highlight.js**: Syntax highlighting for code examples

### Process Management
- **subprocess**: Launch each tool or workflow as an isolated process
- **tempfile**: Stage per-job JSON/YAML config passed to jsonargparse CLIs
- **psutil**: Process monitoring (already in domainator deps)
- **asyncio** (optional): Drive non-blocking status polling if needed

### File Handling
- **pathlib**: Path manipulation
- **tempfile**: Secure temporary directories
- **shutil**: File operations

---

## Installation and Deployment

### Package Structure (if separate repository)

```
domainator_server/
├── setup.py / pyproject.toml
├── conda_recipe/
│   └── meta.yaml
├── src/domainator_server/
│   ├── __init__.py
│   ├── app.py              # Main Flask app
│   ├── cli.py              # CLI entry point
│   ├── config.py           # Configuration
│   ├── api/
│   │   ├── files.py        # File management endpoints
│   │   ├── tools.py        # Tool execution endpoints
│   │   └── workflows.py    # Workflow endpoints
│   ├── core/
│   │   ├── file_manager.py
│   │   ├── tool_registry.py
│   │   ├── tool_executor.py
│   │   └── workflow_registry.py
│   ├── workflows/         # Bundled workflow command-line scripts
│   │   ├── genome_annotation.py
│   │   ├── protein_clustering.py
│   │   └── ...
│   ├── templates/          # Jinja2 templates
│   │   ├── base.html
│   │   ├── index.html
│   │   ├── tools.html
│   │   ├── workflows.html
│   │   └── help.html
│   ├── static/
│   │   ├── css/
│   │   ├── js/
│   │   └── docs/          # Markdown help files
│   └── schemas/           # Tool/workflow definition JSONs
│       ├── tools/
│       │   ├── domainate.json
│       │   └── domain_search.json
│       └── workflows/
│           ├── genome_annotation.json
│           └── protein_clustering.json
├── tests/
│   └── test_*.py
└── README.md
```

### Conda Recipe (meta.yaml)

```yaml
{% set name = "domainator_server" %}
{% set version = "0.1.0" %}

package:
  name: {{ name|lower }}
  version: {{ version }}

source:
  url: https://github.com/nebiolabs/domainator_server/archive/v{{ version }}.tar.gz
  sha256: ...

build:
  number: 0
  noarch: python
  script: {{ PYTHON }} -m pip install . --no-deps -vv
  entry_points:
    - domainator_server = domainator_server.cli:main

requirements:
  host:
    - python >=3.10
    - pip
    - setuptools
  run:
    - python >=3.10
    - flask >=2.3
    - werkzeug >=2.3
    - jsonschema >=4.0
    - domainator >=0.7.2
    - psutil >=5.9

test:
  imports:
    - domainator_server
  commands:
    - domainator_server --help

about:
  home: https://github.com/nebiolabs/domainator_server
  license: MIT
  license_file: LICENSE
  summary: 'Web-based GUI for Domainator workflows'
```

### CLI Entry Point (cli.py)

```python
#!/usr/bin/env python3
"""
Command-line interface for domainator_server
"""
import argparse
import sys
import webbrowser
import time
from pathlib import Path
from domainator_server.app import create_app
from domainator_server.config import Config

def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Start Domainator web server',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--port', type=int, default=8080,
                        help='Port to run server on')
    parser.add_argument('--host', default='127.0.0.1',
                        help='Host to bind to (use 0.0.0.0 for network access)')
    parser.add_argument('--data-dir', type=Path, default=None,
                        help='Directory for file storage (default: ~/domainator_server)')
    parser.add_argument('--no-browser', action='store_true',
                        help='Do not automatically open browser')
    parser.add_argument('--debug', action='store_true',
                        help='Run in debug mode')
    
    args = parser.parse_args(argv)
    
    # Setup config
    config = Config(
        data_dir=args.data_dir,
        debug=args.debug
    )
    
    # Create Flask app
    app = create_app(config)
    
    # Print startup message
    url = f"http://{args.host}:{args.port}"
    print(f"Starting Domainator Server...")
    print(f"Data directory: {config.data_dir}")
    print(f"Server running at: {url}")
    print(f"Press Ctrl+C to stop")
    
    # Open browser
    if not args.no_browser and args.host in ['127.0.0.1', 'localhost']:
        def open_browser():
            time.sleep(1.5)  # Give server time to start
            webbrowser.open(url)
        import threading
        threading.Thread(target=open_browser, daemon=True).start()
    
    # Run server
    app.run(host=args.host, port=args.port, debug=args.debug)

if __name__ == '__main__':
    sys.exit(main())
```

---

## Core Features

### 1. File Management

**Upload**
- Drag-and-drop or file picker interface
- Support for `.gb`, `.gbk`, `.fasta`, `.fa`, `.hmm`, `.hdf5`, `.tsv` formats
- File validation (format checking)
- Automatic organization by file type
- Size limits (configurable, default 5GB per file)

**Download**
- Single file download
- Batch download as ZIP
- Direct links to outputs

**File Browser**
- Hierarchical view: uploads/, outputs/, by timestamp
- File metadata: name, size, type, upload date
- Preview for text files (GenBank, FASTA first 100 lines)
- Delete/rename capabilities

### 2. Tool Execution

**Tool Discovery**
- Auto-discover all Domainator tools from `pyproject.toml` entry points
- JSON schema definitions for each tool's parameters
- Group tools by category (Editors, Reports, Comparisons, etc.)

**Dynamic Form Generation**
- Generate web forms from tool schemas
- Input types: file selector, text, number, boolean, dropdown
- Parameter validation (required vs optional, type checking)
- Help text for each parameter (extracted from argparse help)

**Execution Interface**
- Start button to launch tool
- Real-time progress indicator (spinning/progress bar)
- Log streaming (optional)
- Cancel button to terminate
- Results display upon completion

**Output Handling**
- List of generated output files
- In-browser viewer for HTML outputs
- Download links for all outputs
- Link outputs to inputs (provenance tracking)

### 3. HTML Viewer

**In-Browser Display**
- Iframe-based viewer for HTML reports
- Full-screen mode
- Responsive layout
- Support for interactive elements (d3.js visualizations)

**Supported HTML Types**
- `summary_report.py` HTML output
- `enum_report.py` HTML tables
- `matrix_report.py` histograms
- `plot_contigs.py` interactive plots
- Any other HTML generated by tools

### 4. Help System

**Markdown Documentation**
- Convert existing docs from `/docs` to web-friendly format
- Navigation sidebar
- Search functionality
- Code syntax highlighting
- Links to tool-specific help

**Contextual Help**
- Per-tool documentation
- Parameter descriptions
- Example workflows
- Common error messages and solutions

---

## User Interface Design

### Layout

```
┌────────────────────────────────────────────────────────┐
│  Domainator Server            [User] [Settings] [Help] │
├──────────┬─────────────────────────────────────────────┤
│          │                                             │
│  Files   │  Main Content Area                          │
│  Tools   │                                             │
│  Jobs    │                                             │
│  Help    │                                             │
│          │                                             │
│          │                                             │
│          │                                             │
└──────────┴─────────────────────────────────────────────┘
```

### Page Views

#### Home/Dashboard
- Quick stats: # files uploaded, # jobs run
- Recent jobs with status
- Quick start guide
- Example workflows

#### Files Page
```
┌────────────────────────────────────────────────────┐
│  Files                            [Upload Files]   │
├────────────────────────────────────────────────────┤
│  ┌──────────────────────────────────────────────┐  │
│  │  Drag files here or click to browse          │  │
│  └──────────────────────────────────────────────┘  │
│                                                    │
│  Uploaded Files                     [Delete All]   │
│  ┌──────────────────────────────────────────────┐  │
│  │ 📄 genome.gb       1.2 MB    2025-11-24      │  │
│  │ 📄 pfam.hmm       15.3 MB    2025-11-24      │  │
│  │ ...                                          │  │
│  └──────────────────────────────────────────────┘  │
│                                                    │
│  Output Files                                      │
│  ┌──────────────────────────────────────────────┐  │
│  │ 📊 annotated.gb    1.5 MB   [View][Download] │  │
│  │ 📈 report.html     150 KB   [View][Download] │  │
│  └──────────────────────────────────────────────┘  │
└────────────────────────────────────────────────────┘
```

#### Tool Execution Page
```
┌────────────────────────────────────────────────────┐
│  Tools > domainate.py                              │
├────────────────────────────────────────────────────┤
│  Annotate sequences with HMM profiles              │
│                                                    │
│  Input File (-i)        [Select] genome.gb         │
│  Reference File (-r)    [Select] pfam.hmm          │
│  Output File (-o)       annotated.gb               │
│  E-value threshold      10            [?]          │
│  CPU threads            4             [?]          │
│  Max overlap            0.6           [?]          │
│                                                    │
│  [Show Advanced Options]                           │
│                                                    │
│                     [Run Tool]                     │
└────────────────────────────────────────────────────┘
```

---

## Backend API Design

### RESTful Endpoints

#### File Management

```
GET    /api/files                    # List all files
POST   /api/files/upload             # Upload file(s)
GET    /api/files/<file_id>          # Get file metadata
GET    /api/files/<file_id>/content  # Download file
DELETE /api/files/<file_id>          # Delete file
GET    /api/files/<file_id>/preview  # Preview text file
```

#### Tool Management

```
GET    /api/tools                    # List all tools
GET    /api/tools/<tool_name>        # Get tool schema
POST   /api/tools/<tool_name>/execute # Execute tool
GET    /api/tools/<tool_name>/help   # Get tool help text
```

#### Job Management

```
GET    /api/jobs                     # List all jobs
GET    /api/jobs/<job_id>            # Get job status
POST   /api/jobs/<job_id>/cancel     # Cancel running job
GET    /api/jobs/<job_id>/logs       # Get job logs
DELETE /api/jobs/<job_id>            # Delete job record
```

> Implementation detail: the job listing endpoint walks `data_dir/jobs/`, reading each `job.json` manifest to assemble the response. Missing or malformed JSON files are skipped with a warning.

#### Workflow Management

```
GET    /api/workflows                 # List available workflows
GET    /api/workflows/<workflow_id>   # Get workflow schema/info
POST   /api/workflows/<workflow_id>/execute # Execute workflow
```

### Example API Responses

#### Tool Execution Request
```json
POST /api/tools/domainate/execute
{
  "parameters": {
    "input_file": "file_123",
    "reference_files": ["file_456"],
    "output": "annotated.gb",
    "cpu": 4,
    "evalue": 10,
    "max_overlap": 0.6
  }
}
```

#### Tool Execution Response
```json
{
  "job_id": "job_789",
  "status": "running",
  "started_at": "2025-11-24T10:30:00Z",
  "estimated_duration": null
}
```

#### Job Status Response
```json
{
  "job_id": "job_789",
  "tool": "domainate",
  "status": "completed",  // queued, running, completed, failed, cancelled
  "started_at": "2025-11-24T10:30:00Z",
  "completed_at": "2025-11-24T10:35:00Z",
  "progress": 100,
  "output_files": ["file_999"],
  "logs": "/api/jobs/job_789/logs",
  "error": null
}
```

---

## File Management

### Storage Structure

```
~/domainator_server/
├── config.json              # Server configuration
├── uploads/                 # User uploaded files
│   ├── <file_id>.gb
│   ├── <file_id>.hmm
│   └── ...
├── outputs/                 # Tool output files
│   ├── <job_id>/
│   │   ├── annotated.gb
│   │   ├── report.html
│   │   └── metadata.json
│   └── ...
├── jobs/                   # Job metadata
│   ├── <job_id>.json
│   └── ...
└── logs/                   # Execution logs
    ├── <job_id>.log
    └── ...
```

### Metadata Files (JSON)

- **Per-artifact metadata**: Every upload, generated file, or workflow artifact gets a sibling JSON file (e.g., `annotated.gb` + `annotated.gb.json`) describing name, size, provenance, and tags.
- **Job manifests**: Each job directory contains a `job.json` file summarizing status, parameters, timestamps, and output file references.
- **Discovery strategy**: When a user opens the Files or Jobs page, the server walks the corresponding directories, reads available JSON files, and builds the tables on-the-fly.
- **Error handling**: Invalid or unreadable JSON files are skipped with a warning in the server logs but never block page rendering.
- **Advantages**: Avoids the complexity of an external database while remaining consistent with the file-centric workflows Domainator already uses.

### File Metadata

Each file has associated metadata:
```json
{
  "file_id": "file_123",
  "original_name": "genome.gb",
  "path": "/uploads/file_123/genome.gb",
    "metadata_path": "/uploads/file_123.json",
  "size": 1234567,
  "type": "genbank",
  "uploaded_at": "2025-11-24T10:00:00Z",
  "checksum": "sha256:...",
  "description": "User-provided annotation"
}

> The server never caches these records in a separate database; instead it re-parses the JSON file whenever a list view is requested.
```

### Cleanup Policy

- Keep uploads until explicitly deleted
- Auto-delete job outputs after 7 days (configurable)
- Warn before deleting files referenced by active jobs
- Option to archive outputs as ZIP before cleanup

---

## Workflow Execution

### Tool Execution Flow

```
1. User submits form → Validate parameters against schema
2. Create job record → Generate job_id and working directory
3. Stage input files → Resolve file_ids to absolute paths within job workspace
4. Write args file → Serialize parameters to JSON/YAML understood by jsonargparse
5. Spawn subprocess → Execute CLI with `--config` (or equivalent) pointing to args file
6. Monitor process → Stream logs, capture exit status, enforce timeouts
7. Update manifest → Persist job state transitions to `jobs/<job_id>/job.json`
8. Discover outputs → Record artifacts produced in the job workspace
9. Notify UI → Surface completion, log location, and output downloads
```

### Implementation (tool_executor.py)

```python
import json
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Any
from dataclasses import dataclass
from enum import Enum
from tempfile import NamedTemporaryFile

class JobStatus(Enum):
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"

@dataclass
class Job:
    job_id: str
    tool_name: str
    parameters: Dict[str, Any]
    status: JobStatus
    process: subprocess.Popen | None = None
    started_at: float | None = None
    completed_at: float | None = None
    output_files: List[str] | None = None
    error: str | None = None
    log_file: Path | None = None
    work_dir: Path | None = None

class ToolExecutor:
    def __init__(self, config, registry):
        self.config = config
        self.registry = registry  # ToolRegistry instance responsible for schema lookup
        self.jobs: Dict[str, Job] = {}
    
    def execute_tool(self, tool_name: str, parameters: Dict[str, Any]) -> str:
        job_id = self._generate_job_id()
        work_dir = self._prepare_job_dir(job_id)
        log_file = self.config.data_dir / "logs" / f"{job_id}.log"

        job = Job(
            job_id=job_id,
            tool_name=tool_name,
            parameters=parameters,
            status=JobStatus.QUEUED,
            output_files=[],
            log_file=log_file,
            work_dir=work_dir,
        )
        self.jobs[job_id] = job
        self._write_job_manifest(job)

        cfg_path = self._write_args_file(parameters, work_dir)
        cmd = self._build_command(tool_name, cfg_path)

        log_file.parent.mkdir(parents=True, exist_ok=True)
        log_handle = log_file.open("w", encoding="utf-8")

        proc = subprocess.Popen(
            cmd,
            cwd=work_dir,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            text=True,
        )

        job.process = proc
        job.status = JobStatus.RUNNING
        job.started_at = time.time()
        self._write_job_manifest(job)

        return job_id

    def poll_jobs(self):
        for job in self.jobs.values():
            if job.status != JobStatus.RUNNING or job.process is None:
                continue
            retcode = job.process.poll()
            if retcode is None:
                continue
            job.completed_at = time.time()
            if retcode == 0:
                job.status = JobStatus.COMPLETED
                job.output_files = self._discover_outputs(job)
            else:
                job.status = JobStatus.FAILED
                job.error = f"Process exited with code {retcode}"
            if job.log_file:
                job.log_file.close()
            self._write_job_manifest(job)

    def cancel_job(self, job_id: str):
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
        self._write_job_manifest(job)

    def _write_args_file(self, parameters: Dict[str, Any], work_dir: Path) -> Path:
        with NamedTemporaryFile("w", delete=False, suffix=".json", dir=work_dir) as fh:
            json.dump(parameters, fh, indent=2)
            return Path(fh.name)

    def _build_command(self, tool_name: str, cfg_path: Path) -> List[str]:
        schema = self._schema_for(tool_name)
        runner = schema.get("runner", "binary")
        entry = schema["entry_point"]

        if runner == "python":
            return [self.config.python_executable, entry, f"--config={cfg_path}"]
        elif runner == "module":
            return [self.config.python_executable, "-m", entry, f"--config={cfg_path}"]
        else:
            return [entry, f"--config={cfg_path}"]

    def _schema_for(self, tool_name: str) -> Dict[str, Any]:
        """Fetch the schema describing how to invoke a tool/workflow."""
        schema = self.registry.get(tool_name)
        if schema is None:
            raise KeyError(f"Unknown tool id: {tool_name}")
        return schema

    def _write_job_manifest(self, job: Job):
        job_dir = self.config.data_dir / "jobs" / job.job_id
        job_dir.mkdir(parents=True, exist_ok=True)
        manifest_path = job_dir / "job.json"
        payload = {
            "job_id": job.job_id,
            "tool": job.tool_name,
            "status": job.status.value,
            "parameters": job.parameters,
            "started_at": job.started_at,
            "completed_at": job.completed_at,
            "output_files": job.output_files,
            "error": job.error,
            "log_file": str(job.log_file) if job.log_file else None,
            "work_dir": str(job.work_dir) if job.work_dir else None,
        }
        with manifest_path.open("w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2)
```

---

### Workflow System

### Overview

Workflows are packaged as standalone command-line programs, discovered alongside built-in tools. Each workflow is simply another jsonargparse-enabled script that accepts a `--config` file and emits outputs in its working directory. The server treats workflows exactly like any other tool: it writes a config file, calls the executable, and surfaces outputs and logs.

**Key principles:**
- Workflows live in a dedicated directory (e.g., `share/domainator_server/workflows/`).
- Each workflow is a Python script or module entry point invoked via CLI.
- Authors compose Domainator tools by chaining subprocess calls internally; no server-side base class is required.
- Registration metadata (id, display name, schema, help link) is expressed in the same JSON schema format used for core tools.

### Workflow Discovery

The registry walks configured tool directories (core CLIs, workflow scripts, optional external binaries) and registers any schema files it finds:

```python
# domainator_server/core/tool_registry.py
from pathlib import Path
import json

class ToolRegistry:
    def __init__(self, schema_dirs: list[Path]):
        self.schemas = {}
        for directory in schema_dirs:
            for schema_path in directory.glob("*.json"):
                with schema_path.open("r", encoding="utf-8") as handle:
                    schema = json.load(handle)
                tool_id = schema["id"]
                self.schemas[tool_id] = schema | {"schema_path": schema_path}

    def list_tools(self):
        return [schema for schema in self.schemas.values()]

    def get(self, tool_id: str):
        return self.schemas.get(tool_id)
```

Each schema declares the command to run and any runtime expectations, decoupling the server from tool implementation details.

### Workflow Packaging Guidelines

- Provide a top-level `main()` generated by jsonargparse.
- Accept `--config /path/to/config.{json,yaml}` for batch execution.
- Write outputs into the current working directory or to explicitly declared paths supplied in the config.
- Document expected outputs in the schema so the server can surface them post-run.
- Version scripts through filenames or schema metadata (`version` field) to help users distinguish revisions.

---

## Developer Guide

### Tool and Workflow Schemas

Every executable exposed through the web UI is described by a JSON schema. The schema captures metadata (id, label, category), invocation details (runner, entry point), default configuration, and parameter definitions so the UI can build forms automatically.

Schemas live under `schemas/` and can be split into logical subdirectories (`tools/`, `workflows/`, `external/`) without affecting discovery. The registry reads all `*.json` files and trusts the `id` field to be unique.

### Adding a New Tool

1. **Author or identify the CLI.** Ensure it supports jsonargparse-style config input (Domainator CLIs already do).
2. **Create schema file** (e.g., `schemas/tools/new_tool.json`):

```json
{
    "id": "new_tool",
    "display_name": "New Tool",
    "category": "Editors",
    "description": "Annotate sequences with a custom rule set",
    "help_url": "/help/tools/new_tool",
    "runner": "module",
    "entry_point": "domainator.new_tool",
    "default_config": {
        "output": "output.gb",
        "cpu": 4
    },
    "parameters": [
        {
            "name": "input_file",
            "type": "file",
            "file_types": ["genbank", "fasta"],
            "required": true,
            "description": "Input sequence file"
        },
        {
            "name": "output",
            "type": "output",
            "required": true,
            "description": "Output file path"
        },
        {
            "name": "threshold",
            "type": "number",
            "default": 10,
            "min": 0,
            "description": "Threshold value"
        },
        {
            "name": "cpu",
            "type": "integer",
            "default": 4,
            "min": 1,
            "max": 64,
            "description": "Number of CPU threads"
        },
        {
            "name": "verbose",
            "type": "boolean",
            "default": false,
            "description": "Enable verbose output"
        }
    ]
}
```

3. **Document usage** in `static/docs/tools/new_tool.md` for in-app help.
4. **Optional:** ship presets by dropping YAML files into `presets/new_tool/` that the UI can offer as one-click configurations.

### Adding a Workflow

Workflows are simply CLIs stored in a dedicated directory (e.g., `workflows/`). To register one:

1. Implement the script with jsonargparse, accepting `--config` for batch execution.
2. Place the executable in the workflows directory or expose it as a console entry point via `pyproject.toml`.
3. Add a schema (e.g., `schemas/workflows/genome_annotation.json`) pointing at the script:

```json
{
    "id": "genome_annotation",
    "display_name": "Genome Annotation Workflow",
    "category": "Workflows",
    "description": "Search, annotate, and visualise neighborhoods",
    "help_url": "/help/workflows/genome_annotation",
    "runner": "python",
    "entry_point": "workflows/genome_annotation.py",
    "parameters": [
        {"name": "genome", "type": "file", "file_types": ["genbank", "fasta"], "required": true},
        {"name": "query", "type": "file", "file_types": ["hmm", "fasta"], "required": true},
        {"name": "database", "type": "file", "file_types": ["hmm"], "required": true},
        {"name": "evalue", "type": "number", "default": 1e-10},
        {"name": "context_range", "type": "integer", "default": 10}
    ],
    "expected_outputs": [
        {"name": "annotated", "type": "genbank", "description": "Annotated neighborhoods"},
        {"name": "report", "type": "html", "description": "Interactive visualization"}
    ]
}
```

4. Provide workflow-specific docs in `static/docs/workflows/<id>.md`.
5. Test locally with `python workflows/genome_annotation.py --config examples/genome_annotation.yaml` so the same invocation the server uses is validated.

### Runner Types

The `runner` field tells the executor how to call the tool:

- `binary`: invoke `entry_point` directly (tool already on `$PATH`).
- `module`: run `python -m entry_point ...`.
- `python`: run `python entry_point ...`, with `entry_point` pointing to a script path relative to the install root or absolute.
- `shell`: opt-in for commands that must run via shell; disable by default for safety.

Extended runners (e.g., `container`, `conda`) can be added later without changing the execution contract.

### Testing Checklist

- Validate schemas with `jsonschema` to catch typos.
- Run each tool/workflow manually using the generated config file.
- Confirm outputs land in the job workspace and match schema expectations.
- Ensure help docs and UI descriptions reference the same parameter names.
- If the tool needs external binaries, document detection logic in the schema (`requires_binary`).

---

## Security Considerations

### Single-User Local Deployment
- Bind to `127.0.0.1` by default and display a warning before allowing `0.0.0.0`.
- Treat the active OS user as fully trusted; no authentication UI is planned for v1.
- Store all data under a user-controlled directory (`--data-dir`) to keep permissions predictable.

### Process Isolation
- Execute tools and workflows in subprocesses to prevent shared state leakage.
- For cancellation, send `SIGTERM` followed by `SIGKILL` if the process fails to exit within a grace period.
- Record the command line and config path in job manifests to aid post-mortem analysis.

### File Safety
- Validate uploads with extension and basic magic checks before placing them in `uploads/`.
- Use `werkzeug.utils.secure_filename` to avoid path traversal.
- Enforce configurable size limits (default 5 GB) and refuse archives containing symlinks outside the workspace.

### Data Privacy
- Never transmit telemetry; any future metrics collection must be opt-in.
- Provide a one-click “Purge Data” action that removes uploads, jobs, outputs, and logs.
- Encourage users to run within encrypted home directories when handling sensitive genomes.

---

## Testing Strategy

### Unit Tests
- `test_file_manager.py`: filesystem create/read/delete flows and quota handling.
- `test_tool_registry.py`: schema loading, runner dispatch, and validation errors.
- `test_tool_executor.py`: config file generation, subprocess lifecycle, cancellation, and log capture.
- `test_workflow_registry.py`: multi-directory discovery and schema precedence.

### Integration Tests
- `test_api_endpoints.py`: REST coverage for files, tools, workflows, and jobs.
- `test_full_workflow.py`: Submit representative payloads and assert outputs appear in `outputs/<job_id>/`.
- `test_workflow_execution.py`: Smoke-test bundled workflows via real subprocess execution.

### UI/End-to-End (Optional)
- Browser automation (Playwright/Selenium) for upload, execution, and result preview flows.
- Visual regression on critical pages (Jobs dashboard, Tool form, HTML viewer).

### Test Data
- Reuse fixtures under `test/data/` and add lightweight JSON/YAML configs mirroring UI defaults.
- Keep large datasets out of the repo; prefer Git LFS or downloadable archives referenced in docs.

Continuous integration should run unit and integration suites on pull requests, reserving UI tests for nightly builds or tagged releases.

---

## Future Enhancements

### Phase 2
- Job history persisted in SQLite (optional feature flag) with filtering and search.
- Server-Sent Events for real-time job updates without polling.
- Parameter preset library curated by developers and shareable between users.
- Workflow dependency checks with proactive binary detection UI.

### Phase 3
- Visual workflow authoring that exports jsonargparse scripts.
- Batch execution against multiple inputs with aggregated reporting.
- Desktop notifications (Electron/Toast) when long jobs finish.
- Optional container runners (Apptainer/Podman) for sandboxing heavy tools.

### Stretch Ideas
- Remote execution adapters (SSH, SLURM) while retaining local UI.
- Pluggable authentication layer for institutions that need multi-user access.
- Result comparison dashboards to diff outputs across job runs.

---

## Open Questions

### 1. Repository Structure

**Question**: Should the web server remain in the main Domainator repository or ship as a companion project?

**Option A: Separate Repository (`domainator_server`)**
- ✅ Clear separation of dependencies and release cadence.
- ✅ Easier for CLI-only users to avoid web requirements.
- ✅ Dedicated issue tracker and documentation site.
- ❌ Requires careful version compatibility matrix.
- ❌ Slightly more complex contributor workflow.

**Option B: Integrated (`domainator.server` subpackage)**
- ✅ Single installation and version number.
- ✅ Shared utility modules without duplication.
- ✅ Tight coupling guarantees parity with CLI behaviour.
- ❌ Web dependencies imposed on all users.
- ❌ Larger distribution footprint.
- ❌ Mixed concerns within a single repository.

**Recommendation**: Start with a separate repo while the UI stabilises; revisit consolidation once stable APIs emerge.

### 2. Configuration Options

**Question**: What configuration options should be exposed beyond `--port` and `--data-dir`?

**Possible options**:
- `--max-upload-size`: Max file upload size
- `--max-concurrent-jobs`: Limit parallel executions
- `--job-retention-days`: Auto-delete outputs after N days
- `--default-cpu`: Default CPU count for tools
- `--temp-dir`: Location for temporary files
- `--log-level`: Verbosity (DEBUG, INFO, WARNING, ERROR)
- `--theme`: UI color scheme

**Recommendation**: Start minimal (just port/data-dir), add others based on user feedback.

### 3. Job Persistence

**Question**: Should job history persist across server restarts?

**Option A: In-memory only**
- ✅ Simple implementation
- ✅ No database needed
- ❌ Lose history on restart
- ❌ Can't resume interrupted jobs

**Option B: Persistent (JSON files)**
- ✅ Survives restarts
- ✅ No external dependencies
- ❌ Potential file locking issues
- ❌ Limited query capabilities

**Option C: SQLite database**
- ✅ Survives restarts
- ✅ Efficient queries
- ✅ ACID guarantees
- ❌ Adds dependency
- ❌ More complex

**Recommendation**: Start with JSON files (Option B) for simplicity, migrate to SQLite in Phase 2.

### 4. Real-time Updates

**Question**: How should the UI get updates on job progress?

**Option A: Polling**
- ✅ Simple to implement
- ✅ Works with any web server
- ❌ Delayed updates
- ❌ Wasteful (many empty responses)

**Option B: WebSockets**
- ✅ Real-time updates
- ✅ Server can push updates
- ❌ More complex implementation
- ❌ Requires websocket support

**Option C: Server-Sent Events (SSE)**
- ✅ Real-time updates
- ✅ Simpler than WebSockets
- ✅ Built into Flask
- ❌ One-way only

**Recommendation**: Start with polling (every 2-3 seconds), migrate to SSE in Phase 2 if needed.

### 5. Progress Indicators

**Question**: Can we show meaningful progress (% complete) for tool execution?

**Challenge**: Most Domainator tools run as subprocesses without emitting structured progress events.

**Possible solutions**:
- **Option A**: Indeterminate spinner (no % complete)
- **Option B**: Estimate based on input file size and historical data
- **Option C**: Modify tools to support progress callbacks (requires Domainator changes)
- **Option D**: Parse log output for clues about progress

**Recommendation**: Start with Option A (spinner), evaluate Option B in Phase 2.

### 6. Multi-user Support

**Question**: Even for local use, should we support multiple browser sessions?

**Considerations**:
- User might have multiple workflows in different tabs
- Family members might use same machine
- File conflicts and race conditions

**Recommendation**: Design to be multi-session safe (use proper locking), but don't implement authentication. Each tab can submit jobs independently.

### 7. External Tool Dependencies

**Question**: How to handle optional dependencies (cd-hit, diamond, usearch)?

**Options**:
- Auto-detect and disable tools that require missing binaries
- Show warning in UI when tool is unavailable
- Provide installation instructions
- Bundle binaries with conda package (where possible)

**Recommendation**: Auto-detect and show clear error messages with installation instructions when tool requires missing binary.

### 8. Error Handling

**Question**: How detailed should error messages be for end users?

**Balance**:
- Too technical: Confusing for non-programmers
- Too simple: Not enough info to debug

**Recommendation**: 
- Show simple error message in UI ("Tool failed: Invalid input file format")
- Provide "Show Details" button with full stack trace
- Link to relevant help documentation
- Save full error to log file

---

## Conclusion

This design document outlines a comprehensive plan for a Flask-based Domainator web server that makes domain-based genome analysis accessible to users without command-line experience. The architecture prioritizes simplicity, ease of installation, and extensibility while maintaining the power and flexibility of the underlying Domainator toolkit.

The modular design allows for incremental development:
1. **Phase 1**: Basic file management and single-tool execution
2. **Phase 2**: Workflow catalog and JSON-based job persistence
3. **Phase 3**: Advanced features and visual workflow authoring

**Next Steps**:
1. Resolve open questions (especially repository structure)
2. Create initial project structure
3. Implement core file management and tool registry
4. Build basic web interface
5. Test with representative workflows
6. Gather user feedback and iterate

---

## Appendix: Technology Alternatives Considered

### Web Frameworks
- **Flask** ✓: Lightweight, flexible, good for small-scale
- **FastAPI**: Modern, async, but overkill for local server
- **Django**: Too heavy, includes features we don't need
- **Streamlit**: Great for data apps, but less flexible for custom workflows

### Frontend Frameworks
- **Vanilla JS + Bootstrap** ✓: Simple, no build step, fast
- **React/Vue**: Powerful but adds complexity and build step
- **Svelte**: Interesting but less mature ecosystem
- **Alpine.js** ✓: Good middle ground, considered as optional enhancement

### Job Execution
- **Subprocess** ✓: Matches CLI ergonomics, isolates state, simplest to ship
- **Process Pools**: Useful for controlled concurrency limits, but adds orchestration overhead
- **Celery/RQ**: Provide distributed queues, but introduce external brokers (overkill for single-user)
- **Container Runners**: Strong isolation, higher setup cost (suitable for later phases)

### File Storage
- **Filesystem** ✓: Direct, simple, sufficient
- **Object storage**: S3-compatible, but unnecessary for local
- **Database BLOBs**: Poor performance for large files
