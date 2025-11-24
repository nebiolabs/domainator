# Overview

- Domainator is a Python-first toolkit plus CLI suite for domain-based genome neighborhood analysis, protein annotation, clustering, and reporting.
- Core workflows revolve around GenBank records augmented with Domainator-specific features (e.g., `qualifier="Domainator"`); HMM profiles are treated as first-class inputs alongside sequence files.
- Every CLI tool has a pure-Python backend reachable through `from domainator.<tool> import <tool>` and many can stream via stdin/stdout to support shell-style chaining or in-process pipelines.
- External tooling such as HMMER, CD-HIT, DIAMOND, UMAP, and Foldseek is optional; most day-to-day functionality runs with the Python dependencies declared in `pyproject.toml` / `conda_env.yml`.
- Use GitHub Issues and Discussions (linked in `README.md`) for bugs, feature requests, and roadmap coordination.

## Repository Structure

- `src/domainator/`: Package code. Each module corresponds to a CLI defined in `pyproject.toml` and exposes `_entrypoint`, `main(argv)`, and core functions. Key helpers include `utils.py` (SeqRecord streaming I/O, CLI scaffolding), `data_matrix.py` (dense/sparse HDF5 matrices), `cytoscape.py`, `foldseek.py`, and the vendored BioPython fork under `Bio/`.
- `test/`: Pytest suite with fixture data in `test/data/`. Tests usually call module-level `main([])` functions to exercise CLI behavior without spawning subprocesses.
- `docs/`: User and contributor docs. `docs/README.md` is the index; other files cover examples, file formats, Foldseek/ESM integration, limitations/FAQ, and developer guidance.
- `conda_env.yml`: Reference environment for end users and contributors. Keep aligned with `project.dependencies` and optional extras in `pyproject.toml`.
- `domainator*.def` / `domainator*.sif`: Apptainer/Singularity definitions and images for containerized execution. Update alongside dependency or CLI changes that affect runtime behavior.

## Architecture

- **CLI entry points**: Scripts like `domainate.py` map to `_entrypoint` functions registered in `pyproject.toml`. `_entrypoint()` reads `sys.argv[1:]`, while `main(argv)` accepts an explicit list for easier testing.
- **SeqRecord extensions**: Domainator wraps a patched BioPython subset in `domainator.Bio`. Always use `domainator.utils.parse_seqfiles` and `write_genbank` to preserve Domainator-specific qualifiers (`feature.type == "Domainator"`, color metadata, database provenance, etc.).
- **Matrix handling**: `data_matrix.DataMatrix` and `SparseDataMatrix` abstract HDF5-backed dense matrices and SciPy sparse matrices, maintaining row/column labels, metadata, and conversions used by `compare_contigs`, `seq_dist`, `build_projection`, and reporting commands.
- **Streaming pipeline**: Many CLIs support `-i -` / `-o -` for stdin/stdout streaming and internally operate on generators to avoid loading entire GenBank files into memory. Favor these helpers when composing new tools.
- **Integrations**: Modules under `seq_dist`, `deduplicate_genbank`, `foldseek`, and `domainator_db_download` coordinate with external binaries or web resources. They perform availability checks and raise descriptive errors; ensure binaries are discoverable on `PATH` in development or document skipped functionality.
- **Visualization & styling**: `cytoscape.py`, `color_genbank.py`, `color_table_to_legend.py`, and plotting utilities generate XGMML, SVG, PNG, and HTML outputs while honoring color palettes from `utils.get_palette` and user-supplied tables.

## Testing

- Run the full suite with `pytest test` from the project root (documented in `docs/developing_domainator.md`). Use `pytest -k <pattern>` for targeted runs.
- Optional coverage workflow: `coverage run -m pytest test`, `coverage report -m`, `coverage html`, then open `htmlcov/index.html`.
- Tests rely on `pytest-datadir` to stage fixtures from `test/data/`; avoid mutating fixture files directly and prefer temporary directories for outputs.
- Most tests are pure Python. A subset exercise integrations via optional binaries—when introducing new external requirements, document them in `conda_env.yml`, `docs/developing_domainator.md`, and ensure tests skip gracefully if the executable is missing.

## Documentation

- `docs/README.md` provides the documentation table of contents referenced by the root `README.md`.
- `docs/examples.md` and the external `domainator_examples` repository showcase step-by-step workflows. Mirror that structure when contributing new tutorials.
- `docs/file_formats.md` and `docs/limitations_and_FAQ.md` define on-disk contracts (GenBank feature schema, matrix encodings, known caveats). Review before changing serialization or metadata conventions.
- `docs/esm_3b_foldseek.md` captures optional Foldseek/ESM integrations; update in tandem with modifications to `domainator.foldseek` or related dependencies.
- `docs/developing_domainator.md` houses contributor-focused tips (testing, PDF builds). Extend this guide when introducing new development tooling, linting, or release steps.