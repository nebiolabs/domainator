# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

All work happens inside the `domainator` conda environment:

```bash
conda activate domainator
```

- Run the full test suite: `pytest test`
- Run one test file: `pytest test/test_domainate.py`
- Run by name pattern: `pytest -k <pattern>`
- Coverage: `coverage run -m pytest test && coverage report -m && coverage html` (open `htmlcov/index.html`)
- Browser/SSN-viewer tests (optional, separate extra): `pip install -e ".[test,browser]" && playwright install chromium` — these skip automatically when `playwright` is not importable.
- Launch the web GUI: `domainator_server --port 8080` (then visit `http://127.0.0.1:8080/`)

Tests use `pytest-datadir`: fixtures live in `test/data/` and are staged via the `shared_datadir` fixture. Do not mutate fixture files; write outputs to `tempfile.TemporaryDirectory()`. Most tests call a tool's `main([...args])` directly rather than spawning subprocesses, and compare outputs with `helpers.compare_seqfiles` / `compare_seqrecords`.

## Architecture

Domainator is a suite of several dozen CLI tools that compose into genome-neighborhood / protein-domain analysis workflows. Each tool is also a Python library function. Three things tie the suite together:

**Every tool follows the same module shape.** A module `src/domainator/<tool>.py` exposes:
- `<tool>(...)` — the core function operating on/yielding `SeqRecord` objects (the reusable library API; import via `from domainator.<tool> import <tool>`).
- `main(argv)` — accepts an explicit arg list (used by tests).
- `_entrypoint()` — calls `main(sys.argv[1:])`; this is what `[project.scripts]` in `pyproject.toml` registers as the installed CLI (e.g. `domainate.py`).

When adding a tool, register its `_entrypoint` in `pyproject.toml` and mirror this shape.

**The CLI argument convention is inverted from typical tools.** `-i` / `--input` is the *file being edited*; editing *criteria* come from other flags (notably `-r` for reference sequences/HMMs). E.g. `domain_search.py -i uniprot.fasta -r query.hmm -o hits.gb` searches the references against the input. Tools are grouped by role: **editors** (input format == output format, chainable), **reports** (TSV/HTML summaries), **comparison** (produce score/distance matrices), and **plotting** (XGMML/SVG/HTML/trees). See the program table in `README.md` for the full catalog.

**Streaming pipeline.** Many editors support `-i -` / `-o -` for stdin/stdout GenBank streaming and operate on generators to avoid loading whole files into memory; you can shell-pipe a streaming-output tool into a streaming-input tool. Chaining the core Python functions directly is faster (avoids GenBank↔SeqRecord round-trips) but loses some I/O validation. Tools that need multiple passes over input cannot stream.

### Key shared infrastructure

- **Vendored BioPython fork** at `src/domainator/Bio/`. Domainator's `SeqRecord` is *not* interchangeable with upstream BioPython. **Always read/write sequence files via `domainator.utils.parse_seqfiles` and `domainator.utils.write_genbank`** to preserve Domainator-specific qualifiers.
- **GenBank as the universal carrier.** Annotations are stored as features with `feature.type == "Domainator"` (constant `DOMAIN_FEATURE_NAME` in `__init__.py`; best-hit features use `DOMAIN_SEARCH_BEST_HIT_NAME`). Changing these constants breaks file compatibility across Domainator versions.
- **Matrices.** `data_matrix.DataMatrix` (HDF5-backed dense) and `SparseDataMatrix` (SciPy sparse) carry row/column labels + metadata and underpin `compare_contigs`, `seq_dist`, `build_projection`, `transform_matrix`, and reporting. The on-disk format is versioned via `_MATRIX_FILE_VERSION` in `data_matrix.py`.
- **External binaries** (HMMER, CD-HIT, DIAMOND, Foldseek, etc.) are optional. Tools that use them check availability and raise descriptive errors; tests that need them should skip gracefully when the binary is absent.
- **Output guardrails.** Tools with potentially large output accept `--max_output_gb` (default 25 GB; `0` disables) to fail early — see `output_guardrails.py`.

### Web server (`src/domainator/server/`)

A Flask app (`app.py`, launched by `cli.py` → `domainator_server`) that surfaces CLIs in a browser UI. Each tool is described by a JSON schema under `server/schemas/`; the server's `tool_registry` / `tool_executor` launch tools by writing parameters to a JSON file and invoking the tool with `--config=<file>` (Domainator CLIs parse this via `jsonargparse`).

To add a tool to the UI: generate a schema with `python scripts/server/generate_tool_schemas.py --pyproject pyproject.toml --output-dir src/domainator/server/schemas/generated/`, then hand-edit it (move niche params into `advanced_parameters`). Full guide: `docs/server/adding_tools.md`. Server design: `docs/server/technical_overview.md`.

## Versioning

The version string lives in `src/domainator/__init__.py` (`__version__`) and **must be updated in `README.md` too**. Breaking changes to the matrix file format require bumping `_MATRIX_FILE_VERSION` in `data_matrix.py`.

## Docs map

`docs/README.md` is the index. `docs/file_formats.md` defines on-disk contracts (review before changing serialization). `docs/developing_domainator.md` covers testing + building the PDF manual. `docs/esm_3b_foldseek.md` covers the optional Foldseek/ESM-2 integration.
