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
- `<tool>(input, *, ...)` — the core function operating on/yielding `SeqRecord` objects (the reusable library API; import via `from domainator.<tool> import <tool>`). Only the main input is positional; every other parameter is keyword-only and named like the command's option (its argparse dest). Matrix producers return a `DataMatrix`; `main` writes it.
- `build_parser()` — builds the command's options on `domainator.cli.common.tool_parser(__doc__)`, without parsing. The dispatcher, shell completion, and the schema generator all use it.
- `main(argv)` — `build_parser().parse_args(argv)`, then runs the tool; accepts an explicit arg list (used by tests).

Users run tools through one dispatcher, `domainator <command> [options]` (`src/domainator/cli/`). Its command list is the static table in `cli/registry.py`, so dispatching imports only the requested tool's module. When adding a tool, add it to `cli/registry.py` and mirror this shape. `domainator completion bash|zsh|tcsh` prints a static shell completion script generated from every `build_parser()`.

**The CLI argument convention is inverted from typical tools.** `-i` / `--input` is the *file being edited*; editing *criteria* come from other flags (notably `-r` for reference sequences/HMMs). E.g. `domainator domain_search -i uniprot.fasta -r query.hmm -o hits.gb` searches the references against the input. Tools are grouped by role: **editors** (input format == output format, chainable), **reports** (TSV/HTML summaries), **comparison** (produce score/distance matrices), and **plotting** (XGMML/SVG/HTML/trees). See the program table in `README.md` for the full catalog.

**Streaming pipeline.** Many editors support `-i -` / `-o -` for stdin/stdout GenBank streaming and operate on generators to avoid loading whole files into memory; you can shell-pipe a streaming-output tool into a streaming-input tool. Chaining the core Python functions directly is faster (avoids GenBank↔SeqRecord round-trips) but loses some I/O validation. Tools that need multiple passes over input cannot stream.

### Key shared infrastructure

- **Vendored BioPython fork** at `src/domainator/Bio/`. Domainator's `SeqRecord` is *not* interchangeable with upstream BioPython. **Always read/write sequence files via `domainator.seqio.parse_seqfiles` and `domainator.seqio.write_genbank`** (also exported from the package root) to preserve Domainator-specific qualifiers.
- **Shared modules** (internal; the package root re-exports the public pieces): `seqio` (sequence file I/O and parser backends), `hmmio` (HMM/CM/reference files, alphabets), `model` (`DomainatorCDS`, CDS naming, `ANNOTATION_TYPES`), `locations` (slicing and padding), `partial` (fragment status), `taxonomy` (NCBI taxonomy, and `taxonomy.records` for per-record filtering), `expr` (`BooleanEvaluator`), `colors` (palettes, color tables), `offsets` (record offsets and partitioning), `external` (`resolve_binary`/`run_command` for every external program), `sparse_scores` (chunked sparse score matrices for the matrix tools), `table_writers` (TSV/NDJSON/HTML tables for the reports), `_mp` (multiprocessing). CLI-only helpers live in `cli/common.py`. Qualifier values are read with `model.feature_evalue`/`feature_score`, and every `-e` filter uses `model.evalue_passes` (keeps e-values <= the cutoff).
- **GenBank as the universal carrier.** Annotations are stored as features with `feature.type == "Domainator"` (constant `DOMAIN_FEATURE_NAME` in `__init__.py`; best-hit features use `DOMAIN_SEARCH_BEST_HIT_NAME`). Changing these constants breaks file compatibility across Domainator versions.
- **Matrices.** `data_matrix.DataMatrix` (a SciPy CSR array, read from hdf5 or dense text) carries row/column labels + metadata and underpin `compare_contigs`, `seq_dist`, `kmer_dist`, `structure_dist`, `build_projection`, `transform_matrix`, and reporting. The on-disk format is versioned via `_MATRIX_FILE_VERSION` in `data_matrix.py`; exactly symmetric matrices are stored as their upper triangle. Matrix tools share their `--hdf_out`/`--dense_text_out`/`--symmetrize` options via `matrix_output.py`.
- **External binaries** (HMMER, CD-HIT, DIAMOND, Foldseek, etc.) are optional. Tools that use them check availability and raise descriptive errors; tests that need them use the `@pytest.mark.requires_binary("name")` marker (see `test/conftest.py`) so they skip when the binary is absent. Tests that download data are marked `@pytest.mark.network`.
- **Output guardrails.** Tools with potentially large output accept `--max_output_gb` (default 25 GB; `0` disables) to fail early — see `output_guardrails.py`.
- **Errors and logging.** Raise a `domainator.errors` exception (`InputError`, `ExternalToolError`, `OutputSizeLimitExceeded`, `FormatVersionError`, all `DomainatorError`) for problems the user can fix; the dispatcher prints them as one line with exit status 1. Library code never calls `sys.exit` or prints progress: use `logger = logging.getLogger(__name__)` (progress at INFO, problems at WARNING) and `warnings.warn` for data warnings. `tool_parser` adds `-v/--verbose`/`-q/--quiet` and configures the `domainator` logger when the arguments are parsed.

### Web server (`src/domainator/server/`)

A Flask app (`app.py`, launched by `cli.py` → `domainator_server`) that surfaces CLIs in a browser UI. Each tool is described by a JSON schema under `server/schemas/`; the server's `tool_registry` / `tool_executor` launch tools by writing parameters to a JSON file and invoking the tool with `--config=<file>` (Domainator CLIs parse this via `jsonargparse`).

Schemas use `"runner": "domainator"` with the command name as `entry_point`, and `test/server/test_domainator_schemas.py` checks that every schema parameter is still an option of its tool, so a CLI rename must update the schemas. To add a tool to the UI: generate a schema with `python scripts/server/generate_tool_schemas.py --output-dir src/domainator/server/schemas/generated/ --include <command>`, then hand-edit it (move niche params into `advanced_parameters`). Full guide: `docs/server/adding_tools.md`. Server design: `docs/server/technical_overview.md`.

## Spelling

Use **American English spelling** everywhere: code identifiers, string literals, comments, docstrings, docs, and commit messages. Prefer `color` / `behavior` / `neighbor` / `center` / `labeled` / `modeled` / `canceled` / `analog` / `catalog` / `-ize` / `-ization` / `toward` over their British forms. Two exceptions: the vendored BioPython fork under `src/domainator/Bio/` stays as upstream wrote it, and standard identifiers from external specs keep their own spelling (e.g. the HTML attribute `aria-labelledby`).

## Versioning

The version string lives in `src/domainator/__init__.py` (`__version__`) and **must be updated in `README.md` too**. Breaking changes to the matrix file format require bumping `_MATRIX_FILE_VERSION` in `data_matrix.py`. Changes to the `.dsnv` Domainator Similarity Network Viewer bundle schema require bumping `SSN_VIEWER_BUNDLE_VERSION` in `ssn_bundle.py`; if the change is additive, add the previous version to `SUPPORTED_SSN_VIEWER_BUNDLE_VERSIONS` so existing bundles keep loading (see the `.dsnv` section of `docs/file_formats.md`). Don't bump versions of anything without first confirming with the user.

## Docs map

`docs/README.md` is the index. `docs/file_formats.md` defines on-disk contracts (review before changing serialization). `docs/developing_domainator.md` covers testing and CI. `docs/esm_3b_foldseek.md` covers the optional Foldseek/ESM-2 integration.
