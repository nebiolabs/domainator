# Domainator 1.0 readiness: audit findings and release plan

## Context

Domainator is at 0.12.0 and has grown to 46 CLI tools without a stability commitment. The goal is a stable, consistent, well-designed 1.0, and breaking changes are acceptable. This plan comes from a read-only audit of:
- the CLI surface (every parser);
- the library internals;
- packaging, file formats, docs, tests and the server;
- two follow-up investigations (Biopython migration, GUI interface).

The plan lists the issues found and orders them into release phases.

**Decisions already made with the maintainer**
- **CLI:** one `domainator <subcommand>` dispatcher with shell tab completion. No plain-name aliases. The `foo.py` scripts go away; this is a clean break.
- **Deprecations:** clean break plus `MIGRATION.md`. No alias shims for renamed flags, and the WIP `--sparse`/`--dense` aliases are dropped too.
- **GUI:** moves to a separate project. Domainator's job is a stable, versioned, machine-readable tool manifest, generated from the parsers, from which the GUI derives its forms with zero hand edits.
- **Biopython:** investigate migrating to upstream for 1.0 (findings in §5). Expect that it may not be feasible.
- **Round-two decisions:**
  - Remove `partition_seqids` and `partition_seqfile` (Phase 0).
  - Bump the matrix file version to 2.0, but not `__version__` (Phase 0).
  - `_range` means both sides (Phase 3).
  - Every output path is an explicit `--<content>_out` channel, and channels can be combined, e.g. GenBank and FASTA together (Phase 3).
  - Nothing reads stdin or writes stdout by default; `-i -` is required for stdin (Phase 2d).
  - Remove `--offset`/`--recs_to_read`; users build an index instead (Phase 0).
  - The Python floor is 3.12 (Phase 8).
  - stdin is supported only for single-pass tools (Phase 3).
  - Python 3.15 is a target (Phase 8).

**The working tree has uncommitted WIP:** the matrix-output refactor (`matrix_output.py`, the `data_matrix.py` rewrite, matrix version 1.1). Phase 0 lands it before anything else.

## Progress

Work happens on the `version_one` branch.

- **Phase 0: done.**
  - Matrix file version 2.0; readers accept major versions 1 and 2. The `--sparse`/`--dense` aliases and the deprecated `DataMatrix.write` output types are removed.
  - CI test workflow (`.github/workflows/tests.yml`) and pytest config, with `requires_binary`/`network` markers and `tmp_path_retention_policy = "failed"`.
  - Stale files deleted; `partition_seqids`, `partition_seqfile`, `--offset` and `--recs_to_read` removed.
  - Deviations:
    - `AGENTS.md` and `CLAUDE.md` serve different readers (AI agents using the CLI vs. developers), so both stay.
    - `package.sh` and the README.pdf build are removed. Users install from a clone of the public GitHub repository.
- **Phase 1: done.**
  - `--aligner_arg` is replaced by `--foldseek_args`/`--reseek_args`. Each is one shell-split string, and options Domainator sets itself are refused.
  - `allow_abbrev=False` was done in Phase 2; it exposed abbreviated options in three test files, which now spell options out.
  - Each fix has a regression test.
  - The GenBank/FASTA editors now write compressed output through `open_writable_seqfile`, and `filter_domains` closes its output.
- **Phase 2a–2c: done.**
  - Every tool has `build_parser()` on `cli.common.tool_parser`, which sets `description=`, `allow_abbrev=False` and a version epilog.
  - The `domainator` dispatcher (`src/domainator/cli/`) and `python -m domainator` work, along with `domainator completion {bash,zsh,tcsh}` (static shtab scripts).
  - Deviations:
    - The command metadata lives only in the static `cli/registry.py` table, instead of also as a `TOOL` constant in each module, so it isn't duplicated.
    - `domainator --help` lists the commands, so there is no separate `list` command.
  - Two module renames: `domainator_db_download` runs as `db_download` and `domainator_format_db` as `format_db`.
  - The `<tool>.py` scripts and every `_entrypoint`/`__main__` block are gone. `[project.scripts]` has only `domainator` (and `domainator_server` until Phase 9).
- **Phase 2d: done.**
  - `cli.common` holds the shared helpers: `add_input`, `add_output`, `add_seq_outputs`, `add_switch`, `add_annotation_type`, `check_outputs`, `reject_stdin`, `resolve_cpu`, `parse_keyword_args`.
  - Nothing reads stdin or writes stdout by default. `-`, `/dev/stdin` and `/dev/stdout` are recognized, at most one output may go to stdout, and a tool given no output fails with "No output specified".
  - `outputs.py`: atomic temp-file outputs, compression inferred from the extension, and `SeqOutputs`, which writes GenBank and FASTA in one pass and serializes each record once for the size guardrail.
  - `domainator.errors`: `DomainatorError` with `InputError`, `ExternalToolError`, `OutputSizeLimitExceeded` and `FormatVersionError` (each also a `ValueError` or `RuntimeError`). The dispatcher prints them as `domainator <command>: error: ...` with exit status 1 (`-v` keeps the traceback). The twelve `raise SystemExit(str(exc))` conversions and seq_dist's `exit(1)` calls are gone.
  - Logging: `tool_parser` adds `-v/--verbose` and `-q/--quiet` to every command (so they can follow the command name, and go in config files), and `ToolParser.parse_args` configures the `domainator` logger: INFO by default, DEBUG with `-v`, errors only with `-q` (which also hides Domainator's Python warnings). The remaining progress `print`s are `logger.info`/`logger.warning`; the debug dumps in exception handlers are exception notes. A library import logs nothing below WARNING. build_tree's own `--quiet` is replaced by the shared one.
  - Help: the formatter leaves out empty defaults (None, False, empty lists), so help text says what leaving an option out means; the hand-written phrases that repeated a printed default are removed. Phrases that explain a `None` default ("default: all contigs") stay.
- **Phase 3: done (see MIGRATION.md).**
  - The `_out` output naming, separate GenBank/FASTA channels, `--reference`, `--min_score`/`--threshold`, `--max_hits_per_cds`, `--annotation_type`, `_range` meaning both sides, `--anchor`, `--<program>_args` (shell-split, reserved options refused), `--three_d`, `--co_occurrence_*`, lowercase `--esm2_3di_*`.
  - `--fasta_type` is `--alphabet {nucleotide,protein}` (`cli.common.add_alphabet`). hmmer_build and hmmer_select keep their own `--alphabet {amino,dna,rna}`, the HMM alphabets.
  - Every `-e/--evalue` is a jsonargparse `PositiveFloat` and every `--min_evalue` a `NonNegativeFloat`, so the check happens at parse time and the bound is visible to the Phase 9 manifest. HMMER-backed tools also cap `-e` at 10.
  - Binary paths: deduplicate_genbank's `--bin_path` became `--cdhit_path`/`--cdhit_est_path`/`--usearch_path`/`--diamond_path`; seq_dist gained `--diamond_path` and domainate `--foldseek_path`.
  - Gaps filled: `-e`/`--databases` for trim_contigs (which also fixes `--unannotated` doing nothing on its own, and an `IndexError` with `--cds_range` on a contig without CDSs), and `--color_table` for plot_contigs.
  - Deviations:
    - `domain_search --deduplicate` stays off by default. Turning it on to match `select_by_cds` silently dropped distinct hits, because peptide hits can share a record ID.
    - domainate keeps two devices (`--esm2_3di_device`, default `cuda:0`, and `--foldseek_device`, default `cpu`). ESM2 inference and the Foldseek search usually run on different devices, because GPU Foldseek needs a special build. Tools with one device use `--device`.
    - `--domain_expr` stays in the CDS- and contig-level selectors (select_by_cds, select_by_contig, trim_contigs). It is an expression over the set of domains in a CDS or contig, so it has no meaning for tools that act on single domain annotations (extract_*, filter_domains, rotate_contigs) or rank by score (sort_contigs).
    - No `--keep_name` for domain_search and select_by_cds. They output new regions, and the region names keep several regions from one contig distinct (`--deduplicate` relies on it).
    - Report columns behind one `--columns` option moves to Phase 9, with the enum_report column catalog.
- **Server (interim, until Phase 9).** Schemas use a new `"runner": "domainator"` (`python -m domainator <command>`); the old `python -m domainator.<module>` runner silently did nothing once the `__main__` blocks were removed. The schemas were remapped to the 1.0 option names, `generate_tool_schemas.py` reads `cli/registry.py`, and `test_schema_parameters_match_the_parser` fails whenever a schema names an option its tool no longer has.
- **Phase 11 (early): docs check.** `test/test_docs_commands.py` checks every option shown with a command in a code block or inline code span of README.md, AGENTS.md and docs/. It also parses each complete command in a code block and fails when the command names none of its outputs. It found stale options: `seq_dist -o`, `plot_contigs -o`, and the `--aligner_arg` example. All docs now use `domainator <command>` and the 1.0 option names.

---

## Phase 0: Foundation (do first)

1. **Land the matrix WIP.**
   - Bump `_MATRIX_FILE_VERSION` to **2.0** (approved). `SPARSE_CSR_UPPER_TRIANGLE` can't be read by 0.12 readers, and `file_formats.md:238,272` says such a change bumps the major version.
   - Readers accept 1.x and 2.x and reject a newer major version with a clear error.
   - Update the `file_formats.md` matrix section.
   - **The Domainator `__version__` is not touched.** You will bump it yourself.
   - Remove the `--sparse`/`--dense` aliases from `matrix_output.py` (clean break).
   - Remove the stale TODO at `data_matrix.py:275`.
   - Make `_get_file_type` public, or reuse `utils.get_file_type`, rather than having `matrix_output` import a private.
2. **Add a CI test workflow** (`.github/workflows/tests.yml`):
   - pytest on Linux and macOS for every supported Python version, 3.12 through 3.15 (see Phase 8);
   - one job with the conda env, for the binary-dependent tests;
   - a docs-reference check (Phase 2).
   This is the safety net for every later phase.
3. **Test hygiene:**
   - Add `[tool.pytest.ini_options]` with `testpaths`, registered markers (`network`, `requires_binary`, `browser`, `slow`) and `filterwarnings = error::FutureWarning`.
   - Add skips for missing cd-hit, usearch and diamond (`test_deduplicate_genbank.py`, the diamond tests in `test_seq_dist.py`).
   - Mark the network tests in `test_domainator_db_download.py:77-121`.
4. **Repo hygiene:**
   - Add `/TMbed/` to `.gitignore`.
   - Delete stale tracked files: `.github/copilot-instructions.md` (or regenerate it from CLAUDE.md), `docs/web_server_design.md` (it describes files that don't exist), `Bio/SeqIO/InsdcIO_new.py`, `server/schemas/archived/`, `server/workflows/` (a mock).
   - Drop `README.pdf` and `package.sh` in favor of wheels/sdist and bioconda.
   - Consolidate `AGENTS.md`/`CLAUDE.md`.
5. **Remove `partition_seqids` and `partition_seqfile`.** Both are superseded by the partitioning built into domainate and domain_search.
   - Delete both modules, their tests, the `partition_seqfile.py` pyproject entry, the README row (`README.md:123`), `docs/examples.md:204` and the archived schema.
   - `domain_search` still imports `partition_seqfile.i_partition_seqfiles` (`domain_search.py:26,218,542,552`). Move that function, and the `partition_seqfile()` it calls, into the internal `index` module (Phase 4 split) as private helpers, keeping their tests.
   - Also remove the obsolete `--offset`/`--recs_to_read` options from domainate (`domainate.py:1356-1358`) and domain_search, along with the `seek_to`/`max_recs` plumbing that only they use (keep whatever `parse_seqfiles` callers still need). Docs point users to building an index (`domainator_format_db --index`) instead.

## Phase 1: Bugs found by the audit (fix regardless of 1.0)

| Bug | Location |
|---|---|
| `extract_peptides -e 0.1` is parsed as str, causing a `TypeError` in a float comparison | `extract_peptides.py:219` (add `type=float`) |
| `hmmer_compare --cpu 0` crashes on `Pool(0)`; the default of 8 is inconsistent with other tools | `hmmer_compare.py:83,120` |
| `get_domain_names(evalue=)` drops hits *better* than the threshold | `utils.py:799` |
| A `KeyError` path prints the record to **stdout**, corrupting piped GenBank | `utils.py:946` |
| `domainate(hits_only)`: the last batch is yielded in hit order, not input order; the membership test is O(n·h) on a list | `domainate.py:1272,1282` |
| A bare `except:` around `overlay` silently skips annotations (and catches `KeyboardInterrupt`) | `domainate.py:876,890` |
| `_merge_intervals` is defined twice; the second shadows the first | `utils.py:1827,2394` |
| `summary_report -o -` creates a file named `-` | `summary_report.py:730` |
| `select_by_cds -o x.gb.gz` writes uncompressed output; compression only works when a path reaches `write_genbank` | `select_by_cds`, `trim_contigs`, … |
| The documented example `--aligner_arg --prefilter-mode 2` fails ("expected at least one argument"), and a repeated `--aligner_arg=` keeps only the last value (verified) | `structure_lib.py:1181` (fixed by the Phase 3 `--<program>_args` design) |
| Help text says the wrong default | `filter_domains --annotation_type`, `color_genbank --mode`, cpu help ("all cores" when physical cores are used) |
| `max_recs` is applied per file, not overall as documented | `utils.py:1208` |
| The `BooleanEvaluator` `)` tokenizing edge case | `utils.py:832-835` |
| `ssn_navigator.json` uses param `id`, but the CLI dest is `cluster_id`, so it's broken in the GUI | moot after Phase 9; fix if released before then |
| The README headline example uses `enum_report.py --name`, which doesn't exist (it prefix-matches `--named_domain_count`) | `README.md:274` |
| Prefix matching of abbreviated options causes the bug above | set `allow_abbrev=False` on all parsers (Phase 2) |

## Phase 2: CLI platform

**2a. A tool-module contract.** Every tool module exposes:
```python
TOOL = ToolInfo(name="domainate", category="annotate", summary="…")   # static metadata
def build_parser() -> ArgumentParser: ...   # no heavy imports here
def <tool>(...): ...                        # library function (Phase 4)
def main(argv=None) -> int: ...             # parse, call <tool>, return exit code
```

**2b. The `domainator` dispatcher** (new `src/domainator/cli/__main__.py`, plus `domainator/__main__.py` so that `python -m domainator` works):
- The registry is a static table `{name: "domainator.<module>"}` with summaries, so no tool module is imported to list or dispatch commands. The registry is the single source for the dispatcher, completion, the manifest and the docs.
- `domainator <tool> …` imports exactly one module.
- `domainator --version`, `domainator --help` (grouped by category) and `domainator list`.
- `pyproject` `[project.scripts]` shrinks to `domainator = "domainator.cli:main"`. This also fixes the broken `domainator:X._entrypoint` form, which `importlib.metadata` can't load.
- The command name is the module name without the suffix. The ones to rename are decided in Phase 3 (e.g. `domainator_db_download` → `db_download`, `domainator_format_db` → `format_db`).

**2c. Tab completion:**
- `domainator completion {bash,zsh,tcsh}` prints a static shtab script built from `build_parser()` of every registered tool. It uses jsonargparse 4.49's `get_completion_script`/shtab support, and `shtab` becomes a dependency.
- The script is static, so pressing Tab costs nothing even though some tool imports are slow. Document a one-line install (`domainator completion bash > ~/.local/share/bash-completion/completions/domainator`).
- File arguments complete to files, choices complete to their values.

**2d. Shared CLI conventions:**
- Enforced in a `domainator.cli.common` helper module that all parsers use: `add_input`, `add_output`, `add_cpu`, `add_evalue`, `add_domain_selection`, `add_contig_selection`, `add_fasta_out`, `add_max_output_gb`, plus the existing `matrix_output.add_matrix_output_arguments`. Each helper also attaches the manifest metadata (Phase 9).
- Rules:
  - `description=` instead of the docstring as `prog`; 42 of 46 parsers pass it as `prog` today, which mangles the usage line.
  - `allow_abbrev=False`.
  - All parsers use jsonargparse.
  - The formatter shows defaults; remove the ~140 hand-written "default:" phrases in help text.
  - **I/O channels are always explicit** (see the "I/O channel rules" in Phase 3):
    - No tool reads stdin or writes stdout by default. Today 23 tools read stdin when `-i` is omitted and 28 write to stdout when no output is given.
    - The user names the channel to pipe: `-i -` for input; `--fasta_out -` or `-o /dev/stdout` for output. `-` and `/dev/stdout`/`/dev/stdin` are recognized as the same stream. The shared output helper must special-case them, because the guardrail writers' temp-file + `os.replace` (`output_guardrails.py`) fails on a device path.
    - At most one output channel may target stdout. A tool given no output channel errors out ("no output specified").
    - Compression is inferred from the extension for every GenBank/FASTA output.
  - One `--cpu` semantics: `0` = all *available* cores, via `os.sched_getaffinity` with a fallback, in one `resolve_cpu()` helper that replaces 7 copies.
- **Errors and logging:**
  - `main()` catches `DomainatorError` and prints a one-line message to stderr. The exit code is 1 for runtime errors and 2 for usage errors via `parser.error`; uncaught exceptions keep their traceback.
  - Add a `logging` setup with `-v/--verbose` and `-q/--quiet` on the dispatcher. Replace the ~80 `print(..., file=sys.stderr)` progress messages with `logging`. Keep `warnings.warn` for data warnings.

## Phase 3: CLI vocabulary normalization (one breaking pass, recorded in MIGRATION.md)

Proposed canonical spellings, based on the audit's inventory of name drift:

| Concept | Current variants | 1.0 |
|---|---|---|
| Reference | `--references` (5 tools) / `--reference` (5 tools) | `-r/--reference` (repeatable) |
| E-value cutoff | `-e` with defaults 0.001 / None / 1e8; validation `<=10` vs `>0` | `-e/--evalue`. Filter tools default to None (no filter); search tools default to 0.001. One validator |
| Score cutoff | `--score`, `--score_cutoff`, `--search_score_lb/ub`, `--lb`, `--threshold` | `--min_score` / `--max_score`. For matrices, `--min_value`, replacing the overloaded `--lb` |
| Length bounds | `--lb` (extract_unannotated), `--length_lb/ub`, `--cds_count_lb/ub` | `--min_length/--max_length`, `--min_cds/--max_cds` |
| Hit caps | `--max_hits` (4 meanings), `--max_domains`, `--max_hits_per_contig`, `--first`, `--num_recs`, `--recs_to_read`; "0 = unlimited" handled inconsistently | `--max_hits_per_cds`, `--max_hits_per_contig`, `--max_hits_per_query`, `--max_records`; None = unlimited, 0 is an error |
| Annotation source | `--search_hits`, `--domain_type`, `--annotation_type` | `--annotation_type {domainator,domain_search,all}` |
| Domain lists | `--domains` + `--domains_file` (sort_contigs lacks the file form); `--domain_expr` only in 3 tools | the trio `--domains/--domains_file/--domain_expr` in every domain-filtering tool |
| Databases | `--databases`, `--databases_keep/--databases_filter`; trim_contigs has none | `--databases` everywhere a filter applies |
| Unannotated | `--unannotated` / `--no_domain` | `--unannotated` |
| CDS lists | `--cds/--cds_file`, `--cdss/--cdss_file` | `--cds/--cds_file` |
| Neighborhood / trimming | `--kb_range/up/down`, `--cds_range/up/down`, `--kb_both/left/right`, `--cds_both/left/right`, `--pad_up/down`, `--where` / `--pad_where` | `_range` always means "both sides". Strand-relative tools (select_by_cds, domain_search, extract_domains padding) use `--kb_up/--kb_down/--kb_range` and `--cds_up/--cds_down/--cds_range`. Contig-coordinate trimming (trim_contigs) uses `--kb_left/--kb_right/--kb_range` and `--cds_left/--cds_right/--cds_range`. The anchor option is `--anchor` |
| Dedup | `--deduplicate` (opt-in) / `--skip_deduplicate` (opt-out); `--both_strands` / `--dedup_both_strands` | `--deduplicate/--no_deduplicate` with one default; `--both_strands` |
| Taxonomy update | `--taxonomy_update` / `--skip_taxonomy_update` | `--taxonomy_update/--no_taxonomy_update` (jsonargparse/`BooleanOptionalAction`) |
| Output format | `--fasta_out` is a *boolean* that switches the one `-o` stream from GenBank to FASTA (clean_sequences, deduplicate_genbank, extract_domains, extract_peptides, extract_unannotated, rotate_contigs, select_by_contig) | Separate channels: `-o/--genbank_out PATH` and `--fasta_out PATH`, usable together. One pass writes to both (a shared `write_outputs(records, genbank=…, fasta=…)` helper in `seqio`). Every GenBank editor gets `--fasta_out` |
| Input format / molecule type | `--fasta_type`, `--alphabet`, `--input_type` (a file format in seq_dist and kmer_dist, but the matrix data type in transform_matrix) | `--alphabet {protein,nucleotide}`; `--input_format` (file-format override only); transform_matrix gets `--value_type` |
| Device | `--esm2_3Di_device`, `--foldseek_device`, `--gpu_device`, `--device` | `--device` |
| Binary paths | `--bin_path`, `--cmbuild_path`, `--foldseek_path`, … (diamond has none) | `--<binary>_path` for every external binary |
| Extra program args | `--params '"-c":4, "-b":0.4'` (a JSON body without braces; seq_dist, deduplicate_genbank, build_projection) / `--aligner_arg` (`nargs='+'`; broken, see Phase 1) | `--<program>_args` (see below) |
| Output directory | `--outdir` / `--output_dir` / `--output_prefix` | `--output_dir` |
| Output-path naming | Mixed: `--png_out`, `--coords_out`, `--color_table_out` vs `--hdf`, `--html`, `--json`, `--xgmml`, `--newick`, `--svg`, `--cluster_tsv`, `--hits_tsv`, `--domains_table`, `--cluster_table` | Every output-path flag ends in `_out` (see below) |
| Misc | `--3d` (needs `vars(p)['3d']`), `--co_occurence_*` (misspelled), `--esm2_3Di_*` (mixed case) | `--three_d`, `--co_occurrence_*`, lowercase throughout |

**Output-path naming rule (no exceptions):**
- Every output-path flag is named `--<content>_out`: `--genbank_out`, `--fasta_out`, `--hdf_out`, `--dense_text_out`, `--tsv_out`, `--html_out`, `--json_out`, `--xgmml_out`, `--newick_out`, `--svg_out`, `--png_out`, `--coords_out`, `--color_table_out`, `--cluster_tsv_out`, `--hits_tsv_out`.
- Input flags never end in `_out`; they are named for their content (`--color_table`, `--metadata`, `--coords`).
- This settles the `color_genbank` case: `--color_table` is read and `--color_table_out` is written. My earlier proposal to rename `--color_table_out` to `--color_table` would have collided with the input flag.
- Each tool's primary channel also gets the short alias `-o`: `--genbank_out` in editors, `--tsv_out` in reports, `--hdf_out` in matrix tools. So `-o` stays easy to type, and its meaning is visible in `--help` and in the manifest.
- The manifest lint (Phase 9) enforces the rule: a flag ends in `_out` if and only if its role is output.

**Passing extra arguments to external programs:**
- One string option per program: `--diamond_args`, `--foldseek_args`, `--reseek_args`, `--cdhit_args`, `--umap_args`, …
- The whole value is a single shell token, split with `shlex.split`. That keeps the program's flags apart from Domainator's: `--foldseek_args="--prefilter-mode 2 --exhaustive-search"`.
- Recommend the `=` form in docs and examples. argparse rejects a following token that starts with `-` and contains no space, so `--foldseek_args --exhaustive-search` would fail; the error message should suggest the `=` form.
- In a JSON config the value is a string (or a list of strings).
- Each backend declares the options that Domainator sets itself (seq_dist already keeps such a list for diamond: `--outfmt`, `-p`, `--max-target-seqs`, …). Passing a reserved option is an error, not a silent conflict.
- UMAP keyword arguments (build_projection `--params`) aren't command-line flags. They become `--umap_args="n_neighbors=30 min_dist=0.1"`, parsed as `key=value` pairs with `ast.literal_eval` values.

**Long flags reused with different meanings** (`--domains`, `--score`, `--strand`, `--length`, `--append`, `--acc`, `--desc`, `--alphabet`, `--id`, `--mode`): mostly in the report tools, where column-selection flags clash with filter names. Fix by moving report columns behind one `--columns name,score,…` list (this also simplifies the GUI; see Phase 9).

**I/O channel rules:**
- `-i -` is allowed only for tools that read their input exactly once. Tools that need several passes reject `-` with a clear message ("this tool reads its input more than once; pass a file"), and the manifest declares `stdin: false` for them.
- I checked the tools I had earlier listed for stdin support:
  - `domainate`: multi-pass. It peeks each file's molecule type (`domainate.py:99`) and counts peptides for Z (`:1389`) before the main parse.
  - `seq_dist`: multi-pass (parses the input at `seq_dist.py:567` and again at `:683`).
  - `kmer_dist`: multi-pass (a peek at `kmer_dist.py:557`, then the parse at `:233`).
  - `hmmer_compare`: two passes over the query today (`common_hmm_alphabet` at `hmmer_compare.py:63`, then iteration at `:82`). It can become single-pass by switching to `iter_hmms_with_alphabet`, as hmmer_search already does. Then it accepts `-i -`.
  - `find_features`: single pass (`find_features.py:853`), so it gets `-i -`.
  - **Dropped** from the stdin proposal: domainate, seq_dist, kmer_dist.

**Also fill the symmetry gaps:**
- `--fasta_out` for every GenBank editor;
- `--max_output_gb` for domainate, select_by_cds, structure_to_genbank;
- `--databases`/`-e` for trim_contigs;
- `--keep_name` for domain_search and select_by_cds;
- `--color_table` for plot_contigs;
- `--tsv_out` for matrix_report.

## Phase 4: Library API

- **Define the public API.** `domainator/__init__.py` re-exports the core functions lazily (PEP 562 `__getattr__`) and has an `__all__`. Everything else in the package is private (`_`-prefixed) or documented as internal. Move `RawAndDefaultsFormatter` out of the package root.
- **Consistent core-function signatures:**
  - The first parameter is `records: Iterable[SeqRecord]` for every record transformer, which returns `Iterator[SeqRecord]`.
  - All options are keyword-only (`*,`), because `domain_search` and `domainate` have 27 or 28 positional parameters today.
  - Validate eagerly (a wrapper around the inner generator, as `rotate_contigs.py:266` does).
  - Uniform names: `cpu` (not `cpus`/`threads`), `evalue`, `databases`. Fix the `unanntotated` typo (`select_by_contig.py:90`).
  - Add a record-based core where one is missing: `genbank_to_fasta`, `deduplicate_genbank` (today a class method that takes file paths), `domainator_format_db` (`format_db`), `plot_contigs`, `build_ssn_viewer`. Move the dedup and padding logic out of `select_by_cds.main` (`:624-656`) into the library, shared with `domain_search`.
  - Matrix producers (`seq_dist`, `kmer_dist`, `structure_dist`, `compare_contigs`) return a `DataMatrix` instead of writing files. Writing goes through `matrix_output`.
- **Keep CLI concerns out of the library.** No `sys.exit`/`exit(1)` (`seq_dist.py:475-560`, `deduplicate_genbank.py:53,62`, `structure_lib.py:1385`), no stdout printing, no progress `print`s; use `logging`.
- **Exceptions:** a `domainator.errors` module with `DomainatorError` and subclasses `InputError`, `ExternalToolError`, `OutputSizeLimitExceeded`, `FormatVersionError`. Re-parent the five existing ad-hoc classes onto it.
- **Split `utils.py` (2,746 lines)** along the seams the audit found: `_mp`, `hmmio`, `seqio`, `index`, `model` (CDS naming, `DomainatorCDS`), `expr` (`BooleanEvaluator`: compile once to an AST instead of per-record `eval`), `locations`, `partial`, `taxonomy`, `colors` (also move `read_color_table` out of the `color_genbank` tool).
- **Qualifier schema as code:** a `domainator.schema` module with constants for feature types and qualifier keys (`name`, `database`, `evalue`, `score`, `cds_id`, `rstart`, …), plus typed accessors (`domain_evalue(f)`, `domain_score(f)`). These replace about 11 copies of `float(f.qualifiers['evalue'][0])` and their inconsistent `<` vs `<=` comparisons.
- **Unify overlap filtering:**
  - Today `domainate.filter_by_overlap` rejects at `>=` and `filter_domains.domain_overlap` at `>`.
  - `find_features` builds the domain qualifiers by hand instead of calling `domain_feature_qualifiers`.
- **Shared helpers:**
  - one external-binary runner (`resolve_binary` / `run_command` / `ExternalToolError`), replacing four styles;
  - one report-writer base, replacing the three `Summary*Writer` trios and `enum_report`'s TODOs;
  - one color system (drop seaborn husl or route it through `colors`).
- **Import hygiene:**
  - Remove the 14 module-level `warnings.filterwarnings("ignore", module='numpy')` calls and the global `BiopythonParserWarning` suppression (`utils.py:137`). Library imports must not change process-wide warning state.
  - Make heavy imports lazy: torch/esmologs in `foldseek.py` (and fail loudly instead of `except ImportError: pass`), umap in `build_projection` (`build_tree` takes 7.1 s just to import `write_cytoscape_xgmml`; move that to `cytoscape.py`), seaborn in `cytoscape.py` (unused) and `summary_report`, and numba in `hmmer_search`.
  - Target: under 0.5 s for any `domainator <tool> --help`.
- **Dead code** (delete, per the remove-dead-code policy): `parse_simple_list`, `DEFAULT_BUFFER_SIZE`, `_native_lean_usable`, `color_genbank.auto_color_map`, `pyhmmer_decode` (only needed for pyhmmer <0.12; 41 call sites), the legacy `--color_domains`/`--color_both` flags, the `data_matrix` re-exports of `ssn_edges`, about 70 lines of commented-out argparse blocks, the dead `except` at `enum_report.py:839`, and `jsonschema` in deps.
- **Parser backend:**
  - Replace the global `DOMAINATOR_GB_PARSER` env override with an explicit per-call argument. The override currently hands lean objects to tools that can't handle them.
  - Define a `Protocol` for the record types `parse_seqfiles` can yield.
  - Fix the Rust/Python drift: `/transl_table` is ignored in `rust/src/parser.rs:161,356`; FASTA whitespace is handled differently (`parser.rs:938-945`); the FASTA fast path has no fallback (`lean_record.py:476`); the native parsers build whole partitions as lists.
  - Add a parity test that runs every fixture through both parsers.

## Phase 5: Vendored Biopython → upstream plus a thin layer (feasible for 1.0)

**Finding: a hybrid migration is feasible, with an estimated 1–1.5 weeks of work.**
- I prototyped it in memory: upstream Biopython 1.87 plus about 150 lines of GenBank subclasses.
- Results on 41 fixture files (194 records):
  - every record parsed identically to the fork;
  - the written output was byte-identical for all 40 files;
  - parse speed was the same as the fork or better (`pDONR` ×300: fork 0.56 s, layer 0.60 s, plain upstream 0.73 s; `MT_nbs` ×25: fork 0.44 s, layer 0.33 s).
- The fork is 23,199 lines. About 750 lines would remain.

**What the fork actually changes, and how the layer reproduces each change:**

| Fork behavior | Why it matters | Layer |
|---|---|---|
| A 6-token GenPept LOCUS branch (`Scanner.py:1006-1012`) | Upstream misparses long-name **protein** records (`FeSOD_20*.gb`): it drops `molecule_type` and `date`, sets strand +1, and then fails to write. Long *nucleotide* names are fine upstream. | Scanner subclass |
| No qualifier wrapping except `/translation` (`InsdcIO.py:280-284`) | Upstream wraps at 80 columns, so a long value with no spaces (fixture values reach 388 characters) **comes back with spaces inserted**. This is a real data-integrity problem. | Writer subclass |
| The NCBI quote-escape regex is removed from `feature_qualifier` (`GenBank/utils.py:1228`) | This is the main reason upstream parses 1.3–1.45x slower | Consumer subclass |
| Footer checks and the LOCUS-length warnings are disabled | Noise and speed | Subclass + warning filter |
| Location algebra on `SimpleLocation`/`CompoundLocation` (`SeqFeature.py:1388-1590,1970-2310`) | About 45 call sites in 12 files, plus 34 in `test_SeqFeature.py` | Free functions in `domainator/locations.py` that use only public attributes |
| SeqRecord hit-annotation getters/setters and comparison operators keyed on `hit_best_score` | About 35 getter/setter sites. The operators are used only by 3 heaps (`domain_search.py:149,308`, `sort_contigs.py:82`) | Annotation helpers in `domainator.schema`; heaps use `(score, counter, rec)` tuples, as `structure_search.py` already does |

**Steps:**
1. Add `domainator/genbank.py`: a scanner, consumer and writer subclass, plus `parse_genbank`/`write_genbank`. Pin `biopython>=1.85,<1.89`, because the layer touches private internals (`_FeatureConsumer`, `_feed_first_line`, `_write_feature_qualifier`). A round-trip test over every fixture in CI guards the upper bound.
2. Add `domainator/locations.py` with the location algebra as free functions, and move `test_SeqFeature.py` to `test_locations.py`.
3. Rewrite about 124 `domainator.Bio` imports in `src/`, 25 test files and `lean_record.py`. Fix `summary_report.py:543`, which uses `feature.strand` (removed upstream).
4. Filter upstream's extra `BiopythonParserWarning`s narrowly, inside the layer and not process-wide, so they don't flood users.
5. Delete `src/domainator/Bio/`. The license then becomes plain MIT, and the Phase 8 BSD-3 packaging issue goes away.
6. Update CLAUDE.md and `docs/file_formats.md`. Domainator records are now ordinary Biopython objects, so the warning that they are "not interchangeable" can go.

**Expected behavior changes** (record them in MIGRATION.md; some test expectations may need updating):
- Heap tie order may differ.
- Nested `join(order(...))` warns and sets the location to None, instead of aborting the parse.
- A **fork bug** goes away: `GenBank/utils.py:1128` passes `seq_type.lower()` as `is_circular`, which wraps reversed joins on linear records as if they crossed the origin.
- Other fork bugs that disappear with the fork: a missing `return` at `SeqFeature.py:1511`, and `raise f"..."` at `SeqRecord.py:1311,1333,1342`.

**Fallback:** if the private-API coupling proves too fragile, keep a trimmed fork. Delete `InsdcIO_new.py`, the unused `SeqUtils`/`Data` modules and EMBL/IMGT; ship `Bio/LICENSE.rst`; and fix the three bugs above.

## Phase 6: File formats and stability contract

- Write `docs/stability.md`, stating what 1.x promises: the CLI flags/semantics, the public Python API (`__all__`), the GenBank qualifier schema, and versioned on-disk formats. Anything else is internal.
- Every format gets a version marker and a reader policy (read N and N−1, reject future majors with `FormatVersionError`):
  - GenBank features: add a `domainator_schema_version` record annotation or a COMMENT structured block, written by `write_genbank`.
  - JSON outputs (`summary_report`, `matrix_report`, `ssn_navigator`, `enum_report` NDJSON): add a top-level `schema_version`.
  - Matrix (Phase 0), `.dsnv` (already versioned), `.didx` (already versioned).
- Document the undocumented formats: XGMML attributes, Newick naming, `hmmer_compare` TSV, and the color table. Add the undocumented qualifiers (`accession`, `Color` vs `color`, pick one; `source_contig`, `threedi`, `gene_id`) and decide what `find_features` writes for `rstart/rend/rlen`, which is `'.'` today.
- Resolve the coordinate off-by-one TODOs before freezing the schema (`domainate.py:222-225,252-255`, `utils.py:1302,1716`, `enum_report.py:870`). Add explicit coordinate tests for both strands and for compound/circular locations.
- Taxonomy filtering vs database size Z: there are known-wrong TODOs at `domain_search.py:365,535` and `domainate.py:1319`. Fix them or document them before 1.0, because they affect reported e-values.

## Phase 7: Performance

- `write_records_with_limit` serializes every record twice (`output_guardrails.py:158-189`), and the guardrail is on by default. Measure with a counting writer instead.
- `write_genbank` mutates its input (id/name swap, adds qualifiers, sorts features), and `preserve_original` only shallow-copies (`utils.py:953-1012`). Make it non-mutating, and set up the warning filters once.
- pyhmmer: build `OptimizedProfile`s once per database and one `DigitalSequenceBlock` per batch (`domainate.py:479-486`); reuse the `pyrodigal.GeneFinder` (`domainate.py:1095`).
- Quadratic spots:
  - `features.index(feature)` inside a loop (`utils.py:774`);
  - `slice_record_from_location` scans all features for every slice (`utils.py:1950+`); use a sorted interval index;
  - `seq_dist` parses its input three times (`seq_dist.py:562-799`);
  - the recursive Newick builder (`build_tree.py:64,95`).
- Add a small benchmark script (`scripts/bench/`) to catch regressions: parse, domainate, domain_search on `Staph_phages.gb`.

## Phase 8: Packaging and release engineering

- **Dependencies:**
  - Drop Flask, Werkzeug and jsonschema (the GUI moves out).
  - Declare what the code imports directly: numpy, matplotlib, scikit-learn, numba.
  - Put umap-learn behind a `[projection]` extra with a clear error if it's missing. Drop bashplotlib (unmaintained), or move it to an extra.
  - Add `[structure]` (torch/esmologs) and `[dev]` (pytest, coverage, playwright, shtab) extras.
  - Remove the duplicate `core` extra.
  - Revisit the tight pins: pyhmmer `<0.13` (make `<1` and test), `pytest-datadir~=1.4.1`.
- **Metadata:**
  - Require `setuptools>=77` (needed for the PEP 639 `license` string).
  - Declare `license = "MIT AND BSD-3-Clause"` with `license-files`, and ship `Bio/LICENSE.rst`; it is currently missing from wheels, which is a compliance issue.
  - Add classifiers, keywords, and URLs (Homepage, Documentation, Changelog, citation DOI).
  - Build universal2 or x86_64 macOS wheels.
- **Python 3.15 support** (checked against PyPI and conda-forge on 2026-10-21; it looks achievable):

  | Dependency | 3.15 status |
  |---|---|
  | CPython 3.15.0 | on conda-forge |
  | numpy 2.5, scipy 1.18, pandas 3.0, scikit-learn 1.9, matplotlib 3.11, numba 0.68, llvmlite 0.50 | cp315 wheels. numba has no macOS x86_64 build |
  | pyhmmer 0.12.3, pyinfernal 0.1.1, psutil | `cp312-abi3`/`cp36-abi3` wheels install on 3.15. pyinfernal ships an abi3 wheel only for Linux x86_64, as it already does for 3.12+ |
  | umap-learn, pynndescent, seaborn, jsonargparse, tqdm, requests | pure Python |
  | PyO3 0.29 (`_gbfast`) | knows 3.15 (it has `MINIMUM_SUPPORTED_VERSION_ABI3T = 3.15`) |
  | **h5py 3.16** | no cp315 wheel on PyPI (a pip install would build from sdist against system HDF5), but conda-forge has py315 builds |
  | **pyrodigal 3.7.1** | no cp315 or abi3 wheel, so it builds from the Cython sdist. Verify that it compiles, or ask upstream for wheels |
  | **torch** | no cp315 wheel. It affects only the optional `[structure]` extra (ESM2/foldseek in domainate), which must stay optional and fail with a clear message |
  | **bashplotlib** | sdist only and unmaintained. Already slated for removal |
  | bioconda pyhmmer | builds stop at py313, so the conda env pip-installs the abi3 wheel (or request a rebuild) |

  **Steps:**
  1. Set `requires-python = ">=3.12"` (decided), with no upper cap. Raise dependency floors to match (numpy ≥2, etc.) and update `conda_env*.yml` and the `.def` files.
  2. Upgrade cibuildwheel to 3.x and build cp312–cp315 wheels. Consider building `_gbfast` as abi3 (PyO3's `abi3-py312` feature): one wheel per platform, and it works on future Python versions without rebuilding.
  3. Create a conda env on Python 3.15 and run the full suite. Expect to fix:
     - **UTF-8 mode is the default in 3.15** (PEP 686). Audit `open()` calls with no `encoding=`; mostly harmless, but test outputs and the Windows-agnostic code could change.
     - Anything removed in 3.14/3.15 that the vendored `Bio/` fork still uses (the fork goes away in Phase 5).
     - numba's eager `@jit(cache=True)` in `hmmer_search.py:105-179`. Make it lazy or optional (Phase 4) so a missing numba build can't break an unrelated tool.
  4. The multiprocessing start method is already pinned (`DOMAINATOR_MP_START_METHOD`, default spawn), so the 3.14 change to `forkserver` doesn't affect it.
- **Environments:**
  - `conda_env.yml`: stop doing an editable `.[test]` install for end users. Pin TMbed and move it to the structure extra. Add `infernal` (needed for `cmbuild`), drop `hmmer`, which is unused. Add usearch to arm64 or make it optional.
  - Pin esmologs in the `.def` files, and remove the duplicated verification blocks.
- **Release:**
  - A `CHANGELOG.md` (Keep a Changelog), `MIGRATION.md` (0.x→1.0) and `CITATION.cff`.
  - Single-source the version so the README doesn't hard-code it; render it in docs instead.
  - Tag every release.
  - Add a PyPI trusted-publishing job to `wheels.yml`, and a bioconda recipe.
- **External binaries:** one documented table (tool → binary → which commands need it → install hint) in the docs and in `domainator doctor`, a new subcommand that reports which optional binaries and extras are available.

## Phase 9: GUI split and the tool manifest

**Why the schemas need hand edits today.** I regenerated all 40 schemas and diffed them against the committed versions:
- 608 parameters received about 1,000 hand edits. The largest groups: 197 moves to advanced, 127 help rewrites, 73 `file_types` additions, 48 relabels, 49 `required` changes driven by stdin/stdout, 47 input-vs-output type fixes, 31 hidden parameters, 28 numeric bounds, and 18 choice reorders caused by Python `set` choices.
- The schema format has no field for several things: conditional visibility, mutually exclusive groups, at-least-one-output rules, directory vs. file, and stdin/stdout.
- Drift has already broken things:
  - ssn_navigator `id` vs. `cluster_id`;
  - transform_matrix's default `efi_score`, which the GUI silently submits;
  - compare_contigs `ji`/`ai` typed as strings;
  - advanced `file` parameters that `tool_executor.py:183` never resolves.

**Design: Domainator owns a versioned tool manifest, and the GUI only renders it.**

1. **Metadata channel** (new `src/domainator/cli/meta.py`).
   - The Phase 2d helpers record a `Param` on the argparse action they return, a pattern I checked works with jsonargparse 4.49, or carry it through `Annotated[str, Meta(...)]`, which still parses as a plain `str`.
   - Each `Param` holds:
     - `kind` (path, value, flag, `ordered_list`);
     - `role` (input, output, cache, dir) and `formats` (keys into a formats registry built from `utils.EXTENSION_TO_TYPE` plus hdf5/hmm/cm/tsv/…);
     - `multiple`, `stdin`, `stdout`;
     - `one_of` (output sets);
     - `level` (basic, advanced, hidden), `label`, `short_help`;
     - `depends_on` / `conflicts`;
     - `default_dynamic` (environment-dependent defaults such as the taxonomy path, which must not be baked in);
     - `min` and `max`, read from jsonargparse restricted types (`PositiveInt`, `ClosedUnitInterval`).
   - Use `Literal`/sorted choices everywhere so the output is deterministic.
2. **Tool metadata:** `TOOL = ToolInfo(id, title, category, summary, priority, requires=[binaries], extras=[...])` per module, plus tool-level `constraints` (`at_least_one`, mutually exclusive groups, which are read from the parser).
3. **enum_report and hmmer_report columns:**
   - Replace `DynamicArg` and the roughly 34 `append_const` flags that share `output_cols` with one declared `ordered_list` parameter. Its catalog of column kinds gives each column typed, named arguments (e.g. `taxid(rank: str+)`, `qualifier(feature, qualifier)`, `net_charge(pH: float)`, `append(name, type, value)`).
   - The CLI keeps its ordered-flags syntax, generated from the catalog.
   - The config encoding (`[key | [key, [args…]]]`) already round-trips through `--config` (checked).
   - This also fixes the reused long-flag meanings from Phase 3 (`--domains`, `--score`, `--strand`, … as columns).
4. **Export:**
   - `domainator describe [--tool X] --json` emits `{manifest_version, domainator_version, formats, categories, invocation: {config_flag: "--config", keys: "dest"}, tools: […]}`.
   - The invocation contract is that config keys equal argparse dest names. Use `capture_parser` only as a fallback.
   - Document the manifest schema in `docs/file_formats.md` and version it like the other formats.
5. **Enforcement (CI):**
   - Every path-like argument declares a role.
   - Manifest defaults round-trip through `parse_args(["--config", json])`.
   - The manifest snapshot is deterministic.
   - `allow_abbrev=False` is set.
   This class of drift then becomes a test failure.
6. **Removal from this repo** (in 1.0rc1):
   - Delete `src/domainator/server/`, `scripts/server/`, `test/server/`, `docs/server/`, the server schemas, and the `domainator_server` entry point and Flask deps.
   - The new GUI project depends on `domainator>=1.0` and reads `domainator describe`.
   - GUI-side policy stays in the GUI: category colors, presets that replace GUI-only default overrides (domainate max_overlap 0.6, build_tree quiet), sandboxing, resource limits, and whether binary-path parameters are exposed.
7. **Hand the new GUI project the server security findings**, so they are not repeated:
   - Unknown parameters pass straight through to the CLI (`tool_executor.py:194-198`), which together with exposed binary-path parameters allows arbitrary command execution.
   - Any host path is accepted as a file parameter (`:362-365`).
   - Uploaded HTML is served inline from the server's own origin (stored XSS).
   - No auth, no Host/Origin check, and `--debug` turns on the Werkzeug debugger.
   - `file_id` is not validated (`file_manager.py:141`).
   - Filter unknown keys using the manifest and allowlist `role=input` paths to uploads.

## Phase 10: SSN viewer source layout

- `ssn_viewer_html.py` (13k lines, 600 KB) is HTML, CSS and JavaScript held in Python strings, with 1,323 doubled braces. Move the assets to `src/domainator/ssn_viewer/{viewer.html,viewer.css,js/*.js}` as package data, loaded with `importlib.resources`.
- Inject the runtime values in one `<script type="application/json" id="dsnv-config">` block, and inline the assets at render time so the output stays a single self-contained file.
- Add ESLint/Prettier and Node unit tests for the pure layout math. Keep the Playwright suite.
- The other HTML reports load d3, plotly and tabulator from CDNs with floating versions and no SRI, so they break offline:
  - `plot_contigs.py:269`
  - `matrix_report.py:378`
  - `summary_report.py:221`
  - `enum_report.py:274`
  Vendor these libraries inline, as the viewer does, or at least pin exact versions and add SRI hashes.

## Phase 11: Docs

- Generate the **CLI reference** from the manifest (one page per command), so it can't drift. Add a CI check that every `domainator <tool> --flag` written in the docs and README parses; the audit found several stale ones (`README.md:44,112,221,274`, `docs/ssn_viewer.md:53`, `docs/structure_search.md:91`, `docs/server/adding_tools.md:7`).
- Build the **API reference** with mkdocs-material + mkdocstrings, from `__all__` docstrings, and publish it via GitHub Pages or Read the Docs. Bring docstrings and type hints for the public API to 100%; today 57% of functions have a docstring and 45% have any type hint.
- **README:** make Quickstart and Installation agree (bioconda vs zip), list all commands (6 are missing today), and add tutorial snippets for the commands that have no usage example (about 17).

## Phase 12: Tests

- Stop copying all of `test/data` (12 MB) into every test's `tmp_path`. `shared_datadir` does this about 1,600 times, so a full run writes about 11 GB of temp files; on a 16 GB tmpfs `/tmp`, two runs filled it. Phase 0 added `tmp_path_retention_policy = "failed"`, which deletes them as tests pass. The real fix is pytest-datadir ≥1.5's `lazy_shared_datadir` (copies only the files a test opens), or read-only fixtures opened in place.
- Test the public library API directly. Tests call `X.main(` 564 times but call the core functions almost never.
- Fill the gaps:
  - no tests: `cytoscape`, `edge_types`, `matrix_output`;
  - only indirect coverage: `output_guardrails`, `genbank_to_fasta`, `ssn_navigator`, `ssn_bundle`;
  - thin: `ssn_edges`, `summary_report`, `plot_contigs`.
- Add tests for:
  - the dispatcher (every registered tool resolves, `--help` works, and help stays under the import-time budget);
  - completion-script generation;
  - manifest snapshots;
  - Rust-vs-Python parser parity;
  - round-tripping format versions (old fixtures still read).

---

## Suggested release sequence
1. **0.13**: Phases 0, 1, 2a/2b (the dispatcher alongside the old scripts), and the import-hygiene part of 4. Non-breaking, and gives early feedback.
2. **0.14**: Phase 5, the move to upstream Biopython. Doing it before the 1.0 renames isolates its output changes from the interface changes.
3. **1.0.0rc1**: Phase 3 renames, the Phase 4 API, Phase 6 versioning, Phase 8 packaging, the Phase 9 manifest, and removal of the server and the `.py` scripts. Ship MIGRATION.md.
4. **1.0.0**: after the rc has been used on real workflows. Phases 7, 10, 11 and 12 can land before or after 1.0, since they don't change interfaces.

## Verification
- `pytest test` green in CI on all supported Python versions; binary-dependent jobs run in the conda job.
- `domainator --help`, `domainator <tool> --help` for every registered tool, and `domainator completion bash | bash -n`.
- An import-time budget test: `python -X importtime -c "import domainator.<tool>"` stays under its threshold.
- A manifest snapshot test, plus a test that regenerating the manifest is deterministic (choices sorted).
- A docs check: every command line in `README.md`/`docs/` parses.
- Parser parity: Rust lean vs Biopython on all fixtures.
- Build the wheel and sdist, install in a clean venv, and run `domainator doctor` and a smoke workflow (domainate → select_by_cds → enum_report).
