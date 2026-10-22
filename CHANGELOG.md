# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). From 1.0 on,
Domainator follows [Semantic Versioning](https://semver.org/): what stays stable is listed in
[docs/stability.md](docs/stability.md).

## [Unreleased]

The 1.0 release: one `domainator` command, consistent options, a documented Python API and
documented, versioned file formats. Many options were renamed; [MIGRATION.md](MIGRATION.md)
lists every change from 0.x.

### Added
- The `domainator <command>` dispatcher, with `domainator --help`, `domainator completion
  {bash,zsh,tcsh}` (shell tab completion) and `domainator doctor` (which optional programs
  and packages are installed).
- `-v/--verbose` and `-q/--quiet` on every command; errors you can fix are one-line messages.
- `--fasta_out` on every command that writes GenBank, alongside `-o/--genbank_out`.
- `docs/stability.md`, and a reference for every GenBank qualifier Domainator writes in
  `docs/file_formats.md`.
- A `schema_version` key in the `--json_out` of summary_report, matrix_report and ssn_navigator.
- hmmer_compare `--alignments_out`.
- `domainator describe`: a JSON description of every command and its options (the tool
  manifest, docs/tool_manifest.md), for GUIs and workflow systems.
- A command-line reference generated from the parsers (docs/cli/), and a documentation site
  with the Python API reference, published at https://nebiolabs.github.io/domainator/.
- Python extras: `projection` (UMAP), `structure` (torch), `dev`.
- `scripts/bench/bench.py`, a benchmark of common steps.

### Changed
- Commands run as `domainator <command>`; the `<command>.py` scripts are gone.
- Every output is an explicit `--<content>_out` option, and nothing reads standard input or
  writes standard output unless given `-`.
- Options with the same meaning have the same name everywhere (see MIGRATION.md).
- enum_report and hmmer_report choose their columns with one ordered `--columns` list
  (`--columns domains length taxid:genus`) instead of one flag per column.
- Domainator uses upstream Biopython (1.85 to 1.88) instead of a bundled fork.
- Python 3.12 or newer is required; wheels use the stable ABI, so one wheel serves every
  Python from 3.12 on.
- Matrix files are version 2.0 (symmetric matrices are stored as their upper triangle).
- The conda environment no longer installs TMbed and PyTorch; the README's "Optional: TMbed"
  section says how to add them for `find_features --tmbed`.
- The command functions take only their main input positionally and are importable as
  `from domainator.<command> import <command>`; their keyword parameters are named like the
  command's options and default like them (domainate's e-value 0.001, not 10; `cpu=0`, all
  cores). Every one has a docstring and type hints.
- Faster: neighborhood extraction (select_by_cds, domain_search) and extract_domains look only
  at nearby features; HMM databases are prepared once per search, not once per batch.

### Removed
- `partition_seqfile` and `partition_seqids` (domainate and domain_search partition their
  input themselves), and `--offset`/`--recs_to_read` (index the database with
  `domainator format_db --index` instead).
- The `bashplotlib` dependency; `umap-learn` is now an optional extra.
- The bundled web GUI (`domainator_server`) and its Flask dependency; a separate GUI project is
  being built on `domainator describe`.

### Fixed
- Many bugs found in the 1.0 review, including: `-e` cutoffs that dropped hits equal to the
  cutoff; XGMML values escaped twice; Newick branch lengths rounded to two decimals;
  build_tree crashing on deep trees; extract_peptides `-e` parsed as text; compressed output
  ignored by some commands; a missing-qualifier error that printed a record to standard output;
  matrix_report's kNN edge counts with long neighbor lists; build_projection with
  `--color_table_out` and no colors, and `--pca_args` on sparse matrices.
