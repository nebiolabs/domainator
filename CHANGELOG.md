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
- Python extras: `projection` (UMAP), `structure` (torch), `dev`.
- `scripts/bench/bench.py`, a benchmark of common steps.

### Changed
- Commands run as `domainator <command>`; the `<command>.py` scripts are gone.
- Every output is an explicit `--<content>_out` option, and nothing reads standard input or
  writes standard output unless given `-`.
- Options with the same meaning have the same name everywhere (see MIGRATION.md).
- Domainator uses upstream Biopython (1.85 to 1.88) instead of a bundled fork.
- Python 3.12 or newer is required; wheels use the stable ABI, so one wheel serves every
  Python from 3.12 on.
- Matrix files are version 2.0 (symmetric matrices are stored as their upper triangle).
- The command functions take only their main input positionally and are importable as
  `from domainator.<command> import <command>`.
- Faster: neighborhood extraction (select_by_cds, domain_search) and extract_domains look only
  at nearby features; HMM databases are prepared once per search, not once per batch.

### Removed
- `partition_seqfile` and `partition_seqids` (domainate and domain_search partition their
  input themselves), and `--offset`/`--recs_to_read` (index the database with
  `domainator format_db --index` instead).
- The `bashplotlib` dependency; `umap-learn` is now an optional extra.

### Fixed
- Many bugs found in the 1.0 review, including: `-e` cutoffs that dropped hits equal to the
  cutoff; XGMML values escaped twice; Newick branch lengths rounded to two decimals;
  build_tree crashing on deep trees; extract_peptides `-e` parsed as text; compressed output
  ignored by some commands; a missing-qualifier error that printed a record to standard output.
