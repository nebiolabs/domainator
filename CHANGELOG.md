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
- `--palette_offset` on build_ssn, build_projection and plot_contigs, and a matching start
  offset in the SSN viewer's color picker,
  for coloring a second column without reusing the palette colors a first one used.
- A Zoom box in the SSN viewer's View Settings that shows and sets the zoom, so a zoom can
  be copied from one viewer into another.
- `domainator describe`: a JSON description of every command and its options (the tool
  manifest, docs/tool_manifest.md), for GUIs and workflow systems.
- A command-line reference generated from the parsers (docs/cli/), and a documentation site
  with the Python API reference, published at https://nebiolabs.github.io/domainator/.
- Python extras: `projection` (UMAP), `structure` (torch), `dev`.
- `scripts/bench/bench.py`, a benchmark of common steps.
- The `domain_search` function searches SeqRecords as well as file partitions, reading a stream
  of records only a few batches ahead of the search; it has a `batch_size` parameter.
- On free-threaded Python with the GIL off, parallel work runs in threads instead of processes
  (`DOMAINATOR_MP_START_METHOD` overrides this). Domainator builds on free-threaded Python.

### Changed
- Commands run as `domainator <command>`; the `<command>.py` scripts are gone.
- Every output is an explicit `--<content>_out` option, and nothing reads standard input or
  writes standard output unless given `-`.
- Options with the same meaning have the same name everywhere (see MIGRATION.md).
- enum_report and hmmer_report choose their columns with one ordered `--columns` list
  (`--columns domains length taxid:genus`) instead of one flag per column.
- Domainator uses upstream Biopython (1.85 to 1.88) instead of a bundled fork.
- Python 3.12 or newer is required; wheels use the stable ABI, so one wheel serves every
  Python from 3.12 on (free-threaded Python builds its own).
- Matrix files are version 2.0 (symmetric matrices are stored as their upper triangle).
- The conda environment no longer installs TMbed and PyTorch; the README's "Optional: TMbed"
  section says how to add them for `find_features --tmbed`.
- The command functions take only their main input positionally and are importable as
  `from domainator.<command> import <command>`; their keyword parameters are named like the
  command's options and default like them (domainate's e-value 0.001, not 10; `cpu=0`, all
  cores). Every one has a docstring and type hints.
- Faster: neighborhood extraction (select_by_cds, domain_search) and extract_domains look only
  at nearby features; HMM databases are prepared once per search, not once per batch.
- enum_report and hmmer_report `--html_out` pages work offline: the table (sort, per-column
  filters, paging with first/last and go-to-page, movable and resizable columns, a right-click
  "Copy cell value", CSV/TSV/JSON downloads) is Domainator's own JavaScript
  instead of Tabulator from a CDN. The page carries the same rows as `--json_out`, so the
  downloads keep full precision; the table still shows floats to one decimal place.
- summary_report's `--html_out` page works offline and no longer needs matplotlib or seaborn:
  the contig length histogram (with hover counts and SVG/PNG downloads) and the zoomable
  taxonomy icicle are drawn by Domainator's own JavaScript instead of a seaborn PNG and d3
  from a CDN, and the domain frequency and co-occurrence tables are sortable, filterable
  tables like enum_report's. The page carries the `--json_out` document, so it shows the
  same numbers. The taxonomy tree's children are now in a fixed order (most contigs first).
- matrix_report's `--html_out` page works offline: its charts are Domainator's own JavaScript
  instead of Plotly and d3 from a CDN. The clusters-by-threshold and clusters-by-edges charts
  zoom by dragging across them; the split chart lines up with the threshold slider under it,
  and clicking it moves the threshold there; the arrow keys step the slider from stop to stop;
  and each chart, the 100-largest-clusters bubble chart included, downloads as SVG or PNG. The
  page carries the `--json_out` document, so its summary shows the same numbers.
- plot_contigs works offline: the plot is Domainator's own JavaScript instead of d3 from a CDN,
  with zoom buttons, SVG/PNG downloads and a details box that lists every value in full. Each
  part of a multi-part feature is now drawn where it is, not all at the feature's start. A
  contig's details list a fixed set of record annotations (references and structured comments
  are left out), and its metadata as the text in the file. `--color_by` needs `--metadata`
  and must name one of its columns; a contig with no value is drawn gray instead of failing.

### Removed
- `partition_seqfile` and `partition_seqids` (domainate and domain_search partition their
  input themselves), and `--offset`/`--recs_to_read` (index the database with
  `domainator format_db --index` instead).
- The `bashplotlib` dependency; `umap-learn` is now an optional extra.
- The matplotlib and seaborn dependencies. No HTML page needs them; `build_projection
  --png_out` still does and takes them from the `projection` extra
  (`pip install "domainator[projection]"`).
- The bundled web GUI (`domainator_server`) and its Flask dependency; a separate GUI project is
  being built on `domainator describe`.

### Fixed
- A taxid missing from the NCBI taxonomy database (one newer than the database) crashed
  summary_report and enum_report's taxonomy columns. It is now treated as unknown, with
  blank names and ranks and one warning per taxid suggesting `--taxonomy_update`.
- summary_report listed domains with equal counts in an arbitrary order once there were more
  than 16 domains; they are now in the order they were first seen.
- Many bugs found in the 1.0 review, including: `-e` cutoffs that dropped hits equal to the
  cutoff; XGMML values escaped twice; Newick branch lengths rounded to two decimals;
  build_tree crashing on deep trees; extract_peptides `-e` parsed as text; compressed output
  ignored by some commands; a missing-qualifier error that printed a record to standard output;
  matrix_report's kNN edge counts with long neighbor lists; build_projection with
  `--color_table_out` and no colors, and `--pca_args` on sparse matrices.
- seq_dist's hmmer modes read the target profiles again for every query profile, and
  hmmer_search and hmmer_compare sent every reference profile with every query; each worker
  now gets them once.
- domain_search's `--batch_size` set the partition size but not how many CDSs were searched at
  once.
- Output files were created readable only by their owner (mode 0600); they now get the usual
  permissions from the umask.
- format_db `--shards` failed on compressed (gzip or BGZF) input, and format_db left an empty
  temporary file next to the output whenever its input was compressed.
- db_download `--db uniprot` and `--db uniprot_gb` wrote only TrEMBL: the TrEMBL download
  overwrote the SwissProt records already written. `--max_records` now caps the total across
  both rather than each.
