# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). From 1.0 on,
Domainator follows [Semantic Versioning](https://semver.org/): what stays stable is listed in
[docs/stability.md](docs/stability.md).

## [Unreleased]

The 1.0 release: one `domainator` command, consistent options, a documented Python API and
documented, versioned file formats. Many options were renamed; [MIGRATION.md](MIGRATION.md)
lists every change from 0.x.

Extraction options and matching Python keywords now use `keep_parent_id` instead
of `keep_name` on all three extract tools. `extract_unannotated` also renames
`min_length` to `min_region_length` and `min_score` to `min_annotation_score`.
Defaults and extraction behavior are unchanged; old names are not aliases.

### Added
- `extract_peptides --domain_expr` for Boolean domain selection and
  `select_by_cds --name_field` for explicit CDS name resolution, with matching
  keyword-only library parameters. Existing default name precedence is preserved.
- Shared taxonomy filters on extract_domains, extract_peptides, extract_unannotated,
  and select_by_cds, consistent with select_by_contig and domain_search. CDS/domain
  extraction filters original anchors, not expanded neighborhoods.
- `--taxonomy_unknown {exclude,include,evaluate,error}` on these selection, extraction,
  and sequence-search tools. The default excludes unresolved taxonomy even for
  negated criteria such as `~2`; empty-lineage evaluation is an explicit opt-in.
- `--has_search_hit` on select_by_cds, select_by_contig, extract_peptides, extract_domains,
  and rotate_contigs: an independent search-hit selection criterion that bypasses domain
  names, expressions, and database restrictions, honors quality cutoffs, and participates
  in `--invert` where available. `--annotation_type {domainator,domain_search,all}` now
  consistently selects annotation sources rather than implicitly adding search hits.
- The `domainator <command>` dispatcher, with `domainator --help`, `domainator completion
  {bash,zsh,tcsh}` (shell tab completion) and `domainator doctor` (which optional programs
  and packages are installed).
- `-v/--verbose` and `-q/--quiet` on every command; errors you can fix are one-line messages.
- `--fasta_out` on every command that writes GenBank, alongside `-o/--genbank_out`.
- `docs/stability.md`, and a reference for every GenBank qualifier Domainator writes in
  `docs/file_formats.md`.
- A `schema_version` key in the `--json_out` of summary_report, matrix_report and ssn_navigator.
- hmmer_compare `--alignments_out`.
- The Domainator Viewer: build_ssn, build_tree and build_projection now write one HTML page
  (`--html_out`) and one bundle format (`--dvb_out`, `.dvb`, docs/file_formats.md) for networks,
  trees and projections alike, with the same selection-driven metadata table, editing, colors,
  presets, column charts and saved sessions for all three. `domainator viewer` writes the page
  empty, or opening a bundle (`-i`).
- build_ssn `--html_out`, `--dvb_out`, `--merge_impact_metric`, `--label_by`, `--categorical` and
  `--name`, replacing build_ssn_viewer; `--xgmml_out` is no longer required.
- build_tree `--html_out` and `--dvb_out`: the tree in the Domainator Viewer, with a multiple
  sequence alignment beside it (`--msa`, aligned FASTA); clade selection, collapsing,
  ladderizing, a cut slider that numbers the clades below it, and Newick and FASTA exports of
  the selection. `--newick` reads a tree built elsewhere instead of a matrix. Also
  `--color_by`, `--color_table`, `--color_table_out`, `--palette_offset`, `--max_color_groups`,
  `--label_by`, `--categorical`, `--name`, `--subset`, `--subset_file` and `--max_output_gb`.
- build_projection `--dvb_out`, `--label_by`, `--categorical` and `--name`.
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
- Peptide and CDS-neighborhood selection now invert only the combined OR selectors,
  never strand or other independent restrictions. With no selectors, both library
  and CLI calls select all eligible CDSs; inversion warns and is ignored.
- domain_search and select_by_cds replace deduplication switches with the single-choice
  `--dedup {none,id,sequence,sequence_rc}` policy. Their defaults remain
  none and id, respectively. Ranked searches retain the top `--max_hits` winners
  incrementally in each worker and the global merge, replacing weaker copies
  without buffering all candidates. Deduplication relies on Domainator's
  invariant that the same record ID identifies the same sequence.
  The pre-release multi-value syntax is removed; configuration and library APIs
  also use one policy string.
- Commands run as `domainator <command>`; the `<command>.py` scripts are gone.
- build_projection `--html_out` writes the Domainator Viewer, which selects points into a metadata
  table; its on-screen legend is now the viewer's legend export.
- build_tree's xgmml output is laid out as a rectangular tree, colored with `--color_by`; its
  internal nodes are named by their place in the tree in preorder.
- The SSN viewer is now the Domainator Viewer, and its sources moved from `ssn_viewer/` to
  `viewer_app/` (`ssn_viewer_html.py` is now `viewer_html.py`).
- Networks for the viewer are `.dvb` bundles, and the viewer saves sessions and extractions as
  `.dvb`. The viewer, `domainator viewer -i` and ssn_navigator still read `.dsnv` bundles.
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
  the contig length histogram (with hover counts, wheel zoom and drag to pan, a bin width box,
  and SVG/PNG downloads) and the zoomable
  taxonomy icicle are drawn by Domainator's own JavaScript instead of a seaborn PNG and d3
  from a CDN, and the domain frequency and co-occurrence tables are sortable, filterable
  tables like enum_report's. The page carries the `--json_out` document, so it shows the
  same numbers. The taxonomy tree's children are now in a fixed order (most contigs first).
- summary_report's text report (`-o`) lays out its tables as plain aligned columns instead of
  pandas' `to_string`, so the layout no longer changes with the pandas version: numbers are
  right-aligned, `avg score` has one decimal place, and the co-occurrence and taxonomy tables
  lose pandas' row numbers.
- matrix_report's `--html_out` page works offline: its charts are Domainator's own JavaScript
  instead of Plotly and d3 from a CDN. The clusters-by-threshold and clusters-by-edges charts
  zoom by dragging across them; the split chart lines up with the threshold slider under it,
  and clicking it moves the threshold there; the arrow keys step the slider from stop to stop;
  and each chart, the 100-largest-clusters bubble chart included, downloads as SVG or PNG. The
  page carries the `--json_out` document, so its summary shows the same numbers, in a compact
  block of figures beside each other. Dragging across a chart never selects the page's text.
- plot_contigs works offline: the plot is Domainator's own JavaScript instead of d3 from a CDN,
  with zoom buttons, SVG/PNG downloads and a details box that lists every value in full. Each
  part of a multi-part feature is now drawn where it is, not all at the feature's start. A
  contig's details list a fixed set of record annotations (references and structured comments
  are left out), and its metadata as the text in the file. `--color_by` needs `--metadata`
  and must name one of its columns; a contig with no value is drawn gray instead of failing.
- build_projection `--html_out` replaces `--png_out`: an offline, interactive plot of the
  projection, colored as the xgmml is, with each point's metadata on hover, legend entries that
  hide and show their group, zoom and pan, a search by name, X/Y, X/Z and Y/Z views of a 3D
  projection, and SVG/PNG downloads of the current view. A value missing from `--color_table`
  is gray instead of failing.

### Removed
- Taxonomy filters, `--hits_only`, and `--no_annotations` from domainate;
  structure_domainate also drops the latter two. Their public Python functions
  drop the corresponding keywords. Annotation tools retain every input record
  and annotate hits; use the search tools for hit selection.
- `partition_seqfile` and `partition_seqids` (domainate and domain_search partition their
  input themselves), and `--offset`/`--recs_to_read` (index the database with
  `domainator format_db --index` instead).
- The `bashplotlib` dependency; `umap-learn` is now an optional extra.
- The matplotlib and seaborn dependencies: no command needs them. `build_projection
  --png_out` is gone; its `--html_out` page downloads the plot as PNG or SVG.
- The bundled web GUI (`domainator_server`) and its Flask dependency; a separate GUI project is
  being built on `domainator describe`.
- build_ssn_viewer: build_ssn `--html_out` and `--dvb_out` write the viewer and its bundle, and
  `domainator viewer` writes the page on its own (see MIGRATION.md).

### Fixed
- Strand-independent exact deduplication no longer collapses complement-only or
  reverse-only nucleotide sequences; only reverse-complement copies are equivalent.
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
