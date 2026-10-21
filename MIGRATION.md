# Migrating from Domainator 0.x to 1.0

Domainator 1.0 makes the command line consistent across all programs. Old option names are
**not** kept as aliases: a script that uses one fails with an "unrecognized arguments" error
that names it, so nothing changes silently. This guide lists every change.

## Running programs

| 0.x | 1.0 |
|---|---|
| `domainate.py -i in.gb ...` | `domainator domainate -i in.gb ...` |
| `domainator_db_download.py ...` | `domainator db_download ...` |
| `domainator_format_db.py ...` | `domainator format_db ...` |
| `partition_seqfile.py`, `partition_seqids.py` | removed; index databases with `domainator format_db --index` instead |
| `domainator_server` | moved to a separate project |

`domainator --help` lists all commands. `domainator completion bash` prints a tab-completion
script (see the README).

Options must be spelled out in full. 0.x accepted any unambiguous prefix (`--alignment` for
`--alignments`); 1.0 does not.

Every program takes `-v/--verbose` (debug output) and `-q/--quiet` (errors only: no progress
messages or warnings). `build_tree --quiet` is now this shared option.

Problems with inputs, options, or external programs print one line,
`domainator <command>: error: ...`, and exit with status 1; usage errors exit with status 2.
`-v` shows the traceback. In Python, these errors are subclasses of
`domainator.errors.DomainatorError` (and of `ValueError` or `RuntimeError`, as before). Tools that
used to `exit(1)` from library code (seq_dist's diamond errors) or raise `SystemExit` for
`--max_output_gb` now raise `ExternalToolError` or `OutputSizeLimitExceeded`.

## Inputs and outputs

**Nothing reads stdin or writes stdout unless asked to.** Use `-` as the path:
`-i -` reads stdin and `-o -` writes stdout. `/dev/stdin` and `/dev/stdout` work too. Only one
output can go to stdout at a time. A program given no output at all is an error.

```bash
# 0.x
domainate.py -i in.gb -r pfam.hmm | select_by_cds.py --domains CcdB > out.gb
# 1.0
domainator domainate -i in.gb -r pfam.hmm -o - | domainator select_by_cds -i - --domains CcdB -o out.gb
```

Programs that read their input more than once refuse `-i -`: domainate, domain_search,
seq_dist, kmer_dist, deduplicate_genbank, split_by_cluster, and format_db. So do the programs
whose input is not a stream: the structure programs (structure files or directories) and the
matrix and bundle readers (build_projection, build_ssn, build_tree, build_ssn_viewer,
matrix_report, transform_matrix, ssn_navigator).

**Every output file has its own option, ending in `_out`.** No input option ends in `_out`.
Each program's main output also has the short option `-o`:

| Output | Option |
|---|---|
| GenBank | `-o/--genbank_out` |
| FASTA | `--fasta_out` (in every program that writes GenBank; both can be written at once) |
| HDF5 matrix | `-o/--hdf_out` |
| tab-separated matrix | `--dense_text_out` |
| TSV table | `--tsv_out` |
| HTML / JSON / text report | `--html_out` / `--json_out` / `--text_out` |
| HMM / CM profiles | `-o/--hmm_out` / `-o/--cm_out` |
| Cytoscape XGMML / Newick tree / SVG / PNG | `--xgmml_out` / `--newick_out` / `--svg_out` / `--png_out` |
| log file | `--log_out` |

`--fasta_out` used to be a switch that made `-o` write FASTA instead of GenBank. It is now a
path:

```bash
# 0.x
extract_peptides.py -i in.gb -o peptides.fasta --fasta_out
# 1.0
domainator extract_peptides -i in.gb --fasta_out peptides.fasta
```

Output ending in `.gz` is compressed for every GenBank, FASTA, and TSV output.

## Options renamed in every program that has them

| 0.x | 1.0 | Notes |
|---|---|---|
| `-r/--references` | `-r/--reference` | |
| `--lb` (matrix and network programs) | `--threshold` | keep values strictly above this |
| `--score_cutoff`, `--score` | `--min_score` | |
| `--search_hits` | `--annotation_type domain_search` | `--annotation_type {domainator,domain_search}` |
| `--keep_db` | `--db_out` | |
| `--hits_tsv` | `--hits_tsv_out` | |
| `--input_type` / `--reference_type` (seq_dist, kmer_dist) | `--input_format` / `--reference_format` | |
| `--pad_where` (select_by_cds, domain_search), `--where` (rotate_contigs) | `--anchor` | |
| `--max_domains` (domainate, structure_domainate) | `--max_hits_per_cds` | unlimited by default; `0` is an error |
| `--aligner_arg` | `--foldseek_args`, `--reseek_args` | one quoted string: `--foldseek_args="--prefilter-mode 2"` |
| `--params` | `--diamond_args` (seq_dist), `--<algorithm>_args` (deduplicate_genbank, build_projection) | one quoted string |
| `--fasta_type` (domainate, domain_search, clean_sequences, deduplicate_genbank, find_features, kmer_dist) | `--alphabet {nucleotide,protein}` | the molecule type of FASTA input |
| `-e/--evalue` in filtering programs | unchanged, but the default is now "no cutoff" | was 1e8 in some programs and none in others |
| `-e/--evalue`, `--min_evalue` | unchanged | must be greater than 0 (`--min_evalue`: at least 0); checked when the command line is parsed |
| `--cpu` | unchanged | `0` (the default) now means one thread per physical core everywhere |

## Program-specific changes

Options not listed are unchanged.

| Program | 0.x | 1.0 |
|---|---|---|
| build_projection | `--png_out`, `--xgmml`, `--coords_out` | `--png_out`, `-o/--xgmml_out`, `--coords_out` |
| | `--3d` | `--three_d` |
| | `--params` | `--umap_args`, `--tsne_args`, `--pca_args`: `key=value` pairs, e.g. `--umap_args="n_neighbors=30 min_dist=0.1"` |
| build_ssn | `--xgmml`, `--cluster_tsv` | `-o/--xgmml_out`, `--cluster_tsv_out` |
| | `--lb` | `--threshold` |
| build_ssn_viewer | `-o/--output`, `--html` | `-o/--dsnv_out`, `--html_out` |
| build_tree | `--xgmml`, `--newick` | `--xgmml_out`, `-o/--newick_out` |
| clean_sequences | `-o/--output` + `--fasta_out` switch, `--log` | `-o/--genbank_out`, `--fasta_out`, `--log_out` |
| color_genbank | `-o/--output` | `-o/--genbank_out` |
| | `--mode` | `--target` |
| color_table_to_legend | `--svg` | `-o/--svg_out` |
| compare_contigs | `--hdf`, `--dense_text`, `-o/--output` (sorted GenBank) | `-o/--hdf_out`, `--dense_text_out`, `--genbank_out` |
| | `--lb` | `--threshold` |
| db_download | `-o/--output` | `--genbank_out` or `--fasta_out`, matching the format of `--db` |
| | `--skip_taxonomy_update` | `--no_taxonomy_update` |
| | `--num_recs` | `--max_records` |
| | `--success_rec_log` | `--success_log_out` |
| deduplicate_genbank | `-o/--output` + `--fasta_out` switch, `--log`, `--cluster_table` | `-o/--genbank_out`, `--fasta_out`, `--log_out`, `--cluster_tsv_out` |
| | `--id` | `--identity` |
| | `--params` | `--cdhit_args`, `--usearch_args`, `--diamond_args` |
| | `--bin_path` | `--cdhit_path`, `--cdhit_est_path`, `--usearch_path`, `--diamond_path` |
| domain_search | `-o/--output` | `-o/--genbank_out`, `--fasta_out` |
| | `-r/--references` | `-r/--reference` |
| | `--pad_where` | `--anchor` |
| | `--deduplicate` | `--deduplicate` / `--no_deduplicate` (still off by default) |
| | `--dedup_both_strands` | `--both_strands` |
| domainate | `-o/--output` | `-o/--genbank_out`, `--fasta_out` |
| | `-r/--references` | `-r/--reference` |
| | `--max_domains` | `--max_hits_per_cds` |
| | `--esm2_3Di_weights`, `--esm2_3Di_device` | `--esm2_3di_weights`, `--esm2_3di_device` |
| | `--offset`, `--recs_to_read` | removed |
| | | new: `--foldseek_path` |
| enum_report | `-o/--output`, `--html`, `--json` | `-o/--tsv_out`, `--html_out`, `--json_out` |
| extract_domains | `-o/--output` + `--fasta_out` switch | `-o/--genbank_out`, `--fasta_out` |
| | `--score` | `--min_score` |
| | `--search_hits` | `--annotation_type domain_search` |
| | `--pad_up`, `--pad_down` | `--extend_up`, `--extend_down`, and new `--extend_range` |
| extract_peptides | `-o/--output` + `--fasta_out` switch | `-o/--genbank_out`, `--fasta_out` |
| | `--cdss`, `--cdss_file` | `--cds`, `--cds_file` |
| | `--search_hits` | `--annotation_type domain_search` |
| extract_unannotated | `-o/--output` + `--fasta_out` switch | `-o/--genbank_out`, `--fasta_out` |
| | `--score` | `--min_score` |
| | `--lb` | `--min_length` |
| | `--search_hits` | `--annotation_type domain_search` |
| filter_domains | `-o/--output` | `-o/--genbank_out`, `--fasta_out` |
| | `-m/--max_overlap` | `--max_overlap` |
| | `--databases_filter` | `--databases` |
| | `--databases_keep` | `--keep_databases` |
| find_features | `-o/--output` | `-o/--genbank_out`, `--fasta_out` |
| | `--gpu_device` | `--device` |
| genbank_to_fasta | `-o/--output` | `-o/--fasta_out` |
| hmmer_build, hmmer_select, hmmer_search | `-o/--output` | `-o/--hmm_out` |
| hmmer_compare | `-o/--output` | `-o/--tsv_out` |
| | `--score_cutoff` | `--min_score` |
| hmmer_report | `-o/--output`, `--html`, `--json` | `-o/--tsv_out`, `--html_out`, `--json_out` |
| hmmer_search | `--score_cutoff` | `--min_score` |
| infernal_build | `-o/--output` | `-o/--cm_out` |
| | `--informat` | `--input_format` |
| kmer_dist | `--hdf`, `--dense_text` | `-o/--hdf_out`, `--dense_text_out` |
| | `--lb` | `--threshold` |
| matrix_report | `-o/--output`, `--html`, `--json` | `-o/--text_out`, `--html_out`, `--json_out` |
| plot_contigs | `--html` | `-o/--html_out` |
| | | new: `--color_table`, colors for `--color_by` values |
| rotate_contigs | `-o/--output` + `--fasta_out` switch | `-o/--genbank_out`, `--fasta_out` |
| | `-r/--references` | `-r/--reference` |
| | `--where` | `--anchor` |
| | `--search_hits` | `--annotation_type domain_search` |
| select_by_cds | `-o/--output` | `-o/--genbank_out`, `--fasta_out` |
| | `--pad_where` | `--anchor` |
| | `--skip_deduplicate` | `--no_deduplicate` |
| | `--search_hits` | `--annotation_type domain_search` |
| select_by_contig | `-o/--output` + `--fasta_out` switch | `-o/--genbank_out`, `--fasta_out` |
| | `--domain_type {domain,search,both}` | `--annotation_type {domainator,domain_search,all}` |
| | `--first` | `--max_records` |
| | `--length_lb`, `--length_ub` | `--min_length`, `--max_length` |
| | `--cds_count_lb`, `--cds_count_ub` | `--min_cds`, `--max_cds` |
| | `--search_score_lb`, `--search_score_ub` | `--min_search_score`, `--max_search_score` |
| select_by_coord | `-o/--output`, `-c/--coords` | `-o/--genbank_out`, `--fasta_out`, `--coords` |
| seq_dist | `--hdf`, `--dense_text` | `-o/--hdf_out`, `--dense_text_out` |
| | `--lb` | `--threshold` |
| | `--input_type`, `--reference_type` | `--input_format`, `--reference_format` |
| | `--params '"-c":4'` | `--diamond_args="-c 4"` |
| | | new: `--diamond_path` |
| sort_contigs | `-o/--output` | `-o/--genbank_out`, `--fasta_out` |
| split_by_cluster | `--outdir` | `--output_dir` |
| | `--lb` | `--threshold` |
| ssn_navigator | `-o/--output` | `-o/--json_out` |
| | `--mode` | `--query` |
| | `--id` | `--cluster` |
| structure_dist | `--hdf`, `--dense_text`, `--hits_tsv`, `--keep_db` | `-o/--hdf_out`, `--dense_text_out`, `--hits_tsv_out`, `--db_out` |
| | `--lb` | `--threshold` |
| structure_domainate | `-o/--output`, `--hits_tsv`, `--keep_db` | `-o/--genbank_out`, `--fasta_out`, `--hits_tsv_out`, `--db_out` |
| | `-r/--references` | `-r/--reference` |
| | `--max_domains` | `--max_hits_per_cds` |
| structure_search | `-o/--output`, `--hits_tsv`, `--keep_db` | `-o/--genbank_out`, `--fasta_out`, `--hits_tsv_out`, `--db_out` |
| | `-r/--references` | `-r/--reference` |
| structure_to_genbank | `-o/--output`, `--keep_db` | `-o/--genbank_out`, `--fasta_out`, `--db_out` |
| summary_report | `-o/--output`, `--html`, `--json` | `-o/--text_out`, `--html_out`, `--json_out` |
| | `--domains_table` | `--domains_tsv_out` |
| | `--co_occurence_fraction`, `--co_occurence_count` | `--co_occurrence_fraction_out`, `--co_occurrence_count_out` |
| transform_matrix | `--hdf`, `--dense_text` | `-o/--hdf_out`, `--dense_text_out` |
| | `--lb` | `--threshold` |
| | `--input_type` | `--value_type` |
| trim_contigs | `-o/--output` | `-o/--genbank_out`, `--fasta_out` |
| | `--kb_both`, `--cds_both` | `--kb_range`, `--cds_range` |
| | `--no_domain` | `--unannotated` |
| | | new: `-e/--evalue` and `--databases` limit which domain annotations count when trimming by domain |
| | `--unannotated` on its own did nothing | trims unannotated CDSs from both ends |
| format_db | `--log` | `--log_out` |
| | `--compress` / `--bgzf` | `--compress` |
| | `--cpu` (default 1) | `--cpu` (default 0: all physical cores) |

## Behavior changes

- `-e/--evalue` keeps annotations with an e-value *at or below* the cutoff in every program.
  select_by_contig, select_by_cds and extract_peptides used to drop an annotation whose
  e-value equaled the cutoff.
- domainate `--max_overlap` rejects a hit only when it overlaps a better hit by *more* than the
  allowed fraction, as filter_domains `--max_overlap` does. A hit overlapping by exactly that
  fraction used to be rejected.
- `--domain_expr` and `--taxonomy_expr` are checked when the program starts: a malformed
  expression is an error. `~` binds tighter than `&`, which binds tighter than `|`.
  Expressions like `A~B` or `A(B`, which used to be split into the wrong names, are errors.
- The `DOMAINATOR_GB_PARSER` environment variable no longer exists.
- Missing external programs (foldseek, cmbuild, cd-hit, usearch, diamond, datasets) give one
  error message that names the `--<program>_path` option.

## File formats

- Matrix hdf5 files are now version 2.0. Domainator 1.0 reads 0.x matrix files. 0.x reads 1.0
  files unless they store a symmetric matrix as its upper triangle, which 0.x rejects with an
  error.
- `DataMatrix.write()` no longer accepts the `"sparse"` and `"dense"` output types; use `"hdf"`.

## Python API

The library functions changed along with the command line, so that a function's parameters
have the same names as its command's options.

**Imports.** The package root exports the shared pieces:
`from domainator import parse_seqfiles, write_genbank, DataMatrix, DomainatorError`.
Each command's function is `from domainator.<command> import <command>`. Everything else is
internal; see `docs/stability.md`.

| 0.x | 1.0 |
|---|---|
| `domainator.utils` | split by topic; see the next table |
| `domainator.Taxonomy` | `domainator.taxonomy` (record-level helpers in `domainator.taxonomy.records`) |
| `domainator.domainator_db_download.domainator_db_download` | `domainator.db_download.db_download` |
| `domainator.domainator_format_db` | `domainator.format_db` |
| `domainator.plot_contigs.create_d3_html(output_file, contigs, ...)` | `domainator.plot_contigs.plot_contigs(records, html_out=..., ...)` |
| `DeduplicateGenbank().deduplicate_genbank(...)` and `.cluster_members` | `deduplicate_genbank(inputs, identity=..., clusters=[])`, which appends the clusters to the list |
| `MaxTree` and the edge functions imported from `domainator.data_matrix` | import them from `domainator.ssn_edges` (`SortedUndirectedEdges`, `CompactNeighborRankings` from `domainator.edge_types`) |
| `utils.pyhmmer_decode(x)` | `x` (pyhmmer 0.12 returns `str`) |

Where `domainator.utils` went:

| Module | Contents |
|---|---|
| `domainator.seqio` | `parse_seqfiles`, `write_genbank`, `get_file_type`, `detect_compression`, `open_seqfile`, `open_writable_seqfile`, `copy_SeqRecord`, ... |
| `domainator.hmmio` | `read_hmms`, `iter_hmms`, `open_hmm_file`, alphabets, `read_pyhmmer_fastas`, `read_infernal_cms`, ... |
| `domainator.model` | `DomainatorCDS`, `get_cds_unique_name`, `get_cds_name`, `region_id`, `copy_feature`, `sequence_hashes`, ... |
| `domainator.locations` | `slice_record`, `slice_record_from_location`, `pad_location`, `regions_overlap`, ... |
| `domainator.partial` | `get_fragment_status`, `filter_by_partial`, `PARTIAL_CHOICES`, ... |
| `domainator.taxonomy.records` | `get_taxid`, `filter_by_taxonomy`, `TaxonomyData`, ... |
| `domainator.expr` | `BooleanEvaluator` |
| `domainator.colors` | `get_palette`, `DISTINCT_COLORS`, ... |
| `domainator.offsets` | `get_offsets`, `partition_seqfile`, `count_peptides`, ... |
| `domainator.cli.common` | `list_and_file_to_dict_keys`, `partial_help` |

**Signatures.** Each command's function takes its main input as the only positional argument;
every other argument is keyword-only. Parameters were renamed to match the options:

| Function | 0.x | 1.0 |
|---|---|---|
| all record editors | `seq_iterator`, `sequence_iterator`, `contigs` (the input records) | `records` |
| all | `fasta_type` | `alphabet` |
| all | `threads`, `cpus` | `cpu` |
| all matrix tools | `lb` | `threshold` |
| extract_domains, extract_peptides, extract_unannotated, rotate_contigs, select_by_cds | `search_hits=True` | `annotation_type="domain_search"` |
| select_by_contig | `domain_type="domain"/"search"/"both"` | `annotation_type="domainator"/"domain_search"/"all"` |
| select_by_contig | `target_domains`, `target_contigs`, `domain_evalue`, `first`, `length_lb/ub`, `cds_count_lb/ub`, `search_score_lb/ub` | `domains`, `contigs`, `evalue`, `max_records`, `min_length`/`max_length`, `min_cds`/`max_cds`, `min_search_score`/`max_search_score` |
| select_by_cds, extract_peptides | `target_cdss`, `target_domains` | `cds_ids`, `domains` |
| filter_domains | `databases_keep`, `databases_filter`, `selected_domains` | `keep_databases`, `databases`, `domains` |
| extract_domains | `score`, `pad_up`, `pad_down` | `min_score`, `extend_up`, `extend_down` |
| extract_unannotated | `score`, `lb` | `min_score`, `min_length` |
| trim_contigs | `strip_domains`, `no_domain` | `domains`, `unannotated` |
| rotate_contigs | `where` | `anchor` |
| domainate | `max_hits`, `esm2_3Di_weights`, `esm2_3Di_device` | `max_hits_per_cds`, `esm2_3di_weights`, `esm2_3di_device` |
| domain_search | `decoy_names`, `dedup_both_strands` | `decoys`, `both_strands` |
| find_features | `gpu_device` | `device` |
| hmmer_search, hmmer_compare | `reference_files`, `score_cutoff` | `references`, `min_score` |
| hmmer_select | `query_regex`, `query_exact`, `query_contains` | `regex`, `exact`, `contains` |
| seq_dist, kmer_dist | `input_type`, `reference_path`, `reference_type` | `input_format`, `reference`, `reference_format` |
| seq_dist | `diamond_params` | `diamond_args` |
| structure_* | `reference_values`, `keep_db`, `hits_tsv` | `references`, `db_out`, `hits_tsv_out` |
| structure_domainate | `max_domains` | `max_hits_per_cds` |
| build_ssn, build_tree, build_projection | `metadata_files`, `xgmml`/`xgmml_file`, `newick_file`, `cluster_tsv`, `three_dee` | `metadata`, `xgmml_out`, `newick_out`, `cluster_tsv_out`, `three_d` |
| split_by_cluster | `outdir` | `output_dir` |
| color_genbank | `color_targets` | `targets` |
| color_table_to_legend | `svg` | `svg_out` |

**Matrices.** `seq_dist`, `kmer_dist`, `structure_dist` and `compare_contigs` return a
`DataMatrix` instead of writing files, so they no longer take `hdf`, `dense_text` or
`max_output_bytes`. Write the result with `matrix.write(path, "hdf")`.
`compare_contigs` takes records instead of file names, and no longer takes `reports`,
`contigs`, `name_by_order` or `databases` (give the metrics `databases`); its sorted GenBank
output is `compare_contigs.write_sorted_genbank(records, matrix, path, ...)`.

**Warnings and errors.** Importing Domainator no longer hides numpy's warnings for the whole
process. Missing optional structure-search packages (torch, esmologs) raise
`ExternalToolError` when they are needed, instead of a `NameError` later.
