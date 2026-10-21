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
| `-e/--evalue` in filtering programs | unchanged, but the default is now "no cutoff" | was 1e8 in some programs and none in others |
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
| format_db | `--log` | `--log_out` |
| | `--compress` / `--bgzf` | `--compress` |
| | `--cpu` (default 1) | `--cpu` (default 0: all physical cores) |

## File formats

- Matrix hdf5 files are now version 2.0. Domainator 1.0 reads 0.x matrix files. 0.x reads 1.0
  files unless they store a symmetric matrix as its upper triangle, which 0.x rejects with an
  error.
- `DataMatrix.write()` no longer accepts the `"sparse"` and `"dense"` output types; use `"hdf"`.
