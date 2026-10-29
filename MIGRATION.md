# Migrating from Domainator 0.x to 1.0

Domainator 1.0 makes the command line consistent across all programs. Old option names are
**not** kept as aliases: a script that uses one fails with an "unrecognized arguments" error
that names it, so nothing changes silently. This guide lists every change.

## Installation

- Python 3.12 or newer is required. Recreate the conda environment (`conda env create --yes -f
  conda_env.yml`), which now also installs Infernal and no longer installs HMMER or TMbed.
  TMbed (for `find_features --tmbed`) is installed separately; see "Optional: TMbed" in the
  README.
- UMAP is an optional extra: `pip install "domainator[projection]"` for
  `build_projection --algorithm umap`. matplotlib and seaborn are no longer needed. `torch`
  (for the ESM-2 and TMbed features) is the `structure` extra. The `core` extra is gone (its
  packages are ordinary dependencies), and so is the `bashplotlib` dependency.
- `domainator doctor` lists the optional programs and packages and which of them are installed.

## Running programs

| 0.x | 1.0 |
|---|---|
| `domainate.py -i in.gb ...` | `domainator domainate -i in.gb ...` |
| `domainator_db_download.py ...` | `domainator db_download ...` |
| `domainator_format_db.py ...` | `domainator format_db ...` |
| `partition_seqfile.py`, `partition_seqids.py` | removed; index databases with `domainator format_db --index` instead |
| `domainator_server` | removed; the browser GUI is becoming a separate project, built on `domainator describe` (see docs/tool_manifest.md) |
| `build_ssn_viewer.py` | removed; `domainator build_ssn --html_out` / `--dvb_out` write the viewer and its bundle, and `domainator viewer` writes an empty viewer page (see below) |

`domainator --help` lists all commands. `domainator completion bash` prints a tab-completion
script (see the README).

Options must be spelled out in full. 0.x accepted any unambiguous prefix (`--max_hits_per` for
`--max_hits_per_query`); 1.0 does not.

Every program takes `-v/--verbose` (debug output) and `-q/--quiet` (errors only: no progress
messages or warnings). `build_tree --quiet` is now this shared option.

Problems with inputs, options, or external programs print one line,
`domainator <command>: error: ...`, and exit with status 1; usage errors exit with status 2.
`-v` shows the traceback. In Python, these errors are subclasses of
`domainator.errors.DomainatorError` (and of `ValueError` or `RuntimeError`, as before). Tools that
used to `exit(1)` from library code (seq_dist's diamond errors) or raise `SystemExit` for
`--max_output_gb` now raise `ExternalToolError` or `OutputSizeLimitExceeded`.

## Inputs and outputs

### Annotation versus selection/search

`domainate` and `structure_domainate` always annotate and retain every input record,
including records without hits. Both remove `--hits_only` and `--no_annotations`,
and the corresponding public Python keywords. `domainate` also removes
`--include_taxids`, `--exclude_taxids`, `--taxonomy_expr`, `--ncbi_taxonomy_path`,
and `--taxonomy_update`, and its public taxonomy-related Python keywords.

Use `domain_search` or `structure_search` for hit-only output. Search output includes
`Domain_Search` best-hit features by default; `--add_annotations` also adds `Domainator`
features. Search defaults therefore are not byte-for-byte replacements for old
annotation-tool hit-only output. For whole-contig sequence hits, also consider
`domain_search --whole_contig`.

To annotate only records from a clade, filter first:

```bash
domainator select_by_contig -i genomes.gb --include_taxids 2 -o bacteria.gb
domainator domainate -i bacteria.gb -r pfam.hmm -o annotated.gb
```

### Taxonomy filtering and unresolved labels

`extract_domains`, `extract_peptides`, `extract_unannotated`, `select_by_cds`,
`select_by_contig`, and `domain_search` share include/exclude taxids, Boolean
`--taxonomy_expr`, database setup options, and
`--taxonomy_unknown {exclude,include,evaluate,error}`.

The default **exclude** drops unresolved taxonomy, even for exclusion-only filters
and negation. Missing labels, uncovered regions, unidentified taxid 32644, deleted
taxids, and IDs absent from the local database are unresolved. Merged IDs use their
valid replacement's taxonomy. This intentionally changes older filtering that could
retain unidentified/deleted taxa under negation.

For `--taxonomy_expr "~2"` (not Bacteria):

| Candidate | exclude (default) | evaluate | include | error |
|---|---|---|---|---|
| Known bacterial protein | drop | drop | drop | drop |
| Known archaeal protein | keep | keep | keep | keep |
| Unlabeled protein, possibly bacterial | drop | keep | keep | error |

`evaluate` evaluates unresolved taxonomy against an empty lineage: positive inclusion
fails, while exclusion-only criteria/negation can pass. `include` always passes the
taxonomy gate for unresolved candidates, even with positive inclusion. Neither option
proves those candidates lie outside an excluded clade. `error` reports the unresolved
candidate. Other selection criteria still apply with every policy. With no taxonomy
criteria, the policy is inactive and no database is loaded.

CDS/domain extraction resolves taxonomy on the original anchor before translation,
extension, combination, or splicing. The longest source fully covering the anchor
wins, not the narrowest source; a gappy source must cover every anchor part.
Neighborhoods and extended domains may contain sequence from other taxa.
`select_by_contig` and protein-only `extract_unannotated` filter whole records.
Taxonomy restrictions are independent of `--invert`.

`db_download` retains its existing compiled-membership filtering behavior; this change
does not add an unknown-taxonomy option to downloads.

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
matrix and bundle readers (build_projection, build_ssn, build_tree, viewer,
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
| Cytoscape XGMML / Newick tree / SVG | `--xgmml_out` / `--newick_out` / `--svg_out` |
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
| `--search_hits` | `--has_search_hit` | additive search-hit selection; extract_unannotated uses `--annotation_type all` instead |
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

### Annotation sources versus search-hit selection

`--annotation_type {domainator,domain_search,all}` consistently chooses the annotation
sources used by normal domain criteria. Names from both sources are combined for `all`.
It does not independently select search-hit CDSs or records. Use `--has_search_hit` for
the additive behavior of the old `--search_hits`: normal selection OR search-hit presence.
The combined selection is inverted by `--invert` where available. Search hits must pass
quality cutoffs, but this independent branch ignores domain names, expressions, and
database restrictions. To filter search hits by database, use `--annotation_type
domain_search` or `all` with normal criteria, without `--has_search_hit`.

This also replaces the additive use of `--annotation_type domain_search` in earlier
1.0 development versions. With no other selectors, `--has_search_hit` selects only search
hits; the source option alone does not restrict CDS/contig selection or peptide extraction.
For contig-level nucleotide hits in `select_by_cds`, `--has_search_hit` enables them
automatically; ordinary source-based selection uses `--include_nucleic_acids`.

`extract_unannotated` instead uses `--annotation_type all` to count both annotation types
as covered regions; it has no `--has_search_hit`. `filter_domains` also has no independent
search-hit selector: it filters the selected annotation types and preserves the others.
See [annotation selection](docs/file_formats.md#selecting-annotation-sources-and-search-hits).

### CDS and peptide selection

`extract_peptides` now accepts `--domain_expr`, additive with CDS names, domain
names, unannotated, and search-hit presence, just like `select_by_cds`. For example:

```bash
domainator extract_peptides -i annotated.gb --domain_expr "A & B & ~C" \
  --fasta_out proteins.faa
```

Both tools now apply `--invert` only to the combined OR selectors. Strand is an
independent restriction, like contig IDs, fragment status, and taxonomy. Previously,
`--domains A --strand f --invert` could reinclude reverse-strand CDSs; it now selects
only forward-strand CDSs without A. This is an intentional pre-v1 behavior correction.

With no positive selectors, both CLI and library functions select all otherwise
eligible CDSs. Inversion without selectors warns once and is ignored; strand,
fragment-status, and taxonomy restrictions still apply. In particular,
`extract_peptides(records)` now extracts eligible peptides without requiring the
library-only `extract_all=True`. `select_by_cds(records, invert=True)` no longer
silently selects nothing. Explicit empty selector sets remain different from `None`:
an empty set selects nothing, or all eligible CDSs when inverted.

`select_by_cds` now accepts `--name_field`, like `extract_peptides`. This field affects
CDS name matching; on peptide extraction it also determines peptide IDs. Default
name precedence has not changed: peptides prefer `protein_id`, while neighborhood
selection prefers `gene_id`, then `locus_tag`. Specify the same name field in both
tools when necessary. Neighborhood output IDs remain coordinate-based.
See [the full selection contract](docs/file_formats.md#cds-selection-algebra-and-names).

## Program-specific changes

Options not listed are unchanged.

| Program | 0.x | 1.0 |
|---|---|---|
| build_projection | `--png_out`, `--xgmml`, `--coords_out` | `--html_out` (an interactive plot that saves itself as PNG or SVG; replaces `--png_out`), `-o/--xgmml_out`, `--coords_out` |
| | `--3d` | `--three_d` |
| | `--params` | `--umap_args`, `--tsne_args`, `--pca_args`: `key=value` pairs, e.g. `--umap_args="n_neighbors=30 min_dist=0.1"` |
| build_ssn | `--xgmml`, `--cluster_tsv` | `-o/--xgmml_out`, `--cluster_tsv_out` |
| | `--lb` | `--threshold` |
| build_ssn_viewer | `-o/--output` (a `.dsnv` bundle) | `build_ssn --dvb_out` (a `.dvb` bundle) |
| | `--html --embed_data` | `build_ssn --html_out` |
| | `--html` alone (an empty page) | `viewer --html_out` |
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
| | `--deduplicate`, `--no_deduplicate` | `--dedup id`, `--dedup none` (still off by default) |
| | `--deduplicate_by_sequence` | `--dedup sequence` |
| | `--dedup_both_strands`, `--both_strands` (with sequence deduplication) | `--dedup sequence_rc` |
| domainate | `-o/--output` | `-o/--genbank_out`, `--fasta_out` |
| | `-r/--references` | `-r/--reference` |
| | `--max_domains` | `--max_hits_per_cds` |
| | `--esm2_3Di_weights`, `--esm2_3Di_device` | `--esm2_3di_weights`, `--esm2_3di_device` |
| | `--offset`, `--recs_to_read` | removed |
| | | new: `--foldseek_path` |
| enum_report | `-o/--output`, `--html`, `--json` | `-o/--tsv_out`, `--html_out`, `--json_out` |
| | one flag per column: `--domains --length --taxname genus species --qualifier CDS product` | `--columns domains length taxname:genus:species qualifier:CDS:product` (arguments follow the column name after `:`; `domainator enum_report --help` lists the columns) |
| extract_domains | `-o/--output` + `--fasta_out` switch | `-o/--genbank_out`, `--fasta_out` |
| | `--score` | `--min_score` |
| | `--search_hits` | `--has_search_hit` |
| | `--pad_up`, `--pad_down` | `--extend_up`, `--extend_down`, and new `--extend_range` |
| extract_peptides | `-o/--output` + `--fasta_out` switch | `-o/--genbank_out`, `--fasta_out` |
| | `--cdss`, `--cdss_file` | `--cds`, `--cds_file` |
| | `--search_hits` | `--has_search_hit` |
| extract_unannotated | `-o/--output` + `--fasta_out` switch | `-o/--genbank_out`, `--fasta_out` |
| | `--score` | `--min_score` |
| | `--lb` | `--min_length` |
| | `--search_hits` | `--annotation_type all` |
| filter_domains | `-o/--output` | `-o/--genbank_out`, `--fasta_out` |
| | `-m/--max_overlap` | `--max_overlap` |
| | `--databases_filter` | `--databases` |
| | `--databases_keep` | `--keep_databases` |
| find_features | `-o/--output` | `-o/--genbank_out`, `--fasta_out` |
| | `--gpu_device` | `--device` |
| genbank_to_fasta | `-o/--output` | `-o/--fasta_out` |
| hmmer_build, hmmer_select, hmmer_search | `-o/--output` | `-o/--hmm_out` |
| hmmer_compare | `-o/--output` | `-o/--tsv_out` |
| | `--alignments` (alignments mixed into the TSV) | `--alignments_out PATH`, a separate file; the TSV is now always a plain table |
| | `--score_cutoff` | `--min_score` |
| hmmer_report | `-o/--output`, `--html`, `--json` | `-o/--tsv_out`, `--html_out`, `--json_out` |
| | one flag per column: `--acc --length --append name str value` | `--columns acc length append:name:str:value` |
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
| | `--search_hits` | `--has_search_hit` |
| select_by_cds | `-o/--output` | `-o/--genbank_out`, `--fasta_out` |
| | `--pad_where` | `--anchor` |
| | `--deduplicate`, `--skip_deduplicate`, `--no_deduplicate` | `--dedup id`, `--dedup none`, `--dedup none` (ID deduplication remains the default) |
| | `--search_hits` | `--has_search_hit` |
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
| ssn_navigator | `-o/--output`; config key `cluster_id` | `-o/--json_out`; config key `cluster` (the option is still `--cluster`) |
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

### Output deduplication policies

domain_search and select_by_cds use one single-choice option,
`--dedup {none,id,sequence,sequence_rc}`, with these policies:

| Policy | Meaning |
|---|---|
| `none` | Keep all output records. |
| `id` | Keep one representative per output record ID. |
| `sequence` | Keep unique output sequences, ignoring case and terminal `*` characters. |
| `sequence_rc` | Like `sequence`, also collapsing nucleotide reverse complements. Not valid for protein output. |

Choose exactly one policy. The old switches and the pre-release multi-value
`--dedup` syntax are removed, without aliases.

Domainator assumes that **the same record ID identifies the same sequence**.
Different IDs can identify identical sequences; annotations and focus-hit scores
can differ between copies of the same sequence. Consequently, sequence uniqueness
also enforces ID uniqueness, so composable selectors are unnecessary.
Reusing an ID for different
sequences is unsupported; give those records distinct IDs before using them
together. Deduplication relies on this invariant rather than globally validating it.

Omitting the option preserves each command's default: no deduplication for
domain_search, and ID deduplication for select_by_cds. JSON/YAML configurations use
a string, for example `"dedup": "sequence_rc"`, not an array.
`--dedup none` overrides deduplication from a configuration file.

Deduplication compares extracted/translated records before padding. Without
`--max_hits`, the first arriving representative is kept; parallel search arrival
order is arbitrary. With `domain_search --max_hits N`, each worker and the global
merge retain at most N winning records, replacing a representative when a
better-scoring copy arrives. Only the final winners are sorted, so duplicates do
not consume result slots and no all-hit buffer is needed. Each policy's
retained records and heap metadata are O(N), in
addition to the normal bounded task/result queues, input batches, and search state.
Equal scores favor earlier arrivals; representatives can depend on parallel
batch completion order.

`deduplicate_genbank --both_strands` is unchanged as an option: it belongs to the
separate clustering interface. Exact hash deduplication now correctly equates only
a nucleotide sequence and its reverse complement, not its plain complement or
plain reversal.

## Behavior changes

- `-e/--evalue` keeps annotations with an e-value *at or below* the cutoff in every program.
  select_by_contig, select_by_cds and extract_peptides used to drop an annotation whose
  e-value equaled the cutoff, and so did the searches in domainate, domain_search and the
  structure tools.
- domainate and domain_search refuse `-Z 0` together with `--gene_call`: the number of
  targets isn't known until the genes are called, so give the database size. `-Z 0` counts the
  targets before any taxonomy filter (as before; now documented), so filtering doesn't change
  E-values.
- domainate `--max_overlap` rejects a hit only when it overlaps a better hit by *more* than the
  allowed fraction, as filter_domains `--max_overlap` does. A hit overlapping by exactly that
  fraction used to be rejected.
- `--domain_expr` and `--taxonomy_expr` are checked when the program starts: a malformed
  expression is an error. `~` binds tighter than `&`, which binds tighter than `|`.
  Expressions like `A~B` or `A(B`, which used to be split into the wrong names, are errors.
- The `DOMAINATOR_GB_PARSER` environment variable no longer exists.
- Tabs inside FASTA sequence lines are dropped, like spaces (Biopython's FASTA reader does this).
- domain_search `--max_hits` breaks equal-score ties by result arrival order (which
  depends on parallel batch completion); sort_contigs uses input order.
- Missing external programs (foldseek, cmbuild, cd-hit, usearch, diamond, datasets) give one
  error message that names the `--<program>_path` option.

## File formats

- Matrix hdf5 files are now version 2.0. Domainator 1.0 reads 0.x matrix files. 0.x reads 1.0
  files unless they store a symmetric matrix as its upper triangle, which 0.x rejects with an
  error.
- `DataMatrix.write()` no longer accepts the `"sparse"` and `"dense"` output types; use `"hdf"`.
- The `--json_out` documents of summary_report, matrix_report and ssn_navigator start with a
  `schema_version` key; ssn_navigator's also has `query`, the mode that answered.
- Networks for the viewer are `.dvb` bundles (`build_ssn --dvb_out`), the format trees and
  projections use too, instead of `.dsnv` bundles. The viewer, `domainator viewer -i` and
  ssn_navigator still read `.dsnv` files, and the viewer saves their sessions as `.dvb`.
- build_tree writes Newick branch lengths to six significant digits instead of two decimals.
- XGMML attribute values are no longer escaped twice: a node named `A&B` used to show up in
  Cytoscape as `A&amp;B`.
- The GenBank qualifiers are documented in `docs/file_formats.md`, which is the 1.x contract.

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
| `build_ssn_viewer.build_ssn_viewer_bundle(matrix, metadata_files=..., categorical_columns=...)` | `build_ssn.build_ssn(matrix, dvb_out=..., metadata=..., categorical=...)`, which writes the bundle (`ssn_bundle.load_bundle` reads it back) |
| (none) | `genbank_to_fasta.genbank_to_fasta(input, fasta_out=...)` |
| `summary_report(records, out_text_handle=..., out_html_handle=..., out_json_handle=..., focus_domains=..., report_taxonomy=..., domains_table=..., co_occurrence_fraction_outfile=..., co_occurrence_count_outfile=...)` | `summary_report(records, text_out_handle=..., html_out_handle=..., json_out_handle=..., domains=..., taxonomy=..., domains_tsv_out=..., co_occurrence_fraction_out=..., co_occurrence_count_out=...)` |
| `matrix_report(matrix, out_text_handle=..., out_html_handle=..., out_json_handle=..., profile_stages=...)` | `matrix_report(matrix, text_out_handle=..., html_out_handle=..., json_out_handle=...)` (pass `stage_timings=[]` to collect timings) |
| `ssn_navigator(bundle, mode=..., cluster_id=...)` | `ssn_navigator(bundle, query=..., cluster=...)` |
| `format_db(inputs, log=...)`, `deduplicate_genbank(inputs, log_handle=...)` | `log_out_handle=...` |
| `db_download(..., api_key=...)` | `db_download(..., ncbi_api_key=...)` |
| `hmmer_search(input_files, hmmer_handle=...)` | `hmmer_search(input_files, hmm_out_handle=...)` |
| Keyword defaults that differed from the command's (`domainate(evalue=10)`, `cpu=1`, required keywords) | Every keyword named like an option defaults like the option: notably `domainate` evalue 0.001, not 10; `cpu` 0 (all physical cores) everywhere; e-value, length and hit-count limits default to None (no limit), as on the command line |
| `MaxTree` and the edge functions imported from `domainator.data_matrix` | import them from `domainator.ssn_edges` (`SortedUndirectedEdges`, `CompactNeighborRankings` from `domainator.edge_types`) |
| `utils.pyhmmer_decode(x)` | `x` (pyhmmer 0.12 returns `str`) |
| `write_genbank(records, path, preserve_original=..., add_domain_cds_annotations=..., replace_space_with=...)` | `write_genbank(records, path)`; it never changes the records, so `preserve_original` is not needed |

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
| extract_domains, extract_peptides, rotate_contigs, select_by_cds | `search_hits=True` | `has_search_hit=True` |
| extract_unannotated | `search_hits=True` | `annotation_type="all"` |
| select_by_contig | `domain_type="domain"/"search"/"both"` | `annotation_type="domainator"/"domain_search"/"all"` |
| select_by_contig | `target_domains`, `target_contigs`, `domain_evalue`, `first`, `length_lb/ub`, `cds_count_lb/ub`, `search_score_lb/ub` | `domains`, `contigs`, `evalue`, `max_records`, `min_length`/`max_length`, `min_cds`/`max_cds`, `min_search_score`/`max_search_score` |
| select_by_cds, extract_peptides | `target_cdss`, `target_domains` | `cds_ids`, `domains` |
| filter_domains | `databases_keep`, `databases_filter`, `selected_domains` | `keep_databases`, `databases`, `domains` |
| extract_domains | `score`, `pad_up`, `pad_down` | `min_score`, `extend_up`, `extend_down` |
| extract_unannotated | `score`, `lb` | `min_score`, `min_length` |
| trim_contigs | `strip_domains`, `no_domain`, `kb_truncation`, `cds_truncation` | `domains`, `unannotated`, `kb_range`, `cds_range` |
| db_download | `output_file_path`, `num_recs`, `success_rec_log` | `output`, `max_records`, `success_log_out` |
| split_by_cluster | `matrix_path` | `matrix` |
| infernal_build | `informat` | `input_format` |
| rotate_contigs | `where` | `anchor` |
| domainate | `max_hits`, `esm2_3Di_weights`, `esm2_3Di_device` | `max_hits_per_cds`, `esm2_3di_weights`, `esm2_3di_device` |
| domain_search | `partitions` (the input), `decoy_names` | `targets` (SeqRecords, or the same partitions), `decoys` |
| domain_search | `deduplicate_by_sequence`, `both_strands` / `dedup_both_strands` | `dedup="sequence"`, `dedup="sequence_rc"` |
| select_by_cds.deduplicate_and_pad | `deduplicate=True` / `False` | `dedup="id"` / `"none"` |
| find_features | `gpu_device` | `device` |
| hmmer_search, hmmer_compare | `reference_files`, `score_cutoff` | `references`, `min_score` |
| hmmer_compare | `out_handle`, `alignments` | `tsv_out_handle`, `alignments_out_handle` (either can be None) |
| hmmer_select | `query_regex`, `query_exact`, `query_contains` | `regex`, `exact`, `contains` |
| seq_dist, kmer_dist | `input_type`, `reference_path`, `reference_type` | `input_format`, `reference`, `reference_format` |
| seq_dist | `diamond_params` | `diamond_args` |
| structure_* | `reference_values`, `keep_db`, `hits_tsv` | `references`, `db_out`, `hits_tsv_out` |
| structure_domainate | `max_domains` | `max_hits_per_cds` |
| build_ssn, build_tree, build_projection | `metadata_files`, `xgmml`/`xgmml_file`, `newick_file`, `cluster_tsv`, `three_dee` | `metadata`, `xgmml_out`, `newick_out`, `cluster_tsv_out`, `three_d` |
| split_by_cluster | `outdir` | `output_dir` |
| color_genbank | `color_targets` | `targets` |
| color_table_to_legend | `svg` | `svg_out` |

The pre-release `dedup` collection parameter is now one string in Python as well
as in configuration files. For example, use `dedup="sequence_rc"` instead of
`dedup=["id", "sequence_rc"]`. `domain_search` defaults to `dedup="none"`.

**Matrices.** `seq_dist`, `kmer_dist`, `structure_dist` and `compare_contigs` return a
`DataMatrix` instead of writing files, so they no longer take `hdf`, `dense_text` or
`max_output_bytes`. Write the result with `matrix.write(path, "hdf")`.
`compare_contigs` takes records instead of file names, and no longer takes `reports`,
`contigs`, `name_by_order` or `databases` (give the metrics `databases`); its sorted GenBank
output is `compare_contigs.write_sorted_genbank(records, matrix, path, ...)`.

**Warnings and errors.** Importing Domainator no longer hides numpy's warnings for the whole
process. Missing optional structure-search packages (torch, esmologs) raise
`ExternalToolError` when they are needed, instead of a `NameError` later.

**Biopython.** Domainator now uses upstream Biopython (1.85 to 1.88) instead of its own copy,
`domainator.Bio`, which is gone. Import `SeqRecord`, `Seq`, `SeqFeature` and so on from `Bio`.
Domainator's records are ordinary Biopython records, so they work with other Biopython code.
What the copy added is now in functions:

| 0.x | 1.0 |
|---|---|
| `loc.stranded_start`, `loc.stranded_end`, `loc.stranded_start_human_readable`, `loc.stranded_end_human_readable` | `stranded_start(loc)`, ... in `domainator.location_ops` |
| `loc.overlay(offset, length)`, `loc.circular_envelope(n, strand)`, `loc.inverted(n)`, `loc.overlaps(other)`, `loc.overlap_size(other)`, `loc.contains(other)`, `loc.touches(other)`, `loc.copy()`, `CompoundLocation.merge_parts(parts)` | `overlay(loc, offset, length)`, ... in `domainator.location_ops` (`copy_location(loc)`) |
| `record.set_hit_best_score(x)`, `record.get_hit_scores()`, `record.set_dist_to_start(x)`, ... | `set_hit_best_score(record, x)`, `get_hit_scores(record)`, ... in `domainator.model` |
| comparing records with `<` (by best hit score) | `sorted(records, key=get_hit_best_score)` |
| `feature.strand` | `feature.location.strand` |
