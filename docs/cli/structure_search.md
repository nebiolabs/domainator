[CLI reference](README.md)
# domainator structure_search

Select matching structure-derived protein records and mark their best-hit regions.

```
Select matching structure-derived protein records and mark their best-hit regions.

The structure-input counterpart of domainator domain_search, optimized for large
structure databases (for example a foldseek database built from AlphaFold DB) with
relatively few references. The -i/--input values identify the database being searched;
-r/--reference supplies the reference structures to search for.

Only matching records are written, using each target's full observed amino-acid
sequence, not a crop to the alignment. Search best-hit annotations (Domain_Search
features) mark the regions homologous to the best-matching references.
--add_annotations also writes domain annotations (Domainator features) for retained
hits; the best-hit annotations are already present without this option.

References are aligned as the query and the searched database as the target, which is both
what foldseek is optimized for and what makes --evalue mean the same thing it means
elsewhere in domainator: it is computed against the database being searched. The hit
coordinates are positions on the database record; the matching positions on the reference
are recorded in the /rstart, /rend and /rlen qualifiers.

Hit records are built from the sequences the backend reports, so the searched database is
never read end to end and its size does not bound memory. Note that --max_seqs bounds how
many hits the backend keeps per reference; the default is low for a large database, and a
warning is issued whenever a reference saturates it.

Use domainator structure_domainate instead to annotate every input structure rather than select
hits, and domainator structure_to_genbank to convert structures without searching at all.
```

Give at least one of these outputs: `--genbank_out`, `--fasta_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | the structure database to search: the prefix of a prebuilt foldseek database, or pdb/cif files (optionally gzipped), a directory to search recursively, or a glob. |
| `-r`, `--reference` | FILE ... | required | Reference structures to search for in the --input database. Same accepted forms as --input. |
| `-o`, `--genbank_out` | FILE |  | Write the records to this GenBank file. Use - for standard output. |
| `--fasta_out` | FILE |  | Write the records to this FASTA file. Use - for standard output. |
| `-e`, `--evalue` | NUMBER | 0.001 | threshold E value for a hit. E values are computed against the searched database. |
| `--min_evalue` | NUMBER | 0.0 | minimum E value for a hit. Useful for excluding a reference's match to its own copy in the database. |
| `--max_hits` | INTEGER |  | return at most this many records, the best-scoring ones. Applies to the whole search, not per reference. When set, output is sorted best-first and written after the search completes. |
| `--max_overlap` | NUMBER | 1 | maximum fraction of overlap allowed between annotations on one record. 1 disables overlap filtering. |
| `--overlap_by_db` |  |  | filter overlapping annotations within each database independently, instead of across all databases. |
| `--add_annotations` |  |  | Add domain annotations (Domainator features) for retained hits, in addition to the search best-hit annotations (Domain_Search features) already present. |
| `--database_name` | STRING |  | value for the /database qualifier on the annotations. Defaults to the reference database or directory name. |
| `--db_out` | FILE |  | when the input is structure files rather than a prebuilt database, also write the database built from them, using this path as its prefix. |
| `--hits_tsv_out` | FILE |  | write the aligner's own hit table here, for debugging or for cross-checking against a hand-run search. With reseek this is the table before the per-reference --max_seqs cap is applied, so it may hold more rows than the output. |
| `--algorithm` | {foldseek,reseek} | foldseek | structural alignment backend to use. |
| `--foldseek_path` | PROGRAM |  | path to the foldseek executable. If not supplied, foldseek is looked up on PATH. |
| `--reseek_path` | PROGRAM |  | path to the reseek executable. If not supplied, reseek is looked up on PATH. |
| `--foldseek_args` | STRING |  | extra options for the foldseek search command, for options domainator does not wrap, as one quoted string split like a shell command line. Use the = form: --foldseek_args="--option value". Options that domainator sets itself are refused. |
| `--reseek_args` | STRING |  | extra options for the reseek search command, for options domainator does not wrap, as one quoted string split like a shell command line. Use the = form: --reseek_args="--option value". Options that domainator sets itself are refused. |
| `--cpu` | INTEGER | 0 | number of threads to give the backend. 0 lets the backend use all available cores. |
| `--tmp_dir` | DIR |  | parent directory for temporary databases and alignment files. Structural searches can need a lot of scratch space, so point this at a large filesystem for big inputs. |
| `--reseek_stats` | {family,superfamily,fold} | superfamily | reseek only: SCOP level that reseek's p-value model is calibrated against. This shifts p-values by orders of magnitude, so it changes which hits pass --evalue. Ignored by other backends. |
| `--reseek_sensitivity` | {fast,sensitive} | sensitive | reseek only: search sensitivity. Ignored by other backends. |
| `--alignment_type` | {0,1,2} |  | foldseek only: how to compute the alignment: 0 = 3Di only, 1 = TM-align, 2 = 3Di+AA (the default when unset). Type 1 requires C-alpha coordinates in both databases. reseek has no alignment-type option and rejects this. |
| `--metrics` | {tmscore,lddt,rmsd,prob} ... |  | additional per-hit structural metrics to compute and store as feature qualifiers. tmscore, lddt and rmsd require C-alpha coordinates in both databases. |
| `--max_seqs` | INTEGER |  | maximum number of hits to keep per reference, keeping the most significant ones. Defaults to the backend's own default: 1000 for foldseek, no cap for reseek. For foldseek this bounds the search itself, so raising it costs time but improves coverage; for reseek it bounds only the output, since reseek has no native cap and does the full search either way. A warning is printed whenever a reference saturates the cap. |
| `--device` | STRING | cpu | where to run the search. 'cpu' for a CPU-only search, or a CUDA device string ('cuda', 'cuda:0', 'cuda:1', ...) for a GPU-accelerated search. GPU search requires a recent foldseek. |
| `--max_output_gb` | NUMBER | 25.0 | Fail if a projected output file would exceed this many GB. Set to 0 to disable the guardrail. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
