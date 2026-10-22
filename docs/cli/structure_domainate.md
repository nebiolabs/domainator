[CLI reference](README.md)
# domainator structure_domainate

Annotate protein structures with hits to reference structures.

```
Annotate protein structures with hits to reference structures.

The structure-input counterpart of domainator domainate. Reads structure files (pdb/cif,
optionally gzipped, given as files, globs, or a directory searched recursively) or a
prebuilt foldseek database, aligns reference structures against them, and writes a protein
genbank file where every input chain is one record carrying Domainator features for the
references that hit it.

Every input record is written, whether or not it had hits, unless --hits_only is given.
Use domainator structure_search instead when searching a large database with a few references.

References are aligned as the query and the input database as the target, so --evalue is
always computed against the database being annotated, exactly as it is for the rest of
domainator. The annotated coordinates are positions on the input record; the matching
positions on the reference are recorded in the /rstart, /rend and /rlen qualifiers.

Because the input database is built from real structures it carries C-alpha coordinates,
so --alignment_type 1 (TM-align) and the coordinate-based --metrics are available here.

Record coordinates are indices into the *observed* residue sequence extracted from the
structure, not original PDB residue numbers.
```

Give at least one of these outputs: `--genbank_out`, `--fasta_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | structures to annotate: pdb/cif files (optionally gzipped), a directory to search recursively, a glob, or the prefix of a prebuilt foldseek database. |
| `-r`, `--reference` | FILE ... | required | reference structures to annotate the input with. Same accepted forms as --input. |
| `-o`, `--genbank_out` | FILE |  | Write the records to this GenBank file. Use - for standard output. |
| `--fasta_out` | FILE |  | Write the records to this FASTA file. Use - for standard output. |
| `-e`, `--evalue` | NUMBER | 0.001 | threshold E value for a hit. E values are computed against the input database. |
| `--min_evalue` | NUMBER | 0.0 | minimum E value for a hit. Useful for excluding near-identical matches, such as a structure's hit to itself. |
| `--max_hits_per_cds` | INTEGER |  | maximum number of annotations to add per record. If not supplied, there is no limit. |
| `--max_overlap` | NUMBER | 1 | maximum fraction of overlap allowed between annotations on one record. 1 disables overlap filtering. |
| `--overlap_by_db` |  |  | filter overlapping annotations within each database independently, instead of across all databases. |
| `--hits_only` |  |  | write only records that had at least one hit. |
| `--no_annotations` |  |  | do not write Domainator features. Only meaningful together with --hits_only, where it selects hit records without annotating them. |
| `--database_name` | STRING |  | value for the /database qualifier on the annotations. Defaults to the reference database or directory name. |
| `--db_out` | FILE |  | also write the aligner database built from the input structures, using this path as its prefix, so it can be reused as --input later. |
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
