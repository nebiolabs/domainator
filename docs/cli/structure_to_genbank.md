[CLI reference](README.md)
# domainator structure_to_genbank

Convert protein structures into protein GenBank records.

```
Convert protein structures into a Domainator protein genbank file.

Reads structure files (pdb/cif, optionally gzipped, given as files, globs, or a directory
searched recursively) or a prebuilt foldseek database, and writes one protein genbank
record per structure chain. This makes structures ordinary domainator inputs, so the
sequence-based tools (domainator domainate, domainator domain_search, domainator seq_dist, domainator build_ssn, the
reports) work on them with no further conversion.

Chain splitting, structure parsing and naming are all done by the structural aligner, so
records are named exactly as they are in a search against the same structures.

Use --db_out to also persist the aligner database that was built along the way. That
database carries C-alpha coordinates, so it can be reused as the -i or -r argument of
domainator structure_domainate / domainator structure_search without paying to rebuild it, and it is how
you turn a directory of reference structures into a reusable reference database.

Note that record coordinates are indices into the *observed* residue sequence extracted
from the structure, not original PDB residue numbers: a chain with a disordered loop has a
sequence shorter than its residue-number span.
```

Give at least one of these outputs: `--genbank_out`, `--fasta_out`, `--db_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | structure files (pdb/cif, optionally gzipped), a directory to search recursively, a glob, or the prefix of a prebuilt foldseek database. |
| `-o`, `--genbank_out` | FILE |  | Write the protein records to this GenBank file. Use - for standard output. |
| `--fasta_out` | FILE |  | Write the protein sequences to this FASTA file. Use - for standard output. |
| `--db_out` | FILE |  | also write the aligner database built from the input structures, using this path as its prefix. The database keeps C-alpha coordinates, so it can be reused as an input or reference for domainator structure_domainate and domainator structure_search. |
| `--store_3di` |  |  | store each chain's 3Di structural alphabet string in a /threedi qualifier on a full-length misc_feature. Note that this qualifier is a single opaque string, so it does not survive tools that slice records (domainator extract_domains, domainator select_by_coord). |
| `--algorithm` | {foldseek,reseek} | foldseek | structural alignment backend to use. |
| `--foldseek_path` | PROGRAM |  | path to the foldseek executable. If not supplied, foldseek is looked up on PATH. |
| `--reseek_path` | PROGRAM |  | path to the reseek executable. If not supplied, reseek is looked up on PATH. |
| `--foldseek_args` | STRING |  | extra options for the foldseek search command, for options domainator does not wrap, as one quoted string split like a shell command line. Use the = form: --foldseek_args="--option value". Options that domainator sets itself are refused. |
| `--reseek_args` | STRING |  | extra options for the reseek search command, for options domainator does not wrap, as one quoted string split like a shell command line. Use the = form: --reseek_args="--option value". Options that domainator sets itself are refused. |
| `--cpu` | INTEGER | 0 | number of threads to give the backend. 0 lets the backend use all available cores. |
| `--tmp_dir` | DIR |  | parent directory for temporary databases and alignment files. Structural searches can need a lot of scratch space, so point this at a large filesystem for big inputs. |
| `--reseek_stats` | {family,superfamily,fold} | superfamily | reseek only: SCOP level that reseek's p-value model is calibrated against. This shifts p-values by orders of magnitude, so it changes which hits pass --evalue. Ignored by other backends. |
| `--reseek_sensitivity` | {fast,sensitive} | sensitive | reseek only: search sensitivity. Ignored by other backends. |
| `--max_output_gb` | NUMBER | 25.0 | Fail if a projected output file would exceed this many GB. Set to 0 to disable the guardrail. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
