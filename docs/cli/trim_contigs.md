[CLI reference](README.md)
# domainator trim_contigs

Trim the ends of contigs.

```
Trim contigs based on CDS or domain names or presence of domains in CDSs

Circular contigs are linearized before trimming.

Typically used for nucleotide sequences, but can be used for protein sequences as well, if you use the 'kb' options to mean 1000 amino acids. CDS and domain options are not supported for protein sequences.
```

Give at least one of these outputs: `--genbank_out`, `--fasta_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | GenBank files. Use - to read standard input. |
| `-o`, `--genbank_out` | FILE |  | Write the records to this GenBank file. Use - for standard output. |
| `--fasta_out` | FILE |  | Write the records to this FASTA file. Use - for standard output. |
| `--kb_range` | NUMBER |  | How many kb to trim from both sides of each contig. |
| `--kb_left` | NUMBER |  | How many kb to trim from the left side of each contig. |
| `--kb_right` | NUMBER |  | How many kb to trim from the right side of each contigs. |
| `--cds_range` | INTEGER |  | How many CDSs to trim from both sides of each contig. |
| `--cds_left` | INTEGER |  | How many CDSs to trim from the left side of each contig. |
| `--cds_right` | INTEGER |  | How many CDSs to trim from the right side of each contig. |
| `--contigs` | STRING ... |  | Only return contigs with these names. |
| `--contigs_file` | FILE |  | text file containing names of contigs to return. |
| `--domains` | STRING ... |  | Trim CDSs containing domains with these names from the ends of each contig. |
| `--unannotated` |  |  | Trim CDSs containing no domain annotations from the ends of each contig. |
| `--domains_file` | FILE |  | text file containing names of domains to trim contigs containing. |
| `--domain_expr` | STRING |  | a boolean expression using operators & (AND), \| (OR), and ~(NOT), to specify a desired combination of domains for CDSs to trim |
| `-e`, `--evalue` | NUMBER |  | the evalue cutoff (max) for domain annotations to be considered when trimming by domain. Default: no cutoff. |
| `--databases` | STRING ... |  | Only consider domain annotations from these databases when trimming by domain. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
