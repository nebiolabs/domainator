[CLI reference](README.md)
# domainator color_genbank

Add color qualifiers to CDS and domain features.

```
Add colors to domains or CDS annotations in genbank files

Takes a genbank file which has been annotated using the domainator program, adds color annotations to the sequence features.

Domains are colored based on domain name.

CDSs are colored based on the associated domains, prioritized by domain score.

If a color_table is provided, then the colors in the table will be used, otherwise a random color will be assigned to each domain.
```

Give at least one of these outputs: `--genbank_out`, `--fasta_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | Genbank filenames. Use - to read standard input. |
| `-o`, `--genbank_out` | FILE |  | Write the records to this GenBank file. Use - for standard output. |
| `--fasta_out` | FILE |  | Write the records to this FASTA file. Use - for standard output. |
| `--contigs` | STRING ... |  | only color contigs with ids in this list (other contigs will be excluded from output). Additive with --contigs_file. default: all contigs. |
| `--contigs_file` | FILE |  | only color contigs with ids listed in this file (one per line). Additive with --contigs. default: all contigs. |
| `--color_table` | FILE |  | tab separated file with two columns and no header, columns are: domain id, hex color. For example: CCDB   cc0000 |
| `--color_table_out` | FILE |  | tab separated file with two columns and no header, columns are: annotation, hex color. Written after the color table is updated with new colors, for example if using --color_by, but not supplying an external color table. |
| `--search_hit_color` | STRING |  | color to use for search hit annotations. |
| `--clear` |  |  | if set, then remove all color annotations from the genbank file before adding any new ones. |
| `--databases` | STRING ... |  | Domain databases to use for coloring. default: all databases. |
| `--target` | {domains,cdss,domains_and_cdss} | cdss | what to color. domains: color only domain annotations. cdss: color only CDS annotations. domains_and_cdss: color both domains and CDS annotations. |
| `--domain_database` |  |  | if set, then color domains based on their database instead of their name. If used with --color_table, the first column must contain database ids. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
