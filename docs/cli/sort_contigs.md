[CLI reference](README.md)
# domainator sort_contigs

Sort contigs by their Domain_Search score.

```
Sort contigs by Domain_Search score.

This matches the output ordering used by domainator domain_search when --max_hits is set.
It is primarily intended for post-processing domain_search output written without
--max_hits, where the temporary hit_best_score annotation is not preserved in the
written GenBank file.
```

Give at least one of these outputs: `--genbank_out`, `--fasta_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | GenBank files. Use - to read standard input. |
| `-o`, `--genbank_out` | FILE |  | Write the records to this GenBank file. Use - for standard output. |
| `--fasta_out` | FILE |  | Write the records to this FASTA file. Use - for standard output. |
| `--max_hits` | INTEGER |  | Maximum number of contigs to return, keeping the best-scoring ones. If not supplied, return all contigs. |
| `--domains` | STRING ... |  | If supplied, sort by the best Domainator annotation score from domains with these names instead of Domain_Search annotation scores. Additive with --domains_file. |
| `--domains_file` | FILE |  | Like --domains, but read the domain names from this file, one per line. Additive with --domains. |
| `--databases` | STRING ... |  | When used with --domains, only consider Domainator annotations from these databases. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
