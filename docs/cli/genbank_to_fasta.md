[CLI reference](README.md)
# domainator genbank_to_fasta

Convert GenBank files to FASTA.

```
Convert a GenBank file to a FASTA file
```

Give at least one of these outputs: `--fasta_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | Genbank filenames. Use - to read standard input. |
| `-o`, `--fasta_out` | FILE |  | Write the records to this FASTA file. Use - for standard output. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
