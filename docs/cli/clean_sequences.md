[CLI reference](README.md)
# domainator clean_sequences

Clean sequence names and residues, and filter out sequences with non-canonical residues.

```
Clean and filter sequences in genbank or fasta files

This module provides functionality to:
- Clean sequence names for compatibility with newick format
- Deduplicate sequence names to ensure uniqueness
- Strip non-canonical residues from sequence ends
- Filter out sequences containing non-canonical residues
```

Give at least one of these outputs: `--genbank_out`, `--fasta_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | Names of input genbank files. Use - to read standard input. |
| `--alphabet` | {nucleotide,protein} | protein | Whether the sequences in FASTA inputs are protein or nucleotide. GenBank records declare their own molecule type. |
| `-o`, `--genbank_out` | FILE |  | Write the records to this GenBank file. Use - for standard output. |
| `--fasta_out` | FILE |  | Write the records to this FASTA file. Use - for standard output. |
| `--clean_name` |  |  | if true then make the name compatible with newick format. |
| `--deduplicate_names` |  |  | if true then rename sequences such that output sequences have unique names. |
| `--strip` |  |  | if true then strip non-canonical residues or nucleotides from the ends of sequences (depending on input type). Including * characters. |
| `--filter_by_aa` |  |  | if true then throw out sequences containing non-canonical residues or nucleotides (depending on input type). Including * characters. (applied after strip) |
| `--log_out` | FILE |  | Write the cleaning statistics to this file. If not supplied, writes them to stderr. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
