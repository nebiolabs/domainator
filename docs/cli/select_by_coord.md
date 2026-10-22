[CLI reference](README.md)
# domainator select_by_coord

Extract contig regions listed in a coordinate table.

```
Extract contig regions using explicit coordinate ranges.

Reads a labeled coordinate table into memory, then iterates through the input
contigs and extracts one region for each matching row.

Required table columns:
    contig, start, end

Optional table columns:
    strand, file

Coordinates are 1-based and inclusive.

If end < start, the region is assumed to wrap around the origin of a circular
contig. The same coordinates on a linear contig raise an exception.

Extracted regions are named like those of the other Domainator tools
(see utils.region_id): contig_first:last, where first and last are the
coordinates of the first and last bases of the region as written, followed by
'rc' for reverse-strand regions. For example, start 10, end 20, strand - is
named contig_20:10rc. Whole contigs keep their names (contig_rc if reverse).

If a file column is supplied, it must contain basenames only (with extension).
Rows without a file value match any input file.
```

Give at least one of these outputs: `--genbank_out`, `--fasta_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | Names of input genbank files. Use - to read standard input. |
| `-o`, `--genbank_out` | FILE |  | Write the records to this GenBank file. Use - for standard output. |
| `--fasta_out` | FILE |  | Write the records to this FASTA file. Use - for standard output. |
| `--coords` | FILE | required | labeled coordinate table with columns contig, start, end, strand (optional), file (optional). |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
