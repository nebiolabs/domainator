[CLI reference](README.md)
# domainator hmmer_report

Tabulate properties of each HMM profile.

```
Writes record-wise information from hmm files into tab-separated or html tables

Amino acid, DNA, and RNA profiles are all supported; use --alphabet to report which
alphabet each profile uses.
```

Give at least one of these outputs: `--tsv_out`, `--html_out`, `--json_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | hmm filenames. Use - to read standard input. |
| `-o`, `--tsv_out` | FILE |  | Write the table to this tab-separated file. Use - for standard output. |
| `--html_out` | FILE |  | Write the table in html format to this file. Use - for standard output. |
| `--json_out` | FILE |  | Write streaming NDJSON (one compact JSON object per profile) to this file. Use - for standard output. |
| `--html_max_height` | INTEGER |  | Max height in pixels to set the html output to. |
| `--column_names` | STRING ... |  | If supplied, then this list will be used instead of the default column names. |
| `--columns` | columns... |  | The columns to report, in order. A column with arguments takes them after ':', as in 'qualifier:CDS:product'. Can be given more than once. |

## Columns

| Column | Description |
|---|---|
| `source` | the name of the hmm file. |
| `acc` | the profile's accession. |
| `desc` | the profile's description. |
| `length` | the length of the profile (the number of nodes). |
| `consensus` | the profile's consensus residue line, if set. |
| `alphabet` | the profile's alphabet (amino, DNA, or RNA). |
| `append:name:type:value` | a column with this name and type (str, int or float) holding this value in every row. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
