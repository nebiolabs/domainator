[CLI reference](README.md)
# domainator build_ssn_viewer

Build a Domainator Similarity Network Viewer bundle and HTML page.

```
Build a compact MST-based bundle for the standalone Domainator Similarity Network Viewer (DSNV).
```

Give at least one of these outputs: `--dsnv_out`, `--html_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE |  | Symmetric similarity matrix to convert into a Domainator Similarity Network Viewer bundle. Format can be tab-separated text or Domainator hdf5. |
| `-o`, `--dsnv_out` | FILE |  | Write the gzip-compressed JSON Domainator Similarity Network Viewer bundle (.dsnv) to this file. |
| `--html_out` | FILE |  | Write a standalone static HTML viewer page to this file. Without --embed_data, the page loads local bundle files. |
| `--embed_data` |  |  | Embed the input bundle JSON into --html_out so the viewer autoloads it on open. |
| `--name` | STRING |  | Optional display name recorded in the bundle. Defaults to the output file stem. |
| `--metadata` | FILE ... |  | Tab-separated metadata tables keyed by node id. |
| `--subset` | STRING ... |  | Only consider matrix labels in this list. Additive with --subset_file. |
| `--subset_file` | FILE |  | Text file containing matrix labels to retain. |
| `--color_by` | STRING |  | Record a default metadata field for node coloring in the viewer bundle. |
| `--label_by` | STRING |  | Record a default metadata field for node labels in the viewer bundle. |
| `--categorical` | STRING ... |  | Numeric metadata columns (for example integer cluster numbers) that the viewer should color as discrete categories instead of a gradient. Can be toggled per column in the viewer's color picker. |
| `--merge_impact_metric` | {product,min_child} | min_child | Metric the viewer uses for split events: min_child (node count of the smaller piece) or product (the two piece sizes multiplied). Recorded in the bundle as the default; switchable in the viewer. |
| `--max_output_gb` | NUMBER | 25.0 | Fail if a projected output file would exceed this many GB. Set to 0 to disable the guardrail. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
