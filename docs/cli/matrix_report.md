[CLI reference](README.md)
# domainator matrix_report

Summarize a similarity matrix.

```
Reads a DataMatrix and generates a report suitable for helping with selecting edge score thresholds for similarity networks
```

Give at least one of these outputs: `--text_out`, `--html_out`, `--json_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE | required | name of input matrix file. |
| `-o`, `--text_out` | FILE |  | Write the report as text to this file. Use - for standard output. |
| `--html_out` | FILE |  | Write the report as html to this file. Use - for standard output. |
| `--json_out` | FILE |  | Write a compact JSON summary (edge stats, connected components, split-event series) to this file. Use - for standard output. |
| `--knn_projection` |  |  | Project how many edges 'build_ssn --threshold T --mst --knn K' would keep, for every threshold T and K from 1 to 25, and add a K slider to the HTML report. Disabled by default to reduce memory use. |
| `--merge_impact_metric` | {product,min_child} | min_child | Metric used for split-event plots and tables: 'min_child' emphasizes the smaller cluster breaking away, while 'product' emphasizes balanced large splits. |
| `--max_merge_events` | INTEGER | 500 | Maximum number of strongest merge events to embed in the interactive HTML threshold slider and split plot. A few more may be added on top: after the strongest are taken, the strongest event in each otherwise-empty 5% band of the threshold axis is added back, so the plot never leaves a stretch of its own axis blank. Use 0 to include all merge events. |
| `--progress` |  |  | Print live progress updates to stderr during long-running matrix_report stages. |
| `--profile_stages` |  |  | Print per-stage runtime timings to stderr for profiling large matrix_report runs. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
