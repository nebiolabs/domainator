[CLI reference](README.md)
# domainator format_systems

Write system models as a model file, a table, or the HMM profiles they need.

```
Inspect system models, and write them as a model file, a table, or the profiles they need

Reads system models (--system_models, and/or one model defined on the command line with
--domain_expr and a window, as in find_systems) and writes:

- --yaml_out: a native system model file. Models that fit the short form are written in it, so
  this is how a command-line model becomes a file:

    domainator format_systems --domain_expr "X & Y & (Q | R)" --window_cds 10 --system_name XYZ_QR --yaml_out xyz.yaml

- --tsv_out: one row per model, component and reference, with its thresholds, to see what a
  model (for example a converted one) means;
- --hmm_out: the HMM profiles that the selected models name, from their reference_files, ready
  for domainate -r.

--systems/--systems_file and --system_categories select models, so --systems with --yaml_out
copies one model out of a collection into a file of its own, to edit.
```

Give at least one of these outputs: `--yaml_out`, `--tsv_out`, `--hmm_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `--system_models` | FILE ... |  | System model files (YAML with a domainator_systems: key), or directories holding a systems.yaml. See docs/file_formats.md. |
| `--domain_expr` | STRING |  | A system defined on the command line: a boolean expression over domain (reference) names that must hold within one window, using & (and), \| (or), ~ (not, within the window), parentheses and < (upstream of, in system orientation). database:name limits a name to one database. Needs --window_cds or --window_kb. |
| `--window_cds` | INTEGER |  | Window size for --domain_expr, in CDSs (see --window_mode). |
| `--window_kb` | NUMBER |  | Window size for --domain_expr, in kb (see --window_mode). With --window_cds, both must hold. |
| `--window_mode` | {span,chain,anchor} | span | How --window_cds and --window_kb are measured. span: the whole system, first member to last. chain: the gap between consecutive members. anchor: the distance of every member from one anchor member. |
| `--system_name` | STRING | domain_expr | Name of the --domain_expr system (its database is command_line). |
| `--same_strand` |  |  | Members of the --domain_expr system must all be on one strand. |
| `--systems` | STRING ... |  | Only the models with these names (name or database:name). Additive with --systems_file; duplicate names are ignored. If neither source is supplied, no name criterion is applied. |
| `--systems_file` | FILE |  | Only the models with these names (name or database:name). Read names from a plain-text file (one per line), FASTA/GenBank record IDs, or HMM profile names. Additive with --systems; duplicate names are ignored. Plain-text lines are stripped and blank lines ignored. An empty file adds no names; if no names are supplied directly either, the resulting list is explicitly empty, not an omitted criterion. |
| `--system_categories` | STRING ... |  | Only the models in these categories (e.g. defense or PKS). With --systems, a model must match both. |
| `--yaml_out` | FILE |  | Write the selected models as a native system model file. They must come from one collection; a command-line model is written in a collection named after this file. |
| `--tsv_out` | FILE |  | Write one row per model, component and reference, with its thresholds. |
| `--hmm_out` | FILE |  | Write the HMM profiles the selected models name, from their reference files, for domainate -r. Components that only accept hits from one database match hits labeled with that database, which domainate takes from the file name, so name the file after it. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
