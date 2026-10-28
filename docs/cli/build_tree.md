[CLI reference](README.md)
# domainator build_tree

Build a UPGMA tree from a distance or similarity matrix.

```
Generates trees from distance or similarity matrices, or reads them from Newick files

Distance matrices (data types ending in _dist, or dense text with no data type) are used
as is. Similarity matrices (e.g. score, efi_score) are first converted to distances with
1 - (score / min(row_max, col_max)), so a pair with no score has distance 1.
The two directions of each pair are reconciled by taking the smaller distance.

Instead of a matrix, --newick reads a tree built elsewhere (for example a maximum-likelihood
tree), to write it with metadata, colors, and an alignment.

Outputs are a Newick tree, a Cytoscape xgmml file (laid out as a rectangular tree), and an
interactive Domainator Viewer page, which can show a multiple sequence alignment beside the
tree.
```

Give at least one of these outputs: `--newick_out`, `--xgmml_out`, `--html_out`, `--dvb_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE |  | Square matrix of pairwise distances or similarity scores. Format can be tab-separated text or Domainator hdf5. Give this or --newick. |
| `--newick` | FILE |  | Read the tree from this Newick file instead of building one from a matrix. Numeric internal node labels are read as support values. |
| `--distance` | {auto,precomputed,score_dist} | auto | How to turn matrix values into distances. precomputed: the values are distances. score_dist: the values are similarity scores, distance = 1 - (score / min(row_max, col_max)), so a pair with no score has distance 1. auto: precomputed if the matrix data type is a distance (ends in _dist) or not recorded (e.g. tab-separated text), otherwise score_dist. |
| `--algorithm` | {upgma} | upgma | Which clustering algorithm to use. |
| `--metadata` | FILE ... |  | Tab separated files of sequence metadata, keyed by their first column. These data are added to the xgmml and html outputs. |
| `--subset` | STRING ... |  | Only include these leaves. Additive with --subset_file. |
| `--subset_file` | FILE |  | Text file of leaf names to include. |
| `--msa` | FILE |  | An aligned FASTA file of the leaves' sequences, keyed by name, shown beside the tree in the html output. |
| `--color_by` | STRING |  | Color the leaves in the xgmml, and at first in the html viewer, by this column of the metadata table. |
| `--color_table` | FILE |  | tab separated file with two columns and no header, columns are: annotation, hex color. For example: CCDB   cc0000. Requires --color_by. |
| `--color_table_out` | FILE |  | tab separated file with two columns and no header, columns are: annotation, hex color. Written after the color table is updated with new colors, for example if using --color_by, but not supplying an external color table. |
| `--palette_offset` | INTEGER | 0 | Skip this many colors at the start of the palette before assigning colors to --color_by values. Ignored with --color_table. |
| `--max_color_groups` | INTEGER |  | Give distinct colors only to this many of the largest --color_by groups, all remaining groups share a neutral gray. |
| `--label_by` | STRING |  | Label the leaves in the html viewer by this metadata column at first. |
| `--categorical` | STRING ... |  | Numeric metadata columns (for example integer cluster numbers) that the html viewer should color as discrete categories instead of a gradient. |
| `--name` | STRING |  | The tree's name, shown in the html viewer. Defaults to the html (or dvb) output's file name. |
| `-o`, `--newick_out` | FILE |  | Write the tree to this Newick file. Use - for standard output. |
| `--xgmml_out` | FILE |  | Write the tree to this Cytoscape xgmml file, laid out as a rectangular tree. |
| `--html_out` | FILE |  | Write the Domainator Viewer, showing this tree, to this HTML file. |
| `--dvb_out` | FILE |  | Write the gzip-compressed JSON Domainator Viewer bundle (.dvb) to this file. The viewer opens it, and saves sessions in the same format. |
| `--max_output_gb` | NUMBER | 25.0 | Fail if a projected output file would exceed this many GB. Set to 0 to disable the guardrail. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
