[CLI reference](README.md)
# domainator build_tree

Build a UPGMA tree from a distance or similarity matrix.

```
Generates trees from distance or similarity matrices

Distance matrices (data types ending in _dist, or dense text with no data type) are used
as is. Similarity matrices (e.g. score, efi_score) are first converted to distances with
1 - (score / min(row_max, col_max)), so a pair with no score has distance 1.
The two directions of each pair are reconciled by taking the smaller distance.
```

Give at least one of these outputs: `--newick_out`, `--xgmml_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE | required | Square matrix of pairwise distances or similarity scores. Format can be tab-separated text or Domainator hdf5. |
| `--distance` | {auto,precomputed,score_dist} | auto | How to turn matrix values into distances. precomputed: the values are distances. score_dist: the values are similarity scores, distance = 1 - (score / min(row_max, col_max)), so a pair with no score has distance 1. auto: precomputed if the matrix data type is a distance (ends in _dist) or not recorded (e.g. tab-separated text), otherwise score_dist. |
| `--algorithm` | {upgma} | upgma | Which clustering algorithm to use. |
| `--metadata` | FILE ... |  | Tab separated files of sequence metadata. These data will be added to the xgmml output. |
| `-o`, `--newick_out` | FILE |  | Write the tree to this Newick file. Use - for standard output. |
| `--xgmml_out` | FILE |  | Write the tree to this Cytoscape xgmml file. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
