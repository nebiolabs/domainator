[CLI reference](README.md)
# domainator split_by_cluster

Split a sequence or profile file into one file per cluster of a similarity matrix.

```
Split a sequence or profile file into connected components derived from a matrix.

Given one input file and one symmetric matrix, threshold the matrix at the
requested lower bound, find connected components, and write one output file per
cluster. Clusters smaller than the requested minimum size are combined into a
single small_clusters output file.
```

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE | required | Input genbank, fasta, hmm, or cm file to split by cluster. |
| `--matrix` | FILE | required | Matrix input in dense HDF5, sparse HDF5, or dense text format. |
| `--output_dir` | DIR | required | Directory to write cluster files into. Existing target files cause an error. |
| `--threshold` | NUMBER | 0 | Keep only edges with weights strictly above this value when finding connected components. |
| `--min_cluster_size` | INTEGER | 1 | Minimum cluster size required to write a dedicated cluster file. |
| `--write_small_clusters` |  |  | If set, write clusters smaller than --min_cluster_size into small_clusters.<ext>. Otherwise those records or models are skipped. |
| `--pad_on_search_hits` |  |  | For GenBank inputs, center each contig on its first Domain_Search annotation and pad records within each output file so the search hits align. Raises an error if any record lacks Domain_Search. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
