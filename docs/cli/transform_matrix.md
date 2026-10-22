[CLI reference](README.md)
# domainator transform_matrix

Convert score matrices into normalized scores or distances.

```
Transforms a score matrix into various kinds of normalized matrices

For matrices with any data type, can convert between hdf5 and dense text, or change the datatype.
```

Give at least one of these outputs: `--hdf_out`, `--dense_text_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE | required | name of input matrix file. |
| `--value_type` | {bool,efi_score,efi_score_dist,norm_score,row_norm_score,score,score_dist} |  | What the values of the input matrix are, if the file does not say. Overrides the type stored in the file. |
| `-o`, `--hdf_out` | FILE |  | Write the transformed matrix to this Domainator hdf5 file. Values are stored sparse; a matrix with exactly symmetric values is stored as its upper triangle. |
| `--dense_text_out` | FILE |  | Write the transformed matrix to this tab-separated text file, every cell included. Use - for standard output. |
| `--mode` | {none,score,bool,norm_score,row_norm_score,score_dist,efi_score,efi_score_dist} | none | what kind of values should be in the output matrix. none: preserve the input values and data type, score: raw score, bool: 1 if a hit otherwise 0, score_dist: 1 - (score / min(row_max, col_max)), norm_score: score/min(row_max, col_max), row_norm_score: score/row_max, efi_score: -log10[2^(-score) * (input_seq_length * reference_seq_length)], efi_score_dist: 1 - (efi_score / min(row_max, col_max)). |
| `--threshold` | NUMBER |  | Zero out all values less than or equal to this threshold after any mode transformation. |
| `--mst` |  |  | Keep only the edges of the maximum spanning tree (forest). The clusters (connected components) are the same as the full graph's, with the connections inside each cluster pruned to the minimum needed to keep it connected. Combine with --knn to also keep each node's nearest neighbors. Uses the post-transform values, except that --mode bool selects on the scores and then keeps the selected cells as 1; requires a square, non-distance matrix. |
| `--knn` | INTEGER |  | Keep only OR-symmetric K-nearest-neighbor edges (K >= 1): an edge is kept if it is among the K best for either of its nodes. On its own this can split a cluster of the full graph into several; add --mst to keep the clusters connected. Uses the post-transform values, except that --mode bool selects on the scores and then keeps the selected cells as 1; requires a square, non-distance matrix. |
| `--symmetrize` | {none,min,max,mean} | none | reconcile the two directions of each pair of a square self-comparison, on the raw values before any --mode transform. none: leave the matrix as computed, which may be slightly asymmetric. max: the better of the two directions in both cells, the OR-symmetric convention used by --mst/--knn and by build_ssn. min: the worse of the two -- note a pair found in only one direction has a stored zero in the other, so min DELETES it, making this a mutual-hit filter. mean: the average of the two -- note a pair found in only one direction is HALVED, because the missing direction contributes zero. A symmetric result is stored as its upper triangle by --hdf_out, about half the size. Applied before --mode; requires a square matrix with the same row and column labels, holding similarity (not *_dist) values. Cannot be combined with --mst/--knn, whose output is always max-symmetric. |
| `--max_output_gb` | NUMBER | 25.0 | Fail if a projected output file would exceed this many GB. Set to 0 to disable the guardrail. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
