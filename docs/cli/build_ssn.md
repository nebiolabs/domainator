[CLI reference](README.md)
# domainator build_ssn

Build a sequence similarity network for Cytoscape or the Domainator Viewer.

```
Generates Sequence Similarity Networks

Build a sequence similarity network and do analysis related to that: write it for Cytoscape
(--xgmml_out) or for the Domainator Viewer (--html_out, --dvb_out).
```

Give at least one of these outputs: `--xgmml_out`, `--html_out`, `--dvb_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE | required | pairwise similarity (or distance) scores. In the network, nodes will be derived from row and column names, and edges will be added between pairs of rows and columns meeting the threshold criteria. Format can be tab-separated text, or Domainator hdf5. |
| `--threshold` | NUMBER | 0 | keep only edges with weights strictly above this value. This should be >= 0. |
| `--metadata` | FILE ... |  | tab separated files of sequence metadata. |
| `--subset` | STRING ... |  | Only consider matrix labels in this list. Additive with --subset_file. |
| `--subset_file` | FILE |  | text file containing matrix labels to retain. |
| `--color_by` | STRING |  | Color the points in the output image based on this column of the metadata table. |
| `--color_table` | FILE |  | tab separated file with two columns and no header, columns are: annotation, hex color. For example: CCDB   cc0000 |
| `--max_color_groups` | INTEGER |  | Give distinct colors only to this many of the largest --color_by groups, all remaining groups share a neutral gray. By default every group gets a color, cycling over a palette of 64 distinguishable colors. |
| `--palette_offset` | INTEGER | 0 | Skip this many colors at the start of the palette before assigning colors to --color_by groups. For example, if a network colored by one column used the first 5 colors, --palette_offset 5 colors another column without reusing them. Ignored with --color_table. |
| `--mark_mst` |  |  | Give every edge a boolean 'MST' property saying whether it belongs to the maximum spanning tree. Useful for styling or filtering the MST backbone within a larger network in Cytoscape. |
| `--color_table_out` | FILE |  | tab separated file with two columns and no header, columns are: annotation, hex color. Written after the color table is updated with new colors, for example if using --color_by, but not supplying an external color table. |
| `-o`, `--xgmml_out` | FILE |  | Write the network to this Cytoscape xgmml file. |
| `--html_out` | FILE |  | Write the Domainator Viewer, showing this network, to this HTML file. The viewer shows the network's maximum spanning tree above --threshold, whatever --mst and --knn say. |
| `--dvb_out` | FILE |  | Write the gzip-compressed JSON Domainator Viewer bundle (.dvb) to this file. The viewer opens it (domainator viewer), and saves sessions in the same format. |
| `--merge_impact_metric` | {product,min_child} | min_child | How the viewer's split chart measures a merge: min_child (node count of the smaller piece) or product (the two piece sizes multiplied). Switchable in the viewer. |
| `--label_by` | STRING |  | Label the nodes in the html viewer by this metadata column at first. |
| `--categorical` | STRING ... |  | Numeric metadata columns (for example integer cluster numbers) that the html viewer should color as discrete categories instead of a gradient. |
| `--name` | STRING |  | The network's name, shown in the html viewer. Defaults to the html (or dvb) output's file name. |
| `--cluster` |  |  | Creates a new metadata column called 'SSN_cluster', replacing that column from the metadata file if already present. Finds clusters by following connectivity, assigns each cluster a distinct integer, based on size rank, bigger clusters get smaller numbers. |
| `--cluster_tsv_out` | FILE |  | Path to write a tsv file mapping each sequence to its cluster. If set, then --cluster is automatically activated. If --no_cluster_header is not set, then the file will have the header: contig cluster. |
| `--no_cluster_header` |  |  | If set, then the tsv file will not have a header. Only relevant if --cluster_tsv_out is set. |
| `--max_output_gb` | NUMBER | 25.0 | Fail if a projected output file would exceed this many GB. Set to 0 to disable the guardrail. |
| `--mst` |  |  | Keep only the edges of the maximum spanning tree (forest). The clusters (connected components) are the same as the full graph's, with the connections inside each cluster pruned to the minimum needed to keep it connected. Combine with --knn to also keep each node's nearest neighbors. |
| `--knn` | INTEGER |  | Keep only OR-symmetric K-nearest-neighbor edges (K >= 1): an edge is kept if it is among the K best for either of its nodes. On its own this can split a cluster of the full graph into several; add --mst to keep the clusters connected. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
