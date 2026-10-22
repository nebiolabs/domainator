[CLI reference](README.md)
# domainator compare_contigs

Compute similarities between gene neighborhoods from their domain content and order.

```
Calculates similarity metrics for gene neighborhoods

Takes a genbank file annotated with domainator and calculates similarities between contigs and writes reports.

ji, ai, and dss metrics are similar to and inspired by the metrics used in the BigScape method:
Navarro-Muñoz, Jorge C., Nelly Selem-Mojica, Michael W. Mullowney, Satria A. Kautsar, James H. Tryon, Elizabeth I. Parkinson, Emmanuel L. C. De Los Santos, et al. "A Computational Framework to Explore Large-Scale Biosynthetic Diversity from Large-Scale Genomic Data." Nature Chemical Biology 16, no. 1 (January 2020): 60–68. https://doi.org/10.1038/s41589-019-0400-9.
```

Give at least one of these outputs: `--hdf_out`, `--dense_text_out`, `--genbank_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | Genbank filenames. Use - to read standard input. |
| `--genbank_out` | FILE |  | Write the contigs to this GenBank file, sorted hierarchically by similarity. Use - for standard output. |
| `--name_by_order` |  |  | Will rename contigs in the genbank file such that they have a prefix that can be sorted to restore the sorted order. |
| `--linkage` | {average,ward,complete,single,weighted} | average | Hierarchical clustering linkage method used to sort contigs for the genbank output (-o). Distances are 1 - score/max(score); pairs with no score have distance 1. |
| `--leaf_ordering` | {greedy,optimal,none} | greedy | How to order the leaves of the clustering tree for the genbank output (-o). greedy: flip subtrees so that the most similar contigs meet at each junction; fast and usually about as good as optimal. optimal: scipy's optimal leaf ordering; very slow (minutes to hours) beyond a few thousand contigs. none: plain dendrogram order. |
| `--skip_sort_memory_check` |  |  | Don't fail early when sorting the contigs for the genbank output (-o) is estimated to need more memory than is available. |
| `--contigs` | STRING ... |  | only annotate contigs with ids in this list. Additive with --contigs_file. default: all contigs. |
| `--contigs_file` | FILE |  | only annotate contigs with ids listed in this file (one per line). Additive with --contigs. default: all contigs. |
| `--max_hits_per_query` | INTEGER |  | After combining all metric scores, keep this many of the top hits for each query and consider all other hits to be zero. Applied before --mst and --knn, so a cap can split clusters that --mst would otherwise keep connected. Default: keep all hits. |
| `--threshold` | NUMBER | 0 | Round any final scores less than or equal to this down to zero. This is primarily useful for making sparse outputs smaller. |
| `--mst` |  |  | Keep only the edges of the maximum spanning tree (forest). The clusters (connected components) are the same as the full graph's, with the connections inside each cluster pruned to the minimum needed to keep it connected. Combine with --knn to also keep each node's nearest neighbors. Computed while streaming, so the full matrix is never held in memory. Best paired with --hdf_out. |
| `--knn` | INTEGER |  | Keep only OR-symmetric K-nearest-neighbor edges (K >= 1): an edge is kept if it is among the K best for either of its nodes. On its own this can split a cluster of the full graph into several; add --mst to keep the clusters connected. Computed while streaming, so the full matrix is never held in memory. Best paired with --hdf_out. |
| `--ji` | NUMBER | 0.5 | weighting for jaccard index metric, a comparison of all distinct types of domains in the contigs. When two contigs have no domains, their ji is 0.0. |
| `--ai` | NUMBER | 0.5 | weighting for adjacency index metric, a comparison of shared domain pairs in the contigs. When two contigs have 0 or 1 domains, their ai is 0.0. |
| `--databases` | STRING ... |  | Domain databases to use for domain content comparison. default: all databases. |
| `-o`, `--hdf_out` | FILE |  | Write the contig similarity score matrix to this Domainator hdf5 file. Values are stored sparse; a matrix with exactly symmetric values is stored as its upper triangle. |
| `--dense_text_out` | FILE |  | Write the contig similarity score matrix to this tab-separated text file, every cell included. Use - for standard output. |
| `--cpu` | INTEGER | 0 | How many cpu workers to use. Default: use all available physical cores. |
| `--progress` |  |  | Show a progress bar for long-running steps. |
| `--max_output_gb` | NUMBER | 25.0 | Fail if a projected output file would exceed this many GB. Set to 0 to disable the guardrail. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
