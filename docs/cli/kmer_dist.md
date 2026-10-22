[CLI reference](README.md)
# domainator kmer_dist

Compute exact k-mer similarity scores between sequence records.

```
Calculates similarity scores between sequence records using exact k-mer set comparisons.

Typically the input and reference will be the same sequence file, creating a pairwise similarity matrix.
Protein inputs use protein k-mers. Nucleotide inputs use translated CDS protein k-mers by default,
or strand-invariant nucleotide k-mers when --nucleotide_kmers is supplied.
```

Give at least one of these outputs: `--hdf_out`, `--dense_text_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE | required | Input sequence file. This can be a fasta or genbank file. |
| `--input_format` | {fasta,genbank} |  | File format of the input, if it cannot be told from the file extension. |
| `-r`, `--reference` | FILE |  | Reference file. If not provided, the input file will be used as the reference. |
| `--reference_format` | {fasta,genbank} |  | File format of the reference, if it cannot be told from the file extension. |
| `--alphabet` | {nucleotide,protein} | protein | Whether the sequences in FASTA inputs are protein or nucleotide. GenBank records declare their own molecule type. |
| `--gene_call` | {all,unannotated} |  | When activated, new CDS annotations will be added with Prodigal in Metagenomic mode before translated-CDS k-mer extraction. If 'all', then any existing CDS annotations will be deleted and all contigs will be re-annotated. If 'unannotated', then only contigs without CDS annotations will be annotated. |
| `--nucleotide_kmers` |  |  | Use strand-invariant nucleotide k-mers instead of translated CDS protein k-mers for nucleotide inputs. |
| `--ksize` | INTEGER |  | K-mer size. Defaults to 5 for protein/translated-CDS comparisons and 15 for nucleotide-kmer comparisons. |
| `--metric` | {max_containment,min_containment,jaccard} | max_containment | Similarity metric to compute from shared k-mer counts. |
| `--mode` | {score,bool,score_dist} | score | Output matrix values. score: raw similarity, bool: 1 if a similarity is non-zero, score_dist: 1 - (score / min(row_max, col_max)). |
| `-o`, `--hdf_out` | FILE |  | Write the similarity (or, for score_dist, distance) matrix to this Domainator hdf5 file. Values are stored sparse; a matrix with exactly symmetric values is stored as its upper triangle. |
| `--dense_text_out` | FILE |  | Write the similarity (or, for score_dist, distance) matrix to this tab-separated text file, every cell included. Use - for standard output. |
| `--threshold` | NUMBER | 0 | Round any similarities lower than or equal to this down to zero before mode transforms. |
| `--max_hits_per_query` | INTEGER |  | Include at most this many non-zero entries in the matrix for each input sequence, keeping the highest-valued ones. Applied before --mst and --knn, so a cap can split clusters that --mst would otherwise keep connected. |
| `--mst` |  |  | Keep only the edges of the maximum spanning tree (forest). The clusters (connected components) are the same as the full graph's, with the connections inside each cluster pruned to the minimum needed to keep it connected. Combine with --knn to also keep each node's nearest neighbors. Computed while streaming, so the full matrix is never held in memory. Selects on the k-mer similarity. Requires the reference to be the same file as the input (a square, symmetric comparison) and --mode in {score, bool}. Best paired with --hdf_out. |
| `--knn` | INTEGER |  | Keep only OR-symmetric K-nearest-neighbor edges (K >= 1): an edge is kept if it is among the K best for either of its nodes. On its own this can split a cluster of the full graph into several; add --mst to keep the clusters connected. Computed while streaming, so the full matrix is never held in memory. Selects on the k-mer similarity. Requires the reference to be the same file as the input (a square, symmetric comparison) and --mode in {score, bool}. Best paired with --hdf_out. |
| `--cpu` | INTEGER | 0 | How many cpu threads to use. Default: use all physical cores. |
| `--progress` |  |  | Show a progress bar for long-running steps. |
| `--max_output_gb` | NUMBER | 25.0 | Fail if a projected output file would exceed this many GB. Set to 0 to disable the guardrail. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
