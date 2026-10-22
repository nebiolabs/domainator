[CLI reference](README.md)
# domainator seq_dist

Compute similarity scores between sequences or profiles.

```
Calculates similarity scores between protein sequences or HMM profiles

Typically the input and reference will be the same sequence file, creating a pairwise similarity matrix. 
An hmm file can also be used as a reference to make a table of profile scores for each input peptide.
An hmm file can also be used as both input and reference to compare profile similarity.
```

Give at least one of these outputs: `--hdf_out`, `--dense_text_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE | required | Input query file. Queries in the input file will be compared to the reference sequences/profiles. This can be a fasta, genbank, or hmm file. |
| `--input_format` | {fasta,genbank,hmm} |  | File format of the input, if it cannot be told from the file extension. |
| `-r`, `--reference` | FILE |  | Reference file. If not provided, the input file will be used as the reference, and the output will be a pairwise similarity matrix between the input sequences/profiles. This can be a fasta, genbank, or hmm file. |
| `--reference_format` | {fasta,genbank,hmm} |  | File format of the reference, if it cannot be told from the file extension. Ignored for some kinds of algorithms, for example, hmmer expects a .hmm reference file. |
| `-o`, `--hdf_out` | FILE |  | Write the similarity (or, for *_dist modes, distance) matrix to this Domainator hdf5 file. Values are stored sparse; a matrix with exactly symmetric values is stored as its upper triangle. |
| `--dense_text_out` | FILE |  | Write the similarity (or, for *_dist modes, distance) matrix to this tab-separated text file, every cell included. Use - for standard output. |
| `--threshold` | NUMBER | 0 | Round any scores lower than this down to zero (this applies to raw scores not processed/normalized). |
| `--max_hits_per_query` | INTEGER |  | Keep at most this many hits for each input sequence, the highest-scoring ones. This caps the search itself (diamond --max-target-seqs, or the hmmer hit list). Applied before --mst and --knn, so a cap can split clusters that --mst would otherwise keep connected. Default: keep all hits. |
| `--algorithm` | {diamond_f,diamond_d,diamond_mids,diamond_s,diamond_mors,diamond_vs,diamond_us,hmmer,hmmer_compare} | diamond_us | Which distance metric to use, diamond_*; f: fast, d: default, mids: mid-sensitive, s: sensitive, mors: more-sensitive, vs: very-sensitive, us: ultra-sensitive. default: diamond_us (ultra sensitive) |
| `--diamond_args` | STRING |  | Extra options for diamond blastp, as one quoted string split like a shell command line. Use the = form: --diamond_args="-c 4 -b 0.4 --target-indexed". Options that seq_dist sets itself (query, database, --outfmt, -p, --max-target-seqs, --max-hsps, --min-score, --verbose, and the sensitivity flags) are rejected; use the corresponding seq_dist option instead. seq_dist defaults --algo to 0 and --tmpdir to its own managed temporary directory; both can be replaced here. Useful for large all-vs-all runs, where diamond's seed-search scratch can exceed available disk: a smaller --block-size (-b) and/or more --index-chunks (-c) bound its peak scratch usage. |
| `--diamond_path` | PROGRAM |  | path to the diamond executable, if it is not on the PATH. |
| `--mst` |  |  | Keep only the edges of the maximum spanning tree (forest). The clusters (connected components) are the same as the full graph's, with the connections inside each cluster pruned to the minimum needed to keep it connected. Combine with --knn to also keep each node's nearest neighbors. Computed while streaming, so the full matrix is never held in memory. Requires the reference to be the same file as the input (a square, symmetric comparison) and --mode in {score, bool, efi_score}; bool selects on the scores and then keeps the selected cells as 1. Best paired with --hdf_out. |
| `--knn` | INTEGER |  | Keep only OR-symmetric K-nearest-neighbor edges (K >= 1): an edge is kept if it is among the K best for either of its nodes. On its own this can split a cluster of the full graph into several; add --mst to keep the clusters connected. Computed while streaming, so the full matrix is never held in memory. Requires the reference to be the same file as the input (a square, symmetric comparison) and --mode in {score, bool, efi_score}; bool selects on the scores and then keeps the selected cells as 1. Best paired with --hdf_out. |
| `--symmetrize` | {none,min,max,mean} | none | reconcile the two directions of each pair of a square self-comparison, on the raw values before any --mode transform. none: leave the matrix as computed, which may be slightly asymmetric. max: the better of the two directions in both cells, the OR-symmetric convention used by --mst/--knn and by build_ssn. min: the worse of the two -- note a pair found in only one direction has a stored zero in the other, so min DELETES it, making this a mutual-hit filter. mean: the average of the two -- note a pair found in only one direction is HALVED, because the missing direction contributes zero. A symmetric result is stored as its upper triangle by --hdf_out, about half the size. Applied to the raw scores, after --threshold and --max_hits_per_query. Requires the reference to be the same file as the input. Cannot be combined with --mst/--knn, whose output is always max-symmetric. |
| `--mode` | {score,bool,norm_score,row_norm_score,score_dist,efi_score,efi_score_dist} | score | what kind of values should be in the matrix. score: raw score, bool: 1 if a hit otherwise 0, score_dist: 1 - (score / min(row_max, col_max)), norm_score: score/min(row_max, col_max), efi_score: -log10\[2^(-score) * (input_seq_length * reference_seq_length)\], efi_score_dist: 1 - (efi_score / min(row_max, col_max)). |
| `--cpu` | INTEGER | 0 | how many cpu threads to use. 0: one per physical core. |
| `--progress` |  |  | Show a progress bar for long-running steps. |
| `--max_output_gb` | NUMBER | 25.0 | Fail if a projected output file would exceed this many GB. Set to 0 to disable the guardrail. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
