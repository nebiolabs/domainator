[CLI reference](README.md)
# domainator structure_dist

Compute similarity or distance matrices between protein structures.

```
Calculates similarity or distance matrices between protein structures.

The structure-input counterpart of domainator seq_dist. Typically only -i is supplied, giving an
all-vs-all matrix of the input structures; a separate -r produces a rectangular matrix of
input structures against reference structures.

Every chain of every input is a row, whether or not it had a hit, so the row and column
labels are exactly the record ids that domainator structure_to_genbank writes for the same
structures, and the matrix feeds domainator build_ssn, domainator build_tree or domainator build_projection
directly. That does mean the input database is read end to end, unlike domainator structure_search.

Alignment direction is the opposite of domainator structure_search and domainator structure_domainate: the
input is the QUERY and the reference is the TARGET, so rows are inputs and columns are
references, and --evalue is computed against the reference set. This matches domainator seq_dist,
where diamond is run with the input as -q and the reference as -d.

--mode selects what goes in the cells. The score-based modes (score, bool, norm_score,
row_norm_score, score_dist, efi_score, efi_score_dist) are computed from the aligner's
alignment score and behave exactly as they do in domainator seq_dist. The structural modes put a
similarity that is already bounded in [0, 1] in the cell -- tmscore, lddt, fident -- and
their *_dist forms are simply 1 - that value, needing none of the min(row_max, col_max)
normalization that score_dist does. tmscore and lddt require --algorithm foldseek and
C-alpha coordinates in both databases; fident works with either backend.

Note that with --algorithm reseek there is no bit score: structure_lib records -log10(p)
in its place. That is monotonic in significance, so score, bool and the normalized modes
remain meaningful, but the values are not comparable to foldseek's and efi_score, which
assumes a bit score in its 2^(-score) term, is rejected.
```

Give at least one of these outputs: `--hdf_out`, `--dense_text_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | structures to compare: pdb/cif files (optionally gzipped), a directory searched recursively, a glob, or the prefix of a prebuilt aligner database. |
| `-r`, `--reference` | FILE ... |  | structures to compare the input against. Same accepted forms as --input. If not supplied, the input is compared against itself, giving a square all-vs-all matrix. Rows are inputs and columns are references. |
| `-o`, `--hdf_out` | FILE |  | Write the similarity (or, for *_dist modes, distance) matrix to this Domainator hdf5 file. Values are stored sparse; a matrix with exactly symmetric values is stored as its upper triangle. |
| `--dense_text_out` | FILE |  | Write the similarity (or, for *_dist modes, distance) matrix to this tab-separated text file, every cell included. Use - for standard output. |
| `-e`, `--evalue` | NUMBER | 0.001 | only pairs with an E-value strictly below this become non-zero cells. E-values are computed against the reference set, as they are in seq_dist. Raise it to fill in weaker relationships; lower it to shrink the matrix. |
| `--mode` | {bool,efi_score,efi_score_dist,fident,fident_dist,lddt,lddt_dist,norm_score,row_norm_score,score,score_dist,tmscore,tmscore_dist} | score | what kind of values should be in the matrix. Score-based, from the aligner's alignment score, exactly as in seq_dist: score: raw score, bool: 1 if a hit otherwise 0, norm_score: score/min(row_max, col_max), row_norm_score: score/row_max, score_dist: 1 - norm_score, efi_score: -log10\[2^(-score) * (row_len * col_len)\], efi_score_dist: 1 - (efi_score / min(row_max, col_max)). Structural, already in \[0,1\] so the distance forms are just 1 - value: tmscore, lddt, fident, tmscore_dist, lddt_dist, fident_dist. tmscore and lddt require --algorithm foldseek and C-alpha coordinates in both databases; fident works with either backend. efi_score and efi_score_dist require a real bit score, so they are rejected with --algorithm reseek. |
| `--symmetrize` | {none,min,max,mean} | none | reconcile the two directions of each pair of a square self-comparison, on the raw values before any --mode transform. none: leave the matrix as computed, which may be slightly asymmetric. max: the better of the two directions in both cells, the OR-symmetric convention used by --mst/--knn and by build_ssn. min: the worse of the two -- note a pair found in only one direction has a stored zero in the other, so min DELETES it, making this a mutual-hit filter. mean: the average of the two -- note a pair found in only one direction is HALVED, because the missing direction contributes zero. A symmetric result is stored as its upper triangle by --hdf_out, about half the size. Requires a square self-comparison (no -r). Cannot be combined with --mst/--knn, whose output is always max-symmetric. |
| `--threshold` | NUMBER | 0 | Round any values lower than or equal to this down to zero. This applies to the raw per-hit value before the --mode transform: the aligner's alignment score for the score-based modes, and the metric itself (0-1) for the tmscore/lddt/fident modes. |
| `--max_hits_per_query` | INTEGER |  | Include at most this many non-zero entries in the matrix for each input structure, keeping the highest-valued ones. Applied after --symmetrize. Applied before --mst and --knn, so a cap can split clusters that --mst would otherwise keep connected. Use --max_seqs to bound the search itself. |
| `--mst` |  |  | Keep only the edges of the maximum spanning tree (forest). The clusters (connected components) are the same as the full graph's, with the connections inside each cluster pruned to the minimum needed to keep it connected. Combine with --knn to also keep each node's nearest neighbors. Computed while streaming, so the full matrix is never held in memory. Requires a square self-comparison and --mode in {score, bool, efi_score, tmscore, lddt, fident}; bool selects on the alignment scores and then keeps the selected cells as 1. Best paired with --hdf_out. |
| `--knn` | INTEGER |  | Keep only OR-symmetric K-nearest-neighbor edges (K >= 1): an edge is kept if it is among the K best for either of its nodes. On its own this can split a cluster of the full graph into several; add --mst to keep the clusters connected. Computed while streaming, so the full matrix is never held in memory. Requires a square self-comparison and --mode in {score, bool, efi_score, tmscore, lddt, fident}; bool selects on the alignment scores and then keeps the selected cells as 1. Best paired with --hdf_out. |
| `--db_out` | FILE |  | when the input is structure files rather than a prebuilt database, also write the database built from them, using this path as its prefix. |
| `--hits_tsv_out` | FILE |  | write the aligner's own hit table here, for debugging or for cross-checking against a hand-run search. |
| `--progress` |  |  | Show a progress bar for long-running steps. |
| `--algorithm` | {foldseek,reseek} | foldseek | structural alignment backend to use. |
| `--foldseek_path` | PROGRAM |  | path to the foldseek executable. If not supplied, foldseek is looked up on PATH. |
| `--reseek_path` | PROGRAM |  | path to the reseek executable. If not supplied, reseek is looked up on PATH. |
| `--foldseek_args` | STRING |  | extra options for the foldseek search command, for options domainator does not wrap, as one quoted string split like a shell command line. Use the = form: --foldseek_args="--option value". Options that domainator sets itself are refused. |
| `--reseek_args` | STRING |  | extra options for the reseek search command, for options domainator does not wrap, as one quoted string split like a shell command line. Use the = form: --reseek_args="--option value". Options that domainator sets itself are refused. |
| `--cpu` | INTEGER | 0 | number of threads to give the backend. 0 lets the backend use all available cores. |
| `--tmp_dir` | DIR |  | parent directory for temporary databases and alignment files. Structural searches can need a lot of scratch space, so point this at a large filesystem for big inputs. |
| `--reseek_stats` | {family,superfamily,fold} | superfamily | reseek only: SCOP level that reseek's p-value model is calibrated against. This shifts p-values by orders of magnitude, so it changes which hits pass --evalue. Ignored by other backends. |
| `--reseek_sensitivity` | {fast,sensitive} | sensitive | reseek only: search sensitivity. Ignored by other backends. |
| `--alignment_type` | {0,1,2} |  | foldseek only: how to compute the alignment: 0 = 3Di only, 1 = TM-align, 2 = 3Di+AA (the default when unset). Type 1 requires C-alpha coordinates in both databases. reseek has no alignment-type option and rejects this. |
| `--metrics` | {tmscore,lddt,rmsd,prob} ... |  | additional per-hit structural metrics to compute and store as feature qualifiers. tmscore, lddt and rmsd require C-alpha coordinates in both databases. |
| `--max_seqs` | INTEGER |  | maximum number of hits to keep per reference, keeping the most significant ones. Defaults to the backend's own default: 1000 for foldseek, no cap for reseek. For foldseek this bounds the search itself, so raising it costs time but improves coverage; for reseek it bounds only the output, since reseek has no native cap and does the full search either way. A warning is printed whenever a reference saturates the cap. |
| `--device` | STRING | cpu | where to run the search. 'cpu' for a CPU-only search, or a CUDA device string ('cuda', 'cuda:0', 'cuda:1', ...) for a GPU-accelerated search. GPU search requires a recent foldseek. |
| `--max_output_gb` | NUMBER | 25.0 | Fail if a projected output file would exceed this many GB. Set to 0 to disable the guardrail. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
