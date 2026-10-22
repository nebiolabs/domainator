[CLI reference](README.md)
# domainator deduplicate_genbank

Reduce redundancy by clustering sequences (cd-hit, usearch, diamond, or exact hashes).

```
Remove redundant sequences from a genbank file

Runs a clustering algorithm, such as cdhit, usearch, or diamond on sequences from a genbank (or fasta) file to reduce redundancy.
The "hash" algorithm is a faster way to deduplicate exact duplicates, by hashing the sequences and eliminating hash duplicates.
```

Give at least one of these outputs: `--genbank_out`, `--fasta_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | GenBank or FASTA files. |
| `-o`, `--genbank_out` | FILE |  | Write the records to this GenBank file. Use - for standard output. |
| `--fasta_out` | FILE |  | Write the records to this FASTA file. Use - for standard output. |
| `--log_out` | FILE |  | Write the log to this file. If not supplied, writes it to stderr. |
| `--alphabet` | {nucleotide,protein} | protein | Whether the sequences in FASTA inputs are protein or nucleotide. GenBank records declare their own molecule type. |
| `--algorithm` | {cd-hit,usearch,hash,diamond} | cd-hit | Which clustering algorithm to use. |
| `--cdhit_path` | PROGRAM |  | path to the cd-hit executable, if it is not on the PATH. |
| `--cdhit_est_path` | PROGRAM |  | path to the cd-hit-est executable, if it is not on the PATH. |
| `--usearch_path` | PROGRAM |  | path to the usearch executable, if it is not on the PATH. |
| `--diamond_path` | PROGRAM |  | path to the diamond executable, if it is not on the PATH. |
| `--identity` | NUMBER | required | Identity threshold (between 0 and 1), passed to cd-hit as the -c parameter, to usearch as the -id parameter, or to diamond as --approx-id after converting to percent. |
| `--prefix_count` |  |  | In the output genbank file, will prepend X- to the names of contigs, where X is the number of contigs from the original file collapsed into this centroid. |
| `--suffix_count` |  |  | In the output genbank file, will append -X to the names of contigs, where X is the number of contigs from the original file collapsed into this centroid. |
| `--cpu` | INTEGER | 0 | The number of threads to use for deduplication. 0: one per physical core. |
| `--cdhit_args` | STRING |  | extra options for cd-hit or cd-hit-est, as one quoted string split like a shell command line. Use the = form: --cdhit_args="-s 0.9". Options that domainator sets itself are refused. |
| `--usearch_args` | STRING |  | extra options for usearch, as one quoted string split like a shell command line. Use the = form: --usearch_args="-s 0.9". Options that domainator sets itself are refused. |
| `--diamond_args` | STRING |  | extra options for diamond cluster, as one quoted string split like a shell command line. Use the = form: --diamond_args="-s 0.9". Options that domainator sets itself are refused. |
| `--both_strands` |  |  | When set for nucleotide inputs, consider both strands for clustering/hash deduplication. |
| `--cluster_tsv_out` | FILE |  | If supplied, then write a tab separated table with columns: representative, contigs. |
| `--cluster_sep` | STRING | ; | Separator to separate individual contig names in the second column of the cluster table. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
