[CLI reference](README.md)
# domainator hmmer_compare

Score and align pairs of HMM profiles.

```
Aligns and calculates alignment scores between hmmer3 profiles 

Kind of like hhsearch except much slower and for hmmer3 profiles instead of hhsuite profiles.

Adapted from pseudocode in:
Steinegger, Martin, Markus Meier, Milot Mirdita, Harald Vöhringer, Stephan J. Haunsberger, and Johannes Söding. “HH-Suite3 for Fast Remote Homology Detection and Deep Protein Annotation.” BMC Bioinformatics 20, no. 1 (September 14, 2019): 473. https://doi.org/10.1186/s12859-019-3019-7.

and

Söding, Johannes. "Protein Homology Detection by HMM–HMM Comparison." Bioinformatics 21, no. 7 (April 1, 2005): 951–60. https://doi.org/10.1093/bioinformatics/bti125.

Amino acid, DNA, and RNA profiles are all supported, but the -i and -r sides must use
the same alphabet (DNA and RNA count as different alphabets).
```

Give at least one of these outputs: `--tsv_out`, `--alignments_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | Input query files. One or more hmm text files with one or more hmmer3 profiles. Use - to read standard input. |
| `-r`, `--reference` | FILE ... | required | Reference files. One or more hmm text files with one or more hmmer3 profiles. |
| `--min_score` | NUMBER | 0 | Report alignments with scores greater than or equal to this. Scores are on the score scale of the input alphabet: DNA and RNA scores run at roughly half the amino acid scale, so a threshold tuned on proteins is about twice as strict on nucleotide profiles. |
| `--max_hits_per_query` | INTEGER |  | Include at most this many of the top hits for each query. Default: Include all hits. |
| `-o`, `--tsv_out` | FILE |  | Write the scores to this tab-separated file. Use - for standard output. |
| `--alignments_out` | FILE |  | Also write each hit's alignment to this text file. Use - for standard output. |
| `--cpu` | INTEGER | 0 | how many cpu threads to use. 0: one per physical core. |
| `--max_output_gb` | NUMBER | 25.0 | Fail if a projected output file would exceed this many GB. Set to 0 to disable the guardrail. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
