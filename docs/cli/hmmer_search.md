[CLI reference](README.md)
# domainator hmmer_search

Select HMM profiles by their similarity to reference profiles.

```
Filter an hmm file based on similarity to reference hmm profiles

If you want a distance matrix, use domainator seq_dist with .hmm files for both inputs and outputs.

If you want something more like hhsearch or blast, use domainator hmmer_compare, which uses 
the same algorithm, but a different interface (can output the same input profile multiple times).


Adapted from pseudocode in:
Steinegger, Martin, Markus Meier, Milot Mirdita, Harald Vöhringer, Stephan J. Haunsberger, and Johannes Söding. “HH-Suite3 for Fast Remote Homology Detection and Deep Protein Annotation.” BMC Bioinformatics 20, no. 1 (September 14, 2019): 473. https://doi.org/10.1186/s12859-019-3019-7.

and

Söding, Johannes. “Protein Homology Detection by HMM–HMM Comparison.” Bioinformatics 21, no. 7 (April 1, 2005): 951–60. https://doi.org/10.1093/bioinformatics/bti125.

Amino acid, DNA, and RNA profiles are all supported, but the -i and -r sides must use
the same alphabet (DNA and RNA count as different alphabets). Profile-profile scores are
on the scale of whichever alphabet is in use, so the default --min_score differs
between amino acid and nucleotide profiles.
```

Give at least one of these outputs: `--hmm_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | Input query files. One or more hmm text files with one or more hmmer3 profiles. If -i and -r are different, -i will typically be the bigger of the two, as it is what is being filtered. Use - to read standard input. |
| `-r`, `--reference` | FILE ... | required | Reference files. One or more hmm text files with one or more hmmer3 profiles. If -i and -r are different, -r will typically be the smaller of the two, as it is used as a filtering criterion. |
| `-o`, `--hmm_out` | FILE |  | Write the hit profiles to this hmm file. Writes gzip-compressed output if the path ends in '.gz'; '.bgz' is refused because HMMER cannot read BGZF. Output to stdout is never compressed. Use - for standard output. |
| `--min_score` | NUMBER |  | Report alignments with scores greater than or equal to this. Default: 10 for amino acid profiles, and the equivalent point on the nucleotide score scale for DNA and RNA profiles (about 4.63). A value supplied here is used as-is, without rescaling. |
| `--max_hits` | INTEGER |  | The maximum number of hmms returned by the search, keeping the best-scoring ones by the bitscore of their best scoring profile. If not supplied, return all hits passing the score cutoff. |
| `--cpu` | INTEGER | 0 | how many cpu threads to use. 0: one per physical core. |
| `--max_output_gb` | NUMBER | 25.0 | Fail if a projected output file would exceed this many GB. Set to 0 to disable the guardrail. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
