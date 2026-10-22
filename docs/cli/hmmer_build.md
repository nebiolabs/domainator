[CLI reference](README.md)
# domainator hmmer_build

Build an HMM profile from a multiple sequence alignment.

```
Build an HMM profile from a multiple sequence alignment (MSA)
    
    Allows the user to specify the ACC, NAME, and DESC fields of the HMM profile.
    Some of these options are not available in the standard hmmbuild tool from HMMER.
```

Give at least one of these outputs: `--hmm_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE | required | Path of input msa. Acceptable formats are the same as for hmmbuild. Use - to read standard input. |
| `-o`, `--hmm_out` | FILE |  | Write the profile to this hmm file. Writes gzip-compressed output if the path ends in '.gz'; '.bgz' is refused because HMMER cannot read BGZF. Output to stdout is never compressed. Use - for standard output. |
| `--name` | STRING | required | Name of the HMM profile. |
| `--acc` | STRING |  | Accession of the HMM profile. |
| `--desc` | STRING |  | Description of the HMM profile. |
| `--alphabet` | {amino,dna,rna} |  | Alphabet of the input MSA. If not supplied, inferred from the input, which can guess wrong on short or ambiguous alignments. |
| `--window_length` | INTEGER |  | Max expected hit length (the MAXL field of the profile), used by nhmmer for windowed E-values. Nucleotide alphabets only. If not supplied, hmmer picks one. Takes precedence over --window_beta, so only one of the two may be given. |
| `--window_beta` | NUMBER |  | Tail mass from which to derive the window length. Nucleotide alphabets only. Only one of --window_length and --window_beta may be given. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
