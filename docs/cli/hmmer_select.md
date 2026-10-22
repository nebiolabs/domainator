[CLI reference](README.md)
# domainator hmmer_select

Extract a subset of profiles from HMM files.

```
Find profiles in an hmm file using text-based searches

Search through an hmm file to find profiles with names or descriptions matching a specified regex or string, or containing a specified string.

All specifications are treated as "OR" with relation to each other. So any profile that matches any of the specifications will be returned.

Amino acid, DNA, and RNA profiles are all supported. Because HMMER locks an hmm file to
the alphabet of its first profile, every profile written to one output file must share a
single alphabet, so all --input files must agree. --alphabet restricts the search to input
files of one alphabet, which is how a mixed set of files can be split into per-alphabet
outputs. Unlike the text criteria, --alphabet is a restriction rather than another "OR"
term: it narrows whatever the text criteria matched.
```

Give at least one of these outputs: `--hmm_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | names of input hmm files. Use - to read standard input. |
| `-o`, `--hmm_out` | FILE |  | Write the selected profiles to this hmm file. Writes gzip-compressed output if the path ends in '.gz'; '.bgz' is refused because HMMER cannot read BGZF. Output to stdout is never compressed. Use - for standard output. |
| `--regex` | STRING ... |  | What regex(es) to search. |
| `--exact` | STRING ... |  | What exact string(s) to search. Target must match exactly with no additions at either end. |
| `--contains` | STRING ... |  | What string(s) to search. Target must contain this as a substring. |
| `--field` | {name,desc,acc,all} | all | whether to search the NAME field, the DESC field, the ACC field, or all three. |
| `--case_sensitive` |  |  | If set, then the search is case sensitive. Otherwise it is case insensitive. |
| `--alphabet` | {amino,dna,rna} |  | Only read input files whose profiles use this alphabet, skipping the rest. Unlike the text criteria this narrows the selection rather than adding to it, and it can be used on its own to split a mixed set of input files by alphabet. Because an hmm file can only hold one alphabet, this selects whole files. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
