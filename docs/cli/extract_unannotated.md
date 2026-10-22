[CLI reference](README.md)
# domainator extract_unannotated

Extract the regions of proteins not covered by any domain annotation.

```
Extract unannotated regions

Takes a genbank file which has been annotated using the domainator program, looks for the specified domains, and extracts regions not containing them.

This tool is currently only supports protein sequences, support for nucleic acids may be added someday.

can write a fasta or a genbank file
```

Give at least one of these outputs: `--genbank_out`, `--fasta_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | Genbank filenames. Use - to read standard input. |
| `-o`, `--genbank_out` | FILE |  | Write the records to this GenBank file. Use - for standard output. |
| `--fasta_out` | FILE |  | Write the records to this FASTA file. Use - for standard output. |
| `--contigs` | STRING ... |  | only extract regions from contigs with ids in this list. Additive with --contigs_file. default: all contigs. |
| `--contigs_file` | FILE |  | only extract regions from contigs with ids listed in this file (one per line). Additive with --contigs. default: all contigs. |
| `-e`, `--evalue` | NUMBER |  | the evalue cutoff (max) for domain annotations to be considered. default: no cutoff. |
| `--min_score` | NUMBER |  | the score cutoff (min) for domain annotations to be considered. default: no cutoff. |
| `--databases` | STRING ... |  | Only consider domains from these databases. default: all databases. |
| `--annotation_type` | {domainator,domain_search} | domainator | Which annotations count as annotated. domainator: domain annotations (from domainate). domain_search: best-hit annotations (from domain_search). |
| `--keep_direction` |  |  | by default extracted regions will be flipped so that the source CDS will be on the forward strand. Set this option to not reverse complement anything. In the rare case of a split CDS with parts on different strands, the strand of the first part will be considered. |
| `--keep_name` |  |  | by default, extracted regions will be renamed to note their original coordinates and whether they were reverse complemented. Set this option to keep contigs with their original names. WARNING: if you aren't careful, this can easily result in duplicate names. |
| `--largest` |  |  | Extract only the largest unannotated region from each protein, providing it passes the length cutoff specified by --min_length |
| `--min_length` | INTEGER |  | Extract unannotated regions of at least this length from each protein. default: no cutoff. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
