[CLI reference](README.md)
# domainator filter_domains

Filter domain annotations by e-value, overlap, database, or name.

```
Filter domain anntations

Takes a genbank file which has been annotated through the domainator program and filters domainator annotations,
based on percentage overlap between domains and evalues of domains.
```

Give at least one of these outputs: `--genbank_out`, `--fasta_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | The genbank filenames to be filtered. Use - to read standard input. |
| `-o`, `--genbank_out` | FILE |  | Write the records to this GenBank file. Use - for standard output. |
| `--fasta_out` | FILE |  | Write the records to this FASTA file. Use - for standard output. |
| `-e`, `--evalue` | NUMBER |  | threshold (max) E value for the domain hit. Default: no cutoff. |
| `--max_overlap` | NUMBER | 1.0 | the maximum fraction of overlap between domains to be included in the annotated genbank. |
| `--contigs` | STRING ... |  | only annotate contigs with ids in this list. Additive with --contigs_file. default: annotate all contigs. |
| `--contigs_file` | FILE |  | only annotate contigs with ids listed in this file (one per line). Additive with --contigs. default: annotate all contigs. |
| `--domains` | STRING ... |  | Remove all domains NOT matching these names |
| `--domains_file` | FILE |  | text file containing names of domains to include. (all domains not matching this list will be removed). |
| `--keep_databases` | STRING ... |  | Remove annotations from all databases except these. Applied before other filters, such as domains, evalue, and max_overlap. By default will keep all annotations from all databases. |
| `--databases` | STRING ... |  | Only filter annotations from these databases, according to evalue, domains, and max_overlap; annotations from other databases are left as they are. By default, annotations from all databases are filtered. |
| `--annotation_type` | {domainator,domain_search,all} | all | Which type of annotations to filter (domainator, domain_search, all). |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
