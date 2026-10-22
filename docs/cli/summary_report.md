[CLI reference](README.md)
# domainator summary_report

Summarize the domain and taxonomy content of a dataset.

```
From a domainator-annotated genbank file, calculates contig stats and domain frequency stats
```

Give at least one of these outputs: `--text_out`, `--html_out`, `--json_out`, `--domains_tsv_out`, `--co_occurrence_fraction_out`, `--co_occurrence_count_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | name of input genbank files. Use - to read standard input. |
| `-o`, `--text_out` | FILE |  | Write the report as text to this file. Use - for standard output. |
| `--html_out` | FILE |  | Write the report as html to this file. Use - for standard output. |
| `--json_out` | FILE |  | Write the full summary as a single compact JSON object to this file. Use - for standard output. |
| `--domains_tsv_out` | FILE |  | Write a tab-separated table of domain frequencies to this file. Headers are: domain, database, count, description, avg score. |
| `--co_occurrence_fraction_out` | FILE |  | Write a dense text matrix of domain co-occurrence frequencies to this file. Every row is a domain, the values are the fraction of contigs containing the row-domain that also include the column-domain (so the matrix is asymmetric). |
| `--co_occurrence_count_out` | FILE |  | Write a dense text matrix of domain co-occurrence frequencies to this file. Every row is a domain, the values are the number of contigs containing both the row-domain and the column-domain (so the matrix is symmetric). |
| `--taxonomy` |  |  | Write a summary of the taxonomic origins of the contigs. |
| `--taxonomy_update` |  |  | If taxonomy database exists, check it against the version on the ncbi server and update if there is a newer version. |
| `--ncbi_taxonomy_path` | DIR | /tmp/ncbi_taxonomy_sean | Path to NCBI taxonomy database directory. Will be created and downloaded if it does not exist. |
| `--contigs` | STRING ... |  | Only consider contigs with these names. |
| `--contigs_file` | FILE |  | text file containing names of contigs to select. |
| `--domains` | STRING ... |  | calculate distance graphs for domains with these names |
| `--domains_file` | FILE |  | text file containing names of domains to calculate distance graphs for. |
| `--databases` | STRING ... |  | Ignore domain annotations not from these databases. default: consider all databases. |
| `--partial` | {include,exclude,only} | include | Which CDSs/proteins to count, by fragment status. include: all. exclude: skip fragments. only: only fragments. A CDS is a fragment if its location has '<' or '>', or it has a /partial qualifier. A protein record is a fragment if it has a UniProt ' (Fragment)' or MGnify 'FL=0' header, a 'Flags: Fragment' definition or NON_TER features (UniProt text entries), or a CDS/Protein feature with '<' or '>' spanning it. Protein records that are filtered out are skipped; on nucleotide contigs, filtered CDSs and their domains are not counted, but the contig is. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
