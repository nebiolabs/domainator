[CLI reference](README.md)
# domainator extract_unannotated

Extract the regions of proteins not covered by any domain annotation.

```
Extract unannotated regions

Takes a genbank file which has been annotated using the domainator program, looks for the specified domains, and extracts regions not containing them.

This tool supports only protein sequences.

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
| `--include_taxids` | INTEGER ... |  | Include these NCBI taxids and their descendants (OR). |
| `--exclude_taxids` | INTEGER ... |  | Exclude these NCBI taxids and their descendants; exclusions take precedence. |
| `--taxonomy_expr` | STRING |  | Boolean expression over taxid lineages using &, \|, ~, and parentheses, e.g. "2 & ~1224". Mutually exclusive with include/exclude_taxids. |
| `--taxonomy_unknown` | {exclude,include,evaluate,error} | exclude | Unresolved taxonomy: exclude (default), include, evaluate against an empty lineage, or error. With ~2, exclude rejects unlabeled candidates; evaluate and include allow them. Only active with taxonomy criteria. Not inverted. |
| `--ncbi_taxonomy_path` | DIR | /tmp/ncbi_taxonomy_sean | NCBI taxonomy database directory; downloaded if missing and filtering is active. |
| `--taxonomy_update` |  |  | Check for an updated taxonomy database when filtering is active. |
| `--contigs` | STRING ... |  | only extract regions from contigs with ids in this list. Additive with --contigs_file. default: all contigs. |
| `--contigs_file` | FILE |  | only extract regions from contigs with ids listed in this file (one per line). Additive with --contigs. default: all contigs. |
| `-e`, `--evalue` | NUMBER |  | Count only domain/search-hit annotations with e-value <= this when finding unannotated gaps. Default: no cutoff. |
| `--min_annotation_score` | NUMBER |  | Count only domain/search-hit annotations with score >= this when finding unannotated gaps. Raising this can enlarge the gaps. Default: no cutoff. |
| `--databases` | STRING ... |  | Only consider domains from these databases. default: all databases. |
| `--annotation_type` | {domainator,domain_search,all} | domainator | Which annotations count as annotated. domainator: domain annotations (from domainate). domain_search: best-hit annotations (from domain_search). all: both. |
| `--keep_direction` |  |  | Has no effect: this tool supports only protein sequences and never reverse-complements them. Retained for compatibility with extract_domains. |
| `--keep_parent_id` |  |  | Use the parent protein ID instead of a coordinate-based region ID. This can produce duplicate IDs for different extracted sequences. |
| `--largest` |  |  | Extract only the largest unannotated region from each protein, providing it passes the length cutoff specified by --min_region_length |
| `--min_region_length` | INTEGER |  | Extract only unannotated gaps of at least this many amino acids, before choosing --largest. Not a cutoff on the input protein length. Default: no cutoff. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
