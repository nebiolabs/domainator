[CLI reference](README.md)
# domainator extract_domains

Extract annotated domains as their own records.

```
Extract domain-annotated regions

Takes a genbank file which has been annotated using the domainator program, looks for the specified domains, and extracts them from the contigs.

This tool is typically used with protein sequences, but can also be used with nucleic acid sequences.
Handling of nucleic acids is less stable and will may change in future updates, so use with caution.

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
| `--contigs` | STRING ... |  | only extract domains from contigs with ids in this list. Additive with --contigs_file. default: all contigs. |
| `--contigs_file` | FILE |  | only extract domains from contigs with ids listed in this file (one per line). Additive with --contigs. default: all contigs. |
| `-e`, `--evalue` | NUMBER |  | Extract only domain/search-hit annotations with e-value <= this, before extension, combination, or splicing. Default: no cutoff. |
| `--min_score` | NUMBER |  | Extract only domain/search-hit annotations with score >= this, before extension, combination, or splicing. Default: no cutoff. |
| `--domains` | STRING ... |  | Include domains with these names |
| `--domains_file` | FILE |  | text file containing names of domains to include. |
| `--databases` | STRING ... |  | Domain databases to use for domain extraction. default: all databases. |
| `--annotation_type` | {domainator,domain_search,all} | domainator | Which annotations to extract. domainator: domain annotations (from domainate). domain_search: best-hit annotations (from domain_search). all: both. |
| `--has_search_hit` |  |  | Also select existing domain_search hits, regardless of annotation_type, domain names, domain expressions, or databases. Hits must pass the e-value and score cutoffs where applicable. Participates in --invert where available. |
| `--partial` | {include,exclude,only} | include | Which domains to extract, by the CDS or protein they are on, by fragment status. include: all. exclude: skip fragments. only: only fragments. A CDS is a fragment if its location has '<' or '>', or it has a /partial qualifier. A protein record is a fragment if it has a UniProt ' (Fragment)' or MGnify 'FL=0' header, a 'Flags: Fragment' definition or NON_TER features (UniProt text entries), or a CDS/Protein feature with '<' or '>' spanning it. Contig-level nucleotide hits are never on fragments. |
| `--extend_range` | INTEGER |  | extract this many additional residues (or bases) on each side of the hit region. Sets the default of --extend_up and --extend_down. |
| `--extend_up` | INTEGER |  | extract this many additional residues (or bases) upstream of the hit region. Default: --extend_range, or 0. |
| `--extend_down` | INTEGER |  | extract this many additional residues (or bases) downstream of the hit region. Default: --extend_range, or 0. |
| `--combine` |  |  | if set, then when multiple target domains occur in the same contig, an envelope containing all of the hit domains will be extracted. By default, each hit domain will be extracted separately. This option will raise an exception on circular contigs, because it would be ambiguous. Cannot be used with --splice. |
| `--splice` |  |  | if set, then when multiple target domains occur in the same contig, the intermediate sequence will be deleted and the extracted domains from the contig glued together, in the order they are read on the strand of the domain with the lowest start coordinate. By default, each hit domain will be extracted separately. Cannot be used with --combine. |
| `--keep_direction` |  |  | by default extracted regions will be flipped so that the domain (or the domain with the lowest start coordinate when --combine or --splice is active) is on the forward strand. Set this option to not reverse complement anything. In the rare case of a split domain with parts on different strands, the strand of the first part will be considered. |
| `--keep_parent_id` |  |  | Use the parent nucleotide or protein record ID instead of a coordinate-based region ID. This can produce duplicate IDs for different extracted sequences. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
