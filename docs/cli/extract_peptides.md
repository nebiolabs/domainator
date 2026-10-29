[CLI reference](README.md)
# domainator extract_peptides

Extract translated CDSs as protein records, keeping their domain annotations.

```
Converts annotated CDSs into annotated protein sequences

Takes a domainator-annotated genbank file 
Creates peptide genbank or fasta files with the domainator annotations from the CDSs found in the original
genbank file

By default, extracts peptides from all CDSs in the genbank file.

If any criteria are supplied, then only peptides from CDSs that meet the criteria will be extracted.

CDS names, domain names, domain expressions, unannotated, and search-hit presence
are OR selectors. With no selectors, all eligible CDSs are selected.
Contig IDs, strand, fragment status, and taxonomy are independent AND restrictions.
Invert negates only the combined selectors, never these restrictions. With no
selectors, invert warns and is ignored.
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
| `--contigs` | STRING ... |  | only look for peptides on contigs with ids in this list. Additive with --contigs_file. default: all contigs. |
| `--contigs_file` | FILE |  | only look for peptides on contigs with ids listed in this file (one per line). Additive with --contigs. default: all contigs. |
| `--domains` | STRING ... |  | Only extract peptides containing these domains |
| `--domains_file` | FILE |  | text file containing names of domains, where peptides containing these domains will be extracted. |
| `--domain_expr` | STRING |  | Boolean domain expression using &, \|, ~, and parentheses. Additive with other CDS selectors. |
| `--databases` | STRING ... |  | Consider only domains from these databases when filtering by domain. default: all databases. |
| `--unannotated` |  |  | extract peptides from CDSs with no domain annotations. |
| `--annotation_type` | {domainator,domain_search,all} | domainator | Which annotations to select peptides by. domainator: domain annotations (from domainate). domain_search: best-hit annotations (from domain_search). all: both. |
| `--has_search_hit` |  |  | Also select existing domain_search hits, regardless of annotation_type, domain names, domain expressions, or databases. Hits must pass the e-value and score cutoffs where applicable. Participates in --invert where available. |
| `--strand` | {f,r} |  | Only extract peptides from CDSs on the specified strand. |
| `--partial` | {include,exclude,only} | include | Which CDSs to extract peptides from, by fragment status. include: all. exclude: skip fragments. only: only fragments. A CDS is a fragment if its location has '<' or '>', or it has a /partial qualifier. A protein record is a fragment if it has a UniProt ' (Fragment)' or MGnify 'FL=0' header, a 'Flags: Fragment' definition or NON_TER features (UniProt text entries), or a CDS/Protein feature with '<' or '>' spanning it. Not affected by --invert. |
| `--invert` |  |  | Invert the combined CDS selectors, not contigs, strand, partial, or taxonomy. With no selectors, warn and ignore. |
| `--cds` | STRING ... |  | only look for peptides on cdss with these names |
| `--cds_file` | FILE |  | only look for peptides on cdss with names listed in this text file. |
| `--keep_parent_id` |  |  | Use the original contig ID instead of the resolved CDS name for peptide IDs. This can produce duplicate IDs. |
| `--name_field` | STRING |  | The CDS qualifier to prefer for CDS name matching and peptide IDs. If not supplied, then qualifiers will be checked in order of: \["protein_id", "gene_id", "locus_tag", "gene"\] |
| `-e`, `--evalue` | NUMBER |  | Maximum e-value of domain/search-hit annotations used for CDS selection and carried onto peptide records. Default: no cutoff. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
