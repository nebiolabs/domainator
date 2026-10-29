[CLI reference](README.md)
# domainator select_by_contig

Select whole contigs by their domain content and other properties.

```
Create a new genbank file with a subset of contigs from the input

From a domainator-annotated genbank file, extracts a subset of contigs based on selection criteria.

Selection criteria are combined as follows:

(contig length within specified bounds [--min_length, --max_length])
AND 
(contig CDS count within specified bounds [--min_cds, --max_cds])
AND 
(contig matches specified taxonomy filters)
AND
(contig passes the fragment filter [partial])
AND # INVERT applies to this list of criteria
(
   (contig name in specified list [contigs, contigs_file])
OR (contig name matches specified regex [contigs_regex])
OR (contig definition matches specified regex [definition_regex]) 
OR (contig sequence matches specified regex [sequence_regex])
OR (contig contains at least one domain with name in specified list [domains, domains_file])
OR (contig contains a combination of domains that match the specified boolean expression [domain_expr])
OR (contig contains no domains [unannotated])
OR (contig contains a qualifying search hit [has_search_hit])
OR (contig has a search-hit score within specified bounds [min_search_score, max_search_score])
)
AND
(total number of selected contigs is less than or equal to specified maximum [max_records])

If --invert is specified, the selection criteria in the middle expression (the stuff linked by OR) are inverted. Length, CDS count, taxonomy, fragment status, and max_records still apply,
but the contigs that would have been selected by the middle expression are instead discarded, and contigs that would have been discarded are selected.
```

Give at least one of these outputs: `--genbank_out`, `--fasta_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | Names of input genbank files. Use - to read standard input. |
| `-o`, `--genbank_out` | FILE |  | Write the records to this GenBank file. Use - for standard output. |
| `--fasta_out` | FILE |  | Write the records to this FASTA file. Use - for standard output. |
| `--contigs` | STRING ... |  | Select records with these IDs, OR with other record selectors; affected by --invert. Additive with --contigs_file; duplicate names are ignored. If neither source is supplied, no name criterion is applied. |
| `--contigs_file` | FILE |  | Select records with these IDs, OR with other record selectors; affected by --invert. Read names from a plain-text file (one per line), FASTA/GenBank record IDs, or HMM profile names. Additive with --contigs; duplicate names are ignored. Plain-text lines are stripped and blank lines ignored. An empty file adds no names; if no names are supplied directly either, the resulting list is explicitly empty, not an omitted criterion. |
| `--contigs_regex` | STRING |  | Only consider contigs with names matching this regex. |
| `--definition_regex` | STRING |  | only keep contigs with a definition that matches this regex. |
| `--sequence_regex` | STRING |  | only keep contigs with a sequence that matches this regex. |
| `--domains` | STRING ... |  | Select records containing these domain names, OR with other record selectors. Additive with --domains_file; duplicate names are ignored. If neither source is supplied, no name criterion is applied. |
| `--domains_file` | FILE |  | Select records containing these domain names, OR with other record selectors. Read names from a plain-text file (one per line), FASTA/GenBank record IDs, or HMM profile names. Additive with --domains; duplicate names are ignored. Plain-text lines are stripped and blank lines ignored. An empty file adds no names; if no names are supplied directly either, the resulting list is explicitly empty, not an omitted criterion. |
| `--domain_expr` | STRING |  | a boolean expression using operators & (AND), \| (OR), and ~(NOT), to specify a desired combination of domains |
| `--partial` | {include,exclude,only} | include | Select contigs by fragment status. include: all. exclude: drop fragments. only: keep only fragments. A protein record is a fragment if it has a UniProt ' (Fragment)' or MGnify 'FL=0' header, a 'Flags: Fragment' definition, NON_TER features, or a CDS/Protein feature with '<' or '>' spanning it. A nucleotide contig is a fragment if any of its CDSs has '<' or '>' in its location, or a /partial qualifier. Not affected by --invert. |
| `--annotation_type` | {domainator,domain_search,all} | domainator | Which annotations to select contigs by. domainator: domain annotations (from domainate). domain_search: best-hit annotations (from domain_search). all: both. |
| `--has_search_hit` |  |  | Also select existing domain_search hits, regardless of annotation_type, domain names, domain expressions, or databases. Hits must pass the e-value and score cutoffs where applicable. Participates in --invert where available. |
| `--databases` | STRING ... |  | Consider only domains from these databases when filtering by domain. default: all databases. |
| `--unannotated` |  |  | Select contigs with no domain annotations. |
| `-e`, `--evalue` | NUMBER |  | the evalue cutoff (max) for domain annotations to be considered when selecting by domain. Default: no cutoff. |
| `--include_taxids` | INTEGER ... |  | Include these NCBI taxids and their descendants (OR). |
| `--exclude_taxids` | INTEGER ... |  | Exclude these NCBI taxids and their descendants; exclusions take precedence. |
| `--taxonomy_expr` | STRING |  | Boolean expression over taxid lineages using &, \|, ~, and parentheses, e.g. "2 & ~1224". Mutually exclusive with include/exclude_taxids. |
| `--taxonomy_unknown` | {exclude,include,evaluate,error} | exclude | Unresolved taxonomy: exclude (default), include, evaluate against an empty lineage, or error. With ~2, exclude rejects unlabeled candidates; evaluate and include allow them. Only active with taxonomy criteria. Not inverted. |
| `--ncbi_taxonomy_path` | DIR | /tmp/ncbi_taxonomy_sean | NCBI taxonomy database directory; downloaded if missing and filtering is active. |
| `--taxonomy_update` |  |  | Check for an updated taxonomy database when filtering is active. |
| `--max_records` | INTEGER |  | If set then stop after returning this many contigs. |
| `--min_length` | INTEGER |  | Keep only input records with length >= this, in bp for nucleotide input or amino acids for protein input. Not affected by --invert. |
| `--max_length` | INTEGER |  | Keep only input records with length <= this, in bp for nucleotide input or amino acids for protein input. Not affected by --invert. |
| `--min_cds` | INTEGER |  | skip contigs with fewer CDSs than this. Pseudogenes are not counted. |
| `--max_cds` | INTEGER |  | skip contigs with more CDSs than this. Pseudogenes are not counted. |
| `--min_search_score` | NUMBER |  | Also select records with an existing Domain_Search annotation scoring >= this and passing e-value/database filters. Paired with --max_search_score as an inclusive range; OR with other selectors and affected by --invert. Also restricts --has_search_hit scoring, which ignores databases. |
| `--max_search_score` | NUMBER |  | Upper bound (inclusive) of the search-score selector range. OR with other selectors, not a restriction on records selected by other criteria. Also restricts --has_search_hit scoring, which ignores databases. |
| `--invert` |  |  | Invert the selection criteria. i.e. return contigs that do not match the selection criteria. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
