[CLI reference](README.md)
# domainator select_by_cds

Extract contigs or neighborhoods around selected CDSs.

```
Extract contig regions around selected CDSs

Takes a domain-annotated genbank file and extracts regions of contigs based on contig or CDS name or presence of specified domains within individual CDSs.

Specify a range to extract around selected CDSs, if no range is specified, then only the selected CDSs will be returned and not any of their neighbors.

CDS names, domain names, domain expressions, unannotated, and search-hit presence
are OR selectors. With no selectors, all eligible CDSs are selected.
Contig IDs, strand, fragment status, and taxonomy are independent AND restrictions.
Invert negates only the combined selectors, never these restrictions. With no
selectors, invert warns and is ignored. Restrictions apply to the focus CDS,
not every feature in an expanded neighborhood.
```

Give at least one of these outputs: `--genbank_out`, `--fasta_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | GenBank files. Use - to read standard input. |
| `-o`, `--genbank_out` | FILE |  | Write the records to this GenBank file. Use - for standard output. |
| `--fasta_out` | FILE |  | Write the records to this FASTA file. Use - for standard output. |
| `--include_taxids` | INTEGER ... |  | Include these NCBI taxids and their descendants (OR). |
| `--exclude_taxids` | INTEGER ... |  | Exclude these NCBI taxids and their descendants; exclusions take precedence. |
| `--taxonomy_expr` | STRING |  | Boolean expression over taxid lineages using &, \|, ~, and parentheses, e.g. "2 & ~1224". Mutually exclusive with include/exclude_taxids. |
| `--taxonomy_unknown` | {exclude,include,evaluate,error} | exclude | Unresolved taxonomy: exclude (default), include, evaluate against an empty lineage, or error. With ~2, exclude rejects unlabeled candidates; evaluate and include allow them. Only active with taxonomy criteria. Not inverted. |
| `--ncbi_taxonomy_path` | DIR | /tmp/ncbi_taxonomy_sean | NCBI taxonomy database directory; downloaded if missing and filtering is active. |
| `--taxonomy_update` |  |  | Check for an updated taxonomy database when filtering is active. |
| `--kb_range` | NUMBER |  | How many kb upstream and downstream of the selected CDS to extract. Starting from the ends of the CDS. |
| `--kb_up` | NUMBER |  | How many kb upstream of the selected CDS to extract. Starting from the end of the CDS. |
| `--kb_down` | NUMBER |  | How many kb downstream of the selected CDS to extract. Starting from the end of the CDS. |
| `--cds_range` | INTEGER |  | How many CDSs upstream and downstream of the selected CDS to extract. Example, 1 would mean to return contigs including the selected CDSs and the immediate upstream and downstream CDSs. |
| `--cds_up` | INTEGER |  | How many CDSs upstream of the selected CDS to extract. Example, 1 would mean to return contigs including the selected CDSs and the immediate upstream CDS. |
| `--cds_down` | INTEGER |  | How many CDSs downstream of the selected CDS to extract. Example, 1 would mean to return contigs including the selected CDSs and the immediate downstream CDS. |
| `--whole_contig` |  |  | If set, then return the entire matching contigs. |
| `--contigs` | STRING ... |  | Only read input records with these IDs; not inverted by --invert. Additive with --contigs_file; duplicate names are ignored. If neither source is supplied, no name criterion is applied. |
| `--contigs_file` | FILE |  | Only read input records with these IDs; not inverted by --invert. Read names from a plain-text file (one per line), FASTA/GenBank record IDs, or HMM profile names. Additive with --contigs; duplicate names are ignored. Plain-text lines are stripped and blank lines ignored. An empty file adds no names; if no names are supplied directly either, the resulting list is explicitly empty, not an omitted criterion. |
| `--cds` | STRING ... |  | Include CDSs with these names (additive with --domains) |
| `--cds_file` | FILE |  | text file containing names of CDSs to include. |
| `--name_field` | STRING |  | CDS qualifier to prefer for name matching, falling back to gene_id, locus_tag, gene, then protein_id. Does not change coordinate-based region IDs. |
| `--domains` | STRING ... |  | Select CDSs containing these domain names, OR with other CDS selectors. Additive with --domains_file; duplicate names are ignored. If neither source is supplied, no name criterion is applied. |
| `--domains_file` | FILE |  | Select CDSs containing these domain names, OR with other CDS selectors. Read names from a plain-text file (one per line), FASTA/GenBank record IDs, or HMM profile names. Additive with --domains; duplicate names are ignored. Plain-text lines are stripped and blank lines ignored. An empty file adds no names; if no names are supplied directly either, the resulting list is explicitly empty, not an omitted criterion. |
| `--domain_expr` | STRING |  | a boolean expression using operators & (AND), \| (OR), and ~(NOT), to specify a desired combination of domains |
| `--databases` | STRING ... |  | Consider only domains from these databases when filtering by domain. default: all databases. |
| `--strand` | {f,r} |  | Only extract regions around CDSs on the specified strand. |
| `--partial` | {include,exclude,only} | include | Which CDSs to extract regions around, by fragment status. include: all. exclude: skip fragments. only: only fragments. A CDS is a fragment if its location has '<' or '>', or it has a /partial qualifier. A protein record is a fragment if it has a UniProt ' (Fragment)' or MGnify 'FL=0' header, a 'Flags: Fragment' definition or NON_TER features (UniProt text entries), or a CDS/Protein feature with '<' or '>' spanning it. Not affected by --invert. |
| `--unannotated` |  |  | Select CDSs with no domain annotations. |
| `--annotation_type` | {domainator,domain_search,all} | domainator | Which annotations to select CDSs by. domainator: domain annotations (from domainate). domain_search: best-hit annotations (from domain_search). all: both. |
| `--has_search_hit` |  |  | Also select existing domain_search hits, regardless of annotation_type, domain names, domain expressions, or databases. Hits must pass the e-value and score cutoffs where applicable. Participates in --invert where available. |
| `--pad` |  |  | If set, then ends of sequences will be padded so that the focus CDS aligns. See --anchor. |
| `--anchor` | {start,end,center} | start | Where to align the focus CDSs when --pad is set, relative to the focus CDS in its own orientation. start: align the first base of the focus CDSs. end: align their last base. center: align their midpoints. |
| `--keep_direction` |  |  | by default extracted regions will be flipped so that the focus cds is on the forward strand. Setting this option will keep the focus cds on whatever strand it started on. |
| `--dedup` | {none,id,sequence,sequence_rc} | id | Output uniqueness policy: none keeps all records; id drops repeated output IDs; sequence drops identical output sequences (ignoring case and terminal '*'); sequence_rc also collapses nucleotide reverse complements. The same ID must identify the same sequence, so sequence uniqueness also ensures unique IDs. Deduplication happens before padding. |
| `--max_region_overlap` | NUMBER | 1.0 | the maximum fractional of overlap between any two output regions. If >= 1, then no overlap filtering will be done. Regions are output in a greedy fashion based on CDS start site. New regions are output if less than this fraction of them overlaps with any previously output region. |
| `-e`, `--evalue` | NUMBER |  | the evalue cutoff (max) for domain annotations to be considered when selecting CDSs. Default: no cutoff. |
| `--invert` |  |  | Invert the combined CDS selectors, not contigs, strand, partial, or taxonomy. With no selectors, warn and ignore. |
| `--include_nucleic_acids` |  |  | Include nucleic-acid annotations with no associated CDS (those with cds_id='.') when applying CDS selection criteria. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
