[CLI reference](README.md)
# domainator extract_domains

Extract annotated domains as their own records.

```
Extract domain-annotated regions from protein or nucleotide records.

Reads annotated GenBank records and writes the selected regions as GenBank or FASTA.
Output retains the input sequence type: nucleotide input produces nucleotide regions,
not translated peptides or automatically complete CDSs.

Extension is measured in amino acids for protein input and bases for nucleotide input.
Nucleotide extension is not rounded to codons or confined to the source CDS.
Individual extractions follow annotation locations, including compound locations.

--combine and --splice operate across selected annotations in each input record, not
separately per CDS. They can span different genes or strands and do not guarantee a
valid coding sequence. --combine includes intervening sequence and rejects circular
records; --splice joins selected regions rather than reconstructing a gene.

For protein-domain output from nucleotide records, use extract_peptides first, then
extract_domains. For whole CDSs or nucleotide neighborhoods, use select_by_cds.
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
| `--contigs` | STRING ... |  | Only read input records with these IDs; restricts where extraction occurs. Additive with --contigs_file; duplicate names are ignored. If neither source is supplied, no name criterion is applied. |
| `--contigs_file` | FILE |  | Only read input records with these IDs; restricts where extraction occurs. Read names from a plain-text file (one per line), FASTA/GenBank record IDs, or HMM profile names. Additive with --contigs; duplicate names are ignored. Plain-text lines are stripped and blank lines ignored. An empty file adds no names; if no names are supplied directly either, the resulting list is explicitly empty, not an omitted criterion. |
| `-e`, `--evalue` | NUMBER |  | Extract only domain/search-hit annotations with e-value <= this, before extension, combination, or splicing. Default: no cutoff. |
| `--min_score` | NUMBER |  | Extract only domain/search-hit annotations with score >= this, before extension, combination, or splicing. Default: no cutoff. |
| `--domains` | STRING ... |  | Extract annotations with these domain names; search hits can also be selected by --has_search_hit. Additive with --domains_file; duplicate names are ignored. If neither source is supplied, no name criterion is applied. |
| `--domains_file` | FILE |  | Extract annotations with these domain names; search hits can also be selected by --has_search_hit. Read names from a plain-text file (one per line), FASTA/GenBank record IDs, or HMM profile names. Additive with --domains; duplicate names are ignored. Plain-text lines are stripped and blank lines ignored. An empty file adds no names; if no names are supplied directly either, the resulting list is explicitly empty, not an omitted criterion. |
| `--databases` | STRING ... |  | Domain databases to use for domain extraction. default: all databases. |
| `--annotation_type` | {domainator,domain_search,all} | domainator | Which annotations to extract. domainator: domain annotations (from domainate). domain_search: best-hit annotations (from domain_search). all: both. |
| `--has_search_hit` |  |  | Also select existing domain_search hits, regardless of annotation_type, domain names, domain expressions, or databases. Hits must pass the e-value and score cutoffs where applicable. Participates in --invert where available. |
| `--partial` | {include,exclude,only} | include | Which domains to extract, by the CDS or protein they are on, by fragment status. include: all. exclude: skip fragments. only: only fragments. A CDS is a fragment if its location has '<' or '>', or it has a /partial qualifier. A protein record is a fragment if it has a UniProt ' (Fragment)' or MGnify 'FL=0' header, a 'Flags: Fragment' definition or NON_TER features (UniProt text entries), or a CDS/Protein feature with '<' or '>' spanning it. Contig-level nucleotide hits are never on fragments. |
| `--extend_range` | INTEGER |  | Extend each side by this many amino acids for protein input or bases for nucleotide input. Not rounded to codons or restricted to CDS boundaries. Sets the default of --extend_up and --extend_down. |
| `--extend_up` | INTEGER |  | Extend upstream by this many amino acids for protein input or bases for nucleotide input. With --combine, extends the enclosing span. Default: --extend_range, or 0. |
| `--extend_down` | INTEGER |  | Extend downstream by this many amino acids for protein input or bases for nucleotide input. With --combine, extends the enclosing span. Default: --extend_range, or 0. |
| `--combine` |  |  | Extract one enclosing span of all selected annotations per input record, including intervening sequence. Not grouped by CDS; may span different genes or strands. Does not guarantee a coding sequence. Rejects circular records. By default, extracts each annotation separately. Cannot be used with --splice. |
| `--splice` |  |  | if set, then when multiple target domains occur in the same contig, the intermediate sequence will be deleted and the extracted domains from the contig glued together, in the order they are read on the strand of the domain with the lowest start coordinate. Applies across selected annotations per input record, not per CDS; does not guarantee biological exon splicing or a coding sequence. By default, each hit domain will be extracted separately. Cannot be used with --combine. |
| `--keep_direction` |  |  | by default extracted regions will be flipped so that the domain (or the domain with the lowest start coordinate when --combine or --splice is active) is on the forward strand. Set this option to not reverse complement anything. In the rare case of a split domain with parts on different strands, the strand of the first part will be considered. |
| `--keep_parent_id` |  |  | Use the parent nucleotide or protein record ID instead of a coordinate-based region ID. This can produce duplicate IDs for different extracted sequences. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
