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
OR (contig contains at least one domain with name in specified list [target_domains, target_domains_file]) 
OR (contig contains a combination of domains that match the specified boolean expression [domain_expr])
OR (contig contains no domains [unannotated])
OR (contig contains a qualifying search hit [has_search_hit])
)
AND
(total number of selected contigs is less than or equal to specified maximum [first])

If --invert is specified, the selection criteria in the middle expression (the stuff linked by OR) are inverted. That is, length and first are still applied, 
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
| `--contigs` | STRING ... |  | Only consider contigs with these names. |
| `--contigs_file` | FILE |  | text file containing names of contigs to select. |
| `--contigs_regex` | STRING |  | Only consider contigs with names matching this regex. |
| `--definition_regex` | STRING |  | only keep contigs with a definition that matches this regex. |
| `--sequence_regex` | STRING |  | only keep contigs with a sequence that matches this regex. |
| `--domains` | STRING ... |  | Include domains with these names |
| `--domains_file` | FILE |  | text file containing names of domains to include. |
| `--domain_expr` | STRING |  | a boolean expression using operators & (AND), \| (OR), and ~(NOT), to specify a desired combination of domains |
| `--partial` | {include,exclude,only} | include | Select contigs by fragment status. include: all. exclude: drop fragments. only: keep only fragments. A protein record is a fragment if it has a UniProt ' (Fragment)' or MGnify 'FL=0' header, a 'Flags: Fragment' definition, NON_TER features, or a CDS/Protein feature with '<' or '>' spanning it. A nucleotide contig is a fragment if any of its CDSs has '<' or '>' in its location, or a /partial qualifier. Not affected by --invert. |
| `--annotation_type` | {domainator,domain_search,all} | domainator | Which annotations to select contigs by. domainator: domain annotations (from domainate). domain_search: best-hit annotations (from domain_search). all: both. |
| `--has_search_hit` |  |  | Also select existing domain_search hits, regardless of annotation_type, domain names, domain expressions, or databases. Hits must pass the e-value and score cutoffs where applicable. Participates in --invert where available. |
| `--databases` | STRING ... |  | Consider only domains from these databases when filtering by domain. default: all databases. |
| `--unannotated` |  |  | Select contigs with no domain annotations. |
| `-e`, `--evalue` | NUMBER |  | the evalue cutoff (max) for domain annotations to be considered when selecting by domain. Default: no cutoff. |
| `--include_taxids` | INTEGER ... |  | Space separated list of taxids to include |
| `--exclude_taxids` | INTEGER ... |  | Space separated list of taxids to exclude |
| `--taxonomy_expr` | STRING |  | A boolean expression over taxids using operators & (AND), \| (OR), ~ (NOT), and parentheses, e.g. "2 & ~1224" (within Bacteria but not Proteobacteria). A taxid is true for a contig when it is in the contig's lineage. Mutually exclusive with --include_taxids/--exclude_taxids. |
| `--ncbi_taxonomy_path` | DIR | /tmp/ncbi_taxonomy_sean | Path to NCBI taxonomy database directory. Will be created and downloaded if it does not exist. |
| `--taxonomy_update` |  |  | If taxonomy database exists, check it against the version on the ncbi server and update if there is a newer version. |
| `--max_records` | INTEGER |  | If set then stop after returning this many contigs. |
| `--min_length` | INTEGER |  | skip contigs smaller than this. In units of bp (aa for protein input). |
| `--max_length` | INTEGER |  | skip contigs larger than this. In units of bp (aa for protein input). |
| `--min_cds` | INTEGER |  | skip contigs with fewer CDSs than this. Pseudogenes are not counted. |
| `--max_cds` | INTEGER |  | skip contigs with more CDSs than this. Pseudogenes are not counted. |
| `--min_search_score` | NUMBER |  | Skip contigs with a domain_search annotation score smaller than this. |
| `--max_search_score` | NUMBER |  | Skip contigs with a domain_search annotation score larger than this. |
| `--invert` |  |  | Invert the selection criteria. i.e. return contigs that do not match the selection criteria. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
