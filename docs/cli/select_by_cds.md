[CLI reference](README.md)
# domainator select_by_cds

Extract contigs or neighborhoods around selected CDSs.

```
Extract contig regions around selected CDSs

Takes a domain-annotated genbank file and extracts regions of contigs based on contig or CDS name or presence of specified domains within individual CDSs.

Specify a range to extract around selected CDSs, if no range is specified, then only the selected CDSs will be returned and not any of their neighbors.

If CDS names, domain names, a domain expression, unannotated, or has_search_hit is specified,
selection logic is:

    contigs & (cds_ids or domains or domain_expr or unannotated or has_search_hit) & strand & partial

otherwise:
    return every cds in contigs
```

Give at least one of these outputs: `--genbank_out`, `--fasta_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | GenBank files. Use - to read standard input. |
| `-o`, `--genbank_out` | FILE |  | Write the records to this GenBank file. Use - for standard output. |
| `--fasta_out` | FILE |  | Write the records to this FASTA file. Use - for standard output. |
| `--kb_range` | NUMBER |  | How many kb upstream and downstream of the selected CDS to extract. Starting from the ends of the CDS. |
| `--kb_up` | NUMBER |  | How many kb upstream of the selected CDS to extract. Starting from the end of the CDS. |
| `--kb_down` | NUMBER |  | How many kb downstream of the selected CDS to extract. Starting from the end of the CDS. |
| `--cds_range` | INTEGER |  | How many CDSs upstream and downstream of the selected CDS to extract. Example, 1 would mean to return contigs including the selected CDSs and the immediate upstream and downstream CDSs. |
| `--cds_up` | INTEGER |  | How many CDSs upstream of the selected CDS to extract. Example, 1 would mean to return contigs including the selected CDSs and the immediate upstream CDS. |
| `--cds_down` | INTEGER |  | How many CDSs downstream of the selected CDS to extract. Example, 1 would mean to return contigs including the selected CDSs and the immediate downstream CDS. |
| `--whole_contig` |  |  | If set, then return the entire matching contigs. |
| `--contigs` | STRING ... |  | Only extract regions from contigs with these names. |
| `--contigs_file` | FILE |  | text file containing names of contigs to extract from. |
| `--cds` | STRING ... |  | Include CDSs with these names (additive with --domains) |
| `--cds_file` | FILE |  | text file containing names of CDSs to include. |
| `--domains` | STRING ... |  | Include domains with these names |
| `--domains_file` | FILE |  | text file containing names of domains to include. |
| `--domain_expr` | STRING |  | a boolean expression using operators & (AND), \| (OR), and ~(NOT), to specify a desired combination of domains |
| `--databases` | STRING ... |  | Consider only domains from these databases when filtering by domain. default: all databases. |
| `--strand` | {f,r} |  | Only extract regions around CDSs on the specified strand. |
| `--partial` | {include,exclude,only} | include | Which CDSs to extract regions around, by fragment status. include: all. exclude: skip fragments. only: only fragments. A CDS is a fragment if its location has '<' or '>', or it has a /partial qualifier. A protein record is a fragment if it has a UniProt ' (Fragment)' or MGnify 'FL=0' header, a 'Flags: Fragment' definition or NON_TER features (UniProt text entries), or a CDS/Protein feature with '<' or '>' spanning it. Not affected by --invert. |
| `--unannotated` |  |  | Select contigs with no domain annotations. |
| `--annotation_type` | {domainator,domain_search,all} | domainator | Which annotations to select CDSs by. domainator: domain annotations (from domainate). domain_search: best-hit annotations (from domain_search). all: both. |
| `--has_search_hit` |  |  | Also select existing domain_search hits, regardless of annotation_type, domain names, domain expressions, or databases. Hits must pass the e-value and score cutoffs where applicable. Participates in --invert where available. |
| `--pad` |  |  | If set, then ends of sequences will be padded so that the focus CDS aligns. See --anchor. |
| `--anchor` | {start,end,center} | start | Where to align the focus CDSs when --pad is set, relative to the focus CDS in its own orientation. start: align the first base of the focus CDSs. end: align their last base. center: align their midpoints. |
| `--keep_direction` |  |  | by default extracted regions will be flipped so that the focus cds is on the forward strand. Setting this option will keep the focus cds on whatever strand it started on. |
| `--deduplicate`, `--no_deduplicate` |  |  | if the same region is extracted for multiple hits, then keep only one copy of it. |
| `--max_region_overlap` | NUMBER | 1.0 | the maximum fractional of overlap between any two output regions. If >= 1, then no overlap filtering will be done. Regions are output in a greedy fashion based on CDS start site. New regions are output if less than this fraction of them overlaps with any previously output region. |
| `-e`, `--evalue` | NUMBER |  | the evalue cutoff (max) for domain annotations to be considered when selecting CDSs. Default: no cutoff. |
| `--invert` |  |  | Invert the CDS selection criteria. i.e. return CDSs that don't match the CDS selection criteria. (only applies to CDS selection, not neighborhoods or contigs). |
| `--include_nucleic_acids` |  |  | Include nucleic-acid annotations with no associated CDS (those with cds_id='.') when applying CDS selection criteria. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
