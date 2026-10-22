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

Criteria are combined with the following logic:
    (all OR target_CDSs OR target_domains OR unannotated OR domain_search hits) AND (target_contigs OR strand) AND partial
```

Give at least one of these outputs: `--genbank_out`, `--fasta_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | Genbank filenames. Use - to read standard input. |
| `-o`, `--genbank_out` | FILE |  | Write the records to this GenBank file. Use - for standard output. |
| `--fasta_out` | FILE |  | Write the records to this FASTA file. Use - for standard output. |
| `--contigs` | STRING ... |  | only look for peptides on contigs with ids in this list. Additive with --contigs_file. default: all contigs. |
| `--contigs_file` | FILE |  | only look for peptides on contigs with ids listed in this file (one per line). Additive with --contigs. default: all contigs. |
| `--domains` | STRING ... |  | Only extract peptides containing these domains |
| `--domains_file` | FILE |  | text file containing names of domains, where peptides containing these domains will be extracted. |
| `--databases` | STRING ... |  | Consider only domains from these databases when filtering by domain. default: all databases. |
| `--unannotated` |  |  | extract peptides from CDSs with no domain annotations. |
| `--annotation_type` | {domainator,domain_search} | domainator | Which annotations to select peptides by. With domain_search, the search hits are extracted. domainator: domain annotations (from domainate). domain_search: best-hit annotations (from domain_search). |
| `--strand` | {f,r} |  | Only extract peptides from CDSs on the specified strand. |
| `--partial` | {include,exclude,only} | include | Which CDSs to extract peptides from, by fragment status. include: all. exclude: skip fragments. only: only fragments. A CDS is a fragment if its location has '<' or '>', or it has a /partial qualifier. A protein record is a fragment if it has a UniProt ' (Fragment)' or MGnify 'FL=0' header, a 'Flags: Fragment' definition or NON_TER features (UniProt text entries), or a CDS/Protein feature with '<' or '>' spanning it. Not affected by --invert. |
| `--invert` |  |  | Invert the CDS selection criteria. i.e. return peptides from CDSs that don't match the CDS selection criteria. (Only applies to CDS selection, not contig selection) |
| `--cds` | STRING ... |  | only look for peptides on cdss with these names |
| `--cds_file` | FILE |  | only look for peptides on cdss with names listed in this text file. |
| `--keep_name` |  |  | by default, extracted regions will be renamed to note their original coordinates and whether they were reverse complemented. Set this option to keep contigs with their original names. WARNING: if you aren't careful, this can easily result in duplicate names. |
| `--name_field` | STRING |  | The CDS qualifier to preferentially use as the name of the peptide. If not supplied, then qualifiers will be checked in order of: ["protein_id", "gene_id", "locus_tag", "gene"] |
| `-e`, `--evalue` | NUMBER |  | the evalue cutoff (max) for domain annotations to be put onto the new genbank file. Default: no cutoff. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
