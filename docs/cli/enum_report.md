[CLI reference](README.md)
# domainator enum_report

Tabulate properties of each contig, CDS, or domain.

```
Writes record-wise information from Genbank files into tab-separated or html tables

Takes a genbank file which has been annotated using Domainate, write tabulated descriptions of each contig, CDS, domain
or system.

With --by system, each row is the region of one system (a Domain_System feature, from find_systems), so columns
such as domains, architecture and length describe the system. The systems columns report, for each row, the
systems on the contig (--by contig), the systems a CDS is a member of (--by cds and --by domain), or the system
itself (--by system).
```

Give at least one of these outputs: `--tsv_out`, `--html_out`, `--json_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | Genbank filenames. Use - to read standard input. |
| `-o`, `--tsv_out` | FILE |  | Write the table to this tab-separated file. Use - for standard output. |
| `--html_out` | FILE |  | Write the table in html format to this file. Use - for standard output. |
| `--json_out` | FILE |  | Write streaming NDJSON (one compact JSON object per record) to this file. Use - for standard output. |
| `--html_max_height` | INTEGER |  | Max height in pixels to set the html output to. |
| `--contigs` | STRING ... |  | only examine contigs with ids in this list. Additive with --contigs_file. default: all contigs. |
| `--contigs_file` | FILE |  | only examine contigs with ids listed in this file (one per line). Additive with --contigs. default: all contigs. |
| `--by` | {contig,cds,domain,system} | contig | One line in output for every contig, cds, domain or system (the region of each Domain_System feature, from find_systems). For protein genbanks, contig and cds are treated the same. |
| `--partial` | {include,exclude,only} | include | Which rows to report, by fragment status. include: all. exclude: skip fragments. only: only fragments. A CDS is a fragment if its location has '<' or '>', or it has a /partial qualifier. A protein record is a fragment if it has a UniProt ' (Fragment)' or MGnify 'FL=0' header, a 'Flags: Fragment' definition or NON_TER features (UniProt text entries), or a CDS/Protein feature with '<' or '>' spanning it. By contig, a nucleotide contig is a fragment if any of its CDSs is. By domain, the CDS or protein the domain is on decides. |
| `--databases` | STRING ... |  | Ignore domain annotations not from these databases. default: consider all databases. |
| `--columns` | columns... |  | The columns to report, in order. A column with arguments takes them after ':', as in 'qualifier:CDS:product'. Can be given more than once. |
| `--taxonomy_update` |  |  | If taxonomy database exists, check it against the version on the ncbi server and update if there is a newer version. |
| `--ncbi_taxonomy_path` | DIR | /tmp/ncbi_taxonomy_sean | Path to NCBI taxonomy database directory. Will be created and downloaded if it does not exist. |
| `--column_names` | STRING ... |  | If supplied, then this list will be used instead of the default column names. |

## Columns

| Column | Description |
|---|---|
| `domains` | the domain names, in alphabetical order and separated by '; ', each listed once. |
| `domain_descriptions` | the domain descriptions, in alphabetical order of domain name and separated by '; ', each listed once. |
| `score` | the sum of the domain scores. |
| `identity` | the sum of the domain identities. |
| `architecture` | the domain names in order of start position, separated by '; '. |
| `architecture_detailed` | the domain architecture in order of start position. Within each CDS, domains are separated by '/'; CDSs are separated by ';'. CDSs with no domains are 'no_domains'. |
| `start` | the start coordinate of the feature. |
| `end` | the end coordinate of the feature. |
| `strand` | the strand of the feature. |
| `length` | the length of the contig, CDS or domain. |
| `topology` | the topology (usually only meaningful with --by contig). |
| `sequence` | the sequence of the contig, CDS or domain. |
| `cds_count` | the number of CDSs. |
| `domain_count` | the total number of domains. |
| `definition` | the record's definition line. |
| `filename` | the name of the input file the record came from. |
| `domain_search` | several columns describing the record's best domain_search hit: best_query, query_description, query_start, query_end, query_length, match_start, match_end, match_strand, match_score, match_identity. |
| `taxid_lineage` | the taxids from the root to the record's taxid, separated by ';'. |
| `taxname_lineage` | the taxon names from the root to the record's taxid, separated by ';'. |
| `rank_lineage` | the ranks from the root to the record's taxid, separated by ';'. |
| `taxid:rank[:...]` | the taxid at each given rank (e.g. domain, genus, species; 'self' for the record's own taxid), in columns named taxid_\[rank\]. |
| `taxname:rank[:...]` | the taxon name at each given rank (e.g. domain, genus, species; 'self' for the record's own taxon), in columns named taxname_\[rank\]. |
| `qualifier:feature:qualifier` | the values of this qualifier on features of this type, separated by ';'. |
| `feature_count:feature` | the number of features of this type. |
| `motif_count:motif` | the number of matches of this PROSITE motif in the sequence (nucleotide motifs for nucleotide sequences, protein motifs for proteins; see https://emboss.sourceforge.net/apps/release/6.5/emboss/apps/fuzzpro.html). |
| `isoelectric_point` | the isoelectric point of the protein. Only defined for contigs with exactly one protein (a nucleotide contig with one CDS, or a protein record). |
| `net_charge:pH` | the net charge of the protein at this pH. Only defined for contigs with exactly one protein. |
| `repeat_count:length` | the number of repeats of at least this length; for example AGGKKK has repeat_count_1 = 3, repeat_count_2 = 2 and repeat_count_3 = 1. |
| `named_domain_count:domain` | the number of Domainator domains with this name, in a column named named_domain_count_\[domain\]. |
| `starts_with` | the first character of the sequence. |
| `molecular_weight` | the molecular weight of the sequence in Daltons, for protein, DNA or RNA according to its molecule type. |
| `partial_count` | the number of partial CDSs: by contig, the CDSs with '<' or '>' in their location (for a protein record, 1 if it is a fragment); by cds, 1 if the CDS or protein is a fragment, otherwise 0; by domain, 1 if its CDS or protein is a fragment, whether or not the domain reaches the missing end. |
| `fragment` | whether the protein or CDS is a fragment: complete, N (missing the N-terminus), C (missing the C-terminus), NC (missing both), or ? (missing ends unknown, e.g. from a UniProt ' (Fragment)' or MGnify 'FL=0' header). Only for rows that are a single protein or CDS; otherwise empty, with a warning. |
| `append:name:type:value` | a column with this name and type (str, int or float) holding this value in every row. |
| `systems` | the names of the row's systems (see the description above), in alphabetical order and separated by '; ', each listed once. |
| `system_categories` | the categories of the row's systems, in alphabetical order and separated by '; ', each listed once. |
| `system_components` | role:component for each member of the row's systems, in contig order, separated by '; '. By cds and domain, the components the CDS fills. |
| `system_count:name` | the number of the row's systems with this name (name or database:name), in a column named system_count_\[name\]. |
| `system_search` | several columns describing the record's System_Search feature (from select_by_system): system_search_name, system_search_database, system_search_score, system_search_evalue, system_search_completeness, system_search_core_count, system_search_total_count. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
