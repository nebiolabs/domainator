[CLI reference](README.md)
# domainator find_systems

Annotate systems: domains found together in one neighborhood, in already-annotated records.

```
Annotate systems: sets of domains found together in one genomic neighborhood.

Reads records that already carry domain annotations (from domainate), finds the systems that
the models describe, and adds a Domain_System feature for each one. A system is defined in a
system model file (--system_models), or on the command line with --domain_expr and a window:

    domainator find_systems -i annotated.gb --domain_expr "X & Y & (Q | R)" --window_cds 10 -o systems.gb

Detection uses only the annotations already in the records, so run domainate permissively
enough (a loose -e, no --max_hits_per_cds, --max_overlap 1) that the hits the models need are
there. See docs/file_formats.md for the model format and the Domain_System qualifiers.
```

Give at least one of these outputs: `--genbank_out`, `--fasta_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | GenBank files with domain annotations. Use - to read standard input. |
| `-o`, `--genbank_out` | FILE |  | Write the records to this GenBank file. Use - for standard output. |
| `--fasta_out` | FILE |  | Write the records to this FASTA file. Use - for standard output. |
| `--system_models` | FILE ... |  | System model files (YAML with a domainator_systems: key), or directories holding a systems.yaml. See docs/file_formats.md. |
| `--domain_expr` | STRING |  | A system defined on the command line: a boolean expression over domain (reference) names that must hold within one window, using & (and), \| (or), ~ (not, within the window), parentheses and < (upstream of, in system orientation), and the functions cds(), minimum(), repeatable() and minscore() (see the system model files section of docs/file_formats.md). database:name limits a name to one database. Needs --window_cds or --window_kb. |
| `--window_cds` | INTEGER |  | Window size for --domain_expr, in CDSs (see --window_mode). |
| `--window_kb` | NUMBER |  | Window size for --domain_expr, in kb (see --window_mode). With --window_cds, both must hold. |
| `--window_mode` | {span,chain,anchor} | span | How --window_cds and --window_kb are measured. span: the whole system, first member to last. chain: the gap between consecutive members. anchor: the distance of every member from one anchor member. |
| `--system_name` | STRING | domain_expr | Name of the --domain_expr system (its database is command_line). |
| `--same_strand` |  |  | Members of the --domain_expr system must all be on one strand. |
| `--systems` | STRING ... |  | Only look for the models with these names (name or database:name). Additive with --systems_file; duplicate names are ignored. If neither source is supplied, no name criterion is applied. |
| `--systems_file` | FILE |  | Only look for the models with these names (name or database:name). Read names from a plain-text file (one per line), FASTA/GenBank record IDs, or HMM profile names. Additive with --systems; duplicate names are ignored. Plain-text lines are stripped and blank lines ignored. An empty file adds no names; if no names are supplied directly either, the resulting list is explicitly empty, not an omitted criterion. |
| `-e`, `--evalue` | NUMBER |  | Only hits with e-values at or below this count as evidence, in addition to each model's own thresholds. Default: no cutoff. |
| `--annotation_type` | {domainator,domain_search,all} | domainator | Which annotations count as evidence. domainator: domain annotations (from domainate). domain_search: best-hit annotations (from domain_search). all: both. |
| `--databases` | STRING ... |  | Only hits from these databases count as evidence. Default: all databases. |
| `--min_score` | NUMBER |  | Drop systems with a score below this (the sum of the best bit scores of their counted components). |
| `--min_completeness` | NUMBER |  | Drop systems with a completeness below this (counted components found / model components). Models without a quorum (condition only, like --domain_expr) have no completeness and are kept. |
| `--max_system_overlap` | NUMBER | 1.0 | Drop a system when more than this fraction of its counted members belong to better-scoring systems. 1 keeps every system; 0 makes calls exclusive. |
| `--max_systems_per_contig` | INTEGER |  | Keep at most this many systems per record, best score first. Default: no limit. |
| `--truncated_systems` | {exclude,include} | exclude | include: also report systems that the end of a linear contig may have cut: windows at an end that lack components, but would be systems if the missing ones lay just beyond the end. They have a system_partial qualifier, as do complete systems near an end. Default: exclude. |
| `--clear`, `--no_clear` |  |  | Remove existing Domain_System features from the models' databases first (for --domain_expr, those with the same --system_name). |
| `--partial` | {include,exclude,only} | include | Which CDSs can be system members, by fragment status. include: all. exclude: skip fragments. only: only fragments. A CDS is a fragment if its location has '<' or '>', or it has a /partial qualifier. A protein record is a fragment if it has a UniProt ' (Fragment)' or MGnify 'FL=0' header, a 'Flags: Fragment' definition or NON_TER features (UniProt text entries), or a CDS/Protein feature with '<' or '>' spanning it. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
