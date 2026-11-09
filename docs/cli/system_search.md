[CLI reference](README.md)
# domainator system_search

Search sequence databases for systems and select the neighborhoods around them.

```
Find systems in sequence databases and select the regions around them.

A system is a set of domains found together in one genomic neighborhood (see find_systems).
system_search annotates the contigs of -i/--input with the profiles the system models use (and
any -r/--reference), finds the systems, and returns a region around each one, as domain_search
does for single hits:

    domainator system_search -i genomes.gb --system_models padloc_db --systems gabija --cds_range 5 -o gabija.gb
    domainator system_search -i genomes.gb -r toxins.hmm --domain_expr "TA_toxin & TA_antitoxin" --window_cds 2 -o ta.gb

Systems are defined in model files (--system_models), or on the command line with --domain_expr
and a window, whose names must be references given with -r. Inputs must be nucleotide contigs
with CDSs, or nucleotide sequences with --gene_call.

Padding (--cds_range, --kb_range, ...) is counted from the system's outermost members, in system
orientation, and --model_neighborhood pads each system by its model's own neighborhood.
--min_context_* skips systems without enough context beyond them, such as systems cut by the
end of a linear contig.

Each region has one System_Search feature, marking the system it was selected for, and the
Domainator hits of the system's members. --add_annotations keeps every hit in the region, and
the region's other systems (Domain_System features).

--systems/--systems_file, --system_expr and --system_categories choose which models to look
for. --system_expr is evaluated over the names of the systems found on each contig, and only
systems of the models it names positively become regions: "gabija & ~thoeris_I" returns Gabija
systems on contigs without a Thoeris I system.

Every model profile is searched with the larger of -e and its model's own threshold, and a hit
counts as evidence for a system only if it also passes -e (as in find_systems).

--max_hits applies to the whole search, keeping the best-scoring systems. When --max_hits or
--pad is given, regions are written after the search completes, best first; otherwise they are
written as they are found, in input order.
```

Give at least one of these outputs: `--genbank_out`, `--fasta_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | The GenBank or FASTA files to search: nucleotide contigs with CDSs, or nucleotide sequences with --gene_call. |
| `-r`, `--reference` | FILE ... |  | References to search besides the models' profiles (HMM, FASTA or CM files): the names --domain_expr uses, and other annotations, which are kept with --add_annotations. |
| `--system_models` | FILE ... |  | System model files (YAML with a domainator_systems: key), or directories holding a systems.yaml. See docs/file_formats.md. |
| `--domain_expr` | STRING |  | A system defined on the command line: a boolean expression over domain (reference) names that must hold within one window, using & (and), \| (or), ~ (not, within the window), parentheses and < (upstream of, in system orientation), and the functions cds(), minimum(), repeatable() and minscore() (see the system model files section of docs/file_formats.md). database:name limits a name to one database. Needs --window_cds or --window_kb. |
| `--window_cds` | INTEGER |  | Window size for --domain_expr, in CDSs (see --window_mode). |
| `--window_kb` | NUMBER |  | Window size for --domain_expr, in kb (see --window_mode). With --window_cds, both must hold. |
| `--window_mode` | {span,chain,anchor} | span | How --window_cds and --window_kb are measured. span: the whole system, first member to last. chain: the gap between consecutive members. anchor: the distance of every member from one anchor member. |
| `--system_name` | STRING | domain_expr | Name of the --domain_expr system (its database is command_line). |
| `--same_strand` |  |  | Members of the --domain_expr system must all be on one strand. |
| `-o`, `--genbank_out` | FILE |  | Write the regions around the systems to this GenBank file. Use - for standard output. |
| `--fasta_out` | FILE |  | Write the regions around the systems to this FASTA file. Use - for standard output. |
| `--systems` | STRING ... |  | Select systems with these names (name, or database:name for one collection's model). Additive with --systems_file; duplicate names are ignored. If neither source is supplied, no name criterion is applied. |
| `--systems_file` | FILE |  | Select systems with these names (name, or database:name for one collection's model). Read names from a plain-text file (one per line), FASTA/GenBank record IDs, or HMM profile names. Additive with --systems; duplicate names are ignored. Plain-text lines are stripped and blank lines ignored. An empty file adds no names; if no names are supplied directly either, the resulting list is explicitly empty, not an omitted criterion. |
| `--system_expr` | STRING |  | A boolean expression over the names of the systems in a record, using & (and), \| (or), ~ (not) and parentheses, e.g. "Gabija & ~Thoeris". A name is name or database:name. |
| `--system_categories` | STRING ... |  | Select systems in these categories (the category qualifier, e.g. defense or PKS). With --systems, a system must match both. |
| `-Z` | INTEGER | 1000 | The database size for e-values (hmmsearch's -Z). 0: the number of CDSs in the inputs, which takes an extra pass to count unless they have .didx indexes. |
| `-e`, `--evalue` | NUMBER | 0.001 | Count only hits with an e-value at most this as evidence. Model profiles are searched with this or their model's own threshold, whichever is larger. |
| `--min_evalue` | NUMBER | 0.0 | Hits with an e-value below this are not reported. |
| `--max_mode` |  |  | Turn off hmmsearch's heuristic filters: slower, but finds weak hits. |
| `--min_score` | NUMBER |  | Only select systems with at least this score. |
| `--min_completeness` | NUMBER |  | Only select systems with at least this completeness (0 to 1). Systems of models without a quorum have no completeness, and are not filtered by it. |
| `--max_system_overlap` | NUMBER | 1.0 | Drop a system when more than this fraction of its counted members belong to better-scoring systems. 1 keeps every system; 0 makes calls exclusive. |
| `--max_systems_per_contig` | INTEGER |  | Keep at most this many systems per record, best score first. Default: no limit. |
| `--truncated_systems` | {exclude,include} | exclude | include: also report systems that the end of a linear contig may have cut: windows at an end that lack components, but would be systems if the missing ones lay just beyond the end. They have a system_partial qualifier, as do complete systems near an end. Default: exclude. |
| `--max_hits` | INTEGER |  | The most systems to return, best score first. Default: all. |
| `--include_taxids` | INTEGER ... |  | Include these NCBI taxids and their descendants (OR). |
| `--exclude_taxids` | INTEGER ... |  | Exclude these NCBI taxids and their descendants; exclusions take precedence. |
| `--taxonomy_expr` | STRING |  | Boolean expression over taxid lineages using &, \|, ~, and parentheses, e.g. "2 & ~1224". Mutually exclusive with include/exclude_taxids. |
| `--taxonomy_unknown` | {exclude,include,evaluate,error} | exclude | Unresolved taxonomy: exclude (default), include, evaluate against an empty lineage, or error. With ~2, exclude rejects unlabeled candidates; evaluate and include allow them. Only active with taxonomy criteria. Not inverted. |
| `--ncbi_taxonomy_path` | DIR | /tmp/ncbi_taxonomy_sean | NCBI taxonomy database directory; downloaded if missing and filtering is active. |
| `--taxonomy_update` |  |  | Check for an updated taxonomy database when filtering is active. |
| `--partial` | {include,exclude,only} | include | Which CDSs to search and count as system members, by fragment status. include: all. exclude: skip fragments. only: only fragments. A CDS is a fragment if its location has '<' or '>', or it has a /partial qualifier. A protein record is a fragment if it has a UniProt ' (Fragment)' or MGnify 'FL=0' header, a 'Flags: Fragment' definition or NON_TER features (UniProt text entries), or a CDS/Protein feature with '<' or '>' spanning it. |
| `--cpu` | INTEGER | 0 | The number of processes to use (0: one per physical core; at least 2). |
| `--add_annotations` |  |  | Keep every hit in each region, and its other systems. Default: keep only the hits of the system's members. |
| `--dedup` | {none,id,sequence,sequence_rc} | none | Output uniqueness policy: none keeps all records; id drops repeated output IDs; sequence drops identical output sequences (ignoring case and terminal '*'); sequence_rc also collapses nucleotide reverse complements. The same ID must identify the same sequence, so sequence uniqueness also ensures unique IDs. Deduplication happens before padding. |
| `--max_region_overlap` | NUMBER | 1.0 | The largest fraction of a region that may overlap a region already output from the same contig. If >= 1, then no overlap filtering will be done. |
| `--translate` |  |  | Output one protein per member CDS of each system, instead of the regions. Each protein's CDS feature has a system_membership qualifier. |
| `--kb_range` | NUMBER |  | How many kb upstream and downstream of the system to extract, from its outermost members. |
| `--kb_up` | NUMBER |  | How many kb upstream of the system to extract, in system orientation. |
| `--kb_down` | NUMBER |  | How many kb downstream of the system to extract, in system orientation. |
| `--cds_range` | INTEGER |  | How many CDSs upstream and downstream of the system's outermost members to extract. Only protein-coding CDSs are counted, never hits on a nucleotide reference. |
| `--cds_up` | INTEGER |  | How many CDSs upstream of the system to extract, in system orientation. |
| `--cds_down` | INTEGER |  | How many CDSs downstream of the system to extract, in system orientation. |
| `--whole_contig` |  |  | Extract the whole contig around each system. |
| `--model_neighborhood` |  |  | Pad each system by its model's neighborhood (the neighborhood_kb qualifier, else neighborhood_cds), instead of by --kb_* or --cds_*. Systems whose model has no neighborhood are not padded. |
| `--min_context_cds` | INTEGER |  | Only select systems with at least this many CDSs upstream and downstream of them on the contig, in system orientation, for example to skip systems cut by a contig end. On a circular contig, the CDSs outside the system must add up to at least twice this. Independent of the padding options. Pseudo-CDSs for nucleotide hits are not counted. |
| `--min_context_cds_up` | INTEGER |  | Only select systems with at least this many CDSs upstream of them. Can't be combined with --min_context_cds. |
| `--min_context_cds_down` | INTEGER |  | Only select systems with at least this many CDSs downstream of them. Can't be combined with --min_context_cds. |
| `--min_context_kb` | NUMBER |  | Only select systems with at least this many kb upstream and downstream of them on the contig, in system orientation, for example to skip systems cut by a contig end. On a circular contig, the kb outside the system must add up to at least twice this. Independent of the padding options. |
| `--min_context_kb_up` | NUMBER |  | Only select systems with at least this many kb upstream of them. Can't be combined with --min_context_kb. |
| `--min_context_kb_down` | NUMBER |  | Only select systems with at least this many kb downstream of them. Can't be combined with --min_context_kb. |
| `--strand` | {f,r} |  | Only select systems in this orientation. |
| `--pad` |  |  | Pad the ends of the regions so that their systems align. See --anchor. |
| `--anchor` | {start,end,center} | start | Where to align the systems when --pad is set, in system orientation. start: their first base. end: their last base. center: their midpoints. |
| `--keep_direction` |  |  | By default, regions are flipped so that their system is on the forward strand. Set this to keep each region in the direction of its contig. |
| `--batch_size` | INTEGER | 10000 | About how many CDSs each process searches at once. |
| `--gene_call` | {all,unannotated} |  | Call genes with Prodigal (metagenomic mode) before searching. all: replace every CDS. unannotated: only contigs without CDSs. |
| `--max_output_gb` | NUMBER | 25.0 | Fail if a projected output file would exceed this many GB. Set to 0 to disable the guardrail. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
