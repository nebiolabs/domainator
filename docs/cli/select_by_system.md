[CLI reference](README.md)
# domainator select_by_system

Extract neighborhoods around systems already annotated by find_systems.

```
Extract contig regions around systems already annotated in GenBank records

Takes GenBank records with systems annotated by find_systems (Domain_System features) and
extracts a region around each selected system, the way select_by_cds extracts regions around
CDSs. Padding (--cds_range, --kb_range, ...) is counted from the system's outermost members, in
system orientation, and --model_neighborhood pads each system by its model's own neighborhood.

Each region gets a System_Search feature marking the system it was selected for, so it can be
handled like system_search output.

--systems/--systems_file and --system_categories select systems (both must match when both
are given). --invert negates that selection only. The other criteria (--databases, --min_score,
--min_completeness, --partial_systems, --strand and --min_context_*) restrict which systems can
be selected, and are not inverted.

--min_context_* skips systems without enough context beyond them, for example systems cut by
the end of a linear contig: --cds_range 10 --min_context_cds 10 selects only systems with 10
CDSs on each side.
```

Give at least one of these outputs: `--genbank_out`, `--fasta_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE ... | required | GenBank files with system annotations (from find_systems). Use - to read standard input. |
| `-o`, `--genbank_out` | FILE |  | Write the records to this GenBank file. Use - for standard output. |
| `--fasta_out` | FILE |  | Write the records to this FASTA file. Use - for standard output. |
| `--systems` | STRING ... |  | Select systems with these names (name, or database:name for one collection's model). Additive with --systems_file; duplicate names are ignored. If neither source is supplied, no name criterion is applied. |
| `--systems_file` | FILE |  | Select systems with these names (name, or database:name for one collection's model). Read names from a plain-text file (one per line), FASTA/GenBank record IDs, or HMM profile names. Additive with --systems; duplicate names are ignored. Plain-text lines are stripped and blank lines ignored. An empty file adds no names; if no names are supplied directly either, the resulting list is explicitly empty, not an omitted criterion. |
| `--system_categories` | STRING ... |  | Select systems in these categories (the category qualifier, e.g. defense or PKS). With --systems, a system must match both. |
| `--databases` | STRING ... |  | Only select systems from these model collections (the database qualifier). Default: all. |
| `--system_type` | {domain_system,system_search,all} | domain_system | Which systems to select. domain_system: systems found by find_systems (Domain_System features). system_search: systems selected by a search (System_Search features). all: both. |
| `--min_score` | NUMBER |  | Only select systems with at least this score. |
| `--min_completeness` | NUMBER |  | Only select systems with at least this completeness (0 to 1). Systems of models without a quorum have no completeness, and are not filtered by it. |
| `--partial_systems` | {include,exclude,only} | include | Which systems to select by their system_partial flag, set when a system touches a linear contig end and a missing component could lie beyond it. include: all. exclude: skip flagged systems. only: only flagged systems. |
| `--invert` |  |  | Invert --systems/--systems_file and --system_categories, not the other criteria. With neither, warn and ignore. |
| `--strand` | {f,r} |  | Only select systems in this orientation. |
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
| `--pad` |  |  | Pad the ends of the regions so that their systems align. See --anchor. |
| `--anchor` | {start,end,center} | start | Where to align the systems when --pad is set, in system orientation. start: their first base. end: their last base. center: their midpoints. |
| `--keep_direction` |  |  | By default, regions are flipped so that their system is on the forward strand. Set this to keep each region in the direction of its contig. |
| `--dedup` | {none,id,sequence,sequence_rc} | id | Output uniqueness policy: none keeps all records; id drops repeated output IDs; sequence drops identical output sequences (ignoring case and terminal '*'); sequence_rc also collapses nucleotide reverse complements. The same ID must identify the same sequence, so sequence uniqueness also ensures unique IDs. Deduplication happens before padding. |
| `--max_region_overlap` | NUMBER | 1.0 | The largest fraction of a region that may overlap a region already output from the same contig. If >= 1, then no overlap filtering will be done. |
| `--contigs` | STRING ... |  | Only read input records with these IDs; not inverted by --invert. Additive with --contigs_file; duplicate names are ignored. If neither source is supplied, no name criterion is applied. |
| `--contigs_file` | FILE |  | Only read input records with these IDs; not inverted by --invert. Read names from a plain-text file (one per line), FASTA/GenBank record IDs, or HMM profile names. Additive with --contigs; duplicate names are ignored. Plain-text lines are stripped and blank lines ignored. An empty file adds no names; if no names are supplied directly either, the resulting list is explicitly empty, not an omitted criterion. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
