# Systems: commands and tool integration

## Three routes

Detection is a pure function of an annotated contig ([detection.md](detection.md)), so the
real question is *where users call it*.

### Route A: one new post-hoc editor

Add a single `find_systems` editor. Users run `domainate` with the right profiles, then
`find_systems --system_models ...`. Every other tool learns to *read* `Domain_System`
features.

- **For:** the smallest change. No search tool changes. It streams (`-i -`/`-o -`), chains
  like every editor, and lets users try several definitions against one expensive
  annotation run.
- **Against:** always two steps. Users must know to run domainate permissively (see
  [definition_formats.md](definition_formats.md#thresholds-are-applied-at-detection-time)).
  And it gives no way to *search* a genome database for systems and cut out their
  neighborhoods, which is the domain_search use case in the request.

### Route B: build systems into the existing search commands

Add `--system_models` to `domainate` (annotate) and to `domain_search` (select). In
`domain_search`, the flag switches the focus unit from "best hit in a CDS" to "a system":

- `--cds_range`, `--kb_range`, `_up` and `_down` are measured from the CDSs at either end of
  the system;
- `--max_hits` counts systems;
- the output gets `System_Search` features.

- **For:** one step, and it matches how users already think ("domain_search, but for
  systems"). No new command names.
- **Against:** `domain_search` is already 40 options. In system mode, `-r` becomes optional,
  `--decoys` and `--max_overlap` mean something different or nothing, the output feature
  type changes, and so does the record-annotation protocol between the worker and
  `select_by_cds` (`hit_scores`/`hit_names` are keyed by single CDS). Every option's help
  then needs a "(in system mode: ...)" clause, and `domainator describe` can't express the
  mode-dependent semantics for GUIs.

### Route C: a parallel family of system commands

Add `system_search` and `select_by_system`, as `structure_search` and `structure_domainate`
sit beside `domain_search` and `domainate`. They share helpers, not modules.

- **For:** each command has one meaning per option, help text stays focused, and the
  manifest stays honest.
- **Against:** more command names, and option lists that must be kept in step (solved with
  shared `cli/common.py` helpers, as for every other shared option).

### Chosen: Route A + `domainate --system_models` + Route C for search

| Need | Command |
|---|---|
| Detect systems on records that already have domain hits | **`find_systems`** (new editor) |
| Search profiles *and* detect systems in one pass | **`domainate --system_models`** (new options) |
| Find systems in a database and cut out their neighborhoods | **`system_search`** (new; the domain_search analog) |
| Cut out neighborhoods around systems already annotated | **`select_by_system`** (new editor; the select_by_cds analog) |
| Download or convert DefenseFinder, PADLOC, MacSyFinder and antiSMASH models into a native directory | **`system_db_download`** (new; the only command that reads foreign formats) |
| Inspect model files; write a command-line model as a file; export the profiles a set of models needs; check models against their examples | **`format_systems`** (new; optional) |
| Turn system calls made by antiSMASH (later DefenseFinder, PADLOC) into `Domain_System` features | **`import_systems`** (new; see [antismash.md](antismash.md#importing-antismash-results)) |

`domainate` gets flags rather than a sibling command because annotation is additive: the
output is still "every record, with annotations". Nothing about domainate's existing options
changes meaning. Selection is where the meanings fork, so it gets its own commands.

### Route B vs C for search

**Decided: Route C.** System search is the new `system_search` command, and
`domain_search` is unchanged apart from the shared neighborhood refactor. The deciding
points were:

1. a mode switch would change `--max_hits` and the output feature type inside a 1.x command;
2. the manifest would have to describe one command with conditional semantics.

The 1.0 plan's answer to similar cases was consistent names, each with one meaning (Phase 3,
"Long flags reused with different meanings").

## Shared CLI vocabulary

New helpers in `domainator/cli/common.py`, so every command spells these the same way:

| Helper | Options | Notes |
|---|---|---|
| `add_system_models(parser, required=...)` | `--system_models PATH [PATH ...]` | Native YAML files, or directories written by `system_db_download`. A foreign format gives an `InputError` with the conversion command. The manifest classifies `system_models` as an input path (add it to `INPUT_PATHS`). |
| `add_window_model(parser)` | `--domain_expr`, `--window_cds`, `--window_kb`, `--window_mode`, `--system_name`, `--same_strand` | Defines one short-form model on the command line. See [Command-line models](#command-line-models). Add `system_name` to `TEXT_OPTIONS` (`domain_expr` is already there). |
| `add_system_selection(parser, expr=True)` | `--systems`, `--systems_file` (via `add_list_and_file`), `--system_expr`, `--system_categories` | Mirrors `--domains/--domains_file/--domain_expr`. Names match the `name` qualifier; `database:name` limits a name to one collection. `system_expr` uses `BooleanEvaluator`. `--system_categories` (e.g. `PKS RiPP`, `defense`) matches the `category` qualifier and is ANDed with the name criteria. Add `systems`, `system_expr` and `system_categories` to `TEXT_OPTIONS`. |
| `add_system_type(parser)` | `--system_type {domain_system,system_search,all}` | Mirrors `--annotation_type`: `model.SYSTEM_TYPES`, with `system_feature_types()` beside `annotation_feature_types()`. |
| `add_system_detection(parser)` | `--max_system_overlap` (1.0), `--max_systems_per_contig` (None) | For the three commands that detect: `find_systems`, `domainate`, `system_search`. `--partial_systems {exclude,include}` joins it in Phase 4, when `include` is implemented. |
| `add_neighborhood(parser)` | `--kb_range/--kb_up/--kb_down`, `--cds_range/--cds_up/--cds_down`, `--whole_contig`, `--strand`, `--pad`, `--anchor`, `--keep_direction`, `--max_region_overlap` | A refactor that extracts the existing domain_search/select_by_cds options. No CLI change for those tools; manifest snapshot unchanged. |
| `add_model_neighborhood(parser)` | `--model_neighborhood` | System selectors only. Pad each system by its model's `neighborhood_kb`/`neighborhood_cds` (antiSMASH's protocluster extent). It can't be combined with explicit `--kb_*`/`--cds_*` options. Systems whose model has no neighborhood get none. |
| `add_min_context(parser)` | `--min_context_cds`, `--min_context_cds_up`, `--min_context_cds_down`, `--min_context_kb`, `--min_context_kb_up`, `--min_context_kb_down` | System selectors only. See [Requiring context](#requiring-context). |

### Requiring context

A system near the end of a linear contig is selected with whatever context is there, so
`--cds_range 10` can give a region with 2 CDSs on one side. For domain hits the existing
recipe handles this: a hit is always one CDS, so `domain_search --cds_range 10` followed by
`select_by_contig --min_cds 21` keeps only hits with 10 CDSs on each side. Systems vary in
size, so a region-size filter can't do the same, and `/system_partial` only flags systems
that may be missing a component. `system_search` and `select_by_system` therefore take
context requirements:

- `--min_context_cds N` is short for `--min_context_cds_up N --min_context_cds_down N`, and
  `--min_context_kb X` for the `_up`/`_down` pair. Each `_up`/`_down` option requires at least
  that much context beyond the system's upstream or downstream end, in system orientation,
  so one side alone can be required: `--min_context_cds_down 5`.
- The both-sides form and its `_up`/`_down` forms are mutually exclusive, as `--cds_range` is
  with `--cds_up`/`--cds_down` (shared code with `select_by_cds.validate_range`). CDS and kb
  requirements may be combined, and both must hold.
- CDSs are counted as `--cds_range` counts them: pseudo-CDSs for ncRNA hits never count. So
  `--cds_range N --min_context_cds N` is the system version of the domain_search recipe, and
  `--cds_up 3 --min_context_cds_up 3` does the same for one side. Kb are measured from the
  system's location to the contig end.
- On a circular contig, the test passes when the CDSs (or kb) outside the system add up to
  at least up + down.
- The requirement is independent of the padding options. Systems that fail it are not
  selected, and don't count toward `--max_hits`.
- Library parameters: `min_context_cds: Optional[Tuple[int, int]]` and
  `min_context_kb: Optional[Tuple[float, float]]`, each (up, down), like `cds_range` and
  `kb_range`.

Score and completeness filters: in system-only commands, `--min_score` and
`--min_completeness` refer to the system. In commands where `--min_score` already means a
domain score, the system filters are `--min_system_score` and `--min_system_completeness`.
This follows the 1.0 rule that a long flag has one meaning per command.

## Command-line models

The common case, "windows of N CDSs containing these domains in this combination", needs
no model file. `find_systems`, `system_search` and `format_systems` take one model from
these options:

| Option | Meaning |
|---|---|
| `--domain_expr EXPR` | A boolean expression over domain (reference) names that must hold **within one window**. `db:name` limits a name to one database. `<` gives order, e.g. `"HicB < RelE"` ([Order](definition_formats.md#order)). |
| `--window_cds N`, `--window_kb N` | The window size. At least one is required with `--domain_expr`. |
| `--window_mode {span,chain,anchor}` | Default `span`, the literal reading of "within N CDSs". |
| `--system_name NAME` | Default `domain_expr`. It becomes the `name` qualifier; `database` is `command_line`, and `description` holds the expression. |
| `--same_strand` | All members on one strand. |

They build exactly the short-form model of
[definition_formats.md](definition_formats.md#the-short-form): every name is its own
component, there is no quorum, a CDS carrying both X and Y satisfies both, and `~W` means
"no W in this window". The model is added to any `--system_models`, so a command needs at
least one of the two.

```bash
# search and select in one step
domainator system_search -i genomes/*.gb -r xyzqr.hmm -e 1e-20 \
    --domain_expr "X & Y & Z & (Q | R)" --window_cds 10 --cds_range 5 -o hits.gb

# or annotate once, then try windows cheaply
domainator domainate -i genomes.gb -r xyzqr.hmm -o annotated.gb
domainator find_systems -i annotated.gb -e 1e-20 --domain_expr "X & Y & Z & (Q | R)" --window_cds 10 -o systems.gb
domainator select_by_system -i systems.gb --cds_range 5 -o regions.gb

# keep the model as a file once it grows
domainator format_systems --domain_expr "X & Y & Z & (Q | R)" --window_cds 10 --system_name XYZ_QR --yaml_out xyz.yaml
```

**Why `--domain_expr` (decided).** The option already has one meaning across commands: "a
boolean over domain names, evaluated over the unit this command works on". The unit is the
CDS in `select_by_cds`, `extract_peptides` and `trim_contigs`, and the contig in
`select_by_contig`. In the system commands it is the window. `--system_expr` is taken: it
selects systems by system name. A separate `--window_expr` was considered and rejected. In
the existing commands `--domain_expr` keeps its plain boolean grammar; `<` is accepted only
where the unit is a window.

**Thresholds.** In `system_search`, `-e` is the search e-value, as in `domain_search`. In
`find_systems`, `-e` is an evidence filter: a hit counts only if it passes both `-e` and the
model's own threshold.

## New commands

### `find_systems` (category `annotate`, editor)

> Annotate systems of co-located domains on records that already carry domain annotations.

| Option | Meaning |
|---|---|
| `-i/--input` | GenBank, `-` allowed (single pass) |
| `-o/--genbank_out`, `--fasta_out` | outputs (`add_seq_outputs`) |
| `--system_models` | model files or converted directories |
| `add_window_model` options | a command-line model ([above](#command-line-models)). At least one of `--system_models` and `--domain_expr` is required. |
| `--systems/--systems_file` | evaluate only these models (no expression: this restricts the models, it doesn't select records) |
| `-e/--evalue` | evidence filter: only hits with e-values at or below this count, in addition to each model's own thresholds |
| `--annotation_type` | which hits count as evidence; default `domainator` |
| `--databases` | only hits from these databases count as evidence |
| `--min_score`, `--min_completeness` | drop system hits below these |
| `add_system_detection` options | overlap, cap, partial systems |
| `--clear/--no_clear` | replace existing `Domain_System` features from the same model databases (default on). For a command-line model (database `command_line`), only systems with the same name are replaced. |
| `--partial` | the usual CDS-fragment filter, applied to members |

Library: `find_systems(records, *, system_models=None, domain_expr=None, window_cds=None,
window_kb=None, window_mode="span", system_name="domain_expr", same_strand=False,
systems=None, evalue=None, annotation_type="domainator",
databases=None, min_score=None, min_completeness=None, max_system_overlap=1.0,
max_systems_per_contig=None, clear=True, partial="include")
-> Iterator[SeqRecord]`. `system_models` also accepts already-read `ModelCollection`s (as
`domainate` accepts `pre_parsed_references`). `partial_systems` is added in Phase 4.

### `domainate --system_models`

New options: `add_system_models`, `--systems/--systems_file` (restricts models),
`add_system_detection`.

- domainate already requires "at least one of `-r` or `--foldseek`" (a `parser.error` in
  `main`). That check widens to accept `--system_models` when the models supply references. The profiles named
  by the models' `reference_files` are loaded through the existing `read_references`, keyed
  by their own database names. `-r` files are searched too, and users can combine them, for
  example adding Pfam for context.
- Profiles from models are searched with the per-profile e-value and are exempt from
  `--max_hits_per_cds`/`--max_overlap`, as described in
  [definition_formats.md](definition_formats.md#thresholds-are-applied-at-detection-time).
- Detection runs in the same worker, right after `domainator_inner` annotates a contig. It
  adds no pass, and `-i -` support is unchanged.
- Protein records (single proteins) get no systems, and their hits are written as usual.
  That is correct, not an error.

### `system_search` (category `annotate`, search)

> Find systems in sequence databases and select the regions around them.

It mirrors `domain_search`, sharing its worker structure (`_domainate(hits_only=...)`, then
detection, then neighborhood selection) and its option helpers.

Unlike domain_search, which annotates with `no_annotations=not add_annotations` and
`best_annotation=True`, the worker must annotate every hit (`no_annotations=False`, no
`Domain_Search` best hits), with model profiles exempt from the overlap and per-CDS caps,
because detection needs all the evidence. Hits that are not members of the selected system
are removed afterwards, unless `--add_annotations` is given.

| Option | In system_search |
|---|---|
| `-i`, `--alphabet`, `-o`, `--fasta_out`, `-Z`, `-e`, `--min_evalue`, `--max_mode`, `--cpu`, `--batch_size`, `--gene_call`, `--partial`, `--max_output_gb`, taxonomy filters, `--dedup` | as in domain_search |
| `--system_models` | model files or converted directories |
| `add_window_model` options | a command-line model ([above](#command-line-models)). At least one of `--system_models` and `--domain_expr` is required. The expression's names must be references given with `-r`, or profiles of the `--system_models`. |
| `-r/--reference` | references for a command-line model, and extra context annotations |
| `--systems/--systems_file/--system_expr/--system_categories` | which models to search for. The expression is evaluated over the names of systems found on the contig, so `"Gabija & ~Thoeris_I"` keeps Gabija systems on contigs without Thoeris. Only systems of models named positively in the expression become focus regions (here Gabija, not Thoeris_I). `--system_categories PKS` restricts an antiSMASH rule set to polyketide rules. |
| `--model_neighborhood` | pad by each model's own neighborhood instead of `--kb_*`/`--cds_*`; this reproduces antiSMASH protocluster extents |
| `add_min_context` options | select only systems with at least this much context on each side, or on one side ([Requiring context](#requiring-context)) |
| `add_system_detection` | as in find_systems |
| `--max_hits` | caps the number of **systems** returned, best score first (`dedup.TopHits` over `get_hit_best_score`) |
| `--max_hits_per_contig` | systems per contig |
| neighborhood options | **measured from the system's end CDSs**: `--cds_up 3` adds 3 CDSs upstream of the first member, in system orientation; `--kb_range 5` adds 5 kb beyond each end. `--anchor {start,end,center}` places a fixed-size window relative to the envelope. `--whole_contig` returns the whole contig, centered on the system when it is circular. `--strand` filters by system orientation. `--keep_direction` turns off reorienting the region to the system's strand. |
| `--max_region_overlap` | as in domain_search, between system regions |
| `--add_annotations` | without it, the output carries the `System_Search` feature and the member `Domainator` hits (so the system is visible). With it, every hit in the region and every other `Domain_System` feature is written too. |
| `--translate` | output one protein record per member CDS, carrying the member's hits and the `system_membership` summary qualifier ([data_model.md](data_model.md#optional-per-cds-summary)) |
| `--decoys` | not offered. Decoys are about which reference a CDS hits best, which doesn't carry over to systems. |

Exactly one `System_Search` feature per output region marks the focus system, and any
`System_Search` features already in the input are cleared, as for `Domain_Search`.

### `select_by_system` (category `edit`, editor)

> Select regions around systems already annotated in GenBank records.

Options: `-i` (`-` allowed), outputs, `--systems/--systems_file/--system_categories`,
`--databases`, `--system_type`, `--min_score`, `--min_completeness`, `--partial_systems`,
`--invert`, the neighborhood options, `--model_neighborhood`, the `add_min_context` options,
and `--dedup`. It is to `system_search` what `select_by_cds` is to
`domain_search`. Every selected system becomes a region named by `model.region_id`, and gets
a `System_Search` feature as its focus, so downstream tools (`sort_contigs`,
`split_by_cluster --pad_on_search_hits`) treat the result like system_search output.

### `system_db_download` (category `database`)

> Download system models, or convert a local copy, into a native directory that other commands read.

| Option | Meaning |
|---|---|
| `--source {padloc,defensefinder,macsyfinder,antismash}` | the format to read |
| `--release` | a release tag to download (padloc and defensefinder only; default: the latest, which is then recorded) |
| `-i/--input` | a local copy instead of a download: a directory or a rule file (required for macsyfinder and antismash) |
| `--output_dir` | where to write `systems.yaml`, `profiles.hmm`, `LICENSE` and `conversion.tsv` |
| `--antismash_strictness {strict,relaxed,loose}` | antiSMASH only; default `relaxed` |
| `--coverage_profile` | MacSyFinder only; default 0.4, DefenseFinder's value |
| `--name_style {stem,qualified}` | MacSyFinder only: model names from file stems or fully qualified names |
| `--substitutes_file` | a native file whose components replace components of the same name, such as antiSMASH dynamic profiles ([antismash.md](antismash.md#dynamic-profiles)) |

It prints a summary of what was approximated or dropped, and writes the details to
`conversion.tsv`. The mappings are in
[definition_formats.md](definition_formats.md#converting-foreign-formats). Tests that download
are marked `network`; conversion tests use small fixtures with `--input`.

Why not `db_download` targets: `db_download`'s options (taxonomy filters, gene calling,
GenBank output) don't apply to model files, and conversion of local copies is half of this
command's job.

### `format_systems` (category `profile`, optional)

> Inspect native system models, and write them as YAML, a table, or the profiles they need.

Options:

- `--system_models`, `add_window_model` options, `--systems/--systems_file/--system_categories`;
- `--yaml_out`: the native format. Models that fit the short form are written in it, so this
  is how a command-line model becomes a file;
- `--tsv_out`: one row per model × component × reference, with thresholds, window, category
  and antiSMASH `RELATED` profiles. A way to *see* what an imported model means.
- `--hmm_out`: the concatenated profiles the selected models reference, ready for
  `domainate -r`;
- `--check_examples`: download each model's `examples` (NCBI accessions), run detection, and
  report models that don't find their own examples; `--check_examples_tsv_out` writes the
  results.

`--tsv_out` makes converted models transparent, and `--yaml_out` with `--systems` gives users a
path to fork one DefenseFinder, PADLOC or antiSMASH model into an editable file.

### `import_systems` (category `convert`, editor)

> Convert system calls made by other tools into `Domain_System` features.

| Option | Meaning |
|---|---|
| `-i`, `-o/--genbank_out`, `--fasta_out` | GenBank in and out; `-` allowed (single pass) |
| `--source {antismash}` | the tool that made the calls. `defense_finder` and `padloc` are planned. (Not `--from`: the library parameter would be a Python keyword.) |
| `--level {protocluster,candidate_cluster,region}` | antiSMASH only: which antiSMASH feature becomes a `Domain_System`. Default `protocluster`. |
| `--calls_file` | for DefenseFinder `*_systems.tsv` / PADLOC `.csv`: the call table, matched to CDSs by `locus_tag`/`protein_id` (`--name_field`) |
| `--clear_source` | remove the source tool's features after converting |

Detailed mapping in [antismash.md](antismash.md#importing-antismash-results).

## Changes to existing commands

Priority: **P1** ships with the first release of systems, **P2** with the second, **P3**
later.

| Command | Change | Pri |
|---|---|---|
| `domainate` | `--system_models` and detection options (above) | P2 |
| `domain_search` | none under Route C. The neighborhood code it shares with select_by_cds is generalized from "one focus CDS" to "an envelope of CDSs" (see [implementation.md](implementation.md)), with no CLI change. | P2 |
| `structure_search`, `structure_domainate` | none; see [Structure tools](#structure-tools) | - |
| `find_features` | none. Its `MOTIF`/`TMbed` hits are `Domainator` features, so a model can already use them as components (`{name: transmembrane_alpha_helix, database: TMbed}`, e.g. "a membrane anchor next to the effector"). | - |
| `enum_report` | New columns: `systems` (names of the systems the row overlaps: for `--by cds`, the systems the CDS belongs to), `system_categories`, `system_count:NAME`, `system_components` (role:component per member), and a `system_search` column family mirroring the existing `domain_search` sub-columns (name, database, score, evalue, completeness, core/total counts). A new row unit `--by system` slices each `Domain_System` envelope, so existing columns (`domains`, `architecture`, `length`, `taxname:...`) describe the system's region. | P2 |
| `summary_report` | A system count table (by name and database), systems by taxonomy, and a system co-occurrence matrix reusing the domain co-occurrence code. New `--json_out` keys only. | P2 |
| `select_by_contig` | `add_system_selection` + `--system_type` + `--min_system_score`: select records by the set of systems they contain. Joins the OR'd selectors and `--invert`, as domains do. Optionally (P3), `--window_cds`/`--window_kb` make the existing `--domain_expr` hold within one window instead of over the whole contig. That answers "which contigs have it" without writing system features. | P2 |
| `select_by_cds` | `--systems/--systems_file/--system_expr` select member CDSs. The CDS's system names come from `component` qualifiers by `cds_id`. Useful for "each GajA, with 2 CDSs either side". | P3 |
| `extract_peptides` | `--systems/--system_expr` select member CDSs. Proteins get the `system_membership` summary. | P2 |
| `extract_domains` | `--systems` keeps only domains that served as a member hit (the hit named in `component`), e.g. the GajA domain from Gabija systems only. | P3 |
| `filter_domains` | none (systems are self-contained; see [data_model.md](data_model.md#lifecycle)) | - |
| `trim_contigs` | `--systems/--system_expr`: trim from the ends until reaching a member of a selected system. | P3 |
| `rotate_contigs` | `--systems`: anchor the new origin on a system's envelope (start, end or center, through the existing anchor logic). | P3 |
| `sort_contigs` | sort by the `System_Search` score when present (`--system_type system_search`), as it sorts by `Domain_Search` today. | P3 |
| `split_by_cluster` | `--pad_on_search_hits` also centers on `System_Search`. | P3 |
| `color_genbank` | `--target` gains `systems` and `systems_and_cdss`: color `Domain_System` features, and optionally member CDSs, by `--system_color_by {name,category,database}` (default `name`; `category` gives the familiar antiSMASH-style PKS/NRPS/RiPP coloring). `--search_hit_color` also applies to `System_Search`. The color table is keyed by system name, as for domains. | P1 |
| `plot_contigs` | Add both new types to `DEFAULT_TYPES` with default colors. `plot_contigs.js` already draws any feature generically (parts, color, qualifiers), so systems appear as a bar spanning their members. Worth adding: a separate lane for system features, and a tooltip listing members. | P1 |
| `compare_contigs` | Optional `--sji` weight: Jaccard over system names, beside `--ji` (domains) and `--ai` (adjacency). Lets system content drive neighborhood networks. | P3 |
| `deduplicate_genbank`, `seq_dist`, `kmer_dist`, `structure_dist`, `hmmer_*`, `genbank_to_fasta`, `clean_sequences`, `extract_unannotated` | none | - |
| `build_ssn`, `build_tree`, `build_projection`, `viewer` | none. System information reaches the viewer as node metadata through `enum_report --by contig --columns systems ... --tsv_out` passed as `--metadata`, the same bridge domains use today. | - |
| `db_download` | none; model downloads are `system_db_download` | - |

## Structure tools

`structure_search` and `structure_domainate` work on **protein-only records**: one structure
or predicted structure per record, with `cds_id=0_1_<len>`. A system needs an ordered contig
of several CDSs, so there is nothing for these tools to cluster. Systems with structural
components are still supported, through the other routes:

1. **Structural references on contigs.** `domainate` already runs Foldseek on the proteins
   of CDSs in contigs (`--foldseek`, with ESM-2/ProstT5 3Di prediction). Those hits are
   ordinary `Domainator` features with `program=foldseek` and `database=<foldseek db>`. A
   component rule can name them (`{name: GajA_AF2, database: defense_structures}`). So
   `domainate --foldseek defense_structures --system_models models.yaml` detects systems
   whose components are defined by structure, and `system_search` would accept `--foldseek`
   the same way. `ReferenceRule` gains optional `min_tmscore`, `min_lddt` and `min_prob`
   thresholds, read from the qualifiers foldseek hits already carry.
2. **Mixed evidence.** One model can mix HMM components and structural components, because
   matching is by `(database, name)` on any `Domainator` feature.
3. **Not recommended: `--system_models` on structure_search or structure_domainate.** It
   would only be meaningful if their inputs were ordered proteins from one contig. That would
   mean adding contig parsing and CDS translation, which would duplicate `domainate
   --foldseek`. If a reseek-backed contig search is wanted later, the better path is a
   `--structure_backend` option on domainate.

So for structure tools the design is: **no new options, documented routes.**
`docs/structure_search.md` gets a short "Systems with structural components" section that
points at `domainate --foldseek --system_models`.

## Manifest and generated docs

- `cli/manifest.py`: add `system_models` to `INPUT_PATHS`, and `systems`, `system_expr`,
  `system_categories`, `system_name` and `release` to `TEXT_OPTIONS`. `--systems_file`,
  `--calls_file` and `--substitutes_file` are classified as inputs by their `_file` suffix,
  and `--output_dir` as a directory (as in `format_db`). New outputs end in `_out`
  (`--yaml_out`, `--hmm_out`).
- New commands go in `cli/registry.py` (categories above) and follow the module shape in
  `CLAUDE.md`: `<tool>()`, `build_parser()` ending in `declare_outputs`, `main(argv)`.
- Regenerate `test/data/manifest_snapshot.json` (`python test/test_manifest.py`) and
  `docs/cli/` (`python scripts/docs/cli_reference.py`). `domainate`'s snapshot changes
  because of its new options.
