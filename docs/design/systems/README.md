# Design: multidomain systems in Domainator

Status: **proposal, decisions recorded** (see [Decisions](#decisions) and
[Review amendments](#review-amendments)). Phase 1 is being implemented on the `systems`
branch.

A *system* is a set of domains that must be found together in one genomic neighborhood:
for example a defense system (Gabija = GajA + GajB), a CRISPR-Cas subtype, a toxin/antitoxin
pair, or a biosynthetic operon. Today Domainator annotates domains (`Domainator` features)
and selects single-hit neighborhoods (`Domain_Search` features). It cannot express or detect
"A and B within 5 CDSs of each other", and it cannot select a neighborhood around a group of
CDSs.

These documents explore how to add systems. Where there is a real choice they lay out the
options and recommend one. They are written to fit the 1.0 conventions in
[`docs/stability.md`](../../stability.md), [`docs/v1_release_plan.md`](../../v1_release_plan.md)
(Phases 2-4 and 9) and `CLAUDE.md`.

| Document | Covers |
|---|---|
| [data_model.md](data_model.md) | The new GenBank feature types (`Domain_System`, `System_Search`), their qualifiers and location, and how they behave under slicing and editing |
| [definition_formats.md](definition_formats.md) | The internal `SystemModel`, the three window modes (`chain`, `anchor`, `span`), the native YAML format and its short form, the condition grammar, nucleotide references, and converting DefenseFinder/MacSyFinder, PADLOC and antiSMASH models |
| [detection.md](detection.md) | The detection algorithm: competing hits, hit-to-component matching, candidate windows per mode, the requirement (quorum and condition), extenders, superiors and overlap resolution, scoring, and fidelity to DefenseFinder, PADLOC and antiSMASH |
| [antismash.md](antismash.md) | How antiSMASH's rule-based detection works, what the design borrows from it, converting its rule files, and importing its results |
| [tool_integration.md](tool_integration.md) | Three architectural routes, command-line models, the recommended set of new and changed commands, shared CLI vocabulary, and per-tool changes |
| [implementation.md](implementation.md) | Modules and library API, reuse of existing helpers, tests, docs, stability and a phased roadmap |
| [multi_locus.md](multi_locus.md) | **Future add-on:** systems made of distant loci, and MacSyFinder loners |

## Terms

| Term | Meaning |
|---|---|
| **system model** | A rule saying which components make up a system and how close together they must be. One file can hold many models. |
| **component** | A named role in a model (e.g. `GajA`), satisfied by a hit to any of one or more references (HMMs, sequences, structures, nucleotide profiles). Equivalent to a MacSyFinder *gene* (with its *exchangeables*) or a PADLOC *protein name*. |
| **role** | `core`, `accessory`, `neutral` or `forbidden` (MacSyFinder: mandatory/accessory/neutral/forbidden; PADLOC: core/secondary/neutral/prohibited). |
| **member** | A CDS, or a contig-level nucleotide hit, assigned to a component in one detected system. |
| **system hit** | One detected occurrence of a model on a contig, written as a `Domain_System` feature. |
| **envelope** | The contig interval from the first to the last member of a system hit, in the system's orientation. This is the system's *core*. |
| **window mode** | How "close together" is measured. `chain`: a maximum gap between consecutive members (MacSyFinder, PADLOC). `anchor`: a maximum distance from one anchor member (antiSMASH). `span`: a maximum length of the whole system (the default). Each is in CDSs and/or kb. |
| **reach** | In `span` mode, the members that could join a system with the result still fitting the span limit. Negation and `forbidden` look at the reach. |
| **condition** | A boolean expression over component names, with antiSMASH-style functions (`cds()`, `minimum()`, `repeatable()`, `minscore()`) and an order operator (`a < b`), evaluated over the members in one window. |
| **short form** | A model written as just a name and a condition, such as `XYZ_QR: "X & Y & Z & (Q \| R)"`. Every name in the condition is a domain and its own component. |
| **command-line model** | A short-form model given as options (`--domain_expr`, `--window_cds`), with no file. |
| **converted directory** | What `system_db_download` writes from a foreign format: `systems.yaml`, `profiles.hmm`, the source's license and a conversion report. |
| **neighborhood** | A model's suggested padding around the core (antiSMASH `NEIGHBOURHOOD`). It is recorded, not applied, unless a selector is given `--model_neighborhood`. |

## Recommendation at a glance

1. **Data model.** Add two feature types, `Domain_System` (every detected system, the analog
   of `Domainator`) and `System_Search` (the focus system of a selected region, the analog of
   `Domain_Search`). Both are 13 characters, because the INSDC feature-key limit is 15 and
   Biopython truncates longer keys. A system feature spans its members' CDSs in one
   contiguous region and lists them in a repeated `/component` qualifier. Systems are derived
   data: domains are never rewritten to point at systems.
2. **Definitions.** One internal `SystemModel` and one native YAML format, which is the only
   format commands read.
   - **Simple things are simple.** "10-CDS windows with X & Y & Z & (Q | R)" is one
     command-line option (`--domain_expr "X & Y & Z & (Q | R)" --window_cds 10`) or one line
     of a short-form file. Every key besides a name and a condition is optional.
   - **Hard things are possible.** The full form has an explicit **window mode** (`chain`,
     `anchor` or `span`), a quorum and/or a **condition** in a grammar that adopts
     antiSMASH's semantics and adds an order operator (`a < b`), and optional `category`, `superiors`, `extenders`,
     `neighborhood` and `examples`. Competing-profile groups resolve overlapping hits from
     related profiles.
   - **Foreign formats are converted once.** `system_db_download` downloads PADLOC and
     DefenseFinder releases, or converts local copies of those and of antiSMASH rules and
     other MacSyFinder packages, into a native directory with the profiles, the license and
     a report of every approximation.
   - **Any reference type can be a component**, including nhmmer and Infernal hits on the
     contig.
3. **Detection is a pure function of an annotated contig.** It reads `Domainator` features
   and applies per-reference thresholds (e-value, score, reference coverage, target
   coverage) from qualifiers already written, so it needs no new search.
4. **Commands.**
   - New `find_systems` (an editor: annotate systems on already-annotated records).
   - `domainate --system_models` (search the models' profiles and detect in one pass).
   - New `system_search` (the `domain_search` analog: `--cds_range`/`--kb_range`/`_up`/`_down`
     measured from the CDSs at either end of the system, and `--min_context_*` to require
     context on one or both sides).
   - New `select_by_system` (the editor counterpart, for records that already carry systems).
   - New `system_db_download` (download and convert foreign model formats).
   - New `import_systems` (turn antiSMASH's own output, and later DefenseFinder and PADLOC
     call tables, into `Domain_System` features).
   - `find_systems`, `system_search` and `format_systems` accept command-line models.
   - System-aware options in the reports, selectors, extractors and plotting tools.
   - `structure_search`/`structure_domainate` stay protein-only. Structural components work
     through `domainate --foldseek` on contigs (see [tool_integration.md](tool_integration.md#structure-tools)).
5. **Phasing.**
   1. The data model (including all three window modes), the native format with its short
      form and the order operator, `find_systems` with command-line models, and plotting.
   2. `system_db_download` for PADLOC, `domainate --system_models`, `select_by_system`, the
      reports, and importing antiSMASH results.
   3. `system_search`, and `system_db_download` for DefenseFinder and MacSyFinder.
   4. The full condition grammar and `system_db_download` for antiSMASH rules.
   - Future add-on: [multi-locus systems and loners](multi_locus.md), as linked per-locus
     features.

## Example pipeline (recommended route)

```bash
# the simple case: 10-CDS windows with X, Y, Z and Q or R, plus 5 CDSs of context each side
domainator system_search -i genomes.gb -r xyzqr.hmm -e 1e-20 \
    --domain_expr "X & Y & Z & (Q | R)" --window_cds 10 --cds_range 5 -o xyz_regions.gb

# published models: convert once, then use the directory anywhere
domainator system_db_download --source padloc --output_dir padloc_db/
domainator system_db_download --source defensefinder --output_dir df_models/

# one step: search the models' profiles and detect systems
domainator domainate -i genomes.gb --system_models padloc_db/ -o annotated.gb --cpu 0

# or two steps: annotate once with a big profile set, then try several definitions cheaply
domainator domainate -i genomes.gb -r defense_profiles.hmm -o annotated.gb
domainator find_systems -i annotated.gb --system_models my_systems.yaml -o systems.gb

# find system neighborhoods in a genome database, 5 CDSs beyond each end of the system
domainator system_search -i genomes.gb --system_models df_models/ \
    --systems Gabija Lamassu-Fam --cds_range 5 -o gabija_neighborhoods.gb

# biosynthetic gene clusters with antiSMASH's own rules, padded as antiSMASH would
domainator system_db_download --source antismash --input ~/antismash/antismash/detection/hmm_detection \
    --output_dir antismash_rules/
domainator system_search -i genomes.gb --system_models antismash_rules/ \
    --system_categories PKS NRPS --model_neighborhood -o bgc_regions.gb

# or bring in calls antiSMASH already made
domainator import_systems -i antismash_output.gbk --source antismash -o bgcs.gb

# reports and plots
domainator enum_report -i systems.gb --by system --columns systems system_components score -o systems.tsv
domainator select_by_contig -i systems.gb --system_expr "Gabija & ~Thoeris_I" -o subset.gb
domainator color_genbank -i systems.gb --target systems -o colored.gb
domainator plot_contigs -i colored.gb --html_out systems.html
```

## Decisions

The maintainer settled the open questions of the earlier drafts (November 2026):

1. **Feature-type names:** `Domain_System` and `System_Search`
   ([data_model.md](data_model.md#feature-type-names)).
2. **Search command:** a new `system_search` command (Route C), not a mode of
   `domain_search` ([tool_integration.md](tool_integration.md#route-b-vs-c-for-search)).
3. **Staleness:** systems are self-contained snapshots. `filter_domains` and other domain
   editors leave them alone, and `find_systems` recomputes them
   ([data_model.md](data_model.md#lifecycle)).
4. **Fidelity:** converted DefenseFinder, PADLOC and antiSMASH models are evaluated with the
   same semantics, and the differences are documented. Bit-for-bit reproduction is not a goal
   ([detection.md](detection.md#fidelity)).
5. **Licensing:** `system_db_download` may download model collections, including
   defense-finder-models, to the user's machine. Domainator never distributes model content
   ([definition_formats.md](definition_formats.md#licensing)).
6. **Grammar scope:** antiSMASH's condition semantics, plus a Domainator **order operator**
   `a < b` (a upstream of b in system orientation; N- to C-terminal inside `cds()`)
   ([definition_formats.md](definition_formats.md#order)).
7. **Window keyword:** `span` ([definition_formats.md](definition_formats.md#windows)).
8. **antiSMASH scope:** import antiSMASH results in Phase 2, convert its rules in Phase 4
   ([antismash.md](antismash.md#is-it-worth-it)).
9. **Command-line model option:** `--domain_expr`
   ([tool_integration.md](tool_integration.md#command-line-models)).
10. **Foreign formats on `--system_models`:** refused, with an error that gives the
    `system_db_download` command ([definition_formats.md](definition_formats.md#reading-models)).
11. **Nucleotide hits inside CDSs:** separate members unless the reference sets
    `attach_to_cds` ([definition_formats.md](definition_formats.md#nucleotide-references)).
12. **Span mode:** a negated name or forbidden component vetoes a system when it is within
    the system's *reach* (the system plus it would still fit the span limit), and the system
    is the best-scoring passing run of members, shortest first among equal scores ([detection.md](detection.md#span-domainator-the-whole-system-within-x)).
13. **Partial flag:** the qualifier is `/system_partial`, not `/partial`
    ([data_model.md](data_model.md#qualifiers)).
14. **Context requirement:** `system_search` and `select_by_system` take
    `--min_context_cds`/`--min_context_kb` with `_up`/`_down` forms
    ([tool_integration.md](tool_integration.md#requiring-context)).

## Review amendments

A review against the code (November 2026) found these problems in the earlier drafts. The
documents above have been corrected.

1. The per-CDS summary qualifier can't be `domainator_systems`: `seqio._with_domain_qualifiers`
   removes unknown `domainator_*` CDS qualifiers on write. It is `system_membership`.
2. Span mode tested maximal windows, which missed valid systems under negation and let stray
   duplicates stretch systems. Decision 12 replaces it.
3. Converted directories named their profiles `profiles.hmm`, so hits were labeled
   `database=profiles` when the file was given to `domainate -r`, but by the collection's
   database under `domainate --system_models`. The file is now `<database>.hmm`.
4. A fusion CDS filling several components of a condition-only model listed only one of them.
   `/component` now has one value per (member, component).
5. A name used only under `~` in the short form had no role. It is a `forbidden` component.
6. `nuc_` member ids are coordinates and go stale after slicing. Readers resolve nucleotide
   members by name inside the system's location.
7. `-e` could not tighten model profiles in `domainate --system_models`, unlike
   `find_systems`. It is now an evidence filter in detection on every route.
8. `extenders` appeared in a format-1 example but ships with format 1.1.
9. The format version is parsed as `major.minor` integers, not a YAML float, and unknown keys
   are refused.
10. `system_models` takes paths or `ModelCollection`s (not `SystemModel`s), and the condition
    evaluator lives in `systems/conditions.py`, sharing `expr`'s tokenizer.
11. `superiors` bare names resolve within the model's own collection.
12. Anchor mode now requires the anchor to contribute to the condition, as antiSMASH does.
13. `find_systems --clear` clears command-line models by name, not by the shared
    `command_line` database.
14. `--partial_systems` moved to Phase 4, with the behavior it selects.
15. Smaller corrections: where detection runs (domainate is in-process; only system_search
    uses a worker pool), system_search must annotate every hit before detection, which systems
    of a `--system_expr` become focus regions, the writer checks `|` in hit databases, and
    `program`/`semantics:` records provenance.

## Implementation status

Phases 1 to 4 are implemented, except the optional `merge_systems`. Where the code departs
from these documents:

1. **PADLOC model names.** Converted models are named by padloc's `yaml.name` (e.g.
   `cas_type_I-E`), not `sys_meta.system`: padloc reports systems by `yaml.name`, and names
   such as `CBASS (other)` would break `--system_expr`. `sys_meta.system` is the description.
2. **`domainate --system_models`** searches the model profiles in a second pass over the
   stream, rather than through exemption keys in the reference dict: the second pass has no
   `--max_hits_per_cds` or `--max_overlap`, and `run_search` takes per-profile e-values.
3. **`select_by_system --partial_systems`** is a selector on the `/system_partial` flag
   (`include`, `exclude`, `only`, like `--partial`), not the Phase 4 detection option of the
   same name. The detection commands still don't have it.
4. **`add_neighborhood`** was not extracted: the padding options' help differs between
   `domain_search`, `select_by_cds` and `select_by_system`, and extracting it would change
   existing help text. `select_by_cds.get_neighborhood` is shared instead.
5. **`/system_membership`** (on the CDS features of `extract_peptides` proteins) is
   `role|component|database|name`, mirroring `/component`, rather than `Gabija:GajA:core`.
6. **`select_by_contig --system_expr`** and **`extract_peptides --system_expr`** were added
   beside `--systems`; a name in an expression is `name` or `database:name`.
7. **`import_systems`** names candidate clusters and regions by their protoclusters' products
   joined with `+` (e.g. `hglE-KS+T1PKS`), and replaces systems it imported before.
8. **`system_search`** searches `-r` and the model profiles in one pass, since it has no
   `--max_overlap` or `--max_hits_per_cds` to exempt the profiles from: overlapping hits are
   for detection to resolve. It also drops `--alphabet` (inputs are nucleotide; FASTA needs
   `--gene_call`) and `--max_hits_per_contig` (`--max_systems_per_contig` covers it). With
   `--systems` and `--system_expr` together, the models named by either are searched, the
   expression must hold on the contig, and systems of the `--systems` models and of the
   expression's positive names are selected. With `--translate`, `--max_hits` still counts
   systems.
9. **`trim_contigs --systems`** has no `--system_expr`: it trims to the outermost selected
   system's location, so that the systems stay whole (nucleotide members included), which a
   per-CDS expression wouldn't guarantee. `extract_domains --systems` and `rotate_contigs
   --systems` take `--system_type`, and `rotate_contigs` doesn't filter systems by `-e`.
10. **`sort_contigs --system_type`** defaults to none (sort by `Domain_Search`, as before);
    any value sorts by the best score of those systems.
11. **`select_by_contig --window_*`**: with a window, protein records aren't selected by
    `--domain_expr`, since they have no windows.
12. **MacSyFinder conversion** (checked against defense-finder-models 3.1.0, all 666 models):
    `--source defensefinder` names the collection `defensefinder`, and `--source macsyfinder`
    names it after the package directory. Component names drop the gene's `<family>__` prefix
    (`Gabija__GajA` → `GajA`). The category is the model's top-level folder (`DefenseFinder`,
    `RM`, `ADF`), not `defense`. Two models with the same file name get their path as their
    name (with a note), instead of an error suggesting `--name_style`. A gene listed twice in a
    model (five DefenseFinder models do this) is merged, with a note. The download uses the
    release's tarball asset.
13. **`import_systems --source defensefinder/padloc`** reads DefenseFinder's
    `*_defense_finder_genes.tsv` (not `*_systems.tsv`, which lacks roles and profiles) and
    PADLOC's `*_padloc.csv`, given with `--calls`, into the `defensefinder` and `padloc`
    databases. PADLOC's table has no scores or roles: its systems score 0 and every member is
    core. Members without a CDS (CRISPR arrays, ncRNAs) are skipped.
14. **`--truncated_systems {exclude,include}`** (`find_systems`, `domainate`, `system_search`)
    is the detection option the design calls `--partial_systems include`: that name is
    `select_by_system`'s flag selector (deviation 3). A window passes as truncated when it
    passes with a stand-in member just beyond the contig end that fills every component. In
    span mode, the window at each end (as far as the limit reaches) is tested, because span
    selection would always pick a stand-in alone (it passes, and is the shortest run). The span radius behind `/system_partial`
    now leaves a CDS for the missing member.
15. **Short-form roles.** A name used only negated as a bare name (`~W`) is forbidden. One used
    only negated inside `cds()`, a function or `<`, or only in `extenders`, is neutral: in
    Phase 1, `~(W < V)` made W and V forbidden, so any W vetoed the system.
16. **Format 1.1 additions beyond the design:** `related` (antiSMASH `RELATED`, informational),
    and feature components (`{feature: repeat_region, qualifiers: {rpt_family: CRISPR}}`),
    which take no thresholds and add 0 to the score. The PADLOC converter keeps `CRISPR_array`
    as such a component instead of dropping it and lowering `minimum_total`, which is what
    padloc does without CRISPRDetect output.
17. **Extenders** grow a system in every window mode by the window's `max_cds`/`max_kb` from
    the last member on that side; members that don't meet the extenders are passed over while
    within that distance. Condition contributors are names, so in anchor mode every member
    filling a contributing name is a core member, as are `cds()` contributors on other CDSs.
18. **antiSMASH conversion** reads an antiSMASH checkout, installation or `hmm_detection`
    directory (not a lone rule file), records antiSMASH's version from its `main.py`, writes
    a `NOTICE` on the AGPL beside the license, and approximates `CUTOFF 0` (one loose rule) as
    a 0 kb window. Substitutes files are `domainator_substitutes: 1` with a `components` list.
    Checked against antiSMASH 8.0.4: all 103 rules convert and read back unchanged, and seven
    rule examples (benzoxazole, T2PKS, transAT-PKS, lanthipeptide class I) are found.
19. **`--bit_cutoffs`** is on `system_search` too. With it, `-e` doesn't filter system evidence,
    because the cutoffs replace it.
20. **`format_systems --check_examples_out <tsv>`** (outputs end in `_out`) is the design's
    `--check_examples`. It annotates each example with domainate's defaults.
21. **`merge_systems`** is not built.
