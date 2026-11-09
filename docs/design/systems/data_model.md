# Systems: GenBank data model

GenBank is Domainator's universal carrier. Every tool reads and writes `SeqRecord`s, and the
annotations travel as features. Systems should travel the same way, so that systems
annotated by one command can be selected, reported, colored and plotted by the next, and
streamed through `-i -` / `-o -`.

## Requirements

- Survive a GenBank round trip with no side files.
- Be distinguishable from domain hits, so tools that count or filter `Domainator` features
  are not affected (`enum_report`'s domain columns, `filter_domains`, `compare_contigs`, ...).
- Carry enough to report a system without its definition file: name, source, score, members
  and each member's role.
- Have a contig location, so `plot_contigs` can draw it and slicing tools can tell whether
  it was cut.
- Stay within the 1.x stability contract. New feature types and qualifiers are additions,
  which [`stability.md`](../../stability.md) already allows: "a reader must ignore qualifiers
  it doesn't know" and "feature types".

## Feature types

Two new feature types, defined next to the existing constants in
`src/domainator/__init__.py` and exported in `__all__`:

```python
DOMAIN_FEATURE_NAME = "Domainator"              # existing
DOMAIN_SEARCH_BEST_HIT_NAME = "Domain_Search"   # existing
SYSTEM_FEATURE_NAME = "Domain_System"           # new: every detected system
SYSTEM_SEARCH_BEST_HIT_NAME = "System_Search"   # new: the focus system of a system_search region
```

The two pairs line up like this:

| | domains | systems |
|---|---|---|
| Annotation, as many as found | `Domainator` (domainate, find_features) | `Domain_System` (find_systems, domainate `--system_models`, import_systems) |
| The hit a region was selected for | `Domain_Search` (domain_search, one per region) | `System_Search` (system_search, one per region) |

### Feature-type names

| Option | Length | Verdict |
|---|---|---|
| `Domain_System` + `System_Search` | 13 + 13 | **Chosen.** Parallels `Domain_Search`, says what the feature is, and leaves `Domainator` meaning "a domain hit". |
| `System` + `System_Search` | 6 + 13 | Workable, but `System` is generic enough to collide with feature keys from other tools. |
| `Domainator_System` | 17 | **Not possible.** The INSDC limit for feature keys is 15 characters. Biopython warns, and `_InsdcWriter._write_feature` slices the key to the qualifier indent, so the key is truncated on write. |
| Reuse `Domain_Search` for system searches | - | Rejected. `Domain_Search` means one hit in one CDS (`cds_id`), `DomainatorCDS.domain_search_feature` holds one per CDS, and `extract_peptides` copies it onto the protein. A multi-CDS feature of the same type breaks those readers. |

## Location

The location is the **envelope of the member CDSs**: from the outer edge of the first member
to the outer edge of the last, covering whole CDSs, not just the domain envelopes.

- **Strand** is the system's orientation, which `system_search`'s `--cds_up`/`--cds_down`,
  `--anchor` and `--keep_direction` follow. Rule:
  1. if all counted CDS members (core and accessory) are on one strand, that strand;
  2. otherwise the strand of the majority of counted CDS members;
  3. on a tie, the strand of the highest-scoring core member.

  Nucleotide members (ncRNA pseudo-CDSs) don't vote, because their strand is often
  uninformative. A model may declare `orientation_from: <component>` to make one component
  define "upstream" (for example the effector in a CRISPR-Cas locus, or an ncRNA whose
  strand matters).
- **Circular contigs:** a system that wraps the origin gets a `CompoundLocation` split at
  the origin. This is how `domainate.build_contig_hit_location` already writes wrapped
  contig-level hits, and `location_ops.circular_envelope` already computes such envelopes.
- **Members that are contig-level nucleotide hits** (nhmmer or Infernal hits with
  `cds_id="."`, such as a retron ncRNA) are part of the envelope like CDSs. A nucleotide hit
  whose reference sets `attach_to_cds` is a member through the CDS it overlaps (see
  [definition_formats.md](definition_formats.md#nucleotide-references)).
- **One region.** A system location is always one contiguous region, or two parts when it
  wraps the origin. Systems made of distant loci, and MacSyFinder loners, are a
  [future add-on](multi_locus.md) that adds linked per-locus features without changing this
  rule for single-locus systems.

Why the CDS envelope and not the domain envelope: windows (`--cds_range`, `--kb_range`) are
measured from CDSs, and slicing must keep whole member CDSs for the system to stay valid.

**Core, not neighborhood.** antiSMASH distinguishes a protocluster's *core* (the defining
CDSs) from its *neighborhood* (the core padded by a rule-specific distance). The
`Domain_System` location is always the core, including CDSs added by `extenders`. A model's
suggested neighborhood is recorded in qualifiers and applied only when a selection tool is
asked to (`--model_neighborhood`). This keeps "what the system is" separate from "how much
context to cut out", which is a per-analysis choice in Domainator.

## Qualifiers

Where a qualifier means the same thing as on a `Domainator` feature, it keeps the same name
and formatting (`domainate.domain_feature_qualifiers`). Readers like `model.feature_score`
then work unchanged.

| Qualifier | Value |
|---|---|
| `name` | The model name (`Gabija`, `cas_type_I-A`). |
| `database` | The model collection: the native file's `database:` key, or its stem without extensions (the same rule as `seqio.db_name_from_path` for references). For converted collections it is the source's name (the MacSyFinder package name `DefenseFinder`, or `padloc`, or `antismash`). For a model defined on the command line it is `command_line`, and `description` holds the expression. |
| `program` | Where the model's definition came from: `domainator`, `padloc`, `macsyfinder` or `antismash`, from the collection's `semantics:` key (set by `system_db_download`). Detection is the same for every value (see [detection.md](detection.md#fidelity)), so this is provenance: it tells a reader which tool's conventions the model's thresholds and windows were translated from. |
| `category` | The model's category (`defense`, `PKS`, `RiPP`, `toxin_antitoxin`, ...), or `.`. Used for grouping in reports and coloring. |
| `description` | The model's description, or `.`. |
| `accession` | A stable model identifier when the source has one (the MacSyFinder fully qualified name `DefenseFinder/Gabija/Gabija`, the PADLOC `yaml.name`), or `.`. |
| `score` | System score, `%.1f`: the sum of the best bit scores of the counted members (see [detection.md](detection.md#scoring)). |
| `evalue` | The largest (worst) e-value among the counted members, `%.1e`. A `-e` filter then means "every counted member passes". |
| `completeness` | Distinct counted components found ÷ the model's component total (MacSyFinder's `max_nb_genes` when given), `%.3f`. `.` for condition-only models (antiSMASH rules), where a component total is undefined. |
| `window` | The window that was applied, e.g. `chain:cds=5`, `anchor:kb=19.999`, `span:kb=2.5`. Informational: it lets a reader judge a call without the definition file. |
| `neighborhood_kb`, `neighborhood_cds` | Optional. The model's suggested padding (antiSMASH `NEIGHBOURHOOD`). It is **not** part of the location; `--model_neighborhood` applies it. |
| `kind` | Optional. For systems imported from antiSMASH candidate clusters: `single`, `neighbouring`, `interleaved` or `chemical_hybrid`. These are antiSMASH's own identifiers, so they keep its spelling. |
| `core_count`, `total_count` | Distinct core components, and distinct core + accessory components, found. These are the quantities the quorum tests. |
| `member_count` | Number of member CDSs, including neutral members and duplicates. |
| `component` | **Repeated**, one per (member, component) pair, in contig order. See below. |
| `system_partial` | Optional. `<`, `>` or `<>` when the system touches a linear contig end and a missing component could lie beyond it (see [detection.md](detection.md#contig-ends-and-partial-systems)). The values follow the `<`/`>` convention of `domainator.partial`. It is not named `partial`, because `/partial` is the INSDC qualifier that Domainator reads as CDS fragment status. |

### The `component` qualifier

Each (member, component) pair is one value of a repeated qualifier, with five fields separated
by `|`. A member usually fills one component. In a model with no quorum (a condition-only
model, such as the short form or an antiSMASH rule), a CDS that carries hits to several
components named positively in the condition fills all of them (a fusion), and each gets a
value. In a quorum model, a member lists only its primary component (see
[detection.md](detection.md#4-assign-components-per-member-per-model)).

```
/component="<role>|<component>|<cds_id>|<database>|<reference name>"
```

- `role` is `core`, `accessory`, `neutral` or `extender` (a CDS added by the model's
  `extenders` condition). Forbidden components never appear, because their presence vetoes
  the system.
- `cds_id` is the member CDS's `cds_id`, the same key that links `Domainator` features to
  CDSs (`model.get_cds_unique_name`). `cds_id`s are preserved when records are sliced, so the
  reference survives `select_by_cds`, `domain_search`, `trim_contigs` and so on. A
  contig-level nucleotide member writes its own coordinates in the same
  `start_strand_end` form, prefixed `nuc_`, because those hits share `cds_id="."`. Unlike a
  `cds_id`, these coordinates are not stored on the hit, so they describe the contig as it
  was at detection and go stale when the record is sliced, rotated or reverse complemented.
  Readers find a nucleotide member by its `database` and `reference name` among the
  contig-level hits inside the system's location. A nucleotide hit attached to a CDS
  (`attach_to_cds`) uses that CDS's `cds_id`.
- `database` and `reference name` identify the hit that satisfied the component. Together
  with `cds_id` they match the `database`, `name` and `cds_id` of a `Domainator` feature.
- Parse with `value.split("|", 4)`. The first four fields cannot contain `|`: roles and
  `cds_id`s by construction, and component names at load time with an `InputError`. The
  fourth field is the *hit's* database, which comes from a reference file name rather than
  from the model, so the writer checks it too and raises an `InputError` naming the file.
  The last field may contain `|`.

Example (Gabija on the minus strand):

```
     Domain_System   complement(12011..18930)
                     /name="Gabija"
                     /database="DefenseFinder"
                     /program="macsyfinder"
                     /description="."
                     /accession="DefenseFinder/Gabija/Gabija"
                     /category="defense"
                     /window="chain:cds=5"
                     /score="1203.4"
                     /evalue="3.2e-41"
                     /completeness="1.000"
                     /core_count="2"
                     /total_count="2"
                     /member_count="2"
                     /component="core|GajA|18930_-1_16102|DefenseFinder|Gabija__GajA"
                     /component="core|GajB|15998_-1_12011|DefenseFinder|Gabija__GajB_1"
```

### Alternatives considered for linking members

| Alternative | Problem |
|---|---|
| Parallel repeated qualifiers (`/member_cds_id`, `/member_role`, ...) | Readers must keep the lists aligned, and a hand edit or a tool that drops one value corrupts every later member. |
| A `/system=` qualifier on each member `Domainator` feature | Rewrites existing features. One domain may belong to several systems (MacSyFinder `multi_model`). It gives no envelope, score or quorum counts, and goes stale when systems are recomputed. |
| A per-record JSON blob in `record.annotations` | Not written to GenBank, so it is lost after one command. |
| `misc_feature` with `/note` | Can't be told apart from users' own notes. Breaks the "one feature type per annotation kind" rule. |

### Optional per-CDS summary

`seqio._genbank_output_record` already adds `domainator_<db>` summary qualifiers to CDS
features on write. A `/system_membership="Gabija:GajA:core"` qualifier on member CDSs would
let `extract_peptides` and protein-level tools keep system membership after the contig is
gone, because protein records can't carry a multi-CDS feature. It is optional and can be
added later without breaking anything.

It must not be named `domainator_systems`. `seqio._with_domain_qualifiers` owns the
`domainator_` prefix: on every write it removes each `domainator_*` CDS qualifier that isn't a
summary of a current database, and a reference database named `systems` would write the
same qualifier.

## Lifecycle

**Input clearing.** `domain_search` clears `Domain_Search` from its input
(`clean_rec(clear_best_hit=...)`), and `system_search` clears `System_Search` the same way.
`find_systems` replaces existing `Domain_System` features from the same `database` by default
(`--clear`/`--no_clear`), so re-running it with edited definitions doesn't duplicate systems.
Command-line models all have the database `command_line`, so for them only systems with the
same `database` and `name` are replaced; a second `--domain_expr` run with another
`--system_name` keeps the first run's systems.

**Slicing.** `locations.slice_record_from_location` drops hit features that a slice cuts
through (`locations.py:437`). `Domain_System` and `System_Search` get the same rule: a
system cut by a slice is dropped, not truncated, because a truncated system no longer meets
its quorum. `system_search` and `select_by_system` always pad from the system's ends, so the
system they select is never cut. Neighbor systems at the edge of a region are dropped, which
matches how domain hits behave today.

**Editing domains.** A `Domain_System` is a snapshot of the hits it was detected from.

- **Decided:** systems are self-contained. Editors that remove domains (`filter_domains`,
  `extract_domains`) leave systems alone, and `find_systems` recomputes them. This needs no
  changes to existing editors, and the `component` qualifier still reports which hits were
  used.
- **Rejected:** `filter_domains` dropping (or `--revalidate_systems` re-checking) systems
  whose member hits it removed. That is more consistent, but it couples domain editors to
  system semantics and needs the model definitions at hand to re-check a quorum.

**Protein records.** A protein record has a single CDS, so a protein-level `Domain_System`
can't exist. The membership summary above is how system information reaches proteins.

## Stability and documentation

- `docs/file_formats.md` gets a "Systems" section: the two feature types, the qualifier
  table above, the `component` grammar, the location rule and slicing behavior.
- `docs/stability.md`: the two feature types and their qualifiers join the GenBank contract
  (changes only by addition). The two constants join the Python API list.
- No format version is needed. Old Domainator versions ignore unknown feature types; check
  that `plot_contigs`' `--types` default and any `feature.type` switch tolerate them (they do
  today).
