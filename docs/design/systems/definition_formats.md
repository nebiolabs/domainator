# Systems: definition formats

Four formats are in scope:

- a new Domainator-native format;
- DefenseFinder (MacSyFinder models);
- PADLOC (`padloc-db`);
- antiSMASH cluster rules (`hmm_detection`).

They all describe the same thing: named components, each satisfied by one of several
references, plus a requirement (a quorum, a boolean condition, or both) and a window that
says how close together the members must be. So the design is **one internal model and one
native format**, not four engines.

- **Commands read only the native format.** The three foreign formats are converted once,
  by `system_db_download`, into a native directory (see
  [Converting foreign formats](#converting-foreign-formats)).
- **The native format has a short form** for the common case: "these domains within N CDSs"
  is one line (see [The short form](#the-short-form)), or a command-line option with no file
  at all ([tool_integration.md](tool_integration.md#command-line-models)).

The three foreign formats measure "close together" in three different ways, so the internal
model makes the **window mode** explicit (see [Windows](#windows)).

## The internal model

Internal module `domainator.systems.model`, as frozen dataclasses so models can be shared
read-only with workers through `_mp.Shared`:

```python
Role = Literal["core", "accessory", "neutral", "forbidden"]
WindowMode = Literal["chain", "anchor", "span"]

@dataclass(frozen=True)
class ReferenceRule:
    """One way to satisfy a component: a hit to this reference that passes these thresholds."""
    name: str                                     # the hit's `name` qualifier (HMM NAME, sequence id, structure name)
    database: Optional[str] = None                # the hit's `database` qualifier; None matches any database
    evalue: Optional[float] = None                # keep hits with evalue <= this (model.evalue_passes)
    min_score: Optional[float] = None             # bit score; MacSyFinder GA and antiSMASH cutoffs land here
    min_reference_coverage: Optional[float] = None  # (rend - rstart + 1) / rlen
    min_target_coverage: Optional[float] = None     # hit length / protein length
    min_target_length: Optional[int] = None       # protein length in aa (antiSMASH dynamic profiles, precursor peptides)
    max_target_length: Optional[int] = None
    exchangeable: bool = False                    # MacSyFinder: satisfies the component as a stand-in
    attach_to_cds: bool = False                   # nucleotide hits: count as the CDS they overlap (see "Nucleotide references")

@dataclass(frozen=True)
class Component:
    name: str
    role: Role
    references: Tuple[ReferenceRule, ...]
    loner: bool = False            # MacSyFinder: may sit outside the cluster; future add-on (multi_locus.md), treated as ordinary until then
    multi_system: bool = False     # one hit may serve several occurrences of the model
    multi_model: bool = False      # one hit may be shared with other models
    max_gap_cds: Optional[int] = None  # per-component chain override (MacSyFinder)

@dataclass(frozen=True)
class Window:
    mode: WindowMode
    max_cds: Optional[int] = None        # meaning depends on mode; see "Windows"
    max_kb: Optional[float] = None       # meaning depends on mode; see "Windows"
    max_span_cds: Optional[int] = None   # optional extra cap on the whole system (chain and anchor only)
    max_span_kb: Optional[float] = None

@dataclass(frozen=True)
class SystemModel:
    name: str
    database: str
    program: str                    # "domainator" | "padloc" | "macsyfinder" | "antismash": which semantics apply (the collection's `semantics:` key)
    components: Tuple[Component, ...]
    window: Window
    min_core: int = 0
    min_total: int = 0
    condition: Optional[str] = None      # expression over component names; see "Conditions"
    extenders: Optional[str] = None      # condition a neighboring CDS must meet to extend the system (antiSMASH)
    superiors: Tuple[str, ...] = ()      # models that win when they overlap this one (antiSMASH SUPERIORS; PADLOC *_other)
    category: str = "."                  # e.g. "defense", "PKS", "RiPP" (antiSMASH CATEGORY)
    neighborhood_kb: Optional[float] = None   # default padding for regions around this system (antiSMASH NEIGHBOURHOOD)
    neighborhood_cds: Optional[int] = None
    description: str = "."
    accession: str = "."
    same_strand: bool = False            # PADLOC force_strand
    max_components: Optional[int] = None # completeness denominator (MacSyFinder max_nb_genes)
    multi_loci: bool = False             # MacSyFinder; future add-on (multi_locus.md), evaluated as single-locus until then
    orientation_from: Optional[str] = None
    examples: Tuple["Example", ...] = () # known occurrences: accession, start, end, note
    literature: Tuple[str, ...] = ()     # DOIs / PubMed ids
    reference_files: Tuple[str, ...] = ()  # profiles the model needs, for domainate --system_models

@dataclass(frozen=True)
class ModelCollection:
    """What one --system_models input yields."""
    database: str
    models: Tuple[SystemModel, ...]
    semantics: str = "domainator"               # "padloc" | "macsyfinder" | "antismash" for converted collections; sets SystemModel.program
    source: Optional["Source"] = None           # provenance of a converted collection: source, URL, release, date, converter version
    competing: Tuple[FrozenSet[str], ...] = ()  # groups of reference names whose overlapping hits compete (antiSMASH equivalence groups)
    competing_overlap_aa: int = 20
```

A model needs a `window` (its own, or the collection's `defaults`) and at least one
requirement: `min_core > 0`, `min_total > 0`, or a `condition`. When both a quorum and a condition are given, both must hold. A `condition`
must contain at least one positive (non-negated) term, because something has to anchor the
system. antiSMASH enforces the same rule.

The role vocabularies line up like this:

| Domainator | MacSyFinder `presence` | PADLOC key | antiSMASH | Counts toward | Joins the window | Vetoes |
|---|---|---|---|---|---|---|
| `core` | `mandatory` | `core_genes` | a profile named positively in `CONDITIONS` | `min_core`, `min_total`, score | yes | - |
| `accessory` | `accessory` | `secondary_genes` | - | `min_total`, score | yes | - |
| `neutral` | `neutral` | `neutral_genes` | `EXTENDERS` profiles; `RELATED` (informational) | nothing | yes | - |
| `forbidden` | `forbidden` | `prohibited_genes` | (negation in a condition, see below) | nothing | yes | yes |

## Windows

"Close together" means three different things in the formats we read. The model names which
one applies with `window.mode`:

| Mode | `max_cds` means | `max_kb` means | Used by |
|---|---|---|---|
| **`chain`** | at most this many CDSs **between consecutive members** | at most this many kb between the end of one member and the start of the next | MacSyFinder (`inter_gene_max_space`), PADLOC (`maximum_separation`) |
| **`anchor`** | every member has at most this many CDSs **between it and one anchor member** | every member is within this many kb of the anchor | antiSMASH (`CUTOFF`) |
| **`span`** | the whole system, **first member to last inclusive**, covers at most this many CDSs | the whole system envelope is at most this many kb long | Domainator-native "all within N kb / N CDSs" |

What negation and `forbidden` look at depends on the mode, because each mode has its own idea
of "the window":

- `chain`: the chain cluster. Forbidden members chain like any other, so a forbidden gene
  within the gap of a member joins the cluster and vetoes it.
- `anchor`: everything within the radius of the anchor, as in antiSMASH.
- `span`: the system's **reach** (decided). A member is within reach when the system plus
  that member would still fit the span limit. So with `max_cds: 10`, a 3-CDS system reaches
  7 CDSs beyond either end. `~W` means "no W within reach", and a forbidden member within
  reach vetoes the system. The system itself is the best-scoring run of members that passes,
  and the shortest of equal scores: it takes every member that fills another component, but no
  stray duplicates.

- When both `max_cds` and `max_kb` are given, both must hold.
- `mode` defaults to `span`, the literal reading of "within N CDSs", so
  `window: {max_cds: 10}` is a complete window.
- Contig-level nucleotide members (ncRNA pseudo-CDSs) never count as CDSs, but they do count
  toward kb distances. A nucleotide hit whose reference sets `attach_to_cds` counts as the
  CDS it overlaps instead (see [Nucleotide references](#nucleotide-references)).
- On records with no CDSs (nucleotide FASTA contigs with only nhmmer or cmsearch hits), only
  `max_kb` can be measured. A model whose window sets only `max_cds` is skipped there, with
  one warning per model.
- In `chain` and `anchor` mode, `max_span_cds`/`max_span_kb` optionally add a cap on the
  whole system. Example: chain with gaps of up to 3 CDSs, but never more than 40 kb in total.
  In `span` mode the span is the window, so the caps are not allowed.

Worked example: components A, B, C and D sit at CDS positions 1, 4, 7 and 10. The model
requires all four (`condition: "A & B & C & D"`).

| Window | Result | Why |
|---|---|---|
| `chain`, `max_cds: 2` | **found** | every consecutive pair has 2 CDSs between |
| `anchor`, `max_cds: 2` | not found | no member has all three others within 2 intervening CDSs (B reaches A and C but not D) |
| `anchor`, `max_cds: 5` | **found** | from B or C, every other member is within 5 intervening CDSs |
| `span`, `max_cds: 10` | **found** | positions 1-10 inclusive is 10 CDSs |
| `span`, `max_cds: 9` | not found | the system needs 10 |

Span mode, negation and duplicates: members X, Y and W at CDS positions 1, 2 and 9, with
`span`, `max_cds: 10`.

| Condition | Result | Why |
|---|---|---|
| `X & Y` | **found**, CDSs 1-2 | X, Y is the shortest of the best-scoring runs; W is not part of the system |
| `X & Y & ~W` | not found | W is within reach: the system plus W covers 9 CDSs, within 10 |
| `X & Y & ~W`, `max_cds: 5` | **found**, CDSs 1-2 | W is out of reach: the system plus W would cover 9 CDSs, more than 5 |

With members A, B and A at positions 1, 3 and 9 and `A & B`, the system is CDSs 1-3. The stray
copy of A at 9 doesn't stretch it.

`chain` can drift: a long run of members can stretch well beyond the gap. `anchor` bounds
every member's distance from one center, so it never covers more than about twice the
radius. `span` bounds the whole system directly, which is the most literal reading of "these
domains within x kb".

**Keyword (decided).** The third mode is `span`, because it names the quantity being
limited and matches the `max_span_*` caps. `strict` was rejected: antiSMASH already uses it
for a rule tier (its strict/relaxed/loose rule files), and it suggests stricter thresholds
rather than a different window shape. `fixed`, `bounded` and `total` were also considered.

## Reading models

```python
def read_system_models(paths: Sequence[str | os.PathLike]) -> List[ModelCollection]: ...
```

It also accepts already-read `ModelCollection`s in place of paths, so `system_models`
parameters of the command functions take paths, collections or a mix. Every command reads
models through this one function, and it reads **only the native format**:

| Path | Read as |
|---|---|
| a `.yaml`/`.yml`/`.json` file with a top-level `domainator_systems:` key | one collection |
| a directory holding `systems.yaml` (what `system_db_download` writes) | that file; its `reference_files` resolve inside the directory |

A path that looks like a foreign format (a padloc-db directory, a MacSyFinder package, an
antiSMASH `hmm_detection` directory or rule file) is refused with an `InputError` that gives
the exact `system_db_download` command that converts it. Python users can call the
converters directly (`domainator.systems.formats.convert(source, path, ...) ->
ModelCollection`) and pass the result to any command function, because `system_models`
parameters also accept already-read collections.

`--system_models` takes any number of native paths. Model `(database, name)` pairs must be
unique across all inputs, as reference `db_name`s are today (`hmmio.read_hmms`). Errors use
`InputError`, and point at the file, line and column, as antiSMASH's rule parser does:

- a malformed file;
- an unknown key, at any level. A misspelled key would otherwise silently change a model.
  Keys added later come with a new format version;
- an unknown component in a `condition`;
- a missing profile;
- a duplicate model;
- an unknown `superiors` name. A bare name in `superiors` refers to a model in the same
  collection, because model names are unique only per `(database, name)`. `database:name`
  refers to a model in another collection.

Why only the native format: the converters for foreign formats are the part of the design most
likely to change with upstream releases, and the part where translation is approximate.
Running them once, at conversion time, means:

- approximations are reported once, in a conversion summary, not as warnings on every run;
- what runs is visible: thresholds, windows and cutoffs are written into a file the user can
  read, diff and edit;
- every command's tests and options deal with one format, and format options such as
  `antismash_strictness` belong to the converter alone;
- a converted directory pins one release of the source, so analyses are reproducible.

The cost is one extra step for users who already have a copy of a foreign format. The error
message above makes that step obvious.

## The native format

YAML. PyYAML is already installed (jsonargparse requires it) and would become an explicit
dependency. PADLOC already uses YAML, and since JSON is valid YAML, programs can emit JSON.
TOML (`tomllib`) was considered, but it has no stdlib writer and is clumsy for lists of
mappings.

```yaml
domainator_systems: 1            # format version, "major" or "major.minor"; required, refused if newer (FormatVersionError)
database: my_systems             # optional; default is the file stem
semantics: domainator            # optional; converters write padloc, macsyfinder or antismash
source:                          # optional; written by system_db_download, informational
  name: padloc
  url: https://github.com/padlocbio/padloc-db
  release: v2.0.0
  date: 2026-11-09
  converter: domainator 1.1.0
reference_files:                 # optional; profile files, relative to this file
  - profiles/defense.hmm
defaults:                        # optional; apply to every model / ReferenceRule unless overridden
  evalue: 1.0e-5
  min_reference_coverage: 0.4
  window: {mode: chain, max_cds: 3}
competing:                       # optional; overlapping hits within a group compete, best score wins
  - [Gabija__GajB_1, Gabija__GajB_2, Gabija__GajB_3]
competing_overlap_aa: 20

systems:
  - name: Gabija
    category: defense
    description: Gabija anti-phage system
    window: {mode: chain, max_cds: 5}
    min_core: 2                  # default when no condition: number of core components
    min_total: 2                 # default: min_core
    components:
      - name: GajA
        role: core
        domains: [Gabija__GajA]
      - name: GajB
        role: core
        domains:                 # any one of these satisfies GajB
          - Gabija__GajB_1
          - {name: Gabija__GajB_2, min_score: 40}
          - {name: Gabija__GajB_3, database: defense, evalue: 1.0e-10}
    examples:
      - {accession: NZ_CP000000.1, start: 120110, end: 125930, note: B. cereus VD045}

  - name: RelE_HicB_like_TA
    category: toxin_antitoxin
    window: {mode: span, max_kb: 2.5}      # the whole system within 2.5 kb
    same_strand: true
    components:
      - {name: toxin, role: core, domains: [RelE, ParE_toxin]}
      - {name: antitoxin, role: core, domains: [HicB, PhdYeFM_antitox]}
      - {name: transposase, role: forbidden, domains: [DDE_Tnp_1, Transposase_11]}
    orientation_from: antitoxin
    superiors: [Gabija]          # drop this call where a Gabija system overlaps it
    neighborhood_kb: 5           # system_search --model_neighborhood pads by 5 kb

  - name: lanthipeptide_like
    category: RiPP
    window: {mode: anchor, max_kb: 20}
    neighborhood_kb: 20
    components:
      - {name: LanB, role: core, domains: [Lant_dehydr_N, Lant_dehydr_C]}
      - {name: LanC, role: core, domains: [LANC_like]}
      - {name: LanM, role: core, domains: [DUF4135]}
      - {name: transporter, role: neutral, domains: [ABC_tran]}
    condition: "(LanB & LanC) | LanM"
```

`extenders` (for example `extenders: "transporter"` on the last model) arrives with format
1.1, with the condition functions and `define`. A file that uses them must say
`domainator_systems: "1.1"`, and a reader that supports only 1 refuses it with an error
naming the version it needs.

Rules:

- A domain entry is either a bare reference name, or a mapping with `name` plus optional
  `database`, threshold keys (`evalue`, `min_score`, `min_reference_coverage`,
  `min_target_coverage`, `min_target_length`, `max_target_length`) and `attach_to_cds`.
  Keys use the CLI vocabulary (`evalue`, `min_score`), as config files keyed by option dest
  do.
- Every key of a model is optional except its name and a requirement. Window, thresholds and
  `same_strand` can all come from `defaults`.
- Names are matched exactly. Wildcards or regex are a possible later addition; exact
  matches keep the lookup a dict hit (see [detection.md](detection.md#cost)).
- `reference_files` lets `domainate --system_models my_systems.yaml` find the profiles to
  search without `-r`.
- `examples` are known occurrences. `format_systems --check_examples` (a `network`-marked
  test in CI) downloads each one, runs detection and reports models that no longer find
  their own examples. antiSMASH's `EXAMPLE` lines serve the same purpose.
- A writer (`write_system_models`) emits this format from any `ModelCollection`. It powers
  `system_db_download` and `format_systems` (see [tool_integration.md](tool_integration.md)),
  and the tests use it for round trips. A user forks a converted DefenseFinder, PADLOC or
  antiSMASH model by copying it into their own file. The writer uses the short form for
  models that fit it.
- The format is a file-format contract once released: it is documented in
  `docs/file_formats.md`, versioned by the `domainator_systems:` key, and the current and
  previous major versions are read, as for `.dvb`.
- The version is `major` or `major.minor`. Readers parse it into a pair of integers, never a
  float, because YAML reads `1.10` as the number 1.1. The writer writes a bare `1` for 1.0
  and quotes minor versions (`"1.1"`).

### The short form

Most hand-written models are "these domains, in this combination, within N CDSs". The full
form makes that verbose: every domain needs a component wrapper and a role that doesn't
change the call. Two rules make it short:

1. **Implicit components.** A model with no `components:` list takes its components from its
   `condition`: every name in the condition is a reference name and becomes a component of
   the same name, with that one reference. `db:name` limits a name to one database, as in
   `--systems`.
2. **The mapping shorthand.** `systems:` may be a mapping from model name to condition,
   instead of a list of models. Every other key then comes from `defaults`.

"10-CDS windows with X, Y, Z and Q or R, at e-value 1e-20" is then:

```yaml
domainator_systems: 1
defaults: {evalue: 1.0e-20, window: {max_cds: 10}}
systems:
  XYZ_QR: "X & Y & Z & (Q | R)"
  XYW: "X & Y & W & ~transposase"
  HicB_RelE: "HicB < RelE"          # antitoxin upstream of toxin
```

The same model in the full form, for comparison:

```yaml
domainator_systems: 1
defaults: {evalue: 1.0e-20}
systems:
  - name: XYZ_QR
    window: {mode: span, max_cds: 10}
    components:
      - {name: X, role: core, domains: [X]}
      - {name: Y, role: core, domains: [Y]}
      - {name: Z, role: core, domains: [Z]}
      - {name: Q, role: core, domains: [Q]}
      - {name: R, role: core, domains: [R]}
    condition: "X & Y & Z & (Q | R)"
```

What the short form implies:

- There is no quorum. The condition is the whole requirement, and `completeness` is `.`, as
  for antiSMASH rules.
- A name used positively anywhere in the condition is a `core` component. A name used only
  under `~` (`transposase` in `XYW`) is a `forbidden` component: it is never a member and
  never counted, and it vetoes the system as described under [Windows](#windows).
- Members are every CDS in the window that hits a name used positively in the condition.
- A CDS that carries both X and Y satisfies both, because conditions are evaluated over
  every hit in the window (see [detection.md](detection.md#4-assign-components-per-member-per-model)).
  So names really are domains, and fusions work.
- `~W` means "no W in this window", not "no W on the contig".

A `systems:` block is either a list or a mapping, but not both. A file can mix short and full
models by using the list form, where a short entry is `{name: XYW, condition: "..."}` with no
`components`.

The command-line options `--domain_expr` and `--window_cds`
([tool_integration.md](tool_integration.md#command-line-models)) build exactly this model, so
the command line, the short form and the full form share one semantics. `format_systems`
writes a command-line model out as a short-form file.

### Conditions

The quorum (`min_core`, `min_total`) covers PADLOC fully and most of MacSyFinder.
`condition` covers antiSMASH and anything else. **The condition grammar adopts antiSMASH's
semantics**, a design that has been tested on more than 100 published rules, spelled with
Domainator's operators so it reads like `--domain_expr`:

| Construct | Meaning | antiSMASH spelling |
|---|---|---|
| `a`, `a & b`, `a \| b`, `~a`, `( )` | component present / both / either / absent; precedence `<` > `~` > `&` > `\|` | `a`, `and`, `or`, `not` |
| `a < b`, `a < b < c` | a upstream of b (of c) in system orientation; see [Order](#order) | none (Domainator addition) |
| `cds(a & b)` | a and b on **one** CDS (a fusion). `cds(a & ~b)` is a CDS with a but not b. | `cds(a and b)` |
| `minimum(n, [a, b, c])` | at least n **distinct** listed components present | same |
| `repeatable(n, [a, b, c])` | at least n hits to listed components, counting each (CDS, component) pair, so repeats count | same |
| `minscore(a, 150)` | a present with a bit score of at least 150 | same |
| `DEFINE` aliases | named sub-expressions | `DEFINE X AS ...` |

Aliases are written as a `define:` mapping at collection level:
`define: {ANY_KS: "PKS_KS | ene_KS | mod_KS"}`, used as `cds(PKS_AT & ANY_KS)`.

The condition is evaluated over **the members inside one window**: the chain cluster, the
anchor's radius, or the span window. Negation therefore means "absent from this window".
That is antiSMASH's meaning of `not` around an anchor, and it generalizes the `forbidden`
role, which is simply `~x` with x listed as a component.

`BooleanEvaluator` evaluates over a flat name set, so conditions get their own evaluator,
in `domainator.systems.conditions`. It shares `expr`'s tokenizer, so the two read names the
same way. A condition is parsed once into an AST and evaluated over a window context:
per-CDS component sets, scores and positions. Plain boolean expressions mean the same thing
in both.

Evaluation also returns the components that made the condition true. Anchor mode needs
this, because, as in antiSMASH, an anchor counts only if it contributed (see
[detection.md](detection.md#anchor-antismash)).

#### Order

`<` is a Domainator addition beyond antiSMASH. It covers operon layouts ("the antitoxin
upstream of the toxin") and, inside `cds()`, domain architectures in the style of InterPro
IDA, which the `select_by_contig` TODO asks for.

- **Between members:** `a < b` is true when some member filling `a` lies upstream of some
  member filling `b`, in the window's orientation. Members are compared by CDS position, so
  a and b on the same CDS do not satisfy `a < b` outside `cds()`.
- **Inside `cds()`:** `cds(a < b)` is true when one CDS has a hit to a that starts before a
  hit to b in protein coordinates, N- to C-terminal, whatever the CDS's strand.
  `cds(a < b < c)` is a three-domain architecture. Hits are compared by their start; ties
  are not ordered.
- **Chains:** `a < b < c` needs three distinct members (or hits) in that order, not three
  separate pairs.
- **Operands** are component names, or `define:` aliases that are disjunctions of names
  (`ANY_KS < PKS_AT`). An operand can't contain `&`, `~` or another `<`.
- **Precedence:** `<` binds tightest, so `~a < b` is `~(a < b)`, "no a upstream of b in this
  window", and `a < b & c` is `(a < b) & c`.
- **Orientation:** the window's orientation is computed before the condition is evaluated,
  with the same rule as a system's strand ([data_model.md](data_model.md#location)): the
  strand of the CDS members, by majority, or the `orientation_from` component's strand.
  Windows on circular contigs are compared in unwrapped coordinates.
- **Membership:** `<` doesn't change which CDSs are members. Its names count as positive
  terms, so the short form makes them components.
- **Names:** reference and component names can't contain `<`. That is checked when models
  are read, as for `|`.

Only `<` is offered: `b > a` is the same as `a < b`. An "immediately upstream" operator is a
possible later addition.

**Plan:** the full grammar is fixed now, in the format spec, and implemented in two steps.

1. Plain booleans and `<` between members, in `domainator_systems: 1`.
2. The functions (`cds`, `minimum`, `repeatable`, `minscore`), `define`, and therefore
   `cds(a < b)`, in the release that adds the antiSMASH converter, as
   `domainator_systems: 1.1`.

Readers refuse unknown functions with an error that names the needed format version.

## Nucleotide references

A component is satisfied by a hit to a reference name, and detection doesn't care which
program made the hit. So every reference type Domainator searches can be a component:

| Reference | Hit lands on | How detection sees it |
|---|---|---|
| Protein HMM (hmmsearch), protein sequence (phmmer), Foldseek structure | a CDS (`cds_id` = the CDS) | a member CDS |
| DNA HMM or DNA sequence (nhmmer), covariance model (Infernal cmsearch) | the contig (`cds_id="."`) | a pseudo-CDS member (`DomainatorCDS.list_from_contig(..., include_nucleic_acid_annotations=True)`), which doesn't count as a CDS in windows |
| `find_features` motifs | as written by `find_features` | like any other hit |
| CRISPR arrays and other non-hit features | (no `Domainator` feature) | not supported; a later "feature component" could match a feature type and qualifier |

Rules for nucleotide hits:

- **A nucleotide hit inside a CDS** stays a separate pseudo-CDS member by default. That is
  right for ncRNAs in untranslated regions and riboswitches that overlap a gene start. A
  reference with `attach_to_cds: true` assigns its hits to a CDS they overlap on the same
  strand by at least half the hit's length. The hit then counts as that CDS in `max_cds`
  windows and in `cds()`, and the member's `cds_id` is the CDS's. This is for DNA profiles of
  coding genes, for example on contigs that were searched without gene calls and annotated
  with genes later.
- **Records without CDSs** can only be measured in kb (see [Windows](#windows)).
- **Orientation** is voted by CDS members only, because ncRNA strand is often uninformative
  (see [data_model.md](data_model.md#location)). `orientation_from` can name an ncRNA
  component when its strand should decide.
- **Scores** from hmmsearch and cmsearch are summed into one system score although they
  aren't on the same scale. That is documented, not corrected: the score ranks calls of one
  model, and every call of a model mixes the same kinds of hits.
- **Inputs that can't match**, such as protein FASTA searched with a model that needs an
  ncRNA, warn once per model, not once per record.

## Converting foreign formats

`system_db_download` ([tool_integration.md](tool_integration.md#system_db_download-category-database))
either downloads a pinned release of a foreign format, or reads a copy the user already has,
and writes a native directory:

```
padloc_db/
  systems.yaml     # native format, with database: padloc, semantics:, source: and reference_files: [padloc.hmm]
  padloc.hmm       # every profile the models use, named as systems.yaml names them
  LICENSE          # copied from the source
  conversion.tsv   # every approximation made: model, feature, what was done
```

It recognizes these inputs:

| `--source` | Download | Local `--input` |
|---|---|---|
| `padloc` | a padloc-db GitHub release | a directory with `sys/` and `hmm_meta.txt` |
| `defensefinder` | a defense-finder-models GitHub release | a MacSyFinder package directory (`metadata.yml` and `definitions/`), e.g. one installed by `defense-finder update` |
| `macsyfinder` | (none; packages come from many places) | any MacSyFinder package directory |
| `antismash` | (none; read from an installation, see [licensing](#licensing)) | an antiSMASH `detection/hmm_detection` directory, or a rule `.txt` file with a `hmmdetails.txt` |

What conversion does beyond translating models:

- it collects every profile the models use into `<database>.hmm`, and resolves name clashes
  between collections once, by keeping each collection's `database`. The file is named after
  the collection's database, not `profiles.hmm`, because a hit's `database` qualifier comes
  from the reference file name (`seqio.db_name_from_path`). Then `domainate -r
  padloc_db/padloc.hmm` followed by `find_systems`, and `domainate --system_models
  padloc_db/`, label the hits the same way and find the same systems;
- it writes every threshold explicitly: PADLOC's `hmm_meta.txt` values with their NA
  defaults, MacSyFinder's defaults and GA cutoffs, antiSMASH cutoffs and kb conversions;
- it records the source, release and converter version in `source:`. Commands warn when a
  directory was made by a converter older than the running Domainator's, because a later
  converter may translate more faithfully.

The sections below give each mapping. In each, "unsupported" means the converter drops or
approximates the feature, and lists it in `conversion.tsv` and the summary it prints.

### PADLOC (`padloc-db`)

padloc-db layout: `sys/*.yaml` (about 385 models), `sys_meta.txt`, `hmm_meta.txt` (about
5,000 HMM rows), `hmm/`, and `cm/` + `cm_meta.txt` (retron ncRNA covariance models). It is
MIT licensed.

| PADLOC | `SystemModel` |
|---|---|
| model file `sys/<x>.yaml`, `sys_meta.txt` row with `search == TRUE` | one model; `name` = `sys_meta.system` (or the yaml stem); `accession` = `yaml.name`; `database` = `padloc`; `program` = `padloc`; `category` = `defense` |
| `maximum_separation` | `window = Window("chain", max_cds=maximum_separation)` (intervening genes; PADLOC starts a new cluster when the position difference is greater than `max_sep + 1`) |
| `minimum_core`, `minimum_total` | `min_core`, `min_total` (counted as **distinct** protein names) |
| `force_strand: TRUE` | `same_strand=True` (PADLOC clusters each strand separately, which is equivalent) |
| `core_genes` / `secondary_genes` / `neutral_genes` / `prohibited_genes` | components with roles core / accessory / neutral / forbidden. Overlaps are removed in PADLOC's order: secondary minus (core ∪ prohibited), prohibited minus core, neutral minus the rest. `NA` means an empty list. |
| a protein name, e.g. `GajB` | one component; its references are every `hmm_meta` row whose pipe-split `protein.name` includes `GajB` (e.g. `GajB_00001` … `GajB_00010`) |
| a group name in a gene list, e.g. `cas_adaptation` | expanded to **one component per protein name** in the group (rows whose pipe-split `secondary.name` includes it), all with the list's role. Each protein name counts separately toward `minimum_total`, as in PADLOC. |
| `hmm_meta` `e.val.threshold`, `hmm.coverage.threshold`, `target.coverage.threshold` (NA → 1e-5, 0.3, 0.3) | `ReferenceRule.evalue`, `min_reference_coverage`, `min_target_coverage` |
| `hmm/<hmm.name>.hmm` (or the combined file) | `reference_files` |
| `*_other` catch-all systems | `superiors` = every canonical (non-`_other`) model sharing a component with it. This replaces the earlier `yield_to_overlapping` flag with antiSMASH's more general mechanism. |
| `CRISPR_array` | **Unsupported.** Domainator has no CRISPR array detector. The component is dropped and listed in `conversion.tsv`, and `min_total` is lowered if it can no longer be reached. A later option: let a component match a non-Domainator feature, such as a `repeat_region` with `/rpt_family="CRISPR"` from an external tool. |
| `ncRNA` (Infernal CMs from `cm/`) | Supported as contig-level nucleotide members: Domainator already runs Infernal CMs and writes `cds_id="."` hits. Like PADLOC, these don't count as genes in gap counting. |

Notes:

- `hmm_meta.txt` is read by column **position**, as padloc does, because the README and the
  file header spell the e-value column differently.
- PADLOC computes coverage as `(to - from) / len`, without the `+1`. The converter can store a
  `coverage_offset` so imported models reproduce PADLOC exactly, or accept the one-residue
  difference and document it. Recommended: document it.

### MacSyFinder / DefenseFinder packages

Package layout: `metadata.yml`, `definitions/[subfamily/...]/<Model>.xml`,
`profiles/<gene>.hmm[.gz]` (one profile per gene, file name = gene name), and an optional
`model_conf.xml`. defense-finder-models has about 666 models under
`definitions/{DefenseFinder,RM,ADF}/<Family>/`.

| MacSyFinder | `SystemModel` |
|---|---|
| model `definitions/DefenseFinder/Gabija/Gabija.xml` | `name` = file stem (`Gabija`); `accession` = fully qualified name (`DefenseFinder/Gabija/Gabija`); `database` = package name; `program` = `macsyfinder`; `category` = `defense`. If two stems collide in one package, the converter raises an error suggesting `--name_style qualified`. |
| `inter_gene_max_space` (required; counts genes) | `window = Window("chain", max_cds=inter_gene_max_space)` |
| per-gene `inter_gene_max_space` | `Component.max_gap_cds`. Between two members, use the model value if neither overrides it, the one override if only one does, and **`min(d1, d2)`** if both do. This is what `cluster.py::_colocates` does, not what the docs say. |
| `min_mandatory_genes_required` (default: number of mandatory genes) | `min_core` |
| `min_genes_required` (default: number of mandatory genes) | `min_total` |
| `max_nb_genes` | `max_components` |
| `<gene presence=...>` | `Component.role` |
| `<exchangeables><gene name=.../>` | extra `ReferenceRule`s on the parent component, with `exchangeable=True`. They inherit the parent's role and flags. |
| `loner`, `multi_system`, `multi_model` | the matching `Component` flags. `loner` is kept, but evaluated as an ordinary component until the [multi-locus add-on](multi_locus.md); `conversion.tsv` lists it. |
| `multi_loci="1"` | `multi_loci=True`, evaluated as single-locus until the [multi-locus add-on](multi_locus.md); `conversion.tsv` lists it. |
| gene name → `profiles/<name>.hmm` | `ReferenceRule.name` = the profile's `NAME` line (read from the file, usually equal to the gene name); the file goes into `reference_files` |
| GA cutoffs (`cut_ga`, default on) | `ReferenceRule.min_score` = the profile's domain GA. This approximates `hmmsearch --cut_ga`, which also applies the sequence GA. A later `domainate --bit_cutoffs {gathering,trusted,noise}` (pyhmmer's `bit_cutoffs`) would make it exact. |
| `i_evalue_sel` (0.001), `coverage_profile` (0.5; DefenseFinder passes 0.4) | `ReferenceRule.evalue`, `min_reference_coverage`. The converter takes a `--coverage_profile` option defaulting to DefenseFinder's 0.4, and honors `model_conf.xml`. |
| `vers="1.0"` / `<system>` root (MacSyFinder v1) | refused with an `InputError` that names the version |

### antiSMASH cluster rules

Covered in detail in [antismash.md](antismash.md#converting-antismash-rules). In short:

- every profile named in a rule becomes a component;
- `CONDITIONS` becomes `condition`;
- `CUTOFF` becomes `Window("anchor", max_kb=...)`;
- `NEIGHBOURHOOD` becomes `neighborhood_kb`;
- `EXTENDERS`, `SUPERIORS` and `CATEGORY` map one-to-one;
- `hmmdetails.txt` cutoffs become `min_score`;
- `filterhmmdetails.txt` becomes the collection's `competing` groups.

## Thresholds are applied at detection time

Every per-reference threshold above can be evaluated from qualifiers that domainate already
writes:

- e-value and score from `evalue` and `score` (`model.feature_evalue`, `feature_score`);
- reference coverage from `rstart`, `rend` and `rlen`;
- target coverage and length from the hit length (its location) and the CDS length.

So detection needs no new search and is reproducible from a GenBank file. The catch is that
the search must not have discarded the hits a model needs:

- run domainate with an `-e` at least as permissive as the models' loosest threshold;
- leave `--max_hits_per_cds` unset and `--max_overlap 1`, because exchangeable profiles
  (`GajB_1`, `GajB_2`) and competing groups (antiSMASH's KS profiles) hit the same protein,
  and the system logic, not the search, must choose between them.

`domainate --system_models` takes care of this for the profiles it loads from the models: it
searches each one with the larger of `-e` and that profile's loosest model threshold, and
exempts those hits from `--max_hits_per_cds` and `--max_overlap`. `-r` references keep the
usual options, so the CLI defaults don't change.

Searching with the looser value would make `-e` unable to tighten a model profile. So in
`domainate --system_models` and `system_search`, as in `find_systems`, `-e` is also an evidence
filter in detection: a hit counts only if it passes both `-e` and the model's own threshold.
Then the one-step and two-step routes find the same systems for the same `-e`. `find_systems` warns when the input's hits
look filtered (for example, a model references an HMM that never appears in a file that has
other hits from the same database).

## Licensing

| Source | License | Recommendation |
|---|---|---|
| padloc-db | MIT ("Copyright (c) 2020 Leighton Payne") | `system_db_download` downloads it. Small excerpts may be vendored as test fixtures with the license notice. |
| defense-finder-models | Ambiguous: `LICENSE` is GPL-3.0, but `metadata.yml` has two `license:` keys, GPL-3.0 and CC BY-NC-SA 4.0 | Don't vendor. Downloading to the user's machine is fine (decided). `system_db_download` downloads it to the user's machine, or converts an installed copy (`macsydata install defense-finder-models`, `defense-finder update`, or a git clone). Test fixtures are synthetic packages written for Domainator. |
| CasFinder | CC BY-NC-SA 4.0 | Same as defense-finder-models. |
| MacSyFinder/MacSyLib code | GPL-3.0 | Don't import or copy. The converter is written from the documented format and checked against the observed behavior noted above. |
| antiSMASH code, rule files and bundled profiles | AGPL-3.0 | Read from the user's antiSMASH installation or checkout; never vendor. Write the parser from the documented grammar (the module docstring of `rule_parser.py`), not by porting code. Test fixtures are synthetic rules. See [antismash.md](antismash.md#licensing). |

`system_db_download` fits all of these. It downloads to the user's machine, copies the
source's license next to the converted files, and never ships model content with
Domainator. A converted antiSMASH directory is derived from AGPL files, so it is for the
user's own use; the command says so when it converts antiSMASH rules.
