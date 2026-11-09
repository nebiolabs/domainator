# Systems: detection algorithm

Detection is a **pure function of one annotated contig and a list of model collections**:

```python
def detect_systems(record: SeqRecord, collections: Sequence[ModelCollection], *,
                   evalue: Optional[float] = None,
                   databases: Optional[Collection[str]] = None,
                   annotation_type: str = "domainator",
                   partial: str = "include",
                   max_system_overlap: float = 1.0,
                   max_systems_per_contig: Optional[int] = None) -> List[SystemHit]: ...
```

It reads `Domainator` features (and optionally `Domain_Search` features, through the usual
`--annotation_type`) and returns `SystemHit` objects. `evalue` and `databases` limit which
hits count as evidence, and `partial` which CDSs can be members. `partial_systems` joins
the signature in Phase 4, with `--partial_systems include`. A separate function turns them into
`Domain_System` features. It never searches. `find_systems`, `domainate --system_models` and
`system_search` all call it, so a system found by one is found by the others. A model
defined on the command line (`--domain_expr` with `--window_cds`/`--window_kb`) is an
ordinary short-form model, so it goes through the same steps.

## Steps

### 1. Order the contig

Use `DomainatorCDS.list_from_contig(record, include_nucleic_acid_annotations=True)`, which
already groups hits by `cds_id` and turns contig-level nucleotide hits into pseudo-CDSs.

- Sort by position.
- Give each real CDS an integer **gene index**. Give pseudo-CDSs (ncRNA hits) a fractional
  index between their neighbors, so they never count as CDSs in any window. This matches
  PADLOC, and the `cds_range_members` helper in `select_by_cds.py` already leaves
  pseudo-CDSs out of CDS counts.
- A nucleotide hit whose reference sets `attach_to_cds` is moved onto the CDS it overlaps
  on the same strand by at least half the hit's length, before indexing. It then counts as
  that CDS. With no such CDS it stays a pseudo-CDS (see
  [definition_formats.md](definition_formats.md#nucleotide-references)).
- A record with no real CDSs has no gene indices. Models whose window sets only `max_cds`
  are skipped on it, with one warning per model for the run.
- Pseudo genes (`skip_pseudo`) and partial CDSs follow the caller's `--partial` choice, as
  in other tools.

### 2. Resolve competing references

For each collection with `competing` groups, and within each CDS: when two hits whose
reference names are in the same group overlap by more than `competing_overlap_aa` residues,
keep only the best-scoring one. Hits within 1e-8 bits of the best are all kept, as antiSMASH
does. Hits outside every group are untouched.

This is antiSMASH's equivalence-group filter (`filterhmmdetails.txt`). It is how a KS domain
counts as `t2ks` *or* `PKS_KS` but never both. The same filter is a cleaner answer for
families of exchangeable profiles, such as `GajB_1`/`GajB_2`/`GajB_3`.

The filter works on a view of the hits for detection only. The `Domainator` features in the
record are not changed.

### 3. Match hits to components

Before detection, build a lookup `reference name → [(model, component, rule)]` over all
models, so each hit costs one dict lookup.

For each surviving hit on each CDS:

1. look up the hit's `name`;
2. keep the matches whose `rule.database` is unset or equal to the hit's `database`;
3. apply the rule's thresholds (`model.evalue_passes`, score, reference coverage, target
   coverage, target length).

The result is a list of candidate `(cds, model, component, hit)` tuples. Only models with at
least one positive candidate on the contig go on. With 600+ DefenseFinder models or 100+
antiSMASH rules, that is usually a handful per contig.

### 4. Assign components per member, per model

A CDS can match several components of one model, for example an HMM shared between `LmuB`
and `LmuC` variants. Within one model, each CDS takes one **primary** component, chosen by:

1. role priority: core, then accessory, then neutral, then forbidden (PADLOC's order);
2. then the higher bit score;
3. then the lower e-value.

The primary component is what counts toward the quorum and score. The CDS also keeps the
full set of components it matched, because `cds(a & b)` conditions ask about everything on
one CDS (fusions).

A CDS may still be a member of *different* models; step 8 resolves that.

### 5. Form candidate windows

The model's `window.mode` decides how members are grouped. Each mode yields **candidate
windows** (sets of members) for step 6 to test. See
[definition_formats.md](definition_formats.md#windows) for the definitions and a worked
example.

#### `chain` (MacSyFinder, PADLOC)

Walk the assigned members in gene order. Chain consecutive members into one cluster while
the gap is within `max_cds` (CDSs strictly between them) and `max_kb` (bp from the end of
one to the start of the next).

- Per-component `max_gap_cds` overrides follow MacSyFinder: use the override, or the `min`
  of two overrides.
- **Every role chains, forbidden included**, as in PADLOC. A forbidden gene sitting between
  two core genes still joins them, and then vetoes the cluster in step 6.
- Each maximal cluster is one candidate window.

#### `anchor` (antiSMASH)

1. Every member that fills a positive component is a possible **anchor**.
2. Its candidate window is every member within `max_cds` intervening CDSs and `max_kb` of
   it (whichever are set), measured edge to edge.
3. An anchor whose window passes step 6 counts only if it **contributed**: one of the
   components it fills is among those that made the condition true. This is antiSMASH's
   rule. For example, with `(A & B) | C` and a window holding an A anchor and a C member,
   C alone satisfies the condition, so the A anchor doesn't count (the C member, as an
   anchor, does). Contributing anchors mark themselves, and the members that helped satisfy
   the requirement, as **core members** of the model.
4. Core members are then joined: consecutive core members within the same radius of each
   other form one system, as antiSMASH joins neighboring cores. So a system can be longer
   than one window, but every core member was justified by some anchor.

antiSMASH tests distance with a strict `<`. Domainator uses `<=` everywhere, so the
antiSMASH reader converts `CUTOFF c` to `max_kb = (1000·c − 1) / 1000`. The two are then
identical on integer bp coordinates.

#### `span` (Domainator "the whole system within x")

A candidate is any contiguous run of members, in gene order, that covers at most `max_cds`
CDSs (first to last inclusive) and spans at most `max_kb` kb, whichever are set.

- Positive terms of the requirement are tested on the run's own members.
- Negated terms and `forbidden` components are tested on the run's **reach**: every member
  that could be added to the run with the result still fitting the limit. A maximal window
  around the run never extends past its reach, so "no W within reach" means "no W in any
  window of the allowed size that contains the system" (decided).
- A passing run is a candidate only if no shorter run inside it also passes. So the system
  is the tightest set of members that satisfies the model, and stray duplicates within the
  limit don't stretch it.

Then pick non-overlapping candidates greedily by score, breaking ties toward the smaller
span, then the leftmost.

Testing every run is O(H × w²) for H members and windows of up to w members, which is small
in practice. Testing only maximal windows would be cheaper, but it gets negation wrong.
With X, Y and W at CDSs 1, 2 and 9 and `X & Y & ~W` in a 5-CDS span, every maximal window
containing X and Y also contains W, so a valid 2-CDS system would be missed.

#### Constraints common to every mode

- **What negation sees:** the cluster in `chain` mode, the anchor's radius in `anchor` mode,
  the reach in `span` mode (see [definition_formats.md](definition_formats.md#windows)).
- **`same_strand`:** form windows on each strand separately.
- **`max_span_cds` / `max_span_kb` caps** (chain and anchor only): a system longer than the
  cap is split with the `span` procedure inside the over-long cluster.
- **Circular contigs:** handled in every mode by extending the member list past the origin
  (members repeated with positions shifted by the contig length) and folding results back.
  `locations.circular_dist` gives wrapped distances. A system that wraps gets a two-part
  location (see [data_model.md](data_model.md#location)), as antiSMASH's
  `merge_over_origin` does.

### 6. Test the requirement

A candidate window is a system when all of these hold:

- distinct core components ≥ `min_core`;
- distinct core + accessory components ≥ `min_total`;
- there are no forbidden members in the window;
- `condition`, if any, is true when evaluated over the window. The window supplies each
  member CDS's matched components, scores and hit positions, which is what `cds()`,
  `minimum()`, `repeatable()`, `minscore()` and `<` need (see
  [definition_formats.md](definition_formats.md#conditions)). A negated term means "absent
  from this window". For `<`, the window's orientation is computed first, with the strand
  rule of [data_model.md](data_model.md#location) applied to the window's CDS members.

**Loners and multi-locus systems** (MacSyFinder `loner` and `multi_loci`) are a
[future add-on](multi_locus.md). Until then, a `loner` component is an ordinary component
that must be inside the window, and `multi_loci` models are evaluated as single-locus.
`system_db_download` lists both in its conversion summary. Every system location is
therefore one contiguous region, or two parts when it wraps the origin.

### 7. Trim and extend

**Trim.** Neutral members at either end of a system are trimmed. They may join members
together, but they don't make a system longer.

**Extend.** If the model has `extenders`, grow the system outward from each end, one CDS at a
time. A CDS joins when:

- it meets the extender condition, which is evaluated on that CDS alone, as antiSMASH does;
- it lies within the model's window radius (`max_cds`/`max_kb`) of the previous member in
  that direction.

Extended members are listed with role `extender`, and they are not trimmed. On circular
contigs, extension wraps.

### 8. Resolve overlapping systems

Several system hits can share members, both across models (Lamassu and a `*_other`
catch-all; T1PKS and transAT-PKS) and within one (tandem copies). Two mechanisms apply, in
this order:

1. **Superiors.** A system is dropped when a system of one of its model's `superiors` shares
   at least one CDS with it, or contains its envelope. A bare name in `superiors` means a
   model of the same collection; `database:name` reaches another collection. This is antiSMASH's
   `remove_redundant_protoclusters`, and PADLOC's `*_other` filter imports as superiors.
2. **Greedy overlap**, following the existing domain overlap filter
   (`domainate.filter_by_overlap`, `--max_overlap`):
   - sort the remaining hits by score, descending (ties: `database`, then `name`, so the
     result is deterministic);
   - accept greedily, rejecting a hit when the fraction of its counted members already used
     by accepted systems is greater than `--max_system_overlap`;
   - members whose component is `multi_model` (across models) or `multi_system` (within a
     model) are not counted as used.

`--max_system_overlap 1` (the default) keeps everything. That is PADLOC- and antiSMASH-like:
antiSMASH lets a CDS belong to several protoclusters, and later groups overlapping ones into
"candidate clusters". The default lets overlapping interpretations through to reports, as
domainate's default `--max_overlap 1` does for domains. `--max_system_overlap 0` gives
MacSyFinder-like exclusive calls.

`--max_systems_per_contig` caps the count after overlap filtering, by score, matching
`--max_hits_per_contig` for domains.

The greedy pass approximates MacSyFinder's best solution (a maximum-score clique of
compatible systems). On real genomes they rarely differ. The differences are noted under
[Fidelity](#fidelity).

### 9. Build the system hit

```python
@dataclass(frozen=True)
class SystemMember:
    cds: DomainatorCDS
    component: Component     # the primary component (step 4)
    role: str                # core | accessory | neutral | extender
    hit: SeqFeature          # the Domainator feature that satisfied it

@dataclass(frozen=True)
class SystemHit:
    model: SystemModel
    members: Tuple[SystemMember, ...]   # contig order
    location: SeqFeature location       # the core envelope (see data_model.md)
    score: float
    evalue: float
    completeness: Optional[float]       # None when the model has no quorum (condition-only models)
    core_count: int
    total_count: int
    partial: str                        # "", "<", ">", "<>"
```

The model's `neighborhood_kb`/`neighborhood_cds` is **not** added to the location. The
`Domain_System` location is always the core (antiSMASH's `proto_core`). The neighborhood is
recorded as a qualifier, and `system_search`/`select_by_system` apply it on request (see
[tool_integration.md](tool_integration.md)).

## Scoring

The system `score` is the **sum, over distinct counted components, of the best bit score of
a member filling that component**. Counted means core and accessory, or for condition-only
models, the components the condition matched positively. It has the same units and
direction as domain scores, so:

- `sort_contigs`, `dedup.TopHits` (through `model.set_hit_best_score`) and `--min_score`
  work unchanged;
- a more complete system scores higher than a partial one built from the same hits;
- duplicate copies of one component add nothing, so tandem paralogs don't inflate the score.

Bit scores from different programs (hmmsearch, cmsearch, Foldseek) are added together
although their scales differ. The score ranks calls of one model, and every call of a model
mixes the same kinds of hits, so this is documented rather than corrected.

Alternatives:

- **MacSyFinder's score:** 1 per mandatory, 0.5 per accessory, ×0.8 for exchangeables,
  −1.5 per redundant gene. Comparable across models, but blind to hit quality. It could be
  written as an extra `/msf_score` qualifier on models with `program="macsyfinder"`, if users
  need it to match DefenseFinder's tie-breaking.
- **Completeness alone:** already its own qualifier. It is a poor ranking key because many
  systems tie at 1.0.

## Contig ends and partial systems

On a linear contig, a window that ends within its radius of a contig end may be missing
components that lie past the end. That is common with metagenome-assembled genomes and with
neighborhoods that domain_search has already cut out. Two policies:

1. Report only systems that meet the requirement (the default), but set
   `/system_partial="<"`, `">"` or `"<>"` when a missing positive component *could* lie past
   an end, because that end is within the window's reach of the system. Readers then know that
   `completeness < 1` near an edge may be a cut, not an absence.
2. With `--partial_systems include`, also report windows that fail the requirement but touch
   an end and would pass if every missing positive component were beyond it. They get a
   `partial` flag and `completeness < 1`. This is useful for fragmented assemblies. The
   default is `exclude`.

The flags follow the `<`/`>` convention of `domainator.partial`, so reports can treat
partial systems like partial CDSs. The qualifier is `system_partial` rather than `partial`,
which is the INSDC qualifier Domainator reads as CDS fragment status.

## Cost

Each contig is visited once, after the search. For H passing hits on a contig and M models
with a positive candidate there:

- matching is O(H) dict lookups;
- `chain` grouping is O(H log H) per model;
- `anchor` and `span` are O(H × w), where w is the number of members inside one window
  (small in practice).

That is negligible next to hmmsearch. Detection runs on each record right after it is
annotated: in `domainate` that is in-process, after `domainator_inner`; in `system_search` it is
inside the pool task, as domain_search's selection is. So it needs no extra pass, and keeps
`-i -` streaming wherever the host tool supports it. In `system_search`, collections reach
the pool's workers once through `_mp.Shared(build=...)`, as references do.

## Fidelity

Converted models are evaluated with **the same semantics, and the differences documented**
(decided). They don't aim to reproduce every call of the original tools. Known differences:

| Area | DefenseFinder (MacSyFinder) | PADLOC | antiSMASH | Domainator |
|---|---|---|---|---|
| Search | `hmmsearch --cut_ga`, else `-E 0.1` | `hmmsearch`, then per-HMM thresholds | `hmmsearch`; domain bit score `>` per-profile cutoff | per-reference thresholds on domainate hits; `>=` instead of `>` for antiSMASH; GA approximated by domain GA |
| Competing hits | - | best rank-1 hit per protein | equivalence groups, more than 20 aa overlap | `competing` groups (step 2) |
| Coverage | `(to - from + 1) / len` | `(to - from) / len` | not used | `(rend - rstart + 1) / rlen`, so off by one residue from PADLOC |
| Distance | genes between consecutive members | genes between consecutive members | bp from an anchor (`<`); the anchor must contribute | `chain` / `anchor` (with contributing anchors) / `span` modes; antiSMASH `<` converted to `<=` |
| Gene positions | replicon order | GFF `relative.position` | CDS locations | CDS order in the GenBank record; ncRNA gets a fractional index |
| CRISPR arrays | (CasFinder only) | CRISPRDetect | - | not detected; components dropped with a warning |
| Code-based profiles | - | - | dynamic profiles (Python) | not importable; substitute `find_features` MOTIF hits ([antismash.md](antismash.md#dynamic-profiles)) |
| Overlapping systems | maximum-score clique | keep all, minus `*_other` | keep all, minus inferior rules; later merged into candidate clusters and regions | superiors, then greedy by `--max_system_overlap` |
| Loners / multi-locus | yes / yes | no / no | no / no | not yet ([future add-on](multi_locus.md)): loners must be in the window; multi-locus models are evaluated as single-locus |

A validation note, not a CI test: run DefenseFinder, PADLOC and antiSMASH (detection stage
only) on a few dozen complete genomes, run the imported models through `domainate
--system_models`, and compare the calls. Expect agreement on clear systems and a short,
explainable list of differences. `network`-marked tests can automate a smaller version of
this, and models' own `examples` give a per-model regression check.
