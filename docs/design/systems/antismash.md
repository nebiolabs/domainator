# Systems: antiSMASH cluster rules

antiSMASH detects biosynthetic gene clusters (BGCs) with a rule language over HMM profile
hits. It is the standard for BGC classification and the most mature rule language of the
formats in scope. This document covers:

- how antiSMASH's detection stage works;
- what the Domainator design borrows from it;
- how to convert antiSMASH rule files to the native format;
- how to import antiSMASH *results*.

Sources read (antiSMASH `master`, November 2026):

- `antismash/common/hmm_rule_parser/rule_parser.py`: the grammar, documented in the module
  docstring, and condition evaluation;
- `antismash/common/hmm_rule_parser/cluster_prediction.py`: hit filtering, protocluster
  formation, extenders and superiors;
- `antismash/detection/hmm_detection/`: `__init__.py`, `cluster_rules/{strict,relaxed,loose}.txt`,
  `data/hmmdetails.txt`, `filterhmmdetails.txt`, `dynamic_profiles/`.

## How antiSMASH's detection stage works

1. **Search and filter hits.**
   - `hmmsearch` runs every profile in `hmmdetails.txt` (about 400 profiles: name,
     description, bit-score cutoff, file). A hit is kept when its domain bit score is
     greater than the profile's cutoff.
   - Within one CDS, overlapping hits (more than 20 aa) from profiles in the same
     **equivalence group** (`filterhmmdetails.txt`, e.g. all ketosynthase profiles) compete,
     and only the best is kept.
   - Then only the best hit per profile per CDS is kept.
   - *Dynamic profiles* (Python code, e.g. a regex for darobactin precursors) add hits that
     have no HMM.
2. **Evaluate rules around anchors.** Each CDS with a hit is tried as an anchor for each
   rule. The rule's `CONDITIONS` are evaluated using the anchor and every CDS whose distance
   to it is less than `CUTOFF` kb. A rule is met only when the anchor itself contributes a
   matching profile. The anchor, plus the CDSs that contributed, are marked as core CDSs for
   that rule.
3. **Form protoclusters.**
   - Core CDSs of a rule are joined when they are within `CUTOFF` of each other, giving the
     *core*.
   - The core grows over adjacent CDSs that meet `EXTENDERS`.
   - The core is padded by `NEIGHBOURHOOD` kb, giving the *protocluster*.
   - Cores that span the origin are merged.
4. **Drop inferior rules.** A protocluster is dropped when a protocluster of one of its
   rule's `SUPERIORS` contains or overlaps its core.
5. **Later stages, which are not part of rule detection.**
   - Overlapping protoclusters are grouped into *candidate clusters* (kinds: `single`,
     `neighbouring`, `interleaved`, `chemical_hybrid`), and those into *regions*.
   - Downstream modules analyze the regions (NRPS/PKS domains, ClusterBlast, and so on).

The rule language, from the module docstring:

```
DEFINE ANY_KS AS (PKS_KS or ene_KS or mod_KS or hyb_KS or itr_KS or tra_KS)

RULE transAT-PKS
    CATEGORY PKS
    DESCRIPTION Trans-AT polyketide
    EXAMPLE NCBI CP006259.1 382565-499230 kirromycin
    REFERENCE DOI 10.1002/ange.201709056 toblerol with tra_KS and PP-binding
    CUTOFF 20            # kb
    NEIGHBOURHOOD 20     # kb
    CONDITIONS cds(PKS_AT and not ANY_KS)
               and (cds(ATd and ANY_KS) or cds(tra_KS and PP-binding))
    EXTENDERS cds(ANY_KS and not PKS_AT)
```

The condition constructs:

- `and`, `or`, `not`, with precedence `not` > `and` > `or`, and parentheses;
- `cds(...)`: the condition holds within a single CDS;
- `minimum(n, [...])`: n distinct profiles;
- `repeatable(n, [...])`: n hits, so repeats count;
- `minscore(p, s)`: p with a bit score of at least s;
- `DEFINE` aliases.

A rule must contain at least one positive term. `RELATED` lists profiles that are relevant
but not required, and `SUPERIORS` lists the rules that win when they overlap this one. There
are about 105 rules across the strict, relaxed and loose files. antiSMASH's
`--hmmdetection-strictness` loads them cumulatively: strict ⊂ relaxed ⊂ loose.

## What the Domainator design borrows

These changes are already folded into the other design documents.

| antiSMASH | In the Domainator design | Where |
|---|---|---|
| Distance measured from an anchor, not along a chain | `window.mode: anchor`, beside `chain` (MacSyFinder/PADLOC) and the Domainator-native `span` | [definition_formats.md](definition_formats.md#windows), [detection.md](detection.md#5-form-candidate-windows) |
| `cds()`, `minimum()`, `repeatable()`, `minscore()`, `DEFINE` | the condition grammar, spelled with `&`/`\|`/`~`; aliases in a `define:` mapping | [definition_formats.md](definition_formats.md#conditions) |
| Negation means "absent within the cutoff of the anchor" | negation means "absent from the window", in every mode; `forbidden` is `~x` | [detection.md](detection.md#6-test-the-requirement) |
| Equivalence groups (`filterhmmdetails.txt`) | collection-level `competing` groups, applied before matching | [detection.md](detection.md#2-resolve-competing-references) |
| `SUPERIORS` | `superiors` on the model; also replaces PADLOC's `*_other` flag | [detection.md](detection.md#8-resolve-overlapping-systems) |
| `EXTENDERS` | `extenders` condition; extended members get role `extender` | [detection.md](detection.md#7-trim-and-extend) |
| Core vs `NEIGHBOURHOOD` | `Domain_System` location is the core; `neighborhood_kb`/`neighborhood_cds` are recorded and applied on request by `--model_neighborhood` | [data_model.md](data_model.md#qualifiers), [tool_integration.md](tool_integration.md#shared-cli-vocabulary) |
| `CATEGORY` | `category` field, qualifier, `--system_categories` filter, coloring by category | [data_model.md](data_model.md#qualifiers), [tool_integration.md](tool_integration.md) |
| `EXAMPLE` lines | `examples` in the native format; `format_systems --check_examples` | [definition_formats.md](definition_formats.md#the-native-format) |
| Dynamic profiles (code-based hits) | `find_features` MOTIF hits as components, with `min_target_length`/`max_target_length` | [below](#dynamic-profiles) |
| Parse errors with line and column context | the same for the native reader and every converter | [definition_formats.md](definition_formats.md#reading-models) |
| Separate strict/relaxed/loose rule files | `system_db_download --antismash_strictness`; native users simply choose files | [below](#converting-antismash-rules) |

## Converting antiSMASH rules

`system_db_download --source antismash --input <path>` reads an antiSMASH
`detection/hmm_detection` directory, or a rule `.txt` file plus a `hmmdetails.txt`, and
writes one native collection with `database: antismash` and `semantics: antismash`. It copies
the profiles from the same installation into `profiles.hmm`. There is no download mode,
because antiSMASH is AGPL (see [Licensing](#licensing)).

```bash
domainator system_db_download --source antismash --input ~/antismash/antismash/detection/hmm_detection \
    --antismash_strictness relaxed --output_dir antismash_rules/
```

| antiSMASH | `SystemModel` / `ModelCollection` |
|---|---|
| `RULE name` | `name`; `accession` = `antismash:<name>` |
| `CATEGORY` | `category` |
| `DESCRIPTION`, `REFERENCE` | `description`, `literature` |
| `EXAMPLE NCBI acc start-end note` | `examples` |
| `CUTOFF c` (kb) | `window = Window("anchor", max_kb=(1000·c − 1)/1000)`, converting antiSMASH's `<` to Domainator's `<=` |
| `NEIGHBOURHOOD n` (kb) | `neighborhood_kb = n` |
| `CONDITIONS` | `condition`, translated `and`/`or`/`not` → `&`/`\|`/`~`. Functions keep their names. `DEFINE` aliases → the collection's `define:` mapping. Quorum `min_core = min_total = 0`. |
| a profile identifier in `CONDITIONS` | a component of the same name, role `core`, with one `ReferenceRule` |
| a profile identifier only in `EXTENDERS` | a component, role `neutral` |
| `EXTENDERS` | `extenders` |
| `SUPERIORS` | `superiors` |
| `RELATED` | kept as metadata (`format_systems --tsv_out` lists it). Not a component. |
| `hmmdetails.txt` row: name, description, cutoff, file | the `ReferenceRule` for that identifier: `name` = the profile's `NAME` line (antiSMASH also matches by accession, e.g. `PF01979.hmm` holding `Amidohydro_1`); `min_score` = cutoff, with `>` vs `>=` documented; the file goes into `reference_files` |
| `filterhmmdetails.txt` (one comma-separated group per line) | `competing` groups; `competing_overlap_aa = 20` |
| `--hmmdetection-strictness` | `--antismash_strictness {strict,relaxed,loose}`, loading rule files cumulatively as antiSMASH does. The default is `relaxed`, antiSMASH's default. |
| fungal cutoff and neighborhood multipliers (`--hmmdetection-fungal-neighbourhood-multiplier`, ...) | converter arguments `cutoff_multiplier` and `neighborhood_multiplier` (default 1). Exposing them on the CLI waits until users need it. |

### Dynamic profiles

antiSMASH has six code-based profiles (`dynamic_profiles/`):

- `atropopeptide_precursor`
- `cyanobactin_precursor`
- `darobactin_precursor`
- `darobactin_rsam`
- `triceptide_precursor`
- `triceptide_rsam`

They are Python functions, so they can't be read as data. Some are just a regex plus a
length range: `darobactin_precursor` is `W.W.K..` in proteins of 10-80 aa. Others also look
at HMMER hits.

- The converter turns each one into a component with **no references**, and lists it in
  `conversion.tsv` with the affected rules. Conditions that need it can never be met, and `~`
  terms that mention it are always true.
- A user can supply a substitute, a normal component whose reference is a `find_features`
  MOTIF hit:
  `{name: darobactin_W_W_K, database: MOTIF, min_target_length: 10, max_target_length: 80}`.
  Substitutes live in a small native file passed to the converter with `--substitutes_file`,
  which merges them by name. Keeping them out of the converted `systems.yaml` means they
  survive reconversion against a newer antiSMASH.
- The regex-only profiles could ship as ready-made native substitutes. That is optional, and
  it means re-expressing them from their description, not copying the code.

### What is not reproduced

- Candidate clusters and regions (antiSMASH's step 5), and every downstream module. A
  Domainator antiSMASH-rule run gives **protocluster cores**, one `Domain_System` per core.
  The regions antiSMASH would report can be approximated with `system_search
  --model_neighborhood --max_region_overlap 0` or by merging overlapping systems in a report.
  A `merge_systems` editor that groups overlapping systems into candidate-cluster-like
  `Domain_System` features (`/kind="chemical_hybrid"`, ...) is a possible later addition.
- antiSMASH's own CDS calling and record preprocessing. Domainator uses the CDSs in the
  input, or `--gene_call`.
- Exact hit sets. antiSMASH uses a domain bit score `>` the cutoff, and Domainator applies
  `>=`. Profile versions also differ unless the user points Domainator at the profiles in the
  same antiSMASH installation, which the converter does.

## Importing antiSMASH results

The cheapest way to put antiSMASH's BGC calls into Domainator workflows is to read
antiSMASH's **output** rather than re-run its rules. antiSMASH writes GenBank records with:

- `protocluster` and `proto_core` features (product, category, cutoff, `neighbourhood`, the
  detection rule);
- `cand_cluster` and `region` features;
- per-CDS qualifiers such as `/gene_kind` and `/sec_met_domain` (profile name, e-value and
  bit score).

Domainator keeps unknown features when it parses a file, so these records already pass
through every tool. They just aren't *systems*.

A small converter turns them into systems: a new `import_systems` command (category
`convert`), or `find_systems --source antismash`. **Recommended: `import_systems`**, because it
neither reads models nor detects anything.

| antiSMASH feature | Domainator |
|---|---|
| `proto_core` (with its `protocluster`) | one `Domain_System`. Location = the core. `name` = product. `category`, and `neighborhood_kb` from the protocluster. `database` = `antismash`. `program` = `antismash`. `score` = the sum of the core CDSs' best `sec_met_domain` bit scores. `component` values built from core CDSs and their `sec_met_domain` profiles, with role `core`. |
| `cand_cluster`, `region` | optional, with `--level {protocluster,candidate_cluster,region}`. One `Domain_System` per feature at the chosen level, with `/kind` from the candidate cluster. |
| the original antiSMASH features | kept by default, or removed with `--clear_source` |

Notes:

- The exact qualifier names should be checked against current antiSMASH output (7.x/8.x)
  before implementing. The fixture should be a small antiSMASH GenBank file produced by the
  maintainers from a public genome. antiSMASH output is data, not AGPL code.
- `sec_met_domain` gives profile names and scores but no positions within the protein. So
  imported systems reference profiles by name only; there are no matching `Domainator`
  features unless the user also runs domainate with the antiSMASH profiles.
- The same command is the natural home for importing **DefenseFinder** (`*_systems.tsv`) and
  **PADLOC** (`.csv`) results. Those files name proteins, so the converter matches them to
  CDSs by `locus_tag`/`protein_id` (`model.cds_name_precedence`). Users can then compare
  Domainator's own calls with the reference tools' calls on the same records.

## Licensing

antiSMASH's code, rule files and bundled profiles are **AGPL-3.0**. Domainator is MIT.

- **Don't vendor** rule files, `hmmdetails.txt`, `filterhmmdetails.txt` or the bundled HMMs.
  Read them from the user's antiSMASH installation or checkout, or from a path the user
  gives. (Some bundled profiles derive from Pfam or TIGRFAM, which have their own licenses,
  but the files as distributed come with antiSMASH.)
- **Don't port code.** Write the parser and evaluator from the documented grammar (the
  `rule_parser.py` module docstring) and from the behavior described here. The design
  borrows ideas, which are not covered by copyright, not code.
- **Test fixtures are synthetic:** hand-written rules and a tiny `hmmdetails.txt` that
  reference small Domainator-made HMMs. A `network`-marked test may clone antiSMASH at a
  pinned tag and check that the converter converts every shipped rule without errors.
- **Converted directories** contain antiSMASH rules and profiles in another form, so they
  carry antiSMASH's license. The converter copies antiSMASH's license file into the directory
  and says that it is for the user's own use; Domainator never distributes one.
- **Importing antiSMASH results** (`import_systems`) touches only antiSMASH's output files,
  which raises no license issue.

## Is it worth it?

Yes, in two stages.

1. **`import_systems --source antismash`** (Phase 2). Low cost and exact fidelity, because
   antiSMASH made the calls. Users get Domainator's selection, reporting, coloring, plotting
   and network tools for BGCs right away.
2. **The antiSMASH rule converter** (Phase 4). Medium cost: it needs the `anchor` window,
   the full condition grammar, competing groups, extenders and superiors. Most of that is
   already in the core design, because it improves the native format too, so the converter
   itself is small. The payoff is scanning very large genome databases for BGC neighborhoods
   (`system_search --system_models antismash_rules/ --system_categories PKS`)
   without running the full antiSMASH pipeline on every genome. The cost is that only
   protocluster cores are reproduced, as described above.
