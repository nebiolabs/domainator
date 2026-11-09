# Systems: implementation plan and roadmap

## Modules

A new internal package, like `taxonomy/`. It is not part of the public API, per
[`stability.md`](../../stability.md): only the command functions and root exports are public.

```
src/domainator/systems/
    __init__.py        # read_system_models (native only), detect_systems, window_model, re-exported for the commands
    model.py           # ReferenceRule, Component, Window, SystemModel, ModelCollection, SystemMember, SystemHit
    conditions.py      # condition grammar, including `<`: parse once to an AST; evaluate over a window context (shares expr's tokenizer)
    native.py          # read/write the domainator_systems YAML, full and short form; window_model() for command-line models
    formats/           # converters, used only by system_db_download (and by Python users)
        __init__.py    # convert(source, path, ...) -> ModelCollection; download helpers; conversion report
        padloc.py      # padloc-db converter
        macsyfinder.py # MacSyFinder package converter (xml.etree; no MacSyLib import)
        antismash.py   # antiSMASH rule converter (written from the documented grammar; no antiSMASH import)
    windows.py         # chain / anchor / span candidate windows, circular wrapping
    detect.py          # the pipeline in detection.md: competing → match → assign → windows → requirement → extend → overlap
    features.py        # SystemHit <-> Domain_System/System_Search SeqFeature; component qualifier codec
    records.py         # per-record accessors for selectors: systems_in(record), member_systems_by_cds_id(record)
    imports/
        antismash.py   # antiSMASH GenBank output → SystemHit (import_systems)
```

New commands, one module each, with the usual shape: `find_systems.py`, `system_search.py`,
`select_by_system.py`, `system_db_download.py`, `format_systems.py`, `import_systems.py`.

The condition grammar in `conditions.py` is a superset of `expr.BooleanEvaluator`'s:

- plain `&`/`|`/`~` expressions mean the same thing in both;
- `BooleanEvaluator` stays the evaluator for `--system_expr`, and for `--domain_expr` in the
  existing commands, which work on flat name sets. In the system commands, `--domain_expr`
  becomes a short-form model, so it is parsed by `conditions.py`;
- the tokenizer is shared (`expr.tokenize(expression, operators)`), so the two read names the
  same way. `BooleanEvaluator` keeps its own precedence, without `<`.

## Reuse, not reinvention

| Need | Existing code to use |
|---|---|
| Ordered CDSs with their hits, ncRNA pseudo-CDSs | `model.DomainatorCDS.list_from_contig(..., include_nucleic_acid_annotations=True)` |
| Hit thresholds | `model.feature_evalue`, `feature_score`, `evalue_passes` |
| CDS identity across slicing | `model.get_cds_unique_name` / the `cds_id` qualifier |
| CDS-count windows that skip pseudo-CDSs | `select_by_cds.cds_range_members` |
| Envelope across the origin | `location_ops.circular_envelope`; `domainate.build_contig_hit_location` for wrapped locations |
| Wrapped distances | `locations.circular_dist` |
| Padding, slicing, region names | `locations.pad_location`, `slice_record_from_location`, `model.region_id` |
| Neighborhood selection | `select_by_cds.get_cds_neighborhood`, `validate_range`, `deduplicate_and_pad` (generalized below) |
| Greedy overlap filtering, competing-hit overlap | the pattern of `domainate.filter_by_overlap`; `location_ops.overlap_size` |
| Boolean conditions and selection expressions | `expr.BooleanEvaluator` |
| Ranking and `--max_hits` | `dedup.TopHits`, `model.set_hit_best_score` / `get_hit_best_score` |
| Loading references for workers | `domainate.read_references`, `_mp.Shared(build=...)` |
| Name lists | `cli.common.add_list_and_file` |
| Report columns | `report_columns.Column`, `enum_report`'s `COLUMNS` and analysis registries |
| Errors and warnings | `errors.InputError`, `FormatVersionError`; `warnings.warn` for unsupported foreign-format features |

## The two refactors that carry the weight

### 1. Neighborhoods around an envelope, not one CDS

`select_by_cds.get_cds_neighborhood(contig, cds_list, cds_idx, ...)` pads around one focus
CDS. Systems need padding around a run of CDSs: `cds_up` counted from the first member,
`cds_down` from the last, in system orientation.

Generalize the core to `get_neighborhood(contig, cds_list, first_idx, last_idx, strand, ...)`.
The existing function becomes `first_idx == last_idx` with the CDS's own strand. Both
`select_by_cds` and `domain_search` keep their behavior, and their existing tests protect
the refactor. The circular even-split logic (`pad_circular_cds`) works on a location, so it
carries over.

The worker-to-selector protocol (`model.set_hit_scores`/`set_hit_names`, keyed by CDS
feature index) needs a system equivalent: `model.set_system_hits(record, [(first_idx,
last_idx, strand, score, name)])`. It is set by the `system_search` worker and read by the
generalized selector.

**Fix first:** the exploration found that in `domainate.domainator_inner`
(`domainate.py:1055-1058`), a nucleotide contig with both contig-level and CDS hits has its
`hit_scores`/`hit_names` overwritten by `add_nucleic_acid_annotations`. A system protocol
built on the same pattern would inherit the bug, and contigs with both nucleotide and CDS
members are exactly the case systems need. Fix it, with a test, before generalizing.

### 2. A detection hook in the domainate worker

`domainate --system_models` and `system_search` both run detection on each annotated
contig. The hook goes after `domainator_inner` in `_domainate` (which is private, so its
parameter needs no `_` prefix). `domainate` searches in-process; `system_search` calls
`_domainate` from its pool tasks, as domain_search does. It also needs:

- reference loading that merges the models' `reference_files` with `-r` (a converted
  directory names its profile file after the collection's database, so hits are labeled the
  same on both routes);
- `-e` applied as an evidence filter in detection, on top of the looser search threshold of
  model profiles (see [definition_formats.md](definition_formats.md#thresholds-are-applied-at-detection-time));
- per-reference e-value overrides and exemptions from `--max_hits_per_cds`/`--max_overlap`
  for model profiles. These are carried as a set of `(database, name)` keys in the reference
  dict, so `filter_by_overlap` can skip them.

## Public surface

| Addition | Where |
|---|---|
| `SYSTEM_FEATURE_NAME = "Domain_System"`, `SYSTEM_SEARCH_BEST_HIT_NAME = "System_Search"` | `domainator/__init__.py`, `__all__` |
| `find_systems`, `system_search`, `select_by_system`, `system_db_download`, `format_systems`, `import_systems` | `from domainator.<cmd> import <cmd>`. The first parameter is positional, all others keyword-only, named and defaulted like the options. Listed in `test/test_library_api.py`. |
| New keyword parameters on `domainate` and the selector, report, extractor and plotting functions | keyword-only, named like the options, default "off", so existing callers are unaffected |
| `model.SYSTEM_TYPES`, `model.system_feature_types` | internal, beside `ANNOTATION_TYPES` (Phase 2, with `--system_type`) |
| Feature types, qualifiers, native definition format | `docs/file_formats.md`, and the GenBank and versioned-format contracts in `stability.md` |

`system_models` parameters take paths or already-read `ModelCollection`s, so library users
can build models in code, or convert a foreign format in memory with
`domainator.systems.formats.convert` and skip writing a directory:

```python
from domainator import parse_seqfiles, write_genbank
from domainator.find_systems import find_systems

records = find_systems(parse_seqfiles(["annotated.gb"]), system_models=["my_systems.yaml"])
write_genbank(records, "systems.gb")
```

## Tests

All fixtures go in `test/data/systems/`, are staged with `shared_datadir`, and are never
mutated.

- **Synthetic contigs with hand-written `Domainator` features**, so detection tests need no
  HMMER:
  - linear and circular cases, both strands;
  - a system across the origin;
  - gaps exactly at the limit and one over, in CDS and kb;
  - an ncRNA pseudo-CDS between members (it must not count as a gap gene);
  - a forbidden gene inside a cluster;
  - exchangeable references;
  - a CDS matching two components;
  - a CDS carrying two domains named in a short-form condition (fusions satisfy both);
  - a nucleotide hit inside a CDS, with and without `attach_to_cds`;
  - a record with no CDSs: `max_kb` models evaluated, `max_cds`-only models skipped with one
    warning;
  - orientation voted by CDS members only, with an opposite-strand ncRNA member;
  - tandem duplicate systems;
  - two models sharing members (overlap policy);
  - a system at a linear contig end (`system_partial`);
  - span mode: a negated name within reach vetoes and one out of reach doesn't; a stray
    duplicate doesn't stretch the system;
  - a fusion CDS listed once per component it fills;
  - anchor mode: an anchor that didn't contribute doesn't count.
- **Window modes**, using the worked example in
  [definition_formats.md](definition_formats.md#windows) (members at CDS positions 1, 4, 7
  and 10) as a shared fixture:
  - `chain`, `anchor` and `span` each give exactly the documented result;
  - limits exactly at the boundary and one over, in CDSs and in kb;
  - `max_span_*` caps splitting a long chain;
  - a span window choosing the smaller of two passing windows;
  - the anchor join, where two anchors' cores merge;
  - wrap across the origin in each mode;
  - antiSMASH's `<` cutoff converted to `<=`, tested at the exact bp boundary.
- **Requirement and resolution:**
  - each condition construct (`cds`, `minimum`, `repeatable`, `minscore`, `define`
    aliases), including negation inside `cds()`;
  - `<` between members on both strands, with a mixed-strand window, with
    `orientation_from`, across the origin, with an alias operand, and as a three-member
    chain; `~a < b`; a and b on one CDS (false outside `cds()`); `cds(a < b)` on a
    minus-strand CDS (N- to C-terminal); a name containing `<` refused;
  - a rule with no positive term is refused;
  - competing groups (overlap above and below `competing_overlap_aa`, near-ties kept);
  - extenders growing each end, with circular wrap;
  - superiors dropping an overlapping inferior system and keeping a separate one.
- **Native format:**
  - round trip (`write_system_models` then `read_system_models`), including the writer
    choosing the short form;
  - the short form (implicit components, the mapping shorthand, `defaults.window`) gives the
    same `SystemModel` as the equivalent full form, and the same as the command-line options
    (`--domain_expr`, `--window_cds`): one fixture, three spellings, identical systems;
  - a version newer than supported gives a `FormatVersionError`;
  - malformed files and unknown condition names give an `InputError`;
  - a foreign-format path given to `--system_models` gives an `InputError` naming the
    `system_db_download` command.
- **Converters** (`system_db_download --input`), each checked for its `systems.yaml`,
  `profiles.hmm`, `source:` block and `conversion.tsv`:
  - **PADLOC:** a trimmed padloc-db excerpt (MIT: vendor a few `sys/*.yaml` files and
    `hmm_meta.txt` rows with the license notice), covering pipe aliases, group expansion,
    role-overlap resolution, NA thresholds and `CRISPR_array` dropping.
  - **MacSyFinder:** a synthetic package written for Domainator (no DefenseFinder content),
    covering exchangeables, per-gene gap overrides with the `min` rule, `vers` 2.0/2.1
    booleans, `multi_loci` and `loner` listed in `conversion.tsv`, and a v1 `<system>`
    refusal.
  - **antiSMASH rules:** synthetic rule files and a tiny `hmmdetails.txt` and
    `filterhmmdetails.txt`, referencing small Domainator-made HMMs. No antiSMASH content,
    because it is AGPL. They cover:
    - every grammar construct and comments;
    - `DEFINE`;
    - strictness tiers loaded cumulatively;
    - parse errors with line and column;
    - a dynamic profile becoming an empty component listed in `conversion.tsv`, and its
      substitution with `--substitutes_file`.

    A `network`-marked test clones antiSMASH at a pinned tag and checks that every shipped
    rule converts.
  - **`import_systems`:** a small antiSMASH GenBank output from a public genome, at each
    `--level`.
- **GenBank:**
  - `Domain_System`/`System_Search` round trip through `write_genbank`/`parse_seqfiles`,
    including `|` in a reference name;
  - slicing drops cut systems;
  - `cds_id` references survive `select_by_cds`.
- **Commands:**
  - `find_systems` streaming (`-i -`);
  - `find_systems -e` as an evidence filter combined with model thresholds;
  - `format_systems --yaml_out` from command-line options gives a file that reproduces the
    same systems;
  - `domainate --system_models` with tiny HMMs, giving the same systems as `domainate`
    followed by `find_systems`;
  - `system_search` windows measured from the system's ends on both strands, `--anchor`,
    `--whole_contig` on circular contigs, and `--max_hits` ranking;
  - `--min_context_cds`/`--min_context_kb` and their `_up`/`_down` forms on both strands, at
    linear contig ends and on circular contigs;
  - `find_systems --clear` keeping another command-line model's systems;
  - `select_by_system` gives the same regions as `system_search` on the same input;
  - the `enum_report --by system` columns.
- **Contracts:** `test_library_api.py`, `test_defaults_match_cli.py`, and `test_manifest.py`
  with a regenerated snapshot.
- **`network`-marked tests:** `system_db_download` padloc-db and defense-finder-models, run
  on one complete genome, and assert that the known systems are found (e.g. a genome with a
  documented Gabija).
- **Viewer and plot:** a Playwright smoke test that a `Domain_System` feature renders in
  `plot_contigs` HTML. It is optional and skipped without playwright.

## Documentation

- `docs/file_formats.md`: "Systems" (feature types, qualifiers, `component` grammar, location
  and slicing rules) and "System definition files" (native format, version key).
- `docs/stability.md`: the new feature types, constants, commands, and the definition-format
  version.
- `docs/systems.md`: a user guide. It opens with the simple case (`--domain_expr` with
  `--window_cds`, then the short form), and only then covers the full format, converting
  DefenseFinder, PADLOC and antiSMASH models, the recommended pipelines (from the
  [README](README.md#example-pipeline-recommended-route)), fidelity notes, and licensing of
  model collections.
- `docs/structure_search.md`: a "Systems with structural components" pointer.
- `docs/examples.md`: a systems section. `README.md`: the program table.
- `CHANGELOG.md`: a 1.x feature release entry. Bumping `__version__` is the maintainer's
  call, per `CLAUDE.md`.
- `docs/README.md`: link the user guide, and this design folder while it is under
  discussion.

## Roadmap

| Phase | Delivers | Depends on |
|---|---|---|
| **1. Foundations** | Feature types and constants; `systems/` model with `Window` (all three modes), `category`, `superiors`, `neighborhood_*`, `examples`, `attach_to_cds`; native reader and writer (format 1), full and **short form**, with plain boolean conditions and `<` between members; `windows.py`, `detect.py` (including competing groups, superiors and the `system_partial` flag), `features.py`; `find_systems` with **command-line models** (`--domain_expr`, `--window_*`) and `-e`; slicing rule; `color_genbank --target systems`; `plot_contigs` default types; file-format docs | - |
| **2. Everyday use** | **`system_db_download --source padloc`** (the converter framework, `conversion.tsv`, download and `--input`); `domainate --system_models` (detection hook, model-profile thresholds); the hit-score bug fix and neighborhood generalization; `select_by_system` with `--model_neighborhood` and `--min_context_*`; `enum_report` columns and `--by system`; `summary_report` tables; `select_by_contig` and `extract_peptides` selectors; `format_systems`; **`import_systems --source antismash`** | 1 |
| **3. Search and DefenseFinder** | `system_search`, with command-line models and `--min_context_*`; `system_db_download --source defensefinder/macsyfinder` (exchangeables, gap overrides, GA as `min_score`); `select_by_cds`, `extract_domains`, `trim_contigs`, `rotate_contigs`, `sort_contigs`, `split_by_cluster` options; `select_by_contig --window_*`; structural thresholds (`min_tmscore`, ...); `docs/systems.md`; `import_systems` for DefenseFinder and PADLOC call tables | 2 |
| **4. Depth and antiSMASH rules** | Full condition grammar (`cds`, `minimum`, `repeatable`, `minscore`, `define`, and `cds(a < b)`; format 1.1); `extenders`; **`system_db_download --source antismash`** with `--substitutes_file`; `format_systems --check_examples`; `--partial_systems include`; `compare_contigs --sji`; `domainate --bit_cutoffs`; non-Domainator feature components (CRISPR arrays); possibly `merge_systems` (candidate-cluster-like grouping) | 3 |
| **Future add-on** | [Multi-locus systems and loners](multi_locus.md): linked per-locus `Domain_System` features (`/system_id`, `/locus`), per-locus region selection, the bounded union search | 4, and demand from MacSyFinder users |

Phase 1 is useful on its own: anyone with a domainate-annotated file can find, color and plot
systems with any of the three window modes, from one `--domain_expr` or a short YAML file. Each later phase adds
commands or options without changing earlier ones, which keeps every step inside the 1.x
addition-only contract.

Fixing the `Window` structure and the condition grammar in Phase 1, even though the grammar's
functions ship in Phase 4, is what keeps the antiSMASH converter cheap. The converter then
only translates syntax, and needs no changes to the data model.
