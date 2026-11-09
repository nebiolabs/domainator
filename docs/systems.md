# Systems

A **system** is a set of domains found together in one genomic neighborhood: a defense system
such as Gabija (GajA next to GajB), a toxin-antitoxin pair, an operon, or a biosynthetic gene
cluster. Domainator finds systems on annotated contigs and records each one as a
`Domain_System` feature that spans its members and lists them, so the other tools can select,
extract, color, plot and report by system.

This guide starts with the simplest case, a system defined on the command line, and then
covers model files, published model collections (PADLOC, DefenseFinder and other MacSyFinder
packages, antiSMASH results), the recommended pipelines, and the commands that work with
systems. The file formats are specified in [file_formats.md](file_formats.md#systems), and the
design in [design/systems](design/systems/README.md).

## The simple case: a domain expression and a window

`find_systems` reads records that already have domain annotations (from `domainate`) and finds
the systems a model describes. The quickest model is a boolean expression over domain names
plus a window:

```bash
domainator domainate -i genomes.gb -r ta_profiles.hmm -o annotated.gb
domainator find_systems -i annotated.gb --domain_expr "RelE & (HicB | HigA)" --window_cds 2 \
    --system_name RelE_TA -o systems.gb
```

This finds every place where a CDS with a RelE hit and a CDS with a HicB or HigA hit lie within
a window of 2 CDSs. The expression uses `&` (and), `|` (or), `~` (not), parentheses, and `<`
(upstream of, in the system's orientation): `"RelE < HicB"` requires the RelE gene first.
`database:name` limits a name to one database, such as `pfam:RelE`.

`system_search` does the search and the selection in one step, for large databases, the way
`domain_search` does for single domains:

```bash
domainator system_search -i genomes.gb -r ta_profiles.hmm --domain_expr "RelE & HicB" \
    --window_cds 2 --cds_range 5 -o ta_neighborhoods.gb
```

Each output region has the system plus 5 CDSs on each side, with a `System_Search` feature
marking the system and the hits of its members.

### Windows

`--window_mode` says how "close together" is measured:

| Mode | `--window_cds 10` means | Used by |
|---|---|---|
| `span` (the default) | the whole system, first member to last, covers at most 10 CDSs | Domainator |
| `chain` | at most 10 CDSs between consecutive members | PADLOC, MacSyFinder |
| `anchor` | every member within 10 CDSs of one anchor member | antiSMASH |

`--window_kb` sets the same limit in kb, and with `--window_cds` both must hold. CDS counts
include only protein-coding CDSs, not hits to nucleotide references. `--same_strand` requires
all members on one strand.

A negated name vetoes a system when it is close enough to have joined it:
`"GajA & GajB & ~Tnp"` rejects a Gabija with a transposase gene that would fit in the same
window.

## Model files

When a definition grows, keep it in a file. The short form is one line per system:

```yaml
domainator_systems: 1
database: my_systems
defaults:
  window: {mode: chain, max_cds: 3}
  evalue: 1e-10
systems:
  RelE_TA: "RelE & (HicB | HigA)"
  Gabija: "GajA & GajB & ~Tnp"
```

```bash
domainator find_systems -i annotated.gb --system_models my_systems.yaml -o systems.gb
```

The full form names each component, its role and its references, with per-reference
thresholds, and supports quorums ("at least 3 of these 5"), as PADLOC and MacSyFinder models
do:

```yaml
domainator_systems: 1
database: my_systems
systems:
  - name: Zorya_like
    window: {mode: chain, max_cds: 4}
    min_core: 2          # at least 2 core components
    min_total: 3         # at least 3 core or accessory components
    components:
      - {name: ZorA, role: core, domains: [ZorA_1, ZorA_2]}
      - {name: ZorB, role: core, domains: [{name: ZorB, min_score: 30}]}
      - {name: ZorC, role: accessory, domains: [ZorC]}
      - {name: ZorD, role: accessory, domains: [ZorD]}
      - {name: Tnp, role: forbidden, domains: [Transposase]}
```

Roles: `core` and `accessory` components count toward the quorum, the score and completeness;
`neutral` ones may join a system without counting; a `forbidden` one vetoes any system it is
near. A domain can be a mapping with a `database` and thresholds (`evalue`, `min_score`,
`min_reference_coverage`, `min_target_coverage`, `min_target_length`, `max_target_length`, and
for structural hits `min_tmscore`, `min_lddt`, `min_prob`). `format_systems --yaml_out` writes
a `--domain_expr` model as a file to start from, and `--tsv_out` lists every component and
threshold of a model file as a table. The format is specified in
[file_formats.md](file_formats.md#system-model-files).

## Published model collections

`system_db_download` converts other tools' model collections, once, into a directory that
`--system_models` reads anywhere:

```bash
domainator system_db_download --source padloc --output_dir padloc_db/
domainator system_db_download --source defensefinder --output_dir defensefinder_db/
domainator system_db_download --source macsyfinder --input CasFinder/ --output_dir casfinder_db/
```

The directory holds `systems.yaml`, the profiles the models use (`padloc.hmm`,
`defensefinder.hmm`, ...), the source's `LICENSE`, and `conversion.tsv`, which lists everything
the conversion approximated or dropped. With `--input`, a local copy is converted instead of
downloading one (padloc-db, a defense-finder-models directory as installed by `defense-finder
update`, or any MacSyFinder model package).

Converted models keep their source's names (`gabija`, `cas_type_I-E` for PADLOC; `Gabija`,
`RM_Type_II` for DefenseFinder), and their descriptions, categories and accessions:
PADLOC's human-readable names (`Gabija`, `CBASS (other)`) are in the `description` qualifier,
and DefenseFinder's model paths (`DefenseFinder/Gabija/Gabija`) in `accession`.

antiSMASH's calls can be brought in from its own GenBank output, without converting its rules:

```bash
domainator import_systems -i antismash_output.gbk --source antismash --level protocluster -o bgcs.gb
```

## Pipelines

```bash
# one step: search the models' profiles and annotate the systems
domainator domainate -i genomes.gb --system_models padloc_db/ -o annotated.gb

# two steps: annotate once, permissively, then try definitions cheaply
domainator domainate -i genomes.gb -r padloc_db/padloc.hmm -e 1e-3 -o annotated.gb
domainator find_systems -i annotated.gb --system_models padloc_db/ -o systems.gb
domainator find_systems -i annotated.gb --system_models my_systems.yaml -o systems2.gb

# a large database: regions around one kind of system, 5 CDSs beyond each end
domainator system_search -i genomes.gb --system_models defensefinder_db/ --systems Gabija \
    --cds_range 5 --min_context_cds 5 -o gabija_regions.gb
```

`domainate --system_models` searches each model profile with its model's own threshold (or
`-e`, if that is looser), and keeps every hit of those profiles, whatever `--max_hits_per_cds`
and `--max_overlap` say, so that detection can choose between them. A hit then counts as
evidence only if it also passes `-e`. The two-step route finds the same systems when the first
step keeps the hits the models need: use an `-e` at least as loose as the models' thresholds,
and leave out `--max_hits_per_cds` and `--max_overlap`.

## Working with systems

| Task | Command |
|---|---|
| regions around systems, padded by CDSs or kb from their outermost members | `select_by_system --cds_range 5` (or `--model_neighborhood`, the model's own neighborhood) |
| skip systems cut by a contig end | `--min_context_cds 5` (or `_up`/`_down`, `--min_context_kb`) on `select_by_system` and `system_search` |
| contigs with a combination of systems | `select_by_contig --system_expr "Gabija & ~Thoeris_I"` |
| contigs where domains co-occur within a window, without writing systems | `select_by_contig --domain_expr "A & B" --window_cds 5` |
| the members of systems, with neighbors | `select_by_cds --systems Gabija --cds_range 2` |
| the proteins of members | `extract_peptides --systems Gabija` (each gets a `system_membership` qualifier) |
| the domains that made a system | `extract_domains --systems Gabija` |
| trim contigs to their systems | `trim_contigs --systems Gabija` |
| rotate circular contigs to a system | `rotate_contigs --systems Gabija --anchor start` |
| sort `system_search` output by score | `sort_contigs --system_type system_search` |
| color and plot | `color_genbank --target systems_and_cdss --system_color_by category`, then `plot_contigs` |
| tables | `enum_report --by system --columns systems system_components ...`, `summary_report` (system frequency, co-occurrence, and systems by taxon) |

Systems are snapshots: `filter_domains` and other domain editors leave them as they are.
Re-run `find_systems` (which replaces the systems of the same models) after changing the
annotations. Slicing tools drop a system they cut through.

## How converted models behave

Converted models are evaluated with Domainator's semantics, which follow their sources closely
but not exactly. The differences are listed in each directory's `conversion.tsv`; in general:

- **PADLOC.** CRISPR arrays can't be detected, so they are dropped from the models (and their
  `minimum_total` lowered if needed). PADLOC measures reference coverage as `(to - from) /
  length`, one residue less than Domainator.
- **DefenseFinder and MacSyFinder.** Each profile's domain gathering threshold becomes a
  `min_score`, which approximates `hmmsearch --cut_ga` (that also checks the sequence
  threshold). Loners and multi-locus models are found as single loci: their genes must lie
  within the cluster.
- **antiSMASH** (imported calls): members are the core CDSs whose profiles triggered each
  product, and protocluster systems are padded by antiSMASH's neighborhood with
  `--model_neighborhood`.

## Licensing of model collections

Domainator ships no model content. `system_db_download` downloads a collection to your machine
and copies its license next to the converted files: padloc-db is MIT licensed;
defense-finder-models' metadata gives both GPL-3.0 and CC BY-NC-SA 4.0 (non-commercial). Check
the `LICENSE` file before redistributing a converted directory.
