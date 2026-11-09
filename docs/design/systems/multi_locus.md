# Systems: multi-locus systems (future add-on)

Status: **future add-on**. It is not part of phases 1-4 in the
[roadmap](implementation.md#roadmap). Until it lands, every system has one contiguous
location, apart from wrapping the origin of a circular contig.

- `multi_loci` models (MacSyFinder) are evaluated as single-locus.
- `loner` components (MacSyFinder) are treated as ordinary components, so they must be
  inside the window.
- `system_db_download` reports both approximations in its conversion summary.

This document records the analysis so the add-on can be built without changing the 1.x data
model.

## What is in scope

- **Multi-locus systems** ("scattered", in DefenseFinder's terms): one system made from
  clusters that sit far apart on one replicon. One example is a secretion system whose
  genes are split between two operons.
- **Loners:** a component that may be found anywhere on the replicon, outside the window,
  such as a protein encoded far from the rest of the system. A loner is a one-CDS locus, so
  it needs the same mechanism.

Only MacSyFinder models use either feature. PADLOC and antiSMASH have nothing comparable.

## Representation: one feature per locus

Write one `Domain_System` per locus, linked by a group id:

```
Domain_System   1200..5400
                /name="T4P"  /system_id="T4P_3"  /locus="1/2"  ...
Domain_System   88200..91050
                /name="T4P"  /system_id="T4P_3"  /locus="2/2"  ...
```

| Qualifier | Value |
|---|---|
| `system_id` | Unique within the record. It is shared by every locus of one system. |
| `locus` | `<i>/<n>`: this locus's position among the system's loci, in system orientation. |
| the others | `score`, `evalue`, `completeness`, `core_count` and `total_count` describe **the whole system** and are repeated on every locus. `component` lists only this locus's members. `member_count` counts only this locus's members. |

Why:

- Every location stays contiguous, so each feature's envelope is its locus. Neighborhood
  selection, slicing, CDS counting, report rows and plotting work unchanged: each locus is a
  region.
- Slicing keeps the local part of a system and drops the rest, which suits neighborhood
  work. A sliced record keeps `locus="1/2"`, so readers can tell that the system continues
  elsewhere.
- Readers that don't know `/locus` still see sensible, if incomplete, features.

The costs:

- One feature no longer passes its requirement on its own. That bends "systems are
  self-contained snapshots" ([data_model.md](data_model.md#lifecycle)).
- Anything that counts systems (`enum_report`, `summary_report`, `select_by_contig`, the
  `--max_hits` ranking) must group by `system_id`. A helper in `systems/records.py`
  (`systems_in(record)` yielding grouped systems) keeps that in one place.

## Alternative: one feature with an `order()` location

A system could also be one feature with a compound location. INSDC has two operators:

- `join(a..b,c..d)` means "the parts are joined to form the feature", as in a spliced CDS.
- `order(a..b,c..d)` means "found in this order; nothing is implied about joining".

`order()` is the correct operator for a scattered system. Domainator already writes `join`
for features that wrap the origin, so `order` also keeps the two meanings apart. Biopython
reads and writes it.

### Compared with multi-part CDSs

Domainator already supports CDSs with introns, so the location code copes with multi-part
locations:

- `overlaps`, `contains`, `overlap_size` and the feature mapping in
  `slice_record_from_location` work part by part;
- the merge of adjacent parts after slicing is limited to `join` (`locations.py:513`);
- `plot_contigs` passes every part of a feature to the drawing code.

What differs is what the gaps mean:

| | Multi-part CDS | Scattered system |
|---|---|---|
| Where the parts come from | the input annotation | inferred by detection |
| Operator | `join`: the parts are spliced together | `order`: nothing is joined |
| What the gaps contain | introns, which rarely contain other genes | ordinary genes, often many, sometimes other systems |
| Is the envelope meaningful? | yes. Gene plus introns is the region of interest, so neighborhoods are padded from it. | no. It mixes the system with everything between its loci. |
| How parts are used | `overlay` maps domain hits from protein coordinates through the exons; extraction concatenates them for translation | no mapping or concatenation; each part is a separate region of interest |
| Counts as | one CDS in every window, whatever its length | as many loci as it has. `--cds_up`/`--cds_down` need a per-locus meaning. |
| Strand and order | one strand, parts in transcription order, gaps bounded by intron length | members on either strand; on a circular contig, two distant loci can be ordered either way around the circle (one order crosses the origin, the other doesn't) |
| A slice through a gap | the CDS is truncated, with `<`/`>` markers, and remains a valid partial gene | the system is dropped, because part of a system doesn't meet its requirement |

Domainator's convention "a multi-part feature is its envelope, for neighborhood purposes"
is right for genes and wrong for systems. Selectors reduce a focus location to its
`circular_envelope` before padding (`select_by_cds.py:190`, `:277-310`). The envelope of two
loci 80 kb apart is an 80 kb region of mostly unrelated genes. With one `order()` feature,
every system consumer (neighborhoods, "the CDSs in this system", `--by system` rows, length
and sort keys) would need its own "envelope or parts" rule. `circular_envelope` would also
need two conventions: every part gets the system's strand, and parts are listed in system
orientation. Otherwise it reads a step "backward" between parts as an origin crossing.

Per-locus features avoid all of that. Their cost is limited to detection, the reports, and
the grouping helper.

## Selecting regions

With per-locus features, `system_search` and `select_by_system` produce **one region per
locus**. Each region carries a `System_Search` for its locus, with the same `system_id` and
`locus`. The alternatives are rejected:

- one region from the first locus to the last may be most of a chromosome;
- concatenating the loci into one record breaks coordinates.

So "one region per system" becomes "one region per locus". `--max_hits` and
`--max_hits_per_contig` still count systems, by `system_id`.

## Detection

MacSyFinder searches over combinations of clusters. A bounded version, added after step 6 of
[detection.md](detection.md#6-test-the-requirement):

1. Form single-locus candidate windows as now. For `multi_loci` models, keep windows that
   fail the requirement but contain at least one counted member.
2. For each `multi_loci` model on one replicon, try unions of its windows, smallest first.
3. Keep the best-scoring union that meets the requirement, remove its windows, and repeat.

Contigs rarely have more than a few clusters per model, so trying every union, with a cap
(for example, 2^10 unions), is cheap. The result need not match MacSyFinder's exact choice.

**Loners:** a `loner` component may be filled by a hit anywhere on the replicon that isn't
already a member of another system of the same model. The loner becomes its own one-CDS
locus.

## Phase

The add-on comes after the four phases, and only if MacSyFinder users ask for it. It needs:

- the two qualifiers;
- the grouping helper and its use in counting tools;
- per-locus region selection;
- the detection step above;
- an update to the fidelity table in [detection.md](detection.md#fidelity).

No 1.x contract changes: the qualifiers are additions, and single-locus systems are
unaffected.
