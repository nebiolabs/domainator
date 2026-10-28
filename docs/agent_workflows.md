[index](README.md)
# Agent-driven Domainator workflows

Copy-paste recipes for AI agents (or scripted pipelines) that need
**low-token, high-signal** output from Domainator. For orientation and the
CLI conventions, read [AGENTS.md](https://github.com/nebiolabs/domainator/blob/main/AGENTS.md) first.

The guiding idea: Domainator's HTML reports and binary artifacts (HDF5
matrices, XGMML, the `.dvb` viewer bundle) are meant for humans and browsers.
For an agent, reach for the `--json_out` output modes and `domainator ssn_navigator`, which
return compact, structured slices instead of whole-artifact dumps.

## Machine-readable output cheat-sheet

| Tool | Flag | Shape |
| ---- | ---- | ----- |
| `domainator enum_report` | `--json_out -` | NDJSON, one object per contig/CDS/domain |
| `domainator hmmer_report` | `--json_out -` | NDJSON, one object per HMM profile (add `--columns alphabet` for amino/DNA/RNA) |
| `domainator summary_report` | `--json_out -` | one JSON object (contig stats, domain frequency + co-occurrence, taxonomy) |
| `domainator matrix_report` | `--json_out -` | one JSON object (edge stats, connected components, split-event series) |
| `domainator ssn_navigator` | (always JSON) | targeted cluster / node / threshold queries |

`--json_out -` writes to stdout; pass a path to write a file. These sit alongside the
existing `-o`/`--html_out` outputs — nothing is removed.

## Output shapes

The `enum_report` / `hmmer_report` NDJSON is **dynamic**: each object's keys are
exactly the columns you requested via flags (e.g. `--architecture --length`),
and values are strings, numbers, or `null`. There is no fixed key set — read the
first line to see the shape. Example (`enum_report --by contig --columns architecture length`):

```json
{"contig":"FeSOD_A0A1F4ZT98|unreviewed|Superoxide","architecture":"Sod_Fe_N; Sod_Fe_C","length":196}
```

The remaining outputs have **stable shapes**, shown below (pretty-printed here;
the tools emit them compact). These are the contract — treat added keys as
backward-compatible, removed/renamed keys as breaking.

### `summary_report --json_out`

```json
{
  "contig_stats": {
    "contigs": 20, "cdss": 0, "cds_per_10kb": 0.0,
    "length": {"min": 189, "max": 215, "mean": 200.15, "median": 201.0, "total": 4003}
  },
  "domain_frequency": [
    {"domain": "Sod_Fe_N", "database": "FeSOD_pfam", "count": 20,
     "description": "Iron/manganese superoxide dismutases, alpha-hairpin domain", "avg_score": 101.1}
  ],
  "domain_cooccurrence": null,
  "taxonomy": [{"taxid": 1234, "name": "Bacillus", "rank": "genus", "count": 5}]
}
```

`length` is `null` for empty input; `domain_cooccurrence` is `null` unless focus
domains were given (then it is `{row_domain: {col_domain: percent}}`); `taxonomy`
is present only with `--taxonomy`.

### `matrix_report --json_out`

```json
{
  "nodes": 20, "non_zero_edges": 52,
  "edge_scores": {"mean": 215.37, "median": 216.0, "min": 135.0, "max": 287.0},
  "connected_components": 2,
  "merge_impact_metric": "min_child",
  "split_events": [
    {"edge_index": 0, "threshold_from": "∞", "threshold_to": "287.00",
     "threshold_value": 287.0, "merge_impact": 1.0,
     "merge_size_counts": {"1": 1}, "largest_merge": 1.0, "merge_count": 1,
     "delta_largest": 1.0, "delta_avg_non_singleton": 2.0}
  ]
}
```

`edge_scores` is `null` for a matrix with no non-zero edges. `split_events` is bounded by
`domainator matrix_report --max_merge_events`; a `threshold_value` of `null` denotes ∞.

`merge_impact` is the **sum** of every merge at that threshold, so it cannot by itself
distinguish one large cluster splitting off from a swarm of small ones. The other three
keys describe its terms: `merge_size_counts` maps each distinct merge size to how many
merges of that size occurred, `merge_count` is their total, and `largest_merge` is the
biggest single one. `sum(size * count) == merge_impact` holds exactly, and
`largest_merge == merge_impact` only when `merge_count` is 1. A threshold reporting
`{"3": 1, "1": 33}` is one 3-node cluster plus thirty-three single nodes falling off at
once -- `merge_impact` 36, but no cluster larger than 3.

Every threshold matrix_report reports means the same thing as `build_ssn --threshold`: the
graph of edges scoring **strictly above** the threshold. So a reported threshold `T`
and its edge/cluster counts are exactly what `domainator build_ssn -i <matrix> --threshold T` emits,
and `split_events[i].edge_index` is the last MST edge scoring above that threshold.
Thresholds equal to a score are therefore exclusive of that score's own ties -- the
merge listed at threshold `T` becomes visible at the next lower threshold.

### `ssn_navigator` (one shape per `--mode`)

```jsonc
// --mode overview
{"name": "net", "domainator_version": "0.10.0", "nodes": 4, "mst_edges": 3,
 "connected_components": 1,
 "metadata_columns": [{"name": "category", "type": "str"}, {"name": "count", "type": "int"}],
 "defaults": {"color_by": null, "label_by": null}, "merge_impact_metric": "min_child"}

// --mode thresholds  -> suggested cut-points, finest first (threshold_value null == ∞)
{"thresholds": [{"threshold_value": 7.0, "threshold_label": "7.00",
                 "clusters": 2, "non_singleton_clusters": 1, "largest_cluster": 3}]}

// --mode clusters --threshold X
{"threshold": 5.0, "clusters": 2, "singletons": 1, "shown": 2,
 "cluster_list": [{"id": 5, "size": 3}, {"id": 3, "size": 1}]}
// "truncated": <n> is added when more than --top_n clusters exist.

// --mode cluster --threshold X (--id N | --node NAME) [--members]
{"threshold": 5.0, "cluster_id": 5, "size": 3,
 "metadata_summary": [
   {"name": "category", "type": "str", "count": 3, "missing": 0, "unique": 2,
    "top": [{"value": "alpha", "count": 2}, {"value": "beta", "count": 1}]},
   {"name": "count", "type": "int", "count": 3, "missing": 0,
    "min": 1.0, "max": 3.0, "mean": 2.0, "p25": 1.5, "median": 2.0, "p75": 2.5}
 ]}
// "members": ["A","B","C"] is present only with --members.
// str columns carry top/unique (+ "truncated" past --top_n); numeric columns carry the quartile stats.

// --mode node --node NAME
{"node": "A", "node_index": 0, "metadata": {"category": "alpha", "count": 1},
 "memberships": [{"threshold_value": 7.0, "threshold_label": "7.00",
                  "cluster_id": 5, "cluster_size": 3}]}
```

## Recipe 1 — annotate, then get a structured summary

```bash
domainator domainate -i genomes.gb -r pfam.hmm -o annotated.gb --max_overlap 0.6
domainator summary_report -i annotated.gb --json_out summary.json
domainator enum_report -i annotated.gb --by contig --columns architecture length --json_out contigs.ndjson
```

`summary.json` gives dataset-level domain frequencies and (with `--taxonomy`)
taxonomic breakdown; `contigs.ndjson` gives one line per contig with its domain
architecture — filter/aggregate it with any JSON tooling.

## Recipe 2 — build a similarity network and navigate it

```bash
# 1. all-vs-all similarity -> dense matrix
domainator seq_dist -i proteins.gb -r proteins.gb --hdf_out dist.hdf5

# 2. (optional) per-node metadata to attach, e.g. domain architecture + taxonomy
domainator enum_report -i proteins.gb --by contig --columns architecture taxname:superkingdom -o meta.tsv

# 3. compact SSN bundle
domainator build_ssn -i dist.hdf5 --dvb_out net.dvb --metadata meta.tsv

# 4. navigate without rendering
domainator ssn_navigator -i net.dvb -o - --query overview
domainator ssn_navigator -i net.dvb -o - --query thresholds                 # pick a cut-point
domainator ssn_navigator -i net.dvb -o - --query clusters --threshold 150   # clusters + sizes
domainator ssn_navigator -i net.dvb -o - --query cluster --threshold 150 --cluster 5   # members + metadata distribution
```

The `metadata.tsv` columns (from `enum_report`) become the per-cluster
distributions in `--mode cluster` output: string columns yield the top values
with counts; numeric columns yield min/max/mean/quartiles.

`--mode thresholds` derives the cuts from the bundle's merge order and reports
`merge_events` / `merge_event_total` alongside them: those cuts are one per *selected*
split event, which `domainator ssn_navigator --max_merge_events` caps (default 500, plus up to 20
more so no 5% band of the threshold range is left unrepresented). On a large network that
is a small sample of the real event list, so treat the cuts as good candidate cut-points
rather than as every threshold at which anything happens — and pass a higher cap, or `0`,
to see them all. The cap is yours, not the file's: it is applied when you ask, so any
bundle answers at any cap without being rebuilt.

`--threshold` follows the viewer's slider: higher is finer. `--threshold inf` is the
slider's ∞ end, where every node is its own cluster (reported as `"threshold": null`,
the same spelling the bundle uses); omitting `--threshold` is the other end, one cluster
per connected component.

Add `--html_out net.html` in step 3 for a self-contained browser viewer
of the same network. Annotations made there (edited cells, new columns) are saved
back into an ordinary `.dvb` by its "Save session…" button, so a session file
feeds straight back into step 4 and the new columns show up in the per-cluster
distributions. The chart glyph in each metadata column header is the interactive
counterpart of those distributions -- a frequency table, bar or pie chart,
histogram, box plot, cumulative distribution or summary statistics over the rows
the table is currently showing, exportable as SVG, PNG or TSV.
[ssn_viewer.md](ssn_viewer.md) is the guide to the viewer; the `.dvb` section of
[file_formats.md](file_formats.md) has the schema.

## Recipe 3 — inspect a matrix's cluster structure

```bash
domainator seq_dist -i proteins.gb -r proteins.gb --hdf_out dist.hdf5
domainator matrix_report -i dist.hdf5 --json_out -
```

Returns edge-score stats, the number of connected components, and the strongest
threshold split events — enough to decide where to cut the network before
building an SSN.

## Notes on token efficiency

- `domainator ssn_navigator` and the report `--json_out` modes bound their lists with
  `--top_n` (default 50). Raise it deliberately when you need the long tail.
- `domainator ssn_navigator --query cluster` omits the member id list unless you pass
  `--members`; the metadata distribution summary is usually what you want.
- NDJSON (`enum_report`/`hmmer_report`) streams, so you can `head`/filter huge
  reports without materializing the whole thing.
