[index](README.md)
# Agent-driven Domainator workflows

Copy-paste recipes for AI agents (or scripted pipelines) that need
**low-token, high-signal** output from Domainator. For orientation and the
CLI conventions, read [AGENTS.md](../AGENTS.md) first.

The guiding idea: Domainator's HTML reports and binary artifacts (HDF5
matrices, XGMML, the `.ssnv` SSN bundle) are meant for humans and browsers.
For an agent, reach for the `--json` output modes and `ssn_navigator.py`, which
return compact, structured slices instead of whole-artifact dumps.

## Machine-readable output cheat-sheet

| Tool | Flag | Shape |
| ---- | ---- | ----- |
| `enum_report.py` | `--json -` | NDJSON, one object per contig/CDS/domain |
| `hmmer_report.py` | `--json -` | NDJSON, one object per HMM profile |
| `summary_report.py` | `--json -` | one JSON object (contig stats, domain frequency + co-occurrence, taxonomy) |
| `matrix_report.py` | `--json -` | one JSON object (edge stats, connected components, split-event series) |
| `ssn_navigator.py` | (always JSON) | targeted cluster / node / threshold queries |

`--json -` writes to stdout; pass a path to write a file. These sit alongside the
existing `-o`/`--html` outputs — nothing is removed.

## Output shapes

The `enum_report` / `hmmer_report` NDJSON is **dynamic**: each object's keys are
exactly the columns you requested via flags (e.g. `--architecture --length`),
and values are strings, numbers, or `null`. There is no fixed key set — read the
first line to see the shape. Example (`enum_report --by contig --architecture --length`):

```json
{"contig":"FeSOD_A0A1F4ZT98|unreviewed|Superoxide","architecture":"Sod_Fe_N; Sod_Fe_C","length":196}
```

The remaining outputs have **stable shapes**, shown below (pretty-printed here;
the tools emit them compact). These are the contract — treat added keys as
backward-compatible, removed/renamed keys as breaking.

### `summary_report --json`

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

### `matrix_report --json`

```json
{
  "nodes": 20, "non_zero_edges": 52,
  "edge_scores": {"mean": 215.37, "median": 216.0, "min": 135.0, "max": 287.0},
  "connected_components": 2,
  "merge_impact_metric": "min_child",
  "split_events": [
    {"edge_index": 0, "threshold_from": "∞", "threshold_to": "287.00",
     "threshold_value": 287.0, "merge_impact": 1.0,
     "delta_largest": 1.0, "delta_avg_non_singleton": 2.0}
  ]
}
```

`edge_scores` is `null` for a matrix with no non-zero edges. `split_events` is
bounded by `--max_merge_events`; a `threshold_value` of `null` denotes ∞.

Every threshold matrix_report reports means the same thing as `build_ssn --lb`: the
graph of edges scoring **strictly above** the threshold. So a reported threshold `T`
and its edge/cluster counts are exactly what `build_ssn.py -i <matrix> --lb T` emits,
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

// --mode thresholds  -> suggested cut-points (threshold_value null == ∞)
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
domainate.py -i genomes.gb -r pfam.hmm -o annotated.gb --max_overlap 0.6
summary_report.py -i annotated.gb --json summary.json
enum_report.py -i annotated.gb --by contig --architecture --length --json contigs.ndjson
```

`summary.json` gives dataset-level domain frequencies and (with `--taxonomy`)
taxonomic breakdown; `contigs.ndjson` gives one line per contig with its domain
architecture — filter/aggregate it with any JSON tooling.

## Recipe 2 — build a similarity network and navigate it

```bash
# 1. all-vs-all similarity -> dense matrix
seq_dist.py -i proteins.gb -r proteins.gb --dense dist.hdf5

# 2. (optional) per-node metadata to attach, e.g. domain architecture + taxonomy
enum_report.py -i proteins.gb --by contig --architecture --taxname superkingdom -o meta.tsv

# 3. compact SSN bundle
build_ssn_viewer.py -i dist.hdf5 -o net.ssnv --metadata meta.tsv

# 4. navigate without rendering
ssn_navigator.py -i net.ssnv --mode overview
ssn_navigator.py -i net.ssnv --mode thresholds                 # pick a cut-point
ssn_navigator.py -i net.ssnv --mode clusters --threshold 150   # clusters + sizes
ssn_navigator.py -i net.ssnv --mode cluster --threshold 150 --id 5   # members + metadata distribution
```

The `metadata.tsv` columns (from `enum_report`) become the per-cluster
distributions in `--mode cluster` output: string columns yield the top values
with counts; numeric columns yield min/max/mean/quartiles.

## Recipe 3 — inspect a matrix's cluster structure

```bash
seq_dist.py -i proteins.gb -r proteins.gb --dense dist.hdf5
matrix_report.py -i dist.hdf5 --json -
```

Returns edge-score stats, the number of connected components, and the strongest
threshold split events — enough to decide where to cut the network before
building an SSN.

## Notes on token efficiency

- `ssn_navigator.py` and the report `--json` modes bound their lists with
  `--top_n` (default 50). Raise it deliberately when you need the long tail.
- `ssn_navigator.py --mode cluster` omits the member id list unless you pass
  `--members`; the metadata distribution summary is usually what you want.
- NDJSON (`enum_report`/`hmmer_report`) streams, so you can `head`/filter huge
  reports without materializing the whole thing.
