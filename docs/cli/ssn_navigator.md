[CLI reference](README.md)
# domainator ssn_navigator

Query a similarity network bundle: clusters, members, and metadata.

```
Navigate the network in a Domainator Viewer bundle (build_ssn --dvb_out) and emit compact JSON summaries.

This is a read-only, agent-friendly companion to the interactive HTML SSN
viewer. Instead of rendering the whole network, it answers targeted questions
about the sequence-similarity network at a chosen similarity threshold and
returns compact JSON: which clusters exist, what a cluster contains, and how a
cluster's metadata is distributed. Higher thresholds give finer clusters (more
splitting), matching the viewer's slider: ``--threshold inf`` is the "∞" stop at
its fine end, where every node is its own cluster, and omitting ``--threshold``
uses the floor at the coarse end, one cluster per connected component.

The cut-points offered by ``thresholds`` and ``node`` are derived from the bundle's
merge order, not read out of it, so ``--max_merge_events`` controls how many are offered
(``0`` for one per split event) whatever the bundle was built with.

Queries (--query):
  overview    network size, metadata columns, component count, defaults.
  thresholds  suggested similarity cut-points with the resulting cluster counts.
  clusters    list the clusters active at --threshold (with sizes).
  cluster     detail one cluster (--cluster or --node): size, members, metadata.
  node        a node's metadata and the cluster it falls in at each cut-point.
```

Give at least one of these outputs: `--json_out`. `-` as an input or output means standard input or output, where the option says it is allowed.

| Option | Value | Default | Description |
|---|---|---|---|
| `-v`, `--verbose` |  |  | show debug output, and a traceback for errors. |
| `-q`, `--quiet` |  |  | show errors only: no progress messages or warnings. |
| `-i`, `--input` | FILE | required | A Domainator Viewer bundle (.dvb) with a network, from build_ssn --dvb_out or the viewer's Save session. Older network bundles (.dsnv) work too. |
| `-o`, `--json_out` | FILE |  | Write the JSON answer to this file. Use - for standard output. |
| `--query` | {overview,thresholds,clusters,cluster,node} | overview | Which navigation query to run. |
| `--threshold` | STRING |  | Similarity threshold to cut the network at (higher = finer clusters), matching the viewer's slider. Pass 'inf' for the slider's fine end, where every node is its own cluster. Omit it for the coarse end: one cluster per connected component. |
| `--cluster` | INTEGER |  | Cluster id to describe (see --query clusters). For --query cluster. |
| `--node` | STRING |  | Node id/name to describe or locate. For --query cluster/node. |
| `--min_size` | INTEGER | 1 | Only report clusters with at least this many members (--query clusters). |
| `--top_n` | INTEGER | 50 | Bound list output: at most this many clusters, category values, etc. |
| `--members` |  |  | Include the full member id list in --query cluster output. |
| `--max_merge_events` | INTEGER | 500 | How many split events the --query thresholds/node cut-points are drawn from: the strongest this many by impact, plus a few so that no 5% band of the threshold range with an event to offer is left without one. 0 offers a cut-point at every split event. The series is derived from the bundle's merge order, so any value works on any bundle. |

Every command also takes `--config FILE` (options from a JSON or YAML file, keyed by option name) and `--print_config`.
