# Domainator for AI agents

Domainator is a suite of ~40 composable CLI tools for genome-neighborhood and
protein-domain analysis. Every tool is also an importable Python function. This
file is the fast orientation for an AI agent driving the suite; see
[docs/agent_workflows.md](docs/agent_workflows.md) for copy-paste recipes and
[README.md](README.md) for the full program table.

## The one convention that trips people up

The `-i` / `--input` flag is **the file being edited**, not the query. Editing
*criteria* come from other flags — usually `-r` / `--references` (HMMs or
sequences to search *for*). Example:

```bash
domain_search.py -i database.fasta -r query.hmm -o hits.gb
#                   ^^ searched through      ^^ searched for
```

## Getting machine-readable output (low-token, high-signal)

Prefer these over the HTML reports, which are large and meant for humans:

| Want | Use |
| ---- | --- |
| Per-contig / per-CDS / per-domain table | `enum_report.py --json -` (streaming NDJSON, one object per record) |
| Dataset-level domain/taxonomy summary | `summary_report.py --json -` (one JSON object) |
| HMM profile table | `hmmer_report.py --json -` (NDJSON) |
| Matrix / network structure summary | `matrix_report.py --json -` (one JSON object) |
| Navigate an SSN (clusters, cluster contents, metadata) | `ssn_navigator.py` (see below) |
| Raw tabular matrices | any comparison tool's `--dense_text` (TSV) |

`--json -` writes to stdout; give a path to write a file. TSV outputs
(`enum_report`/`hmmer_report` default `-o`, `--dense_text` matrices,
`summary_report --domains_table`) are also directly parseable.

## Navigating a sequence-similarity network without rendering it

`build_ssn_viewer.py` produces a compact `.ssnv` bundle (gzip JSON) that the HTML
viewer renders for humans. `ssn_navigator.py` reads the same bundle and answers
targeted questions as compact JSON — no rendering, no whole-network dump:

```bash
ssn_navigator.py -i net.ssnv --mode overview                       # size + metadata columns
ssn_navigator.py -i net.ssnv --mode thresholds                     # suggested similarity cut-points
ssn_navigator.py -i net.ssnv --mode clusters --threshold 50        # clusters (+sizes) at a threshold
ssn_navigator.py -i net.ssnv --mode cluster --threshold 50 --id 5  # one cluster: members + metadata distributions
ssn_navigator.py -i net.ssnv --mode cluster --threshold 50 --node ABC123   # the cluster a node falls in
ssn_navigator.py -i net.ssnv --mode node --node ABC123             # a node's metadata + membership across thresholds
```

Higher `--threshold` = finer clusters (more splitting). Omit `--threshold` (or
pass `inf`) for one cluster per connected component. `--top_n` bounds list sizes;
`--members` includes full member id lists.

## Tool taxonomy (from README.md)

- **editors** — input format == output format, so they chain/stream (GenBank in,
  GenBank out): `domainate`, `domain_search`, `filter_domains`, `select_by_*`,
  `extract_*`, `trim_contigs`, `color_genbank`, `clean_sequences`, `sort_contigs`.
- **reports** — TSV/HTML/JSON summaries: `summary_report`, `enum_report`,
  `hmmer_report`, `matrix_report`, `ssn_navigator`.
- **comparison** — score/distance matrices: `seq_dist`, `kmer_dist`,
  `compare_contigs`, `hmmer_compare`, `transform_matrix`.
- **plotting** — XGMML/SVG/HTML/trees: `build_ssn`, `build_ssn_viewer`,
  `build_projection`, `build_tree`, `plot_contigs`, `color_table_to_legend`.

## Streaming and piping

Editors accept `-i -` (stdin) and `-o -` (stdout) for GenBank and can be
shell-piped: `domainate.py ... -o - | filter_domains.py -i - -o out.gb`. Tools
needing multiple passes (`--pad`, most reports/comparisons) cannot stream.

## Invoking programmatically

Every tool accepts `--config=params.json` (jsonargparse), so parameters can be
supplied as a JSON file instead of flags. Each tool's core function is also
importable: `from domainator.enum_report import enum_report`.
