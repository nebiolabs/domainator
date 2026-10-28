# Domainator for AI agents

Domainator is a suite of ~40 composable CLI tools for genome-neighborhood and
protein-domain analysis. Every tool is also an importable Python function. This
file is the fast orientation for an AI agent driving the suite; see
[docs/agent_workflows.md](docs/agent_workflows.md) for copy-paste recipes and
[README.md](README.md) for the full program table.

## The one convention that trips people up

The `-i` / `--input` flag is **the file being edited**, not the query. Editing
*criteria* come from other flags — usually `-r` / `--reference` (HMMs or
sequences to search *for*). Example:

```bash
domainator domain_search -i database.fasta -r query.hmm -o hits.gb
#                   ^^ searched through      ^^ searched for
```

## Getting machine-readable output (low-token, high-signal)

Prefer these over the HTML reports, which are large and meant for humans:

| Want | Use |
| ---- | --- |
| Per-contig / per-CDS / per-domain table | `domainator enum_report --json_out -` (streaming NDJSON, one object per record) |
| Dataset-level domain/taxonomy summary | `domainator summary_report --json_out -` (one JSON object) |
| HMM profile table | `domainator hmmer_report --json_out -` (NDJSON) |
| Matrix / network structure summary | `domainator matrix_report --json_out -` (one JSON object) |
| Navigate an SSN (clusters, cluster contents, metadata) | `domainator ssn_navigator` (see below) |
| Raw tabular matrices | any comparison tool's `--dense_text_out` (TSV) |

`--json_out -` writes to stdout; give a path to write a file. TSV outputs
(`enum_report`/`hmmer_report` default `-o`, `--dense_text_out` matrices,
`summary_report --domains_tsv_out`) are also directly parseable.

## Navigating a sequence-similarity network without rendering it

`domainator build_ssn --dvb_out` produces a compact `.dvb` bundle (gzip JSON) that the
Domainator Viewer renders for humans. `domainator ssn_navigator` reads
the same bundle and answers targeted questions as compact JSON — no rendering, no
whole-network dump:

```bash
domainator ssn_navigator -i net.dvb -o - --query overview                       # size + metadata columns
domainator ssn_navigator -i net.dvb -o - --query thresholds                     # suggested similarity cut-points
domainator ssn_navigator -i net.dvb -o - --query clusters --threshold 50        # clusters (+sizes) at a threshold
domainator ssn_navigator -i net.dvb -o - --query cluster --threshold 50 --cluster 5  # one cluster: members + metadata distributions
domainator ssn_navigator -i net.dvb -o - --query cluster --threshold 50 --node ABC123   # the cluster a node falls in
domainator ssn_navigator -i net.dvb -o - --query node --node ABC123             # a node's metadata + membership across thresholds
```

Higher `--threshold` = finer clusters (more splitting). Omit `--threshold` (or
pass `inf`) for one cluster per connected component. `--top_n` bounds list sizes;
`--members` includes full member id lists.

## Tool taxonomy (from README.md)

- **editors** — input format == output format, so they chain/stream (GenBank in,
  GenBank out): `domainate`, `domain_search`, `filter_domains`, `select_by_*`,
  `extract_*`, `trim_contigs`, `rotate_contigs`, `color_genbank`, `clean_sequences`, `sort_contigs`.
- **reports** — TSV/HTML/JSON summaries: `summary_report`, `enum_report`,
  `hmmer_report`, `matrix_report`, `ssn_navigator`.
- **comparison** — score/distance matrices: `seq_dist`, `kmer_dist`, `structure_dist`,
  `compare_contigs`, `hmmer_compare`, `transform_matrix`.
- **plotting** — XGMML/SVG/HTML/trees: `build_ssn`, `viewer`,
  `build_projection`, `build_tree`, `plot_contigs`, `color_table_to_legend`.

## Streaming and piping

Editors accept `-i -` (stdin) and `-o -` (stdout) for GenBank and can be
shell-piped: `domainator domainate ... -o - | domainator filter_domains -i - -o out.gb`. Tools
needing multiple passes (`--pad`, most reports/comparisons) cannot stream.

## Invoking programmatically

Every tool accepts `--config=params.json` (jsonargparse), so parameters can be
supplied as a JSON file instead of flags. Each tool's core function is also
importable: `from domainator.enum_report import enum_report`.
