""" Reads a DataMatrix and generates a report suitable for helping with selecting edge score thresholds for similarity networks
"""
from jsonargparse import ActionConfigFile
import sys
import math
from time import perf_counter
from typing import Callable, List, Optional, TextIO

from domainator.cli.common import add_input, add_output, check_outputs, declare_outputs, tool_parser
from domainator.outputs import output_file
from domainator.data_matrix import DataMatrix
from domainator.ssn_edges import MaxTree, build_symmetric_neighbor_rankings, mst_knn_edge_counts_by_threshold
from contextlib import ExitStack
from domainator.text_plot import text_histogram
import numpy as np
import json
from domainator import report_html
from domainator.ssn_hierarchy import (
    DEFAULT_MAX_MERGE_EVENTS,
    MERGE_IMPACT_CHOICES,
    MERGE_IMPACT_MIN_CHILD,
    component_size_summary_by_threshold,
    filter_merge_event_rows,
    format_merge_impact_metric,
    merge_event_moving_sum,
    merge_event_table_rows,
    threshold_merge_event_rows,
    threshold_slider_stops,
)

KNN_PROJECTION_MIN_K = 1
KNN_PROJECTION_MAX_K = 25


def _get_mst_knn_report_config(tree):
    max_neighbors = max(tree.n_nodes - 1, KNN_PROJECTION_MIN_K)
    max_k = min(KNN_PROJECTION_MAX_K, max_neighbors)
    default_k = min(5, max_k)
    return {
        "min_k": KNN_PROJECTION_MIN_K,
        "max_k": max_k,
        "default_k": default_k,
    }


def _json_ready(value):
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    if isinstance(value, np.ndarray):
        return _json_ready(value.tolist())
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        value = float(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            return None
        return value
    return value


def _merge_event_table_rows(component_summary, max_items=25):
    return merge_event_table_rows(component_summary, max_items=max_items)


def _record_stage_timing(stage_timings, stage_name, start_time, **metrics):
    elapsed = perf_counter() - start_time
    stage_timings.append((stage_name, elapsed, metrics))
    return elapsed


def _emit_stage_timings(stage_timings, out_handle=None):
    if out_handle is None:
        out_handle = sys.stderr

    print("matrix_report stage timings:", file=out_handle)
    for stage_name, elapsed, metrics in stage_timings:
        suffix = ""
        if len(metrics) > 0:
            suffix = " (" + ", ".join(f"{key}={value}" for key, value in metrics.items()) + ")"
        print(f"  {stage_name}: {elapsed:.3f}s{suffix}", file=out_handle)


def _make_progress_callback(out_handle=None):
    if out_handle is None:
        out_handle = sys.stderr

    start_time = perf_counter()

    def emit(message):
        elapsed = perf_counter() - start_time
        print(f"[matrix_report +{elapsed:.1f}s] {message}", file=out_handle, flush=True)

    return emit


def _estimate_mst_knn_counts_bytes(n_thresholds: int, max_k: int) -> int:
    if max_k < KNN_PROJECTION_MIN_K:
        return 0
    return int(n_thresholds) * int(max_k - KNN_PROJECTION_MIN_K + 1) * np.dtype(int).itemsize

# The version of the --json_out document; see docs/file_formats.md.
JSON_SCHEMA_VERSION = 1


class SummaryJSONWriter():
    """Emit a compact JSON summary of the matrix/network (agent-friendly).

    Includes edge-score stats, the connected-component count, and the bounded
    strongest-split-event series, but omits the bulky per-node arrays that the
    HTML report embeds for the browser.
    """
    def __init__(self, out_handle, max_merge_events=DEFAULT_MAX_MERGE_EVENTS):
        self.out_handle = out_handle
        self.payload = {"schema_version": JSON_SCHEMA_VERSION}

    def write_header(self, tree, edge_scores):
        n = len(edge_scores)
        self.payload["nodes"] = tree.n_nodes
        self.payload["non_zero_edges"] = n
        if n > 0:
            median = (float(edge_scores[n // 2]) if n % 2 == 1
                      else (float(edge_scores[n // 2 - 1]) + float(edge_scores[n // 2])) / 2)
            self.payload["edge_scores"] = {
                "mean": float(np.mean(edge_scores)),
                "median": median,
                "min": float(edge_scores[-1]),
                "max": float(edge_scores[0]),
            }
        else:
            self.payload["edge_scores"] = None

    def write_plots(self, tree, edge_scores, mst_knn_config, mst_knn_counts, component_summary, merge_impact_metric, max_merge_events=DEFAULT_MAX_MERGE_EVENTS):
        self.payload["connected_components"] = tree.n_nodes - len(tree.mst_edges)
        self.payload["merge_impact_metric"] = merge_impact_metric
        if component_summary is not None and len(component_summary) > 1:
            all_merge_event_rows = threshold_merge_event_rows(component_summary)
            self.payload["split_events"] = filter_merge_event_rows(
                all_merge_event_rows, max_merge_events=max_merge_events
            )
            # split_events is a capped selection; this is what it was selected from.
            self.payload["split_event_total"] = len(all_merge_event_rows)
            self.payload["max_merge_events"] = int(max_merge_events)
        else:
            self.payload["split_events"] = []
            self.payload["split_event_total"] = 0
            self.payload["max_merge_events"] = int(max_merge_events)

    def write_footer(self):
        json.dump(_json_ready(self.payload), self.out_handle, separators=(",", ":"))
        self.out_handle.write("\n")

class SummaryTextWriter():
    def __init__(self, out_handle):
        self.out_handle = out_handle

    def write_header(self, tree, edge_scores):
        non_zero_values = edge_scores
        n_nodes = tree.n_nodes
        n = len(non_zero_values)
        # edge_scores is SortedUndirectedEdges.score — already sorted descending,
        # so min/max/median are O(1) index accesses and mean avoids a Python-level loop.
        mean_val   = float(np.mean(non_zero_values)) if n > 0 else 0
        median_val = float(non_zero_values[n // 2]) if n % 2 == 1 else (float(non_zero_values[n // 2 - 1]) + float(non_zero_values[n // 2])) / 2 if n > 0 else 0
        min_val    = float(non_zero_values[-1]) if n > 0 else 0
        max_val    = float(non_zero_values[0])  if n > 0 else 0

        print(f"""Matrix Report
Nodes: {n_nodes}
Non-zero edges: {n}
Mean: {mean_val:.1f}
Median: {median_val:.1f}
Min: {min_val:.1f}
Max: {max_val:.1f}
        """, file=self.out_handle)
        
    def write_plots(self, tree, edge_scores, mst_knn_config, mst_knn_counts, component_summary, merge_impact_metric, max_merge_events=DEFAULT_MAX_MERGE_EVENTS):
        non_zero_values = edge_scores
        
        # Histogram of all triangular values
        if len(non_zero_values) > 0:
            hist_str = text_histogram(non_zero_values, title="Edge counts by score")

        else:
            hist_str = ""
        print(f"""
{hist_str}""", file=self.out_handle)
        
        # MST analysis summary - print key values for text output
        cluster_by_thresh = tree.cluster_count_by_threshold
        if len(cluster_by_thresh) > 1:
            print(f"\nCluster count at key thresholds:", file=self.out_handle)
            # Show first, middle, and last few points
            show_indices = [0, 1, 2, len(cluster_by_thresh)//2, -3, -2, -1]
            seen = set()
            for i in show_indices:
                if 0 <= i < len(cluster_by_thresh) and i not in seen:
                    seen.add(i)
                    thresh, count = cluster_by_thresh[i]
                    if thresh == float('inf'):
                        print(f"  Threshold: ∞ → Clusters: {int(count)}", file=self.out_handle)
                    else:
                        print(f"  Threshold: {thresh:.2f} → Clusters: {int(count)}", file=self.out_handle)
        
        cluster_by_edges = tree.cluster_count_by_edge_count
        if len(cluster_by_edges) > 1:
            print(f"\nCluster count at key edge counts:", file=self.out_handle)
            show_indices = [0, 1, 2, len(cluster_by_edges)//2, -3, -2, -1]
            seen = set()
            for i in show_indices:
                if 0 <= i < len(cluster_by_edges) and i not in seen:
                    seen.add(i)
                    edges, count = cluster_by_edges[i]
                    print(f"  Edges: {int(edges)} → Clusters: {int(count)}", file=self.out_handle)
        
        edges_by_thresh = tree.edges_by_threshold
        if len(edges_by_thresh) > 0:
            print(f"\nEdges kept at key thresholds (equivalent to build_ssn --threshold <threshold>):", file=self.out_handle)
            show_indices = [0, 1, 2, len(edges_by_thresh)//2, -3, -2, -1]
            seen = set()
            for i in show_indices:
                if 0 <= i < len(edges_by_thresh) and i not in seen:
                    seen.add(i)
                    edges, thresh = edges_by_thresh[i]
                    print(f"  Threshold: {thresh:.2f} → Edges: {int(edges)}", file=self.out_handle)

        if mst_knn_counts is not None and len(mst_knn_counts) > 0:
            report_k = mst_knn_config["default_k"]
            print(f"\nProjected MST + kNN edge counts (k={report_k}):", file=self.out_handle)
            show_indices = [0, 1, 2, len(mst_knn_counts)//2, -3, -2, -1]
            seen = set()
            for i in show_indices:
                if 0 <= i < len(mst_knn_counts) and i not in seen:
                    seen.add(i)
                    threshold = edges_by_thresh[i, 1]
                    projected_edges = mst_knn_counts[i, report_k - mst_knn_config["min_k"]]
                    print(f"  Threshold: {threshold:.2f} → MST + kNN edges: {int(projected_edges)}", file=self.out_handle)

        if component_summary is not None and len(component_summary) > 1:
            merge_event_rows = _merge_event_table_rows(component_summary)
            print(f"\nStrongest MST split events (metric={format_merge_impact_metric(merge_impact_metric)}):", file=self.out_handle)
            if len(merge_event_rows) == 0:
                print("  No split events found.", file=self.out_handle)
            else:
                print(
                    f"  {'From':>8}  {'To':>8}  {'Impact':>8}",
                    file=self.out_handle,
                )
                for merge_row in merge_event_rows:
                    print(
                        f"  {merge_row['threshold_from']:>8}  {merge_row['threshold_to']:>8}  {merge_row['merge_impact']:>8}",
                        file=self.out_handle,
                    )
    
    def write_footer(self):
        pass

class SummaryHTMLWriter(SummaryJSONWriter):
    """The interactive page (web/reports/matrix_report.*): the --json_out document, and what
    the page needs to show the network at any threshold: the maximum spanning tree, the
    threshold tables, the moving sum of split impact, the slider's stops and, with
    --knn_projection, the MST + kNN edge counts."""
    def write_plots(self, tree, edge_scores, mst_knn_config, mst_knn_counts, component_summary, merge_impact_metric, max_merge_events=DEFAULT_MAX_MERGE_EVENTS):
        super().write_plots(tree, edge_scores, mst_knn_config, mst_knn_counts, component_summary, merge_impact_metric, max_merge_events=max_merge_events)
        merge_event_rows = threshold_merge_event_rows(component_summary) if component_summary is not None else []
        self.data = {
            "summary": self.payload,
            "mst_edges": tree.export_for_interactive_viz()["mst_edges"],
            # One row per distinct threshold, the rows of edges_by_threshold and
            # cluster_count_by_threshold, as columns.
            "threshold_table": {
                "threshold": [float(threshold) for threshold in tree.thresholds],
                "edges": [int(edges) for edges in tree.edges_by_threshold[:, 0]],
                "clusters": [int(clusters) for clusters in tree.cluster_count_by_threshold[:, 1]],
            },
            # From the uncapped rows: filtering keeps the top rows by merge_impact, so a moving
            # sum over the filtered series would drop exactly the small events it exists to show.
            "merge_moving_sum": merge_event_moving_sum(merge_event_rows),
            # The page reads each cut's edge counts out of threshold_table, hence the lookup.
            "slider_stops": threshold_slider_stops(
                self.payload["split_events"], tree=tree, threshold_index_lookup=tree.threshold_row_index,
            ),
            "knn_projection": None if mst_knn_counts is None else {
                "min_k": mst_knn_config["min_k"],
                "max_k": mst_knn_config["max_k"],
                "default_k": mst_knn_config["default_k"],
                "counts": mst_knn_counts.tolist(),
            },
        }

    def write_footer(self):
        self.out_handle.write(report_html.render_report(
            "matrix_report", "Matrix Report", "matrix_report", _json_ready(self.data), encoding="gzip+base64",
        ))

def matrix_report(matrix:DataMatrix, *, text_out_handle: Optional[TextIO] = None, html_out_handle: Optional[TextIO] = None,
                  knn_projection: bool = False,
                  merge_impact_metric: str = MERGE_IMPACT_MIN_CHILD,
                  stage_timings: Optional[List[tuple]] = None,
                  progress_callback: Optional[Callable[[str], None]] = None,
                  max_merge_events: int = DEFAULT_MAX_MERGE_EVENTS, json_out_handle: Optional[TextIO] = None) -> None:
    """Write a report on ``matrix``, a similarity matrix, to any of the text, HTML and JSON
    handles (each can be None): its size and value distribution, and for a symmetric matrix,
    how its clusters merge as the threshold falls.

    Args:
        knn_projection: also project the edge counts of ``build_ssn --mst --knn K`` for each K.
        merge_impact_metric: how split events are measured, "min_child" or "product".
        max_merge_events: how many of the strongest merge events to embed (0: all).
        stage_timings: a list to append each stage's (stage, seconds, metrics) timing to.
        progress_callback: called with the name of each stage as it starts.
    """
    if merge_impact_metric not in MERGE_IMPACT_CHOICES:
        raise ValueError(f"merge_impact_metric must be one of {sorted(MERGE_IMPACT_CHOICES)}")
    if max_merge_events < 0:
        raise ValueError("max_merge_events must be >= 0")

    longform_outputs = list()
    if text_out_handle is not None:
        longform_outputs.append(SummaryTextWriter(text_out_handle))
    if html_out_handle is not None:
        longform_outputs.append(SummaryHTMLWriter(html_out_handle))
    if json_out_handle is not None:
        longform_outputs.append(SummaryJSONWriter(json_out_handle, max_merge_events=max_merge_events))

    if stage_timings is None:
        stage_timings = []
    
    if progress_callback is not None:
        progress_callback("building edge table")
    stage_start = perf_counter()
    matrix.require_similarity("matrix_report")
    edge_table = matrix.sorted_undirected_edges(skip_zeros=True, agg=max)
    del matrix  # free O(n²) data array; edge_table holds all we need
    _record_stage_timing(stage_timings, "build_edge_table", stage_start, n_edges=len(edge_table), n_nodes=edge_table.n_nodes)

    if progress_callback is not None:
        progress_callback("building maximum spanning tree")
    stage_start = perf_counter()
    tree = MaxTree(edge_table)
    _record_stage_timing(stage_timings, "build_max_tree", stage_start, n_mst_edges=len(tree.mst_edges))
    component_summary = component_size_summary_by_threshold(tree, merge_impact_metric=merge_impact_metric)
    mst_knn_config = None
    mst_knn_counts = None
    neighbor_rankings = None
    if knn_projection:
        if progress_callback is not None:
            progress_callback("building kNN projection neighbor rankings")
        stage_start = perf_counter()
        mst_knn_config = _get_mst_knn_report_config(tree)
        neighbor_rankings = build_symmetric_neighbor_rankings(edge_table, max_k=mst_knn_config["max_k"])
        _record_stage_timing(
            stage_timings,
            "build_mst_knn_neighbor_rankings",
            stage_start,
            max_k=mst_knn_config["max_k"],
            kept_directed_edges=len(neighbor_rankings.target),
            estimated_counts_bytes=_estimate_mst_knn_counts_bytes(len(tree.thresholds), mst_knn_config["max_k"]),
        )

    if knn_projection:
        if progress_callback is not None:
            progress_callback("building kNN projection threshold counts")
        stage_start = perf_counter()
        mst_knn_counts = mst_knn_edge_counts_by_threshold(None, tree, mst_knn_config["max_k"], neighbor_rankings=neighbor_rankings)
        _record_stage_timing(
            stage_timings,
            "build_mst_knn_counts",
            stage_start,
            thresholds=len(tree.thresholds),
            max_k=mst_knn_config["max_k"],
            output_shape=mst_knn_counts.shape,
            output_bytes=mst_knn_counts.nbytes,
        )

    if progress_callback is not None:
        progress_callback("rendering outputs")
    stage_start = perf_counter()
    for output in longform_outputs:
        output.write_header(tree, edge_table.score)
        output.write_plots(tree, edge_table.score, mst_knn_config, mst_knn_counts, component_summary, merge_impact_metric, max_merge_events=max_merge_events)
        output.write_footer()
    _record_stage_timing(stage_timings, "render_outputs", stage_start, outputs=len(longform_outputs))

def build_parser():
    parser = tool_parser(__doc__)

    add_input(parser, "name of input matrix file.", multiple=False, stdin=False)

    add_output(parser, "text", "Write the report as text to this file.", primary=True)
    add_output(parser, "html", "Write the report as html to this file.")
    add_output(parser, "json", "Write a compact JSON summary (edge stats, connected components, split-event series) to this file.")
    parser.add_argument('--knn_projection', action='store_true', default=False,
                        help=f"Project how many edges 'build_ssn --threshold T --mst --knn K' would keep, for every threshold T and K from {KNN_PROJECTION_MIN_K} to {KNN_PROJECTION_MAX_K}, and add a K slider to the HTML report. Disabled by default to reduce memory use.")
    parser.add_argument('--merge_impact_metric', choices=list(MERGE_IMPACT_CHOICES), default=MERGE_IMPACT_MIN_CHILD,
                        help="Metric used for split-event plots and tables: 'min_child' emphasizes the smaller cluster breaking away, while 'product' emphasizes balanced large splits.")
    parser.add_argument('--max_merge_events', type=int, default=DEFAULT_MAX_MERGE_EVENTS,
                        help="Maximum number of strongest merge events to embed in the interactive HTML threshold slider and split plot. A few more may be added on top: after the strongest are taken, the strongest event in each otherwise-empty 5%% band of the threshold axis is added back, so the plot never leaves a stretch of its own axis blank. Use 0 to include all merge events.")
    parser.add_argument('--progress', action='store_true', default=False,
                        help="Print live progress updates to stderr during long-running matrix_report stages.")
    parser.add_argument('--profile_stages', action='store_true', default=False,
                        help="Print per-stage runtime timings to stderr for profiling large matrix_report runs.")
    
    parser.add_argument('--config', action=ActionConfigFile)
    declare_outputs(parser, ["text_out", "html_out", "json_out"])
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)

    check_outputs(parser, params)

    with ExitStack() as output_stack:
        def open_output(path):
            return None if path is None else output_stack.enter_context(output_file(path))
        out = open_output(params.text_out)
        html_out_handle = open_output(params.html_out)
        json_out_handle = open_output(params.json_out)

        progress_callback = None
        if params.progress or params.profile_stages:
            progress_callback = _make_progress_callback()

        stage_timings = []
        if progress_callback is not None:
            progress_callback(f"loading matrix from {params.input}")
        stage_start = perf_counter()
        matrix = DataMatrix.from_file(params.input)
        if params.profile_stages:
            _record_stage_timing(stage_timings, "load_matrix", stage_start, shape=matrix.shape, matrix_type=type(matrix).__name__)

        ### Run
        if out is not None or html_out_handle is not None or json_out_handle is not None:
            stage_start = perf_counter()
            matrix_report(
                matrix,
                text_out_handle=out,
                html_out_handle=html_out_handle,
                knn_projection=params.knn_projection,
                merge_impact_metric=params.merge_impact_metric,
                stage_timings=stage_timings,
                progress_callback=progress_callback,
                max_merge_events=params.max_merge_events,
                json_out_handle=json_out_handle,
            )
            if params.profile_stages:
                _record_stage_timing(stage_timings, "matrix_report_total", stage_start)
                _emit_stage_timings(stage_timings)
