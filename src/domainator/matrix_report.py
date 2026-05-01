""" Reads a DataMatrix and generates a report suitable for helping with selecting edge score thresholds for similarity networks
"""
from jsonargparse import ArgumentParser, ActionConfigFile
import sys
import math
from time import perf_counter
import base64
import gzip

from domainator.data_matrix import DataMatrix, MaxTree, build_symmetric_neighbor_rankings, mst_knn_edge_counts_by_threshold
from domainator import __version__, RawAndDefaultsFormatter
from bashplotlib.histogram import plot_hist
import io
from contextlib import ExitStack, redirect_stdout
from domainator.summary_report import histplot_base64
import statistics
import numpy as np
import json
from domainator.ssn_hierarchy import (
    DEFAULT_MAX_MERGE_EVENTS,
    MERGE_IMPACT_CHOICES,
    MERGE_IMPACT_MIN_CHILD,
    MERGE_IMPACT_PRODUCT,
    component_size_summary_by_threshold,
    filter_merge_event_rows,
    format_merge_impact_metric,
    merge_event_table_rows,
    threshold_merge_event_rows,
)

MST_KNN_MIN_K = 2
MST_KNN_MAX_K = 25
SLIDER_POSITION_SCALE = 10_000
SLIDER_RIGHT_STOP_PADDING_FRACTION = 0.05
SLIDER_LEFT_STOP_PADDING_FRACTION = 0.05
SPLIT_PLOT_MARGIN_TOP = 72
SPLIT_PLOT_MARGIN_BOTTOM = 50
SPLIT_PLOT_MARGIN_LEFT = 60
SPLIT_PLOT_MARGIN_RIGHT = 20
SLIDER_PANEL_PADDING = 15
THRESHOLD_SLIDER_ALIGN_LEFT = max(0, SPLIT_PLOT_MARGIN_LEFT - SLIDER_PANEL_PADDING)
THRESHOLD_SLIDER_ALIGN_RIGHT = max(0, SPLIT_PLOT_MARGIN_RIGHT - SLIDER_PANEL_PADDING)


def _get_mst_knn_report_config(tree):
    max_neighbors = max(tree.n_nodes - 1, MST_KNN_MIN_K)
    max_k = min(MST_KNN_MAX_K, max_neighbors)
    default_k = min(5, max_k)
    return {
        "min_k": MST_KNN_MIN_K,
        "max_k": max_k,
        "default_k": default_k,
    }


def _slider_stop_rows(merge_event_rows):
    stops = [{
        "edge_index": -1,
        "threshold_label": "∞",
        "threshold_value": None,
        "slider_position": 0,
    }]
    for merge_row in merge_event_rows:
        stops.append({
            "edge_index": int(merge_row["edge_index"]),
            "threshold_label": merge_row["threshold_to"],
            "threshold_value": merge_row["threshold_value"],
            "slider_position": 0,
        })

    if len(stops) <= 1:
        return stops

    finite_stops = stops[1:]
    finite_thresholds = [float(stop["threshold_value"]) for stop in finite_stops]
    min_threshold = min(finite_thresholds)
    max_threshold = max(finite_thresholds)
    min_finite_slider_position = max(1, int(round(SLIDER_POSITION_SCALE * SLIDER_RIGHT_STOP_PADDING_FRACTION)))
    max_finite_slider_position = max(
        min_finite_slider_position + 1,
        int(round(SLIDER_POSITION_SCALE * (1.0 - SLIDER_LEFT_STOP_PADDING_FRACTION))),
    )

    if math.isclose(max_threshold, min_threshold):
        if len(finite_stops) == 1:
            finite_stops[0]["slider_position"] = max_finite_slider_position
            return stops
        for stop_idx, stop in enumerate(finite_stops, start=1):
            stop["slider_position"] = int(round(
                min_finite_slider_position
                + stop_idx * (max_finite_slider_position - min_finite_slider_position) / len(finite_stops)
            ))
        finite_stops[-1]["slider_position"] = max_finite_slider_position
        return stops

    previous_position = 0
    for stop in finite_stops:
        normalized_position = (max_threshold - float(stop["threshold_value"])) / (max_threshold - min_threshold)
        slider_position = min_finite_slider_position + int(round(
            normalized_position * (max_finite_slider_position - min_finite_slider_position)
        ))
        slider_position = max(previous_position + 1, slider_position)
        stop["slider_position"] = slider_position
        previous_position = slider_position

    finite_stops[-1]["slider_position"] = max_finite_slider_position
    return stops


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


def _compressed_json_base64(payload):
    json_bytes = json.dumps(_json_ready(payload), separators=(",", ":")).encode("utf-8")
    return base64.b64encode(gzip.compress(json_bytes, mtime=0)).decode("ascii")


def _merge_event_table_rows(component_summary, max_items=25):
    return merge_event_table_rows(component_summary, max_items=max_items)


def _interactive_report_payload(tree, mst_knn_config, mst_knn_counts, component_summary, max_merge_events=DEFAULT_MAX_MERGE_EVENTS):
    merge_event_rows = threshold_merge_event_rows(component_summary) if component_summary is not None else []
    filtered_merge_event_rows = filter_merge_event_rows(merge_event_rows, max_merge_events=max_merge_events)
    return {
        "viz_data": tree.export_for_interactive_viz(),
        "cluster_by_thresh": tree.cluster_count_by_threshold,
        "cluster_by_edges": tree.cluster_count_by_edge_count,
        "edges_by_thresh": tree.edges_by_threshold,
        "has_mst_knn": mst_knn_counts is not None,
        "mst_knn_counts": mst_knn_counts if mst_knn_counts is not None else [],
        "mst_knn_min_k": mst_knn_config['min_k'] if mst_knn_config is not None else 0,
        "merge_event_series": filtered_merge_event_rows,
        "slider_stops": _slider_stop_rows(filtered_merge_event_rows),
    }


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
    if max_k < MST_KNN_MIN_K:
        return 0
    return int(n_thresholds) * int(max_k - MST_KNN_MIN_K + 1) * np.dtype(int).itemsize

class SummaryTextWriter():
    def __init__(self, out_handle): 
        self.out_handle = out_handle
    
    def write_header(self, tree, edge_scores):
        non_zero_values = edge_scores
        n_nodes = tree.n_nodes
        
        print(f"""Matrix Report
Nodes: {n_nodes}
Non-zero edges: {len(non_zero_values)}
Mean: {statistics.mean(non_zero_values) if len(non_zero_values) > 0 else 0:.1f}
Median: {statistics.median(non_zero_values) if len(non_zero_values) > 0 else 0:.1f}
Min: {min(non_zero_values) if len(non_zero_values) > 0 else 0:.1f}
Max: {max(non_zero_values) if len(non_zero_values) > 0 else 0:.1f}
        """, file=self.out_handle)
        
    def write_plots(self, tree, edge_scores, mst_knn_config, mst_knn_counts, component_summary, merge_impact_metric, max_merge_events=DEFAULT_MAX_MERGE_EVENTS):
        non_zero_values = edge_scores
        
        # Histogram of all triangular values
        if len(non_zero_values) > 0:
            hist_str = io.StringIO()
            with redirect_stdout(hist_str):
                plot_hist(non_zero_values, bincount=50, title="Edge counts by score", height=20.0, xlab=True, regular=True, showSummary=False)
            hist_str = hist_str.getvalue()
            hist_str = hist_str.replace('[39m', '')

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
            print(f"\nEdge counts at key thresholds:", file=self.out_handle)
            show_indices = [0, 1, 2, len(edges_by_thresh)//2, -3, -2, -1]
            seen = set()
            for i in show_indices:
                if 0 <= i < len(edges_by_thresh) and i not in seen:
                    seen.add(i)
                    edges, thresh = edges_by_thresh[i]
                    print(f"  Threshold: {thresh:.2f} → Cumulative edges: {int(edges)}", file=self.out_handle)

        if mst_knn_counts is not None and len(mst_knn_counts) > 0:
            report_k = mst_knn_config["default_k"]
            print(f"\nProjected MST_KNN edge counts (k={report_k}):", file=self.out_handle)
            show_indices = [0, 1, 2, len(mst_knn_counts)//2, -3, -2, -1]
            seen = set()
            for i in show_indices:
                if 0 <= i < len(mst_knn_counts) and i not in seen:
                    seen.add(i)
                    threshold = edges_by_thresh[i, 1]
                    projected_edges = mst_knn_counts[i, report_k - mst_knn_config["min_k"]]
                    print(f"  Threshold: {threshold:.2f} → MST_KNN edges: {int(projected_edges)}", file=self.out_handle)

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

class SummaryHTMLWriter():
    def __init__(self, out_handle): 
        self.out_handle = out_handle
    
    def write_header(self, tree, edge_scores):
        non_zero_values = edge_scores
        n_nodes = tree.n_nodes
        
        print(f"""<!doctype html><html><head><meta charset="UTF-8" /><title>Matrix Report</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/d3@7"></script>
<style>
    body {{
        font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
        margin: 0;
        padding: 20px;
        background: #f5f5f5;
    }}
    .container {{
        max-width: 1400px;
        margin: 0 auto;
        background: white;
        padding: 20px;
        border-radius: 8px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }}
    h1 {{
        margin-top: 0;
        color: #333;
    }}
    table {{
        border-collapse: collapse;
        margin: 10px 0;
    }}
    th, td {{
        padding: 8px 12px;
        text-align: left;
        border-bottom: 1px solid #ddd;
    }}
    th {{
        font-weight: 600;
        color: #666;
    }}
    .dashboard {{
        display: grid;
        grid-template-columns: repeat(2, minmax(0, 1fr));
        gap: 20px;
        margin-top: 20px;
        align-items: start;
    }}
    .chart {{
        border: 1px solid #e0e0e0;
        padding: 15px;
        background: white;
        border-radius: 4px;
        min-height: 400px;
        position: relative;
        min-width: 0;
        overflow: hidden;
    }}
    .chart-with-controls {{
        border: 1px solid #e0e0e0;
        padding: 15px;
        background: white;
        border-radius: 4px;
        min-height: 400px;
        position: relative;
        min-width: 0;
        overflow: hidden;
        grid-column: 1 / -1;
    }}
    .chart-wide {{
        grid-column: 1 / -1;
    }}
    .distribution-panels {{
        margin-top: 15px;
    }}
    #cluster-size-bubble {{
        min-height: 240px;
    }}
    #cluster-size-bubble svg {{
        display: block;
        width: 100%;
        height: 240px;
    }}
    .bubble-title {{
        margin: 0 0 8px 0;
        font-size: 0.95rem;
        color: #444;
    }}
    @media (max-width: 1100px) {{
        .dashboard {{
            grid-template-columns: minmax(0, 1fr);
        }}
        .chart-wide,
        .chart-with-controls {{
            grid-column: auto;
        }}
    }}
    .slider-container {{
        margin: 15px 0 15px 0;
        padding: 15px;
        background: #f9f9f9;
        border-radius: 4px;
    }}
    .slider-container label {{
        display: block;
        margin-bottom: 8px;
        font-weight: 500;
        color: #555;
    }}
    .threshold-slider-track {{
        padding-left: {THRESHOLD_SLIDER_ALIGN_LEFT}px;
        padding-right: {THRESHOLD_SLIDER_ALIGN_RIGHT}px;
    }}
    #threshold-slider {{
        width: 100%;
        height: 6px;
        border-radius: 3px;
        outline: none;
        -webkit-appearance: none;
        direction: rtl;
    }}
    #threshold-slider::-webkit-slider-thumb {{
        -webkit-appearance: none;
        appearance: none;
        width: 18px;
        height: 18px;
        border-radius: 50%;
        background: #d62728;
        cursor: pointer;
    }}
    #threshold-slider::-moz-range-thumb {{
        width: 18px;
        height: 18px;
        border-radius: 50%;
        background: #d62728;
        cursor: pointer;
        border: none;
    }}
    .stats {{
        display: grid;
        grid-template-columns: repeat(2, 1fr);
        gap: 10px;
        margin: 15px 0;
    }}
    .stat-box {{
        padding: 10px;
        background: white;
        border-radius: 4px;
        border-left: 4px solid #d62728;
    }}
    .stat-box strong {{
        display: block;
        color: #666;
        font-size: 0.75em;
        margin-bottom: 3px;
        text-transform: uppercase;
        letter-spacing: 0.5px;
    }}
    .stat-box span {{
        display: block;
        font-size: 1.4em;
        color: #333;
        font-weight: 300;
    }}
    .plotly {{
        width: 100%;
        height: 100%;
    }}
    .js-plotly-plot,
    .plot-container {{
        width: 100% !important;
    }}
</style>
</head><body>
<div class="container">""", 
              file=self.out_handle)
        print(f"""<h1>Matrix Report</h1>""", file=self.out_handle)
        print(f"""<div><h2>Summary</h2>""", file=self.out_handle)
        print(f"""<table>""", file=self.out_handle)
        print(f"""<tr><th>Nodes</th><th>Non-zero edges</th><th>Mean</th><th>Median</th><th>Min</th><th>Max</th></tr>""", file=self.out_handle)
        print(f"""<tr><td>{n_nodes}</td><td>{len(non_zero_values)}</td><td>{statistics.mean(non_zero_values) if len(non_zero_values) > 0 else 0 : .1f}</td><td>{statistics.median(non_zero_values) if len(non_zero_values) > 0 else 0 : .1f}</td><td>{min(non_zero_values) if len(non_zero_values) > 0 else 0 : .1f}</td><td>{max(non_zero_values) if len(non_zero_values) > 0 else 0 : .1f}</td></tr>""", file=self.out_handle)
        print(f"""</table>""", file=self.out_handle)
        print(f"""</div>""", file=self.out_handle)
    
    
    def write_plots(self, tree, edge_scores, mst_knn_config, mst_knn_counts, component_summary, merge_impact_metric, max_merge_events=DEFAULT_MAX_MERGE_EVENTS):
        include_mst_knn = mst_knn_counts is not None
        include_component_summary = component_summary is not None and len(component_summary) > 1
        payload = _interactive_report_payload(tree, mst_knn_config, mst_knn_counts, component_summary, max_merge_events=max_merge_events)
        filtered_merge_event_rows = payload["merge_event_series"]
        slider_stops = payload["slider_stops"]
        compressed_payload = _compressed_json_base64(payload)
        mst_knn_controls = ""
        mst_knn_stat_box = ""
        mst_knn_current_k_expr = "null"
        mst_knn_update_block = ""
        mst_knn_listener_block = ""
        component_signal_chart = ""
        component_signal_plot_block = ""
        if include_component_summary and len(filtered_merge_event_rows) > 0:
            component_signal_chart = """
    <div class=\"chart chart-wide\">
        <div id=\"cluster-discontinuity-by-threshold\"></div>
    </div>"""
            component_signal_plot_block = f"""

        const mergeEventPoints = MERGE_EVENT_SERIES;
        const mergeEventStemX = [];
        const mergeEventStemY = [];
        mergeEventPoints.forEach(d => {{
            mergeEventStemX.push(d.threshold_value, d.threshold_value, null);
            mergeEventStemY.push(0, d.merge_impact, null);
        }});
        Plotly.newPlot('cluster-discontinuity-by-threshold', [
            {{
                x: mergeEventStemX,
                y: mergeEventStemY,
                mode: 'lines',
                name: 'Split event',
                line: {{color: '#72b7b2', width: 1}},
                hoverinfo: 'skip'
            }},
            {{
                x: mergeEventPoints.map(d => d.threshold_value),
                y: mergeEventPoints.map(d => d.merge_impact),
                mode: 'markers',
                name: 'Split size',
                marker: {{color: '#e45756', size: 6, opacity: 0.75}},
                customdata: mergeEventPoints.map(d => [d.threshold_from, d.threshold_to]),
                hovertemplate: 'Threshold: %{{x:.2f}}<br>Smallest new cluster: %{{y:.2f}}<br>From: %{{customdata[0]}}<br>To: %{{customdata[1]}}<extra></extra>'
            }}
        ], {{
            ...chartLayout,
            margin: {{t: {SPLIT_PLOT_MARGIN_TOP}, b: {SPLIT_PLOT_MARGIN_BOTTOM}, l: {SPLIT_PLOT_MARGIN_LEFT}, r: {SPLIT_PLOT_MARGIN_RIGHT}}},
            title: 'Cluster Splits vs Threshold',
            xaxis: {{title: 'Threshold', type: 'linear', autorange: true}},
            yaxis: {{title: 'Size of smallest new cluster'}},
            legend: {{
                orientation: 'h',
                yanchor: 'bottom',
                y: 1.08,
                xanchor: 'left',
                x: 0
            }},
            hovermode: 'closest'
        }}, chartConfig);"""
        if include_mst_knn:
            mst_knn_controls = f"""
            <label for=\"mst-knn-k-slider\">
                MST_KNN k: <span id=\"mst-knn-k-number\">{mst_knn_config['default_k']}</span>
            </label>
            <input type=\"range\" id=\"mst-knn-k-slider\" min=\"{mst_knn_config['min_k']}\" max=\"{mst_knn_config['max_k']}\" value=\"{mst_knn_config['default_k']}\" step=\"1\">"""
            mst_knn_stat_box = """
            <div class=\"stat-box\">
                <strong>Projected MST_KNN Edges</strong>
                <span id=\"mst-knn-edges\">0</span>
            </div>"""
            mst_knn_current_k_expr = "parseInt(document.getElementById('mst-knn-k-slider').value)"
            mst_knn_update_block = """
        let mstKnnEdgeCount = 0;
        if (edgeIndex >= 0 && edgeIndex < MST_KNN_COUNTS.length) {
            mstKnnEdgeCount = MST_KNN_COUNTS[edgeIndex][currentK - MST_KNN_MIN_K];
        }
        document.getElementById('mst-knn-edges').textContent = formatNumber(mstKnnEdgeCount);
        document.getElementById('mst-knn-k-number').textContent = currentK.toString();"""
            mst_knn_listener_block = """
    document.getElementById('mst-knn-k-slider').addEventListener('input', () => {
        clearTimeout(updateTimeout);
        updateTimeout = setTimeout(() => {
            let position = parseInt(document.getElementById('threshold-slider').value);
            updateHistogramFromSliderPosition(position);
        }, 10);
    });"""
        
        print(f"""
<div class="dashboard">
    <div class="chart">
        <div id="cluster-by-threshold"></div>
    </div>
    <div class="chart">
        <div id="cluster-by-edges"></div>
    </div>
    {component_signal_chart}
    <div class="chart-with-controls">
        <div class="slider-container">
            <label for="threshold-slider">
                Threshold: <span id="threshold-number">∞</span>
            </label>
            <div class="threshold-slider-track">
                <input type="range" id="threshold-slider" min="0" max="{SLIDER_POSITION_SCALE if len(slider_stops) > 0 else 0}" value="0" step="1">
            </div>
            {mst_knn_controls}
        </div>
        <div class="stats">
            <div class="stat-box">
                <strong>Number of Clusters</strong>
                <span id="num-clusters">{tree.n_nodes}</span>
            </div>
            <div class="stat-box">
                <strong>Number of Singletons</strong>
                <span id="num-singletons">{tree.n_nodes}</span>
            </div>
            <div class="stat-box">
                <strong>Largest Cluster</strong>
                <span id="largest-cluster">1</span>
            </div>
            <div class="stat-box">
                <strong>Number of Edges</strong>
                <span id="num-edges">0</span>
            </div>
            <div class="stat-box">
                <strong>Avg Cluster Size</strong>
                <span id="avg-cluster-size">1.0</span>
            </div>
            {mst_knn_stat_box}
        </div>
        <div class="distribution-panels">
            <h3 class="bubble-title">100 Largest Clusters</h3>
            <div id="cluster-size-bubble"></div>
        </div>
    </div>
</div>
</div>

<script>
    const EMBEDDED_REPORT_DATA_GZIP_BASE64 = "{compressed_payload}";
    let VIZ_DATA;
    let CLUSTER_BY_THRESH;
    let CLUSTER_BY_EDGES;
    let EDGES_BY_THRESH;
    let HAS_MST_KNN;
    let MST_KNN_COUNTS;
    let MST_KNN_MIN_K;
    let MERGE_EVENT_SERIES;
    let SLIDER_STOPS;
    let CLUSTER_CHECKPOINTS = [];
    const CLUSTER_CHECKPOINT_STRIDE = 50;

    function base64ToBytes(base64Text) {{
        const binaryText = atob(base64Text);
        const bytes = new Uint8Array(binaryText.length);
        for (let i = 0; i < binaryText.length; i++) {{
            bytes[i] = binaryText.charCodeAt(i);
        }}
        return bytes;
    }}

    async function decodeEmbeddedReportData(base64Text) {{
        if (!('DecompressionStream' in window)) {{
            throw new Error('This browser cannot decompress embedded matrix_report data. Please use a current Chromium, Firefox, or Safari browser.');
        }}
        const bytes = base64ToBytes(base64Text);
        const stream = new Response(bytes).body.pipeThrough(new DecompressionStream('gzip'));
        const text = await new Response(stream).text();
        return JSON.parse(text);
    }}

    function installReportData(reportData) {{
        VIZ_DATA = reportData.viz_data;
        CLUSTER_BY_THRESH = reportData.cluster_by_thresh;
        CLUSTER_BY_EDGES = reportData.cluster_by_edges;
        EDGES_BY_THRESH = reportData.edges_by_thresh;
        HAS_MST_KNN = reportData.has_mst_knn;
        MST_KNN_COUNTS = reportData.mst_knn_counts;
        MST_KNN_MIN_K = reportData.mst_knn_min_k;
        MERGE_EVENT_SERIES = reportData.merge_event_series;
        SLIDER_STOPS = reportData.slider_stops;
    }}
    
    // Union-Find implementation for fast cluster computation
    class UnionFind {{
        constructor(n, parent = null, size = null) {{
            this.parent = parent ? parent.slice() : Array.from({{length: n}}, (_, i) => i);
            this.size = size ? size.slice() : Array(n).fill(1);
        }}
        
        find(x) {{
            if (this.parent[x] !== x) {{
                this.parent[x] = this.find(this.parent[x]);
            }}
            return this.parent[x];
        }}
        
        union(x, y) {{
            let rx = this.find(x), ry = this.find(y);
            if (rx !== ry) {{
                this.size[ry] += this.size[rx];
                this.parent[rx] = ry;
            }}
        }}
        
        getClusterSizes() {{
            let seen = new Set();
            let sizes = [];
            for (let i = 0; i < this.parent.length; i++) {{
                let root = this.find(i);
                if (!seen.has(root)) {{
                    seen.add(root);
                    sizes.push(this.size[root]);
                }}
            }}
            return sizes.sort((a, b) => b - a);
        }}

        snapshot(edgeIndex) {{
            return {{
                edgeIndex: edgeIndex,
                parent: this.parent.slice(),
                size: this.size.slice()
            }};
        }}
    }}

    function buildClusterCheckpoints() {{
        CLUSTER_CHECKPOINTS = [];
        let uf = new UnionFind(VIZ_DATA.n_nodes);
        CLUSTER_CHECKPOINTS.push(uf.snapshot(-1));

        let nextEdgeIndex = 0;
        for (let position = 1; position < SLIDER_STOPS.length; position++) {{
            const targetEdgeIndex = SLIDER_STOPS[position].edge_index;
            while (nextEdgeIndex <= targetEdgeIndex && nextEdgeIndex < VIZ_DATA.mst_edges.length) {{
                let [n1, n2, val] = VIZ_DATA.mst_edges[nextEdgeIndex];
                uf.union(n1, n2);
                nextEdgeIndex++;
            }}
            if (position % CLUSTER_CHECKPOINT_STRIDE === 0 || position === SLIDER_STOPS.length - 1) {{
                CLUSTER_CHECKPOINTS.push(uf.snapshot(targetEdgeIndex));
            }}
        }}
    }}

    function checkpointForEdgeIndex(edgeIndex) {{
        let checkpoint = CLUSTER_CHECKPOINTS[0];
        for (let i = 1; i < CLUSTER_CHECKPOINTS.length; i++) {{
            if (CLUSTER_CHECKPOINTS[i].edgeIndex > edgeIndex) {{
                break;
            }}
            checkpoint = CLUSTER_CHECKPOINTS[i];
        }}
        return checkpoint;
    }}
    
    // Compute clusters at specific MST edge index
    function computeClusters(edgeIndex) {{
        if (edgeIndex < 0) {{
            return Array(VIZ_DATA.n_nodes).fill(1);
        }}
        const checkpoint = checkpointForEdgeIndex(edgeIndex);
        let uf = new UnionFind(VIZ_DATA.n_nodes, checkpoint.parent, checkpoint.size);
        for (let i = checkpoint.edgeIndex + 1; i <= edgeIndex && i < VIZ_DATA.mst_edges.length; i++) {{
            let [n1, n2, val] = VIZ_DATA.mst_edges[i];
            uf.union(n1, n2);
        }}
        return uf.getClusterSizes();
    }}
    
    // Format large numbers
    function formatNumber(num) {{
        return num.toLocaleString();
    }}

    function updateClusterBubbleChart(sizes) {{
        const bubbleRoot = d3.select('#cluster-size-bubble');
        bubbleRoot.selectAll('*').remove();

        const topSizes = sizes.slice().sort((a, b) => b - a).slice(0, 100);
        const width = Math.max(320, bubbleRoot.node()?.clientWidth || 320);
        const height = 240;

        if (topSizes.length === 0) {{
            bubbleRoot.append('div')
                .style('color', '#666')
                .style('font-size', '0.95rem')
                .text('No clusters to display.');
            return;
        }}

        const colorScale = d3.scaleSequential(d3.interpolateYlOrRd)
            .domain([0, d3.max(topSizes) || 1]);

        const root = d3.hierarchy({{
            children: topSizes.map((size, index) => ({{
                label: 'Cluster ' + (index + 1),
                value: size,
            }}))
        }})
            .sum(d => d.value || 0)
            .sort((a, b) => (b.value || 0) - (a.value || 0));

        d3.pack()
            .size([width, height])
            .padding(4)(root);

        const svg = bubbleRoot.append('svg')
            .attr('viewBox', '0 0 ' + width + ' ' + height)
            .attr('preserveAspectRatio', 'xMidYMid meet');

        const nodes = svg.selectAll('g')
            .data(root.leaves())
            .enter()
            .append('g')
            .attr('transform', d => 'translate(' + d.x + ',' + d.y + ')');

        nodes.append('circle')
            .attr('r', d => d.r)
            .attr('fill', d => colorScale(d.data.value))
            .attr('fill-opacity', 0.85)
            .attr('stroke', '#9a3412')
            .attr('stroke-width', 1);

        nodes.append('title')
            .text(d => d.data.label + ': ' + formatNumber(d.data.value));

        nodes.filter(d => d.r >= 16)
            .append('text')
            .attr('text-anchor', 'middle')
            .attr('dy', '0.35em')
            .style('font-size', d => Math.max(9, Math.min(14, d.r / 2.4)) + 'px')
            .style('font-weight', 600)
            .style('fill', '#3f1d0f')
            .text(d => d.data.value);
    }}
    
    // Initialize static charts
    function initCharts() {{
        const chartConfig = {{displayModeBar: false, responsive: true}};
        const chartLayout = {{
            autosize: true,
            height: 350,
            margin: {{t: 40, b: 50, l: 60, r: 20}}
        }};
        
        // Clusters and edges vs threshold
        Plotly.newPlot('cluster-by-threshold', [{{
            x: CLUSTER_BY_THRESH.slice(1).map(d => d[0]),
            y: CLUSTER_BY_THRESH.slice(1).map(d => d[1]),
            mode: 'lines',
            name: 'Clusters',
            line: {{color: '#1f77b4', width: 2}},
            hovertemplate: 'Threshold: %{{x:.2f}}<br>Clusters: %{{y}}<extra></extra>'
        }}, {{
            x: EDGES_BY_THRESH.map(d => d[1]),
            y: EDGES_BY_THRESH.map(d => d[0]),
            mode: 'lines',
            name: 'Edges',
            yaxis: 'y2',
            line: {{color: '#2ca02c', width: 2}},
            hovertemplate: 'Threshold: %{{x:.2f}}<br>Edges: %{{y}}<extra></extra>'
        }}], {{
            ...chartLayout,
            title: 'Clusters and Edge Count vs Threshold',
            xaxis: {{title: 'Threshold', type: 'log', autorange: true}},
            yaxis: {{title: 'Number of Clusters'}},
            yaxis2: {{
                title: 'Cumulative Edges',
                overlaying: 'y',
                side: 'right'
            }},
            hovermode: 'closest'
        }}, chartConfig);
        
        // Cluster count vs edge count
        Plotly.newPlot('cluster-by-edges', [{{
            x: CLUSTER_BY_EDGES.map(d => d[0]),
            y: CLUSTER_BY_EDGES.map(d => d[1]),
            mode: 'lines',
            name: 'Clusters',
            line: {{color: '#ff7f0e', width: 2}},
            hovertemplate: 'Edges: %{{x}}<br>Clusters: %{{y}}<extra></extra>'
        }}], {{
            ...chartLayout,
            title: 'Clusters vs Edge Count',
            xaxis: {{title: 'Edge Count'}},
            yaxis: {{title: 'Number of Clusters'}},
            hovermode: 'closest'
        }}, chartConfig);
        
{component_signal_plot_block}
        
    buildClusterCheckpoints();

        // Initial cluster size histogram
        updateHistogramFromSliderPosition(0);
    }}

    function sliderStopAt(position) {{
        if (SLIDER_STOPS.length === 0) {{
            return {{edge_index: -1, threshold_label: '∞', slider_position: 0}};
        }}
        const maxPosition = SLIDER_STOPS[SLIDER_STOPS.length - 1].slider_position;
        const boundedPosition = Math.max(0, Math.min(position, maxPosition));
        let nearestStop = SLIDER_STOPS[0];
        let nearestDistance = Math.abs(boundedPosition - nearestStop.slider_position);
        for (let stopIdx = 1; stopIdx < SLIDER_STOPS.length; stopIdx++) {{
            const stop = SLIDER_STOPS[stopIdx];
            const distance = Math.abs(boundedPosition - stop.slider_position);
            if (distance < nearestDistance) {{
                nearestStop = stop;
                nearestDistance = distance;
            }}
        }}
        return nearestStop;
    }}

    function updateHistogramFromSliderPosition(position, snapSlider = true) {{
        const stop = sliderStopAt(position);
        if (snapSlider) {{
            document.getElementById('threshold-slider').value = stop.slider_position.toString();
        }}
        updateHistogram(stop.edge_index, stop);
    }}
    
    // Update histogram based on slider
    function updateHistogram(edgeIndex, sliderStop) {{
        const chartConfig = {{displayModeBar: false, responsive: true}};
        const currentK = {mst_knn_current_k_expr};
        
        let sizes;
        if (edgeIndex < 0) {{
            sizes = Array(VIZ_DATA.n_nodes).fill(1);
        }} else {{
            sizes = computeClusters(edgeIndex);
        }}
        
        updateClusterBubbleChart(sizes);
        
        // Update stats
        document.getElementById('num-clusters').textContent = formatNumber(sizes.length);
        document.getElementById('num-singletons').textContent = formatNumber(sizes.filter(size => size === 1).length);
        document.getElementById('largest-cluster').textContent = formatNumber(sizes[0] || 0);
        
        // Get the actual edge count from EDGES_BY_THRESH
        let edgeCount = 0;
        if (edgeIndex >= 0 && edgeIndex < EDGES_BY_THRESH.length) {{
            edgeCount = EDGES_BY_THRESH[edgeIndex][0];
        }}
        document.getElementById('num-edges').textContent = formatNumber(edgeCount);
{mst_knn_update_block}
        
        const avgSize = sizes.reduce((a, b) => a + b, 0) / sizes.length;
        document.getElementById('avg-cluster-size').textContent = avgSize.toFixed(1);
        
        let threshVal;
        let threshPos;
        if (sliderStop) {{
            threshVal = sliderStop.threshold_label;
            threshPos = (SLIDER_STOPS.indexOf(sliderStop) + 1).toString();
        }} else if (edgeIndex < 0) {{
            threshVal = '∞';
            threshPos = '0';
        }} else if (edgeIndex >= VIZ_DATA.mst_edges.length) {{
            threshVal = VIZ_DATA.mst_edges[VIZ_DATA.mst_edges.length - 1][2].toFixed(2);
            threshPos = VIZ_DATA.mst_edges.length.toString();
        }} else {{
            threshVal = VIZ_DATA.mst_edges[edgeIndex][2].toFixed(2);
            threshPos = (edgeIndex + 1).toString();
        }}
        document.getElementById('threshold-number').textContent = threshVal;
        
    }}
    
    // Slider event listener with debouncing for better performance
    let updateTimeout;
    document.getElementById('threshold-slider').addEventListener('input', (e) => {{
        clearTimeout(updateTimeout);
        updateTimeout = setTimeout(() => {{
            let position = parseInt(e.target.value);
            updateHistogramFromSliderPosition(position, false);
        }}, 10); // Small delay to batch rapid slider movements
    }});
    document.getElementById('threshold-slider').addEventListener('change', (e) => {{
        let position = parseInt(e.target.value);
        updateHistogramFromSliderPosition(position, true);
    }});
{mst_knn_listener_block}

    decodeEmbeddedReportData(EMBEDDED_REPORT_DATA_GZIP_BASE64)
        .then(reportData => {{
            installReportData(reportData);
            initCharts();
        }})
        .catch(error => {{
            console.error(error);
            document.querySelector('.dashboard').insertAdjacentHTML(
                'afterbegin',
                '<div class="chart chart-wide">Unable to load embedded report data in this browser.</div>'
            );
        }});
</script>""", file=self.out_handle)

    def write_footer(self):
        print("</body></html>", file=self.out_handle)

def matrix_report(matrix:DataMatrix, out_text_handle, out_html_handle, include_mst_knn: bool = False,
                  merge_impact_metric: str = MERGE_IMPACT_MIN_CHILD,
                  profile_stages: bool = False, stage_timings=None, progress_callback=None,
                  max_merge_events: int = DEFAULT_MAX_MERGE_EVENTS):
    """
        Write a report on the matrix to the given handles.
    """
    if merge_impact_metric not in MERGE_IMPACT_CHOICES:
        raise ValueError(f"merge_impact_metric must be one of {sorted(MERGE_IMPACT_CHOICES)}")
    if max_merge_events < 0:
        raise ValueError("max_merge_events must be >= 0")

    longform_outputs = list()
    if out_text_handle is not None:
        longform_outputs.append(SummaryTextWriter(out_text_handle))
    if out_html_handle is not None:
        longform_outputs.append(SummaryHTMLWriter(out_html_handle))

    if stage_timings is None:
        stage_timings = []
    
    if progress_callback is not None:
        progress_callback("building edge table")
    stage_start = perf_counter()
    edge_table = matrix.sorted_undirected_edges(skip_zeros=True, agg=max)
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
    if include_mst_knn:
        if progress_callback is not None:
            progress_callback("building MST_KNN neighbor rankings")
        stage_start = perf_counter()
        mst_knn_config = _get_mst_knn_report_config(tree)
        neighbor_rankings = build_symmetric_neighbor_rankings(edge_table, max_k=mst_knn_config["max_k"])
        _record_stage_timing(
            stage_timings,
            "build_mst_knn_neighbor_rankings",
            stage_start,
            max_k=mst_knn_config["max_k"],
            kept_directed_edges=len(neighbor_rankings.target),
            estimated_counts_bytes=_estimate_mst_knn_counts_bytes(len(tree.mst_edges), mst_knn_config["max_k"]),
        )

    if include_mst_knn:
        if progress_callback is not None:
            progress_callback("building MST_KNN threshold counts")
        stage_start = perf_counter()
        mst_knn_counts = mst_knn_edge_counts_by_threshold(matrix, tree, mst_knn_config["max_k"], neighbor_rankings=neighbor_rankings)
        _record_stage_timing(
            stage_timings,
            "build_mst_knn_counts",
            stage_start,
            thresholds=len(tree.mst_edges),
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

def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument("-i", "--input", default=None, required=True, type=str, help="name of input matrix file.")
    
    parser.add_argument('-o', '--output', default=None, required=False,
                        help="Text file to write output to.")
    
    parser.add_argument('--html', default=None, required=False,
                        help="html file to write output to.")
    parser.add_argument('--include_mst_knn', action='store_true', default=False,
                        help="Include projected MST_KNN edge counts and HTML controls. Disabled by default to reduce memory use.")
    parser.add_argument('--merge_impact_metric', choices=list(MERGE_IMPACT_CHOICES), default=MERGE_IMPACT_MIN_CHILD,
                        help="Metric used for split-event plots and tables: 'min_child' emphasizes the smaller cluster breaking away, while 'product' emphasizes balanced large splits.")
    parser.add_argument('--max_merge_events', type=int, default=DEFAULT_MAX_MERGE_EVENTS,
                        help="Maximum number of strongest merge events to embed in the interactive HTML threshold slider and split plot. Use 0 to include all merge events.")
    parser.add_argument('--progress', action='store_true', default=False,
                        help="Print live progress updates to stderr during long-running matrix_report stages.")
    parser.add_argument('--profile_stages', action='store_true', default=False,
                        help="Print per-stage runtime timings to stderr for profiling large matrix_report runs.")
    
    parser.add_argument('--config', action=ActionConfigFile)

    params = parser.parse_args(argv)

    if params.output is None and params.html is None:
        parser.error("No output file specified. Use --output, and/or --html to specify at least one output file.")

    with ExitStack() as output_stack:
        out = output_stack.enter_context(open(params.output, "w")) if params.output is not None else None
        out_html_handle = output_stack.enter_context(open(params.html, "w")) if params.html is not None else None

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
        if params.output is not None or params.html is not None:
            stage_start = perf_counter()
            matrix_report(
                matrix,
                out,
                out_html_handle,
                include_mst_knn=params.include_mst_knn,
                merge_impact_metric=params.merge_impact_metric,
                profile_stages=params.profile_stages,
                stage_timings=stage_timings,
                progress_callback=progress_callback,
                max_merge_events=params.max_merge_events,
            )
            if params.profile_stages:
                _record_stage_timing(stage_timings, "matrix_report_total", stage_start)
                _emit_stage_timings(stage_timings)

def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':

    main(sys.argv[1:])
  
