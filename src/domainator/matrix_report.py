""" Reads a DataMatrix and generates a report suitable for helping with selecting edge score thresholds for similarity networks
"""
from jsonargparse import ArgumentParser, ActionConfigFile
import sys
import math
from time import perf_counter

from domainator.data_matrix import DataMatrix, MaxTree, build_symmetric_neighbor_rankings, mst_knn_edge_counts_by_threshold, average_closeness_by_threshold
from domainator import __version__, RawAndDefaultsFormatter
from bashplotlib.histogram import plot_hist
from bashplotlib.scatterplot import plot_scatter
import io
from contextlib import redirect_stdout
from domainator.summary_report import histplot_base64
import statistics
import numpy as np
import json

MST_KNN_MIN_K = 2
MST_KNN_MAX_K = 25


def _get_mst_knn_report_config(tree):
    max_neighbors = max(tree.n_nodes - 1, MST_KNN_MIN_K)
    max_k = min(MST_KNN_MAX_K, max_neighbors)
    default_k = min(5, max_k)
    return {
        "min_k": MST_KNN_MIN_K,
        "max_k": max_k,
        "default_k": default_k,
    }


def _format_threshold_value(threshold):
    if math.isinf(float(threshold)):
        return "∞"
    return f"{float(threshold):.2f}"


def _summarize_closeness_curve(closeness_curve, max_items=5):
    if closeness_curve is None or len(closeness_curve.get("points", [])) == 0:
        return [], []

    points = np.asarray(closeness_curve["points"], dtype=float)
    values = points[:, 1]

    if len(values) == 1:
        local_maxima = [0]
    else:
        local_maxima = []
        for idx, value in enumerate(values):
            left = values[idx - 1] if idx > 0 else -np.inf
            right = values[idx + 1] if idx + 1 < len(values) else -np.inf
            if value >= left and value >= right and (value > left or value > right):
                local_maxima.append(idx)

    if len(local_maxima) == 0:
        local_maxima = [int(np.argmax(values))]

    ranked_optima = sorted(
        local_maxima,
        key=lambda idx: (-points[idx, 1], -points[idx, 0]),
    )[:max_items]

    change_rows = []
    for idx in range(len(values) - 1):
        delta = float(values[idx + 1] - values[idx])
        change_rows.append((abs(delta), idx, delta))

    ranked_changes = sorted(change_rows, reverse=True)[:max_items]
    return ranked_optima, ranked_changes


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
        
    def write_plots(self, tree, edge_scores, mst_knn_config, mst_knn_counts, closeness_curve):
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

        if closeness_curve is not None and len(closeness_curve.get("points", [])) > 0:
            optima, changes = _summarize_closeness_curve(closeness_curve)
            points = np.asarray(closeness_curve["points"], dtype=float)
            print(
                f"\nAverage closeness summary ({closeness_curve['mode']}; evaluated {closeness_curve['evaluated_thresholds']}/{closeness_curve['total_thresholds']} thresholds):",
                file=self.out_handle,
            )
            print("  Local optima:", file=self.out_handle)
            for idx in optima:
                threshold, closeness, edge_count, active_nodes, component_count = points[idx]
                print(
                    f"    Threshold: {_format_threshold_value(threshold)} → Avg closeness: {closeness:.4f}, "
                    f"Edges: {int(edge_count)}, Active nodes: {int(active_nodes)}, Components: {int(component_count)}",
                    file=self.out_handle,
                )
            print("  Steepest changes:", file=self.out_handle)
            for _, idx, delta in changes:
                left = points[idx]
                right = points[idx + 1]
                direction = "rise" if delta > 0 else "drop"
                print(
                    f"    {_format_threshold_value(left[0])} → {_format_threshold_value(right[0])}: {direction} {delta:+.4f}",
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
        grid-template-columns: 1fr 1fr;
        gap: 20px;
        margin-top: 20px;
    }}
    .chart {{
        border: 1px solid #e0e0e0;
        padding: 15px;
        background: white;
        border-radius: 4px;
        min-height: 400px;
        position: relative;
    }}
    .chart-with-controls {{
        border: 1px solid #e0e0e0;
        padding: 15px;
        background: white;
        border-radius: 4px;
        min-height: 400px;
        position: relative;
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
    
    
    def write_plots(self, tree, edge_scores, mst_knn_config, mst_knn_counts, closeness_curve):
        # Export data for interactive visualizations
        viz_data = tree.export_for_interactive_viz()
        cluster_by_thresh = tree.cluster_count_by_threshold
        cluster_by_edges = tree.cluster_count_by_edge_count
        edges_by_thresh = tree.edges_by_threshold
        include_mst_knn = mst_knn_counts is not None
        include_closeness = closeness_curve is not None and len(closeness_curve.get("points", [])) > 0
        mst_knn_controls = ""
        mst_knn_stat_box = ""
        mst_knn_current_k_expr = "null"
        mst_knn_update_block = ""
        mst_knn_listener_block = ""
        closeness_chart = ""
        closeness_plot_block = ""
        if include_closeness:
            closeness_chart = """
    <div class=\"chart\">
        <div id=\"closeness-by-threshold\"></div>
    </div>"""
            closeness_plot_block = """

        const closenessPoints = CLOSENESS_CURVE.points;
        Plotly.newPlot('closeness-by-threshold', [{
            x: closenessPoints.map(d => d[0]),
            y: closenessPoints.map(d => d[1]),
            mode: 'lines+markers',
            name: 'Avg closeness',
            line: {color: '#9467bd', width: 2},
            marker: {size: 6},
            customdata: closenessPoints.map(d => [d[2], d[3], d[4]]),
            hovertemplate: 'Threshold: %{x:.2f}<br>Avg closeness: %{y:.4f}<br>Edges: %{customdata[0]}<br>Active nodes: %{customdata[1]}<br>Components: %{customdata[2]}<extra></extra>'
        }], {
            ...chartLayout,
            title: `Average Closeness vs Threshold (${CLOSENESS_CURVE.mode}, ${CLOSENESS_CURVE.evaluated_thresholds}/${CLOSENESS_CURVE.total_thresholds})`,
            xaxis: {title: 'Threshold', type: 'log', autorange: true},
            yaxis: {title: 'Average Closeness', range: [0, 1.05]},
            hovermode: 'closest'
        }, chartConfig);"""
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
            let index = parseInt(document.getElementById('threshold-slider').value);
            updateHistogram(index);
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
    <div class="chart">
        <div id="edges-by-threshold"></div>
    </div>
    {closeness_chart}
    <div class="chart-with-controls">
        <div class="slider-container">
            <label for="threshold-slider">
                Threshold: <span id="threshold-number">∞</span>
            </label>
            <input type="range" id="threshold-slider" min="-1" max="{len(tree.mst) - 1}" value="-1" step="1">
            {mst_knn_controls}
        </div>
        <div class="stats">
            <div class="stat-box">
                <strong>Number of Clusters</strong>
                <span id="num-clusters">{tree.n_nodes}</span>
            </div>
            <div class="stat-box">
                <strong>Largest Cluster</strong>
                <span id="largest-cluster">1</span>
            </div>
            <div class="stat-box">
                <strong>Number of Edges</strong>
                <span id="num-edges">0</span>
            </div>
            {mst_knn_stat_box}
            <div class="stat-box">
                <strong>Avg Cluster Size</strong>
                <span id="avg-cluster-size">1.0</span>
            </div>
        </div>
        <div id="cluster-size-hist"></div>
    </div>
</div>
</div>

<script>
    // Embedded data
    const VIZ_DATA = {json.dumps(viz_data)};
    const CLUSTER_BY_THRESH = {json.dumps(cluster_by_thresh.tolist())};
    const CLUSTER_BY_EDGES = {json.dumps(cluster_by_edges.tolist())};
    const EDGES_BY_THRESH = {json.dumps(edges_by_thresh.tolist())};
    const HAS_MST_KNN = {json.dumps(include_mst_knn)};
    const MST_KNN_COUNTS = {json.dumps(mst_knn_counts.tolist() if mst_knn_counts is not None else [])};
    const MST_KNN_MIN_K = {mst_knn_config['min_k'] if mst_knn_config is not None else 0};
    const HAS_CLOSENESS = {json.dumps(include_closeness)};
    const CLOSENESS_CURVE = {json.dumps(closeness_curve if closeness_curve is not None else {"mode": "exact", "total_thresholds": 0, "evaluated_thresholds": 0, "points": []})};
    
    // Union-Find implementation for fast cluster computation
    class UnionFind {{
        constructor(n) {{
            this.parent = Array.from({{length: n}}, (_, i) => i);
            this.size = Array(n).fill(1);
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
    }}
    
    // Compute clusters at specific MST edge index
    function computeClusters(edgeIndex) {{
        let uf = new UnionFind(VIZ_DATA.n_nodes);
        for (let i = 0; i <= edgeIndex && i < VIZ_DATA.mst_edges.length; i++) {{
            let [n1, n2, val] = VIZ_DATA.mst_edges[i];
            uf.union(n1, n2);
        }}
        return uf.getClusterSizes();
    }}
    
    // Create histogram from cluster sizes
    function createHistogram(sizes) {{
        let hist = {{}};
        sizes.forEach(s => hist[s] = (hist[s] || 0) + 1);
        return hist;
    }}
    
    // Format large numbers
    function formatNumber(num) {{
        return num.toLocaleString();
    }}
    
    // Initialize static charts
    function initCharts() {{
        const chartConfig = {{displayModeBar: false, responsive: true}};
        const chartLayout = {{
            autosize: true,
            height: 350,
            margin: {{t: 40, b: 50, l: 60, r: 20}}
        }};
        
        // Cluster count vs threshold
        Plotly.newPlot('cluster-by-threshold', [{{
            x: CLUSTER_BY_THRESH.slice(1).map(d => d[0]),
            y: CLUSTER_BY_THRESH.slice(1).map(d => d[1]),
            mode: 'lines',
            name: 'Clusters',
            line: {{color: '#1f77b4', width: 2}},
            hovertemplate: 'Threshold: %{{x:.2f}}<br>Clusters: %{{y}}<extra></extra>'
        }}], {{
            ...chartLayout,
            title: 'Clusters vs Threshold',
            xaxis: {{title: 'Threshold', type: 'log', autorange: true}},
            yaxis: {{title: 'Number of Clusters'}},
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
        
        // Edges vs threshold
        Plotly.newPlot('edges-by-threshold', [{{
            x: EDGES_BY_THRESH.map(d => d[1]),
            y: EDGES_BY_THRESH.map(d => d[0]),
            mode: 'lines',
            name: 'Edges',
            line: {{color: '#2ca02c', width: 2}},
            hovertemplate: 'Threshold: %{{x:.2f}}<br>Edges: %{{y}}<extra></extra>'
        }}], {{
            ...chartLayout,
            title: 'Edge Count vs Threshold',
            xaxis: {{title: 'Threshold', type: 'log', autorange: true}},
            yaxis: {{title: 'Cumulative Edges'}},
            hovermode: 'closest'
        }}, chartConfig);

{closeness_plot_block}
        
        // Initial cluster size histogram
        updateHistogram(-1);
    }}
    
    // Update histogram based on slider
    function updateHistogram(edgeIndex) {{
        const chartConfig = {{displayModeBar: false, responsive: true}};
        const currentK = {mst_knn_current_k_expr};
        
        let sizes;
        if (edgeIndex < 0) {{
            sizes = Array(VIZ_DATA.n_nodes).fill(1);
        }} else {{
            sizes = computeClusters(edgeIndex);
        }}
        
        let hist = createHistogram(sizes);
        let x = Object.keys(hist).map(Number).sort((a, b) => a - b);
        let y = x.map(size => hist[size]);
        
        Plotly.newPlot('cluster-size-hist', [{{
            x: x,
            y: y,
            type: 'bar',
            marker: {{color: '#d62728'}},
            hovertemplate: 'Size: %{{x}}<br>Count: %{{y}}<extra></extra>'
        }}], {{
            title: 'Cluster Size Distribution',
            xaxis: {{
                title: 'Cluster Size',
                type: 'log'
            }},
            yaxis: {{title: 'Number of Clusters',
                    type: 'log'
            }},
            autosize: true,
            height: 240,
            margin: {{t: 40, b: 50, l: 60, r: 20}},
            hovermode: 'closest'
        }}, chartConfig);
        
        // Update stats
        document.getElementById('num-clusters').textContent = formatNumber(sizes.length);
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
        if (edgeIndex < 0) {{
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
            let index = parseInt(e.target.value);
            updateHistogram(index);
        }}, 10); // Small delay to batch rapid slider movements
    }});
{mst_knn_listener_block}
    
    // Initialize on load
    initCharts();
</script>""", file=self.out_handle)

    def write_footer(self):
        print("</body></html>", file=self.out_handle)

def matrix_report(matrix:DataMatrix, out_text_handle, out_html_handle, include_mst_knn: bool = False,
                  include_closeness: bool = False, closeness_mode: str = "auto", closeness_max_points: int = 128,
                  profile_stages: bool = False, stage_timings=None, progress_callback=None):
    """
        Write a report on the matrix to the given handles.
    """
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
    mst_knn_config = None
    mst_knn_counts = None
    closeness_curve = None
    if include_mst_knn:
        if progress_callback is not None:
            progress_callback("building MST_KNN neighbor rankings")
        stage_start = perf_counter()
        mst_knn_config = _get_mst_knn_report_config(tree)
        neighbor_rankings = build_symmetric_neighbor_rankings(edge_table, max_k=mst_knn_config["max_k"])
        _record_stage_timing(stage_timings, "build_mst_knn_neighbor_rankings", stage_start, max_k=mst_knn_config["max_k"], kept_edges=len(neighbor_rankings.target))

        if progress_callback is not None:
            progress_callback("building MST_KNN threshold counts")
        stage_start = perf_counter()
        mst_knn_counts = mst_knn_edge_counts_by_threshold(matrix, tree, mst_knn_config["max_k"], neighbor_rankings=neighbor_rankings)
        _record_stage_timing(stage_timings, "build_mst_knn_counts", stage_start, thresholds=len(tree.mst_edges), max_k=mst_knn_config["max_k"])
    if include_closeness:
        if progress_callback is not None:
            progress_callback("building closeness curve")
        stage_start = perf_counter()
        closeness_curve = average_closeness_by_threshold(
            edge_table,
            tree=tree,
            mode=closeness_mode,
            max_points=closeness_max_points,
            progress_callback=progress_callback,
        )
        _record_stage_timing(stage_timings, "build_closeness_curve", stage_start, mode=closeness_curve["mode"], evaluated=closeness_curve["evaluated_thresholds"], total=closeness_curve["total_thresholds"])

    if progress_callback is not None:
        progress_callback("rendering outputs")
    stage_start = perf_counter()
    for output in longform_outputs:
        output.write_header(tree, edge_table.score)
        output.write_plots(tree, edge_table.score, mst_knn_config, mst_knn_counts, closeness_curve)
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
    parser.add_argument('--include_closeness', action='store_true', default=False,
                        help="Include paper-style average closeness by threshold in the text and HTML reports.")
    parser.add_argument('--closeness_mode', choices=['auto', 'exact', 'sampled'], default='auto',
                        help="How to evaluate closeness thresholds when --include_closeness is enabled.")
    parser.add_argument('--closeness_max_points', type=int, default=128,
                        help="Maximum number of thresholds to evaluate when sampled closeness mode is used.")
    parser.add_argument('--progress', action='store_true', default=False,
                        help="Print live progress updates to stderr during long-running matrix_report stages.")
    parser.add_argument('--profile_stages', action='store_true', default=False,
                        help="Print per-stage runtime timings to stderr for profiling large matrix_report runs.")
    
    parser.add_argument('--config', action=ActionConfigFile)

    params = parser.parse_args(argv)

    if params.output is None and params.html is None:
        parser.error("No output file specified. Use --output, and/or --html to specify at least one output file.")

    out = None
    if params.output is not None:
        out = open(params.output, "w")
    
    out_html_handle = None
    if params.html is not None:
        out_html_handle = open(params.html, "w")

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
            include_closeness=params.include_closeness,
            closeness_mode=params.closeness_mode,
            closeness_max_points=params.closeness_max_points,
            profile_stages=params.profile_stages,
            stage_timings=stage_timings,
            progress_callback=progress_callback,
        )
        if params.profile_stages:
            _record_stage_timing(stage_timings, "matrix_report_total", stage_start)
            _emit_stage_timings(stage_timings)

    if params.output is not None:
        out.close()
    
    if params.html is not None:
        out_html_handle.close()

def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':

    main(sys.argv[1:])
  
