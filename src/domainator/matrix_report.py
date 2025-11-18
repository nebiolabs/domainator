""" reads a DataMatrix and generates a report suitable for helping with selecting edge score thresholds for generating similatity networks.
"""
from jsonargparse import ArgumentParser, ActionConfigFile
import sys

from domainator.data_matrix import DataMatrix, MaxTree
from domainator import __version__, RawAndDefaultsFormatter
from bashplotlib.histogram import plot_hist
from bashplotlib.scatterplot import plot_scatter
import io
from contextlib import redirect_stdout
from domainator.summary_report import histplot_base64
import statistics
import numpy as np
import json

class SummaryTextWriter():
    def __init__(self, out_handle): 
        self.out_handle = out_handle
    
    def write_header(self, tree):
        # Get triangular values from the tree
        non_zero_values = tree.edges[:, 2]
        n_nodes = tree.n_nodes
        
        print(f"""Matrix Report
Nodes: {n_nodes}
Non-zero edges: {len(non_zero_values)}
Mean: {statistics.mean(non_zero_values) if len(non_zero_values) > 0 else 0:.1f}
Median: {statistics.median(non_zero_values) if len(non_zero_values) > 0 else 0:.1f}
Min: {min(non_zero_values) if len(non_zero_values) > 0 else 0:.1f}
Max: {max(non_zero_values) if len(non_zero_values) > 0 else 0:.1f}
        """, file=self.out_handle)
        
    def write_plots(self, tree):
        values = tree.edges[:, 2]
        non_zero_values = values[values > 0]
        
        # Histogram of all triangular values
        if len(non_zero_values) > 0:
            hist_str = io.StringIO()
            with redirect_stdout(hist_str):
                plot_hist(list(non_zero_values), bincount=50, title="Edge counts by score", height=20.0, xlab=True, regular=True, showSummary=False)
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
    
    def write_footer(self):
        pass

class SummaryHTMLWriter():
    def __init__(self, out_handle): 
        self.out_handle = out_handle
    
    def write_header(self, tree):
        # Get triangular values from the tree
        non_zero_values = tree.edges[:, 2]
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
    
    
    def write_plots(self, tree):
        # Export data for interactive visualizations
        viz_data = tree.export_for_interactive_viz()
        cluster_by_thresh = tree.cluster_count_by_threshold
        cluster_by_edges = tree.cluster_count_by_edge_count
        edges_by_thresh = tree.edges_by_threshold
        
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
    <div class="chart-with-controls">
        <div class="slider-container">
            <label for="threshold-slider">
                Threshold: <span id="threshold-number">∞</span>
            </label>
            <input type="range" id="threshold-slider" min="-1" max="{len(tree.mst) - 1}" value="-1" step="1">
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
        
        // Initial cluster size histogram
        updateHistogram(-1);
    }}
    
    // Update histogram based on slider
    function updateHistogram(edgeIndex) {{
        const chartConfig = {{displayModeBar: false, responsive: true}};
        
        let sizes;
        if (edgeIndex < 0) {{
            sizes = Array(VIZ_DATA.n_nodes).fill(1);
        }} else {{
            sizes = computeClusters(edgeIndex);
        }}
        
        let hist = createHistogram(sizes);
        let x = Object.keys(hist).map(Number).sort((a, b) => a - b);
        let y = x.map(size => hist[size]);
        
        // Determine if log scale is appropriate
        const maxSize = Math.max(...x);
        const useLogScale = maxSize > 100;
        
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
                type: useLogScale ? 'log' : 'linear'
            }},
            yaxis: {{title: 'Number of Clusters',
                    type: useLogScale ? 'log' : 'linear'
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
    
    // Initialize on load
    initCharts();
</script>""", file=self.out_handle)

    def write_footer(self):
        print("</body></html>", file=self.out_handle)

def matrix_report(matrix:DataMatrix, out_text_handle, out_html_handle):
    """
        Write a report on the matrix to the given handles.
    """
    longform_outputs = list()
    if out_text_handle is not None:
        longform_outputs.append(SummaryTextWriter(out_text_handle))
    if out_html_handle is not None:
        longform_outputs.append(SummaryHTMLWriter(out_html_handle))
    
    tree = MaxTree(matrix)

    for output in longform_outputs:
        output.write_header(tree)
        output.write_plots(tree)
        output.write_footer()

def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument("-i", "--input", default=None, required=True, type=str, help="name of input matrix file.")
    
    parser.add_argument('-o', '--output', default=None, required=False,
                        help="Text file to write output to.")
    
    parser.add_argument('--html', default=None, required=False,
                        help="html file to write output to.")
    
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
 
    matrix = DataMatrix.from_file(params.input)

    ### Run
    if params.output is not None or params.html is not None:
        matrix_report(
            matrix,
            out,
            out_html_handle,
        )

    if params.output is not None:
        out.close()
    
    if params.html is not None:
        out_html_handle.close()

def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':

    main(sys.argv[1:])
  
