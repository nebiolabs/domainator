"""Static HTML shell for the standalone SSN viewer."""

from pathlib import Path

from domainator.output_guardrails import make_temporary_output_path


def ssn_viewer_html(title: str = "Domainator SSN Viewer") -> str:
    escaped_title = title.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width, initial-scale=1.0" />
<title>{escaped_title}</title>
<style>
    :root {{
        color-scheme: light;
        --bg: #f3efe7;
        --panel: #fffaf2;
        --panel-strong: #fffdf9;
        --ink: #1e2a2f;
        --muted: #5c6a70;
        --accent: #c8553d;
        --accent-soft: #f3cdb9;
        --line: #d8cec2;
        --shadow: 0 16px 40px rgba(61, 43, 31, 0.10);
    }}
    * {{ box-sizing: border-box; }}
    body {{
        margin: 0;
        font-family: Georgia, "Iowan Old Style", "Palatino Linotype", serif;
        color: var(--ink);
        background:
            radial-gradient(circle at top left, rgba(200, 85, 61, 0.10), transparent 30%),
            linear-gradient(180deg, #f7f2eb 0%, var(--bg) 100%);
    }}
    .shell {{
        max-width: 1680px;
        margin: 0 auto;
        padding: 24px;
    }}
    .hero {{
        display: flex;
        flex-wrap: wrap;
        justify-content: space-between;
        gap: 18px;
        margin-bottom: 20px;
        padding: 22px 24px;
        background: rgba(255, 250, 242, 0.85);
        border: 1px solid rgba(216, 206, 194, 0.8);
        box-shadow: var(--shadow);
        border-radius: 18px;
        backdrop-filter: blur(8px);
    }}
    .hero h1 {{
        margin: 0;
        font-size: clamp(2rem, 4vw, 3.6rem);
        font-weight: 600;
        letter-spacing: -0.03em;
    }}
    .hero p {{
        margin: 8px 0 0 0;
        max-width: 60ch;
        color: var(--muted);
        line-height: 1.45;
    }}
    .loader {{
        display: flex;
        flex-direction: column;
        gap: 10px;
        align-items: flex-start;
        justify-content: center;
        min-width: 280px;
    }}
    .loader input[type=file] {{
        max-width: 100%;
        font: inherit;
    }}
    .status {{
        font-size: 0.96rem;
        color: var(--muted);
    }}
    .grid {{
        display: grid;
        grid-template-columns: minmax(320px, 0.95fr) minmax(520px, 1.55fr) minmax(320px, 0.9fr);
        gap: 18px;
        align-items: start;
    }}
    .panel {{
        background: var(--panel);
        border: 1px solid var(--line);
        border-radius: 18px;
        box-shadow: var(--shadow);
        overflow: hidden;
        min-width: 0;
    }}
    .panel h2 {{
        margin: 0;
        padding: 18px 20px 0 20px;
        font-size: 1.15rem;
        letter-spacing: 0.01em;
    }}
    .panel-body {{
        padding: 16px 20px 20px 20px;
    }}
    .wide {{ grid-column: 1 / span 2; }}
    .sidebar {{ grid-column: 3; }}
    .full-width {{ grid-column: 1 / -1; }}
    .controls {{
        display: flex;
        flex-direction: column;
        gap: 12px;
        margin-bottom: 16px;
    }}
    .control label {{
        display: block;
        margin-bottom: 5px;
        color: var(--muted);
        font-size: 0.92rem;
    }}
    .control input,
    .control select,
    .toolbar button {{
        width: 100%;
        border: 1px solid var(--line);
        background: var(--panel-strong);
        border-radius: 10px;
        padding: 9px 10px;
        font: inherit;
        color: var(--ink);
    }}
    .checkbox {{
        display: flex;
        align-items: center;
        gap: 8px;
    }}
    .checkbox input {{ width: auto; }}
    .toolbar {{
        display: flex;
        gap: 10px;
        flex-wrap: wrap;
        margin-bottom: 12px;
    }}
    .toolbar button {{
        width: auto;
        min-width: 150px;
        cursor: pointer;
        transition: transform 140ms ease, background 140ms ease;
    }}
    .toolbar button:hover:not(:disabled) {{
        transform: translateY(-1px);
        background: #fff;
    }}
    .toolbar button:disabled {{
        opacity: 0.45;
        cursor: not-allowed;
    }}
    .split-topline {{
        display: block;
        margin-bottom: 12px;
    }}
    .stats {{
        display: grid;
        grid-template-columns: repeat(2, minmax(0, 1fr));
        gap: 10px;
        margin: 0;
    }}
    .slider-stack {{
        display: flex;
        flex-direction: column;
        gap: 8px;
        margin-top: 12px;
    }}
    .slider-track-wrap {{
        padding-left: clamp(24px, 5.1%, 56px);
        padding-right: clamp(14px, 2.8%, 30px);
    }}
    .slider-row {{
        display: grid;
        grid-template-columns: minmax(0, 1fr);
        gap: 8px;
        align-items: center;
    }}
    .slider-row input {{
        width: 100%;
        margin: 0;
        accent-color: var(--accent);
    }}
    .slider-extents {{
        display: flex;
        justify-content: space-between;
        gap: 12px;
        color: var(--muted);
        font-size: 0.88rem;
        letter-spacing: 0.01em;
    }}
    .stat {{
        padding: 12px;
        border-radius: 12px;
        background: rgba(255,255,255,0.65);
        border: 1px solid rgba(216, 206, 194, 0.7);
    }}
    .stat strong {{
        display: block;
        font-size: 0.78rem;
        color: var(--muted);
        text-transform: uppercase;
        letter-spacing: 0.06em;
        margin-bottom: 6px;
    }}
    .stat span {{
        font-size: 1.35rem;
    }}
    .canvas-wrap {{
        border-radius: 16px;
        overflow: hidden;
        background: linear-gradient(180deg, rgba(255,255,255,0.9), rgba(248,242,234,0.92));
        border: 1px solid rgba(216, 206, 194, 0.8);
    }}
    canvas {{
        display: block;
        width: 100%;
        height: auto;
        touch-action: none;
    }}
    .note {{
        margin-top: 10px;
        color: var(--muted);
        font-size: 0.92rem;
        line-height: 1.4;
    }}
    .table-wrap {{
        max-height: 680px;
        overflow: auto;
        border-radius: 12px;
        border: 1px solid rgba(216, 206, 194, 0.8);
        background: rgba(255,255,255,0.72);
    }}
    table {{
        border-collapse: collapse;
        width: 100%;
        font-size: 0.92rem;
    }}
    th, td {{
        padding: 8px 10px;
        border-bottom: 1px solid rgba(216, 206, 194, 0.65);
        text-align: left;
        vertical-align: top;
    }}
    thead th {{
        position: sticky;
        top: 0;
        background: #fffaf4;
        z-index: 1;
    }}
    .pill {{
        display: inline-flex;
        align-items: center;
        gap: 6px;
        padding: 6px 10px;
        border-radius: 999px;
        background: rgba(200, 85, 61, 0.09);
        color: #8a3b2c;
        font-size: 0.9rem;
    }}
    .split-legend {{
        display: flex;
        gap: 12px;
        flex-wrap: wrap;
        margin-top: 10px;
        color: var(--muted);
        font-size: 0.9rem;
        align-items: center;
        justify-content: space-between;
    }}
    .split-legend-items {{
        display: flex;
        gap: 12px;
        flex-wrap: wrap;
        align-items: center;
    }}
    .legend-threshold-readout {{
        display: inline-flex;
        align-items: center;
        gap: 8px;
    }}
    .legend-threshold-value {{
        display: inline-flex;
        align-items: center;
        padding: 4px 10px;
        border-radius: 999px;
        background: rgba(200, 85, 61, 0.09);
        color: #8a3b2c;
        font-weight: 600;
        letter-spacing: 0.01em;
    }}
    .threshold-jump {{
        display: inline-flex;
        align-items: center;
        gap: 8px;
        flex-wrap: wrap;
    }}
    .threshold-jump label {{
        color: var(--muted);
        font-size: 0.9rem;
    }}
    .threshold-jump input {{
        width: 132px;
        border: 1px solid var(--line);
        background: var(--panel-strong);
        border-radius: 10px;
        padding: 6px 9px;
        font: inherit;
        color: var(--ink);
    }}
    .legend-line::before,
    .legend-dot::before {{
        content: "";
        display: inline-block;
        vertical-align: middle;
        margin-right: 8px;
    }}
    .legend-line::before {{
        width: 18px;
        height: 2px;
        background: #c8553d;
    }}
    .legend-dot::before {{
        width: 10px;
        height: 10px;
        border-radius: 50%;
        background: #e29b4b;
    }}
    .network-summary {{
        display: flex;
        gap: 10px;
        flex-wrap: wrap;
        margin-bottom: 12px;
    }}
    .panel-body.tight {{
        padding-top: 12px;
    }}
    @media (max-width: 1380px) {{
        .grid {{ grid-template-columns: minmax(0, 1fr); }}
        .wide, .sidebar, .full-width {{ grid-column: auto; }}
        .split-topline {{ grid-template-columns: minmax(0, 1fr); }}
    }}
</style>
</head>
<body>
<div class="shell">
    <section class="hero">
        <div>
            <h1>{escaped_title}</h1>
            <p>Open a Domainator SSN viewer bundle produced by build_ssn_viewer.py. The viewer runs entirely in the browser, reconstructs threshold-dependent components from the MST hierarchy, and supports cluster selection, metadata export, and basic styling without any backend.</p>
        </div>
        <div class="loader">
            <input id="bundle-file" type="file" accept=".ssnv,.gz,.json,.ssnview" />
            <div id="bundle-status" class="status">No bundle loaded.</div>
            <div id="browser-warning" class="status"></div>
        </div>
    </section>
    <div class="grid">
        <section class="panel wide">
            <h2>Cluster Splits vs Threshold</h2>
            <div class="panel-body">
                <div class="split-topline">
                    <div class="stats">
                        <div class="stat"><strong>Clusters</strong><span id="stat-clusters">0</span></div>
                        <div class="stat"><strong>Split Links</strong><span id="stat-links">0</span></div>
                        <div class="stat"><strong>Shown Nodes</strong><span id="stat-shown-nodes">0</span></div>
                        <div class="stat"><strong>Hidden Nodes</strong><span id="stat-hidden-nodes">0</span></div>
                    </div>
                </div>
                <div class="canvas-wrap"><canvas id="split-chart" width="1100" height="320"></canvas></div>
                <div class="split-legend">
                    <div class="split-legend-items">
                        <span class="legend-line">Split impact trace</span>
                        <span class="legend-dot legend-threshold-readout">Current threshold <span id="threshold-label" class="legend-threshold-value">∞</span></span>
                    </div>
                    <div class="threshold-jump">
                        <label for="threshold-input">Jump to threshold</label>
                        <input id="threshold-input" type="number" step="any" inputmode="decimal" placeholder="Nearest split" disabled />
                    </div>
                </div>
                <div class="slider-stack">
                    <div class="slider-track-wrap">
                        <div class="slider-row">
                            <input id="threshold-slider" type="range" min="0" max="0" value="0" step="1" disabled />
                        </div>
                        <div class="slider-extents">
                            <span id="threshold-min-label">min</span>
                            <span id="threshold-max-label">∞</span>
                        </div>
                    </div>
                </div>
            </div>
        </section>
        <section class="panel wide">
            <h2>Hierarchy View</h2>
            <div class="panel-body">
                <div class="network-summary">
                    <div id="hidden-summary" class="pill">0 nodes hidden</div>
                    <div id="selection-summary" class="pill">0 nodes selected</div>
                </div>
                <div class="canvas-wrap"><canvas id="cluster-view" width="1100" height="700"></canvas></div>
                <div class="note">Wheel to zoom, drag the background to pan, and Shift-drag a box to select clusters. Click a cluster bubble to toggle selection for every node inside it.</div>
            </div>
        </section>
        <section class="panel sidebar">
            <h2>View Settings</h2>
            <div class="panel-body tight">
                <div class="controls">
                    <div class="control">
                        <label for="layout-algorithm">Layout algorithm</label>
                        <select id="layout-algorithm">
                            <option value="tree">Tree</option>
                            <option value="force">Force-directed</option>
                            <option value="organic">Organic</option>
                        </select>
                    </div>
                    <div class="control">
                        <label for="min-cluster-size">Minimum cluster size</label>
                        <input id="min-cluster-size" type="number" min="1" value="1" step="1" />
                    </div>
                    <div class="control">
                        <label for="color-by">Color by</label>
                        <select id="color-by"></select>
                    </div>
                    <div class="control">
                        <label for="label-by">Label by</label>
                        <select id="label-by"></select>
                    </div>
                    <div class="control checkbox">
                        <input id="show-labels" type="checkbox" />
                        <label for="show-labels">Show labels for singletons and small clusters</label>
                    </div>
                    <div class="control checkbox">
                        <input id="show-node-counts" type="checkbox" checked />
                        <label for="show-node-counts">Show node count labels</label>
                    </div>
                    <div class="control checkbox">
                        <input id="show-edge-scores" type="checkbox" />
                        <label for="show-edge-scores">Show edge score labels</label>
                    </div>
                    <div class="control checkbox">
                        <input id="leaf-pruning-only" type="checkbox" checked />
                        <label for="leaf-pruning-only">Minimum cluster size trims leaf clusters only</label>
                    </div>
                </div>
                <div class="toolbar">
                    <button id="reset-view" disabled>Reset view</button>
                    <button id="clear-selection" disabled>Clear selection</button>
                </div>
                <div class="note" id="selection-note">Load a bundle to begin exploring metadata.</div>
            </div>
        </section>
        <section class="panel full-width">
            <h2>Node Metadata</h2>
            <div class="panel-body">
                <div class="toolbar">
                    <button id="export-selected" disabled>Export table TSV</button>
                </div>
                <div class="table-wrap">
                    <table id="metadata-table">
                        <thead></thead>
                        <tbody></tbody>
                    </table>
                </div>
            </div>
        </section>
    </div>
</div>
<script>
    const state = {{
        bundle: null,
        metadataColumns: [],
        metadataRows: [],
        metadataByNodeIndex: [],
        activeClusters: [],
        visibleClusters: [],
        selectedNodeIndices: new Set(),
        visibleLayout: [],
        splitLinks: [],
        sliderModel: null,
        viewTransform: {{scale: 1, offsetX: 0, offsetY: 0, minScale: 0.02, maxScale: 10}},
        selectionBox: null,
        dragState: null,
        suppressClick: false,
        layoutCache: new Map(),
    }};

    const splitCanvas = document.getElementById('split-chart');
    const splitContext = splitCanvas.getContext('2d');
    const clusterCanvas = document.getElementById('cluster-view');
    const clusterContext = clusterCanvas.getContext('2d');

    function setStatus(message) {{
        document.getElementById('bundle-status').textContent = message;
    }}

    function browserSupportsBundleLoading() {{
        return 'DecompressionStream' in window;
    }}

    function warnIfUnsupported() {{
        if (!browserSupportsBundleLoading()) {{
            document.getElementById('browser-warning').textContent = 'This browser lacks DecompressionStream support; gzipped bundles cannot be loaded.';
        }}
    }}

    async function decodeBundleFile(file) {{
        if (browserSupportsBundleLoading()) {{
            try {{
                const stream = file.stream().pipeThrough(new DecompressionStream('gzip'));
                const text = await new Response(stream).text();
                return JSON.parse(text);
            }} catch (error) {{
                console.warn('gzip decode failed, trying plain JSON', error);
            }}
        }}
        const text = await file.text();
        return JSON.parse(text);
    }}

    function formatValue(value) {{
        if (value === null || value === undefined || value === '') {{
            return '—';
        }}
        if (typeof value === 'number') {{
            return Number.isInteger(value) ? value.toLocaleString() : value.toFixed(2);
        }}
        return String(value);
    }}

    function nodeId(nodeIndex) {{
        return state.bundle.graph.nodes[nodeIndex];
    }}

    function installBundle(bundle) {{
        if (bundle.format !== 'domainator_ssn_viewer_bundle') {{
            throw new Error('Unsupported bundle format: ' + bundle.format);
        }}
        state.bundle = bundle;
        state.metadataColumns = bundle.metadata.columns || [];
        state.metadataRows = bundle.metadata.rows || [];
        state.metadataByNodeIndex = state.metadataRows;
        state.selectedNodeIndices = new Set();
        state.sliderModel = buildSliderModel(bundle.graph.slider_stops || []);

        populateMetadataControls();
        const slider = document.getElementById('threshold-slider');
        slider.max = String(state.sliderModel.maxPosition);
        slider.value = String(state.sliderModel.initialPosition);
        slider.disabled = state.sliderModel.stops.length === 0;
        document.getElementById('threshold-min-label').textContent = state.sliderModel.minLabel;
        document.getElementById('threshold-max-label').textContent = state.sliderModel.maxLabel;
        document.getElementById('reset-view').disabled = false;
        document.getElementById('threshold-input').disabled = state.sliderModel.stops.length === 0;
        setStatus('Loaded ' + (bundle.name || 'bundle') + ' with ' + bundle.graph.nodes.length.toLocaleString() + ' nodes.');
        updateThresholdUI();
    }}

    function buildSliderModel(sourceStops) {{
        const stops = (sourceStops || []).map(stop => ({{...stop}}));
        if (stops.length === 0) {{
            return {{stops: [], maxPosition: 0, initialPosition: 0, minLabel: 'min', maxLabel: '∞'}};
        }}

        const infinityStop = stops.find(stop => stop.threshold_value === null) || {{edge_index: -1, threshold_label: '∞', threshold_value: null}};
        const finiteStops = stops
            .filter(stop => stop.threshold_value !== null)
            .sort((left, right) => left.threshold_value - right.threshold_value);

        if (finiteStops.length === 0) {{
            return {{
                stops: [{{...infinityStop, sliderPosition: 1000}}],
                maxPosition: 1000,
                initialPosition: 1000,
                minLabel: '∞',
                maxLabel: '∞',
            }};
        }}

        const minThreshold = finiteStops[0].threshold_value;
        const maxThreshold = finiteStops[finiteStops.length - 1].threshold_value;
        const finiteSpan = maxThreshold - minThreshold;
        const positionedStops = finiteStops.map((stop, stopIndex) => {{
            let sliderPosition = 0;
            if (finiteStops.length === 1) {{
                sliderPosition = 0;
            }} else if (finiteSpan <= 0) {{
                sliderPosition = Math.round((stopIndex / (finiteStops.length - 1)) * 920);
            }} else {{
                sliderPosition = Math.round(((stop.threshold_value - minThreshold) / finiteSpan) * 920);
            }}
            return {{...stop, sliderPosition}};
        }});
        positionedStops.push({{...infinityStop, sliderPosition: 1000}});

        return {{
            stops: positionedStops,
            maxPosition: 1000,
            initialPosition: 1000,
            minLabel: finiteStops[0].threshold_label,
            maxLabel: '∞',
        }};
    }}

    function currentSliderStop() {{
        if (!state.sliderModel || state.sliderModel.stops.length === 0) {{
            return null;
        }}
        const sliderPosition = Number(document.getElementById('threshold-slider').value);
        let nearestStop = state.sliderModel.stops[0];
        let nearestDistance = Math.abs(sliderPosition - nearestStop.sliderPosition);
        for (const stop of state.sliderModel.stops) {{
            const distance = Math.abs(sliderPosition - stop.sliderPosition);
            if (distance < nearestDistance) {{
                nearestStop = stop;
                nearestDistance = distance;
            }}
        }}
        return nearestStop;
    }}

    function snapSliderToStop(stop) {{
        if (!stop) {{
            return;
        }}
        document.getElementById('threshold-slider').value = String(stop.sliderPosition);
    }}

    function nearestStopForThreshold(targetValue) {{
        if (!state.sliderModel || state.sliderModel.stops.length === 0 || !Number.isFinite(targetValue)) {{
            return null;
        }}
        const finiteStops = state.sliderModel.stops.filter(stop => stop.threshold_value !== null);
        if (finiteStops.length === 0) {{
            return currentSliderStop();
        }}
        let nearestStop = finiteStops[0];
        let nearestDistance = Math.abs(finiteStops[0].threshold_value - targetValue);
        finiteStops.forEach(stop => {{
            const distance = Math.abs(stop.threshold_value - targetValue);
            if (distance < nearestDistance) {{
                nearestStop = stop;
                nearestDistance = distance;
            }}
        }});
        return nearestStop;
    }}

    function populateMetadataControls() {{
        const colorBy = document.getElementById('color-by');
        const labelBy = document.getElementById('label-by');
        colorBy.innerHTML = '';
        labelBy.innerHTML = '';

        const emptyOption = document.createElement('option');
        emptyOption.value = '';
        emptyOption.textContent = 'None';
        colorBy.appendChild(emptyOption.cloneNode(true));
        labelBy.appendChild(emptyOption.cloneNode(true));

        state.metadataColumns.forEach(column => {{
            const colorOption = document.createElement('option');
            colorOption.value = column.name;
            colorOption.textContent = column.name;
            colorBy.appendChild(colorOption);

            const labelOption = document.createElement('option');
            labelOption.value = column.name;
            labelOption.textContent = column.name;
            labelBy.appendChild(labelOption);
        }});

        colorBy.value = state.bundle.defaults.color_by || '';
        labelBy.value = state.bundle.defaults.label_by || '';
        document.getElementById('show-labels').checked = Boolean(state.bundle.defaults.label_by);
    }}

    function selectedThresholdValue() {{
        const stop = currentSliderStop();
        if (!stop) {{
            return Infinity;
        }}
        return stop.threshold_value === null ? Infinity : stop.threshold_value;
    }}

    function activeClustersAtThreshold(thresholdValue) {{
        const hierarchy = state.bundle.graph.hierarchy;
        const active = [];
        const stack = [...hierarchy.roots].reverse();
        while (stack.length > 0) {{
            const componentId = stack.pop();
            const component = hierarchy.nodes[componentId];
            if (component.kind === 'leaf') {{
                active.push(componentId);
                continue;
            }}
            if (component.threshold < thresholdValue) {{
                stack.push(component.right);
                stack.push(component.left);
                continue;
            }}
            active.push(componentId);
        }}
        return active;
    }}

    function componentMembers(componentId) {{
        const hierarchy = state.bundle.graph.hierarchy;
        const component = hierarchy.nodes[componentId];
        const start = component.leaf_start;
        const count = component.leaf_count;
        return hierarchy.leaf_order.slice(start, start + count);
    }}

    function leafPruningOnlyEnabled() {{
        return document.getElementById('leaf-pruning-only').checked;
    }}

    function currentLayoutAlgorithm() {{
        return document.getElementById('layout-algorithm').value || 'tree';
    }}

    function showNodeCountsEnabled() {{
        return document.getElementById('show-node-counts').checked;
    }}

    function showEdgeScoresEnabled() {{
        return document.getElementById('show-edge-scores').checked;
    }}

    function currentColorField() {{
        return document.getElementById('color-by').value;
    }}

    function currentLabelField() {{
        return document.getElementById('label-by').value;
    }}

    function metadataColumn(name) {{
        return state.metadataColumns.find(column => column.name === name) || null;
    }}

    function metadataValue(nodeIndex, columnName) {{
        if (!columnName) {{
            return null;
        }}
        const columnIndex = state.metadataColumns.findIndex(column => column.name === columnName);
        if (columnIndex < 0) {{
            return null;
        }}
        return state.metadataByNodeIndex[nodeIndex]?.[columnIndex] ?? null;
    }}

    function categoricalColor(value) {{
        if (value === null || value === undefined || value === '') {{
            return '#b3a89d';
        }}
        const text = String(value);
        let hash = 0;
        for (let i = 0; i < text.length; i++) {{
            hash = ((hash << 5) - hash) + text.charCodeAt(i);
            hash |= 0;
        }}
        const hue = Math.abs(hash) % 360;
        return 'hsl(' + hue + ' 58% 54%)';
    }}

    function numericColor(value, minValue, maxValue) {{
        if (value === null || value === undefined || Number.isNaN(value)) {{
            return '#b3a89d';
        }}
        const fraction = maxValue <= minValue ? 0.5 : (value - minValue) / (maxValue - minValue);
        const clamped = Math.max(0, Math.min(1, fraction));
        const hue = 200 - (160 * clamped);
        const light = 72 - (24 * clamped);
        return 'hsl(' + hue + ' 72% ' + light + '%)';
    }}

    function colorInfo(columnName) {{
        const column = metadataColumn(columnName);
        if (!column) {{
            return null;
        }}
        if (column.type === 'int' || column.type === 'float') {{
            const values = [];
            state.metadataByNodeIndex.forEach(row => {{
                const columnIndex = state.metadataColumns.findIndex(item => item.name === columnName);
                const value = row?.[columnIndex];
                if (typeof value === 'number' && !Number.isNaN(value)) {{
                    values.push(value);
                }}
            }});
            return {{type: 'numeric', min: Math.min(...values), max: Math.max(...values)}};
        }}
        return {{type: 'categorical'}};
    }}

    function nodeColor(nodeIndex) {{
        const columnName = currentColorField();
        const info = colorInfo(columnName);
        const value = metadataValue(nodeIndex, columnName);
        if (!info) {{
            return '#d88f3d';
        }}
        if (info.type === 'numeric') {{
            return numericColor(value, info.min, info.max);
        }}
        return categoricalColor(value);
    }}

    function labelForComponent(componentId) {{
        const hierarchyNode = state.bundle.graph.hierarchy.nodes[componentId];
        if (hierarchyNode.kind === 'leaf') {{
            const field = currentLabelField();
            return formatValue(metadataValue(hierarchyNode.node_index, field)) === '—' ? nodeId(hierarchyNode.node_index) : formatValue(metadataValue(hierarchyNode.node_index, field));
        }}
        return hierarchyNode.size.toLocaleString() + ' nodes';
    }}

    function activeClusterAssignments(activeClusterIds) {{
        const assignments = new Array(state.bundle.graph.nodes.length).fill(-1);
        activeClusterIds.forEach(componentId => {{
            const members = componentMembers(componentId);
            members.forEach(nodeIndex => {{
                assignments[nodeIndex] = componentId;
            }});
        }});
        return assignments;
    }}

    function mstLinksForActiveClusters(activeClusterIds, minClusterSize, leafPruningOnly) {{
        const assignments = activeClusterAssignments(activeClusterIds);
        const linkMap = new Map();

        state.bundle.graph.mst_edges.forEach(edge => {{
            const sourceComponentId = assignments[edge[0]];
            const targetComponentId = assignments[edge[1]];
            if (sourceComponentId < 0 || targetComponentId < 0 || sourceComponentId === targetComponentId) {{
                return;
            }}
            const leftId = sourceComponentId < targetComponentId ? sourceComponentId : targetComponentId;
            const rightId = sourceComponentId < targetComponentId ? targetComponentId : sourceComponentId;
            const key = leftId + ':' + rightId;
            if (!linkMap.has(key)) {{
                linkMap.set(key, {{sourceId: leftId, targetId: rightId, weight: edge[2]}});
            }}
        }});

        if (minClusterSize <= 1) {{
            return {{
                visibleIds: [...activeClusterIds],
                links: Array.from(linkMap.values()),
                hiddenNodes: 0,
            }};
        }}

        const allLinks = Array.from(linkMap.values());
        if (!leafPruningOnly) {{
            const visibleIds = activeClusterIds.filter(componentId => state.bundle.graph.hierarchy.nodes[componentId].size >= minClusterSize);
            const visibleSet = new Set(visibleIds);
            const hiddenNodes = activeClusterIds.reduce((sum, componentId) => {{
                if (visibleSet.has(componentId)) {{
                    return sum;
                }}
                return sum + state.bundle.graph.hierarchy.nodes[componentId].size;
            }}, 0);
            return {{
                visibleIds,
                links: allLinks.filter(link => visibleSet.has(link.sourceId) && visibleSet.has(link.targetId)),
                hiddenNodes,
            }};
        }}

        const visibleSet = new Set(activeClusterIds);
        const adjacency = new Map();
        activeClusterIds.forEach(componentId => adjacency.set(componentId, []));
        allLinks.forEach(link => {{
            adjacency.get(link.sourceId)?.push(link.targetId);
            adjacency.get(link.targetId)?.push(link.sourceId);
        }});

        const degree = new Map();
        activeClusterIds.forEach(componentId => degree.set(componentId, (adjacency.get(componentId) || []).length));
        let removedAny = true;
        while (removedAny) {{
            removedAny = false;
            const pruneIds = [];
            visibleSet.forEach(componentId => {{
                const size = state.bundle.graph.hierarchy.nodes[componentId].size;
                const currentDegree = degree.get(componentId) || 0;
                if (currentDegree <= 1 && size < minClusterSize) {{
                    pruneIds.push(componentId);
                }}
            }});
            if (pruneIds.length === 0) {{
                break;
            }}
            removedAny = true;
            pruneIds.forEach(componentId => {{
                visibleSet.delete(componentId);
            }});
            pruneIds.forEach(componentId => {{
                (adjacency.get(componentId) || []).forEach(neighborId => {{
                    if (!visibleSet.has(neighborId)) {{
                        return;
                    }}
                    degree.set(neighborId, Math.max(0, (degree.get(neighborId) || 0) - 1));
                }});
                degree.set(componentId, 0);
            }});
        }}

        const visibleIds = activeClusterIds.filter(componentId => visibleSet.has(componentId));
        const hiddenNodes = activeClusterIds.reduce((sum, componentId) => {{
            if (visibleSet.has(componentId)) {{
                return sum;
            }}
            return sum + state.bundle.graph.hierarchy.nodes[componentId].size;
        }}, 0);

        return {{
            visibleIds,
            links: allLinks.filter(link => visibleSet.has(link.sourceId) && visibleSet.has(link.targetId)),
            hiddenNodes,
        }};
    }}

    function treeCenter(nodeIds, adjacency) {{
        if (nodeIds.length <= 2) {{
            return nodeIds[0];
        }}
        const degree = new Map();
        nodeIds.forEach(nodeId => {{
            degree.set(nodeId, (adjacency.get(nodeId) || []).length);
        }});
        let leaves = nodeIds.filter(nodeId => (degree.get(nodeId) || 0) <= 1);
        let remaining = nodeIds.length;
        while (remaining > 2 && leaves.length > 0) {{
            remaining -= leaves.length;
            const nextLeaves = [];
            leaves.forEach(leafId => {{
                (adjacency.get(leafId) || []).forEach(neighborId => {{
                    if (!degree.has(neighborId)) {{
                        return;
                    }}
                    degree.set(neighborId, degree.get(neighborId) - 1);
                    if (degree.get(neighborId) === 1) {{
                        nextLeaves.push(neighborId);
                    }}
                }});
                degree.delete(leafId);
            }});
            leaves = nextLeaves;
        }}
        const candidates = degree.size > 0 ? Array.from(degree.keys()) : nodeIds;
        candidates.sort((leftId, rightId) => {{
            const leftSize = state.bundle.graph.hierarchy.nodes[leftId].size;
            const rightSize = state.bundle.graph.hierarchy.nodes[rightId].size;
            return rightSize - leftSize || leftId - rightId;
        }});
        return candidates[0];
    }}

    function rootedTree(rootId, adjacency) {{
        const parent = new Map([[rootId, null]]);
        const order = [rootId];
        for (let index = 0; index < order.length; index++) {{
            const nodeId = order[index];
            (adjacency.get(nodeId) || []).forEach(neighborId => {{
                if (parent.has(neighborId)) {{
                    return;
                }}
                parent.set(neighborId, nodeId);
                order.push(neighborId);
            }});
        }}

        const children = new Map();
        order.forEach(nodeId => children.set(nodeId, []));
        for (let index = 1; index < order.length; index++) {{
            const nodeId = order[index];
            children.get(parent.get(nodeId)).push(nodeId);
        }}
        children.forEach((childIds, nodeId) => {{
            childIds.sort((leftId, rightId) => {{
                const leftNode = state.bundle.graph.hierarchy.nodes[leftId];
                const rightNode = state.bundle.graph.hierarchy.nodes[rightId];
                return leftNode.leaf_start - rightNode.leaf_start || rightNode.size - leftNode.size || leftId - rightId;
            }});
        }});

        return {{parent, order, children}};
    }}

    function componentRadiusForSize(size) {{
        return Math.max(7, Math.min(52, 10 + Math.sqrt(size) * 4.7));
    }}

    function clusterGraphComponents(visibleIds, links) {{
        const adjacency = new Map();
        visibleIds.forEach(nodeId => adjacency.set(nodeId, []));
        links.forEach(link => {{
            adjacency.get(link.sourceId)?.push(link.targetId);
            adjacency.get(link.targetId)?.push(link.sourceId);
        }});

        const components = [];
        const seen = new Set();
        visibleIds.forEach(nodeId => {{
            if (seen.has(nodeId)) {{
                return;
            }}
            const stack = [nodeId];
            const componentIds = [];
            seen.add(nodeId);
            while (stack.length > 0) {{
                const currentId = stack.pop();
                componentIds.push(currentId);
                (adjacency.get(currentId) || []).forEach(neighborId => {{
                    if (seen.has(neighborId)) {{
                        return;
                    }}
                    seen.add(neighborId);
                    stack.push(neighborId);
                }});
            }}
            components.push(componentIds);
        }});

        components.sort((leftIds, rightIds) => {{
            const leftStart = Math.min(...leftIds.map(nodeId => state.bundle.graph.hierarchy.nodes[nodeId].leaf_start));
            const rightStart = Math.min(...rightIds.map(nodeId => state.bundle.graph.hierarchy.nodes[nodeId].leaf_start));
            return leftStart - rightStart || rightIds.length - leftIds.length;
        }});
        return {{adjacency, components}};
    }}

    function packLayouts(componentLayouts, options = {{}}) {{
        const gapX = options.gapX ?? 120;
        const gapY = options.gapY ?? 120;
        const rowTargetWidth = options.rowTargetWidth ?? 2200;
        const outerPadding = options.outerPadding ?? 72;
        const packed = [];
        let cursorX = outerPadding;
        let cursorY = outerPadding;
        let rowHeight = 0;

        componentLayouts.forEach(component => {{
            if (!component || component.items.length === 0) {{
                return;
            }}
            if (cursorX > outerPadding && cursorX + component.width > rowTargetWidth) {{
                cursorX = outerPadding;
                cursorY += rowHeight + gapY;
                rowHeight = 0;
            }}
            component.items.forEach(item => {{
                packed.push({{
                    componentId: item.componentId,
                    x: item.x + cursorX,
                    y: item.y + cursorY,
                    radius: item.radius,
                }});
            }});
            cursorX += component.width + gapX;
            rowHeight = Math.max(rowHeight, component.height);
        }});
        return packed;
    }}

    function normalizeComponentLayout(items, padding = 20) {{
        if (items.length === 0) {{
            return {{items: [], width: 0, height: 0}};
        }}
        let minX = Infinity;
        let minY = Infinity;
        let maxX = -Infinity;
        let maxY = -Infinity;
        items.forEach(item => {{
            minX = Math.min(minX, item.x - item.radius);
            minY = Math.min(minY, item.y - item.radius);
            maxX = Math.max(maxX, item.x + item.radius);
            maxY = Math.max(maxY, item.y + item.radius);
        }});
        const normalizedItems = items.map(item => ({{
            ...item,
            x: item.x - minX + padding,
            y: item.y - minY + padding,
        }}));
        return {{
            items: normalizedItems,
            width: (maxX - minX) + (padding * 2),
            height: (maxY - minY) + (padding * 2),
        }};
    }}

    function seededUnit(componentId, salt) {{
        const raw = Math.sin((componentId + 1) * 12.9898 + salt * 78.233) * 43758.5453;
        return raw - Math.floor(raw);
    }}

    function simulateComponentLayout(componentIds, adjacency, options = {{}}) {{
        const radii = new Map(componentIds.map(nodeId => [nodeId, componentRadiusForSize(state.bundle.graph.hierarchy.nodes[nodeId].size)]));
        const tree = rootedTree(treeCenter(componentIds, adjacency), adjacency);
        const positions = new Map();
        const velocities = new Map();
        const anchors = new Map();
        const baseSpacing = options.baseSpacing ?? 110;
        const angleScale = options.angleScale ?? 0.5;
        const damping = options.damping ?? 0.8;
        const repulsion = options.repulsion ?? 6000;
        const spring = options.spring ?? 0.07;
        const gravity = options.gravity ?? 0.008;
        const anchorStrength = options.anchorStrength ?? 0;
        const iterations = options.iterations ?? 110;
        const preferredEdgeLength = options.preferredEdgeLength ?? 120;
        const collisionStrength = options.collisionStrength ?? 0.28;

        tree.order.forEach((nodeId, index) => {{
            const depth = (() => {{
                let currentDepth = 0;
                let currentId = nodeId;
                while (tree.parent.get(currentId) !== null) {{
                    currentDepth += 1;
                    currentId = tree.parent.get(currentId);
                }}
                return currentDepth;
            }})();
            const angle = ((index / Math.max(1, tree.order.length)) * Math.PI * 2 * angleScale) + (seededUnit(nodeId, 1) - 0.5) * 0.35;
            const radius = depth * baseSpacing;
            const anchor = options.useTreeSeed === false
                ? {{x: Math.cos(angle) * radius, y: Math.sin(angle) * radius}}
                : {{x: depth * baseSpacing, y: ((state.bundle.graph.hierarchy.nodes[nodeId].leaf_start ?? index) - componentIds[0]) * 12 + (seededUnit(nodeId, 2) - 0.5) * 24}};
            anchors.set(nodeId, anchor);
            positions.set(nodeId, {{
                x: anchor.x + (seededUnit(nodeId, 3) - 0.5) * 18,
                y: anchor.y + (seededUnit(nodeId, 4) - 0.5) * 18,
            }});
            velocities.set(nodeId, {{x: 0, y: 0}});
        }});

        const edgePairs = [];
        componentIds.forEach(nodeId => {{
            (adjacency.get(nodeId) || []).forEach(neighborId => {{
                if (nodeId < neighborId) {{
                    edgePairs.push([nodeId, neighborId]);
                }}
            }});
        }});

        for (let iteration = 0; iteration < iterations; iteration++) {{
            const forces = new Map(componentIds.map(nodeId => [nodeId, {{x: 0, y: 0}}]));
            for (let leftIndex = 0; leftIndex < componentIds.length; leftIndex++) {{
                const leftId = componentIds[leftIndex];
                const leftPos = positions.get(leftId);
                for (let rightIndex = leftIndex + 1; rightIndex < componentIds.length; rightIndex++) {{
                    const rightId = componentIds[rightIndex];
                    const rightPos = positions.get(rightId);
                    let dx = rightPos.x - leftPos.x;
                    let dy = rightPos.y - leftPos.y;
                    let distSq = (dx * dx) + (dy * dy);
                    if (distSq < 1e-6) {{
                        dx = (seededUnit(leftId + rightId, iteration + 5) - 0.5) * 0.01;
                        dy = (seededUnit(leftId + rightId, iteration + 6) - 0.5) * 0.01;
                        distSq = (dx * dx) + (dy * dy);
                    }}
                    const dist = Math.sqrt(distSq);
                    const repelForce = repulsion / distSq;
                    const overlap = ((radii.get(leftId) || 0) + (radii.get(rightId) || 0) + 12) - dist;
                    const collisionForce = overlap > 0 ? overlap * collisionStrength : 0;
                    const forceX = (dx / dist) * (repelForce + collisionForce);
                    const forceY = (dy / dist) * (repelForce + collisionForce);
                    forces.get(leftId).x -= forceX;
                    forces.get(leftId).y -= forceY;
                    forces.get(rightId).x += forceX;
                    forces.get(rightId).y += forceY;
                }}
            }}

            edgePairs.forEach(([leftId, rightId]) => {{
                const leftPos = positions.get(leftId);
                const rightPos = positions.get(rightId);
                let dx = rightPos.x - leftPos.x;
                let dy = rightPos.y - leftPos.y;
                let dist = Math.hypot(dx, dy);
                if (dist < 1e-6) {{
                    dist = 1e-6;
                    dx = preferredEdgeLength;
                    dy = 0;
                }}
                const delta = dist - preferredEdgeLength;
                const force = spring * delta;
                const forceX = (dx / dist) * force;
                const forceY = (dy / dist) * force;
                forces.get(leftId).x += forceX;
                forces.get(leftId).y += forceY;
                forces.get(rightId).x -= forceX;
                forces.get(rightId).y -= forceY;
            }});

            componentIds.forEach(nodeId => {{
                const position = positions.get(nodeId);
                const velocity = velocities.get(nodeId);
                const force = forces.get(nodeId);
                const anchor = anchors.get(nodeId);
                force.x += (anchor.x - position.x) * anchorStrength;
                force.y += (anchor.y - position.y) * anchorStrength;
                force.x += -position.x * gravity;
                force.y += -position.y * gravity;
                velocity.x = (velocity.x + force.x) * damping;
                velocity.y = (velocity.y + force.y) * damping;
                position.x += velocity.x;
                position.y += velocity.y;
            }});
        }}

        return normalizeComponentLayout(componentIds.map(nodeId => {{
            const position = positions.get(nodeId);
            return {{componentId: nodeId, x: position.x, y: position.y, radius: radii.get(nodeId) || 10}};
        }}), 24);
    }}

    function forceDirectedForestLayout(visibleIds, links) {{
        if (visibleIds.length === 0) {{
            clusterCanvas.height = 760;
            return [];
        }}
        clusterCanvas.height = Math.max(760, Math.min(1180, Math.round(window.innerHeight * 0.8)));
        const {{adjacency, components}} = clusterGraphComponents(visibleIds, links);
        const componentLayouts = components.map(componentIds => simulateComponentLayout(componentIds, adjacency, {{
            useTreeSeed: false,
            baseSpacing: 82,
            repulsion: 7200,
            spring: 0.08,
            gravity: 0.010,
            damping: 0.82,
            preferredEdgeLength: 118,
            iterations: Math.min(140, 68 + componentIds.length),
            collisionStrength: 0.32,
            anchorStrength: 0.005,
            angleScale: 1,
        }}));
        return packLayouts(componentLayouts, {{gapX: 130, gapY: 130, rowTargetWidth: 2300, outerPadding: 72}});
    }}

    function organicForestLayout(visibleIds, links) {{
        if (visibleIds.length === 0) {{
            clusterCanvas.height = 760;
            return [];
        }}
        clusterCanvas.height = Math.max(760, Math.min(1180, Math.round(window.innerHeight * 0.8)));
        const {{adjacency, components}} = clusterGraphComponents(visibleIds, links);
        const componentLayouts = components.map(componentIds => simulateComponentLayout(componentIds, adjacency, {{
            useTreeSeed: true,
            baseSpacing: 88,
            repulsion: 5600,
            spring: 0.11,
            gravity: 0.012,
            damping: 0.84,
            preferredEdgeLength: 104,
            iterations: Math.min(160, 84 + componentIds.length),
            collisionStrength: 0.34,
            anchorStrength: 0.028,
            angleScale: 0.65,
        }}));
        return packLayouts(componentLayouts, {{gapX: 115, gapY: 115, rowTargetWidth: 2200, outerPadding: 72}});
    }}

    function tidyForestLayout(visibleIds, links) {{
        if (visibleIds.length === 0) {{
            clusterCanvas.height = 760;
            return [];
        }}

        const {{adjacency, components}} = clusterGraphComponents(visibleIds, links);

        clusterCanvas.height = Math.max(760, Math.min(1180, Math.round(window.innerHeight * 0.8)));
        const outerPadding = 58;
        const componentGap = 110;
        const siblingGap = 18;
        const levelGap = 152;
        const layout = [];
        let currentX = outerPadding;

        components.forEach(componentIds => {{
            const rootId = treeCenter(componentIds, adjacency);
            const tree = rootedTree(rootId, adjacency);
            const depthByNode = new Map([[rootId, 0]]);
            tree.order.forEach(nodeId => {{
                (tree.children.get(nodeId) || []).forEach(childId => {{
                    depthByNode.set(childId, (depthByNode.get(nodeId) || 0) + 1);
                }});
            }});
            const radii = new Map(componentIds.map(nodeId => [nodeId, componentRadiusForSize(state.bundle.graph.hierarchy.nodes[nodeId].size)]));
            const subtreeSpan = new Map();
            [...tree.order].reverse().forEach(nodeId => {{
                const childIds = tree.children.get(nodeId) || [];
                const nodeSpan = Math.max(16, Math.min(30, (radii.get(nodeId) || 12) * 0.72 + 8));
                if (childIds.length === 0) {{
                    subtreeSpan.set(nodeId, nodeSpan);
                    return;
                }}
                const childrenSpan = childIds.reduce((sum, childId) => sum + subtreeSpan.get(childId), 0) + (Math.max(0, childIds.length - 1) * siblingGap);
                subtreeSpan.set(nodeId, Math.max(nodeSpan, childrenSpan));
            }});

            const positionById = new Map();
            function placeNode(nodeId, topY) {{
                const nodeSpan = subtreeSpan.get(nodeId) || 26;
                const childIds = tree.children.get(nodeId) || [];
                const x = currentX + (depthByNode.get(nodeId) || 0) * levelGap;
                if (childIds.length === 0) {{
                    positionById.set(nodeId, {{
                        componentId: nodeId,
                        x,
                        y: topY + (nodeSpan / 2),
                        radius: radii.get(nodeId) || 10,
                    }});
                    return;
                }}

                const childrenSpan = childIds.reduce((sum, childId) => sum + subtreeSpan.get(childId), 0) + (Math.max(0, childIds.length - 1) * siblingGap);
                let childCursor = topY + Math.max(0, (nodeSpan - childrenSpan) / 2);
                const childCenters = [];
                childIds.forEach(childId => {{
                    placeNode(childId, childCursor);
                    childCenters.push(positionById.get(childId).y);
                    childCursor += subtreeSpan.get(childId) + siblingGap;
                }});
                const centerY = childCenters.reduce((sum, value) => sum + value, 0) / Math.max(1, childCenters.length);
                positionById.set(nodeId, {{
                    componentId: nodeId,
                    x,
                    y: centerY,
                    radius: radii.get(nodeId) || 10,
                }});
            }}

            placeNode(rootId, outerPadding);
            const componentLayout = tree.order.map(nodeId => positionById.get(nodeId));
            layout.push(...componentLayout);
            const maxDepth = tree.order.reduce((maxValue, nodeId) => Math.max(maxValue, depthByNode.get(nodeId) || 0), 0);
            const componentWidth = (maxDepth * levelGap) + Math.max(...componentLayout.map(item => item.radius), 10) * 2 + outerPadding;
            currentX += componentWidth + componentGap;
        }});

        return layout;
    }}

    function layoutCacheKey(visibleIds, links, algorithm) {{
        const linkKey = links.map(link => link.sourceId + ':' + link.targetId + ':' + Number(link.weight).toFixed(4)).join('|');
        return algorithm + '::' + visibleIds.join(',') + '::' + linkKey;
    }}

    function computeVisibleLayout(visibleIds, links, algorithm) {{
        const key = layoutCacheKey(visibleIds, links, algorithm);
        if (state.layoutCache.has(key)) {{
            return state.layoutCache.get(key).map(item => ({{...item}}));
        }}
        const layout = (() => {{
            if (algorithm === 'force') {{
                return forceDirectedForestLayout(visibleIds, links);
            }}
            if (algorithm === 'organic') {{
                return organicForestLayout(visibleIds, links);
            }}
            return tidyForestLayout(visibleIds, links);
        }})();
        state.layoutCache.set(key, layout.map(item => ({{...item}})));
        if (state.layoutCache.size > 24) {{
            const oldestKey = state.layoutCache.keys().next().value;
            state.layoutCache.delete(oldestKey);
        }}
        return layout;
    }}

    function drawSplitChart() {{
        splitContext.clearRect(0, 0, splitCanvas.width, splitCanvas.height);
        splitContext.fillStyle = '#fffaf4';
        splitContext.fillRect(0, 0, splitCanvas.width, splitCanvas.height);

        if (!state.bundle) {{
            splitContext.fillStyle = '#5c6a70';
            splitContext.font = '18px Georgia';
            splitContext.fillText('Load a bundle to render split events.', 30, 50);
            return;
        }}

        const events = state.bundle.graph.merge_event_series;
        const margin = {{top: 26, right: 30, bottom: 62, left: 76}};
        const width = splitCanvas.width - margin.left - margin.right;
        const height = splitCanvas.height - margin.top - margin.bottom;

        splitContext.strokeStyle = 'rgba(92,106,112,0.35)';
        splitContext.lineWidth = 1;
        splitContext.beginPath();
        splitContext.moveTo(margin.left, margin.top + height);
        splitContext.lineTo(margin.left + width, margin.top + height);
        splitContext.moveTo(margin.left, margin.top);
        splitContext.lineTo(margin.left, margin.top + height);
        splitContext.stroke();

        splitContext.fillStyle = '#5c6a70';
        splitContext.font = '12px Georgia';
        splitContext.textAlign = 'center';
        splitContext.textBaseline = 'alphabetic';
        splitContext.fillText('Threshold', margin.left + (width / 2), splitCanvas.height - 14);
        splitContext.save();
        splitContext.translate(20, margin.top + (height / 2));
        splitContext.rotate(-Math.PI / 2);
        splitContext.textAlign = 'center';
        splitContext.textBaseline = 'alphabetic';
        splitContext.fillText('Split impact', 0, 0);
        splitContext.restore();

        if (events.length === 0) {{
            splitContext.fillStyle = '#5c6a70';
            splitContext.font = '18px Georgia';
            splitContext.fillText('No split events in this bundle.', 30, 50);
            return;
        }}

        const impacts = events.map(event => event.merge_impact);
        const maxImpact = Math.max(...impacts, 1);
        const minThreshold = Math.min(...events.map(event => event.threshold_value));
        const maxThreshold = Math.max(...events.map(event => event.threshold_value));
        const thresholdSpan = Math.max(1e-9, maxThreshold - minThreshold || 1);
        const tickLength = 6;

        const xFor = value => margin.left + ((value - minThreshold) / thresholdSpan) * width;
        const yFor = value => margin.top + height - (value / maxImpact) * height;

        splitContext.fillStyle = '#5c6a70';
        splitContext.font = '11px Georgia';
        splitContext.textBaseline = 'top';
        splitContext.textAlign = 'center';
        const xTickCount = Math.min(6, Math.max(2, events.length > 1 ? 5 : 2));
        for (let tickIndex = 0; tickIndex < xTickCount; tickIndex++) {{
            const fraction = xTickCount === 1 ? 0 : tickIndex / (xTickCount - 1);
            const thresholdValue = minThreshold + (thresholdSpan * fraction);
            const x = xFor(thresholdValue);
            splitContext.strokeStyle = 'rgba(92,106,112,0.4)';
            splitContext.lineWidth = 1;
            splitContext.beginPath();
            splitContext.moveTo(x, margin.top + height);
            splitContext.lineTo(x, margin.top + height + tickLength);
            splitContext.stroke();
            splitContext.fillText(formatValue(thresholdValue), x, margin.top + height + tickLength + 4);
        }}

        splitContext.textAlign = 'right';
        splitContext.textBaseline = 'middle';
        const yTickCount = 5;
        for (let tickIndex = 0; tickIndex < yTickCount; tickIndex++) {{
            const fraction = tickIndex / (yTickCount - 1);
            const impactValue = maxImpact * fraction;
            const y = yFor(impactValue);
            splitContext.strokeStyle = 'rgba(92,106,112,0.4)';
            splitContext.lineWidth = 1;
            splitContext.beginPath();
            splitContext.moveTo(margin.left - tickLength, y);
            splitContext.lineTo(margin.left, y);
            splitContext.stroke();
            splitContext.fillText(formatValue(impactValue), margin.left - tickLength - 6, y);
        }}

        splitContext.strokeStyle = '#c8553d';
        splitContext.lineWidth = 1.5;
        splitContext.beginPath();
        events.forEach((event, index) => {{
            const x = xFor(event.threshold_value);
            const y = yFor(event.merge_impact);
            if (index === 0) {{
                splitContext.moveTo(x, y);
            }} else {{
                splitContext.lineTo(x, y);
            }}
        }});
        splitContext.stroke();

        events.forEach(event => {{
            const x = xFor(event.threshold_value);
            const y = yFor(event.merge_impact);
            splitContext.fillStyle = '#f3cdb9';
            splitContext.beginPath();
            splitContext.arc(x, y, 4, 0, Math.PI * 2);
            splitContext.fill();
            splitContext.strokeStyle = '#8a3b2c';
            splitContext.stroke();
        }});

        const stop = currentSliderStop();
        if (stop) {{
            const x = stop.threshold_value === null ? margin.left + width : xFor(stop.threshold_value);
            splitContext.strokeStyle = '#1e2a2f';
            splitContext.setLineDash([6, 6]);
            splitContext.beginPath();
            splitContext.moveTo(x, margin.top);
            splitContext.lineTo(x, margin.top + height);
            splitContext.stroke();
            splitContext.setLineDash([]);
            splitContext.fillStyle = '#e29b4b';
            splitContext.beginPath();
            splitContext.arc(x, margin.top + height + 10, 5, 0, Math.PI * 2);
            splitContext.fill();
        }}
    }}

    function worldToScreenPoint(x, y) {{
        return {{
            x: (x * state.viewTransform.scale) + state.viewTransform.offsetX,
            y: (y * state.viewTransform.scale) + state.viewTransform.offsetY,
        }};
    }}

    function screenToWorldPoint(x, y) {{
        return {{
            x: (x - state.viewTransform.offsetX) / Math.max(state.viewTransform.scale, 1e-9),
            y: (y - state.viewTransform.offsetY) / Math.max(state.viewTransform.scale, 1e-9),
        }};
    }}

    function canvasCoordinatesFromEvent(event) {{
        const rect = clusterCanvas.getBoundingClientRect();
        const scaleX = clusterCanvas.width / rect.width;
        const scaleY = clusterCanvas.height / rect.height;
        return {{
            x: (event.clientX - rect.left) * scaleX,
            y: (event.clientY - rect.top) * scaleY,
        }};
    }}

    function clusterLayoutBounds() {{
        if (state.visibleLayout.length === 0) {{
            return null;
        }}
        let minX = Infinity;
        let minY = Infinity;
        let maxX = -Infinity;
        let maxY = -Infinity;
        state.visibleLayout.forEach(item => {{
            minX = Math.min(minX, item.x - item.radius);
            minY = Math.min(minY, item.y - item.radius);
            maxX = Math.max(maxX, item.x + item.radius);
            maxY = Math.max(maxY, item.y + item.radius);
        }});
        return {{minX, minY, maxX, maxY}};
    }}

    function fitClusterViewToLayout() {{
        const bounds = clusterLayoutBounds();
        if (!bounds) {{
            state.viewTransform.scale = 1;
            state.viewTransform.offsetX = 0;
            state.viewTransform.offsetY = 0;
            return;
        }}
        const padding = 42;
        const width = Math.max(1, bounds.maxX - bounds.minX);
        const height = Math.max(1, bounds.maxY - bounds.minY);
        const scale = Math.max(
            state.viewTransform.minScale,
            Math.min(
                state.viewTransform.maxScale,
                Math.min((clusterCanvas.width - (padding * 2)) / width, (clusterCanvas.height - (padding * 2)) / height),
            ),
        );
        state.viewTransform.scale = scale;
        state.viewTransform.offsetX = padding + ((clusterCanvas.width - (padding * 2) - (width * scale)) / 2) - (bounds.minX * scale);
        state.viewTransform.offsetY = padding + ((clusterCanvas.height - (padding * 2) - (height * scale)) / 2) - (bounds.minY * scale);
    }}

    function drawBadge(text, x, y) {{
        clusterContext.save();
        clusterContext.font = '11px Georgia';
        const textWidth = clusterContext.measureText(text).width;
        const badgeWidth = textWidth + 12;
        const badgeHeight = 18;
        clusterContext.fillStyle = 'rgba(255, 250, 242, 0.92)';
        clusterContext.strokeStyle = 'rgba(92, 106, 112, 0.32)';
        clusterContext.lineWidth = 1;
        clusterContext.beginPath();
        clusterContext.roundRect(x - (badgeWidth / 2), y - (badgeHeight / 2), badgeWidth, badgeHeight, 8);
        clusterContext.fill();
        clusterContext.stroke();
        clusterContext.fillStyle = '#334147';
        clusterContext.textAlign = 'center';
        clusterContext.textBaseline = 'middle';
        clusterContext.fillText(text, x, y + 0.5);
        clusterContext.restore();
    }}

    function renderClusterView() {{
        clusterContext.clearRect(0, 0, clusterCanvas.width, clusterCanvas.height);
        clusterContext.fillStyle = '#fffaf4';
        clusterContext.fillRect(0, 0, clusterCanvas.width, clusterCanvas.height);

        if (!state.bundle) {{
            clusterContext.fillStyle = '#5c6a70';
            clusterContext.font = '18px Georgia';
            clusterContext.fillText('Load a bundle to render cluster bubbles.', 30, 50);
            return;
        }}

        clusterContext.save();
        clusterContext.translate(state.viewTransform.offsetX, state.viewTransform.offsetY);
        clusterContext.scale(state.viewTransform.scale, state.viewTransform.scale);

        state.splitLinks.forEach(link => {{
            clusterContext.strokeStyle = 'rgba(92, 106, 112, 0.42)';
            clusterContext.lineWidth = 1.6 / Math.max(state.viewTransform.scale, 0.1);
            clusterContext.beginPath();
            clusterContext.moveTo(link.left.x, link.left.y);
            clusterContext.lineTo(link.right.x, link.right.y);
            clusterContext.stroke();
        }});

        state.visibleLayout.forEach(item => {{
            const component = state.bundle.graph.hierarchy.nodes[item.componentId];
            const members = componentMembers(item.componentId);
            const isSelected = members.some(nodeIndex => state.selectedNodeIndices.has(nodeIndex));

            clusterContext.fillStyle = 'rgba(248, 243, 235, 0.95)';
            clusterContext.strokeStyle = isSelected ? '#1e2a2f' : '#c8b8a6';
            clusterContext.lineWidth = (isSelected ? 3 : 1.5) / Math.max(state.viewTransform.scale, 0.1);
            clusterContext.beginPath();
            clusterContext.arc(item.x, item.y, item.radius, 0, Math.PI * 2);
            clusterContext.fill();
            clusterContext.stroke();

            const sampleCount = Math.min(component.size, 140);
            for (let sampleIndex = 0; sampleIndex < sampleCount; sampleIndex++) {{
                const memberIndex = members[Math.floor(sampleIndex * members.length / sampleCount)];
                const angle = sampleIndex * 2.399963229728653;
                const distance = component.size === 1 ? 0 : Math.sqrt((sampleIndex + 0.5) / sampleCount) * Math.max(0, item.radius - 8);
                const dotX = item.x + Math.cos(angle) * distance;
                const dotY = item.y + Math.sin(angle) * distance;
                clusterContext.fillStyle = nodeColor(memberIndex);
                clusterContext.beginPath();
                clusterContext.arc(dotX, dotY, component.size === 1 ? 4.8 / Math.max(state.viewTransform.scale, 0.35) : 2.6 / Math.max(state.viewTransform.scale, 0.35), 0, Math.PI * 2);
                clusterContext.fill();
            }}
        }});
        clusterContext.restore();

        if (showEdgeScoresEnabled() && state.viewTransform.scale >= 0.16) {{
            state.splitLinks.forEach(link => {{
                const left = worldToScreenPoint(link.left.x, link.left.y);
                const right = worldToScreenPoint(link.right.x, link.right.y);
                const dx = right.x - left.x;
                const dy = right.y - left.y;
                if (Math.hypot(dx, dy) < 46) {{
                    return;
                }}
                drawBadge(formatValue(link.threshold), left.x + (dx / 2), left.y + (dy / 2) - 8);
            }});
        }}

        state.visibleLayout.forEach(item => {{
            const component = state.bundle.graph.hierarchy.nodes[item.componentId];
            const screenPoint = worldToScreenPoint(item.x, item.y);
            const screenRadius = item.radius * state.viewTransform.scale;
            let labelText = '';
            let labelY = screenPoint.y + 4;
            let font = '600 12px Georgia';
            if (document.getElementById('show-labels').checked && (component.size <= 8 || component.kind === 'leaf')) {{
                labelText = labelForComponent(item.componentId).slice(0, 20);
                labelY = screenPoint.y + screenRadius + 14;
                font = component.size === 1 ? '12px Georgia' : '11px Georgia';
            }} else if (showNodeCountsEnabled()) {{
                labelText = component.size.toLocaleString();
            }}
            if (!labelText || state.viewTransform.scale < 0.11) {{
                return;
            }}
            clusterContext.fillStyle = labelText === component.size.toLocaleString() ? '#5c6a70' : '#1e2a2f';
            clusterContext.font = font;
            clusterContext.textAlign = 'center';
            clusterContext.textBaseline = 'middle';
            clusterContext.fillText(labelText, screenPoint.x, labelY);
        }});

        if (state.selectionBox) {{
            const box = state.selectionBox;
            clusterContext.save();
            clusterContext.fillStyle = 'rgba(200, 85, 61, 0.10)';
            clusterContext.strokeStyle = 'rgba(200, 85, 61, 0.62)';
            clusterContext.setLineDash([8, 6]);
            clusterContext.lineWidth = 1.5;
            clusterContext.fillRect(box.left, box.top, box.width, box.height);
            clusterContext.strokeRect(box.left, box.top, box.width, box.height);
            clusterContext.restore();
        }}
    }}

    function drawClusterView(resetView = true) {{
        if (!state.bundle) {{
            renderClusterView();
            return;
        }}

        const minClusterSize = Math.max(1, Number(document.getElementById('min-cluster-size').value) || 1);
        const thresholdValue = selectedThresholdValue();
        const activeClusterIds = activeClustersAtThreshold(thresholdValue);
        const visibleGraph = mstLinksForActiveClusters(activeClusterIds, minClusterSize, leafPruningOnlyEnabled());
        const visibleLayout = computeVisibleLayout(visibleGraph.visibleIds, visibleGraph.links, currentLayoutAlgorithm());
        const layoutById = new Map(visibleLayout.map(item => [item.componentId, item]));
        const links = visibleGraph.links
            .map(link => {{
                const left = layoutById.get(link.sourceId);
                const right = layoutById.get(link.targetId);
                if (!left || !right) {{
                    return null;
                }}
                return {{left, right, threshold: link.weight}};
            }})
            .filter(Boolean);

        state.activeClusters = activeClusterIds;
        state.visibleClusters = visibleGraph.visibleIds;
        state.visibleLayout = visibleLayout;
        state.splitLinks = links;

        const hidden = visibleGraph.hiddenNodes;
        const shown = state.bundle.graph.nodes.length - hidden;
        document.getElementById('hidden-summary').textContent = hidden.toLocaleString() + ' nodes hidden by minimum cluster size';
        document.getElementById('stat-clusters').textContent = visibleLayout.length.toLocaleString();
        document.getElementById('stat-links').textContent = links.length.toLocaleString();
        document.getElementById('stat-shown-nodes').textContent = shown.toLocaleString();
        document.getElementById('stat-hidden-nodes').textContent = hidden.toLocaleString();
        if (resetView) {{
            fitClusterViewToLayout();
        }}
        renderClusterView();
    }}

    function updateMetadataTable() {{
        const selected = Array.from(state.selectedNodeIndices).sort((left, right) => left - right);
        const displayedNodeIndices = selected.length > 0
            ? selected
            : state.bundle.graph.nodes.map((_, nodeIndex) => nodeIndex);
        const thead = document.querySelector('#metadata-table thead');
        const tbody = document.querySelector('#metadata-table tbody');
        thead.innerHTML = '';
        tbody.innerHTML = '';

        const headerRow = document.createElement('tr');
        ['node_id', ...state.metadataColumns.map(column => column.name)].forEach(label => {{
            const th = document.createElement('th');
            th.textContent = label;
            headerRow.appendChild(th);
        }});
        thead.appendChild(headerRow);

        displayedNodeIndices.slice(0, 250).forEach(nodeIndex => {{
            const row = document.createElement('tr');
            const idCell = document.createElement('td');
            idCell.textContent = nodeId(nodeIndex);
            row.appendChild(idCell);

            const metadata = state.metadataByNodeIndex[nodeIndex] || [];
            metadata.forEach(value => {{
                const cell = document.createElement('td');
                cell.textContent = formatValue(value);
                row.appendChild(cell);
            }});
            tbody.appendChild(row);
        }});

        document.getElementById('selection-summary').textContent = selected.length.toLocaleString() + ' nodes selected';
        document.getElementById('selection-note').textContent = selected.length === 0
            ? (
                displayedNodeIndices.length > 250
                    ? 'No clusters selected. Showing the first 250 rows from the full network table.'
                    : 'No clusters selected. Showing the full network table.'
            )
            : (
                selected.length > 250
                    ? 'Showing the first 250 selected rows in the table. Export includes the full selection.'
                    : 'Click a cluster to toggle it, or Shift-drag a box to add multiple clusters.'
            );
        document.getElementById('export-selected').disabled = !state.bundle;
        document.getElementById('clear-selection').disabled = selected.length === 0;
    }}

    function exportSelection() {{
        const selected = Array.from(state.selectedNodeIndices).sort((left, right) => left - right);
        if (!state.bundle) {{
            return;
        }}
        const exportedNodeIndices = selected.length > 0
            ? selected
            : state.bundle.graph.nodes.map((_, nodeIndex) => nodeIndex);
        const header = ['node_id', ...state.metadataColumns.map(column => column.name)];
        const lines = [header.join('\\t')];
        exportedNodeIndices.forEach(nodeIndex => {{
            const metadata = state.metadataByNodeIndex[nodeIndex] || [];
            const row = [nodeId(nodeIndex), ...metadata.map(value => value === null || value === undefined ? '' : String(value))];
            lines.push(row.join('\\t'));
        }});
        const blob = new Blob([lines.join('\\n')], {{type: 'text/tab-separated-values'}});
        const link = document.createElement('a');
        link.href = URL.createObjectURL(blob);
        link.download = (
            selected.length > 0
                ? ((state.bundle?.name || 'selection') + '_selected.tsv')
                : ((state.bundle?.name || 'network') + '_table.tsv')
        ).replace(/\\s+/g, '_');
        link.click();
        URL.revokeObjectURL(link.href);
    }}

    function updateThresholdUI(resetView = true) {{
        if (!state.bundle) {{
            return;
        }}
        const stop = currentSliderStop();
        if (!stop) {{
            return;
        }}
        document.getElementById('threshold-label').textContent = stop.threshold_label;
        document.getElementById('threshold-input').value = stop.threshold_value === null ? '' : String(stop.threshold_value);
        drawSplitChart();
        drawClusterView(resetView);
        updateMetadataTable();
    }}

    function jumpToThresholdValue() {{
        if (!state.bundle) {{
            return;
        }}
        const input = document.getElementById('threshold-input');
        const targetValue = Number(input.value);
        if (!Number.isFinite(targetValue)) {{
            const stop = currentSliderStop();
            input.value = stop && stop.threshold_value !== null ? String(stop.threshold_value) : '';
            return;
        }}
        const stop = nearestStopForThreshold(targetValue);
        if (!stop) {{
            return;
        }}
        snapSliderToStop(stop);
        updateThresholdUI();
    }}

    function toggleSelectionForComponent(componentId) {{
        const members = componentMembers(componentId);
        const allSelected = members.every(nodeIndex => state.selectedNodeIndices.has(nodeIndex));
        members.forEach(nodeIndex => {{
            if (allSelected) {{
                state.selectedNodeIndices.delete(nodeIndex);
            }} else {{
                state.selectedNodeIndices.add(nodeIndex);
            }}
        }});
        renderClusterView();
        updateMetadataTable();
    }}

    function hitTestComponentAt(screenX, screenY) {{
        const worldPoint = screenToWorldPoint(screenX, screenY);
        for (let index = state.visibleLayout.length - 1; index >= 0; index--) {{
            const item = state.visibleLayout[index];
            const dx = worldPoint.x - item.x;
            const dy = worldPoint.y - item.y;
            if ((dx * dx) + (dy * dy) <= item.radius * item.radius) {{
                return item;
            }}
        }}
        return null;
    }}

    function normalizedSelectionBox(box) {{
        return {{
            left: Math.min(box.startX, box.endX),
            top: Math.min(box.startY, box.endY),
            width: Math.abs(box.endX - box.startX),
            height: Math.abs(box.endY - box.startY),
        }};
    }}

    function circleIntersectsRect(cx, cy, radius, rect) {{
        const nearestX = Math.max(rect.left, Math.min(cx, rect.left + rect.width));
        const nearestY = Math.max(rect.top, Math.min(cy, rect.top + rect.height));
        const dx = cx - nearestX;
        const dy = cy - nearestY;
        return (dx * dx) + (dy * dy) <= radius * radius;
    }}

    function selectComponentsInBox(rect) {{
        let changed = false;
        state.visibleLayout.forEach(item => {{
            const screenPoint = worldToScreenPoint(item.x, item.y);
            const screenRadius = item.radius * state.viewTransform.scale;
            if (!circleIntersectsRect(screenPoint.x, screenPoint.y, screenRadius, rect)) {{
                return;
            }}
            componentMembers(item.componentId).forEach(nodeIndex => {{
                if (!state.selectedNodeIndices.has(nodeIndex)) {{
                    state.selectedNodeIndices.add(nodeIndex);
                    changed = true;
                }}
            }});
        }});
        if (changed) {{
            renderClusterView();
            updateMetadataTable();
        }}
    }}

    clusterCanvas.addEventListener('click', event => {{
        if (!state.bundle) {{
            return;
        }}
        if (state.suppressClick) {{
            state.suppressClick = false;
            return;
        }}
        const point = canvasCoordinatesFromEvent(event);
        const hit = hitTestComponentAt(point.x, point.y);
        if (hit) {{
            toggleSelectionForComponent(hit.componentId);
        }}
    }});

    clusterCanvas.addEventListener('wheel', event => {{
        if (!state.bundle) {{
            return;
        }}
        event.preventDefault();
        const point = canvasCoordinatesFromEvent(event);
        const anchor = screenToWorldPoint(point.x, point.y);
        const zoomFactor = Math.exp(-event.deltaY * 0.0012);
        const nextScale = Math.max(state.viewTransform.minScale, Math.min(state.viewTransform.maxScale, state.viewTransform.scale * zoomFactor));
        state.viewTransform.scale = nextScale;
        state.viewTransform.offsetX = point.x - (anchor.x * nextScale);
        state.viewTransform.offsetY = point.y - (anchor.y * nextScale);
        renderClusterView();
    }}, {{passive: false}});

    clusterCanvas.addEventListener('mousedown', event => {{
        if (!state.bundle) {{
            return;
        }}
        const point = canvasCoordinatesFromEvent(event);
        if (event.shiftKey) {{
            state.dragState = {{mode: 'select', startX: point.x, startY: point.y, endX: point.x, endY: point.y, moved: false}};
            state.selectionBox = normalizedSelectionBox(state.dragState);
        }} else {{
            state.dragState = {{mode: 'pan', startX: point.x, startY: point.y, originOffsetX: state.viewTransform.offsetX, originOffsetY: state.viewTransform.offsetY, moved: false}};
        }}
        clusterCanvas.style.cursor = state.dragState.mode === 'select' ? 'crosshair' : 'grabbing';
    }});

    window.addEventListener('mousemove', event => {{
        if (!state.dragState) {{
            return;
        }}
        const point = canvasCoordinatesFromEvent(event);
        if (state.dragState.mode === 'pan') {{
            const dx = point.x - state.dragState.startX;
            const dy = point.y - state.dragState.startY;
            state.dragState.moved = state.dragState.moved || Math.abs(dx) > 3 || Math.abs(dy) > 3;
            state.viewTransform.offsetX = state.dragState.originOffsetX + dx;
            state.viewTransform.offsetY = state.dragState.originOffsetY + dy;
            renderClusterView();
            return;
        }}
        state.dragState.endX = point.x;
        state.dragState.endY = point.y;
        state.dragState.moved = state.dragState.moved || Math.abs(point.x - state.dragState.startX) > 3 || Math.abs(point.y - state.dragState.startY) > 3;
        state.selectionBox = normalizedSelectionBox(state.dragState);
        renderClusterView();
    }});

    window.addEventListener('mouseup', () => {{
        if (!state.dragState) {{
            return;
        }}
        state.suppressClick = Boolean(state.dragState.moved);
        if (state.dragState.mode === 'select' && state.selectionBox && (state.selectionBox.width > 4 || state.selectionBox.height > 4)) {{
            selectComponentsInBox(state.selectionBox);
        }}
        state.dragState = null;
        state.selectionBox = null;
        clusterCanvas.style.cursor = 'grab';
        renderClusterView();
    }});

    clusterCanvas.addEventListener('mousemove', event => {{
        if (state.dragState) {{
            return;
        }}
        clusterCanvas.style.cursor = event.shiftKey ? 'crosshair' : 'grab';
    }});

    clusterCanvas.addEventListener('mouseleave', () => {{
        if (!state.dragState) {{
            clusterCanvas.style.cursor = 'grab';
        }}
    }});

    document.getElementById('bundle-file').addEventListener('change', async event => {{
        const file = event.target.files?.[0];
        if (!file) {{
            return;
        }}
        try {{
            const bundle = await decodeBundleFile(file);
            installBundle(bundle);
        }} catch (error) {{
            console.error(error);
            setStatus('Failed to load bundle: ' + error.message);
        }}
    }});
    document.getElementById('threshold-slider').addEventListener('input', () => {{
        updateThresholdUI();
    }});
    document.getElementById('threshold-slider').addEventListener('change', () => {{
        const stop = currentSliderStop();
        snapSliderToStop(stop);
        updateThresholdUI();
    }});
    document.getElementById('threshold-input').addEventListener('change', jumpToThresholdValue);
    document.getElementById('threshold-input').addEventListener('keydown', event => {{
        if (event.key !== 'Enter') {{
            return;
        }}
        event.preventDefault();
        jumpToThresholdValue();
    }});
    document.getElementById('min-cluster-size').addEventListener('input', updateThresholdUI);
    document.getElementById('layout-algorithm').addEventListener('change', () => {{
        if (!state.bundle) {{
            return;
        }}
        drawClusterView(true);
    }});
    document.getElementById('leaf-pruning-only').addEventListener('change', updateThresholdUI);
    document.getElementById('color-by').addEventListener('change', renderClusterView);
    document.getElementById('label-by').addEventListener('change', renderClusterView);
    document.getElementById('show-labels').addEventListener('change', renderClusterView);
    document.getElementById('show-node-counts').addEventListener('change', renderClusterView);
    document.getElementById('show-edge-scores').addEventListener('change', renderClusterView);
    document.getElementById('export-selected').addEventListener('click', exportSelection);
    document.getElementById('reset-view').addEventListener('click', () => {{
        fitClusterViewToLayout();
        renderClusterView();
    }});
    document.getElementById('clear-selection').addEventListener('click', () => {{
        state.selectedNodeIndices = new Set();
        renderClusterView();
        updateMetadataTable();
    }});
    window.addEventListener('resize', () => {{
        if (!state.bundle) {{
            return;
        }}
        drawClusterView(true);
    }});

    warnIfUnsupported();
    drawSplitChart();
    drawClusterView();
</script>
</body>
</html>
"""


def write_ssn_viewer_html(out_path: str, title: str = "Domainator SSN Viewer") -> None:
    temp_path = make_temporary_output_path(out_path)
    try:
        with open(temp_path, "w", encoding="utf-8") as out_handle:
            out_handle.write(ssn_viewer_html(title=title))
        Path(temp_path).replace(out_path)
        temp_path = None
    finally:
        if temp_path is not None and Path(temp_path).exists():
            Path(temp_path).unlink()