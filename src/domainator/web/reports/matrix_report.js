// The matrix report (matrix_report): a similarity matrix's size and scores, and how its
// network falls apart into clusters as the threshold rises, to help choose a threshold for
// build_ssn.
//
// report-data's data is {summary, mst_edges, threshold_table, merge_moving_sum, slider_stops,
// knn_projection}:
// - summary is the document --json_out writes, split_events being the split events to plot.
// - mst_edges is the maximum spanning tree, [node, node, score] strongest first.
// - threshold_table has one row per distinct threshold, as the columns threshold, edges and
//   clusters: what build_ssn --threshold <threshold> keeps, edges scoring strictly above it.
// - merge_moving_sum is {window, x, y}, the split impact within window/2 of each x.
// - slider_stops are the thresholds the slider offers: infinity (threshold_value null), one
//   per split event, then one below every edge. Each has the edge_index of the last MST edge it
//   keeps and the threshold_index of its threshold_table row (-1 for none).
// - knn_projection is null, or {min_k, max_k, default_k, counts}, counts[i][k - min_k] being
//   the edges build_ssn --threshold <threshold_table.threshold[i]> --mst --knn k keeps.

// @include common/payload.js
// @include common/panels.js
// @include common/download.js
// @include common/table.js
// @include common/charts.js
// @include common/line_chart.js
// @include common/pack.js

// What a split's size means under each --merge_impact_metric; ssn_hierarchy.py's
// merge_impact_axis_labels, which test_report_html.py checks these against.
const SPLIT_AXIS_LABELS = {
    min_child: {
        largest: "Largest single split (nodes)",
        moving_sum: "Moving sum of split size (5% window)",
        moving_sum_short: "Moving sum (5% window)",
        impact_amount: "{} nodes",
        moving_sum_hover: "Nodes displaced within",
    },
    product: {
        largest: "Largest single split (size product)",
        moving_sum: "Moving sum of split impact (5% window)",
        moving_sum_short: "Moving sum (5% window)",
        impact_amount: "impact {}",
        moving_sum_hover: "Split impact within",
    },
};
// ssn_hierarchy.MERGE_EVENT_DENSITY_BINS: the bands of the axis --max_merge_events backfills.
const MERGE_EVENT_DENSITY_BINS = 20;
// The slider's positions run 0 to SLIDER_SCALE.
const SLIDER_SCALE = 10000;
// How much of the threshold axis is left empty past the outermost stops, on each side.
const AXIS_PADDING_FRACTION = 0.05;
// The split chart, whose plot the slider underneath it lines up with.
const SPLIT_CHART = {width: 1000, height: 340, margin: {top: 40, right: 72, bottom: 48, left: 72}};
const SPLIT_STEM_COLOR = "#a8d2cf";
const SPLIT_BEAD_COLOR = "#72b7b2";
const MOVING_SUM_COLOR = "#e45756";
// A cluster-size union-find snapshot is kept every this many slider stops.
const CHECKPOINT_STRIDE = 50;

function summaryPanel(root, summary) {
    const section = reportPanel(root, "Summary");
    const scores = summary.edge_scores;
    const columns = [[
        ["nodes", formatChartNumber(summary.nodes)],
        ["non-zero edges", formatChartNumber(summary.non_zero_edges)],
        ["connected components", formatChartNumber(summary.connected_components)],
    ]];
    if (scores) {
        columns.push([
            ["mean score", formatFixed1(scores.mean)],
            ["median score", formatFixed1(scores.median)],
            ["lowest score", formatFixed1(scores.min)],
            ["highest score", formatFixed1(scores.max)],
        ]);
    }
    section.appendChild(statGrid(columns));
}

// Clusters, and the edges kept, by threshold; and clusters by edges kept.
function thresholdCharts(root, table) {
    const section = reportPanel(root, "Clusters by Threshold",
        "Drag across a chart to zoom in, and double-click it to zoom back out.");
    const grid = panelPart(section, null, "chart-grid");
    // Ascending threshold; the table is descending.
    const rows = table.threshold.map((threshold, i) => i).reverse();
    // A log axis needs positive thresholds, and earns its place over a decade or more.
    const positive = rows.filter(i => table.threshold[i] > 0);
    const xLog = positive.length >= 2 && table.threshold[positive[positive.length - 1]] / table.threshold[positive[0]] >= 10;
    const plotted = xLog ? positive : rows;
    const byThreshold = drawLineChart(panelPart(grid, "clusters-by-threshold"), {
        series: [
            {name: "Clusters", color: "#1f77b4", points: plotted.map(i => [table.threshold[i], table.clusters[i]])},
            {name: "Edges kept", color: "#2ca02c", axis: "right", points: plotted.map(i => [table.threshold[i], table.edges[i]])},
        ],
        xLabel: xLog ? "Threshold (log scale)" : "Threshold", xName: "Threshold", xLog,
        xFormat: value => value.toLocaleString("en-US", {maximumFractionDigits: 2}),
        yLabel: "Clusters", y2Label: "Edges kept (build_ssn --threshold)",
        fileStem: "matrix_report_clusters_by_threshold", idPrefix: "threshold-chart-",
    });
    // Edges kept rises as the threshold falls, so this is the table read bottom up.
    const edgeRows = rows.slice().reverse();
    const byEdges = drawLineChart(panelPart(grid, "clusters-by-edges"), {
        series: [{name: "Clusters", color: "#ff7f0e", points: edgeRows.map(i => [table.edges[i], table.clusters[i]])}],
        xLabel: "Edges kept", yLabel: "Clusters",
        fileStem: "matrix_report_clusters_by_edges", idPrefix: "edges-chart-",
    });
    return {byThreshold, byEdges};
}

// The threshold axis the split chart and the slider share: the stops' and the moving sum's
// range, padded, or null when there is no finite stop.
function thresholdAxis(stops, movingSum) {
    const values = stops.filter(stop => stop.threshold_value !== null).map(stop => stop.threshold_value);
    if (movingSum && movingSum.x.length > 0) {
        values.push(movingSum.x[0], movingSum.x[movingSum.x.length - 1]);
    }
    if (values.length === 0) {
        return null;
    }
    let lo = values[0];
    let hi = values[0];
    for (const value of values) {
        lo = Math.min(lo, value);
        hi = Math.max(hi, value);
    }
    const pad = hi > lo ? (hi - lo) * AXIS_PADDING_FRACTION : Math.max(Math.abs(hi) * AXIS_PADDING_FRACTION, 1);
    return [lo - pad, hi + pad];
}

// Each stop's place on the slider. The slider runs right to left, so position 0, infinity, is
// at its right end, and the others sit under their threshold on the axis [lo, hi] above. Every
// stop is kept one position past the last, so none is out of reach.
function sliderPositions(stops, [lo, hi]) {
    const positions = [];
    let previous = -1;
    for (const stop of stops) {
        let position = stop.threshold_value === null
            ? 0 : Math.round((hi - stop.threshold_value) / (hi - lo) * SLIDER_SCALE);
        position = Math.min(SLIDER_SCALE, Math.max(previous + 1, position));
        positions.push(position);
        previous = position;
    }
    return positions;
}

function fillTemplate(template, value) {
    return template.replace("{}", value);
}

// Stems at each split event's threshold as tall as its largest split, a bead for each split
// size, and the moving sum of split impact on the right axis. A tall stem with one bead is one
// large clean split; short stems under a high line are a network crumbling into small pieces.
function drawSplitChart(container, {events, movingSum, labels, axis, onPick}) {
    const {width, height, margin} = SPLIT_CHART;
    const plot = {left: margin.left, right: width - margin.right, top: margin.top, bottom: height - margin.bottom};
    const beads = [];
    let tallest = 1;
    for (const event of events) {
        tallest = Math.max(tallest, event.largest_merge);
        for (const [size, count] of Object.entries(event.merge_size_counts || {})) {
            // JSON object keys are strings.
            beads.push({event, size: Number(size), count});
        }
    }
    const sumPoints = movingSum ? movingSum.x.map((x, i) => [x, movingSum.y[i]]) : [];
    let highestSum = 1;
    for (const [, y] of sumPoints) {
        highestSum = Math.max(highestSum, y);
    }
    const sumWindow = movingSum ? movingSum.window : 0;
    const formatThreshold = value => value.toLocaleString("en-US", {minimumFractionDigits: 2, maximumFractionDigits: 2});

    const draw = (svg, [lo, hi]) => {
        const x = axisScale(lo, hi, plot.left, plot.right);
        const [yTop, yTicks] = zeroBasedAxis(tallest, 5, true);
        const y = axisScale(0, yTop, plot.bottom, plot.top);
        const [y2Top, y2Ticks] = zeroBasedAxis(highestSum, 5);
        const y2 = axisScale(0, y2Top, plot.bottom, plot.top);
        drawAxes(svg, {x, y, xTicks: niceTicks(lo, hi, 10), yTicks, xLabel: "Threshold", yLabel: labels.largest, ...plot});
        drawRightAxis(svg, {y: y2, ticks: y2Ticks, right: plot.right, top: plot.top, bottom: plot.bottom, label: labels.moving_sum});
        drawLegend(svg, [
            {label: "Individual cluster split", color: SPLIT_BEAD_COLOR, mark: "dot"},
            {label: labels.moving_sum_short, color: MOVING_SUM_COLOR, mark: "line"},
        ], {left: plot.left, top: 14});
        const stems = svgElement("g", {stroke: SPLIT_STEM_COLOR, "stroke-width": 1.6}, svg);
        for (const event of events) {
            const sx = x(event.threshold_value);
            svgElement("line", {x1: sx, x2: sx, y1: y(0), y2: y(event.largest_merge)}, stems);
        }
        // No outline: on a stem with several close beads, outlines would merge into a band
        // that hides the stem.
        const dots = svgElement("g", {fill: SPLIT_BEAD_COLOR}, svg);
        for (const bead of beads) {
            svgElement("circle", {cx: x(bead.event.threshold_value), cy: y(bead.size), r: 2.5}, dots);
        }
        svgElement("path", {
            d: linePath(sumPoints, x, y2, true), fill: "none", stroke: MOVING_SUM_COLOR, "stroke-width": 2,
        }, svg);

        // The bead nearest the pointer within 10 pixels; else the stem within 14; else the
        // moving sum under it.
        const hover = (px, py) => {
            let best = null;
            let bestDistance = 10;
            for (const bead of beads) {
                const distance = Math.hypot(x(bead.event.threshold_value) - px, y(bead.size) - py);
                if (distance < bestDistance) {
                    best = bead;
                    bestDistance = distance;
                }
            }
            if (best) {
                const event = best.event;
                return {
                    text: [
                        `Threshold: ${formatThreshold(event.threshold_value)}`,
                        `${formatChartNumber(best.count)} split(s) of ${fillTemplate(labels.impact_amount, formatChartNumber(best.size))}`,
                        `${formatChartNumber(event.merge_count)} split(s) here, ${fillTemplate(labels.impact_amount, formatChartNumber(event.merge_impact))} total`,
                        `From: ${event.threshold_from}`,
                        `To: ${event.threshold_to}`,
                    ].join("\n"),
                    dots: [{x: x(event.threshold_value), y: y(best.size), color: SPLIT_BEAD_COLOR}],
                };
            }
            let stem = null;
            let stemDistance = 14;
            for (const event of events) {
                const distance = Math.abs(x(event.threshold_value) - px);
                if (distance < stemDistance && py >= y(event.largest_merge) - 10) {
                    stem = event;
                    stemDistance = distance;
                }
            }
            if (stem) {
                return {
                    text: [
                        `Threshold: ${formatThreshold(stem.threshold_value)}`,
                        `Largest single split: ${fillTemplate(labels.impact_amount, formatChartNumber(stem.largest_merge))}`,
                        `${formatChartNumber(stem.merge_count)} split(s) here, ${fillTemplate(labels.impact_amount, formatChartNumber(stem.merge_impact))} total`,
                        `From: ${stem.threshold_from}`,
                        `To: ${stem.threshold_to}`,
                    ].join("\n"),
                    guide: x(stem.threshold_value),
                    dots: [],
                };
            }
            if (sumPoints.length === 0) {
                return null;
            }
            const value = x.invert(px);
            // The step line holds each value until the next x.
            let i = nearestIndex(sumPoints, value);
            if (sumPoints[i][0] > value && i > 0) {
                i -= 1;
            }
            return {
                text: `Threshold: ${formatThreshold(value)}\n${labels.moving_sum_hover} ±${formatThreshold(sumWindow / 2)}: ${formatChartNumber(sumPoints[i][1])}`,
                guide: px,
                dots: [{x: px, y: y2(sumPoints[i][1]), color: MOVING_SUM_COLOR}],
            };
        };
        return {plot, x, hover, click: px => onPick(x.invert(px))};
    };
    return interactiveChart(container, {
        width, height, domain: axis, draw, zoom: false,
        ariaLabel: "Cluster splits by threshold", fileStem: "split_events", idPrefix: "split-chart-",
    });
}

// How much of the split series is drawn: the x axis and the moving sum span every split event,
// but --max_merge_events caps the stems, and without this note a stretch of axis with no stems
// reads as a quiet one rather than a filtered one.
function splitEventNote(plotted, total, cap) {
    const plural = count => (count === 1 ? "" : "s");
    if (plotted >= total) {
        return {text: `All ${formatChartNumber(total)} split event${plural(total)} plotted. Click the chart to move the threshold there.`};
    }
    const backfilled = Math.max(0, plotted - cap);
    const how = backfilled > 0
        ? `the strongest ${formatChartNumber(cap)} by impact, plus ${formatChartNumber(backfilled)} so that every ${100 / MERGE_EVENT_DENSITY_BINS}% of the axis with an event to show has one`
        : "the strongest by impact";
    return {
        text: `${formatChartNumber(plotted)} of ${formatChartNumber(total)} split events plotted: ${how}. Click the chart to move the threshold there.`,
        title: "Set by domainator matrix_report --max_merge_events (0 plots them all). A stretch with no stems "
            + "is a stretch whose events were too small to make the cut, not necessarily a quiet one.",
    };
}

// Cluster sizes as the MST's edges are added strongest first, with a union-find snapshot kept
// every CHECKPOINT_STRIDE stops so that moving the slider never replays the whole tree.
class ClusterSizes {
    constructor(nodes, mstEdges, stops) {
        this.nodes = nodes;
        this.edges = mstEdges;
        this.checkpoints = [{edgeIndex: -1, parent: Array.from({length: nodes}, (_, i) => i), size: new Array(nodes).fill(1)}];
        const parent = this.checkpoints[0].parent.slice();
        const size = this.checkpoints[0].size.slice();
        let next = 0;
        stops.forEach((stop, i) => {
            while (next <= stop.edge_index && next < mstEdges.length) {
                ClusterSizes.union(parent, size, mstEdges[next][0], mstEdges[next][1]);
                next += 1;
            }
            if (i > 0 && (i % CHECKPOINT_STRIDE === 0 || i === stops.length - 1)) {
                this.checkpoints.push({edgeIndex: next - 1, parent: parent.slice(), size: size.slice()});
            }
        });
    }

    static find(parent, node) {
        let root = node;
        while (parent[root] !== root) {
            root = parent[root];
        }
        while (parent[node] !== root) {
            const up = parent[node];
            parent[node] = root;
            node = up;
        }
        return root;
    }

    static union(parent, size, a, b) {
        let ra = ClusterSizes.find(parent, a);
        let rb = ClusterSizes.find(parent, b);
        if (ra === rb) {
            return;
        }
        if (size[ra] > size[rb]) {
            [ra, rb] = [rb, ra];
        }
        parent[ra] = rb;
        size[rb] += size[ra];
    }

    // The cluster sizes once MST edges 0 to edgeIndex are in, largest first.
    at(edgeIndex) {
        let checkpoint = this.checkpoints[0];
        for (const candidate of this.checkpoints) {
            if (candidate.edgeIndex <= edgeIndex) {
                checkpoint = candidate;
            }
        }
        const parent = checkpoint.parent.slice();
        const size = checkpoint.size.slice();
        for (let i = checkpoint.edgeIndex + 1; i <= edgeIndex && i < this.edges.length; i += 1) {
            ClusterSizes.union(parent, size, this.edges[i][0], this.edges[i][1]);
        }
        const sizes = [];
        for (let node = 0; node < this.nodes; node += 1) {
            if (parent[node] === node) {
                sizes.push(size[node]);
            }
        }
        return sizes.sort((a, b) => b - a);
    }
}

// The threshold slider, what the network looks like at its threshold, and (with
// --knn_projection) the k slider.
function thresholdPanel(root, data, charts) {
    const summary = data.summary;
    const stops = data.slider_stops;
    const events = summary.split_events;
    const section = reportPanel(root, "Cluster Splits by Threshold");
    const axis = thresholdAxis(stops, data.merge_moving_sum);
    let selectStop = () => {};
    const pickThreshold = value => {
        const position = Math.round((axis[1] - value) / (axis[1] - axis[0]) * SLIDER_SCALE);
        selectStop(nearestIndex(positions, position, p => p));
    };
    let splitChart = null;
    if (events.length > 0 && axis) {
        const labels = SPLIT_AXIS_LABELS[summary.merge_impact_metric] || SPLIT_AXIS_LABELS.min_child;
        splitChart = drawSplitChart(panelPart(section, "split-chart"), {
            events, movingSum: data.merge_moving_sum, labels, axis, onPick: pickThreshold,
        });
        const note = splitEventNote(events.length, summary.split_event_total, summary.max_merge_events);
        const text = document.createElement("p");
        text.className = "note";
        text.id = "split-event-count";
        text.textContent = note.text;
        if (note.title) {
            text.title = note.title;
        }
        section.appendChild(text);
    }

    const positions = axis ? sliderPositions(stops, axis) : [0];
    const controls = panelPart(section, null, "slider-controls");
    // As wide as the split chart, which is at most its own width.
    controls.style.maxWidth = `${SPLIT_CHART.width}px`;
    const label = document.createElement("label");
    label.htmlFor = "threshold-slider";
    label.append("Threshold (build_ssn --threshold): ");
    const number = document.createElement("strong");
    number.id = "threshold-number";
    label.appendChild(number);
    const track = panelPart(controls, null, "slider-track");
    // Line the track up with the split chart's plot, allowing for the thumb's half width.
    const {width, margin} = SPLIT_CHART;
    if (splitChart) {
        track.style.paddingLeft = `calc(${(margin.left / width * 100).toFixed(3)}% - 9px)`;
        track.style.paddingRight = `calc(${(margin.right / width * 100).toFixed(3)}% - 9px)`;
    }
    const slider = document.createElement("input");
    slider.type = "range";
    slider.id = "threshold-slider";
    slider.min = "0";
    slider.max = String(stops.length > 1 ? SLIDER_SCALE : 0);
    slider.step = "1";
    slider.value = "0";
    slider.disabled = stops.length <= 1;
    track.appendChild(slider);
    controls.insertBefore(label, track);

    let kSlider = null;
    const knn = data.knn_projection;
    if (knn) {
        const kLabel = document.createElement("label");
        kLabel.htmlFor = "mst-knn-k-slider";
        kLabel.append("MST + Knn k (build_ssn --mst --knn): ");
        const kNumber = document.createElement("strong");
        kNumber.id = "mst-knn-k-number";
        kLabel.appendChild(kNumber);
        kSlider = document.createElement("input");
        kSlider.type = "range";
        kSlider.id = "mst-knn-k-slider";
        kSlider.min = String(knn.min_k);
        kSlider.max = String(knn.max_k);
        kSlider.step = "1";
        kSlider.value = String(knn.default_k);
        controls.append(kLabel, panelPart(controls, null, "slider-track"));
        controls.lastChild.appendChild(kSlider);
    }

    const columns = [
        [
            ["clusters", formatChartNumber(summary.nodes), "num-clusters"],
            ["largest cluster", "1", "largest-cluster"],
            ["mean cluster size", "1.0", "avg-cluster-size"],
        ],
        [
            ["singletons", formatChartNumber(summary.nodes), "num-singletons"],
            ["edges kept", "0", "num-edges"],
        ],
    ];
    if (knn) {
        columns[1].push(["projected MST + Knn edges", "0", "mst-knn-edges"]);
    }
    section.appendChild(statGrid(columns));

    const heading = document.createElement("h2");
    heading.textContent = "100 Largest Clusters";
    section.appendChild(heading);
    const bubbles = createBubbleChart(panelPart(section, "cluster-bubbles"), {
        itemLabel: "nodes", label: i => `Cluster ${i + 1}`, emptyText: "No clusters to show.",
        fileStem: "matrix_report_largest_clusters", idPrefix: "bubbles-",
    });

    const clusters = new ClusterSizes(summary.nodes, data.mst_edges, stops);
    const table = data.threshold_table;
    let current = 0;
    const show = () => {
        const stop = stops[current];
        const sizes = clusters.at(stop.edge_index);
        const row = stop.threshold_index;
        const edges = row >= 0 && row < table.edges.length ? table.edges[row] : 0;
        number.textContent = stop.threshold_label;
        document.getElementById("num-clusters").textContent = formatChartNumber(sizes.length);
        document.getElementById("num-singletons").textContent = formatChartNumber(sizes.filter(size => size === 1).length);
        document.getElementById("largest-cluster").textContent = formatChartNumber(sizes[0] || 0);
        document.getElementById("num-edges").textContent = formatChartNumber(edges);
        const mean = sizes.length > 0 ? sizes.reduce((a, b) => a + b, 0) / sizes.length : 0;
        document.getElementById("avg-cluster-size").textContent = mean.toFixed(1);
        if (knn) {
            const k = Number(kSlider.value);
            document.getElementById("mst-knn-k-number").textContent = String(k);
            const counts = row >= 0 && row < knn.counts.length ? knn.counts[row][k - knn.min_k] : 0;
            document.getElementById("mst-knn-edges").textContent = formatChartNumber(counts);
        }
        bubbles.update(sizes.slice(0, 100));
        const threshold = stop.threshold_value;
        if (splitChart) {
            splitChart.setMarker(threshold);
        }
        charts.byThreshold.setMarker(row >= 0 ? table.threshold[row] : null);
        charts.byEdges.setMarker(row >= 0 ? edges : null);
    };
    selectStop = index => {
        current = Math.max(0, Math.min(stops.length - 1, index));
        slider.value = String(positions[current]);
        show();
    };
    // The nearest stop to the thumb while it moves, snapping the thumb to it when let go.
    let pending = null;
    slider.addEventListener("input", () => {
        cancelAnimationFrame(pending);
        pending = requestAnimationFrame(() => {
            current = nearestIndex(positions, Number(slider.value), p => p);
            show();
        });
    });
    slider.addEventListener("change", () => selectStop(nearestIndex(positions, Number(slider.value), p => p)));
    // The arrow keys step from stop to stop; the slider's own step would never leave one.
    slider.addEventListener("keydown", event => {
        const moves = {ArrowLeft: 1, ArrowDown: 1, ArrowRight: -1, ArrowUp: -1, PageDown: 10, PageUp: -10};
        if (event.key in moves) {
            selectStop(current + moves[event.key]);
        } else if (event.key === "Home" || event.key === "End") {
            selectStop(event.key === "Home" ? 0 : stops.length - 1);
        } else {
            return;
        }
        event.preventDefault();
    });
    if (kSlider) {
        kSlider.addEventListener("input", show);
    }
    // Start at the lowest threshold, the whole network, at the slider's left end.
    selectStop(stops.length - 1);
    return {select: selectStop, positions};
}

(async () => {
    const root = document.getElementById("report-root");
    try {
        const {data} = await readReportData();
        root.replaceChildren();
        summaryPanel(root, data.summary);
        const charts = thresholdCharts(root, data.threshold_table);
        window.matrixReport = {charts, threshold: thresholdPanel(root, data, charts)};
    } catch (error) {
        showReportError(root, error);
    }
})();
