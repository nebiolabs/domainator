// The projection plot (build_projection): each row of the matrix as a point, colored by a
// metadata column.
//
// report-data's data is {algorithm, axes, ids, coords, color_by, groups, group, metadata}:
// - axes are the coordinate columns, ["projection_X", "projection_Y"] or with "projection_Z";
//   coords holds one array per axis, coords[a][i] being point i's coordinate on axes[a].
// - ids are the points' names.
// - color_by is the metadata column the points are colored by, or null. groups are then the
//   legend's entries {value, color, count} (value null for points with no value), and group
//   each point's index in groups; both are null without color_by.
// - metadata is {columns, values}, values[c][i] being point i's value in columns[c] (null when
//   it has none).
//
// The frame (axes, legend) is an SVG and the points are drawn on a canvas over its plot area, so
// that large projections stay quick to zoom. The downloads draw the points into the SVG, or
// onto the PNG, the same way.

// @include common/payload.js
// @include common/download.js
// @include common/charts.js

// The plot area, in the SVG's units.
const PLOT = {left: 64, top: 16, width: 600, height: 600, bottom: 48, right: 16};
const LEGEND_GAP = 24;
const LEGEND_ROW = 18;
const LEGEND_SWATCH = 12;
// Longest legend label shown in full; the entry's title holds all of it.
const LEGEND_LABEL_LENGTH = 40;
// Longest metadata value a tooltip shows in full.
const TOOLTIP_VALUE_LENGTH = 120;
const DEFAULT_POINT_COLOR = "#1F77B4";
const POINT_OPACITY = 0.75;
const HIGHLIGHT_COLOR = "#000000";
// A press that moves less than this many pixels is a click, not a drag.
const CLICK_SLOP = 3;
// How far, in screen pixels, past a point's edge the pointer still hovers it.
const HOVER_SLOP = 3;
const AXIS_NAMES = ["X", "Y", "Z"];

function defaultPointRadius(count) {
    if (count <= 500) {
        return 4;
    }
    if (count <= 5000) {
        return 3;
    }
    return count <= 50000 ? 2 : 1.5;
}

function shorten(text, length) {
    return text.length > length ? `${text.slice(0, length - 1)}…` : text;
}

function groupLabel(group) {
    return `${group.value === null ? "no value" : group.value} (${formatChartNumber(group.count)})`;
}

// The plot's state: which axes are shown, the view (the data point at the plot's center and
// the data units per SVG unit), the groups hidden from the legend and the highlighted point.
function projectionState(data) {
    const count = data.ids.length;
    const groups = data.groups || [{value: null, color: DEFAULT_POINT_COLOR, count}];
    const group = data.group || new Array(count).fill(0);
    const groupPoints = groups.map(() => []);
    group.forEach((g, i) => groupPoints[g].push(i));
    return {
        data, count, groups, groupPoints,
        pair: [0, 1],
        view: null,
        radius: defaultPointRadius(count),
        hidden: new Set(),
        highlight: null,
    };
}

// The view that shows every point, with the same units on both axes.
function fittedView(state) {
    const [xs, ys] = state.pair.map(axis => state.data.coords[axis]);
    let [x0, x1, y0, y1] = [Infinity, -Infinity, Infinity, -Infinity];
    for (let i = 0; i < state.count; i += 1) {
        x0 = Math.min(x0, xs[i]);
        x1 = Math.max(x1, xs[i]);
        y0 = Math.min(y0, ys[i]);
        y1 = Math.max(y1, ys[i]);
    }
    if (!(x1 >= x0)) {
        return {cx: 0, cy: 0, u: 1};
    }
    const u = Math.max((x1 - x0) / PLOT.width, (y1 - y0) / PLOT.height) * 1.1 || 1;
    return {cx: (x0 + x1) / 2, cy: (y0 + y1) / 2, u};
}

// The data range the view shows, and functions from data to plot-area coordinates.
function plotScales(view) {
    const x0 = view.cx - view.u * PLOT.width / 2;
    const y1 = view.cy + view.u * PLOT.height / 2;
    return {
        x0, x1: x0 + view.u * PLOT.width, y0: y1 - view.u * PLOT.height, y1,
        sx: x => (x - x0) / view.u,
        sy: y => (y1 - y) / view.u,
    };
}

// Calls visit(i, px, py, g) for each point of a shown group that is in the plot area, group by
// group, px and py being plot-area coordinates.
function forEachShownPoint(state, visit) {
    const {sx, sy} = plotScales(state.view);
    const [xs, ys] = state.pair.map(axis => state.data.coords[axis]);
    const r = state.radius;
    state.groupPoints.forEach((points, g) => {
        if (state.hidden.has(g)) {
            return;
        }
        for (const i of points) {
            const px = sx(xs[i]);
            const py = sy(ys[i]);
            if (px >= -r && px <= PLOT.width + r && py >= -r && py <= PLOT.height + r) {
                visit(i, px, py, g);
            }
        }
    });
}

// Draws the points on ctx, whose units are plot-area units. Returns how many it drew.
function drawPoints(ctx, state) {
    let drawn = 0;
    let current = null;
    ctx.globalAlpha = POINT_OPACITY;
    const r = state.radius;
    forEachShownPoint(state, (i, px, py, g) => {
        if (g !== current) {
            if (current !== null) {
                ctx.fill();
            }
            current = g;
            ctx.fillStyle = state.groups[g].color;
            ctx.beginPath();
        }
        ctx.moveTo(px + r, py);
        ctx.arc(px, py, r, 0, 2 * Math.PI);
        drawn += 1;
    });
    if (current !== null) {
        ctx.fill();
    }
    ctx.globalAlpha = 1;
    if (state.highlight !== null) {
        const {sx, sy} = plotScales(state.view);
        const [xs, ys] = state.pair.map(axis => state.data.coords[axis]);
        ctx.strokeStyle = HIGHLIGHT_COLOR;
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.arc(sx(xs[state.highlight]), sy(ys[state.highlight]), r + 3, 0, 2 * Math.PI);
        ctx.stroke();
    }
    return drawn;
}

function legendWidth(groups) {
    const longest = Math.max(...groups.map(group => shorten(groupLabel(group), LEGEND_LABEL_LENGTH).length));
    return Math.max(120, LEGEND_SWATCH + 8 + longest * 7);
}

// The SVG's size: the plot area, its axes and, with color_by, the legend of `groups`.
function frameSize(state, groups) {
    const plotRight = PLOT.left + PLOT.width;
    const plotBottom = PLOT.top + PLOT.height + PLOT.bottom;
    if (state.data.color_by === null) {
        return {width: plotRight + PLOT.right, height: plotBottom};
    }
    return {
        width: plotRight + LEGEND_GAP + legendWidth(groups) + PLOT.right,
        height: Math.max(plotBottom, PLOT.top + (groups.length + 1) * LEGEND_ROW + 8),
    };
}

function drawPlotAxes(layer, state) {
    const {x0, x1, y0, y1, sx, sy} = plotScales(state.view);
    layer.replaceChildren();
    drawAxes(layer, {
        x: value => PLOT.left + sx(value), y: value => PLOT.top + sy(value),
        xTicks: niceTicks(x0, x1, 8), yTicks: niceTicks(y0, y1, 8),
        left: PLOT.left, right: PLOT.left + PLOT.width, top: PLOT.top, bottom: PLOT.top + PLOT.height,
        xLabel: state.data.axes[state.pair[0]], yLabel: state.data.axes[state.pair[1]],
    });
    svgElement("rect", {
        x: PLOT.left, y: PLOT.top, width: PLOT.width, height: PLOT.height, fill: "none", stroke: CHART_INK,
    }, layer);
}

// The legend's entries for `indices` (indices into state.groups); a hidden group is dimmed.
// Interactive entries can be clicked or focused and toggled with Enter or Space.
function drawLegend(layer, state, indices, interactive) {
    layer.replaceChildren();
    const left = PLOT.left + PLOT.width + LEGEND_GAP;
    const legend = svgElement("g", {"font-family": CHART_FONT, "font-size": 12, fill: CHART_INK}, layer);
    const title = svgElement("text", {x: left, y: PLOT.top + 12, "font-weight": 600}, legend);
    title.textContent = state.data.color_by;
    indices.forEach((g, row) => {
        const group = state.groups[g];
        const top = PLOT.top + (row + 1) * LEGEND_ROW + 4;
        const entry = svgElement("g", {"data-group": g}, legend);
        if (interactive) {
            entry.setAttribute("tabindex", 0);
            entry.setAttribute("role", "button");
            entry.setAttribute("aria-pressed", String(!state.hidden.has(g)));
        }
        if (state.hidden.has(g)) {
            entry.setAttribute("opacity", 0.35);
        }
        const label = groupLabel(group);
        svgElement("title", {}, entry).textContent = label;
        svgElement("rect", {x: left, y: top, width: LEGEND_SWATCH, height: LEGEND_SWATCH, fill: group.color}, entry);
        const text = svgElement("text", {x: left + LEGEND_SWATCH + 6, y: top + LEGEND_SWATCH - 2}, entry);
        text.textContent = shorten(label, LEGEND_LABEL_LENGTH);
    });
}

// A standalone SVG of the frame: white background, axes and the legend of the shown groups.
function exportFrame(state) {
    const shown = state.groups.map((group, g) => g).filter(g => !state.hidden.has(g));
    const {width, height} = frameSize(state, shown.map(g => state.groups[g]));
    const svg = svgElement("svg", {viewBox: `0 0 ${width} ${height}`});
    svgElement("rect", {width, height, fill: "#ffffff"}, svg);
    drawPlotAxes(svgElement("g", {}, svg), state);
    if (state.data.color_by !== null) {
        drawLegend(svgElement("g", {}, svg), state, shown, false);
    }
    return svg;
}

function exportSvg(state) {
    const svg = exportFrame(state);
    const defs = svgElement("defs", {}, svg);
    const clip = svgElement("clipPath", {id: "projection-plot-clip"}, defs);
    svgElement("rect", {width: PLOT.width, height: PLOT.height}, clip);
    const points = svgElement("g", {
        transform: `translate(${PLOT.left} ${PLOT.top})`, "clip-path": "url(#projection-plot-clip)",
        "fill-opacity": POINT_OPACITY,
    }, svg);
    const round = value => Math.round(value * 100) / 100;
    forEachShownPoint(state, (i, px, py, g) => {
        svgElement("circle", {cx: round(px), cy: round(py), r: state.radius, fill: state.groups[g].color}, points);
    });
    if (state.highlight !== null) {
        const {sx, sy} = plotScales(state.view);
        const [xs, ys] = state.pair.map(axis => state.data.coords[axis]);
        svgElement("circle", {
            cx: round(sx(xs[state.highlight])), cy: round(sy(ys[state.highlight])), r: state.radius + 3,
            fill: "none", stroke: HIGHLIGHT_COLOR, "stroke-width": 2, "fill-opacity": 1,
        }, points);
    }
    return svg;
}

// The frame drawn from its SVG at `scale`, then the points drawn over it.
function exportPng(state, fileStem, scale = 2) {
    const frame = exportFrame(state);
    const box = frame.viewBox.baseVal;
    const url = URL.createObjectURL(new Blob([svgMarkup(frame)], {type: "image/svg+xml"}));
    const image = new Image();
    image.onload = () => {
        const canvas = document.createElement("canvas");
        canvas.width = Math.round(box.width * scale);
        canvas.height = Math.round(box.height * scale);
        const context = canvas.getContext("2d");
        context.drawImage(image, 0, 0, canvas.width, canvas.height);
        URL.revokeObjectURL(url);
        context.setTransform(scale, 0, 0, scale, PLOT.left * scale, PLOT.top * scale);
        context.beginPath();
        context.rect(0, 0, PLOT.width, PLOT.height);
        context.clip();
        drawPoints(context, state);
        canvas.toBlob(blob => downloadBlob(`${fileStem}.png`, blob), "image/png");
    };
    image.src = url;
}

function toolbarButton(id, label, onClick) {
    const button = document.createElement("button");
    button.type = "button";
    button.id = id;
    button.textContent = label;
    button.addEventListener("click", onClick);
    return button;
}

function tooltipText(state, i) {
    const {data} = state;
    const lines = [data.ids[i]];
    data.axes.forEach((axis, a) => lines.push(`${axis}: ${formatChartNumber(data.coords[a][i])}`));
    data.metadata.columns.forEach((column, c) => {
        const value = data.metadata.values[c][i];
        if (value !== null) {
            lines.push(`${column}: ${shorten(String(value), TOOLTIP_VALUE_LENGTH)}`);
        }
    });
    return lines.join("\n");
}

function drawProjection(container, data) {
    const state = projectionState(data);
    state.view = fittedView(state);
    const {width, height} = frameSize(state, state.groups);

    const svg = svgElement("svg", {
        viewBox: `0 0 ${width} ${height}`, role: "img",
        "aria-label": `Projection of ${formatChartNumber(state.count)} points`,
    });
    const axesLayer = svgElement("g", {}, svg);
    const legendLayer = svgElement("g", {}, svg);
    const canvas = document.createElement("canvas");
    Object.assign(canvas.style, {
        left: `${PLOT.left / width * 100}%`, top: `${PLOT.top / height * 100}%`,
        width: `${PLOT.width / width * 100}%`, height: `${PLOT.height / height * 100}%`,
    });
    const plot = document.createElement("div");
    plot.className = "projection-plot";
    plot.style.maxWidth = `${width}px`;
    plot.append(svg, canvas);

    let drawnPoints = 0;
    let pending = false;
    const redraw = () => {
        pending = false;
        const rect = canvas.getBoundingClientRect();
        const ratio = window.devicePixelRatio || 1;
        const pixelWidth = Math.max(1, Math.round(rect.width * ratio));
        const pixelHeight = Math.max(1, Math.round(rect.height * ratio));
        if (canvas.width !== pixelWidth || canvas.height !== pixelHeight) {
            canvas.width = pixelWidth;
            canvas.height = pixelHeight;
        }
        const context = canvas.getContext("2d");
        context.setTransform(1, 0, 0, 1, 0, 0);
        context.clearRect(0, 0, canvas.width, canvas.height);
        context.setTransform(canvas.width / PLOT.width, 0, 0, canvas.height / PLOT.height, 0, 0);
        drawnPoints = drawPoints(context, state);
        drawPlotAxes(axesLayer, state);
    };
    const scheduleRedraw = () => {
        if (!pending) {
            pending = true;
            requestAnimationFrame(redraw);
        }
    };
    const drawLegendEntries = () => {
        if (data.color_by !== null) {
            drawLegend(legendLayer, state, state.groups.map((group, g) => g), true);
        }
    };

    // The pointer's position in plot-area units.
    const plotPoint = event => {
        const rect = canvas.getBoundingClientRect();
        return {
            x: (event.clientX - rect.left) / rect.width * PLOT.width,
            y: (event.clientY - rect.top) / rect.height * PLOT.height,
            perPixel: PLOT.width / rect.width,
        };
    };
    const zoomAt = (factor, x, y) => {
        const {view} = state;
        const dataX = view.cx + (x - PLOT.width / 2) * view.u;
        const dataY = view.cy - (y - PLOT.height / 2) * view.u;
        const u = view.u / factor;
        state.view = {cx: dataX - (dataX - view.cx) * u / view.u, cy: dataY - (dataY - view.cy) * u / view.u, u};
        scheduleRedraw();
    };
    const nearestPoint = event => {
        const point = plotPoint(event);
        const reach = state.radius + HOVER_SLOP * point.perPixel;
        let best = null;
        let bestDistance = reach * reach;
        forEachShownPoint(state, (i, px, py) => {
            const distance = (px - point.x) ** 2 + (py - point.y) ** 2;
            if (distance <= bestDistance) {
                best = i;
                bestDistance = distance;
            }
        });
        return best;
    };

    // Wheel zooms around the pointer, dragging pans, and double-click zooms in, as in plot_contigs.
    const tooltip = chartTooltip();
    let press = null;
    let dragged = false;
    canvas.addEventListener("wheel", event => {
        event.preventDefault();
        const point = plotPoint(event);
        const delta = -event.deltaY * (event.deltaMode === 1 ? 0.05 : event.deltaMode ? 1 : 0.002);
        zoomAt(2 ** delta, point.x, point.y);
    }, {passive: false});
    canvas.addEventListener("dblclick", event => {
        const point = plotPoint(event);
        zoomAt(event.shiftKey ? 0.5 : 2, point.x, point.y);
    });
    canvas.addEventListener("pointerdown", event => {
        if (event.button !== 0) {
            return;
        }
        press = {id: event.pointerId, clientX: event.clientX, clientY: event.clientY, from: plotPoint(event), view: state.view};
        dragged = false;
    });
    canvas.addEventListener("pointermove", event => {
        if (press && event.pointerId === press.id) {
            if (!dragged && Math.hypot(event.clientX - press.clientX, event.clientY - press.clientY) >= CLICK_SLOP) {
                dragged = true;
                canvas.setPointerCapture(event.pointerId);
                canvas.classList.add("panning");
                tooltip.hide();
            }
            if (dragged) {
                const point = plotPoint(event);
                state.view = {
                    ...press.view,
                    cx: press.view.cx - (point.x - press.from.x) * press.view.u,
                    cy: press.view.cy + (point.y - press.from.y) * press.view.u,
                };
                scheduleRedraw();
                return;
            }
        }
        const i = nearestPoint(event);
        if (i === null) {
            tooltip.hide();
        } else {
            tooltip.show(event, tooltipText(state, i));
        }
    });
    const release = () => {
        press = null;
        dragged = false;
        canvas.classList.remove("panning");
    };
    canvas.addEventListener("pointerup", release);
    canvas.addEventListener("pointercancel", release);
    canvas.addEventListener("pointerleave", () => tooltip.hide());

    const toggleGroup = entry => {
        const g = Number(entry.dataset.group);
        if (state.hidden.has(g)) {
            state.hidden.delete(g);
        } else {
            state.hidden.add(g);
        }
        drawLegendEntries();
        legendLayer.querySelector(`[data-group="${g}"]`).focus();
        scheduleRedraw();
    };
    legendLayer.addEventListener("click", event => {
        const entry = event.target.closest("[data-group]");
        if (entry) {
            toggleGroup(entry);
        }
    });
    legendLayer.addEventListener("keydown", event => {
        const entry = event.target.closest("[data-group]");
        if (entry && (event.key === "Enter" || event.key === " ")) {
            event.preventDefault();
            toggleGroup(entry);
        }
    });

    const fileStem = () => (data.axes.length === 3 ? `projection_${state.pair.map(a => AXIS_NAMES[a]).join("")}` : "projection");
    const toolbar = document.createElement("div");
    toolbar.className = "toolbar";
    if (data.axes.length === 3) {
        const pairButtons = [[0, 1], [0, 2], [1, 2]].map(pair => {
            const names = pair.map(a => AXIS_NAMES[a]);
            const button = toolbarButton(`axes-${names.join("").toLowerCase()}`, names.join(" / "), () => {
                state.pair = pair;
                state.view = fittedView(state);
                pairButtons.forEach(other => other.setAttribute("aria-pressed", String(other === button)));
                scheduleRedraw();
            });
            button.setAttribute("aria-pressed", String(pair[1] === 1));
            return button;
        });
        toolbar.append(...pairButtons);
    }

    const sizeLabel = document.createElement("label");
    sizeLabel.textContent = "Point size ";
    const sizeInput = document.createElement("input");
    Object.assign(sizeInput, {type: "number", id: "point-size", className: "point-size", min: 0.5, max: 20, step: 0.5, value: state.radius});
    sizeInput.addEventListener("input", () => {
        const radius = Number(sizeInput.value);
        if (radius > 0) {
            state.radius = radius;
            scheduleRedraw();
        }
    });
    sizeLabel.appendChild(sizeInput);

    const searchInput = document.createElement("input");
    Object.assign(searchInput, {type: "search", id: "point-search", className: "point-search", placeholder: "Find a point by name"});
    const searchStatus = document.createElement("span");
    searchStatus.id = "search-status";
    searchStatus.className = "status";
    const find = () => {
        const text = searchInput.value.trim();
        if (text === "") {
            state.highlight = null;
            searchStatus.textContent = "";
        } else {
            const lower = text.toLowerCase();
            let i = data.ids.indexOf(text);
            if (i < 0) {
                i = data.ids.findIndex(id => id.toLowerCase().includes(lower));
            }
            if (i < 0) {
                state.highlight = null;
                searchStatus.textContent = "No point matches.";
            } else {
                state.highlight = i;
                const [xs, ys] = state.pair.map(axis => data.coords[axis]);
                state.view = {...state.view, cx: xs[i], cy: ys[i]};
                searchStatus.textContent = data.ids[i];
            }
        }
        scheduleRedraw();
    };
    searchInput.addEventListener("keydown", event => {
        if (event.key === "Enter") {
            find();
        }
    });

    const spacer = document.createElement("span");
    spacer.className = "spacer";
    toolbar.append(
        sizeLabel, searchInput, toolbarButton("find-point", "Find", find), searchStatus, spacer,
        toolbarButton("zoom-in", "Zoom in", () => zoomAt(1.5, PLOT.width / 2, PLOT.height / 2)),
        toolbarButton("zoom-out", "Zoom out", () => zoomAt(1 / 1.5, PLOT.width / 2, PLOT.height / 2)),
        toolbarButton("reset-view", "Reset view", () => {
            state.view = fittedView(state);
            scheduleRedraw();
        }),
        toolbarButton("download-svg", "Download SVG",
            () => downloadText(`${fileStem()}.svg`, "image/svg+xml", svgMarkup(exportSvg(state)))),
        toolbarButton("download-png", "Download PNG", () => exportPng(state, fileStem())),
    );

    container.replaceChildren(toolbar, plot);
    drawLegendEntries();
    redraw();
    window.addEventListener("resize", scheduleRedraw);

    return {
        drawnPoints: () => drawnPoints,
        domain: () => {
            const {x0, x1, y0, y1} = plotScales(state.view);
            return {x0, x1, y0, y1};
        },
        axes: () => state.pair.map(a => data.axes[a]),
        // Where point `id` is on the screen, in client coordinates.
        pointPosition: id => {
            const i = data.ids.indexOf(id);
            const {sx, sy} = plotScales(state.view);
            const rect = canvas.getBoundingClientRect();
            return {
                x: rect.left + sx(data.coords[state.pair[0]][i]) / PLOT.width * rect.width,
                y: rect.top + sy(data.coords[state.pair[1]][i]) / PLOT.height * rect.height,
            };
        },
    };
}

(async () => {
    const root = document.getElementById("report-root");
    try {
        const {data} = await readReportData();
        const section = document.createElement("section");
        section.className = "panel";
        root.replaceChildren(section);
        if (data.ids.length === 0) {
            const status = document.createElement("div");
            status.className = "status";
            status.textContent = "No points to plot.";
            section.appendChild(status);
            return;
        }
        const note = document.createElement("p");
        note.className = "note";
        const colored = data.color_by === null ? "" : `, colored by ${data.color_by} (click a legend entry to hide or show it)`;
        note.textContent = `${formatChartNumber(data.ids.length)} points projected with ${data.algorithm}${colored}. `
            + "Point at a point to see its metadata. Zoom with the mouse wheel or double-click (shift to zoom out), and drag to pan.";
        const plot = document.createElement("div");
        plot.id = "projection";
        section.append(note, plot);
        window.projectionView = drawProjection(plot, data);
    } catch (error) {
        showReportError(root, error);
    }
})();
