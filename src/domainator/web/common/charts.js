// SVG charts drawn without a charting library: axes with round-number ticks, a hover tooltip,
// a histogram of integer values that zooms and pans, and buttons that save a chart as SVG or PNG.
//
// Colors and fonts are set as attributes, not by the page's stylesheet, so a saved SVG looks
// the same outside the page.

const SVG_NS = "http://www.w3.org/2000/svg";
const CHART_FONT = "ui-sans-serif, system-ui, -apple-system, \"Segoe UI\", sans-serif";
const CHART_INK = "#1e2a2f";
const CHART_MUTED = "#5c6a70";
const CHART_GRID = "#e3e6ea";
const CHART_ACCENT = "#c8553d";

function svgElement(tag, attributes = {}, parent = null) {
    const element = document.createElementNS(SVG_NS, tag);
    for (const [name, value] of Object.entries(attributes)) {
        element.setAttribute(name, String(value));
    }
    if (parent) {
        parent.appendChild(element);
    }
    return element;
}

// The step, 1, 2 or 5 times a power of ten, that splits [lo, hi] into about `count` parts.
function niceStep(lo, hi, count) {
    const span = hi - lo;
    if (!(span > 0)) {
        return 1;
    }
    const raw = span / Math.max(1, count);
    const power = 10 ** Math.floor(Math.log10(raw));
    for (const factor of [1, 2, 5, 10]) {
        if (power * factor >= raw) {
            return power * factor;
        }
    }
    return power * 10;
}

// The multiples of niceStep(lo, hi, count) in [lo, hi].
function niceTicks(lo, hi, count) {
    const step = niceStep(lo, hi, count);
    const ticks = [];
    for (let i = Math.ceil(lo / step); i * step <= hi + step * 1e-9; i += 1) {
        // Rounded so 0.1 * 3 is shown as 0.3.
        ticks.push(Number((i * step).toPrecision(12)));
    }
    return ticks;
}

// The axis from 0 up to the first multiple of niceStep(0, top, count) at or above top, and its
// ticks: [axisTop, ticks]. integer: true keeps whole-number ticks only.
function zeroBasedAxis(top, count = 5, integer = false) {
    const step = integer ? Math.max(1, niceStep(0, top, count)) : niceStep(0, top, count);
    const axisTop = Math.max(step, Math.ceil(top / step - 1e-9) * step);
    return [axisTop, niceTicks(0, axisTop, Math.round(axisTop / step))];
}

function formatChartNumber(value) {
    return value.toLocaleString("en-US", {maximumFractionDigits: 3});
}

function linearScale(domainLo, domainHi, rangeLo, rangeHi) {
    const span = domainHi - domainLo || 1;
    return value => rangeLo + (value - domainLo) / span * (rangeHi - rangeLo);
}

// One tooltip per page, shown beside the pointer.
function chartTooltip() {
    let tooltip = document.querySelector(".chart-tooltip");
    if (!tooltip) {
        tooltip = document.createElement("div");
        tooltip.className = "chart-tooltip";
        tooltip.setAttribute("role", "tooltip");
        tooltip.hidden = true;
        document.body.appendChild(tooltip);
    }
    return {
        show(event, text) {
            tooltip.textContent = text;
            tooltip.hidden = false;
            const x = Math.min(event.clientX + 14, window.innerWidth - tooltip.offsetWidth - 4);
            const y = Math.min(event.clientY + 14, window.innerHeight - tooltip.offsetHeight - 4);
            tooltip.style.left = `${Math.max(0, x)}px`;
            tooltip.style.top = `${Math.max(0, y)}px`;
        },
        hide() {
            tooltip.hidden = true;
        },
    };
}

// Shows `text(element)` while the pointer is over an element of `group` that has data-tip.
function attachTooltip(group, text) {
    const tooltip = chartTooltip();
    group.addEventListener("pointermove", event => {
        const target = event.target.closest("[data-tip]");
        if (target && group.contains(target)) {
            tooltip.show(event, text(target));
        } else {
            tooltip.hide();
        }
    });
    group.addEventListener("pointerleave", () => tooltip.hide());
}

function svgMarkup(svg) {
    const copy = svg.cloneNode(true);
    copy.setAttribute("xmlns", SVG_NS);
    const box = svg.viewBox.baseVal;
    copy.setAttribute("width", String(box.width));
    copy.setAttribute("height", String(box.height));
    copy.removeAttribute("style");
    return `<?xml version="1.0" encoding="UTF-8"?>\n${new XMLSerializer().serializeToString(copy)}\n`;
}

function downloadSvg(svg, fileStem) {
    downloadText(`${fileStem}.svg`, "image/svg+xml", svgMarkup(svg));
}

// At twice the chart's size, on white.
function downloadPng(svg, fileStem, scale = 2) {
    const box = svg.viewBox.baseVal;
    const url = URL.createObjectURL(new Blob([svgMarkup(svg)], {type: "image/svg+xml"}));
    const image = new Image();
    image.onload = () => {
        const canvas = document.createElement("canvas");
        canvas.width = Math.round(box.width * scale);
        canvas.height = Math.round(box.height * scale);
        const context = canvas.getContext("2d");
        context.fillStyle = "#ffffff";
        context.fillRect(0, 0, canvas.width, canvas.height);
        context.drawImage(image, 0, 0, canvas.width, canvas.height);
        URL.revokeObjectURL(url);
        canvas.toBlob(blob => downloadBlob(`${fileStem}.png`, blob), "image/png");
    };
    image.src = url;
}

// "Download SVG" and "Download PNG" buttons for the chart `getSvg()` returns when clicked.
function chartDownloadButtons(getSvg, fileStem, idPrefix = "") {
    return [["SVG", downloadSvg], ["PNG", downloadPng]].map(([label, save]) => {
        const button = document.createElement("button");
        button.type = "button";
        button.id = `${idPrefix}download-${label.toLowerCase()}`;
        button.textContent = `Download ${label}`;
        button.addEventListener("click", () => save(getSvg(), fileStem));
        return button;
    });
}

// A chart's toolbar (its download buttons) above the chart.
function chartFigure(container, svg, fileStem, idPrefix) {
    const toolbar = document.createElement("div");
    toolbar.className = "toolbar";
    toolbar.append(...chartDownloadButtons(() => svg, fileStem, idPrefix));
    const figure = document.createElement("div");
    figure.className = "chart";
    figure.appendChild(svg);
    container.replaceChildren(toolbar, figure);
}

function drawAxes(svg, {x, y, xTicks, yTicks, left, right, top, bottom, xLabel, yLabel}) {
    const axes = svgElement("g", {"font-family": CHART_FONT, "font-size": 11, fill: CHART_MUTED}, svg);
    for (const tick of yTicks) {
        const ty = y(tick);
        svgElement("line", {x1: left, x2: right, y1: ty, y2: ty, stroke: CHART_GRID}, axes);
        const label = svgElement("text", {x: left - 6, y: ty, "text-anchor": "end", "dominant-baseline": "middle"}, axes);
        label.textContent = formatChartNumber(tick);
    }
    svgElement("line", {x1: left, x2: right, y1: bottom, y2: bottom, stroke: CHART_INK}, axes);
    for (const tick of xTicks) {
        const tx = x(tick);
        svgElement("line", {x1: tx, x2: tx, y1: bottom, y2: bottom + 4, stroke: CHART_INK}, axes);
        const label = svgElement("text", {x: tx, y: bottom + 16, "text-anchor": "middle"}, axes);
        label.textContent = formatChartNumber(tick);
    }
    const titles = svgElement("g", {"font-family": CHART_FONT, "font-size": 12, fill: CHART_INK}, svg);
    const xTitle = svgElement("text", {x: (left + right) / 2, y: bottom + 36, "text-anchor": "middle"}, titles);
    xTitle.textContent = xLabel;
    const yTitle = svgElement("text", {
        transform: `translate(14 ${(top + bottom) / 2}) rotate(-90)`, "text-anchor": "middle",
    }, titles);
    yTitle.textContent = yLabel;
}

// Histograms of integer values aim for about sqrt(n) bins, at least Sturges' ceil(log2(n)) + 1
// and at most this.
const MAX_HISTOGRAM_BINS = 50;
// A histogram never zooms in to fewer values than this across.
const MIN_HISTOGRAM_SPAN = 10;
// How far the pointer moves before a press is a drag rather than a click.
const CHART_DRAG_SLOP = 3;

// The smallest 1, 2 or 5 times a power of ten that is at least raw (>= 1).
function niceBinWidth(raw) {
    for (let power = 1; ; power *= 10) {
        for (const factor of [1, 2, 5]) {
            if (factor * power >= raw) {
                return factor * power;
            }
        }
    }
}

// The bin width for n integer values from lo to hi: niceBinWidth(ceil((hi - lo + 1) / target)),
// target = min(MAX_HISTOGRAM_BINS, max(bit length of (n - 1) + 1, floor(sqrt(n)))).
function autoBinWidth(lo, hi, n) {
    const bitLength = n > 1 ? (n - 1).toString(2).length : 0;
    const target = Math.min(MAX_HISTOGRAM_BINS, Math.max(bitLength + 1, Math.floor(Math.sqrt(n))));
    return niceBinWidth(Math.ceil((hi - lo + 1) / target));
}

// The first index in sorted values at or above value.
function lowerBound(values, value) {
    let lo = 0;
    let hi = values.length;
    while (lo < hi) {
        const mid = (lo + hi) >> 1;
        if (values[mid] < value) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    return lo;
}

// Bins of the given width over [lo, hi): {edges, counts}, the first edge lo rounded down to a
// multiple of width, bin i counting the values in [edges[i], edges[i + 1]). values are sorted
// distinct integers, each occurring counts[i] times.
function binValues(values, counts, lo, hi, width) {
    const start = Math.floor(lo / width) * width;
    const bins = Math.max(1, Math.ceil((hi - start) / width));
    const binCounts = new Array(bins).fill(0);
    for (let i = lowerBound(values, start); i < values.length && values[i] < start + bins * width; i += 1) {
        binCounts[Math.floor((values[i] - start) / width)] += counts[i];
    }
    return {edges: Array.from({length: bins + 1}, (_, i) => start + i * width), counts: binCounts};
}

// A histogram of integer values: values are sorted and distinct, each occurring counts[i]
// times. The wheel zooms the x axis around the pointer, dragging pans it, and double-click
// zooms in (shift to zoom out). The bin width is chosen for the values in view (autoBinWidth)
// unless one is typed into the toolbar. Returns {svg, view(), bins(), setView(range),
// setBinWidth(width)}, width null for the automatic one.
function drawHistogram(container, {values, counts, xLabel, yLabel, itemLabel = "items",
    fileStem = "histogram", idPrefix = ""}) {
    const width = 760;
    const height = 300;
    const margin = {top: 12, right: 20, bottom: 48, left: 64};
    const plot = {left: margin.left, right: width - margin.right, top: margin.top, bottom: height - margin.bottom};
    const svg = svgElement("svg", {
        viewBox: `0 0 ${width} ${height}`, role: "img", "aria-label": `Histogram of ${xLabel}`,
        style: `width: 100%; max-width: ${width}px; height: auto;`,
    });
    chartFigure(container, svg, fileStem, idPrefix);
    svg.style.cursor = "grab";
    svg.style.touchAction = "none";

    const lo = values[0];
    const hi = values[values.length - 1];
    const total = counts.reduce((a, b) => a + b, 0);
    const fullWidth = autoBinWidth(lo, hi, total);
    const fullStart = Math.floor(lo / fullWidth) * fullWidth;
    // The whole histogram: the automatic bins over every value.
    const domain = [fullStart, fullStart + (Math.floor((hi - fullStart) / fullWidth) + 1) * fullWidth];
    let view = null;
    let chosenWidth = null;
    let shown = null;

    const toolbar = container.querySelector(".toolbar");
    const resetButton = document.createElement("button");
    resetButton.type = "button";
    resetButton.id = `${idPrefix}reset-zoom`;
    resetButton.textContent = "Reset zoom";
    resetButton.addEventListener("click", () => setView(null));
    const widthLabel = document.createElement("label");
    widthLabel.textContent = "Bin width ";
    const widthInput = document.createElement("input");
    widthInput.type = "number";
    widthInput.id = `${idPrefix}bin-width`;
    widthInput.className = "bin-width";
    widthInput.min = "1";
    widthInput.step = "1";
    widthInput.placeholder = "auto";
    widthInput.title = "Leave empty to fit the bins to the values in view.";
    widthInput.addEventListener("change", () => {
        const typed = Math.round(Number(widthInput.value));
        chosenWidth = widthInput.value !== "" && typed >= 1 ? typed : null;
        widthInput.value = chosenWidth === null ? "" : String(chosenWidth);
        render();
    });
    widthLabel.appendChild(widthInput);
    const status = document.createElement("span");
    status.className = "status";
    status.id = `${idPrefix}bin-status`;
    toolbar.append(resetButton, widthLabel, status);

    function render() {
        const [viewLo, viewHi] = view || domain;
        let binWidth = chosenWidth;
        if (binWidth === null) {
            // Fit to the values in view; the whole histogram keeps the bins it opened with.
            const first = Math.ceil(viewLo);
            const last = Math.ceil(viewHi) - 1;
            const a = lowerBound(values, first);
            const b = lowerBound(values, last + 1);
            let n = 0;
            for (let i = a; i < b; i += 1) {
                n += counts[i];
            }
            binWidth = view ? autoBinWidth(first, Math.max(first, last), Math.max(1, n)) : fullWidth;
        }
        shown = binValues(values, counts, viewLo, viewHi, binWidth);
        svg.replaceChildren();
        const x = linearScale(viewLo, viewHi, plot.left, plot.right);
        const maxCount = Math.max(1, ...shown.counts);
        // Counts get whole-number ticks.
        const [yTop, yTicks] = zeroBasedAxis(maxCount, 5, true);
        const y = linearScale(0, yTop, plot.bottom, plot.top);
        const xTicks = niceTicks(viewLo, viewHi, 8).filter(Number.isInteger);
        drawAxes(svg, {x, y, xTicks, yTicks, xLabel, yLabel, ...plot});
        const bars = svgElement("g", {fill: CHART_ACCENT}, svg);
        shown.counts.forEach((count, i) => {
            const start = shown.edges[i];
            const end = shown.edges[i + 1];
            // Bins cut by the edge of the view are drawn up to it.
            const x0 = Math.max(plot.left, x(start));
            const x1 = Math.min(plot.right, x(end));
            if (x1 - x0 < 0.5) {
                return;
            }
            const bar = svgElement("rect", {
                x: x0 + 0.5, y: y(count), width: Math.max(1, x1 - x0 - 1), height: y(0) - y(count),
                "fill-opacity": 0.85,
            }, bars);
            const range = end - start === 1 ? formatChartNumber(start) : `${formatChartNumber(start)}–${formatChartNumber(end - 1)}`;
            bar.dataset.tip = `${range}: ${formatChartNumber(count)} ${itemLabel}`;
            // Empty bins can still be pointed at: a transparent strip the height of the plot.
            if (count === 0) {
                bar.setAttribute("y", plot.top);
                bar.setAttribute("height", y(0) - plot.top);
                bar.setAttribute("fill-opacity", 0);
            }
        });
        attachTooltip(bars, bar => bar.dataset.tip);
        resetButton.disabled = view === null;
        status.textContent = `${formatChartNumber(shown.counts.length)} bin${shown.counts.length === 1 ? "" : "s"}, ${formatChartNumber(binWidth)} wide`;
    }

    // Keeps a view inside the whole histogram and at least MIN_HISTOGRAM_SPAN across; null is
    // the whole histogram.
    function setView(range) {
        if (range !== null) {
            const full = domain[1] - domain[0];
            const span = Math.min(full, Math.max(MIN_HISTOGRAM_SPAN, range[1] - range[0]));
            let start = Math.min(domain[1] - span, Math.max(domain[0], range[0]));
            if (span >= full) {
                start = domain[0];
            }
            range = span >= full ? null : [start, start + span];
        }
        view = range;
        chartTooltip().hide();
        render();
    }

    const toChart = event => new DOMPoint(event.clientX, event.clientY).matrixTransform(svg.getScreenCTM().inverse()).x;
    const valueAt = px => {
        const [viewLo, viewHi] = view || domain;
        return viewLo + (px - plot.left) / (plot.right - plot.left) * (viewHi - viewLo);
    };
    const zoomAt = (factor, px) => {
        const [viewLo, viewHi] = view || domain;
        const at = valueAt(Math.max(plot.left, Math.min(plot.right, px)));
        setView([at - (at - viewLo) / factor, at + (viewHi - at) / factor]);
    };
    svg.addEventListener("wheel", event => {
        event.preventDefault();
        const delta = -event.deltaY * (event.deltaMode === 1 ? 0.05 : event.deltaMode ? 1 : 0.002);
        zoomAt(2 ** delta, toChart(event));
    }, {passive: false});
    svg.addEventListener("dblclick", event => zoomAt(event.shiftKey ? 0.5 : 2, toChart(event)));
    let press = null;
    svg.addEventListener("pointerdown", event => {
        if (event.button === 0) {
            press = {id: event.pointerId, clientX: event.clientX, from: toChart(event), view: view || domain, dragged: false};
        }
    });
    svg.addEventListener("pointermove", event => {
        if (!press || event.pointerId !== press.id) {
            return;
        }
        if (!press.dragged) {
            if (Math.abs(event.clientX - press.clientX) < CHART_DRAG_SLOP) {
                return;
            }
            press.dragged = true;
            svg.setPointerCapture(event.pointerId);
            svg.style.cursor = "grabbing";
        }
        const [viewLo, viewHi] = press.view;
        const shift = (toChart(event) - press.from) / (plot.right - plot.left) * (viewHi - viewLo);
        setView([viewLo - shift, viewHi - shift]);
    });
    const release = () => {
        press = null;
        svg.style.cursor = "grab";
    };
    svg.addEventListener("pointerup", release);
    svg.addEventListener("pointercancel", release);

    render();
    return {
        svg,
        view: () => (view ? view.slice() : null),
        bins: () => ({edges: shown.edges.slice(), counts: shown.counts.slice()}),
        setView,
        setBinWidth(binWidth) {
            widthInput.value = binWidth === null ? "" : String(binWidth);
            widthInput.dispatchEvent(new Event("change"));
        },
    };
}
