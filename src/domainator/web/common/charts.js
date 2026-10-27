// SVG charts drawn without a charting library: axes with round-number ticks, a hover tooltip,
// a histogram of precomputed bins, and buttons that save a chart as SVG or PNG.
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

// binEdges: n + 1 numbers, bin i holding the values in [binEdges[i], binEdges[i + 1]);
// counts: n counts. With integerBins, a bin is labeled by the integers it holds.
function drawHistogram(container, {binEdges, counts, xLabel, yLabel, itemLabel = "items",
    integerBins = false, fileStem = "histogram", idPrefix = ""}) {
    const width = 760;
    const height = 300;
    const margin = {top: 12, right: 20, bottom: 48, left: 64};
    const svg = svgElement("svg", {
        viewBox: `0 0 ${width} ${height}`, role: "img", "aria-label": `Histogram of ${xLabel}`,
        style: `width: 100%; max-width: ${width}px; height: auto;`,
    });
    const lo = binEdges[0];
    const hi = binEdges[binEdges.length - 1];
    const maxCount = Math.max(1, ...counts);
    // Counts get whole-number ticks.
    const yTicks = niceTicks(0, maxCount, 5).filter(Number.isInteger);
    const yTop = Math.max(maxCount, yTicks[yTicks.length - 1]);
    const x = linearScale(lo, hi, margin.left, width - margin.right);
    const y = linearScale(0, yTop, height - margin.bottom, margin.top);
    // Integer bins get integer ticks.
    const xTicks = niceTicks(lo, hi, 8).filter(tick => !integerBins || Number.isInteger(tick));
    drawAxes(svg, {
        x, y, xTicks, yTicks, xLabel, yLabel,
        left: margin.left, right: width - margin.right, top: margin.top, bottom: height - margin.bottom,
    });

    const bars = svgElement("g", {fill: CHART_ACCENT}, svg);
    counts.forEach((count, i) => {
        const x0 = x(binEdges[i]);
        const x1 = x(binEdges[i + 1]);
        const bar = svgElement("rect", {
            x: x0 + 0.5, y: y(count), width: Math.max(1, x1 - x0 - 1), height: y(0) - y(count),
            "fill-opacity": 0.85,
        }, bars);
        const start = binEdges[i];
        const end = binEdges[i + 1];
        const range = integerBins
            ? (end - start === 1 ? formatChartNumber(start) : `${formatChartNumber(start)}–${formatChartNumber(end - 1)}`)
            : `${formatChartNumber(start)}–${formatChartNumber(end)}`;
        bar.dataset.tip = `${range}: ${formatChartNumber(count)} ${itemLabel}`;
        // Empty bins can still be pointed at: a transparent strip the height of the plot.
        if (count === 0) {
            bar.setAttribute("y", margin.top);
            bar.setAttribute("height", y(0) - margin.top);
            bar.setAttribute("fill-opacity", 0);
        }
    });
    attachTooltip(bars, bar => bar.dataset.tip);
    chartFigure(container, svg, fileStem, idPrefix);
    return svg;
}
