// Line charts on SVG, built on charts.js: a chart frame that redraws itself for a range of x
// (zoom by dragging across it), log or linear x axes, a second y axis on the right, and a
// tooltip with each series' value under the pointer.

// Ticks for a log axis over [lo, hi] (lo > 0): the powers of ten, with 2 and 5 times them when
// the range spans few decades, and round linear ticks when it spans less than one.
function logTicks(lo, hi) {
    const first = Math.floor(Math.log10(lo));
    const last = Math.ceil(Math.log10(hi));
    const factors = last - first <= 2 ? [1, 2, 5] : [1];
    const stride = Math.max(1, Math.ceil((last - first) / 8));
    const ticks = [];
    for (let power = first; power <= last; power += stride) {
        for (const factor of factors) {
            const tick = Number((factor * 10 ** power).toPrecision(12));
            if (tick >= lo * (1 - 1e-9) && tick <= hi * (1 + 1e-9)) {
                ticks.push(tick);
            }
        }
    }
    return ticks.length >= 2 ? ticks : niceTicks(lo, hi, 5);
}

// A scale for [lo, hi] onto [rangeLo, rangeHi], linear or by log10, with invert().
function axisScale(lo, hi, rangeLo, rangeHi, log = false) {
    const f = log ? Math.log10 : value => value;
    const a = f(lo);
    const b = f(hi);
    const span = b - a || 1;
    const scale = value => rangeLo + (f(value) - a) / span * (rangeHi - rangeLo);
    scale.invert = position => {
        const value = a + (position - rangeLo) / (rangeHi - rangeLo) * span;
        return log ? 10 ** value : value;
    };
    return scale;
}

// The y axis on the right: its ticks and title, without grid lines.
function drawRightAxis(svg, {y, ticks, right, top, bottom, label}) {
    const axis = svgElement("g", {"font-family": CHART_FONT, "font-size": 11, fill: CHART_MUTED}, svg);
    for (const tick of ticks) {
        const text = svgElement("text", {x: right + 6, y: y(tick), "dominant-baseline": "middle"}, axis);
        text.textContent = formatChartNumber(tick);
    }
    const title = svgElement("text", {
        transform: `translate(${right + 52} ${(top + bottom) / 2}) rotate(90)`, "text-anchor": "middle",
        "font-family": CHART_FONT, "font-size": 12, fill: CHART_INK,
    }, svg);
    title.textContent = label;
}

// A legend in one row above the plot: [{label, color, mark: "line" | "dot"}].
function drawLegend(svg, entries, {left, top}) {
    const legend = svgElement("g", {"font-family": CHART_FONT, "font-size": 11, fill: CHART_INK}, svg);
    let x = left;
    for (const {label, color, mark} of entries) {
        if (mark === "dot") {
            svgElement("circle", {cx: x + 7, cy: top, r: 4, fill: color}, legend);
        } else {
            svgElement("line", {x1: x, x2: x + 14, y1: top, y2: top, stroke: color, "stroke-width": 2}, legend);
        }
        const text = svgElement("text", {x: x + 19, y: top, "dominant-baseline": "middle"}, legend);
        text.textContent = label;
        x += 19 + label.length * 6.2 + 18;
    }
}

// An SVG path through points ([x, y] pairs, sorted by x). step: true holds each y until the
// next x, as a running total over discrete events does.
function linePath(points, x, y, step = false) {
    const parts = [];
    points.forEach(([px, py], i) => {
        const sx = x(px).toFixed(2);
        const sy = y(py).toFixed(2);
        if (i === 0) {
            parts.push(`M${sx},${sy}`);
        } else if (step) {
            parts.push(`H${sx}V${sy}`);
        } else {
            parts.push(`L${sx},${sy}`);
        }
    });
    return parts.join("");
}

// The index of the point in points (sorted by x) nearest to value.
function nearestIndex(points, value, key = point => point[0]) {
    let lo = 0;
    let hi = points.length - 1;
    while (hi - lo > 1) {
        const mid = (lo + hi) >> 1;
        if (key(points[mid]) < value) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return Math.abs(key(points[hi]) - value) < Math.abs(key(points[lo]) - value) ? hi : lo;
}

// The points of a sorted series inside [lo, hi], with one more on each side so a line runs to
// the plot's edges.
function pointsInRange(points, lo, hi) {
    let start = 0;
    while (start < points.length && points[start][0] < lo) {
        start += 1;
    }
    let end = points.length;
    while (end > start && points[end - 1][0] > hi) {
        end -= 1;
    }
    return points.slice(Math.max(0, start - 1), Math.min(points.length, end + 1));
}

// A chart that redraws itself: draw(svg, [lo, hi]) draws it for that x range and returns
// {plot: {left, right, top, bottom}, x, hover(px, py), click(px, py)}, hover giving
// {text, guide, dots: [{x, y, color}]} or null for the point under the pointer (in the
// chart's own coordinates), and click being optional.
//
// With zoom, dragging across the plot zooms to that range of x and double-clicking (or the
// Reset zoom button) zooms back out. setMarker(value) draws a line across the plot at x = value
// (null for none), the value the page is showing elsewhere.
function interactiveChart(container, {width, height, domain, draw, zoom = true, ariaLabel = "Chart",
    fileStem = "chart", idPrefix = ""}) {
    const svg = svgElement("svg", {
        viewBox: `0 0 ${width} ${height}`, role: "img", "aria-label": ariaLabel,
        style: `width: 100%; max-width: ${width}px; height: auto;`,
    });
    chartFigure(container, svg, fileStem, idPrefix);
    const resetButton = document.createElement("button");
    if (zoom) {
        resetButton.type = "button";
        resetButton.id = `${idPrefix}reset-zoom`;
        resetButton.textContent = "Reset zoom";
        resetButton.disabled = true;
        resetButton.title = "Drag across the chart to zoom in; double-click it to zoom back out.";
        resetButton.addEventListener("click", () => setView(null));
        container.querySelector(".toolbar").appendChild(resetButton);
    }
    const tooltip = chartTooltip();
    let view = null;
    let marker = null;
    let layout = null;
    let markerLayer = null;
    let hoverLayer = null;
    let brush = null;

    const toChart = event => {
        const point = new DOMPoint(event.clientX, event.clientY).matrixTransform(svg.getScreenCTM().inverse());
        return [point.x, point.y];
    };
    const inPlot = ([px, py]) => layout && px >= layout.plot.left && px <= layout.plot.right
        && py >= layout.plot.top && py <= layout.plot.bottom;

    function drawMarker() {
        markerLayer.replaceChildren();
        if (marker === null || marker === undefined) {
            return;
        }
        const mx = layout.x(marker);
        if (Number.isFinite(mx) && mx >= layout.plot.left - 0.5 && mx <= layout.plot.right + 0.5) {
            svgElement("line", {
                x1: mx, x2: mx, y1: layout.plot.top, y2: layout.plot.bottom,
                stroke: "#d62728", "stroke-width": 1.5, "stroke-dasharray": "4 3",
            }, markerLayer);
        }
    }

    function render() {
        svg.replaceChildren();
        layout = draw(svg, view || domain);
        const {left, right, top, bottom} = layout.plot;
        markerLayer = svgElement("g", {"pointer-events": "none"}, svg);
        hoverLayer = svgElement("g", {"pointer-events": "none"}, svg);
        drawMarker();
        const overlay = svgElement("rect", {
            x: left, y: top, width: Math.max(0, right - left), height: Math.max(0, bottom - top),
            fill: "transparent", "data-overlay": "",
        }, svg);
        overlay.style.cursor = zoom ? "crosshair" : (layout.click ? "pointer" : "default");
        resetButton.disabled = view === null;
    }

    function setView(range) {
        view = range;
        tooltip.hide();
        render();
    }

    function showHover(point) {
        hoverLayer.replaceChildren();
        const found = inPlot(point) && layout.hover ? layout.hover(point[0], point[1]) : null;
        if (!found) {
            tooltip.hide();
            return null;
        }
        if (found.guide !== undefined && found.guide !== null) {
            svgElement("line", {
                x1: found.guide, x2: found.guide, y1: layout.plot.top, y2: layout.plot.bottom,
                stroke: CHART_MUTED, "stroke-width": 1, "stroke-opacity": 0.5,
            }, hoverLayer);
        }
        for (const dot of found.dots || []) {
            svgElement("circle", {cx: dot.x, cy: dot.y, r: 4, fill: dot.color, stroke: "#ffffff", "stroke-width": 1.5}, hoverLayer);
        }
        return found;
    }

    svg.addEventListener("pointerdown", event => {
        const point = toChart(event);
        if (event.button !== 0 || !inPlot(point)) {
            return;
        }
        brush = {start: point[0], current: point[0], rect: null};
        svg.setPointerCapture(event.pointerId);
    });
    svg.addEventListener("pointermove", event => {
        const point = toChart(event);
        if (brush && zoom) {
            brush.current = Math.max(layout.plot.left, Math.min(layout.plot.right, point[0]));
            if (Math.abs(brush.current - brush.start) > 4) {
                tooltip.hide();
                hoverLayer.replaceChildren();
                if (!brush.rect) {
                    brush.rect = svgElement("rect", {fill: CHART_MUTED, "fill-opacity": 0.15, "pointer-events": "none"}, svg);
                }
                brush.rect.setAttribute("x", Math.min(brush.start, brush.current));
                brush.rect.setAttribute("width", Math.abs(brush.current - brush.start));
                brush.rect.setAttribute("y", layout.plot.top);
                brush.rect.setAttribute("height", layout.plot.bottom - layout.plot.top);
                return;
            }
        }
        const found = showHover(point);
        if (found) {
            tooltip.show(event, found.text);
        }
    });
    svg.addEventListener("pointerup", event => {
        if (!brush) {
            return;
        }
        const {start, current} = brush;
        brush = null;
        if (zoom && Math.abs(current - start) > 4) {
            const a = layout.x.invert(Math.min(start, current));
            const b = layout.x.invert(Math.max(start, current));
            setView([a, b]);
        } else if (layout.click) {
            layout.click(...toChart(event));
        }
    });
    svg.addEventListener("pointerleave", () => {
        if (!brush) {
            tooltip.hide();
            hoverLayer.replaceChildren();
        }
    });
    if (zoom) {
        svg.addEventListener("dblclick", () => setView(null));
    }

    render();
    return {
        svg,
        setMarker(value) {
            marker = value;
            drawMarker();
        },
        zoomTo(range) {
            setView(range);
        },
        layout: () => layout,
    };
}

// A line chart of series [{name, color, points, axis: "left" | "right", step, format}], each
// points a list of [x, y] sorted by x. Every y axis starts at zero. The tooltip shows each
// series' value at the point nearest the pointer.
function drawLineChart(container, {series, xLabel, yLabel, y2Label = null, xLog = false, xName = null,
    xFormat = formatChartNumber, width = 760, height = 320, fileStem = "chart", idPrefix = "",
    integerY = true}) {
    const margin = {top: 34, right: y2Label ? 72 : 20, bottom: 48, left: 72};
    const firsts = series.filter(s => s.points.length > 0).map(s => s.points[0][0]);
    const lasts = series.filter(s => s.points.length > 0).map(s => s.points[s.points.length - 1][0]);
    const domain = firsts.length > 0 ? [Math.min(...firsts), Math.max(...lasts)] : [0, 1];
    const draw = (svg, [lo, hi]) => {
        const plot = {left: margin.left, right: width - margin.right, top: margin.top, bottom: height - margin.bottom};
        const x = axisScale(lo, hi, plot.left, plot.right, xLog);
        const visible = series.map(s => pointsInRange(s.points, lo, hi));
        const yScaleFor = axis => {
            let top = 1;
            visible.forEach((points, i) => {
                if ((series[i].axis || "left") === axis) {
                    for (const point of points) {
                        top = Math.max(top, point[1]);
                    }
                }
            });
            const [yTop, ticks] = zeroBasedAxis(top, 5, integerY);
            return {y: axisScale(0, yTop, plot.bottom, plot.top), ticks};
        };
        const left = yScaleFor("left");
        const right = y2Label ? yScaleFor("right") : null;
        const xTicks = xLog ? logTicks(lo, hi) : niceTicks(lo, hi, 8);
        drawAxes(svg, {
            x, y: left.y, xTicks, yTicks: left.ticks, xLabel, yLabel, ...plot,
        });
        if (right) {
            drawRightAxis(svg, {y: right.y, ticks: right.ticks, right: plot.right, top: plot.top, bottom: plot.bottom, label: y2Label});
        }
        drawLegend(svg, series.map(s => ({label: s.name, color: s.color, mark: "line"})), {left: plot.left, top: 14});
        const clipId = `${idPrefix}clip`;
        const clip = svgElement("clipPath", {id: clipId}, svgElement("defs", {}, svg));
        svgElement("rect", {x: plot.left, y: plot.top - 2, width: plot.right - plot.left, height: plot.bottom - plot.top + 4}, clip);
        const lines = svgElement("g", {"clip-path": `url(#${clipId})`, fill: "none"}, svg);
        series.forEach((s, i) => {
            const y = (s.axis === "right" && right ? right : left).y;
            svgElement("path", {d: linePath(visible[i], x, y, s.step), stroke: s.color, "stroke-width": 2}, lines);
        });
        const hover = px => {
            const value = x.invert(px);
            const dots = [];
            const text = [];
            let guide = null;
            series.forEach((s, i) => {
                const points = visible[i];
                if (points.length === 0) {
                    return;
                }
                const scaled = xLog ? Math.log10 : v => v;
                const point = points[nearestIndex(points, scaled(value), p => scaled(p[0]))];
                const y = (s.axis === "right" && right ? right : left).y;
                if (guide === null) {
                    guide = x(point[0]);
                    text.push(`${xName || xLabel}: ${xFormat(point[0])}`);
                }
                dots.push({x: x(point[0]), y: y(point[1]), color: s.color});
                text.push(`${s.name}: ${(s.format || formatChartNumber)(point[1])}`);
            });
            return dots.length > 0 ? {text: text.join("\n"), guide, dots} : null;
        };
        return {plot, x, hover};
    };
    return interactiveChart(container, {
        width, height, domain, draw, ariaLabel: `${yLabel} by ${xLabel}`, fileStem, idPrefix,
    });
}
