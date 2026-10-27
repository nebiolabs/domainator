// The contig plot (plot_contigs): each contig as a line, with its features drawn above it as
// arrows, stacked so that overlapping features do not cover each other.
//
// report-data's data is {types, color_by, layout, contigs}:
// - types are the feature types drawn, in stacking order: features of earlier types are drawn
//   closer to their contig.
// - color_by is the metadata column the contig lines are colored by, or null.
// - layout holds plot_contigs' size options: width, height, feature_height, arrow_width,
//   contig_thickness, contig_label_font, contig_label_px, feature_label_font,
//   feature_font_size, details_width and details_height.
// - contigs are {name, length, color, details, features}. details are [label, text] pairs: the
//   record's name, description, length, annotations, cross-references and metadata. features
//   are {name, type, color, parts, qualifiers}, parts being [start, end, strand] with 0-based,
//   end-exclusive positions and strand "+" or "-", and qualifiers [key, [values]] pairs.

// @include common/payload.js
// @include common/download.js
// @include common/charts.js

// Feature types that share one level above a contig however they overlap, so that a contig's
// genes stay on one line.
const SQUASHED_TYPES = new Set(["CDS", "exon"]);
// Pixels between the bottom of a feature and its label's baseline.
const FEATURE_LABEL_INSET = 7;
// Pixels between a contig's label and the features of the contig below it.
const CONTIG_GAP = 5;
// Longest qualifier value a tooltip shows in full; the details box shows all of it.
const TOOLTIP_VALUE_LENGTH = 120;
// A press that moves less than this many pixels is a click, not a drag.
const CLICK_SLOP = 3;
const MIN_SCALE = 1e-4;
const MAX_SCALE = 100;

function featureSpan(feature) {
    return {
        start: Math.min(...feature.parts.map(part => part[0])),
        end: Math.max(...feature.parts.map(part => part[1])),
    };
}

// True if the features share a position (touching does not count).
function partsOverlap(a, b) {
    return a.parts.some(([aStart, aEnd]) => b.parts.some(([bStart, bEnd]) => Math.min(aEnd, bEnd) - Math.max(aStart, bStart) > 0));
}

// The features with their span and a level, 0 being next to the contig. Features are placed
// by type (in `types` order), then by start, each on the lowest level that no feature placed
// before it and overlapping it is on; but a SQUASHED_TYPES feature goes on the level of the
// first one placed.
function stackFeatures(features, types) {
    const rank = new Map(types.map((type, i) => [type, i]));
    const order = type => (rank.has(type) ? rank.get(type) : types.length + 1);
    const stacked = features.map(feature => ({...feature, ...featureSpan(feature), level: 0}))
        .sort((a, b) => order(a.type) - order(b.type) || a.start - b.start);
    const placed = [];
    for (const feature of stacked) {
        const occupied = new Set();
        let shared = null;
        for (const other of placed) {
            if (SQUASHED_TYPES.has(feature.type) && SQUASHED_TYPES.has(other.type)) {
                shared = other.level;
                break;
            }
            if (partsOverlap(feature, other)) {
                occupied.add(other.level);
            }
        }
        if (shared !== null) {
            feature.level = shared;
        } else {
            while (occupied.has(feature.level)) {
                feature.level += 1;
            }
        }
        placed.push(feature);
    }
    return stacked;
}

// The contigs with stacked features, and each one's offset: the y of its line, below its
// features and the contig above's label.
function layoutContigs(contigs, types, layout) {
    const gap = layout.contig_label_px + layout.contig_thickness / 2 + CONTIG_GAP;
    let offset = 0;
    return contigs.map(contig => {
        const features = stackFeatures(contig.features, types);
        const levels = Math.max(0, ...features.map(feature => feature.level)) + 1;
        offset += levels * layout.feature_height + gap;
        return {...contig, features, offset};
    });
}

// An arrow `length` long and `height` high from (0, 0) up, pointing right on "+" and left on
// "-". Its head is at most as long as the arrow.
function arrowPath(length, strand, arrowWidth, height) {
    const head = Math.min(arrowWidth, length);
    if (strand === "-") {
        return `M${head},0 L${length},0 L${length},${-height} L${head},${-height} L0,${-height / 2} Z`;
    }
    return `M0,0 L${length - head},0 L${length},${-height / 2} L${length - head},${-height} L0,${-height} Z`;
}

function featureDetails(feature) {
    const strands = new Set(feature.parts.map(part => part[2]));
    return [
        ["start", String(feature.start)],
        ["end", String(feature.end)],
        ["strand", strands.size === 1 ? [...strands][0] : "mixed"],
        ...feature.qualifiers.map(([key, values]) => [key, values.join(";")]),
    ];
}

function tooltipText({heading, lines}) {
    const shorten = text => (text.length > TOOLTIP_VALUE_LENGTH ? `${text.slice(0, TOOLTIP_VALUE_LENGTH - 1)}…` : text);
    return [heading, ...lines.map(([label, text]) => `${label}: ${shorten(text)}`)].join("\n");
}

// Fills the details box with `lines`, the `swatch` label's value marked with its color.
function showDetails(box, {heading, lines, swatch}) {
    const title = document.createElement("div");
    title.className = "details-heading";
    title.textContent = heading;
    const list = document.createElement("dl");
    for (const [label, text] of lines) {
        const term = document.createElement("dt");
        term.textContent = label;
        const value = document.createElement("dd");
        if (swatch && label === swatch.label) {
            const mark = document.createElement("span");
            mark.className = "color-swatch";
            mark.style.background = swatch.color;
            value.appendChild(mark);
        }
        value.appendChild(document.createTextNode(text));
        list.append(term, value);
    }
    box.replaceChildren(title, list);
}

// Wheel zooms around the pointer, dragging pans, and double-click zooms in, as d3-zoom does.
// A drag that ends on a feature is not a click on it.
function panZoom(svg, layer, initialScale, center) {
    let view = {k: initialScale, x: 0, y: 0};
    const apply = () => layer.setAttribute("transform", `translate(${view.x},${view.y}) scale(${view.k})`);
    const toSvg = event => {
        const point = svg.createSVGPoint();
        point.x = event.clientX;
        point.y = event.clientY;
        return point.matrixTransform(svg.getScreenCTM().inverse());
    };
    const zoomAt = (factor, x, y) => {
        const k = Math.min(MAX_SCALE, Math.max(MIN_SCALE, view.k * factor));
        view = {k, x: x - (x - view.x) * k / view.k, y: y - (y - view.y) * k / view.k};
        apply();
    };
    svg.addEventListener("wheel", event => {
        event.preventDefault();
        const point = toSvg(event);
        const delta = -event.deltaY * (event.deltaMode === 1 ? 0.05 : event.deltaMode ? 1 : 0.002);
        zoomAt(2 ** delta, point.x, point.y);
    }, {passive: false});
    svg.addEventListener("dblclick", event => {
        const point = toSvg(event);
        zoomAt(event.shiftKey ? 0.5 : 2, point.x, point.y);
    });

    let press = null;
    let dragged = false;
    svg.addEventListener("pointerdown", event => {
        if (event.button !== 0) {
            return;
        }
        press = {id: event.pointerId, clientX: event.clientX, clientY: event.clientY, from: toSvg(event), view};
        dragged = false;
    });
    svg.addEventListener("pointermove", event => {
        if (!press || event.pointerId !== press.id) {
            return;
        }
        if (!dragged) {
            if (Math.hypot(event.clientX - press.clientX, event.clientY - press.clientY) < CLICK_SLOP) {
                return;
            }
            dragged = true;
            svg.setPointerCapture(event.pointerId);
            svg.classList.add("panning");
        }
        const point = toSvg(event);
        view = {...press.view, x: press.view.x + point.x - press.from.x, y: press.view.y + point.y - press.from.y};
        apply();
    });
    const release = () => {
        press = null;
        svg.classList.remove("panning");
    };
    svg.addEventListener("pointerup", release);
    svg.addEventListener("pointercancel", release);
    svg.addEventListener("click", event => {
        if (dragged) {
            event.stopPropagation();
            dragged = false;
        }
    }, true);

    apply();
    return {
        zoomBy: factor => zoomAt(factor, center.x, center.y),
        reset: () => {
            view = {k: initialScale, x: 0, y: 0};
            apply();
        },
        view: () => ({...view}),
        isPanning: () => dragged && press !== null,
    };
}

function toolbarButton(id, label, onClick) {
    const button = document.createElement("button");
    button.type = "button";
    button.id = id;
    button.textContent = label;
    button.addEventListener("click", onClick);
    return button;
}

function drawContigPlot(container, data, detailsBox) {
    const {layout} = data;
    const contigs = layoutContigs(data.contigs, data.types, layout);
    const widest = Math.max(...contigs.map(contig => contig.length));
    // Fit the widest contig, but never shrink the plot below a tenth of its size.
    const initialScale = Math.max(0.1, Math.min(1, layout.width / widest));

    const svg = svgElement("svg", {
        viewBox: `0 0 ${layout.width} ${layout.height}`, role: "img",
        "aria-label": `Plot of ${contigs.length} contigs`, class: "contig-plot",
        style: `width: 100%; max-width: ${layout.width}px; height: auto;`,
    });
    // Something to press on between the contigs, so that every part of the plot pans.
    svgElement("rect", {width: layout.width, height: layout.height, fill: "#ffffff"}, svg);
    const layer = svgElement("g", {}, svg);
    contigs.forEach((contig, contigIndex) => {
        const group = svgElement("g", {transform: `translate(0,${contig.offset})`, "data-contig": contigIndex}, layer);
        svgElement("line", {
            x1: 0, x2: contig.length, y1: 0, y2: 0, stroke: contig.color, "stroke-width": layout.contig_thickness,
        }, group);
        contig.features.forEach((feature, featureIndex) => {
            const shape = svgElement("g", {"data-feature": featureIndex}, group);
            const y = -layout.feature_height * feature.level;
            for (const [start, end, strand] of feature.parts) {
                svgElement("path", {
                    d: arrowPath(end - start, strand, layout.arrow_width, layout.feature_height),
                    transform: `translate(${start},${y})`, fill: feature.color, stroke: "#000000",
                }, shape);
            }
            const label = svgElement("text", {
                x: feature.start + layout.arrow_width, y: y - FEATURE_LABEL_INSET, fill: "#000000",
                "font-family": layout.feature_label_font, "font-size": layout.feature_font_size,
            }, shape);
            label.textContent = feature.name;
        });
        const label = svgElement("text", {
            x: 0, y: layout.contig_label_px + layout.contig_thickness / 2, fill: "#000000",
            "font-family": layout.contig_label_font, "font-size": layout.contig_label_px,
        }, group);
        label.textContent = contig.name;
    });

    const describe = element => {
        const featureElement = element.closest("[data-feature]");
        const contigElement = element.closest("[data-contig]");
        if (!contigElement) {
            return null;
        }
        const contig = contigs[Number(contigElement.dataset.contig)];
        if (featureElement) {
            const feature = contig.features[Number(featureElement.dataset.feature)];
            return {heading: feature.type, lines: featureDetails(feature)};
        }
        const swatch = data.color_by === null ? null : {label: data.color_by, color: contig.color};
        return {heading: "Contig", lines: contig.details, swatch};
    };

    const view = panZoom(svg, layer, initialScale, {x: layout.width / 2, y: layout.height / 2});
    const tooltip = chartTooltip();
    layer.addEventListener("pointermove", event => {
        const described = view.isPanning() ? null : describe(event.target);
        if (described) {
            tooltip.show(event, tooltipText(described));
        } else {
            tooltip.hide();
        }
    });
    layer.addEventListener("pointerleave", () => tooltip.hide());
    layer.addEventListener("click", event => {
        const described = describe(event.target);
        if (described) {
            showDetails(detailsBox, described);
        }
    });

    chartFigure(container, svg, "contigs", "");
    container.querySelector(".toolbar").prepend(
        toolbarButton("zoom-in", "Zoom in", () => view.zoomBy(1.5)),
        toolbarButton("zoom-out", "Zoom out", () => view.zoomBy(1 / 1.5)),
        toolbarButton("reset-view", "Reset view", () => view.reset()),
    );
    return {contigs, view, svg};
}

(async () => {
    const root = document.getElementById("report-root");
    try {
        const {data} = await readReportData();
        const section = document.createElement("section");
        section.className = "panel";
        root.replaceChildren(section);
        if (data.contigs.length === 0) {
            const status = document.createElement("div");
            status.className = "status";
            status.textContent = "No contigs to plot.";
            section.appendChild(status);
            return;
        }
        const detailsBox = document.createElement("div");
        detailsBox.id = "details";
        detailsBox.className = "contig-details";
        detailsBox.style.width = `${data.layout.details_width}px`;
        detailsBox.style.height = `${data.layout.details_height}px`;
        detailsBox.textContent = "Click a contig or a feature to show its details here.";
        const note = document.createElement("p");
        note.className = "note";
        note.textContent = "Zoom with the mouse wheel or double-click (shift to zoom out), drag to pan, and click a contig or a feature to show its details above.";
        const plot = document.createElement("div");
        plot.id = "contig-plot";
        section.append(detailsBox, note, plot);
        window.contigPlot = drawContigPlot(plot, data, detailsBox);
    } catch (error) {
        showReportError(root, error);
    }
})();
