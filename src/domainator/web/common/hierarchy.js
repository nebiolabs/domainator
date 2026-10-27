// Trees laid out for drawing, and a zoomable icicle chart of one.
//
// The icicle follows the "Zoomable icicle" example in Observable's d3 gallery (Copyright 2018-2023
// Observable, Inc., ISC license), drawn without d3.

// A tree to lay out: {name, count, children}. Each layout node is {data, parent, children,
// depth, height, value}, value being the node's own count plus its children's values.
function hierarchyNodes(tree, depth = 0, parent = null) {
    const node = {data: tree, parent, depth, children: [], value: tree.count || 0, height: 0};
    for (const child of tree.children || []) {
        const childNode = hierarchyNodes(child, depth + 1, node);
        node.children.push(childNode);
        node.value += childNode.value;
        node.height = Math.max(node.height, childNode.height + 1);
    }
    return node;
}

function hierarchyDescendants(root) {
    const nodes = [];
    const stack = [root];
    while (stack.length > 0) {
        const node = stack.pop();
        nodes.push(node);
        for (let i = node.children.length - 1; i >= 0; i -= 1) {
            stack.push(node.children[i]);
        }
    }
    return nodes;
}

// d3's partition layout, turned on its side: a node spans [x0, x1) down the chart in
// proportion to its value, and the column [y0, y1) of its depth across it. A node's own count
// is the part of its span its children leave empty.
function partitionLayout(root, height, columnWidth) {
    root.x0 = 0;
    root.x1 = height;
    for (const node of hierarchyDescendants(root)) {
        node.y0 = node.depth * columnWidth;
        node.y1 = node.y0 + columnWidth;
        const scale = node.value > 0 ? (node.x1 - node.x0) / node.value : 0;
        let x = node.x0;
        for (const child of node.children) {
            child.x0 = x;
            x += child.value * scale;
            child.x1 = x;
        }
    }
    return root;
}

function hierarchyPath(node) {
    const names = [];
    for (let current = node; current; current = current.parent) {
        names.push(current.data.name);
    }
    return names.reverse().join(" / ");
}

// Evenly spaced hues, one per top-level branch.
function branchColors(count) {
    return Array.from({length: count}, (_, i) => `hsl(${Math.round(360 * i / Math.max(1, count))}, 62%, 52%)`);
}

// Click a cell to zoom to it, and click it again to zoom back out to its parent. Three depths
// are shown across the chart's width.
function drawIcicle(container, tree, {width = 928, height = 600, itemLabel = "items",
    fileStem = "icicle", idPrefix = ""} = {}) {
    const root = partitionLayout(hierarchyNodes(tree), height, width / 3);
    const nodes = hierarchyDescendants(root);
    const colors = branchColors(root.children.length);
    const colorOf = node => {
        if (node.depth === 0) {
            return "#cccccc";
        }
        while (node.depth > 1) {
            node = node.parent;
        }
        return colors[node.parent.children.indexOf(node)];
    };

    const svg = svgElement("svg", {
        viewBox: `0 0 ${width} ${height}`, role: "img", "aria-label": "Icicle chart",
        style: `width: 100%; max-width: ${width}px; height: auto;`,
        "font-family": CHART_FONT, "font-size": 10,
    });
    const cells = nodes.map(node => {
        const group = svgElement("g", {}, svg);
        const rect = svgElement("rect", {fill: colorOf(node), "fill-opacity": 0.6}, group);
        rect.dataset.tip = "";
        rect.style.cursor = "pointer";
        const label = svgElement("text", {x: 4, y: 13, fill: CHART_INK, "pointer-events": "none"}, group);
        const name = svgElement("tspan", {}, label);
        const value = svgElement("tspan", {"fill-opacity": 0.7}, label);
        value.textContent = ` ${formatChartNumber(node.value)}`;
        rect.addEventListener("click", () => zoom(node));
        return {node, group, rect, label, name, value};
    });
    const tips = new Map(cells.map(cell => [cell.rect, `${hierarchyPath(cell.node)}\n${formatChartNumber(cell.node.value)} ${itemLabel}`]));
    attachTooltip(svg, rect => tips.get(rect));

    // Where each node is drawn when `focus` fills the chart.
    const targets = focus => {
        const span = focus.x1 - focus.x0 || 1;
        return new Map(nodes.map(node => [node, {
            x0: (node.x0 - focus.x0) / span * height,
            x1: (node.x1 - focus.x0) / span * height,
            y0: node.y0 - focus.y0,
            y1: node.y1 - focus.y0,
        }]));
    };

    const place = positions => {
        for (const {node, group, rect, label, name} of cells) {
            const p = positions.get(node);
            const cellHeight = p.x1 - p.x0;
            const cellWidth = p.y1 - p.y0 - 1;
            group.setAttribute("transform", `translate(${p.y0},${p.x0})`);
            rect.setAttribute("width", Math.max(0, cellWidth));
            rect.setAttribute("height", Math.max(0, cellHeight - Math.min(1, cellHeight / 2)));
            const visible = p.y1 <= width + 1e-6 && p.y0 >= -1e-6 && cellHeight > 16;
            label.setAttribute("visibility", visible ? "visible" : "hidden");
            if (visible) {
                // About 6 pixels a character; the full name is in the tooltip.
                const room = Math.max(1, Math.floor((cellWidth - 8) / 6) - String(node.value).length - 1);
                const text = node.data.name;
                name.textContent = text.length > room ? `${text.slice(0, Math.max(1, room - 1))}…` : text;
            }
        }
    };

    let focus = root;
    let current = targets(root);
    let animation = null;
    function zoom(node) {
        focus = focus === node && node.parent ? node.parent : node;
        const from = current;
        const to = targets(focus);
        const start = performance.now();
        const duration = 600;
        cancelAnimationFrame(animation);
        const step = now => {
            const t = Math.min(1, (now - start) / duration);
            // Cubic ease in and out, as d3's transitions.
            const e = t < 0.5 ? 4 * t * t * t : 1 - (-2 * t + 2) ** 3 / 2;
            current = new Map(nodes.map(n => {
                const a = from.get(n);
                const b = to.get(n);
                return [n, {
                    x0: a.x0 + (b.x0 - a.x0) * e, x1: a.x1 + (b.x1 - a.x1) * e,
                    y0: a.y0 + (b.y0 - a.y0) * e, y1: a.y1 + (b.y1 - a.y1) * e,
                }];
            }));
            place(current);
            if (t < 1) {
                animation = requestAnimationFrame(step);
            }
        };
        animation = requestAnimationFrame(step);
    }

    place(current);
    chartFigure(container, svg, fileStem, idPrefix);
    return {svg, zoom: name => zoom(nodes.find(node => node.data.name === name) || root)};
}
