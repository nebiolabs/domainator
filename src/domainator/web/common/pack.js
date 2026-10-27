// Circles packed side by side, and a bubble chart drawn with them (built on charts.js).
//
// The packing is ported from d3-hierarchy's pack/siblings.js and pack/enclose.js (Copyright
// 2010-2021 Mike Bostock, ISC license).

// d3-hierarchy's linear congruential generator, so a packing is the same on every load.
function packRandom() {
    let state = 1;
    return () => (state = (1664525 * state + 1013904223) % 4294967296) / 4294967296;
}

function shuffleWith(array, random) {
    let m = array.length;
    while (m) {
        const i = random() * m-- | 0;
        [array[m], array[i]] = [array[i], array[m]];
    }
    return array;
}

function enclosesNot(a, b) {
    const dr = a.r - b.r;
    const dx = b.x - a.x;
    const dy = b.y - a.y;
    return dr < 0 || dr * dr < dx * dx + dy * dy;
}

function enclosesWeak(a, b) {
    const dr = a.r - b.r + Math.max(a.r, b.r, 1) * 1e-9;
    const dx = b.x - a.x;
    const dy = b.y - a.y;
    return dr > 0 && dr * dr > dx * dx + dy * dy;
}

function enclosesWeakAll(a, circles) {
    return circles.every(circle => enclosesWeak(a, circle));
}

function encloseBasis2(a, b) {
    const x21 = b.x - a.x;
    const y21 = b.y - a.y;
    const r21 = b.r - a.r;
    const l = Math.sqrt(x21 * x21 + y21 * y21);
    return {x: (a.x + b.x + x21 / l * r21) / 2, y: (a.y + b.y + y21 / l * r21) / 2, r: (l + a.r + b.r) / 2};
}

function encloseBasis3(a, b, c) {
    const {x: x1, y: y1, r: r1} = a;
    const {x: x2, y: y2, r: r2} = b;
    const {x: x3, y: y3, r: r3} = c;
    const a2 = x1 - x2;
    const a3 = x1 - x3;
    const b2 = y1 - y2;
    const b3 = y1 - y3;
    const c2 = r2 - r1;
    const c3 = r3 - r1;
    const d1 = x1 * x1 + y1 * y1 - r1 * r1;
    const d2 = d1 - x2 * x2 - y2 * y2 + r2 * r2;
    const d3 = d1 - x3 * x3 - y3 * y3 + r3 * r3;
    const ab = a3 * b2 - a2 * b3;
    const xa = (b2 * d3 - b3 * d2) / (ab * 2) - x1;
    const xb = (b3 * c2 - b2 * c3) / ab;
    const ya = (a3 * d2 - a2 * d3) / (ab * 2) - y1;
    const yb = (a2 * c3 - a3 * c2) / ab;
    const qa = xb * xb + yb * yb - 1;
    const qb = 2 * (r1 + xa * xb + ya * yb);
    const qc = xa * xa + ya * ya - r1 * r1;
    const r = -(Math.abs(qa) > 1e-6 ? (qb + Math.sqrt(qb * qb - 4 * qa * qc)) / (2 * qa) : qc / qb);
    return {x: x1 + xa + xb * r, y: y1 + ya + yb * r, r};
}

function encloseBasis(basis) {
    if (basis.length === 1) {
        return {x: basis[0].x, y: basis[0].y, r: basis[0].r};
    }
    return basis.length === 2 ? encloseBasis2(basis[0], basis[1]) : encloseBasis3(basis[0], basis[1], basis[2]);
}

function extendBasis(basis, p) {
    if (enclosesWeakAll(p, basis)) {
        return [p];
    }
    for (const b of basis) {
        if (enclosesNot(p, b) && enclosesWeakAll(encloseBasis2(b, p), basis)) {
            return [b, p];
        }
    }
    for (let i = 0; i < basis.length - 1; i += 1) {
        for (let j = i + 1; j < basis.length; j += 1) {
            if (enclosesNot(encloseBasis2(basis[i], basis[j]), p)
                && enclosesNot(encloseBasis2(basis[i], p), basis[j])
                && enclosesNot(encloseBasis2(basis[j], p), basis[i])
                && enclosesWeakAll(encloseBasis3(basis[i], basis[j], p), basis)) {
                return [basis[i], basis[j], p];
            }
        }
    }
    throw new Error("Could not enclose the circles.");
}

// The smallest circle {x, y, r} enclosing every circle {x, y, r}.
function packEnclose(circles, random) {
    const shuffled = shuffleWith(Array.from(circles), random);
    let basis = [];
    let enclosing = null;
    let i = 0;
    while (i < shuffled.length) {
        const p = shuffled[i];
        if (enclosing && enclosesWeak(enclosing, p)) {
            i += 1;
        } else {
            basis = extendBasis(basis, p);
            enclosing = encloseBasis(basis);
            i = 0;
        }
    }
    return enclosing;
}

// Places c tangent to both a and b.
function placeTangent(b, a, c) {
    const dx = b.x - a.x;
    const dy = b.y - a.y;
    const d2 = dx * dx + dy * dy;
    if (d2) {
        const a2 = (a.r + c.r) ** 2;
        const b2 = (b.r + c.r) ** 2;
        if (a2 > b2) {
            const x = (d2 + b2 - a2) / (2 * d2);
            const y = Math.sqrt(Math.max(0, b2 / d2 - x * x));
            c.x = b.x - x * dx - y * dy;
            c.y = b.y - x * dy + y * dx;
        } else {
            const x = (d2 + a2 - b2) / (2 * d2);
            const y = Math.sqrt(Math.max(0, a2 / d2 - x * x));
            c.x = a.x + x * dx - y * dy;
            c.y = a.y + x * dy + y * dx;
        }
    } else {
        c.x = a.x + c.r;
        c.y = a.y;
    }
}

function circlesIntersect(a, b) {
    const dr = a.r + b.r - 1e-6;
    const dx = b.x - a.x;
    const dy = b.y - a.y;
    return dr > 0 && dr * dr > dx * dx + dy * dy;
}

// How far from the origin the point between a front-chain node and the next one is.
function chainScore(node) {
    const a = node.circle;
    const b = node.next.circle;
    const ab = a.r + b.r;
    const dx = (a.x * b.r + b.x * a.r) / ab;
    const dy = (a.y * b.r + b.y * a.r) / ab;
    return dx * dx + dy * dy;
}

// Places circles {r} side by side around the origin, each tangent to two placed before it, by
// setting their x and y. Returns the radius of the circle enclosing them all, centered on the
// origin.
function packSiblings(circles, random) {
    const n = circles.length;
    if (n === 0) {
        return 0;
    }
    let a = circles[0];
    a.x = 0;
    a.y = 0;
    if (n === 1) {
        return a.r;
    }
    let b = circles[1];
    a.x = -b.r;
    b.x = a.r;
    b.y = 0;
    if (n === 2) {
        return a.r + b.r;
    }
    let c = circles[2];
    placeTangent(b, a, c);

    // The front chain: a ring of the outermost circles.
    a = {circle: a};
    b = {circle: b};
    c = {circle: c};
    a.next = c.previous = b;
    b.next = a.previous = c;
    c.next = b.previous = a;

    pack: for (let i = 3; i < n; i += 1) {
        placeTangent(a.circle, b.circle, circles[i]);
        c = {circle: circles[i]};
        // Find the closest circle on the chain that c intersects, if any, by distance along
        // the chain, ahead of b or behind a.
        let j = b.next;
        let k = a.previous;
        let sj = b.circle.r;
        let sk = a.circle.r;
        do {
            if (sj <= sk) {
                if (circlesIntersect(j.circle, c.circle)) {
                    b = j;
                    a.next = b;
                    b.previous = a;
                    i -= 1;
                    continue pack;
                }
                sj += j.circle.r;
                j = j.next;
            } else {
                if (circlesIntersect(k.circle, c.circle)) {
                    a = k;
                    a.next = b;
                    b.previous = a;
                    i -= 1;
                    continue pack;
                }
                sk += k.circle.r;
                k = k.previous;
            }
        } while (j !== k.next);

        // c fits between a and b.
        c.previous = a;
        c.next = b;
        a.next = b.previous = b = c;

        // The next circle goes beside the chain pair closest to the origin.
        let best = chainScore(a);
        while ((c = c.next) !== b) {
            const score = chainScore(c);
            if (score < best) {
                a = c;
                best = score;
            }
        }
        b = a.next;
    }

    const chain = [b.circle];
    for (let node = b.next; node !== b; node = node.next) {
        chain.push(node.circle);
    }
    const enclosing = packEnclose(chain, random);
    for (const circle of circles) {
        circle.x -= enclosing.x;
        circle.y -= enclosing.y;
    }
    return enclosing.r;
}

// d3's pack() of one level: circles with areas in proportion to values, packed into a circle
// that fills the smaller side of width x height, padding pixels apart. Returns [{value, index,
// x, y, r}] in the order of values.
function packValues(values, width, height, padding = 0) {
    const circles = values.map((value, index) => ({value, index, r: Math.sqrt(Math.max(0, value)), x: 0, y: 0}));
    // As d3: pack once to learn the scale, then again with the padding in those units.
    let radius = packSiblings(circles, packRandom());
    const size = Math.min(width, height);
    if (padding > 0 && radius > 0) {
        const pad = padding * radius / size;
        circles.forEach(circle => { circle.r += pad; });
        radius = packSiblings(circles, packRandom()) + pad;
        circles.forEach(circle => { circle.r -= pad; });
    }
    const k = radius > 0 ? size / (2 * radius) : 1;
    for (const circle of circles) {
        circle.x = width / 2 + circle.x * k;
        circle.y = height / 2 + circle.y * k;
        circle.r *= k;
    }
    return circles;
}

// ColorBrewer's YlOrRd (Cynthia Brewer, Apache 2.0), interpolated: t in [0, 1], light to dark.
const YL_OR_RD = ["#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026", "#800026"];

function sequentialColor(t, scheme = YL_OR_RD) {
    const position = Math.max(0, Math.min(1, t)) * (scheme.length - 1);
    const i = Math.min(scheme.length - 2, Math.floor(position));
    const f = position - i;
    const rgb = hex => [1, 3, 5].map(offset => parseInt(hex.slice(offset, offset + 2), 16));
    const [from, to] = [rgb(scheme[i]), rgb(scheme[i + 1])];
    return `rgb(${from.map((value, j) => Math.round(value + (to[j] - value) * f)).join(",")})`;
}

// A bubble chart: one circle per value, its area in proportion to it. Returns {update(values)},
// values being the numbers to draw, largest first; label(i) names circle i in its tooltip.
function createBubbleChart(container, {width = 760, height = 360, padding = 4, itemLabel = "items",
    label = i => `#${i + 1}`, emptyText = "Nothing to show.", fileStem = "bubbles", idPrefix = ""} = {}) {
    const svg = svgElement("svg", {
        viewBox: `0 0 ${width} ${height}`, role: "img", "aria-label": "Bubble chart",
        style: `width: 100%; max-width: ${width}px; height: auto;`,
        "font-family": CHART_FONT,
    });
    chartFigure(container, svg, fileStem, idPrefix);
    const bubbles = svgElement("g", {}, svg);
    attachTooltip(bubbles, circle => circle.dataset.tip);
    return {
        svg,
        update(values) {
            bubbles.replaceChildren();
            if (values.length === 0) {
                const text = svgElement("text", {x: width / 2, y: height / 2, "text-anchor": "middle", fill: CHART_MUTED, "font-size": 13}, bubbles);
                text.textContent = emptyText;
                return;
            }
            const top = Math.max(1, values[0], values[values.length - 1]);
            for (const circle of packValues(values, width, height, padding)) {
                const group = svgElement("g", {transform: `translate(${circle.x.toFixed(2)},${circle.y.toFixed(2)})`}, bubbles);
                const shape = svgElement("circle", {
                    r: circle.r.toFixed(2), fill: sequentialColor(circle.value / top), "fill-opacity": 0.85,
                    stroke: "#9a3412", "stroke-width": 1,
                }, group);
                shape.dataset.tip = `${label(circle.index)}: ${formatChartNumber(circle.value)} ${itemLabel}`;
                if (circle.r >= 16) {
                    const text = svgElement("text", {
                        "text-anchor": "middle", dy: "0.35em", "font-weight": 600, fill: "#3f1d0f",
                        "font-size": Math.max(9, Math.min(14, circle.r / 2.4)).toFixed(1), "pointer-events": "none",
                    }, group);
                    text.textContent = formatChartNumber(circle.value);
                }
            }
        },
    };
}
