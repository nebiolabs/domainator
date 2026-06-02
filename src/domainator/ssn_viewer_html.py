"""Static HTML shell for the standalone SSN viewer."""

import base64
import json
from pathlib import Path

from domainator.output_guardrails import make_temporary_output_path


def _layout_worker_js() -> str:
    """Standalone JS for the layout Web Worker (no f-string escaping needed here)."""
    return """
function seededUnit(componentId, salt) {
    const raw = Math.sin((componentId + 1) * 12.9898 + salt * 78.233) * 43758.5453;
    return raw - Math.floor(raw);
}
function componentRadiusForSize(size) {
    return Math.sqrt(Math.max(1, size) * 48 / Math.PI);
}
function componentLeafOrder(componentIds, hierarchyNodes) {
    return [...componentIds].sort((a, b) => {
        const an = hierarchyNodes[a], bn = hierarchyNodes[b];
        return an.leaf_start - bn.leaf_start || bn.size - an.size || a - b;
    });
}
function treeCenter(nodeIds, adjacency, hierarchyNodes) {
    if (nodeIds.length <= 2) { return nodeIds[0]; }
    const degree = new Map();
    nodeIds.forEach(id => degree.set(id, (adjacency.get(id) || []).length));
    let leaves = nodeIds.filter(id => (degree.get(id) || 0) <= 1), remaining = nodeIds.length;
    while (remaining > 2 && leaves.length > 0) {
        remaining -= leaves.length;
        const nextLeaves = [];
        leaves.forEach(leafId => {
            (adjacency.get(leafId) || []).forEach(nid => {
                if (!degree.has(nid)) { return; }
                degree.set(nid, degree.get(nid) - 1);
                if (degree.get(nid) === 1) { nextLeaves.push(nid); }
            });
            degree.delete(leafId);
        });
        leaves = nextLeaves;
    }
    const candidates = degree.size > 0 ? Array.from(degree.keys()) : nodeIds;
    candidates.sort((a, b) => hierarchyNodes[b].size - hierarchyNodes[a].size || a - b);
    return candidates[0];
}
function rootedTree(rootId, adjacency, hierarchyNodes) {
    const parent = new Map([[rootId, null]]), order = [rootId];
    for (let i = 0; i < order.length; i++) {
        const nodeId = order[i];
        (adjacency.get(nodeId) || []).forEach(nid => { if (!parent.has(nid)) { parent.set(nid, nodeId); order.push(nid); } });
    }
    const children = new Map();
    order.forEach(id => children.set(id, []));
    for (let i = 1; i < order.length; i++) { children.get(parent.get(order[i])).push(order[i]); }
    children.forEach(childIds => {
        childIds.sort((a, b) => {
            const an = hierarchyNodes[a], bn = hierarchyNodes[b];
            return an.leaf_start - bn.leaf_start || bn.size - an.size || a - b;
        });
    });
    return {parent, order, children};
}
function pointSegmentDistance(point, start, end) {
    const dx = end.x - start.x, dy = end.y - start.y, lsq = dx * dx + dy * dy;
    if (lsq < 1e-9) { return {distance: Math.hypot(point.x - start.x, point.y - start.y), t: 0, closestX: start.x, closestY: start.y}; }
    const t = Math.max(0, Math.min(1, ((point.x - start.x) * dx + (point.y - start.y) * dy) / lsq));
    const cx = start.x + dx * t, cy = start.y + dy * t;
    return {distance: Math.hypot(point.x - cx, point.y - cy), t, closestX: cx, closestY: cy};
}
function segmentOrientation(a, b, c) { return (b.y - a.y) * (c.x - b.x) - (b.x - a.x) * (c.y - b.y); }
function segmentsCross(a, b, c, d) {
    return segmentOrientation(a,b,c) * segmentOrientation(a,b,d) < 0 && segmentOrientation(c,d,a) * segmentOrientation(c,d,b) < 0;
}
function refineLayoutGeometry(items, edgePairs, options) {
    options = options || {};
    const refined = items.map(item => Object.assign({}, item));
    if (refined.length <= 1) { return refined; }
    const itemById = new Map(refined.map(item => [item.componentId, item]));
    const links = edgePairs.map(pair => Array.isArray(pair) ? {sourceId: pair[0], targetId: pair[1]} : pair)
        .filter(link => itemById.has(link.sourceId) && itemById.has(link.targetId));
    const bp = options.bubblePadding != null ? options.bubblePadding : 14;
    const ep = options.edgePadding != null ? options.edgePadding : 8;
    const oi = options.overlapIterations != null ? options.overlapIterations : 5;
    const ei = options.edgeIterations != null ? options.edgeIterations : 3;
    const ci = options.crossingIterations != null ? options.crossingIterations : 2;
    const mpc = 180000, menc = 140000, mcc = 90000;
    const pairChecks = refined.length * (refined.length - 1) / 2;
    function overlapPass(iters, strength) {
        if (pairChecks > mpc) { return; }
        for (let it = 0; it < iters; it++) {
            for (let li = 0; li < refined.length; li++) {
                const l = refined[li];
                for (let ri = li + 1; ri < refined.length; ri++) {
                    const r = refined[ri];
                    let dx = r.x - l.x, dy = r.y - l.y, dist = Math.hypot(dx, dy);
                    if (dist < 1e-6) {
                        dx = (seededUnit(l.componentId + r.componentId, it + 21) - 0.5) || 0.01;
                        dy = (seededUnit(l.componentId + r.componentId, it + 22) - 0.5) || 0.01;
                        dist = Math.hypot(dx, dy);
                    }
                    const md = l.radius + r.radius + bp;
                    if (dist >= md) { continue; }
                    const sh = (md - dist) / 2 * strength;
                    l.x -= dx / dist * sh; l.y -= dy / dist * sh;
                    r.x += dx / dist * sh; r.y += dy / dist * sh;
                }
            }
        }
    }
    function edgePass(iters, strength) {
        if (links.length * refined.length > menc) { return; }
        for (let it = 0; it < iters; it++) {
            links.forEach(link => {
                const s = itemById.get(link.sourceId), t = itemById.get(link.targetId);
                if (!s || !t) { return; }
                refined.forEach(item => {
                    if (item.componentId === link.sourceId || item.componentId === link.targetId) { return; }
                    const hit = pointSegmentDistance(item, s, t);
                    if (hit.t <= 0.03 || hit.t >= 0.97) { return; }
                    const md = item.radius + ep;
                    if (hit.distance >= md) { return; }
                    let nx = item.x - hit.closestX, ny = item.y - hit.closestY, nl = Math.hypot(nx, ny);
                    if (nl < 1e-6) { const ex = t.x - s.x, ey = t.y - s.y; nx = -ey || 1; ny = ex || 0; nl = Math.hypot(nx, ny); }
                    const push = (md - hit.distance) * strength;
                    item.x += nx / nl * push; item.y += ny / nl * push;
                });
            });
        }
    }
    overlapPass(oi, 0.72);
    edgePass(ei, 0.68);
    if (links.length * (links.length - 1) / 2 <= mcc) {
        for (let it = 0; it < ci; it++) {
            for (let li = 0; li < links.length; li++) {
                const la = links[li], a = itemById.get(la.sourceId), b = itemById.get(la.targetId);
                if (!a || !b) { continue; }
                for (let ri = li + 1; ri < links.length; ri++) {
                    const rb = links[ri];
                    if (la.sourceId === rb.sourceId || la.sourceId === rb.targetId || la.targetId === rb.sourceId || la.targetId === rb.targetId) { continue; }
                    const c = itemById.get(rb.sourceId), d = itemById.get(rb.targetId);
                    if (!c || !d || !segmentsCross(a, b, c, d)) { continue; }
                    const dx = b.x - a.x, dy = b.y - a.y, len = Math.hypot(dx, dy);
                    if (len < 1e-6) { continue; }
                    const nx = -dy / len, ny = dx / len, push = 4.5 + it * 1.5;
                    a.x += nx * push; a.y += ny * push; b.x += nx * push; b.y += ny * push;
                    c.x -= nx * push; c.y -= ny * push; d.x -= nx * push; d.y -= ny * push;
                }
            }
        }
    }
    edgePass(1, 0.82);
    overlapPass(4, 0.9);
    return refined;
}
function normalizeComponentLayout(items, padding) {
    padding = padding != null ? padding : 20;
    if (!items.length) { return {items: [], width: 0, height: 0}; }
    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    items.forEach(item => {
        minX = Math.min(minX, item.x - item.radius); minY = Math.min(minY, item.y - item.radius);
        maxX = Math.max(maxX, item.x + item.radius); maxY = Math.max(maxY, item.y + item.radius);
    });
    return {
        items: items.map(item => Object.assign({}, item, {x: item.x - minX + padding, y: item.y - minY + padding})),
        width: (maxX - minX) + padding * 2, height: (maxY - minY) + padding * 2,
    };
}
function tidyComponentLayout(componentIds, adjacency, hierarchyNodes, originX) {
    const rootId = treeCenter(componentIds, adjacency, hierarchyNodes);
    const tree = rootedTree(rootId, adjacency, hierarchyNodes);
    const depthByNode = new Map([[rootId, 0]]);
    tree.order.forEach(nodeId => { (tree.children.get(nodeId) || []).forEach(childId => { depthByNode.set(childId, (depthByNode.get(nodeId) || 0) + 1); }); });
    const radii = new Map(componentIds.map(id => [id, componentRadiusForSize(hierarchyNodes[id].size)]));
    const maxRadius = Math.max(...componentIds.map(id => radii.get(id) || 10), 10);
    const siblingGap = Math.max(20, maxRadius * 0.42);
    const levelGap = Math.max(168, maxRadius * 2.2 + 96);
    const subtreeSpan = new Map();
    [...tree.order].reverse().forEach(nodeId => {
        const childIds = tree.children.get(nodeId) || [];
        const nodeSpan = (radii.get(nodeId) || 12) * 2 + 18;
        if (childIds.length === 0) { subtreeSpan.set(nodeId, nodeSpan); return; }
        const childrenSpan = childIds.reduce((s, c) => s + subtreeSpan.get(c), 0) + Math.max(0, childIds.length - 1) * siblingGap;
        subtreeSpan.set(nodeId, Math.max(nodeSpan, childrenSpan));
    });
    const positionById = new Map();
    function placeNode(nodeId, topY) {
        const nodeSpan = subtreeSpan.get(nodeId) || 26;
        const childIds = tree.children.get(nodeId) || [];
        const x = originX + (depthByNode.get(nodeId) || 0) * levelGap;
        if (childIds.length === 0) { positionById.set(nodeId, {componentId: nodeId, x, y: topY + nodeSpan / 2, radius: radii.get(nodeId) || 10}); return; }
        const childrenSpan = childIds.reduce((s, c) => s + subtreeSpan.get(c), 0) + Math.max(0, childIds.length - 1) * siblingGap;
        let childCursor = topY + Math.max(0, (nodeSpan - childrenSpan) / 2);
        const childCenters = [];
        childIds.forEach(childId => { placeNode(childId, childCursor); childCenters.push(positionById.get(childId).y); childCursor += subtreeSpan.get(childId) + siblingGap; });
        const centerY = childCenters.reduce((s, v) => s + v, 0) / Math.max(1, childCenters.length);
        positionById.set(nodeId, {componentId: nodeId, x, y: centerY, radius: radii.get(nodeId) || 10});
    }
    placeNode(rootId, 0);
    return {positionById, order: tree.order};
}
function simulateComponentLayout(componentIds, adjacency, options, hierarchyNodes) {
    options = options || {};
    const radii = new Map(componentIds.map(id => [id, componentRadiusForSize(hierarchyNodes[id].size)]));
    const positions = new Map(), velocities = new Map(), anchors = new Map();
    const damping = options.damping != null ? options.damping : 0.8;
    const repulsion = options.repulsion != null ? options.repulsion : 6000;
    const spring = options.spring != null ? options.spring : 0.07;
    const gravity = options.gravity != null ? options.gravity : 0.008;
    const anchorStrength = options.anchorStrength != null ? options.anchorStrength : 0;
    const iterations = options.iterations != null ? options.iterations : 110;
    const preferredEdgeLength = options.preferredEdgeLength != null ? options.preferredEdgeLength : 120;
    const collisionStrength = options.collisionStrength != null ? options.collisionStrength : 0.28;
    const maxStep = options.maxStep != null ? options.maxStep : 48;
    const maxPhysicsNodes = options.maxPhysicsNodes != null ? options.maxPhysicsNodes : 500;
    // Seed positions and anchors from the crossing-free tidy tree layout, centered on the origin.
    const tidy = tidyComponentLayout(componentIds, adjacency, hierarchyNodes, 0);
    let centerX = 0, centerY = 0;
    componentIds.forEach(id => { const s = tidy.positionById.get(id) || {x: 0, y: 0}; centerX += s.x; centerY += s.y; });
    centerX /= Math.max(1, componentIds.length); centerY /= Math.max(1, componentIds.length);
    componentIds.forEach(id => {
        const s = tidy.positionById.get(id) || {x: 0, y: 0};
        const ax = s.x - centerX, ay = s.y - centerY;
        anchors.set(id, {x: ax, y: ay});
        positions.set(id, {x: ax, y: ay});
        velocities.set(id, {x: 0, y: 0});
    });
    const edgePairs = [];
    componentIds.forEach(id => { (adjacency.get(id) || []).forEach(nid => { if (id < nid) { edgePairs.push([id, nid]); } }); });
    if (componentIds.length <= maxPhysicsNodes) {
    for (let iter = 0; iter < iterations; iter++) {
        const forces = new Map(componentIds.map(id => [id, {x: 0, y: 0}]));
        for (let li = 0; li < componentIds.length; li++) {
            const lid = componentIds[li], lp = positions.get(lid);
            for (let ri = li + 1; ri < componentIds.length; ri++) {
                const rid = componentIds[ri], rp = positions.get(rid);
                let dx = rp.x - lp.x, dy = rp.y - lp.y, dsq = dx * dx + dy * dy;
                if (dsq < 1e-6) {
                    dx = (seededUnit(lid + rid, iter + 5) - 0.5) * 0.01;
                    dy = (seededUnit(lid + rid, iter + 6) - 0.5) * 0.01;
                    dsq = dx * dx + dy * dy;
                }
                const dist = Math.sqrt(dsq), rf = repulsion / Math.max(dsq, 16);
                const ov = (radii.get(lid) || 0) + (radii.get(rid) || 0) + 12 - dist;
                const cf = ov > 0 ? ov * collisionStrength : 0;
                const fx = dx / dist * (rf + cf), fy = dy / dist * (rf + cf);
                forces.get(lid).x -= fx; forces.get(lid).y -= fy;
                forces.get(rid).x += fx; forces.get(rid).y += fy;
            }
        }
        edgePairs.forEach(([lid, rid]) => {
            const lp = positions.get(lid), rp = positions.get(rid);
            let dx = rp.x - lp.x, dy = rp.y - lp.y, dist = Math.hypot(dx, dy);
            if (dist < 1e-6) { dist = 1e-6; dx = preferredEdgeLength; dy = 0; }
            const df = spring * (dist - preferredEdgeLength), fx = dx / dist * df, fy = dy / dist * df;
            forces.get(lid).x += fx; forces.get(lid).y += fy;
            forces.get(rid).x -= fx; forces.get(rid).y -= fy;
        });
        let totalMovement = 0;
        componentIds.forEach(id => {
            const pos = positions.get(id), vel = velocities.get(id), f = forces.get(id), anch = anchors.get(id);
            f.x += (anch.x - pos.x) * anchorStrength - pos.x * gravity;
            f.y += (anch.y - pos.y) * anchorStrength - pos.y * gravity;
            vel.x = (vel.x + f.x) * damping; vel.y = (vel.y + f.y) * damping;
            const speed = Math.hypot(vel.x, vel.y);
            if (speed > maxStep) { vel.x *= maxStep / speed; vel.y *= maxStep / speed; }
            pos.x += vel.x; pos.y += vel.y;
            totalMovement += Math.abs(vel.x) + Math.abs(vel.y);
        });
        if (iter > 12 && totalMovement / componentIds.length < 0.05) { break; }
    }
    }
    const rawItems = componentIds.map(id => { const pos = positions.get(id); return {componentId: id, x: pos.x, y: pos.y, radius: radii.get(id) || 10}; });
    return normalizeComponentLayout(refineLayoutGeometry(rawItems, edgePairs, {
        bubblePadding: options.bubblePadding != null ? options.bubblePadding : 14,
        edgePadding: options.edgePadding != null ? options.edgePadding : 8,
        overlapIterations: options.geometryIterations != null ? options.geometryIterations : 5,
        edgeIterations: options.edgeIterations != null ? options.edgeIterations : 3,
        crossingIterations: options.crossingIterations != null ? options.crossingIterations : 2,
    }), 24);
}
function clusterGraphComponents(visibleIds, links, hierarchyNodes, sortBySizeEnabled) {
    const adjacency = new Map();
    visibleIds.forEach(id => adjacency.set(id, []));
    links.forEach(link => {
        (adjacency.get(link.sourceId) || []).push(link.targetId);
        (adjacency.get(link.targetId) || []).push(link.sourceId);
    });
    const components = [], seen = new Set();
    visibleIds.forEach(id => {
        if (seen.has(id)) { return; }
        const stack = [id], comp = [];
        seen.add(id);
        while (stack.length > 0) { const cur = stack.pop(); comp.push(cur); (adjacency.get(cur) || []).forEach(nid => { if (!seen.has(nid)) { seen.add(nid); stack.push(nid); } }); }
        components.push(comp);
    });
    components.sort((a, b) => {
        const ac = a.reduce((s, id) => s + hierarchyNodes[id].size, 0), bc = b.reduce((s, id) => s + hierarchyNodes[id].size, 0);
        const as_ = Math.min(...a.map(id => hierarchyNodes[id].leaf_start)), bs_ = Math.min(...b.map(id => hierarchyNodes[id].leaf_start));
        return sortBySizeEnabled ? (bc - ac || as_ - bs_ || b.length - a.length) : (as_ - bs_ || b.length - a.length);
    });
    return {adjacency, components};
}
function packLayouts(componentLayouts, options) {
    options = options || {};
    const gapX = options.gapX != null ? options.gapX : 120, gapY = options.gapY != null ? options.gapY : 120;
    const op = options.outerPadding != null ? options.outerPadding : 72;
    let rw;
    if (options.rowTargetWidth != null) {
        rw = options.rowTargetWidth;
    } else {
        const valid = componentLayouts.filter(c => c && c.items.length > 0);
        const totalArea = valid.reduce((s, c) => s + c.width * c.height, 0);
        const widest = valid.reduce((m, c) => Math.max(m, c.width), 0);
        rw = Math.max(widest, Math.sqrt(totalArea) * 1.3) + op;
    }
    const packed = []; let cx = op, cy = op, rh = 0;
    componentLayouts.forEach(comp => {
        if (!comp || !comp.items.length) { return; }
        if (cx > op && cx + comp.width > rw) { cx = op; cy += rh + gapY; rh = 0; }
        comp.items.forEach(item => packed.push({componentId: item.componentId, x: item.x + cx, y: item.y + cy, radius: item.radius}));
        cx += comp.width + gapX; rh = Math.max(rh, comp.height);
    });
    return packed;
}
self.onmessage = function(event) {
    const msg = event.data;
    if (msg.type !== 'computeLayout') { return; }
    const {requestId, key, algorithm, visibleIds, links, hierarchyNodes, sortBySizeEnabled} = msg;
    const {adjacency, components} = clusterGraphComponents(visibleIds, links, hierarchyNodes, sortBySizeEnabled);
    const componentLayouts = components.map(ids => simulateComponentLayout(ids, adjacency, {
        repulsion: 4000, spring: 0.08, gravity: 0.006,
        damping: 0.82, preferredEdgeLength: 124, iterations: Math.min(70, 40 + ids.length),
        collisionStrength: 0.5, bubblePadding: 17, edgePadding: 10, geometryIterations: 7,
        edgeIterations: 4, crossingIterations: 4, anchorStrength: 0.12, maxPhysicsNodes: 500,
    }, hierarchyNodes));
    const layout = packLayouts(componentLayouts, {gapX: 36, gapY: 36, outerPadding: 72});
    self.postMessage({requestId, key, layout});
};
"""


def ssn_viewer_html(
    title: str = "Domainator SSN Viewer",
    embedded_bundle_json: bytes | None = None,
) -> str:
    escaped_title = title.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    embedded_bundle_base64 = None
    if embedded_bundle_json is not None:
        embedded_bundle_base64 = base64.b64encode(embedded_bundle_json).decode("ascii")
    layout_worker_code_json = json.dumps(_layout_worker_js())
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
    .toolbar button,
    .toolbar input,
    .toolbar select {{
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
    .toolbar button[aria-pressed="true"] {{
        background: rgba(200, 85, 61, 0.14);
        border-color: rgba(200, 85, 61, 0.55);
        color: #7b2f22;
    }}
    .toolbar button:disabled {{
        opacity: 0.45;
        cursor: not-allowed;
    }}
    .toolbar input,
    .toolbar select {{
        width: auto;
        min-width: 170px;
        flex: 1 1 200px;
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
    .metadata-pager {{
        display: flex;
        justify-content: space-between;
        align-items: center;
        gap: 12px;
        flex-wrap: wrap;
        margin-top: 12px;
        color: var(--muted);
        font-size: 0.9rem;
    }}
    .metadata-pager-status {{
        min-width: 0;
    }}
    .metadata-pager-controls {{
        display: inline-flex;
        align-items: center;
        gap: 10px;
        flex-wrap: wrap;
    }}
    .metadata-pager-controls label {{
        display: inline-flex;
        align-items: center;
        gap: 8px;
    }}
    .metadata-pager button {{
        width: auto;
        min-width: 0;
    }}
    table {{
        border-collapse: collapse;
        width: 100%;
        table-layout: fixed;
        font-size: 0.92rem;
    }}
    th, td {{
        padding: 8px 10px;
        border-bottom: 1px solid rgba(216, 206, 194, 0.65);
        text-align: left;
        vertical-align: top;
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
    }}
    thead th {{
        position: sticky;
        top: 0;
        background: #fffaf4;
        z-index: 1;
    }}
    th.metadata-header {{
        position: sticky;
        top: 0;
        padding: 0;
        background: #fffaf4;
    }}
    .metadata-header-cell {{
        display: flex;
        align-items: center;
        gap: 8px;
        min-width: 0;
        padding: 8px 14px 8px 10px;
    }}
    tbody tr {{
        transition: background 140ms ease;
    }}
    tbody tr:not(.metadata-empty-row) {{
        cursor: pointer;
    }}
    tbody tr:not(.metadata-empty-row):hover {{
        background: rgba(200, 85, 61, 0.08);
    }}
    .metadata-row-selected {{
        background: rgba(200, 85, 61, 0.14);
    }}
    .metadata-sort-button {{
        display: inline-flex;
        align-items: center;
        gap: 6px;
        min-width: 0;
        flex: 1 1 auto;
        padding: 0;
        border: 0;
        background: transparent;
        color: inherit;
        font: inherit;
        font-weight: 600;
        text-align: left;
        cursor: pointer;
    }}
    .metadata-sort-button:hover {{
        color: #7b2f22;
    }}
    .metadata-sort-indicator {{
        color: var(--muted);
        font-size: 0.82rem;
        letter-spacing: 0.01em;
        flex: 0 0 auto;
    }}
    .metadata-copy-button {{
        flex: 0 0 auto;
        border: 1px solid rgba(216, 206, 194, 0.9);
        border-radius: 999px;
        background: rgba(255, 255, 255, 0.9);
        color: var(--muted);
        padding: 4px 9px;
        font: inherit;
        font-size: 0.78rem;
        letter-spacing: 0.01em;
        cursor: pointer;
        transition: background 140ms ease, color 140ms ease, border-color 140ms ease;
    }}
    .metadata-copy-button:hover {{
        color: #7b2f22;
        border-color: rgba(200, 85, 61, 0.4);
        background: rgba(200, 85, 61, 0.10);
    }}
    .metadata-resize-handle {{
        position: absolute;
        top: 0;
        right: 0;
        width: 12px;
        height: 100%;
        cursor: col-resize;
        touch-action: none;
    }}
    .metadata-resize-handle::before {{
        content: "";
        position: absolute;
        top: 9px;
        bottom: 9px;
        left: 5px;
        width: 2px;
        border-radius: 999px;
        background: rgba(216, 206, 194, 0.95);
        transition: background 140ms ease;
    }}
    .metadata-resize-handle:hover::before {{
        background: rgba(200, 85, 61, 0.58);
    }}
    .metadata-cell-number {{
        text-align: right;
        font-variant-numeric: tabular-nums;
    }}
    .metadata-cell-null {{
        color: var(--muted);
        font-style: italic;
    }}
    .metadata-empty-row td {{
        color: var(--muted);
        font-style: italic;
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
                <div class="note">Wheel to zoom, drag the background to pan, Shift-drag a box to select clusters, and Ctrl-click a node to toggle it individually. Hold Ctrl while Shift-dragging a box to select individual nodes within it. Hold Alt while Shift-dragging to deselect instead (combine with Ctrl to deselect individual nodes). Click a cluster bubble to toggle every node inside it.</div>
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
                            <option value="grid" selected>Grid (no edges)</option>
                            <option value="packed">Packed (no edges)</option>
                        </select>
                    </div>
                    <div class="control">
                        <label for="min-cluster-size">Minimum cluster size</label>
                        <input id="min-cluster-size" type="number" min="1" value="5" step="1" />
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
                        <input id="render-cluster-bounds" type="checkbox" checked />
                        <label for="render-cluster-bounds">Render cluster bounds</label>
                    </div>
                    <div class="control checkbox">
                        <input id="render-nodes" type="checkbox" checked />
                        <label for="render-nodes">Render nodes</label>
                    </div>
                    <div class="control checkbox">
                        <input id="reduce-elongation" type="checkbox" />
                        <label for="reduce-elongation">Reduce subcluster elongation</label>
                    </div>
                    <div class="control checkbox">
                        <input id="leaf-pruning-only" type="checkbox" />
                        <label for="leaf-pruning-only">Minimum cluster size trims leaf clusters only</label>
                    </div>
                </div>
                <div class="toolbar">
                    <button id="sort-components-by-size" type="button" aria-pressed="true" disabled>Sort clusters by size: On</button>
                    <button id="focus-largest-cluster" type="button" disabled>Focus largest cluster</button>
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
                    <button id="metadata-select-nodes" type="button" disabled>Select nodes</button>
                    <button id="metadata-deselect-rows" type="button" disabled>Deselect rows</button>
                    <button id="metadata-reset-sort" type="button" disabled>Reset table sort</button>
                    <input id="metadata-filter" type="search" placeholder="Search node_id and metadata" disabled />
                    <select id="metadata-null-order" disabled>
                        <option value="last" selected>Nulls last</option>
                        <option value="first">Nulls first</option>
                    </select>
                </div>
                <div class="table-wrap">
                    <table id="metadata-table">
                        <colgroup></colgroup>
                        <thead></thead>
                        <tbody></tbody>
                    </table>
                </div>
                <div class="metadata-pager">
                    <div id="metadata-page-status" class="metadata-pager-status">Load a bundle to page through metadata.</div>
                    <div class="metadata-pager-controls">
                        <label for="metadata-rows-per-page">Rows per page
                            <select id="metadata-rows-per-page" disabled>
                                <option value="50">50</option>
                                <option value="100">100</option>
                                <option value="250" selected>250</option>
                                <option value="500">500</option>
                                <option value="1000">1000</option>
                                <option value="2500">2500</option>
                                <option value="5000">5000</option>
                                <option value="all">All rows</option>
                            </select>
                        </label>
                        <button id="metadata-prev-page" type="button" disabled>Previous page</button>
                        <button id="metadata-next-page" type="button" disabled>Next page</button>
                    </div>
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
        metadataColumnByName: new Map(),
        metadataColumnIndexByName: new Map(),
        metadataColorInfoByName: new Map(),
        metadataColumnWidths: new Map(),
        metadataSearchTextByNodeIndex: [],
        metadataPage: 0,
        activeClusters: [],
        visibleClusters: [],
        selectedNodeIndices: new Set(),
        selectedMetadataNodeIndices: new Set(),
        metadataRowSelectionAnchor: null,
        metadataResize: null,
        visibleLayout: [],
        splitLinks: [],
        sliderModel: null,
        metadataSort: {{columnKey: null, direction: 'asc'}},
        viewTransform: {{scale: 1, offsetX: 0, offsetY: 0, minScale: 0.02, maxScale: 10}},
        selectionBox: null,
        dragState: null,
        suppressClick: false,
        layoutCache: new Map(),
        dotLayoutCache: new Map(),
        pendingThresholdUIFrame: null,
        pendingThresholdUIResetView: false,
        pendingMetadataFilterTimer: null,
        pendingClusterRenderFrame: null,
        allNodeIndices: [],
        nodeColorCache: [],
        metadataBaseNoteText: '',
        renderedNodeIndices: [],
        layoutWorker: null,
        pendingLayoutRequestId: 0,
        _pendingVisibleGraph: null,
        _pendingLayoutAlgorithm: null,
        _pendingLayoutResetView: true,
        layoutComputing: false,
    }};

    const splitCanvas = document.getElementById('split-chart');
    const splitContext = splitCanvas.getContext('2d');
    const clusterCanvas = document.getElementById('cluster-view');
    const clusterContext = clusterCanvas.getContext('2d');
    const DEFAULT_METADATA_PAGE_SIZE = 250;
    const EMBEDDED_BUNDLE_BASE64 = {json.dumps(embedded_bundle_base64)};
    const LAYOUT_WORKER_CODE = {layout_worker_code_json};

    function setStatus(message) {{
        document.getElementById('bundle-status').textContent = message;
    }}

    async function writeTextToClipboard(text) {{
        if (navigator.clipboard?.writeText) {{
            await navigator.clipboard.writeText(text);
            return;
        }}
        const helper = document.createElement('textarea');
        helper.value = text;
        helper.setAttribute('readonly', 'readonly');
        helper.style.position = 'fixed';
        helper.style.opacity = '0';
        document.body.appendChild(helper);
        helper.select();
        document.execCommand('copy');
        document.body.removeChild(helper);
    }}


    function browserSupportsBundleLoading() {{
        return 'DecompressionStream' in window;
    }}

    function warnIfUnsupported() {{
        if (!browserSupportsBundleLoading()) {{
            document.getElementById('browser-warning').textContent = EMBEDDED_BUNDLE_BASE64
                ? 'This browser lacks DecompressionStream support; bundled data will still open, but loading external gzipped bundles is unavailable.'
                : 'This browser lacks DecompressionStream support; gzipped bundles cannot be loaded.';
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

    function base64ToBytes(base64Value) {{
        return Uint8Array.from(atob(base64Value), c => c.charCodeAt(0));
    }}

    async function autoloadEmbeddedBundle() {{
        if (!EMBEDDED_BUNDLE_BASE64) {{
            return;
        }}
        setStatus('Loading bundled data...');
        try {{
            const bytes = base64ToBytes(EMBEDDED_BUNDLE_BASE64);
            let text;
            if (browserSupportsBundleLoading()) {{
                const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));
                text = await new Response(stream).text();
            }} else {{
                text = new TextDecoder().decode(bytes);
            }}
            installBundle(JSON.parse(text));
        }} catch (error) {{
            console.error(error);
            setStatus('Failed to load bundled data: ' + error.message);
        }}
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
        state.metadataSort = {{columnKey: null, direction: 'asc'}};
        state.metadataColumnWidths = new Map();
        state.metadataPage = 0;
        state.metadataResize = null;
        state.selectedNodeIndices = new Set();
        state.selectedMetadataNodeIndices = new Set();
        state.metadataRowSelectionAnchor = null;
        state.dotLayoutCache = new Map();
        state.layoutCache = new Map();
        state.sliderModel = buildSliderModel(bundle.graph.slider_stops || []);
        state.allNodeIndices = bundle.graph.nodes.map((_, i) => i);
        rebuildMetadataCaches();
        rebuildNodeColorCache();

        populateMetadataControls();
        const slider = document.getElementById('threshold-slider');
        slider.max = String(state.sliderModel.maxPosition);
        slider.value = String(state.sliderModel.initialPosition);
        slider.disabled = state.sliderModel.stops.length === 0;
        document.getElementById('threshold-min-label').textContent = state.sliderModel.minLabel;
        document.getElementById('threshold-max-label').textContent = state.sliderModel.maxLabel;
        document.getElementById('sort-components-by-size').disabled = false;
        document.getElementById('focus-largest-cluster').disabled = false;
        document.getElementById('reset-view').disabled = false;
        document.getElementById('threshold-input').disabled = state.sliderModel.stops.length === 0;
        document.getElementById('metadata-select-nodes').disabled = true;
        document.getElementById('metadata-deselect-rows').disabled = true;
        document.getElementById('metadata-reset-sort').disabled = true;
        document.getElementById('metadata-filter').disabled = false;
        document.getElementById('metadata-filter').value = '';
        document.getElementById('metadata-null-order').disabled = false;
        document.getElementById('metadata-null-order').value = 'last';
        document.getElementById('metadata-rows-per-page').disabled = false;
        document.getElementById('metadata-prev-page').disabled = true;
        document.getElementById('metadata-next-page').disabled = true;
        setStatus(
            'Loaded ' + (bundle.name || 'bundle') + ' with ' + bundle.graph.nodes.length.toLocaleString() + ' nodes.'
        );
        updateThresholdUI();
        updateMetadataTable();
    }}

    function rebuildMetadataCaches() {{
        state.metadataColumnByName = new Map();
        state.metadataColumnIndexByName = new Map();
        state.metadataColorInfoByName = new Map();

        state.metadataColumns.forEach((column, columnIndex) => {{
            state.metadataColumnByName.set(column.name, column);
            state.metadataColumnIndexByName.set(column.name, columnIndex);
            if (column.type === 'int' || column.type === 'float') {{
                state.metadataColorInfoByName.set(column.name, {{type: 'numeric', min: Infinity, max: -Infinity}});
                return;
            }}
            state.metadataColorInfoByName.set(column.name, {{type: 'categorical'}});
        }});

        state.metadataByNodeIndex.forEach(row => {{
            if (!row) {{
                return;
            }}
            state.metadataColumns.forEach((column, columnIndex) => {{
                if (column.type !== 'int' && column.type !== 'float') {{
                    return;
                }}
                const value = row[columnIndex];
                if (typeof value !== 'number' || Number.isNaN(value)) {{
                    return;
                }}
                const info = state.metadataColorInfoByName.get(column.name);
                if (!info) {{
                    return;
                }}
                info.min = Math.min(info.min, value);
                info.max = Math.max(info.max, value);
            }});
        }});

        state.metadataColorInfoByName.forEach(info => {{
            if (info.type !== 'numeric') {{
                return;
            }}
            if (!Number.isFinite(info.min) || !Number.isFinite(info.max)) {{
                info.min = 0;
                info.max = 0;
            }}
        }});

        state.metadataSearchTextByNodeIndex = state.bundle.graph.nodes.map((nodeName, nodeIndex) => {{
            const row = state.metadataByNodeIndex[nodeIndex] || [];
            const parts = [String(nodeName)];
            row.forEach(value => {{
                if (value === null || value === undefined) {{
                    return;
                }}
                if (Array.isArray(value)) {{
                    parts.push(value.join(' '));
                    return;
                }}
                parts.push(String(value));
            }});
            return parts.join(' ').toLowerCase();
        }});
    }}

    function rebuildNodeColorCache() {{
        if (!state.bundle) {{
            state.nodeColorCache = [];
            return;
        }}
        const nodeCount = state.bundle.graph.nodes.length;
        const cache = new Array(nodeCount);
        for (let i = 0; i < nodeCount; i++) {{
            cache[i] = nodeColor(i);
        }}
        state.nodeColorCache = cache;
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
                stops: [{{...infinityStop, sliderPosition: 0}}],
                maxPosition: 0,
                initialPosition: 0,
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
            initialPosition: 0,
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
        return document.getElementById('layout-algorithm').value || 'grid';
    }}

    function compareVisibleClusterIds(leftId, rightId) {{
        const leftNode = state.bundle.graph.hierarchy.nodes[leftId];
        const rightNode = state.bundle.graph.hierarchy.nodes[rightId];
        if (sortComponentsBySizeEnabled()) {{
            return rightNode.size - leftNode.size || leftNode.leaf_start - rightNode.leaf_start || leftId - rightId;
        }}
        return leftNode.leaf_start - rightNode.leaf_start || rightNode.size - leftNode.size || leftId - rightId;
    }}

    function sortComponentsBySizeEnabled() {{
        return document.getElementById('sort-components-by-size').getAttribute('aria-pressed') === 'true';
    }}

    function updateComponentSortButton() {{
        const button = document.getElementById('sort-components-by-size');
        const enabled = sortComponentsBySizeEnabled();
        const layoutMode = currentLayoutAlgorithm();
        const gridMode = layoutMode === 'grid' || layoutMode === 'packed';
        button.textContent = gridMode
            ? (enabled ? 'Sort clusters by size: On' : 'Sort clusters by size: Off')
            : (enabled ? 'Sort components by size: On' : 'Sort components by size: Off');
        button.title = gridMode
            ? (enabled
                ? 'Visible clusters are ordered by cluster size before grid placement.'
                : 'Visible clusters keep their hierarchy-derived order before grid placement.')
            : (enabled
                ? 'Connected components are ordered by total visible node count.'
                : 'Connected components keep their hierarchy-derived order.');
    }}

    function showNodeCountsEnabled() {{
        return document.getElementById('show-node-counts').checked;
    }}

    function showEdgeScoresEnabled() {{
        return document.getElementById('show-edge-scores').checked;
    }}

    function renderClusterBoundsEnabled() {{
        return document.getElementById('render-cluster-bounds').checked;
    }}

    function renderNodesEnabled() {{
        return document.getElementById('render-nodes').checked;
    }}

    function reduceElongationEnabled() {{
        return document.getElementById('reduce-elongation').checked;
    }}

    function currentColorField() {{
        return document.getElementById('color-by').value;
    }}

    function currentLabelField() {{
        return document.getElementById('label-by').value;
    }}

    function metadataColumn(name) {{
        return state.metadataColumnByName.get(name) || null;
    }}

    function metadataColumnKeys() {{
        return ['node_id', ...state.metadataColumns.map(column => column.name)];
    }}

    function defaultMetadataColumnWidth(columnKey) {{
        if (columnKey === 'node_id') {{
            return 230;
        }}
        const column = metadataColumn(columnKey);
        if (!column) {{
            return 180;
        }}
        if (column.type === 'int' || column.type === 'float') {{
            return 140;
        }}
        if (column.type === 'bool' || column.type === 'boolean') {{
            return 120;
        }}
        return 190;
    }}

    function metadataColumnWidth(columnKey) {{
        return state.metadataColumnWidths.get(columnKey) || defaultMetadataColumnWidth(columnKey);
    }}

    function applyMetadataColumnWidths() {{
        const colgroup = document.querySelector('#metadata-table colgroup');
        if (!colgroup) {{
            return;
        }}
        colgroup.innerHTML = '';
        metadataColumnKeys().forEach(columnKey => {{
            const col = document.createElement('col');
            col.style.width = metadataColumnWidth(columnKey) + 'px';
            colgroup.appendChild(col);
        }});
    }}

    function updateMetadataColumnWidth(columnKey, nextWidth) {{
        const clampedWidth = Math.max(96, Math.min(640, Math.round(nextWidth)));
        state.metadataColumnWidths.set(columnKey, clampedWidth);
        applyMetadataColumnWidths();
    }}

    function startMetadataColumnResize(columnKey, event) {{
        event.preventDefault();
        event.stopPropagation();
        state.metadataResize = {{
            columnKey,
            startX: event.clientX,
            startWidth: metadataColumnWidth(columnKey),
        }};
        document.body.style.cursor = 'col-resize';
        document.body.style.userSelect = 'none';
    }}

    function metadataValue(nodeIndex, columnName) {{
        if (!columnName) {{
            return null;
        }}
        const columnIndex = state.metadataColumnIndexByName.get(columnName);
        if (columnIndex === undefined) {{
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
        return state.metadataColorInfoByName.get(columnName) || null;
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

    function treeCenter(nodeIds, adjacency, hierarchyNodes) {{
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
            const leftSize = hierarchyNodes[leftId].size;
            const rightSize = hierarchyNodes[rightId].size;
            return rightSize - leftSize || leftId - rightId;
        }});
        return candidates[0];
    }}

    function rootedTree(rootId, adjacency, hierarchyNodes) {{
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
                const leftNode = hierarchyNodes[leftId];
                const rightNode = hierarchyNodes[rightId];
                return leftNode.leaf_start - rightNode.leaf_start || rightNode.size - leftNode.size || leftId - rightId;
            }});
        }});

        return {{parent, order, children}};
    }}

    function componentRadiusForSize(size) {{
        // Bubble radius scales so its area is proportional to the node count (every node
        // is always drawn).
        const areaPerNode = 48;
        return Math.sqrt((Math.max(1, size) * areaPerNode) / Math.PI);
    }}

    function componentDotCount(componentSize) {{
        // Every node in the component is drawn as a dot.
        return componentSize;
    }}

    function componentDotRadius(componentSize, bubbleRadius) {{
        return 1.95;
    }}

    function componentDotOffset(sampleIndex, sampleCount, packingRadius) {{
        if (sampleCount <= 1 || packingRadius <= 0) {{
            return {{x: 0, y: 0}};
        }}

        if (sampleCount <= 6) {{
            const ringRadius = sampleCount === 2
                ? packingRadius * 0.58
                : sampleCount <= 4
                    ? packingRadius * 0.72
                    : packingRadius * 0.78;
            const angleOffset = sampleCount === 4 ? Math.PI / 4 : -Math.PI / 2;
            const angle = angleOffset + (sampleIndex * ((Math.PI * 2) / sampleCount));
            return {{
                x: Math.cos(angle) * ringRadius,
                y: Math.sin(angle) * ringRadius,
            }};
        }}

        const angle = sampleIndex * 2.399963229728653;
        const distance = Math.sqrt((sampleIndex + 0.5) / sampleCount) * packingRadius;
        return {{
            x: Math.cos(angle) * distance,
            y: Math.sin(angle) * distance,
        }};
    }}


    function sampledMembersForComponent(componentId, sampleCount) {{
        const hierarchy = state.bundle.graph.hierarchy;
        const component = hierarchy.nodes[componentId];
        const sampledMembers = [];
        for (let sampleIndex = 0; sampleIndex < sampleCount; sampleIndex++) {{
            const leafPosition = component.leaf_start + Math.floor(sampleIndex * component.leaf_count / sampleCount);
            sampledMembers.push({{
                nodeIndex: hierarchy.leaf_order[leafPosition],
                leafPosition,
            }});
        }}
        return sampledMembers;
    }}

    function dotLayoutCacheKey(componentId, sampleCount, arrangement, spacingKey = 'base') {{
        return arrangement + ':' + componentId + ':' + sampleCount + ':' + spacingKey;
    }}

    function radialDotPositions(sampleCount) {{
        const positions = [];
        for (let sampleIndex = 0; sampleIndex < sampleCount; sampleIndex++) {{
            const offset = componentDotOffset(sampleIndex, sampleCount, 1);
            positions.push({{x: offset.x, y: offset.y, sampleIndex}});
        }}
        return positions;
    }}

    function radialDotLayout(sampledMembers) {{
        const positions = radialDotPositions(sampledMembers.length);
        return sampledMembers.map((member, sampleIndex) => {{
            const position = positions[sampleIndex];
            return {{memberIndex: member.nodeIndex, x: position.x, y: position.y}};
        }});
    }}

    const phyllotaxisCoordCache = new Map();

    function phyllotaxisFor(n) {{
        // Flat coordinate arrays for the phyllotaxis (sunflower) positions of size n.
        // Positions for a given n are identical across components, so cache by n. The
        // grouped layout only reads these, never mutates them, so sharing is safe.
        let coords = phyllotaxisCoordCache.get(n);
        if (coords === undefined) {{
            const points = radialDotPositions(n);
            const posX = new Float64Array(n);
            const posY = new Float64Array(n);
            for (let i = 0; i < n; i++) {{
                posX[i] = points[i].x;
                posY[i] = points[i].y;
            }}
            coords = {{posX: posX, posY: posY}};
            phyllotaxisCoordCache.set(n, coords);
        }}
        return coords;
    }}

    function selectByCoord(idx, start, end, leftCount, posX, posY, useY) {{
        // In-place Hoare quickselect: rearrange idx[start, end) so the leftCount entries
        // with the smallest coordinate on the chosen axis occupy [start, start+leftCount).
        // Coordinates are read indirectly via the shared position arrays. Deterministic
        // (median-of-three pivot); no allocation, no closures.
        const target = start + leftCount;
        let lo = start;
        let hi = end;
        while (hi - lo > 1) {{
            const a = lo;
            const b = (lo + hi) >> 1;
            const c = hi - 1;
            const va = useY ? posY[idx[a]] : posX[idx[a]];
            const vb = useY ? posY[idx[b]] : posX[idx[b]];
            const vc = useY ? posY[idx[c]] : posX[idx[c]];
            let pivot;
            if (va < vb) {{
                pivot = vb < vc ? vb : (va < vc ? vc : va);
            }} else {{
                pivot = va < vc ? va : (vb < vc ? vc : vb);
            }}
            let i = lo;
            let j = hi - 1;
            while (i <= j) {{
                while ((useY ? posY[idx[i]] : posX[idx[i]]) < pivot) {{ i++; }}
                while ((useY ? posY[idx[j]] : posX[idx[j]]) > pivot) {{ j--; }}
                if (i <= j) {{
                    const t = idx[i];
                    idx[i] = idx[j];
                    idx[j] = t;
                    i++;
                    j--;
                }}
            }}
            if (target <= j) {{
                hi = j + 1;
            }} else if (target >= i) {{
                lo = i;
            }} else {{
                break;
            }}
        }}
    }}

    function regionLongerAxisIsY(idx, pStart, pEnd, posX, posY) {{
        let minX = Infinity;
        let maxX = -Infinity;
        let minY = Infinity;
        let maxY = -Infinity;
        for (let k = pStart; k < pEnd; k++) {{
            const pi = idx[k];
            const xx = posX[pi];
            const yy = posY[pi];
            if (xx < minX) {{ minX = xx; }}
            if (xx > maxX) {{ maxX = xx; }}
            if (yy < minY) {{ minY = yy; }}
            if (yy > maxY) {{ maxY = yy; }}
        }}
        return (maxY - minY) > (maxX - minX);
    }}

    function regionPrincipalAxis(idx, pStart, pEnd, posX, posY) {{
        // Principal (largest-variance) axis of the points in [pStart, pEnd), via the
        // eigenvector of the 2x2 covariance matrix. Cutting perpendicular to this axis
        // keeps each child region rounder than the bounding-box axis when the region is
        // diagonally elongated. Used only when "reduce elongation" is enabled.
        const m = pEnd - pStart;
        let mx = 0;
        let my = 0;
        for (let k = pStart; k < pEnd; k++) {{
            const pi = idx[k];
            mx += posX[pi];
            my += posY[pi];
        }}
        mx /= m;
        my /= m;
        let sxx = 0;
        let sxy = 0;
        let syy = 0;
        for (let k = pStart; k < pEnd; k++) {{
            const pi = idx[k];
            const dx = posX[pi] - mx;
            const dy = posY[pi] - my;
            sxx += dx * dx;
            sxy += dx * dy;
            syy += dy * dy;
        }}
        const tr = sxx + syy;
        const diff = sxx - syy;
        const lambda = (tr / 2) + Math.sqrt(Math.max(0, (diff * diff) / 4 + sxy * sxy));
        // Two candidate eigenvectors; use whichever is better-conditioned (larger norm).
        const ax1 = lambda - syy;
        const ay1 = sxy;
        const ax2 = sxy;
        const ay2 = lambda - sxx;
        let ax;
        let ay;
        if ((ax1 * ax1 + ay1 * ay1) >= (ax2 * ax2 + ay2 * ay2)) {{ ax = ax1; ay = ay1; }} else {{ ax = ax2; ay = ay2; }}
        const norm = Math.hypot(ax, ay);
        if (norm < 1e-12) {{ return {{ax: 1, ay: 0}}; }}   // isotropic; any axis is fine
        return {{ax: ax / norm, ay: ay / norm}};
    }}

    function selectByProjection(idx, start, end, leftCount, posX, posY, ax, ay) {{
        // Like selectByCoord but partitions by projection onto an arbitrary axis (ax, ay)
        // instead of a single coordinate. In-place Hoare quickselect, median-of-three.
        const target = start + leftCount;
        let lo = start;
        let hi = end;
        while (hi - lo > 1) {{
            const a = lo;
            const b = (lo + hi) >> 1;
            const c = hi - 1;
            const va = posX[idx[a]] * ax + posY[idx[a]] * ay;
            const vb = posX[idx[b]] * ax + posY[idx[b]] * ay;
            const vc = posX[idx[c]] * ax + posY[idx[c]] * ay;
            let pivot;
            if (va < vb) {{
                pivot = vb < vc ? vb : (va < vc ? vc : va);
            }} else {{
                pivot = va < vc ? va : (vb < vc ? vc : vb);
            }}
            let i = lo;
            let j = hi - 1;
            while (i <= j) {{
                while ((posX[idx[i]] * ax + posY[idx[i]] * ay) < pivot) {{ i++; }}
                while ((posX[idx[j]] * ax + posY[idx[j]] * ay) > pivot) {{ j--; }}
                if (i <= j) {{
                    const t = idx[i];
                    idx[i] = idx[j];
                    idx[j] = t;
                    i++;
                    j--;
                }}
            }}
            if (target <= j) {{
                hi = j + 1;
            }} else if (target >= i) {{
                lo = i;
            }} else {{
                break;
            }}
        }}
    }}

    function collectMajorChildren(nodes, sampledMembers, nodeId, lo, hi, frac) {{
        // Contract the lopsided-split chain rooted at nodeId into an ordered list of
        // "major children" that exactly tile [lo, hi) in increasing leaf_start order.
        // Descend the big-child chain: peel each small sibling as a major child, stop at
        // the first balanced split (emit both children) or a leaf. This keeps every cut
        // on a real sibling boundary (nothing straddled) while collapsing a deep chain
        // into one multiway node so the layout stays O(n log n). Flat loop, no recursion.
        const head = [];
        const tail = [];
        let cur = nodeId;
        let curLo = lo;
        let curHi = hi;
        while (true) {{
            const node = cur >= 0 ? nodes[cur] : null;
            const m = curHi - curLo;
            if (!node || node.kind === 'leaf' || m <= 1) {{
                head.push({{node: cur, lo: curLo, hi: curHi}});
                break;
            }}
            const leftNode = nodes[node.left];
            const leftBoundary = leftNode.leaf_start + leftNode.leaf_count;
            // lowerBound in [curLo, curHi): first member whose leafPosition >= leftBoundary.
            let sLo = curLo;
            let sHi = curHi;
            while (sLo < sHi) {{
                const mid = (sLo + sHi) >> 1;
                if (sampledMembers[mid].leafPosition >= leftBoundary) {{ sHi = mid; }} else {{ sLo = mid + 1; }}
            }}
            const splitMember = sLo;
            const leftCount = splitMember - curLo;
            const rightCount = curHi - splitMember;
            if (leftCount <= 0 || rightCount <= 0) {{
                // Sampling can't resolve this boundary; keep cur as one atomic child.
                head.push({{node: cur, lo: curLo, hi: curHi}});
                break;
            }}
            if (Math.min(leftCount, rightCount) >= frac * m) {{
                // Balanced split: emit both children (left precedes right) and stop.
                head.push({{node: node.left, lo: curLo, hi: splitMember}});
                head.push({{node: node.right, lo: splitMember, hi: curHi}});
                break;
            }}
            // Lopsided: peel the smaller sibling, descend into the bigger.
            if (leftCount <= rightCount) {{
                head.push({{node: node.left, lo: curLo, hi: splitMember}});   // left precedes -> head
                cur = node.right; curLo = splitMember;
            }} else {{
                tail.push({{node: node.right, lo: splitMember, hi: curHi}});   // right follows -> tail
                cur = node.left; curHi = splitMember;
            }}
        }}
        // tail holds right-peels in decreasing leaf_start; reverse to increasing and append.
        for (let k = tail.length - 1; k >= 0; k--) {{ head.push(tail[k]); }}
        return head;
    }}

    function partitionAmongChildren(children, idx, posX, posY, pStart0, pEnd0, superStack, usePca) {{
        // Partition the point slice [pStart0, pEnd0) among the ordered child list by
        // recursively bisecting the LIST near its member-count median. Every list cut
        // falls between two children (a real sibling boundary), so no subtree is split.
        // Points are split by the longer bounding-box axis, or (usePca) by the region's
        // principal axis to reduce elongation. Emits a super-node frame per child. Own
        // explicit stack.
        const partStack = [{{ci: 0, cj: children.length, pStart: pStart0, pEnd: pEnd0}}];
        while (partStack.length > 0) {{
            const f = partStack.pop();
            const ci = f.ci;
            const cj = f.cj;
            const ps = f.pStart;
            const pe = f.pEnd;
            if (cj - ci === 1) {{
                const ch = children[ci];
                superStack.push({{lo: ch.lo, hi: ch.hi, pStart: ps, pEnd: pe, nodeId: ch.node}});
                continue;
            }}
            const mlo = children[ci].lo;
            const mhi = children[cj - 1].hi;
            const target = mlo + ((mhi - mlo) >> 1);
            // First list index in (ci, cj) whose child.lo >= target (balanced member split).
            let aLo = ci + 1;
            let aHi = cj;
            while (aLo < aHi) {{
                const mid = (aLo + aHi) >> 1;
                if (children[mid].lo >= target) {{ aHi = mid; }} else {{ aLo = mid + 1; }}
            }}
            let s = aLo;
            if (s <= ci) {{ s = ci + 1; }}
            if (s >= cj) {{ s = cj - 1; }}
            const leftMembers = children[s].lo - mlo;
            const rightMembers = mhi - children[s].lo;
            if (usePca) {{
                const axis = regionPrincipalAxis(idx, ps, pe, posX, posY);
                selectByProjection(idx, ps, pe, leftMembers, posX, posY, axis.ax, axis.ay);
            }} else {{
                const useY = regionLongerAxisIsY(idx, ps, pe, posX, posY);
                selectByCoord(idx, ps, pe, leftMembers, posX, posY, useY);
            }}
            const pSplit = ps + leftMembers;
            // Push the larger group first so the smaller is processed first (bounded stack).
            if (leftMembers >= rightMembers) {{
                partStack.push({{ci: ci, cj: s, pStart: ps, pEnd: pSplit}});
                partStack.push({{ci: s, cj: cj, pStart: pSplit, pEnd: pe}});
            }} else {{
                partStack.push({{ci: s, cj: cj, pStart: pSplit, pEnd: pe}});
                partStack.push({{ci: ci, cj: s, pStart: ps, pEnd: pSplit}});
            }}
        }}
    }}

    function groupedDotLayout(componentId, sampledMembers, usePca) {{
        // Assign each member to a phyllotaxis position (uniform spacing, from
        // radialDotPositions) so that each hierarchy subtree occupies a compact,
        // contiguous blob. Fully faithful: every cut lands on a true sibling boundary, so
        // no subtree -- at any scale -- is ever split across a partition. Deep lopsided
        // "caterpillar" chains are contracted into multiway super-nodes (collectMajorChildren)
        // so cost stays O(n log n) without ever resorting to a hierarchy-independent cut.
        // Iterative (explicit stacks) so deep trees cannot overflow the call stack.
        // usePca: cut along each region's principal axis instead of its bounding-box axis,
        // which reduces elongation of subcluster blobs at some extra per-node cost.
        const n = sampledMembers.length;
        if (n <= 6) {{
            // Keep the tuned ring placement for tiny clusters.
            return radialDotLayout(sampledMembers);
        }}
        const nodes = state.bundle.graph.hierarchy.nodes;
        const phy = phyllotaxisFor(n);
        const posX = phy.posX;
        const posY = phy.posY;
        const FRAC = 0.18;

        const idx = new Int32Array(n);
        for (let i = 0; i < n; i++) {{ idx[i] = i; }}
        const out = new Array(n);

        // Super-node frame invariant: pEnd - pStart === hi - lo (members paired one-to-one
        // with the points allocated to this subtree). idx is partitioned in place.
        const superStack = [{{lo: 0, hi: n, pStart: 0, pEnd: n, nodeId: componentId}}];
        while (superStack.length > 0) {{
            const fr = superStack.pop();
            const lo = fr.lo;
            const hi = fr.hi;
            const pStart = fr.pStart;
            const pEnd = fr.pEnd;
            const m = hi - lo;
            if (m <= 0) {{ continue; }}
            if (m === 1) {{
                const p = idx[pStart];
                out[lo] = {{memberIndex: sampledMembers[lo].nodeIndex, x: posX[p], y: posY[p]}};
                continue;
            }}
            const children = collectMajorChildren(nodes, sampledMembers, fr.nodeId, lo, hi, FRAC);
            if (children.length <= 1) {{
                // Atomic block (leaf with m>1, or unresolved sampling): assign in order.
                for (let k = 0; k < m; k++) {{
                    const p = idx[pStart + k];
                    out[lo + k] = {{memberIndex: sampledMembers[lo + k].nodeIndex, x: posX[p], y: posY[p]}};
                }}
                continue;
            }}
            partitionAmongChildren(children, idx, posX, posY, pStart, pEnd, superStack, usePca);
        }}
        return out;
    }}

    function normalizedComponentDotLayout(componentId, sampleCount, minimumDistance = 0) {{
        // Dots are always laid out with the grouped (chain-contraction) algorithm; the
        // arrangement dropdown was removed. The optional PCA axis reduces elongation.
        const usePca = reduceElongationEnabled();
        const spacingKey = usePca ? 'grouped-pca' : 'grouped-phyllotaxis';
        const cacheKey = dotLayoutCacheKey(componentId, sampleCount, 'grouped', spacingKey);
        const cached = state.dotLayoutCache.get(cacheKey);
        if (cached) {{
            return cached;
        }}

        const sampledMembers = sampledMembersForComponent(componentId, sampleCount);
        const layout = groupedDotLayout(componentId, sampledMembers, usePca);
        state.dotLayoutCache.set(cacheKey, layout);
        return layout;
    }}

    function clusterGraphComponents(visibleIds, links, hierarchyNodes, sortBySizeEnabled) {{
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
            const leftNodeCount = leftIds.reduce((sum, nodeId) => sum + hierarchyNodes[nodeId].size, 0);
            const rightNodeCount = rightIds.reduce((sum, nodeId) => sum + hierarchyNodes[nodeId].size, 0);
            const leftStart = Math.min(...leftIds.map(nodeId => hierarchyNodes[nodeId].leaf_start));
            const rightStart = Math.min(...rightIds.map(nodeId => hierarchyNodes[nodeId].leaf_start));
            if (sortBySizeEnabled) {{
                return rightNodeCount - leftNodeCount || leftStart - rightStart || rightIds.length - leftIds.length;
            }}
            return leftStart - rightStart || rightIds.length - leftIds.length;
        }});
        return {{adjacency, components}};
    }}

    function packLayouts(componentLayouts, options = {{}}) {{
        const gapX = options.gapX ?? 120;
        const gapY = options.gapY ?? 120;
        const outerPadding = options.outerPadding ?? 72;
        const valid = componentLayouts.filter(component => component && component.items.length > 0);
        // Adaptive near-square arrangement when no explicit width is given: keeps many disconnected
        // single-cluster components from spreading into a wide sparse grid (which would make
        // fit-to-view collapse to invisible specks).
        let rowTargetWidth = options.rowTargetWidth;
        if (rowTargetWidth == null) {{
            const totalArea = valid.reduce((sum, component) => sum + (component.width * component.height), 0);
            const widest = valid.reduce((maxWidth, component) => Math.max(maxWidth, component.width), 0);
            rowTargetWidth = Math.max(widest, Math.sqrt(totalArea) * 1.3) + outerPadding;
        }}
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

    function gridClusterLayout(visibleIds) {{
        if (visibleIds.length === 0) {{
            clusterCanvas.height = 760;
            return [];
        }}

        clusterCanvas.height = Math.max(760, Math.min(1180, Math.round(window.innerHeight * 0.8)));
        const orderedIds = [...visibleIds].sort(compareVisibleClusterIds);
        const items = orderedIds.map(componentId => {{
            const node = state.bundle.graph.hierarchy.nodes[componentId];
            return {{
                componentId,
                radius: componentRadiusForSize(node.size),
            }};
        }});
        const radiiDesc = items.map(item => item.radius).sort((a, b) => b - a);
        const r1 = radiiDesc[0] || 11;
        const r2 = radiiDesc.length > 1 ? radiiDesc[1] : r1;
        // Center-to-center = 0.5 * (largest diameter + second-largest diameter) = r1 + r2.
        // This is the tightest square-grid spacing that still guarantees no overlap (the only
        // pair that can touch is the single largest beside the second-largest). +12 visual gap.
        const cellSize = Math.max(22, r1 + r2 + 12);
        const outerPadding = 58;
        const aspectRatio = clusterCanvas.width / Math.max(clusterCanvas.height, 1);
        const columnCount = Math.max(1, Math.ceil(Math.sqrt(items.length * Math.max(0.75, aspectRatio))));

        return items.map((item, index) => {{
            const column = index % columnCount;
            const row = Math.floor(index / columnCount);
            return {{
                componentId: item.componentId,
                radius: item.radius,
                x: outerPadding + (column * cellSize) + (cellSize / 2),
                y: outerPadding + (row * cellSize) + (cellSize / 2),
            }};
        }});
    }}

    // Place circle c externally tangent to already-placed circles a and b (Wang front-chain).
    // Parameter order (b, a, c) mirrors d3's place() so the front-chain insertion stays overlap-free.
    function packPlace(b, a, c) {{
        const dx = b.x - a.x;
        const dy = b.y - a.y;
        const d2 = (dx * dx) + (dy * dy);
        if (d2 > 1e-9) {{
            const a2 = (a.radius + c.radius) * (a.radius + c.radius);
            const b2 = (b.radius + c.radius) * (b.radius + c.radius);
            if (a2 > b2) {{
                const x = (d2 + b2 - a2) / (2 * d2);
                const y = Math.sqrt(Math.max(0, (b2 / d2) - (x * x)));
                c.x = b.x - (x * dx) - (y * dy);
                c.y = b.y - (x * dy) + (y * dx);
            }} else {{
                const x = (d2 + a2 - b2) / (2 * d2);
                const y = Math.sqrt(Math.max(0, (a2 / d2) - (x * x)));
                c.x = a.x + (x * dx) - (y * dy);
                c.y = a.y + (x * dy) + (y * dx);
            }}
        }} else {{
            c.x = a.x + c.radius;
            c.y = a.y;
        }}
    }}

    function packIntersects(a, b) {{
        const dr = a.radius + b.radius - 1e-6;
        const dx = b.x - a.x;
        const dy = b.y - a.y;
        return dr > 0 && (dr * dr) > ((dx * dx) + (dy * dy));
    }}

    function packScore(node) {{
        const a = node._;
        const b = node.next._;
        const ab = a.radius + b.radius;
        const dx = ((a.x * b.radius) + (b.x * a.radius)) / ab;
        const dy = ((a.y * b.radius) + (b.y * a.radius)) / ab;
        return (dx * dx) + (dy * dy);
    }}

    // Front-chain circle packing (the algorithm behind d3's packSiblings). Packs circles in the
    // order given, compactly from the center outward, with no overlaps. Assigns x,y in place.
    function packSiblingsTight(circles) {{
        const n = circles.length;
        if (n === 0) {{ return; }}
        let a = circles[0];
        a.x = 0; a.y = 0;
        if (n === 1) {{ return; }}
        let b = circles[1];
        a.x = -b.radius; b.x = a.radius; b.y = 0;
        if (n === 2) {{ return; }}
        let c = circles[2];
        packPlace(b, a, c);

        let na = {{_: a, next: null, previous: null}};
        let nb = {{_: b, next: null, previous: null}};
        let nc = {{_: c, next: null, previous: null}};
        na.next = nc.previous = nb;
        nb.next = na.previous = nc;
        nc.next = nb.previous = na;

        pack: for (let i = 3; i < n; i++) {{
            c = circles[i];
            packPlace(na._, nb._, c);
            nc = {{_: c, next: null, previous: null}};

            let j = nb.next;
            let k = na.previous;
            let sj = nb._.radius;
            let sk = na._.radius;
            do {{
                if (sj <= sk) {{
                    if (packIntersects(j._, c)) {{
                        nb = j; na.next = nb; nb.previous = na; i--;
                        continue pack;
                    }}
                    sj += j._.radius; j = j.next;
                }} else {{
                    if (packIntersects(k._, c)) {{
                        na = k; na.next = nb; nb.previous = na; i--;
                        continue pack;
                    }}
                    sk += k._.radius; k = k.previous;
                }}
            }} while (j !== k.next);

            nc.previous = na; nc.next = nb; na.next = nb.previous = nb = nc;

            let aa = packScore(na);
            let nn = nc;
            while ((nn = nn.next) !== nb) {{
                const ca = packScore(nn);
                if (ca < aa) {{ na = nn; aa = ca; }}
            }}
            nb = na.next;
        }}
    }}

    function packedClusterLayout(visibleIds) {{
        if (visibleIds.length === 0) {{
            clusterCanvas.height = 760;
            return [];
        }}
        clusterCanvas.height = Math.max(760, Math.min(1180, Math.round(window.innerHeight * 0.8)));
        const orderedIds = [...visibleIds].sort(compareVisibleClusterIds);
        const circles = orderedIds.map(componentId => ({{
            componentId,
            radius: componentRadiusForSize(state.bundle.graph.hierarchy.nodes[componentId].size),
            x: 0,
            y: 0,
        }}));
        packSiblingsTight(circles);
        return normalizeComponentLayout(circles, 40).items;
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

    function componentLeafOrder(componentIds, hierarchyNodes) {{
        return [...componentIds].sort((leftId, rightId) => {{
            const leftNode = hierarchyNodes[leftId];
            const rightNode = hierarchyNodes[rightId];
            return leftNode.leaf_start - rightNode.leaf_start || rightNode.size - leftNode.size || leftId - rightId;
        }});
    }}

    function layoutLinkPairs(linksOrPairs) {{
        return linksOrPairs.map(link => {{
            if (Array.isArray(link)) {{
                return {{sourceId: link[0], targetId: link[1]}};
            }}
            return {{sourceId: link.sourceId, targetId: link.targetId}};
        }});
    }}

    function pointSegmentDistance(point, start, end) {{
        const dx = end.x - start.x;
        const dy = end.y - start.y;
        const lengthSq = (dx * dx) + (dy * dy);
        if (lengthSq < 1e-9) {{
            return {{distance: Math.hypot(point.x - start.x, point.y - start.y), t: 0, closestX: start.x, closestY: start.y}};
        }}
        const t = Math.max(0, Math.min(1, (((point.x - start.x) * dx) + ((point.y - start.y) * dy)) / lengthSq));
        const closestX = start.x + (dx * t);
        const closestY = start.y + (dy * t);
        return {{distance: Math.hypot(point.x - closestX, point.y - closestY), t, closestX, closestY}};
    }}

    function segmentOrientation(a, b, c) {{
        return ((b.y - a.y) * (c.x - b.x)) - ((b.x - a.x) * (c.y - b.y));
    }}

    function segmentsCross(a, b, c, d) {{
        const o1 = segmentOrientation(a, b, c);
        const o2 = segmentOrientation(a, b, d);
        const o3 = segmentOrientation(c, d, a);
        const o4 = segmentOrientation(c, d, b);
        return (o1 * o2 < 0) && (o3 * o4 < 0);
    }}

    function refineLayoutGeometry(items, linksOrPairs, options = {{}}) {{
        const refined = items.map(item => ({{...item}}));
        if (refined.length <= 1) {{
            return refined;
        }}

        const itemById = new Map(refined.map(item => [item.componentId, item]));
        const links = layoutLinkPairs(linksOrPairs).filter(link => itemById.has(link.sourceId) && itemById.has(link.targetId));
        const bubblePadding = options.bubblePadding ?? 14;
        const edgePadding = options.edgePadding ?? 8;
        const overlapIterations = options.overlapIterations ?? 5;
        const edgeIterations = options.edgeIterations ?? 3;
        const crossingIterations = options.crossingIterations ?? 2;
        const maxPairChecks = options.maxPairChecks ?? 180000;
        const maxEdgeNodeChecks = options.maxEdgeNodeChecks ?? 140000;
        const maxCrossingChecks = options.maxCrossingChecks ?? 90000;
        const pairChecks = (refined.length * (refined.length - 1)) / 2;

        if (pairChecks <= maxPairChecks) {{
            for (let iteration = 0; iteration < overlapIterations; iteration++) {{
                for (let leftIndex = 0; leftIndex < refined.length; leftIndex++) {{
                    const left = refined[leftIndex];
                    for (let rightIndex = leftIndex + 1; rightIndex < refined.length; rightIndex++) {{
                        const right = refined[rightIndex];
                        let dx = right.x - left.x;
                        let dy = right.y - left.y;
                        let distance = Math.hypot(dx, dy);
                        if (distance < 1e-6) {{
                            dx = (seededUnit(left.componentId + right.componentId, iteration + 21) - 0.5) || 0.01;
                            dy = (seededUnit(left.componentId + right.componentId, iteration + 22) - 0.5) || 0.01;
                            distance = Math.hypot(dx, dy);
                        }}
                        const minimumDistance = left.radius + right.radius + bubblePadding;
                        if (distance >= minimumDistance) {{
                            continue;
                        }}
                        const shift = ((minimumDistance - distance) / 2) * 0.72;
                        const shiftX = (dx / distance) * shift;
                        const shiftY = (dy / distance) * shift;
                        left.x -= shiftX;
                        left.y -= shiftY;
                        right.x += shiftX;
                        right.y += shiftY;
                    }}
                }}
            }}
        }}

        if (links.length * refined.length <= maxEdgeNodeChecks) {{
            for (let iteration = 0; iteration < edgeIterations; iteration++) {{
                links.forEach(link => {{
                    const source = itemById.get(link.sourceId);
                    const target = itemById.get(link.targetId);
                    if (!source || !target) {{
                        return;
                    }}
                    refined.forEach(item => {{
                        if (item.componentId === link.sourceId || item.componentId === link.targetId) {{
                            return;
                        }}
                        const hit = pointSegmentDistance(item, source, target);
                        if (hit.t <= 0.03 || hit.t >= 0.97) {{
                            return;
                        }}
                        const minimumDistance = item.radius + edgePadding;
                        if (hit.distance >= minimumDistance) {{
                            return;
                        }}
                        let normalX = item.x - hit.closestX;
                        let normalY = item.y - hit.closestY;
                        let normalLength = Math.hypot(normalX, normalY);
                        if (normalLength < 1e-6) {{
                            const edgeDx = target.x - source.x;
                            const edgeDy = target.y - source.y;
                            normalX = -edgeDy || 1;
                            normalY = edgeDx || 0;
                            normalLength = Math.hypot(normalX, normalY);
                        }}
                        const push = (minimumDistance - hit.distance) * 0.68;
                        item.x += (normalX / normalLength) * push;
                        item.y += (normalY / normalLength) * push;
                    }});
                }});
            }}
        }}

        const crossingChecks = (links.length * (links.length - 1)) / 2;
        if (crossingChecks <= maxCrossingChecks) {{
            for (let iteration = 0; iteration < crossingIterations; iteration++) {{
                for (let leftIndex = 0; leftIndex < links.length; leftIndex++) {{
                    const leftLink = links[leftIndex];
                    const a = itemById.get(leftLink.sourceId);
                    const b = itemById.get(leftLink.targetId);
                    if (!a || !b) {{
                        continue;
                    }}
                    for (let rightIndex = leftIndex + 1; rightIndex < links.length; rightIndex++) {{
                        const rightLink = links[rightIndex];
                        if (leftLink.sourceId === rightLink.sourceId || leftLink.sourceId === rightLink.targetId || leftLink.targetId === rightLink.sourceId || leftLink.targetId === rightLink.targetId) {{
                            continue;
                        }}
                        const c = itemById.get(rightLink.sourceId);
                        const d = itemById.get(rightLink.targetId);
                        if (!c || !d || !segmentsCross(a, b, c, d)) {{
                            continue;
                        }}
                        let dx = b.x - a.x;
                        let dy = b.y - a.y;
                        let length = Math.hypot(dx, dy);
                        if (length < 1e-6) {{
                            continue;
                        }}
                        const normalX = -dy / length;
                        const normalY = dx / length;
                        const push = 4.5 + (iteration * 1.5);
                        a.x += normalX * push;
                        a.y += normalY * push;
                        b.x += normalX * push;
                        b.y += normalY * push;
                        c.x -= normalX * push;
                        c.y -= normalY * push;
                        d.x -= normalX * push;
                        d.y -= normalY * push;
                    }}
                }}
            }}
        }}

        if (links.length * refined.length <= maxEdgeNodeChecks) {{
            links.forEach(link => {{
                const source = itemById.get(link.sourceId);
                const target = itemById.get(link.targetId);
                if (!source || !target) {{
                    return;
                }}
                refined.forEach(item => {{
                    if (item.componentId === link.sourceId || item.componentId === link.targetId) {{
                        return;
                    }}
                    const hit = pointSegmentDistance(item, source, target);
                    if (hit.t <= 0.03 || hit.t >= 0.97) {{
                        return;
                    }}
                    const minimumDistance = item.radius + edgePadding;
                    if (hit.distance >= minimumDistance) {{
                        return;
                    }}
                    let normalX = item.x - hit.closestX;
                    let normalY = item.y - hit.closestY;
                    let normalLength = Math.hypot(normalX, normalY);
                    if (normalLength < 1e-6) {{
                        const edgeDx = target.x - source.x;
                        const edgeDy = target.y - source.y;
                        normalX = -edgeDy || 1;
                        normalY = edgeDx || 0;
                        normalLength = Math.hypot(normalX, normalY);
                    }}
                    const push = (minimumDistance - hit.distance) * 0.82;
                    item.x += (normalX / normalLength) * push;
                    item.y += (normalY / normalLength) * push;
                }});
            }});
        }}

        if (pairChecks <= maxPairChecks) {{
            for (let iteration = 0; iteration < 4; iteration++) {{
                for (let leftIndex = 0; leftIndex < refined.length; leftIndex++) {{
                    const left = refined[leftIndex];
                    for (let rightIndex = leftIndex + 1; rightIndex < refined.length; rightIndex++) {{
                        const right = refined[rightIndex];
                        let dx = right.x - left.x;
                        let dy = right.y - left.y;
                        let distance = Math.hypot(dx, dy);
                        if (distance < 1e-6) {{
                            dx = (seededUnit(left.componentId + right.componentId, iteration + 41) - 0.5) || 0.01;
                            dy = (seededUnit(left.componentId + right.componentId, iteration + 42) - 0.5) || 0.01;
                            distance = Math.hypot(dx, dy);
                        }}
                        const minimumDistance = left.radius + right.radius + bubblePadding;
                        if (distance >= minimumDistance) {{
                            continue;
                        }}
                        const shift = ((minimumDistance - distance) / 2) * 0.9;
                        const shiftX = (dx / distance) * shift;
                        const shiftY = (dy / distance) * shift;
                        left.x -= shiftX;
                        left.y -= shiftY;
                        right.x += shiftX;
                        right.y += shiftY;
                    }}
                }}
            }}
        }}

        return refined;
    }}

    function seededUnit(componentId, salt) {{
        const raw = Math.sin((componentId + 1) * 12.9898 + salt * 78.233) * 43758.5453;
        return raw - Math.floor(raw);
    }}

    function simulateComponentLayout(componentIds, adjacency, options = {{}}, hierarchyNodes) {{
        const radii = new Map(componentIds.map(nodeId => [nodeId, componentRadiusForSize(hierarchyNodes[nodeId].size)]));
        const positions = new Map();
        const velocities = new Map();
        const anchors = new Map();
        const damping = options.damping ?? 0.8;
        const repulsion = options.repulsion ?? 6000;
        const spring = options.spring ?? 0.07;
        const gravity = options.gravity ?? 0.008;
        const anchorStrength = options.anchorStrength ?? 0;
        const iterations = options.iterations ?? 110;
        const preferredEdgeLength = options.preferredEdgeLength ?? 120;
        const collisionStrength = options.collisionStrength ?? 0.28;
        const maxStep = options.maxStep ?? 48;
        // Above this many nodes the O(n^2) physics is skipped entirely (the tidy seed is already a
        // good, crossing-free layout) so large forests render fast instead of freezing.
        const maxPhysicsNodes = options.maxPhysicsNodes ?? 500;

        // Seed positions AND anchors from the tidy tree layout: crossing-free for forest topology
        // even with long branches. Centered on the origin so gravity pulls toward the middle.
        const tidy = tidyComponentLayout(componentIds, adjacency, hierarchyNodes, 0);
        let centerX = 0;
        let centerY = 0;
        componentIds.forEach(nodeId => {{
            const seed = tidy.positionById.get(nodeId) || {{x: 0, y: 0}};
            centerX += seed.x;
            centerY += seed.y;
        }});
        centerX /= Math.max(1, componentIds.length);
        centerY /= Math.max(1, componentIds.length);
        componentIds.forEach(nodeId => {{
            const seed = tidy.positionById.get(nodeId) || {{x: 0, y: 0}};
            const anchorX = seed.x - centerX;
            const anchorY = seed.y - centerY;
            anchors.set(nodeId, {{x: anchorX, y: anchorY}});
            positions.set(nodeId, {{x: anchorX, y: anchorY}});
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

        if (componentIds.length <= maxPhysicsNodes) {{
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
                    // Floor the repulsion distance so two near-coincident nodes can't produce an
                    // unbounded force that flings a node to infinity (the numerical blow-up that
                    // collapsed fit-to-view to a single visible bubble).
                    const repelForce = repulsion / Math.max(distSq, 16);
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

            let totalMovement = 0;
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
                // Cap per-step displacement to keep the simulation numerically stable.
                const speed = Math.hypot(velocity.x, velocity.y);
                if (speed > maxStep) {{
                    velocity.x *= maxStep / speed;
                    velocity.y *= maxStep / speed;
                }}
                position.x += velocity.x;
                position.y += velocity.y;
                totalMovement += Math.abs(velocity.x) + Math.abs(velocity.y);
            }});
            // Early-exit once the system has settled (avoids the worker hanging on big forests).
            if (iteration > 12 && (totalMovement / componentIds.length) < 0.05) {{
                break;
            }}
        }}
        }}

        const rawItems = componentIds.map(nodeId => {{
            const position = positions.get(nodeId);
            return {{componentId: nodeId, x: position.x, y: position.y, radius: radii.get(nodeId) || 10}};
        }});
        const refinedItems = refineLayoutGeometry(rawItems, edgePairs, {{
            bubblePadding: options.bubblePadding ?? 14,
            edgePadding: options.edgePadding ?? 8,
            overlapIterations: options.geometryIterations ?? 5,
            edgeIterations: options.edgeIterations ?? 3,
            crossingIterations: options.crossingIterations ?? 2,
        }});
        return normalizeComponentLayout(refinedItems, 24);
    }}

    function forceDirectedForestLayout(visibleIds, links, hierarchyNodes, sortBySizeEnabled) {{
        if (visibleIds.length === 0) {{
            clusterCanvas.height = 760;
            return [];
        }}
        clusterCanvas.height = Math.max(760, Math.min(1180, Math.round(window.innerHeight * 0.8)));
        const {{adjacency, components}} = clusterGraphComponents(visibleIds, links, hierarchyNodes, sortBySizeEnabled);
        const componentLayouts = components.map(componentIds => simulateComponentLayout(componentIds, adjacency, {{
            repulsion: 4000,
            spring: 0.08,
            gravity: 0.006,
            damping: 0.82,
            preferredEdgeLength: 124,
            iterations: Math.min(70, 40 + componentIds.length),
            collisionStrength: 0.5,
            bubblePadding: 17,
            edgePadding: 10,
            geometryIterations: 7,
            edgeIterations: 4,
            crossingIterations: 4,
            anchorStrength: 0.12,
            maxPhysicsNodes: 500,
        }}, hierarchyNodes));
        return packLayouts(componentLayouts, {{gapX: 36, gapY: 36, outerPadding: 72}});
    }}

    // Tidy (Reingold-Tilford style) layout for a single component, rooted at its tree center.
    // Returns crossing-free positions for forest topology even with long branches. Shared by the
    // Tree layout and used to seed the Force-directed simulation. O(n), x offset by originX.
    function tidyComponentLayout(componentIds, adjacency, hierarchyNodes, originX) {{
        const rootId = treeCenter(componentIds, adjacency, hierarchyNodes);
        const tree = rootedTree(rootId, adjacency, hierarchyNodes);
        const depthByNode = new Map([[rootId, 0]]);
        tree.order.forEach(nodeId => {{
            (tree.children.get(nodeId) || []).forEach(childId => {{
                depthByNode.set(childId, (depthByNode.get(nodeId) || 0) + 1);
            }});
        }});
        const radii = new Map(componentIds.map(nodeId => [nodeId, componentRadiusForSize(hierarchyNodes[nodeId].size)]));
        const maxRadius = Math.max(...componentIds.map(nodeId => radii.get(nodeId) || 10), 10);
        const siblingGap = Math.max(20, maxRadius * 0.42);
        const levelGap = Math.max(168, (maxRadius * 2.2) + 96);
        const subtreeSpan = new Map();
        [...tree.order].reverse().forEach(nodeId => {{
            const childIds = tree.children.get(nodeId) || [];
            const nodeSpan = ((radii.get(nodeId) || 12) * 2) + 18;
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
            const x = originX + (depthByNode.get(nodeId) || 0) * levelGap;
            if (childIds.length === 0) {{
                positionById.set(nodeId, {{componentId: nodeId, x, y: topY + (nodeSpan / 2), radius: radii.get(nodeId) || 10}});
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
            positionById.set(nodeId, {{componentId: nodeId, x, y: centerY, radius: radii.get(nodeId) || 10}});
        }}

        placeNode(rootId, 0);
        const maxDepth = tree.order.reduce((maxValue, nodeId) => Math.max(maxValue, depthByNode.get(nodeId) || 0), 0);
        const width = (maxDepth * levelGap) + (Math.max(...tree.order.map(nodeId => radii.get(nodeId) || 10), 10) * 2);
        return {{positionById, order: tree.order, width}};
    }}

    function tidyForestLayout(visibleIds, links, hierarchyNodes, sortBySizeEnabled) {{
        if (visibleIds.length === 0) {{
            clusterCanvas.height = 760;
            return [];
        }}

        const {{adjacency, components}} = clusterGraphComponents(visibleIds, links, hierarchyNodes, sortBySizeEnabled);

        clusterCanvas.height = Math.max(760, Math.min(1180, Math.round(window.innerHeight * 0.8)));
        const outerPadding = 58;
        const componentGap = 110;
        const layout = [];
        let currentX = outerPadding;

        components.forEach(componentIds => {{
            const tidy = tidyComponentLayout(componentIds, adjacency, hierarchyNodes, currentX);
            tidy.order.forEach(nodeId => layout.push(tidy.positionById.get(nodeId)));
            currentX += tidy.width + outerPadding + componentGap;
        }});

        return layout;
    }}

    function layoutCacheKey(visibleIds, links, algorithm) {{
        const flags = sortComponentsBySizeEnabled() ? 1 : 0;
        let h = flags;
        for (let i = 0; i < visibleIds.length; i++) {{
            h = (Math.imul(h, 1664525) + visibleIds[i] + 1013904223) | 0;
        }}
        let lh = 0;
        for (let i = 0; i < links.length; i++) {{
            const link = links[i];
            lh ^= (Math.imul((link.sourceId * 31 + link.targetId) | 0, 2654435761) | 0) + Math.round((link.weight || 0) * 10000);
        }}
        h = (Math.imul(h, 1664525) + lh + 1013904223) | 0;
        return algorithm + ':' + flags + ':' + visibleIds.length + ':' + links.length + ':' + h;
    }}

    function computeVisibleLayout(visibleIds, links, algorithm, hierarchyNodes, sortBySizeEnabled) {{
        const key = layoutCacheKey(visibleIds, links, algorithm);
        const cached = state.layoutCache.get(key);
        if (cached !== undefined) {{
            return {{layout: cached.map(item => ({{...item}})), key, async: false}};
        }}

        if (algorithm === 'force' && state.layoutWorker) {{
            const workerHierarchyNodes = {{}};
            visibleIds.forEach(id => {{
                const node = hierarchyNodes[id];
                workerHierarchyNodes[id] = {{size: node.size, leaf_start: node.leaf_start}};
            }});
            const requestId = ++state.pendingLayoutRequestId;
            state.layoutWorker.postMessage({{
                type: 'computeLayout',
                requestId,
                key,
                algorithm,
                visibleIds,
                links: links.map(link => ({{sourceId: link.sourceId, targetId: link.targetId, weight: link.weight}})),
                hierarchyNodes: workerHierarchyNodes,
                sortBySizeEnabled,
            }});
            return {{layout: null, key, async: true}};
        }}

        const layout = (() => {{
            if (algorithm === 'grid') {{
                return gridClusterLayout(visibleIds);
            }}
            if (algorithm === 'packed') {{
                return packedClusterLayout(visibleIds);
            }}
            if (algorithm === 'force') {{
                return forceDirectedForestLayout(visibleIds, links, hierarchyNodes, sortBySizeEnabled);
            }}
            return tidyForestLayout(visibleIds, links, hierarchyNodes, sortBySizeEnabled);
        }})();
        state.layoutCache.set(key, layout.map(item => ({{...item}})));
        if (state.layoutCache.size > 24) {{
            const oldestKey = state.layoutCache.keys().next().value;
            state.layoutCache.delete(oldestKey);
        }}
        return {{layout, key, async: false}};
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

    function trimmedLinkEndpoints(link) {{
        const dx = link.right.x - link.left.x;
        const dy = link.right.y - link.left.y;
        const distance = Math.hypot(dx, dy);
        if (distance < 1e-6) {{
            return {{startX: link.left.x, startY: link.left.y, endX: link.right.x, endY: link.right.y}};
        }}
        const unitX = dx / distance;
        const unitY = dy / distance;
        const leftOffset = Math.min(distance / 2, link.left.radius + 3);
        const rightOffset = Math.min(distance / 2, link.right.radius + 3);
        return {{
            startX: link.left.x + (unitX * leftOffset),
            startY: link.left.y + (unitY * leftOffset),
            endX: link.right.x - (unitX * rightOffset),
            endY: link.right.y - (unitY * rightOffset),
        }};
    }}

    function renderedLinkSegments(link) {{
        const trimmed = trimmedLinkEndpoints(link);
        if (currentLayoutAlgorithm() !== 'tree') {{
            return [trimmed];
        }}
        const middleX = (trimmed.startX + trimmed.endX) / 2;
        return [
            {{startX: trimmed.startX, startY: trimmed.startY, endX: middleX, endY: trimmed.startY}},
            {{startX: middleX, startY: trimmed.startY, endX: middleX, endY: trimmed.endY}},
            {{startX: middleX, startY: trimmed.endY, endX: trimmed.endX, endY: trimmed.endY}},
        ].filter(segment => Math.hypot(segment.endX - segment.startX, segment.endY - segment.startY) > 1e-6);
    }}

    function componentDotGeometry(component, item) {{
        const sampleCount = componentDotCount(component.size);
        const dotRadius = componentDotRadius(component.size, item.radius);
        const packingRadius = Math.max(0, item.radius - dotRadius - 0.6);
        return {{sampleCount, dotRadius, packingRadius}};
    }}

    function componentDotLayout(component, item) {{
        const {{sampleCount, dotRadius, packingRadius}} = componentDotGeometry(component, item);
        const minimumDistance = packingRadius > 0 ? Math.min(0.22, (dotRadius * 2.12) / packingRadius) : 0;
        const normalizedLayout = normalizedComponentDotLayout(component.id, sampleCount, minimumDistance);
        return normalizedLayout.map(layout => ({{
            memberIndex: layout.memberIndex,
            x: item.x + (layout.x * packingRadius),
            y: item.y + (layout.y * packingRadius),
            radius: dotRadius,
        }}));
    }}

    function componentSelectionState(members) {{
        let selectedCount = 0;
        members.forEach(nodeIndex => {{
            if (state.selectedNodeIndices.has(nodeIndex)) {{
                selectedCount += 1;
            }}
        }});
        return {{
            selectedCount,
            anySelected: selectedCount > 0,
            allSelected: members.length > 0 && selectedCount === members.length,
        }};
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

        const selectedNodeOutlines = [];
        const TAU = Math.PI * 2;
        const drawScale = Math.max(state.viewTransform.scale, 1e-9);

        // World-space viewport bounds for culling (#4)
        const worldMinX = (0 - state.viewTransform.offsetX) / drawScale;
        const worldMaxX = (clusterCanvas.width - state.viewTransform.offsetX) / drawScale;
        const worldMinY = (0 - state.viewTransform.offsetY) / drawScale;
        const worldMaxY = (clusterCanvas.height - state.viewTransform.offsetY) / drawScale;

        clusterContext.save();
        clusterContext.translate(state.viewTransform.offsetX, state.viewTransform.offsetY);
        clusterContext.scale(drawScale, drawScale);

        state.splitLinks.forEach(link => {{
            clusterContext.strokeStyle = 'rgba(92, 106, 112, 0.42)';
            clusterContext.lineWidth = 1.6 / drawScale;
            clusterContext.beginPath();
            const segments = renderedLinkSegments(link);
            segments.forEach((segment, index) => {{
                if (index === 0) {{
                    clusterContext.moveTo(segment.startX, segment.startY);
                }} else {{
                    clusterContext.lineTo(segment.startX, segment.startY);
                }}
                clusterContext.lineTo(segment.endX, segment.endY);
            }});
            clusterContext.stroke();
        }});

        // Collect dots batched by color for batch drawing (#3)
        const dotColorBuckets = new Map();

        const showClusterBounds = renderClusterBoundsEnabled();
        const showNodes = renderNodesEnabled();

        state.visibleLayout.forEach(item => {{
            // Viewport culling (#4)
            if (item.x + item.radius < worldMinX || item.x - item.radius > worldMaxX ||
                item.y + item.radius < worldMinY || item.y - item.radius > worldMaxY) {{
                return;
            }}

            const component = state.bundle.graph.hierarchy.nodes[item.componentId];
            const members = componentMembers(item.componentId);
            const selectionState = componentSelectionState(members);

            if (showClusterBounds) {{
                clusterContext.fillStyle = 'rgba(248, 243, 235, 0.95)';
                clusterContext.strokeStyle = selectionState.allSelected ? '#1e2a2f' : '#c8b8a6';
                clusterContext.lineWidth = (selectionState.allSelected ? 3 : 1.5) / drawScale;
                clusterContext.beginPath();
                clusterContext.arc(item.x, item.y, item.radius, 0, TAU);
                clusterContext.fill();
                clusterContext.stroke();
            }}

            if (!showNodes) {{
                return;
            }}

            const dotLayout = componentDotLayout(component, item);
            for (const dot of dotLayout) {{
                // Use pre-computed color cache (#2)
                const color = state.nodeColorCache[dot.memberIndex] ?? nodeColor(dot.memberIndex);
                let bucket = dotColorBuckets.get(color);
                if (bucket === undefined) {{
                    bucket = [];
                    dotColorBuckets.set(color, bucket);
                }}
                bucket.push(dot);
                if (!selectionState.allSelected && state.selectedNodeIndices.has(dot.memberIndex)) {{
                    selectedNodeOutlines.push({{x: dot.x, y: dot.y, radius: dot.radius}});
                }}
            }}
        }});

        // Batch draw all dots grouped by color (#3)
        dotColorBuckets.forEach((dots, color) => {{
            clusterContext.fillStyle = color;
            clusterContext.beginPath();
            for (const dot of dots) {{
                clusterContext.moveTo(dot.x + dot.radius, dot.y);
                clusterContext.arc(dot.x, dot.y, dot.radius, 0, TAU);
            }}
            clusterContext.fill();
        }});

        clusterContext.restore();

        if (selectedNodeOutlines.length > 0) {{
            clusterContext.save();
            clusterContext.strokeStyle = '#1e2a2f';
            clusterContext.lineWidth = 1.8;
            selectedNodeOutlines.forEach(outline => {{
                const screenPoint = worldToScreenPoint(outline.x, outline.y);
                const screenRadius = Math.max(3, (outline.radius * state.viewTransform.scale) + 1.6);
                clusterContext.beginPath();
                clusterContext.arc(screenPoint.x, screenPoint.y, screenRadius, 0, TAU);
                clusterContext.stroke();
            }});
            clusterContext.restore();
        }}

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
            const screenPoint = worldToScreenPoint(item.x, item.y);
            const screenRadius = item.radius * state.viewTransform.scale;
            // Viewport culling for labels (#4)
            if (screenPoint.x + screenRadius < 0 || screenPoint.x - screenRadius > clusterCanvas.width ||
                screenPoint.y + screenRadius < 0 || screenPoint.y - screenRadius > clusterCanvas.height) {{
                return;
            }}
            const component = state.bundle.graph.hierarchy.nodes[item.componentId];
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

        if (state.layoutComputing) {{
            const cx = clusterCanvas.width / 2;
            const cy = clusterCanvas.height / 2;
            clusterContext.save();
            clusterContext.fillStyle = 'rgba(255, 250, 242, 0.72)';
            clusterContext.fillRect(cx - 120, cy - 22, 240, 44);
            clusterContext.fillStyle = '#5c6a70';
            clusterContext.font = '600 14px Georgia';
            clusterContext.textAlign = 'center';
            clusterContext.textBaseline = 'middle';
            clusterContext.fillText('Computing layout…', cx, cy);
            clusterContext.restore();
        }}
    }}

    function scheduleClusterRender() {{
        if (state.pendingClusterRenderFrame !== null) {{
            return;
        }}
        state.pendingClusterRenderFrame = window.requestAnimationFrame(() => {{
            state.pendingClusterRenderFrame = null;
            renderClusterView();
        }});
    }}

    function applyComputedLayout(layout, visibleGraph, layoutAlgorithm, resetView) {{
        const layoutById = new Map(layout.map(item => [item.componentId, item]));
        const links = (layoutAlgorithm === 'grid' || layoutAlgorithm === 'packed')
            ? []
            : visibleGraph.links
                .map(link => {{
                    const left = layoutById.get(link.sourceId);
                    const right = layoutById.get(link.targetId);
                    if (!left || !right) {{ return null; }}
                    return {{left, right, threshold: link.weight}};
                }})
                .filter(Boolean);
        state.visibleLayout = layout;
        state.splitLinks = links;
        document.getElementById('stat-clusters').textContent = layout.length.toLocaleString();
        document.getElementById('stat-links').textContent = links.length.toLocaleString();
        if (resetView) {{ fitClusterViewToLayout(); }}
        renderClusterView();
    }}

    function drawClusterView(resetView = true) {{
        if (!state.bundle) {{
            renderClusterView();
            return;
        }}

        const minClusterSize = Math.max(1, Number(document.getElementById('min-cluster-size').value) || 1);
        const layoutAlgorithm = currentLayoutAlgorithm();
        const thresholdValue = selectedThresholdValue();
        const activeClusterIds = activeClustersAtThreshold(thresholdValue);
        const visibleGraph = mstLinksForActiveClusters(activeClusterIds, minClusterSize, leafPruningOnlyEnabled());
        const hierarchyNodes = state.bundle.graph.hierarchy.nodes;
        const sortBySizeEnabled = sortComponentsBySizeEnabled();

        state.activeClusters = activeClusterIds;
        state.visibleClusters = visibleGraph.visibleIds;

        const hidden = visibleGraph.hiddenNodes;
        const shown = state.bundle.graph.nodes.length - hidden;
        document.getElementById('hidden-summary').textContent = hidden.toLocaleString() + ' nodes hidden by minimum cluster size';
        document.getElementById('stat-shown-nodes').textContent = shown.toLocaleString();
        document.getElementById('stat-hidden-nodes').textContent = hidden.toLocaleString();

        const result = computeVisibleLayout(visibleGraph.visibleIds, visibleGraph.links, layoutAlgorithm, hierarchyNodes, sortBySizeEnabled);
        if (result.async) {{
            state._pendingVisibleGraph = visibleGraph;
            state._pendingLayoutAlgorithm = layoutAlgorithm;
            state._pendingLayoutResetView = resetView;
            state.layoutComputing = true;
            renderClusterView();
            return;
        }}
        state.layoutComputing = false;
        applyComputedLayout(result.layout, visibleGraph, layoutAlgorithm, resetView);
    }}

    function htmlEscape(value) {{
        return String(value).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
    }}

    function setupLayoutWorker() {{
        try {{
            const blob = new Blob([LAYOUT_WORKER_CODE], {{type: 'application/javascript'}});
            const workerUrl = URL.createObjectURL(blob);
            const worker = new Worker(workerUrl);
            worker.onmessage = function(event) {{
                const {{requestId, key, layout}} = event.data;
                if (requestId !== state.pendingLayoutRequestId) {{ return; }}
                state.layoutCache.set(key, layout.map(item => ({{...item}})));
                if (state.layoutCache.size > 24) {{
                    state.layoutCache.delete(state.layoutCache.keys().next().value);
                }}
                state.layoutComputing = false;
                applyComputedLayout(layout, state._pendingVisibleGraph, state._pendingLayoutAlgorithm, state._pendingLayoutResetView);
            }};
            worker.onerror = function(error) {{
                console.warn('Layout worker error, falling back to main thread:', error);
                state.layoutWorker = null;
                state.layoutComputing = false;
                renderClusterView();
            }};
            state.layoutWorker = worker;
        }} catch (error) {{
            console.warn('Layout worker unavailable, using main thread:', error);
        }}
    }}

    function setupMetadataTableDelegation() {{
        const thead = document.querySelector('#metadata-table thead');
        const tbody = document.querySelector('#metadata-table tbody');

        thead.addEventListener('click', event => {{
            const sortBtn = event.target.closest('[data-column-key]');
            if (sortBtn) {{ toggleMetadataSort(sortBtn.dataset.columnKey); return; }}
            const copyBtn = event.target.closest('[data-copy-column]');
            if (copyBtn) {{ event.preventDefault(); event.stopPropagation(); copyMetadataColumn(copyBtn.dataset.copyColumn); }}
        }});

        thead.addEventListener('pointerdown', event => {{
            const resizeHandle = event.target.closest('[data-resize-column]');
            if (resizeHandle) {{ startMetadataColumnResize(resizeHandle.dataset.resizeColumn, event); }}
        }});

        tbody.addEventListener('mousedown', event => {{
            if (event.shiftKey && !event.target.closest('button, a, input, select, textarea')) {{
                event.preventDefault();
            }}
        }});

        tbody.addEventListener('click', event => {{
            if (event.target.closest('button, a, input, select, textarea')) {{ return; }}
            if (!event.shiftKey && metadataTableHasActiveTextSelection()) {{ return; }}
            const row = event.target.closest('tr[data-node-index]');
            if (!row) {{ return; }}
            toggleMetadataRowSelection(Number(row.dataset.nodeIndex), state.renderedNodeIndices, {{range: event.shiftKey}});
        }});

        tbody.addEventListener('keydown', event => {{
            if (event.key !== 'Enter' && event.key !== ' ') {{ return; }}
            const row = event.target.closest('tr[data-node-index]');
            if (!row) {{ return; }}
            event.preventDefault();
            toggleMetadataRowSelection(Number(row.dataset.nodeIndex), state.renderedNodeIndices, {{range: event.shiftKey}});
        }});
    }}

    function updateMetadataTable() {{
        const selected = Array.from(state.selectedNodeIndices).sort((left, right) => left - right);
        const baseNodeIndices = metadataBaseNodeIndices();
        const filteredNodeIndices = filteredMetadataNodeIndices(baseNodeIndices);
        const sortedNodeIndices = sortedMetadataNodeIndices(filteredNodeIndices);
        const pagination = metadataPagination(sortedNodeIndices.length);
        state.metadataPage = pagination.pageIndex;
        const renderedNodeIndices = sortedNodeIndices.slice(pagination.start, pagination.end);
        state.renderedNodeIndices = renderedNodeIndices;
        pruneMetadataRowSelection(renderedNodeIndices);
        applyMetadataColumnWidths();
        const thead = document.querySelector('#metadata-table thead');
        const tbody = document.querySelector('#metadata-table tbody');

        const columnKeys = metadataColumnKeys();
        const headerCells = columnKeys.map(label => {{
            const escapedLabel = htmlEscape(label);
            const ariasort = state.metadataSort.columnKey === label
                ? (state.metadataSort.direction === 'asc' ? 'ascending' : 'descending')
                : 'none';
            const indicator = htmlEscape(metadataSortIndicator(label));
            return '<th class="metadata-header" aria-sort="' + ariasort + '">'
                + '<div class="metadata-header-cell">'
                + '<button type="button" class="metadata-sort-button" data-column-key="' + escapedLabel + '" title="Sort by ' + escapedLabel + '">'
                + '<span>' + escapedLabel + '</span><span class="metadata-sort-indicator">' + indicator + '</span></button>'
                + '<button type="button" class="metadata-copy-button" data-copy-column="' + escapedLabel + '" title="Copy the currently displayed values from ' + escapedLabel + '">Copy</button>'
                + '</div>'
                + '<div class="metadata-resize-handle" role="separator" aria-orientation="vertical" title="Drag to resize column" data-resize-column="' + escapedLabel + '"></div>'
                + '</th>';
        }});
        thead.innerHTML = '<tr>' + headerCells.join('') + '</tr>';

        let bodyHtml;
        if (sortedNodeIndices.length === 0) {{
            const emptyMsg = htmlEscape(metadataFilterText()
                ? 'No metadata rows match the current filter.'
                : 'No metadata rows to display.');
            bodyHtml = '<tr class="metadata-empty-row"><td colspan="' + (state.metadataColumns.length + 1) + '">' + emptyMsg + '</td></tr>';
        }} else {{
            bodyHtml = renderedNodeIndices.map(nodeIndex => {{
                const isSelected = state.selectedMetadataNodeIndices.has(nodeIndex);
                const rowClass = isSelected ? ' class="metadata-row-selected"' : '';
                const ariaSelected = isSelected ? 'true' : 'false';
                const idText = htmlEscape(nodeId(nodeIndex));
                let cellsHtml = '<td title="' + idText + '">' + idText + '</td>';
                for (const column of state.metadataColumns) {{
                    const value = metadataValue(nodeIndex, column.name);
                    const cellClass = metadataCellClass(column.name, value);
                    const displayText = htmlEscape(formatMetadataDisplayValue(column.name, value));
                    cellsHtml += '<td' + (cellClass ? ' class="' + cellClass + '"' : '') + ' title="' + displayText + '">' + displayText + '</td>';
                }}
                return '<tr' + rowClass + ' tabindex="0" aria-selected="' + ariaSelected + '" data-node-index="' + nodeIndex + '">' + cellsHtml + '</tr>';
            }}).join('');
        }}
        tbody.innerHTML = bodyHtml;

        document.getElementById('selection-summary').textContent = selected.length.toLocaleString() + ' nodes selected';
        const pageStatus = document.getElementById('metadata-page-status');
        if (sortedNodeIndices.length === 0) {{
            pageStatus.textContent = metadataFilterText()
                ? 'No rows match the current filter.'
                : 'No metadata rows available.';
        }} else if (pagination.showAll) {{
            pageStatus.textContent = 'Showing all ' + sortedNodeIndices.length.toLocaleString() + ' rows on one page.';
        }} else {{
            pageStatus.textContent = 'Showing rows ' + (pagination.start + 1).toLocaleString() + ' to ' + pagination.end.toLocaleString() +
                ' of ' + sortedNodeIndices.length.toLocaleString() + ' (page ' + (pagination.pageIndex + 1).toLocaleString() +
                ' of ' + pagination.pageCount.toLocaleString() + ').';
        }}
        document.getElementById('metadata-prev-page').disabled = !state.bundle || sortedNodeIndices.length === 0 || pagination.pageIndex === 0;
        document.getElementById('metadata-next-page').disabled = !state.bundle || sortedNodeIndices.length === 0 || pagination.pageIndex >= pagination.pageCount - 1;
        const sortDescription = metadataSortDescription();
        const filterText = metadataFilterText();
        const filterDescription = filterText ? ' matching filter "' + filterText + '"' : '';
        const pagerActive = pagination.pageCount > 1;
        state.metadataBaseNoteText = (selected.length === 0
            ? (
                pagerActive
                    ? 'No clusters selected. Browse the full network table with the pager' + filterDescription + (sortDescription ? ', sorted by ' + sortDescription : '') + '.'
                    : 'No clusters selected. Showing the full network table' + filterDescription + (sortDescription ? ', sorted by ' + sortDescription : '') + '.'
            )
            : (
                pagerActive
                    ? 'Showing one page of the selected rows in the table' + filterDescription + (sortDescription ? ', sorted by ' + sortDescription : '') + '. Use the pager to browse more rows. Export includes the full filtered selection.'
                    : 'Ctrl-click a node to toggle it individually, click a cluster to toggle it, Shift-drag a box to add multiple clusters, click table rows to stage them, shift-click to select or deselect row ranges, use the search box to filter rows, and click a column header to sort' + (sortDescription ? ' by ' + sortDescription : '') + '.'
            ));
        document.getElementById('export-selected').disabled = !state.bundle;
        document.getElementById('metadata-reset-sort').disabled = !state.bundle || !state.metadataSort.columnKey;
        document.getElementById('clear-selection').disabled = selected.length === 0;
        applyMetadataTableRowHighlights();
    }}

    function applyMetadataTableRowHighlights() {{
        const tbody = document.querySelector('#metadata-table tbody');
        if (tbody) {{
            tbody.querySelectorAll('tr[data-node-index]').forEach(row => {{
                const idx = Number(row.dataset.nodeIndex);
                const sel = state.selectedMetadataNodeIndices.has(idx);
                row.className = sel ? 'metadata-row-selected' : '';
                row.setAttribute('aria-selected', sel ? 'true' : 'false');
            }});
        }}
        const metadataSelectionCount = state.selectedMetadataNodeIndices.size;
        const metadataSelectionDescription = metadataSelectionCount > 0
            ? metadataSelectionCount.toLocaleString() + ' table rows selected. Click Select nodes to promote them into the graph selection. '
            : '';
        document.getElementById('selection-note').textContent = metadataSelectionDescription + state.metadataBaseNoteText;
        document.getElementById('metadata-select-nodes').disabled = !state.bundle || metadataSelectionCount === 0;
        document.getElementById('metadata-deselect-rows').disabled = !state.bundle || metadataSelectionCount === 0;
    }}

    function exportSelection() {{
        const selected = Array.from(state.selectedNodeIndices).sort((left, right) => left - right);
        if (!state.bundle) {{
            return;
        }}
        const activeClusterIds = activeClustersAtThreshold(selectedThresholdValue());
        const componentAssignments = activeClusterAssignments(activeClusterIds);
        const rankedClusterIds = [...activeClusterIds].sort((leftId, rightId) => {{
            const leftNode = state.bundle.graph.hierarchy.nodes[leftId];
            const rightNode = state.bundle.graph.hierarchy.nodes[rightId];
            return rightNode.size - leftNode.size || leftNode.leaf_start - rightNode.leaf_start || leftId - rightId;
        }});
        const clusterByComponentId = new Map(rankedClusterIds.map((componentId, clusterIndex) => [componentId, clusterIndex + 1]));
        const metadataColumns = state.metadataColumns.filter(column => column.name !== 'SSN_cluster');
        const exportedNodeIndices = metadataDisplayNodeIndices(selected.length > 0
            ? selected
            : state.bundle.graph.nodes.map((_, nodeIndex) => nodeIndex));
        const header = ['node_id', 'SSN_cluster', ...metadataColumns.map(column => column.name)];
        const lines = [header.join('\\t')];
        exportedNodeIndices.forEach(nodeIndex => {{
            const row = [
                nodeId(nodeIndex),
                clusterByComponentId.get(componentAssignments[nodeIndex]) ?? '',
                ...metadataColumns.map(column => {{
                    const value = metadataValue(nodeIndex, column.name);
                    return value === null || value === undefined ? '' : String(value);
                }}),
            ];
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
    }}

    function scheduleThresholdUI(resetView = true) {{
        if (!state.bundle) {{
            return;
        }}
        state.pendingThresholdUIResetView = state.pendingThresholdUIResetView || resetView;
        if (state.pendingThresholdUIFrame !== null) {{
            return;
        }}
        state.pendingThresholdUIFrame = window.requestAnimationFrame(() => {{
            const nextResetView = state.pendingThresholdUIResetView;
            state.pendingThresholdUIFrame = null;
            state.pendingThresholdUIResetView = false;
            updateThresholdUI(nextResetView);
        }});
    }}

    function metadataSortIndicator(columnKey) {{
        if (state.metadataSort.columnKey !== columnKey) {{
            return '↕';
        }}
        return state.metadataSort.direction === 'asc' ? '↑' : '↓';
    }}

    function metadataSortDescription() {{
        if (!state.metadataSort.columnKey) {{
            return null;
        }}
        return state.metadataSort.columnKey + ' (' + (state.metadataSort.direction === 'asc' ? 'ascending' : 'descending') + ')';
    }}

    function metadataFilterText() {{
        const input = document.getElementById('metadata-filter');
        return input ? input.value.trim().toLowerCase() : '';
    }}

    function metadataNullPlacement() {{
        const select = document.getElementById('metadata-null-order');
        return select ? (select.value || 'last') : 'last';
    }}

    function isMissingMetadataValue(value) {{
        return value === null || value === undefined || value === '';
    }}

    function metadataSortValue(nodeIndex, columnKey) {{
        if (columnKey === 'node_id') {{
            return nodeId(nodeIndex);
        }}
        const columnIndex = state.metadataColumnIndexByName.get(columnKey);
        if (columnIndex === undefined) {{
            return null;
        }}
        return state.metadataByNodeIndex[nodeIndex]?.[columnIndex] ?? null;
    }}

    function metadataMatchesFilter(nodeIndex, filterText) {{
        if (!filterText) {{
            return true;
        }}
        return (state.metadataSearchTextByNodeIndex[nodeIndex] || '').includes(filterText);
    }}

    function scheduleMetadataFilterUpdate(delayMs = 120) {{
        if (state.pendingMetadataFilterTimer !== null) {{
            window.clearTimeout(state.pendingMetadataFilterTimer);
        }}
        state.pendingMetadataFilterTimer = window.setTimeout(() => {{
            state.pendingMetadataFilterTimer = null;
            resetMetadataPage();
            updateMetadataTable();
        }}, delayMs);
    }}

    function filteredMetadataNodeIndices(nodeIndices) {{
        const filterText = metadataFilterText();
        if (!filterText) {{
            return [...nodeIndices];
        }}
        return nodeIndices.filter(nodeIndex => metadataMatchesFilter(nodeIndex, filterText));
    }}

    function metadataDisplayNodeIndices(nodeIndices) {{
        return sortedMetadataNodeIndices(filteredMetadataNodeIndices(nodeIndices));
    }}

    function metadataRowsPerPageSetting() {{
        const select = document.getElementById('metadata-rows-per-page');
        const value = select ? select.value : String(DEFAULT_METADATA_PAGE_SIZE);
        if (value === 'all') {{
            return {{showAll: true, pageSize: Number.POSITIVE_INFINITY}};
        }}
        const pageSize = Number.parseInt(value, 10);
        if (!Number.isFinite(pageSize) || pageSize <= 0) {{
            return {{showAll: false, pageSize: DEFAULT_METADATA_PAGE_SIZE}};
        }}
        return {{showAll: false, pageSize}};
    }}

    function metadataPagination(totalRowCount) {{
        const rowsPerPage = metadataRowsPerPageSetting();
        if (rowsPerPage.showAll) {{
            return {{pageCount: 1, pageIndex: 0, start: 0, end: totalRowCount, showAll: true}};
        }}
        const pageCount = Math.max(1, Math.ceil(totalRowCount / rowsPerPage.pageSize));
        const pageIndex = Math.max(0, Math.min(state.metadataPage, pageCount - 1));
        const start = totalRowCount === 0 ? 0 : pageIndex * rowsPerPage.pageSize;
        const end = Math.min(totalRowCount, start + rowsPerPage.pageSize);
        return {{pageCount, pageIndex, start, end, showAll: false}};
    }}

    function resetMetadataPage() {{
        state.metadataPage = 0;
    }}

    function stepMetadataPage(delta) {{
        if (!state.bundle) {{
            return;
        }}
        const totalRowCount = metadataDisplayNodeIndices(metadataBaseNodeIndices()).length;
        const pagination = metadataPagination(totalRowCount);
        const nextPage = Math.max(0, Math.min(pagination.pageCount - 1, pagination.pageIndex + delta));
        if (nextPage === pagination.pageIndex) {{
            return;
        }}
        state.metadataPage = nextPage;
        updateMetadataTable();
    }}

    function metadataBaseNodeIndices() {{
        const selected = Array.from(state.selectedNodeIndices).sort((left, right) => left - right);
        if (selected.length > 0) {{
            return selected;
        }}
        return state.allNodeIndices;
    }}

    function metadataColumnValues(nodeIndices, columnKey) {{
        if (columnKey === 'node_id') {{
            return nodeIndices.map(nodeIndex => nodeId(nodeIndex));
        }}
        return nodeIndices.map(nodeIndex => formatMetadataDisplayValue(columnKey, metadataSortValue(nodeIndex, columnKey)));
    }}

    async function copyMetadataColumn(columnKey) {{
        if (!state.bundle) {{
            return;
        }}
        const nodeIndices = metadataDisplayNodeIndices(metadataBaseNodeIndices());
        const values = metadataColumnValues(nodeIndices, columnKey);
        try {{
            await writeTextToClipboard(values.join('\\n'));
            setStatus('Copied column "' + columnKey + '" for ' + values.length.toLocaleString() + ' rows.');
        }} catch (error) {{
            console.error(error);
            setStatus('Failed to copy column "' + columnKey + '": ' + error.message);
        }}
    }}

    function pruneMetadataRowSelection(visibleNodeIndices) {{
        const visibleSet = new Set(visibleNodeIndices);
        if (state.metadataRowSelectionAnchor !== null && !visibleSet.has(state.metadataRowSelectionAnchor)) {{
            state.metadataRowSelectionAnchor = null;
        }}
        if (state.selectedMetadataNodeIndices.size === 0) {{
            return;
        }}
        state.selectedMetadataNodeIndices = new Set(
            Array.from(state.selectedMetadataNodeIndices).filter(nodeIndex => visibleSet.has(nodeIndex))
        );
    }}

    function clearMetadataRowSelection() {{
        state.selectedMetadataNodeIndices = new Set();
        state.metadataRowSelectionAnchor = null;
    }}

    function metadataTableHasActiveTextSelection() {{
        const selection = window.getSelection ? window.getSelection() : null;
        if (!selection || selection.rangeCount === 0 || selection.isCollapsed) {{
            return false;
        }}
        const metadataTable = document.getElementById('metadata-table');
        if (!metadataTable) {{
            return false;
        }}
        for (let rangeIndex = 0; rangeIndex < selection.rangeCount; rangeIndex += 1) {{
            const range = selection.getRangeAt(rangeIndex);
            if (metadataTable.contains(range.commonAncestorContainer)) {{
                return true;
            }}
        }}
        return false;
    }}

    function toggleMetadataRowSelection(nodeIndex, renderedNodeIndices, options = {{}}) {{
        const nextSelection = new Set(state.selectedMetadataNodeIndices);
        if (options.range && state.metadataRowSelectionAnchor !== null) {{
            const anchorPosition = renderedNodeIndices.indexOf(state.metadataRowSelectionAnchor);
            const currentPosition = renderedNodeIndices.indexOf(nodeIndex);
            if (anchorPosition >= 0 && currentPosition >= 0) {{
                const start = Math.min(anchorPosition, currentPosition);
                const end = Math.max(anchorPosition, currentPosition);
                const shouldSelect = !nextSelection.has(nodeIndex);
                for (let position = start; position <= end; position += 1) {{
                    const rangeNodeIndex = renderedNodeIndices[position];
                    if (shouldSelect) {{
                        nextSelection.add(rangeNodeIndex);
                    }} else {{
                        nextSelection.delete(rangeNodeIndex);
                    }}
                }}
            }} else if (nextSelection.has(nodeIndex)) {{
                nextSelection.delete(nodeIndex);
            }} else {{
                nextSelection.add(nodeIndex);
            }}
        }} else if (nextSelection.has(nodeIndex)) {{
            nextSelection.delete(nodeIndex);
        }} else {{
            nextSelection.add(nodeIndex);
        }}
        state.selectedMetadataNodeIndices = nextSelection;
        state.metadataRowSelectionAnchor = nodeIndex;
        applyMetadataTableRowHighlights();
    }}

    function selectNodesFromMetadataRows() {{
        if (state.selectedMetadataNodeIndices.size === 0) {{
            return;
        }}
        state.selectedNodeIndices = new Set(state.selectedMetadataNodeIndices);
        resetMetadataPage();
        clearMetadataRowSelection();
        renderClusterView();
        updateMetadataTable();
    }}

    function metadataColumnType(columnKey) {{
        if (columnKey === 'node_id') {{
            return 'string';
        }}
        return metadataColumn(columnKey)?.type || 'string';
    }}

    function formatMetadataDisplayValue(columnKey, value) {{
        if (isMissingMetadataValue(value)) {{
            return '—';
        }}
        const columnType = metadataColumnType(columnKey);
        if (columnType === 'int') {{
            const numeric = typeof value === 'number' ? value : Number(value);
            return Number.isFinite(numeric) ? Math.trunc(numeric).toLocaleString() : String(value);
        }}
        if (columnType === 'float') {{
            const numeric = typeof value === 'number' ? value : Number(value);
            if (!Number.isFinite(numeric)) {{
                return String(value);
            }}
            const absValue = Math.abs(numeric);
            if (absValue >= 10000 || (absValue > 0 && absValue < 0.001)) {{
                return numeric.toExponential(3);
            }}
            return new Intl.NumberFormat(undefined, {{maximumFractionDigits: 4}}).format(numeric);
        }}
        if (columnType === 'bool' || columnType === 'boolean') {{
            return value ? 'True' : 'False';
        }}
        if (Array.isArray(value)) {{
            return value.join(', ');
        }}
        return String(value);
    }}

    function metadataCellClass(columnKey, value) {{
        if (isMissingMetadataValue(value)) {{
            return 'metadata-cell-null';
        }}
        const columnType = metadataColumnType(columnKey);
        return columnType === 'int' || columnType === 'float' ? 'metadata-cell-number' : '';
    }}

    function compareMetadataValues(leftValue, rightValue) {{
        if (leftValue === rightValue) {{
            return 0;
        }}
        if (typeof leftValue === 'number' && typeof rightValue === 'number') {{
            return leftValue - rightValue;
        }}
        const leftText = String(leftValue);
        const rightText = String(rightValue);
        const leftNumeric = Number(leftText);
        const rightNumeric = Number(rightText);
        if (Number.isFinite(leftNumeric) && Number.isFinite(rightNumeric) && leftText.trim() !== '' && rightText.trim() !== '') {{
            return leftNumeric - rightNumeric;
        }}
        return leftText.localeCompare(rightText, undefined, {{numeric: true, sensitivity: 'base'}});
    }}

    function sortedMetadataNodeIndices(nodeIndices) {{
        if (!state.metadataSort.columnKey) {{
            return [...nodeIndices];
        }}
        const direction = state.metadataSort.direction === 'desc' ? -1 : 1;
        const nullPlacement = metadataNullPlacement();
        return [...nodeIndices].sort((leftIndex, rightIndex) => {{
            const leftValue = metadataSortValue(leftIndex, state.metadataSort.columnKey);
            const rightValue = metadataSortValue(rightIndex, state.metadataSort.columnKey);
            const leftMissing = isMissingMetadataValue(leftValue);
            const rightMissing = isMissingMetadataValue(rightValue);
            if (leftMissing || rightMissing) {{
                if (leftMissing && rightMissing) {{
                    return nodeId(leftIndex).localeCompare(nodeId(rightIndex), undefined, {{numeric: true, sensitivity: 'base'}});
                }}
                return leftMissing
                    ? (nullPlacement === 'first' ? -1 : 1)
                    : (nullPlacement === 'first' ? 1 : -1);
            }}
            const comparison = compareMetadataValues(
                leftValue,
                rightValue,
            );
            if (comparison !== 0) {{
                return comparison * direction;
            }}
            return nodeId(leftIndex).localeCompare(nodeId(rightIndex), undefined, {{numeric: true, sensitivity: 'base'}});
        }});
    }}

    window.addEventListener('pointermove', event => {{
        if (!state.metadataResize) {{
            return;
        }}
        updateMetadataColumnWidth(
            state.metadataResize.columnKey,
            state.metadataResize.startWidth + (event.clientX - state.metadataResize.startX),
        );
    }});

    window.addEventListener('pointerup', () => {{
        if (!state.metadataResize) {{
            return;
        }}
        state.metadataResize = null;
        document.body.style.cursor = '';
        document.body.style.userSelect = '';
    }});

    window.addEventListener('pointercancel', () => {{
        if (!state.metadataResize) {{
            return;
        }}
        state.metadataResize = null;
        document.body.style.cursor = '';
        document.body.style.userSelect = '';
    }});

        document.getElementById('metadata-select-nodes').addEventListener('click', selectNodesFromMetadataRows);

    function toggleMetadataSort(columnKey) {{
        if (state.metadataSort.columnKey === columnKey) {{
            state.metadataSort.direction = state.metadataSort.direction === 'asc' ? 'desc' : 'asc';
        }} else {{
            state.metadataSort = {{columnKey, direction: 'asc'}};
        }}
        resetMetadataPage();
        updateMetadataTable();
    }}

    function resetMetadataSort() {{
        state.metadataSort = {{columnKey: null, direction: 'asc'}};
        resetMetadataPage();
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

    function toggleSelectionForNode(nodeIndex) {{
        if (state.selectedNodeIndices.has(nodeIndex)) {{
            state.selectedNodeIndices.delete(nodeIndex);
        }} else {{
            state.selectedNodeIndices.add(nodeIndex);
        }}
        renderClusterView();
        updateMetadataTable();
    }}

    function hitTestNodeAt(screenX, screenY) {{
        const worldPoint = screenToWorldPoint(screenX, screenY);
        for (let index = state.visibleLayout.length - 1; index >= 0; index--) {{
            const item = state.visibleLayout[index];
            const dx = worldPoint.x - item.x;
            const dy = worldPoint.y - item.y;
            if ((dx * dx) + (dy * dy) > item.radius * item.radius) {{
                continue;
            }}
            const component = state.bundle.graph.hierarchy.nodes[item.componentId];
            const dotLayout = componentDotLayout(component, item);
            const dotRadius = dotLayout.length > 0 ? dotLayout[0].radius : componentDotGeometry(component, item).dotRadius;
            const hitRadius = Math.max(dotRadius, 5 / Math.max(state.viewTransform.scale, 0.1));
            let nearestNodeIndex = null;
            let nearestDistanceSq = Infinity;
            for (const dot of dotLayout) {{
                const nodeDx = worldPoint.x - dot.x;
                const nodeDy = worldPoint.y - dot.y;
                const distanceSq = (nodeDx * nodeDx) + (nodeDy * nodeDy);
                if (distanceSq <= hitRadius * hitRadius && distanceSq < nearestDistanceSq) {{
                    nearestDistanceSq = distanceSq;
                    nearestNodeIndex = dot.memberIndex;
                }}
            }}
            if (nearestNodeIndex !== null) {{
                return nearestNodeIndex;
            }}
        }}
        return null;
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

    function selectComponentsInBox(rect, deselect) {{
        let changed = false;
        state.visibleLayout.forEach(item => {{
            const screenPoint = worldToScreenPoint(item.x, item.y);
            const screenRadius = item.radius * state.viewTransform.scale;
            if (!circleIntersectsRect(screenPoint.x, screenPoint.y, screenRadius, rect)) {{
                return;
            }}
            componentMembers(item.componentId).forEach(nodeIndex => {{
                if (deselect) {{
                    if (state.selectedNodeIndices.delete(nodeIndex)) {{
                        changed = true;
                    }}
                }} else if (!state.selectedNodeIndices.has(nodeIndex)) {{
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

    function selectNodesInBox(rect, deselect) {{
        let changed = false;
        state.visibleLayout.forEach(item => {{
            const screenPoint = worldToScreenPoint(item.x, item.y);
            const screenRadius = item.radius * state.viewTransform.scale;
            if (!circleIntersectsRect(screenPoint.x, screenPoint.y, screenRadius, rect)) {{
                return;
            }}
            const component = state.bundle.graph.hierarchy.nodes[item.componentId];
            const dotLayout = componentDotLayout(component, item);
            dotLayout.forEach(dot => {{
                const dotScreen = worldToScreenPoint(dot.x, dot.y);
                const dotRadius = dot.radius * state.viewTransform.scale;
                if (!circleIntersectsRect(dotScreen.x, dotScreen.y, dotRadius, rect)) {{
                    return;
                }}
                if (deselect) {{
                    if (state.selectedNodeIndices.delete(dot.memberIndex)) {{
                        changed = true;
                    }}
                }} else if (!state.selectedNodeIndices.has(dot.memberIndex)) {{
                    state.selectedNodeIndices.add(dot.memberIndex);
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
        if (event.ctrlKey || event.metaKey) {{
            const nodeIndex = hitTestNodeAt(point.x, point.y);
            if (nodeIndex !== null) {{
                toggleSelectionForNode(nodeIndex);
            }}
            return;
        }}
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
        scheduleClusterRender();
    }}, {{passive: false}});

    clusterCanvas.addEventListener('mousedown', event => {{
        if (!state.bundle) {{
            return;
        }}
        const point = canvasCoordinatesFromEvent(event);
        if (event.shiftKey) {{
            state.dragState = {{mode: 'select', nodeMode: event.ctrlKey || event.metaKey, deselect: event.altKey, startX: point.x, startY: point.y, endX: point.x, endY: point.y, moved: false}};
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
            scheduleClusterRender();
            return;
        }}
        state.dragState.endX = point.x;
        state.dragState.endY = point.y;
        state.dragState.moved = state.dragState.moved || Math.abs(point.x - state.dragState.startX) > 3 || Math.abs(point.y - state.dragState.startY) > 3;
        state.selectionBox = normalizedSelectionBox(state.dragState);
        scheduleClusterRender();
    }});

    window.addEventListener('mouseup', () => {{
        if (!state.dragState) {{
            return;
        }}
        state.suppressClick = Boolean(state.dragState.moved);
        if (state.dragState.mode === 'select' && state.selectionBox && (state.selectionBox.width > 4 || state.selectionBox.height > 4)) {{
            if (state.dragState.nodeMode) {{
                selectNodesInBox(state.selectionBox, state.dragState.deselect);
            }} else {{
                selectComponentsInBox(state.selectionBox, state.dragState.deselect);
            }}
        }}
        state.dragState = null;
        state.selectionBox = null;
        clusterCanvas.style.cursor = 'grab';
        scheduleClusterRender();
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
        scheduleThresholdUI(false);
    }});
    document.getElementById('threshold-slider').addEventListener('change', () => {{
        const stop = currentSliderStop();
        snapSliderToStop(stop);
        // Keep the user's current pan/zoom when the slider is released; only snap to the
        // nearest stop. (Use the Reset view button to re-fit.)
        scheduleThresholdUI(false);
    }});
    document.getElementById('threshold-input').addEventListener('change', jumpToThresholdValue);
    document.getElementById('threshold-input').addEventListener('keydown', event => {{
        if (event.key !== 'Enter') {{
            return;
        }}
        event.preventDefault();
        jumpToThresholdValue();
    }});
    document.getElementById('min-cluster-size').addEventListener('input', () => {{
        scheduleThresholdUI(true);
    }});
    document.getElementById('layout-algorithm').addEventListener('change', () => {{
        updateComponentSortButton();
        if (!state.bundle) {{
            return;
        }}
        drawClusterView(true);
    }});
    document.getElementById('render-cluster-bounds').addEventListener('change', renderClusterView);
    document.getElementById('render-nodes').addEventListener('change', renderClusterView);
    document.getElementById('leaf-pruning-only').addEventListener('change', () => {{
        scheduleThresholdUI(true);
    }});
    document.getElementById('color-by').addEventListener('change', () => {{ rebuildNodeColorCache(); renderClusterView(); }});
    document.getElementById('label-by').addEventListener('change', renderClusterView);
    document.getElementById('show-labels').addEventListener('change', renderClusterView);
    document.getElementById('show-node-counts').addEventListener('change', renderClusterView);
    document.getElementById('show-edge-scores').addEventListener('change', renderClusterView);
    document.getElementById('reduce-elongation').addEventListener('change', () => {{
        if (!state.bundle) {{
            return;
        }}
        renderClusterView();
    }});
    document.getElementById('sort-components-by-size').addEventListener('click', () => {{
        const button = document.getElementById('sort-components-by-size');
        const enabled = button.getAttribute('aria-pressed') === 'true';
        button.setAttribute('aria-pressed', enabled ? 'false' : 'true');
        updateComponentSortButton();
        if (!state.bundle) {{
            return;
        }}
        drawClusterView(true);
    }});
    document.getElementById('export-selected').addEventListener('click', exportSelection);
    document.getElementById('metadata-deselect-rows').addEventListener('click', () => {{
        clearMetadataRowSelection();
        applyMetadataTableRowHighlights();
    }});
    document.getElementById('metadata-reset-sort').addEventListener('click', resetMetadataSort);
    document.getElementById('metadata-filter').addEventListener('input', () => {{
        scheduleMetadataFilterUpdate();
    }});
    document.getElementById('metadata-null-order').addEventListener('change', () => {{
        resetMetadataPage();
        updateMetadataTable();
    }});
    document.getElementById('metadata-rows-per-page').addEventListener('change', () => {{
        resetMetadataPage();
        updateMetadataTable();
    }});
    document.getElementById('metadata-prev-page').addEventListener('click', () => {{
        stepMetadataPage(-1);
    }});
    document.getElementById('metadata-next-page').addEventListener('click', () => {{
        stepMetadataPage(1);
    }});
    document.getElementById('focus-largest-cluster').addEventListener('click', () => {{
        if (!state.bundle || state.visibleLayout.length === 0) {{ return; }}
        const hierarchyNodes = state.bundle.graph.hierarchy.nodes;
        let largestItem = state.visibleLayout[0];
        let largestSize = hierarchyNodes[largestItem.componentId].size;
        for (const item of state.visibleLayout) {{
            const size = hierarchyNodes[item.componentId].size;
            if (size > largestSize) {{ largestSize = size; largestItem = item; }}
        }}
        state.viewTransform.offsetX = clusterCanvas.width / 2 - largestItem.x * state.viewTransform.scale;
        state.viewTransform.offsetY = clusterCanvas.height / 2 - largestItem.y * state.viewTransform.scale;
        scheduleClusterRender();
    }});
    document.getElementById('reset-view').addEventListener('click', () => {{
        fitClusterViewToLayout();
        renderClusterView();
    }});
    document.getElementById('clear-selection').addEventListener('click', () => {{
        state.selectedNodeIndices = new Set();
        resetMetadataPage();
        clearMetadataRowSelection();
        renderClusterView();
        updateMetadataTable();
    }});
    window.addEventListener('resize', () => {{
        if (!state.bundle) {{
            return;
        }}
        drawClusterView(true);
    }});

    setupLayoutWorker();
    setupMetadataTableDelegation();
    warnIfUnsupported();
    updateComponentSortButton();
    drawSplitChart();
    drawClusterView();
    autoloadEmbeddedBundle();
</script>
</body>
</html>
"""


def write_ssn_viewer_html(
    out_path: str,
    title: str = "Domainator SSN Viewer",
    embedded_bundle_json: bytes | None = None,
) -> None:
    temp_path = make_temporary_output_path(out_path)
    try:
        with open(temp_path, "w", encoding="utf-8") as out_handle:
            out_handle.write(ssn_viewer_html(title=title, embedded_bundle_json=embedded_bundle_json))
        Path(temp_path).replace(out_path)
        temp_path = None
    finally:
        if temp_path is not None and Path(temp_path).exists():
            Path(temp_path).unlink()