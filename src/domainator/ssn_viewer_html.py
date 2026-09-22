"""Static HTML shell for the standalone Domainator Similarity Network Viewer."""

import base64
import json
from pathlib import Path

from domainator import __version__
from domainator.output_guardrails import make_temporary_output_path
from domainator.ssn_bundle import (
    SSN_VIEWER_BUNDLE_FORMAT,
    SSN_VIEWER_BUNDLE_VERSION,
    SUPPORTED_SSN_VIEWER_BUNDLE_VERSIONS,
)
from domainator.ssn_hierarchy import (
    MERGE_IMPACT_CHOICES,
    MERGE_IMPACT_MIN_CHILD,
    merge_impact_axis_labels,
)
from domainator.utils import NAMED_CATEGORICAL_PALETTES, OTHER_COLOR


VIEWER_APP_NAME = "Domainator Similarity Network Viewer"
# Titles that name no particular network, so the heading shows the app name alone
# rather than "<app name>: <app name>". The old app name is kept here so bundles
# written before the rename still get the bare heading.
_GENERIC_TITLES = frozenset({"", VIEWER_APP_NAME, "Domainator SSN Viewer"})


def viewer_heading(network_name: str | None) -> str:
    """The page's <h1>: the app name, plus the network name when there is one."""
    name = (network_name or "").strip()
    return VIEWER_APP_NAME if name in _GENERIC_TITLES else f"{VIEWER_APP_NAME}: {name}"


def _layout_core_js() -> str:
    """Layout math shared by the layout Web Worker and the main-thread fallback.

    Pure functions only -- no `state`, no DOM -- because this same source is
    compiled into the Worker, where neither exists. Plain (non-f) string, so the
    JavaScript uses normal single braces; it is interpolated both into the worker
    blob (via `_layout_worker_js`) and into the page's one `<script>` scope (via
    the `{layout_core_js}` placeholder). Function declarations hoist, so the
    order below is for reading, not for dependency.

    This used to exist twice -- once here in minified form and once inline in the
    page f-string -- and the two copies had already drifted. Keep it single-copy.
    """
    return """
function treeCenter(nodeIds, adjacency, hierarchyNodes) {
    if (nodeIds.length <= 2) {
        return nodeIds[0];
    }
    const degree = new Map();
    nodeIds.forEach(nodeId => {
        degree.set(nodeId, (adjacency.get(nodeId) || []).length);
    });
    let leaves = nodeIds.filter(nodeId => (degree.get(nodeId) || 0) <= 1);
    let remaining = nodeIds.length;
    while (remaining > 2 && leaves.length > 0) {
        remaining -= leaves.length;
        const nextLeaves = [];
        leaves.forEach(leafId => {
            (adjacency.get(leafId) || []).forEach(neighborId => {
                if (!degree.has(neighborId)) {
                    return;
                }
                degree.set(neighborId, degree.get(neighborId) - 1);
                if (degree.get(neighborId) === 1) {
                    nextLeaves.push(neighborId);
                }
            });
            degree.delete(leafId);
        });
        leaves = nextLeaves;
    }
    const candidates = degree.size > 0 ? Array.from(degree.keys()) : nodeIds;
    candidates.sort((leftId, rightId) => {
        const leftSize = hierarchyNodes[leftId].size;
        const rightSize = hierarchyNodes[rightId].size;
        return rightSize - leftSize || leftId - rightId;
    });
    return candidates[0];
}

function rootedTree(rootId, adjacency, hierarchyNodes) {
    const parent = new Map([[rootId, null]]);
    const order = [rootId];
    for (let index = 0; index < order.length; index++) {
        const nodeId = order[index];
        (adjacency.get(nodeId) || []).forEach(neighborId => {
            if (parent.has(neighborId)) {
                return;
            }
            parent.set(neighborId, nodeId);
            order.push(neighborId);
        });
    }

    const children = new Map();
    order.forEach(nodeId => children.set(nodeId, []));
    for (let index = 1; index < order.length; index++) {
        const nodeId = order[index];
        children.get(parent.get(nodeId)).push(nodeId);
    }
    children.forEach((childIds, nodeId) => {
        childIds.sort((leftId, rightId) => {
            const leftNode = hierarchyNodes[leftId];
            const rightNode = hierarchyNodes[rightId];
            return leftNode.leaf_start - rightNode.leaf_start || rightNode.size - leftNode.size || leftId - rightId;
        });
    });

    return {parent, order, children};
}

function componentRadiusForSize(size) {
    // Bubble radius scales so its area is proportional to the node count (every node
    // is always drawn).
    const areaPerNode = 48;
    return Math.sqrt((Math.max(1, size) * areaPerNode) / Math.PI);
}

function normalizeComponentLayout(items, padding = 20) {
    if (items.length === 0) {
        return {items: [], width: 0, height: 0};
    }
    let minX = Infinity;
    let minY = Infinity;
    let maxX = -Infinity;
    let maxY = -Infinity;
    items.forEach(item => {
        minX = Math.min(minX, item.x - item.radius);
        minY = Math.min(minY, item.y - item.radius);
        maxX = Math.max(maxX, item.x + item.radius);
        maxY = Math.max(maxY, item.y + item.radius);
    });
    const normalizedItems = items.map(item => ({
        ...item,
        x: item.x - minX + padding,
        y: item.y - minY + padding,
    }));
    return {
        items: normalizedItems,
        width: (maxX - minX) + (padding * 2),
        height: (maxY - minY) + (padding * 2),
    };
}

function componentLeafOrder(componentIds, hierarchyNodes) {
    return [...componentIds].sort((leftId, rightId) => {
        const leftNode = hierarchyNodes[leftId];
        const rightNode = hierarchyNodes[rightId];
        return leftNode.leaf_start - rightNode.leaf_start || rightNode.size - leftNode.size || leftId - rightId;
    });
}

function layoutLinkPairs(linksOrPairs) {
    return linksOrPairs.map(link => {
        if (Array.isArray(link)) {
            return {sourceId: link[0], targetId: link[1]};
        }
        return {sourceId: link.sourceId, targetId: link.targetId};
    });
}

function pointSegmentDistance(point, start, end) {
    const dx = end.x - start.x;
    const dy = end.y - start.y;
    const lengthSq = (dx * dx) + (dy * dy);
    if (lengthSq < 1e-9) {
        return {distance: Math.hypot(point.x - start.x, point.y - start.y), t: 0, closestX: start.x, closestY: start.y};
    }
    const t = Math.max(0, Math.min(1, (((point.x - start.x) * dx) + ((point.y - start.y) * dy)) / lengthSq));
    const closestX = start.x + (dx * t);
    const closestY = start.y + (dy * t);
    return {distance: Math.hypot(point.x - closestX, point.y - closestY), t, closestX, closestY};
}

function segmentOrientation(a, b, c) {
    return ((b.y - a.y) * (c.x - b.x)) - ((b.x - a.x) * (c.y - b.y));
}

function segmentsCross(a, b, c, d) {
    const o1 = segmentOrientation(a, b, c);
    const o2 = segmentOrientation(a, b, d);
    const o3 = segmentOrientation(c, d, a);
    const o4 = segmentOrientation(c, d, b);
    return (o1 * o2 < 0) && (o3 * o4 < 0);
}

function refineLayoutGeometry(items, linksOrPairs, options = {}) {
    const refined = items.map(item => ({...item}));
    if (refined.length <= 1) {
        return refined;
    }

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

    if (pairChecks <= maxPairChecks) {
        for (let iteration = 0; iteration < overlapIterations; iteration++) {
            for (let leftIndex = 0; leftIndex < refined.length; leftIndex++) {
                const left = refined[leftIndex];
                for (let rightIndex = leftIndex + 1; rightIndex < refined.length; rightIndex++) {
                    const right = refined[rightIndex];
                    let dx = right.x - left.x;
                    let dy = right.y - left.y;
                    let distance = Math.hypot(dx, dy);
                    if (distance < 1e-6) {
                        dx = (seededUnit(left.componentId + right.componentId, iteration + 21) - 0.5) || 0.01;
                        dy = (seededUnit(left.componentId + right.componentId, iteration + 22) - 0.5) || 0.01;
                        distance = Math.hypot(dx, dy);
                    }
                    const minimumDistance = left.radius + right.radius + bubblePadding;
                    if (distance >= minimumDistance) {
                        continue;
                    }
                    const shift = ((minimumDistance - distance) / 2) * 0.72;
                    const shiftX = (dx / distance) * shift;
                    const shiftY = (dy / distance) * shift;
                    left.x -= shiftX;
                    left.y -= shiftY;
                    right.x += shiftX;
                    right.y += shiftY;
                }
            }
        }
    }

    if (links.length * refined.length <= maxEdgeNodeChecks) {
        for (let iteration = 0; iteration < edgeIterations; iteration++) {
            links.forEach(link => {
                const source = itemById.get(link.sourceId);
                const target = itemById.get(link.targetId);
                if (!source || !target) {
                    return;
                }
                refined.forEach(item => {
                    if (item.componentId === link.sourceId || item.componentId === link.targetId) {
                        return;
                    }
                    const hit = pointSegmentDistance(item, source, target);
                    if (hit.t <= 0.03 || hit.t >= 0.97) {
                        return;
                    }
                    const minimumDistance = item.radius + edgePadding;
                    if (hit.distance >= minimumDistance) {
                        return;
                    }
                    let normalX = item.x - hit.closestX;
                    let normalY = item.y - hit.closestY;
                    let normalLength = Math.hypot(normalX, normalY);
                    if (normalLength < 1e-6) {
                        const edgeDx = target.x - source.x;
                        const edgeDy = target.y - source.y;
                        normalX = -edgeDy || 1;
                        normalY = edgeDx || 0;
                        normalLength = Math.hypot(normalX, normalY);
                    }
                    const push = (minimumDistance - hit.distance) * 0.68;
                    item.x += (normalX / normalLength) * push;
                    item.y += (normalY / normalLength) * push;
                });
            });
        }
    }

    const crossingChecks = (links.length * (links.length - 1)) / 2;
    if (crossingChecks <= maxCrossingChecks) {
        for (let iteration = 0; iteration < crossingIterations; iteration++) {
            for (let leftIndex = 0; leftIndex < links.length; leftIndex++) {
                const leftLink = links[leftIndex];
                const a = itemById.get(leftLink.sourceId);
                const b = itemById.get(leftLink.targetId);
                if (!a || !b) {
                    continue;
                }
                for (let rightIndex = leftIndex + 1; rightIndex < links.length; rightIndex++) {
                    const rightLink = links[rightIndex];
                    if (leftLink.sourceId === rightLink.sourceId || leftLink.sourceId === rightLink.targetId || leftLink.targetId === rightLink.sourceId || leftLink.targetId === rightLink.targetId) {
                        continue;
                    }
                    const c = itemById.get(rightLink.sourceId);
                    const d = itemById.get(rightLink.targetId);
                    if (!c || !d || !segmentsCross(a, b, c, d)) {
                        continue;
                    }
                    let dx = b.x - a.x;
                    let dy = b.y - a.y;
                    let length = Math.hypot(dx, dy);
                    if (length < 1e-6) {
                        continue;
                    }
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
                }
            }
        }
    }

    if (links.length * refined.length <= maxEdgeNodeChecks) {
        links.forEach(link => {
            const source = itemById.get(link.sourceId);
            const target = itemById.get(link.targetId);
            if (!source || !target) {
                return;
            }
            refined.forEach(item => {
                if (item.componentId === link.sourceId || item.componentId === link.targetId) {
                    return;
                }
                const hit = pointSegmentDistance(item, source, target);
                if (hit.t <= 0.03 || hit.t >= 0.97) {
                    return;
                }
                const minimumDistance = item.radius + edgePadding;
                if (hit.distance >= minimumDistance) {
                    return;
                }
                let normalX = item.x - hit.closestX;
                let normalY = item.y - hit.closestY;
                let normalLength = Math.hypot(normalX, normalY);
                if (normalLength < 1e-6) {
                    const edgeDx = target.x - source.x;
                    const edgeDy = target.y - source.y;
                    normalX = -edgeDy || 1;
                    normalY = edgeDx || 0;
                    normalLength = Math.hypot(normalX, normalY);
                }
                const push = (minimumDistance - hit.distance) * 0.82;
                item.x += (normalX / normalLength) * push;
                item.y += (normalY / normalLength) * push;
            });
        });
    }

    if (pairChecks <= maxPairChecks) {
        for (let iteration = 0; iteration < 4; iteration++) {
            for (let leftIndex = 0; leftIndex < refined.length; leftIndex++) {
                const left = refined[leftIndex];
                for (let rightIndex = leftIndex + 1; rightIndex < refined.length; rightIndex++) {
                    const right = refined[rightIndex];
                    let dx = right.x - left.x;
                    let dy = right.y - left.y;
                    let distance = Math.hypot(dx, dy);
                    if (distance < 1e-6) {
                        dx = (seededUnit(left.componentId + right.componentId, iteration + 41) - 0.5) || 0.01;
                        dy = (seededUnit(left.componentId + right.componentId, iteration + 42) - 0.5) || 0.01;
                        distance = Math.hypot(dx, dy);
                    }
                    const minimumDistance = left.radius + right.radius + bubblePadding;
                    if (distance >= minimumDistance) {
                        continue;
                    }
                    const shift = ((minimumDistance - distance) / 2) * 0.9;
                    const shiftX = (dx / distance) * shift;
                    const shiftY = (dy / distance) * shift;
                    left.x -= shiftX;
                    left.y -= shiftY;
                    right.x += shiftX;
                    right.y += shiftY;
                }
            }
        }
    }

    return refined;
}

function seededUnit(componentId, salt) {
    const raw = Math.sin((componentId + 1) * 12.9898 + salt * 78.233) * 43758.5453;
    return raw - Math.floor(raw);
}

// Radial tree layout seed: root (tree center) at origin, children placed in concentric rings
// with angular wedges allocated per subtree leaf count. Crossing-free and near-circular (compact
// aspect), so it both fixes the tall linear-tidy seed and gives Force its radial branch shape.
// Clearance between a parent bubble's rim and its child's in the radial seed.
const RADIAL_SEED_RING_PAD = 36;

function radialTreeSeed(componentIds, adjacency, hierarchyNodes, ringGap) {
    const rootId = treeCenter(componentIds, adjacency, hierarchyNodes);
    const tree = rootedTree(rootId, adjacency, hierarchyNodes);
    const radii = new Map(componentIds.map(id => [id, componentRadiusForSize(hierarchyNodes[id].size)]));
    // Each ring sits one parent-plus-child apart from the ring inside it, rather than every
    // ring being spaced by the largest bubble anywhere in the component. One big cluster used
    // to push every ring out, including rings joining two singletons, and anchorStrength then
    // held the simulation in that spread-out arrangement.
    const gap = Math.max(ringGap || 0, 60);
    const ringRadius = new Map([[rootId, 0]]);
    tree.order.forEach(nodeId => {
        (tree.children.get(nodeId) || []).forEach(c => {
            const step = Math.max(gap, (radii.get(nodeId) || 10) + (radii.get(c) || 10) + RADIAL_SEED_RING_PAD);
            ringRadius.set(c, (ringRadius.get(nodeId) || 0) + step);
        });
    });
    const leafCount = new Map();
    [...tree.order].reverse().forEach(id => {
        const ch = tree.children.get(id) || [];
        leafCount.set(id, ch.length === 0 ? 1 : ch.reduce((s, c) => s + (leafCount.get(c) || 1), 0));
    });
    const positionById = new Map();
    positionById.set(rootId, {componentId: rootId, x: 0, y: 0, radius: radii.get(rootId) || 10});
    const wedge = new Map([[rootId, [0, Math.PI * 2]]]);
    tree.order.forEach(id => {
        const span = wedge.get(id) || [0, Math.PI * 2];
        const a0 = span[0], a1 = span[1];
        const ch = tree.children.get(id) || [];
        const total = ch.reduce((s, c) => s + (leafCount.get(c) || 1), 0) || 1;
        let cursor = a0;
        ch.forEach(c => {
            const ca0 = cursor, ca1 = cursor + (a1 - a0) * ((leafCount.get(c) || 1) / total);
            cursor = ca1;
            wedge.set(c, [ca0, ca1]);
            const ang = (ca0 + ca1) / 2, r = ringRadius.get(c) || 0;
            positionById.set(c, {componentId: c, x: Math.cos(ang) * r + (seededUnit(c, 11) - 0.5) * 0.5, y: Math.sin(ang) * r + (seededUnit(c, 12) - 0.5) * 0.5, radius: radii.get(c) || 10});
        });
    });
    return {positionById, order: tree.order};
}

// --- Barnes-Hut quadtree: O(n log n) repulsion so physics scales to large forests ---
function bhBuild(ids, positions, radii) {
    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    for (const id of ids) { const p = positions.get(id); if (p.x < minX) minX = p.x; if (p.x > maxX) maxX = p.x; if (p.y < minY) minY = p.y; if (p.y > maxY) maxY = p.y; }
    if (!isFinite(minX)) { minX = 0; minY = 0; maxX = 1; maxY = 1; }
    const size = Math.max(maxX - minX, maxY - minY, 1) * 1.0001;
    const root = {x0: minX, y0: minY, size, mass: 0, comX: 0, comY: 0, maxR: 0, id: null, px: 0, py: 0, r: 0, children: null, bucket: false};
    for (const id of ids) { const p = positions.get(id); bhInsert(root, id, p.x, p.y, radii.get(id) || 10, 0); }
    return root;
}
function bhInsert(cell, id, px, py, r, depth) {
    cell.mass += 1; cell.comX += px; cell.comY += py; if (r > cell.maxR) cell.maxR = r;
    if (cell.children === null) {
        if (cell.id === null) { cell.id = id; cell.px = px; cell.py = py; cell.r = r; return; }
        if (depth >= 48 || cell.size < 1e-6) { cell.bucket = true; return; }
        const eid = cell.id, epx = cell.px, epy = cell.py, er = cell.r;
        cell.id = null; cell.children = [null, null, null, null];
        bhInsertChild(cell, eid, epx, epy, er, depth);
        bhInsertChild(cell, id, px, py, r, depth);
        return;
    }
    bhInsertChild(cell, id, px, py, r, depth);
}
function bhInsertChild(cell, id, px, py, r, depth) {
    const half = cell.size / 2;
    const qx = px >= cell.x0 + half ? 1 : 0, qy = py >= cell.y0 + half ? 1 : 0, qi = qy * 2 + qx;
    let child = cell.children[qi];
    if (!child) { child = {x0: cell.x0 + qx * half, y0: cell.y0 + qy * half, size: half, mass: 0, comX: 0, comY: 0, maxR: 0, id: null, px: 0, py: 0, r: 0, children: null, bucket: false}; cell.children[qi] = child; }
    bhInsert(child, id, px, py, r, depth + 1);
}
function bhAccumulate(cell, selfId, px, py, selfR, acc, p) {
    if (cell.mass === 0) return;
    if (cell.children === null) {
        if (cell.bucket) {
            const comX = cell.comX / cell.mass, comY = cell.comY / cell.mass;
            let dx = px - comX, dy = py - comY, distSq = dx * dx + dy * dy;
            if (distSq < 1) return;
            const dist = Math.sqrt(distSq), rf = p.repulsion * cell.mass / Math.max(distSq, 16);
            acc.x += dx / dist * rf; acc.y += dy / dist * rf;
            return;
        }
        if (cell.id === null || cell.id === selfId) return;
        let dx = px - cell.px, dy = py - cell.py, distSq = dx * dx + dy * dy;
        if (distSq < 1e-6) { dx = (seededUnit(selfId + cell.id, p.iter + 7) - 0.5) * 0.01; dy = (seededUnit(selfId + cell.id, p.iter + 8) - 0.5) * 0.01; distSq = dx * dx + dy * dy; }
        const dist = Math.sqrt(distSq), rf = p.repulsion / Math.max(distSq, 16);
        const ov = selfR + cell.r + 12 - dist;
        const cf = ov > 0 ? ov * p.collisionStrength : 0;
        acc.x += dx / dist * (rf + cf); acc.y += dy / dist * (rf + cf);
        return;
    }
    const comX = cell.comX / cell.mass, comY = cell.comY / cell.mass;
    let dx = px - comX, dy = py - comY, distSq = dx * dx + dy * dy;
    if (distSq < 1e-9) distSq = 1e-9;
    // Always descend cells whose nearest edge is within collision range of self, so the radii-aware
    // collision is exact in the near field even when theta would aggregate.
    const nx = px < cell.x0 ? cell.x0 : (px > cell.x0 + cell.size ? cell.x0 + cell.size : px);
    const ny = py < cell.y0 ? cell.y0 : (py > cell.y0 + cell.size ? cell.y0 + cell.size : py);
    const bdx = px - nx, bdy = py - ny, colR = selfR + cell.maxR + 12;
    const near = (bdx * bdx + bdy * bdy) < colR * colR;
    if (!near && cell.size * cell.size < p.thetaSq * distSq) {
        const dist = Math.sqrt(distSq), rf = p.repulsion * cell.mass / Math.max(distSq, 16);
        acc.x += dx / dist * rf; acc.y += dy / dist * rf;
        return;
    }
    for (let i = 0; i < 4; i++) { const c = cell.children[i]; if (c) bhAccumulate(c, selfId, px, py, selfR, acc, p); }
}
// Collect the separation push needed to lift `self` out of any bubble it overlaps (local descent).
function bhCollect(cell, selfId, px, py, selfR, d, pad) {
    if (cell.mass === 0) return;
    if (cell.children === null) {
        if (cell.bucket || cell.id === null || cell.id === selfId) return;
        let dx = px - cell.px, dy = py - cell.py, distSq = dx * dx + dy * dy;
        const minD = selfR + cell.r + pad;
        if (distSq >= minD * minD) return;
        let dist = Math.sqrt(distSq);
        if (dist < 1e-6) { dx = (seededUnit(selfId + cell.id, 7) - 0.5) || 0.01; dy = (seededUnit(selfId + cell.id, 8) - 0.5) || 0.01; dist = Math.hypot(dx, dy); }
        const push = (minD - dist) * 0.5;
        d.x += dx / dist * push; d.y += dy / dist * push;
        return;
    }
    const nx = px < cell.x0 ? cell.x0 : (px > cell.x0 + cell.size ? cell.x0 + cell.size : px);
    const ny = py < cell.y0 ? cell.y0 : (py > cell.y0 + cell.size ? cell.y0 + cell.size : py);
    const bdx = px - nx, bdy = py - ny, colR = selfR + cell.maxR + pad;
    if (bdx * bdx + bdy * bdy > colR * colR) return;
    for (let i = 0; i < 4; i++) { const c = cell.children[i]; if (c) bhCollect(c, selfId, px, py, selfR, d, pad); }
}
// Iteratively separate overlapping bubbles using the quadtree (replaces the O(n^2) overlap pass on
// large components, where refineLayoutGeometry is skipped). Mutates positions in place.
function bhResolveOverlaps(ids, positions, radii, passes, pad) {
    for (let pass = 0; pass < passes; pass++) {
        const tree = bhBuild(ids, positions, radii);
        const disp = new Map(ids.map(id => [id, {x: 0, y: 0}]));
        ids.forEach(id => { const pos = positions.get(id); bhCollect(tree, id, pos.x, pos.y, radii.get(id) || 10, disp.get(id), pad); });
        let maxMoved = 0;
        ids.forEach(id => { const d = disp.get(id), pos = positions.get(id); pos.x += d.x; pos.y += d.y; maxMoved = Math.max(maxMoved, Math.abs(d.x) + Math.abs(d.y)); });
        if (maxMoved < 1) break;
    }
}

function simulateComponentLayout(componentIds, adjacency, options = {}, hierarchyNodes) {
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
    const collisionStrength = options.collisionStrength ?? 0.28;
    const maxStep = options.maxStep ?? 48;
    // Rim-to-rim clearance a relaxed edge aims for. It is 2x bubblePadding so the spring and
    // refineLayoutGeometry's overlapPass (which separates to rL + rR + bubblePadding) are not
    // pulling against each other, which also guarantees the drawn edge survives that pass.
    const linkGap = options.linkGap ?? 34;
    // Floor, so two singletons stay visually distinct rather than fusing into one blob.
    const minEdgeLength = options.minEdgeLength ?? 46;
    // Minimum spacing between rings of the radial seed. This is the single biggest control on
    // how spread out the finished drawing is, because anchorStrength holds the simulation near
    // the seed. Tightening it compacts the drawing but eventually costs edge crossings: on the
    // reference networks, 60 roughly doubles the crossing count while 90 nearly halves it.
    const seedRingGap = options.seedRingGap ?? 90;
    // Barnes-Hut keeps repulsion O(n log n), so physics now runs on large components too; only
    // truly huge ones fall back to the (already compact) radial seed.
    const maxPhysicsNodes = options.maxPhysicsNodes ?? 50000;
    const theta = options.theta ?? 0.9;
    const thetaSq = theta * theta;

    // Seed positions AND anchors from the crossing-free radial tree layout (compact + circular),
    // centered on the origin so gravity pulls toward the middle.
    const seed = radialTreeSeed(componentIds, adjacency, hierarchyNodes, seedRingGap);
    let centerX = 0;
    let centerY = 0;
    componentIds.forEach(nodeId => {
        const s = seed.positionById.get(nodeId) || {x: 0, y: 0};
        centerX += s.x;
        centerY += s.y;
    });
    centerX /= Math.max(1, componentIds.length);
    centerY /= Math.max(1, componentIds.length);
    componentIds.forEach(nodeId => {
        const s = seed.positionById.get(nodeId) || {x: 0, y: 0};
        const anchorX = s.x - centerX;
        const anchorY = s.y - centerY;
        anchors.set(nodeId, {x: anchorX, y: anchorY});
        positions.set(nodeId, {x: anchorX, y: anchorY});
        velocities.set(nodeId, {x: 0, y: 0});
    });

    const edgePairs = [];
    componentIds.forEach(nodeId => {
        (adjacency.get(nodeId) || []).forEach(neighborId => {
            if (nodeId < neighborId) {
                edgePairs.push([nodeId, neighborId]);
            }
        });
    });

    if (componentIds.length <= maxPhysicsNodes) {
    for (let iteration = 0; iteration < iterations; iteration++) {
        const forces = new Map(componentIds.map(nodeId => [nodeId, {x: 0, y: 0}]));
        // Repulsion + near-field collision via Barnes-Hut (O(n log n) instead of O(n^2)). The
        // distance floor (max(distSq, 16)) keeps a near-coincident pair from exploding.
        const tree = bhBuild(componentIds, positions, radii);
        const bhParams = {repulsion, collisionStrength, thetaSq, iter: iteration};
        componentIds.forEach(nodeId => {
            const pos = positions.get(nodeId);
            bhAccumulate(tree, nodeId, pos.x, pos.y, radii.get(nodeId) || 10, forces.get(nodeId), bhParams);
        });

        edgePairs.forEach(([leftId, rightId]) => {
            const leftPos = positions.get(leftId);
            const rightPos = positions.get(rightId);
            let dx = rightPos.x - leftPos.x;
            let dy = rightPos.y - leftPos.y;
            let dist = Math.hypot(dx, dy);
            if (dist < 1e-6) {
                dist = 1e-6;
                dx = minEdgeLength;
                dy = 0;
            }
            // Rest length keeps connected neighbors outside each other's bubbles, and is set
            // by the two bubbles actually being joined. It used to carry a flat 124px floor,
            // which is what made a pair of 4px singletons settle 124px apart.
            const rest = Math.max(minEdgeLength, (radii.get(leftId) || 0) + (radii.get(rightId) || 0) + linkGap);
            const delta = dist - rest;
            const force = spring * delta;
            const forceX = (dx / dist) * force;
            const forceY = (dy / dist) * force;
            forces.get(leftId).x += forceX;
            forces.get(leftId).y += forceY;
            forces.get(rightId).x -= forceX;
            forces.get(rightId).y -= forceY;
        });

        let totalMovement = 0;
        componentIds.forEach(nodeId => {
            const position = positions.get(nodeId);
            const velocity = velocities.get(nodeId);
            const force = forces.get(nodeId);
            const anchor = anchors.get(nodeId);
            // Anchor and gravity are summed in one expression on purpose. Splitting them
            // into two `+=` steps is the same arithmetic but not the same floating point
            // (addition is not associative), and over ~90 iterations of a chaotic
            // simulation that reshuffles the whole layout. This grouping is the one the
            // Web Worker has always used, i.e. the layout people actually see.
            force.x += (anchor.x - position.x) * anchorStrength - position.x * gravity;
            force.y += (anchor.y - position.y) * anchorStrength - position.y * gravity;
            velocity.x = (velocity.x + force.x) * damping;
            velocity.y = (velocity.y + force.y) * damping;
            // Cap per-step displacement to keep the simulation numerically stable.
            const speed = Math.hypot(velocity.x, velocity.y);
            if (speed > maxStep) {
                velocity.x *= maxStep / speed;
                velocity.y *= maxStep / speed;
            }
            position.x += velocity.x;
            position.y += velocity.y;
            totalMovement += Math.abs(velocity.x) + Math.abs(velocity.y);
        });
        // Early-exit once the system has settled (avoids the worker hanging on big forests).
        if (iteration > 12 && (totalMovement / componentIds.length) < 0.05) {
            break;
        }
    }
    // Large components skip the O(n^2) refine overlap pass; clean residual overlaps via the quadtree.
    if (componentIds.length > 600) { bhResolveOverlaps(componentIds, positions, radii, 14, options.bubblePadding ?? 14); }
    }

    const rawItems = componentIds.map(nodeId => {
        const position = positions.get(nodeId);
        return {componentId: nodeId, x: position.x, y: position.y, radius: radii.get(nodeId) || 10};
    });
    const refinedItems = refineLayoutGeometry(rawItems, edgePairs, {
        bubblePadding: options.bubblePadding ?? 14,
        edgePadding: options.edgePadding ?? 8,
        overlapIterations: options.geometryIterations ?? 5,
        edgeIterations: options.edgeIterations ?? 3,
        crossingIterations: options.crossingIterations ?? 2,
    });
    return normalizeComponentLayout(refinedItems, 24);
}

function clusterGraphComponents(visibleIds, links, hierarchyNodes, sortBySizeEnabled) {
    const adjacency = new Map();
    visibleIds.forEach(nodeId => adjacency.set(nodeId, []));
    links.forEach(link => {
        adjacency.get(link.sourceId)?.push(link.targetId);
        adjacency.get(link.targetId)?.push(link.sourceId);
    });

    const components = [];
    const seen = new Set();
    visibleIds.forEach(nodeId => {
        if (seen.has(nodeId)) {
            return;
        }
        const stack = [nodeId];
        const componentIds = [];
        seen.add(nodeId);
        while (stack.length > 0) {
            const currentId = stack.pop();
            componentIds.push(currentId);
            (adjacency.get(currentId) || []).forEach(neighborId => {
                if (seen.has(neighborId)) {
                    return;
                }
                seen.add(neighborId);
                stack.push(neighborId);
            });
        }
        components.push(componentIds);
    });

    components.sort((leftIds, rightIds) => {
        const leftNodeCount = leftIds.reduce((sum, nodeId) => sum + hierarchyNodes[nodeId].size, 0);
        const rightNodeCount = rightIds.reduce((sum, nodeId) => sum + hierarchyNodes[nodeId].size, 0);
        const leftStart = Math.min(...leftIds.map(nodeId => hierarchyNodes[nodeId].leaf_start));
        const rightStart = Math.min(...rightIds.map(nodeId => hierarchyNodes[nodeId].leaf_start));
        if (sortBySizeEnabled) {
            return rightNodeCount - leftNodeCount || leftStart - rightStart || rightIds.length - leftIds.length;
        }
        return leftStart - rightStart || rightIds.length - leftIds.length;
    });
    return {adjacency, components};
}

function packLayouts(componentLayouts, options = {}) {
    const gapX = options.gapX ?? 120;
    const gapY = options.gapY ?? 120;
    const outerPadding = options.outerPadding ?? 72;
    const valid = componentLayouts.filter(component => component && component.items.length > 0);
    // Adaptive near-square arrangement when no explicit width is given: keeps many disconnected
    // single-cluster components from spreading into a wide sparse grid (which would make
    // fit-to-view collapse to invisible specks).
    let rowTargetWidth = options.rowTargetWidth;
    if (rowTargetWidth == null) {
        const totalArea = valid.reduce((sum, component) => sum + (component.width * component.height), 0);
        const widest = valid.reduce((maxWidth, component) => Math.max(maxWidth, component.width), 0);
        rowTargetWidth = Math.max(widest, Math.sqrt(totalArea) * 1.3) + outerPadding;
    }
    const packed = [];
    let cursorX = outerPadding;
    let cursorY = outerPadding;
    let rowHeight = 0;

    componentLayouts.forEach(component => {
        if (!component || component.items.length === 0) {
            return;
        }
        if (cursorX > outerPadding && cursorX + component.width > rowTargetWidth) {
            cursorX = outerPadding;
            cursorY += rowHeight + gapY;
            rowHeight = 0;
        }
        component.items.forEach(item => {
            packed.push({
                componentId: item.componentId,
                x: item.x + cursorX,
                y: item.y + cursorY,
                radius: item.radius,
            });
        });
        cursorX += component.width + gapX;
        rowHeight = Math.max(rowHeight, component.height);
    });
    return packed;
}

// Tidy (Reingold-Tilford style) layout for a single component, rooted at its tree center.
// Returns crossing-free positions for forest topology even with long branches. Shared by the
// Tree layout and used to seed the Force-directed simulation. O(n), x offset by originX.
// Tree layout spacing. TREE_MIN_LEVEL_GAP is the floor on the horizontal distance between
// adjacent columns; it also decides how many edge score labels have room to be drawn.
const TREE_MIN_LEVEL_GAP = 72;
const TREE_LEVEL_PAD = 48;
const TREE_SIBLING_GAP = 20;
const MIN_TREE_NODE_RADIUS = 10;

function tidyComponentLayout(componentIds, adjacency, hierarchyNodes, originX) {
    const rootId = treeCenter(componentIds, adjacency, hierarchyNodes);
    const tree = rootedTree(rootId, adjacency, hierarchyNodes);
    const depthByNode = new Map([[rootId, 0]]);
    tree.order.forEach(nodeId => {
        (tree.children.get(nodeId) || []).forEach(childId => {
            depthByNode.set(childId, (depthByNode.get(nodeId) || 0) + 1);
        });
    });
    const radii = new Map(componentIds.map(nodeId => [nodeId, componentRadiusForSize(hierarchyNodes[nodeId].size)]));
    const maxDepth = componentIds.reduce((maxValue, nodeId) => Math.max(maxValue, depthByNode.get(nodeId) || 0), 0);
    // Each column is placed relative to the one before it, spaced by the widest bubbles on
    // those two depths. Using one component-wide levelGap meant a single large cluster set the
    // column spacing for the whole tree, including columns joining two 4px singletons -- on the
    // reference networks that was a 288px gap between bubbles 8px across.
    const maxRadiusAtDepth = new Array(maxDepth + 1).fill(MIN_TREE_NODE_RADIUS);
    componentIds.forEach(nodeId => {
        const depth = depthByNode.get(nodeId) || 0;
        maxRadiusAtDepth[depth] = Math.max(maxRadiusAtDepth[depth], radii.get(nodeId) || MIN_TREE_NODE_RADIUS);
    });
    const columnX = [0];
    for (let depth = 0; depth < maxDepth; depth++) {
        columnX.push(columnX[depth] + Math.max(TREE_MIN_LEVEL_GAP,
            maxRadiusAtDepth[depth] + maxRadiusAtDepth[depth + 1] + TREE_LEVEL_PAD));
    }
    // subtreeSpan already reserves each node's own diameter, so this only has to keep adjacent
    // subtrees apart; it does not need to scale with the biggest bubble in the component.
    const siblingGap = TREE_SIBLING_GAP;
    const subtreeSpan = new Map();
    [...tree.order].reverse().forEach(nodeId => {
        const childIds = tree.children.get(nodeId) || [];
        const nodeSpan = ((radii.get(nodeId) || 12) * 2) + 18;
        if (childIds.length === 0) {
            subtreeSpan.set(nodeId, nodeSpan);
            return;
        }
        const childrenSpan = childIds.reduce((sum, childId) => sum + subtreeSpan.get(childId), 0) + (Math.max(0, childIds.length - 1) * siblingGap);
        subtreeSpan.set(nodeId, Math.max(nodeSpan, childrenSpan));
    });

    const positionById = new Map();
    function placeNode(nodeId, topY) {
        const nodeSpan = subtreeSpan.get(nodeId) || 26;
        const childIds = tree.children.get(nodeId) || [];
        const x = originX + columnX[depthByNode.get(nodeId) || 0];
        if (childIds.length === 0) {
            positionById.set(nodeId, {componentId: nodeId, x, y: topY + (nodeSpan / 2), radius: radii.get(nodeId) || 10});
            return;
        }
        const childrenSpan = childIds.reduce((sum, childId) => sum + subtreeSpan.get(childId), 0) + (Math.max(0, childIds.length - 1) * siblingGap);
        let childCursor = topY + Math.max(0, (nodeSpan - childrenSpan) / 2);
        const childCenters = [];
        childIds.forEach(childId => {
            placeNode(childId, childCursor);
            childCenters.push(positionById.get(childId).y);
            childCursor += subtreeSpan.get(childId) + siblingGap;
        });
        const centerY = childCenters.reduce((sum, value) => sum + value, 0) / Math.max(1, childCenters.length);
        positionById.set(nodeId, {componentId: nodeId, x, y: centerY, radius: radii.get(nodeId) || 10});
    }

    placeNode(rootId, 0);
    // tidyForestLayout advances its cursor by this, so it has to track columnX or the gap
    // between components grows to fill the space the columns no longer take.
    const width = columnX[maxDepth] + (maxRadiusAtDepth[maxDepth] * 2);
    return {positionById, order: tree.order, width};
}

// Force-layout tuning, shared by the Web Worker and the main-thread fallback so the
// two paths cannot drift into producing different layouts for the same network.
function forestForceOptions(componentSize) {
    return {
        // repulsion and anchorStrength are the two dominant levers on how spread out the
        // drawing is; both were tuned against real networks. Loosening them further compacts
        // more but starts trading edge crossings for it, sharply so below anchorStrength 0.01.
        repulsion: 1000,
        seedRingGap: 90,
        spring: 0.08,
        gravity: 0.004,
        damping: 0.85,
        iterations: Math.max(30, Math.min(90, Math.round(850000 / componentSize))),
        collisionStrength: 0.5,
        bubblePadding: 17,
        edgePadding: 10,
        geometryIterations: 7,
        edgeIterations: 4,
        crossingIterations: 4,
        anchorStrength: 0.02,
        maxPhysicsNodes: 50000,
        theta: 1.5,
    };
}
const FOREST_PACK_OPTIONS = {gapX: 36, gapY: 36, outerPadding: 72};
"""


def _layout_worker_js() -> str:
    """Standalone JS for the layout Web Worker: the shared core plus its entry point."""
    return _layout_core_js() + """
self.onmessage = function(event) {
    const msg = event.data;
    if (msg.type !== 'computeLayout') { return; }
    const {requestId, key, algorithm, visibleIds, links, hierarchyNodes, sortBySizeEnabled} = msg;
    const {adjacency, components} = clusterGraphComponents(visibleIds, links, hierarchyNodes, sortBySizeEnabled);
    const componentLayouts = components.map(ids =>
        simulateComponentLayout(ids, adjacency, forestForceOptions(ids.length), hierarchyNodes));
    const layout = packLayouts(componentLayouts, FOREST_PACK_OPTIONS);
    self.postMessage({requestId, key, layout});
};
"""


def _session_state_js() -> str:
    """Session save/load: the `app_state` section of a v4 bundle.

    Plain (non-f) string like `_layout_worker_js`, so the JavaScript below uses
    normal single braces. It is interpolated into the viewer's one `<script>`
    scope, so it can call `state`, `renderClusterView()` and friends directly,
    and its function declarations hoist above the listener wiring.
    """
    return r"""
    // ---------------------------------------------------------------------
    // Session state
    //
    // A saved session is an ordinary bundle whose `metadata` carries whatever
    // the user edited plus one extra top-level section, `app_state`. Because
    // the viewer's UI is expected to churn, state is *not* dumped ad hoc:
    // every persisted control is one entry in VIEW_STATE_FIELDS below.
    // Serializing walks the registry; applying looks each saved key up in it.
    // That makes both directions tolerant -- a key this build dropped is
    // skipped and reported, and a key the file omits keeps its default -- so
    // old files open in new viewers and new files open in old ones.
    // ---------------------------------------------------------------------
    const SESSION_STATE_VERSION = 1;
    const PRESET_SLOT_COUNT = 10;

    function domValueField(section, key, elementId, options) {
        const parse = (options && options.parse) || (value => value);
        const format = (options && options.format) || (value => String(value));
        return {
            section,
            key,
            get: () => parse(document.getElementById(elementId).value),
            set: value => {
                const element = document.getElementById(elementId);
                const text = format(value);
                if (element.tagName === 'SELECT' &&
                        !Array.from(element.options).some(option => option.value === text)) {
                    // Restoring a column or option this bundle does not have would
                    // silently blank the control; report it as skipped instead.
                    throw new Error('no option "' + text + '"');
                }
                element.value = text;
            },
        };
    }

    function domCheckboxField(section, key, elementId) {
        return {
            section,
            key,
            get: () => document.getElementById(elementId).checked,
            set: value => { document.getElementById(elementId).checked = Boolean(value); },
        };
    }

    function nodeIndexByIdMap() {
        if (!state.nodeIndexById) {
            state.nodeIndexById = new Map(state.bundle.graph.nodes.map((id, index) => [id, index]));
        }
        return state.nodeIndexById;
    }

    // Selections persist as node ids, never indices: an index is a position in
    // graph.nodes and silently means a different sequence if the bundle is
    // rebuilt or subsetted. Unmatched ids are counted so the load can say so.
    function nodeIndicesFromIds(ids) {
        const map = nodeIndexByIdMap();
        const indices = [];
        let missing = 0;
        (ids || []).forEach(id => {
            const index = map.get(id);
            if (index === undefined) { missing += 1; } else { indices.push(index); }
        });
        state.sessionMissingNodeIds += missing;
        return indices;
    }

    function nodeIdsFromIndices(indices) {
        return Array.from(indices).sort((left, right) => left - right).map(index => nodeId(index));
    }

    function serializeSelectionPresets() {
        const presets = {};
        state.selectionPresets.forEach((preset, slot) => {
            presets[String(slot)] = {
                node_ids: nodeIdsFromIndices(preset.nodeIndices),
                saved_at: preset.savedAt || null,
            };
        });
        return presets;
    }

    function deserializeSelectionPresets(presets) {
        state.selectionPresets = new Map();
        Object.keys(presets || {}).forEach(slotKey => {
            const slot = Number.parseInt(slotKey, 10);
            if (!Number.isInteger(slot) || slot < 0 || slot >= PRESET_SLOT_COUNT) {
                throw new Error('bad preset slot "' + slotKey + '"');
            }
            const entry = presets[slotKey] || {};
            const indices = nodeIndicesFromIds(entry.node_ids);
            if (indices.length === 0) { return; }
            state.selectionPresets.set(slot, {
                nodeIndices: new Set(indices),
                savedAt: entry.saved_at || null,
            });
        });
    }

    const VIEW_STATE_FIELDS = [
        domValueField('view', 'layout_algorithm', 'layout-algorithm'),
        domValueField('view', 'min_cluster_size', 'min-cluster-size'),
        domValueField('view', 'color_by', 'color-by'),
        domValueField('view', 'label_by', 'label-by'),
        domValueField('view', 'export_png_scale', 'export-png-scale'),
        domCheckboxField('view', 'show_node_counts', 'show-node-counts'),
        domCheckboxField('view', 'show_edge_scores', 'show-edge-scores'),
        domCheckboxField('view', 'render_cluster_bounds', 'render-cluster-bounds'),
        domCheckboxField('view', 'render_nodes', 'render-nodes'),
        domCheckboxField('view', 'reduce_elongation', 'reduce-elongation'),
        domCheckboxField('view', 'leaf_pruning_only', 'leaf-pruning-only'),
        domCheckboxField('view', 'collapse_long_paths', 'collapse-long-paths'),
        {
            section: 'view',
            key: 'sort_components_by_size',
            get: () => sortComponentsBySizeEnabled(),
            set: value => {
                document.getElementById('sort-components-by-size')
                    .setAttribute('aria-pressed', value ? 'true' : 'false');
            },
        },
        {
            // The split chart's x window, saved for the same reason the canvas
            // transform below is: reopening a session onto the whole axis would throw
            // away the part of it the session was about. Stored in threshold units, so
            // it survives any change of cap; splitChartLayout clamps it into whatever
            // range the bundle actually has.
            //
            // Restored before max_merge_events, because the window is what that cap is
            // spent on: selecting first and setting the window afterwards would leave a
            // zoomed session showing the unzoomed selection.
            section: 'view',
            key: 'split_chart_zoom',
            get: () => state.splitChartZoom
                ? {min: state.splitChartZoom.min, max: state.splitChartZoom.max}
                : null,
            set: value => {
                if (value === null || value === undefined) {
                    state.splitChartZoom = null;
                    return;
                }
                if (!Number.isFinite(value.min) || !Number.isFinite(value.max) || !(value.max > value.min)) {
                    throw new Error('malformed split chart zoom');
                }
                state.splitChartZoom = {min: Number(value.min), max: Number(value.max)};
            },
        },
        {
            // Restored after split_chart_zoom (the window it is spent on) and before
            // threshold_value, because it decides which stops exist for that field to
            // snap to. Uses rebuildMergeSeries rather than applyMaxMergeEvents for the
            // same reason: snapping here would be undone.
            section: 'view',
            key: 'max_merge_events',
            get: () => state.maxMergeEvents,
            set: value => {
                const cap = Number(value);
                if (!Number.isFinite(cap) || cap < 0) { throw new Error('not a count: ' + value); }
                rebuildMergeSeries(cap);
            },
        },
        {
            // The threshold is stored as its value, not as the slider position:
            // positions are derived from the slider's stops and shift whenever the
            // split-event cap changes.
            section: 'view',
            key: 'threshold_value',
            get: () => {
                const value = selectedThresholdValue();
                return Number.isFinite(value) ? value : null;
            },
            set: value => {
                if (value === null || value === undefined) {
                    const infinityStop = state.sliderModel.stops.find(stop => stop.threshold_value === null);
                    snapSliderToStop(infinityStop);
                    return;
                }
                const target = Number(value);
                if (!Number.isFinite(target)) { throw new Error('not a number: ' + value); }
                snapSliderToStop(nearestStopForThreshold(target));
            },
        },
        {
            section: 'view',
            key: 'view_transform',
            get: () => ({
                scale: state.viewTransform.scale,
                offsetX: state.viewTransform.offsetX,
                offsetY: state.viewTransform.offsetY,
            }),
            set: value => {
                if (!value || !Number.isFinite(value.scale) || value.scale <= 0 ||
                        !Number.isFinite(value.offsetX) || !Number.isFinite(value.offsetY)) {
                    throw new Error('malformed view transform');
                }
                // Consumed by applyComputedLayout once a layout exists, which is
                // also what would otherwise auto-fit the view and discard this.
                state.pendingRestoreViewTransform = {
                    scale: value.scale, offsetX: value.offsetX, offsetY: value.offsetY,
                };
            },
        },
        domValueField('table', 'filter', 'metadata-filter'),
        domValueField('table', 'null_order', 'metadata-null-order'),
        domValueField('table', 'rows_per_page', 'metadata-rows-per-page'),
        {
            section: 'table',
            key: 'sort',
            get: () => ({
                column_key: state.metadataSort.columnKey,
                direction: state.metadataSort.direction,
            }),
            set: value => {
                const columnKey = value && value.column_key;
                if (columnKey && !metadataColumnKeys().includes(columnKey)) {
                    throw new Error('no column "' + columnKey + '"');
                }
                state.metadataSort = {
                    columnKey: columnKey || null,
                    direction: (value && value.direction === 'desc') ? 'desc' : 'asc',
                };
            },
        },
        {
            section: 'table',
            key: 'column_widths',
            get: () => Object.fromEntries(state.metadataColumnWidths),
            set: value => {
                const widths = new Map();
                Object.keys(value || {}).forEach(columnKey => {
                    const width = Number(value[columnKey]);
                    if (!Number.isFinite(width)) { return; }
                    widths.set(columnKey, Math.max(96, Math.min(640, Math.round(width))));
                });
                state.metadataColumnWidths = widths;
            },
        },
        {
            section: 'colors',
            key: 'custom_palettes',
            get: () => state.customPalettes,
            set: value => {
                if (value === null || typeof value !== 'object') { throw new Error('not an object'); }
                // Sessions written before gradients became stop lists carry
                // {low, mid, high}; rewrite them so nothing downstream has to
                // know about that shape.
                state.customPalettes = normalizeCustomPalettes(value);
            },
        },
        {
            section: 'colors',
            key: 'categorical_columns',
            get: () => Array.from(state.categoricalColumns),
            set: value => {
                state.categoricalColumns = new Set(value || []);
                // Whether a numeric column colors as categories is *derived* from this
                // set, and rebuildMetadataCaches -- which derives it -- has already run
                // by the time a session is applied. Without recomputing it here the
                // checkbox would read categorical while the colors and the picker
                // stayed in gradient mode until the box was toggled by hand.
                state.metadataColorInfoByName.forEach((info, columnName) => {
                    if (info.baseType === 'numeric') {
                        info.type = state.categoricalColumns.has(columnName) ? 'categorical' : 'numeric';
                    }
                });
            },
        },
        {
            section: 'selection',
            key: 'node_ids',
            get: () => nodeIdsFromIndices(state.selectedNodeIndices),
            set: value => { state.selectedNodeIndices = new Set(nodeIndicesFromIds(value)); },
        },
        {
            section: 'selection',
            key: 'presets',
            get: () => serializeSelectionPresets(),
            set: value => { deserializeSelectionPresets(value); },
        },
    ];

    function collectSessionState() {
        const sections = {};
        VIEW_STATE_FIELDS.forEach(field => {
            if (!sections[field.section]) { sections[field.section] = {}; }
            sections[field.section][field.key] = field.get();
        });
        return Object.assign({
            state_version: SESSION_STATE_VERSION,
            saved_by: DOMAINATOR_VERSION,
            saved_at: new Date().toISOString(),
        }, sections);
    }

    // Returns a human-readable summary of anything that could not be restored.
    function applySessionState(appState) {
        if (!appState || typeof appState !== 'object') { return ''; }
        state.sessionMissingNodeIds = 0;
        const fieldsByPath = new Map(
            VIEW_STATE_FIELDS.map(field => [field.section + '.' + field.key, field])
        );
        const skipped = [];
        Object.keys(appState).forEach(sectionName => {
            const section = appState[sectionName];
            // Scalars at the top level are metadata about the save itself
            // (state_version, saved_by, saved_at), not restorable fields.
            if (!section || typeof section !== 'object' || Array.isArray(section)) { return; }
            Object.keys(section).forEach(key => {
                const field = fieldsByPath.get(sectionName + '.' + key);
                if (!field) {
                    skipped.push(sectionName + '.' + key + ' (unknown)');
                    return;
                }
                try {
                    field.set(section[key]);
                } catch (error) {
                    skipped.push(sectionName + '.' + key + ' (' + error.message + ')');
                }
            });
        });

        // Controls whose appearance is derived from the values just restored.
        updateComponentSortButton();
        updateCollapseLongPathsControl();
        // Re-select and re-place under the window that was just restored. The registry
        // order already does this for a session carrying view.max_merge_events, so this
        // is idempotent there; it is load-bearing for an older session that carries a
        // window but no cap, whose selection was otherwise made before the window was.
        repositionSliderStops();

        const notes = [];
        if (skipped.length > 0) {
            notes.push('skipped ' + skipped.length + ' saved setting' + (skipped.length === 1 ? '' : 's') +
                ': ' + skipped.join(', '));
        }
        if (state.sessionMissingNodeIds > 0) {
            notes.push(state.sessionMissingNodeIds.toLocaleString() +
                ' saved node id' + (state.sessionMissingNodeIds === 1 ? '' : 's') + ' are not in this bundle');
        }
        return notes.length > 0 ? ' (' + notes.join('; ') + ')' : '';
    }

    // One small dialog serves both "Rename network" and naming an extraction: the
    // network's name is a single field on the bundle, and it drives the heading, the
    // browser tab and the name of every file saved from here.
    function openNameDialog(options) {
        state.pendingNameAction = options.onConfirm;
        document.getElementById('name-dialog-title').textContent = options.title;
        document.getElementById('name-dialog-label').textContent = options.label;
        document.getElementById('name-dialog-note').textContent = options.note || '';
        document.getElementById('name-apply').textContent = options.confirmLabel || 'Apply';
        const input = document.getElementById('name-value');
        input.value = options.value || '';
        document.getElementById('name-overlay').hidden = false;
        input.focus();
        input.select();
    }

    function closeNameDialog() {
        document.getElementById('name-overlay').hidden = true;
        state.pendingNameAction = null;
    }

    function nameDialogIsOpen() {
        return !document.getElementById('name-overlay').hidden;
    }

    function confirmNameDialog() {
        const value = document.getElementById('name-value').value.trim();
        if (value === '') {
            setStatus('Enter a name.');
            return;
        }
        const action = state.pendingNameAction;
        closeNameDialog();
        if (action) {
            action(value);
        }
    }

    function renameNetwork(name) {
        if (!state.bundle) { return; }
        const previous = state.bundle.name || '';
        state.bundle.name = name;
        updateViewerTitle(name);
        setStatus('Renamed "' + previous + '" to "' + name + '". Files saved from here use the new name.');
    }

    function openRenameNetworkDialog() {
        if (!state.bundle) { return; }
        openNameDialog({
            title: 'Rename network',
            label: 'Network name',
            note: 'Shown in the heading and the browser tab, stored in the bundle, and used to name saved files.',
            value: state.bundle.name || '',
            confirmLabel: 'Rename',
            onConfirm: renameNetwork,
        });
    }

    function setupNameDialog() {
        document.getElementById('name-cancel').addEventListener('click', closeNameDialog);
        document.getElementById('name-apply').addEventListener('click', confirmNameDialog);
        document.getElementById('name-overlay').addEventListener('click', event => {
            if (event.target === event.currentTarget) { closeNameDialog(); }
        });
        document.getElementById('name-value').addEventListener('keydown', event => {
            if (event.key === 'Enter') {
                event.preventDefault();
                confirmNameDialog();
            }
        });
        document.addEventListener('keydown', event => {
            if (event.key === 'Escape' && nameDialogIsOpen()) { closeNameDialog(); }
        });
        // The heading shows the name, so it is the obvious thing to double-click.
        document.getElementById('viewer-title').addEventListener('dblclick', openRenameNetworkDialog);
    }

    function browserSupportsBundleSaving() {
        return 'CompressionStream' in window;
    }

    // Serialize a bundle and hand it to the browser as a download. Shared by
    // "Save session" and "Save extraction" so both produce the same file format.
    async function downloadBundleFile(payload, filenameBase, description) {
        try {
            setStatus(description + '...');
            const bytes = new TextEncoder().encode(JSON.stringify(payload));
            let blob;
            let extension;
            if (browserSupportsBundleSaving()) {
                const stream = new Blob([bytes]).stream().pipeThrough(new CompressionStream('gzip'));
                blob = await new Response(stream).blob();
                extension = '.dsnv';
            } else {
                // decodeBundleFile falls back to plain JSON, so this still round-trips.
                blob = new Blob([bytes], {type: 'application/json'});
                extension = '.json';
            }
            triggerDownload(URL.createObjectURL(blob), filenameBase + extension, true);
            setStatus(description + ' (' + Math.round(blob.size / 1024).toLocaleString() + ' KB).');
            return true;
        } catch (error) {
            console.error(error);
            setStatus(description + ' failed: ' + error.message);
            return false;
        }
    }

    async function saveSessionFile() {
        if (!state.bundle) { return; }
        // state.metadataByNodeIndex is bundle.metadata.rows, so table edits are
        // already in here; the file stays a valid bundle that Python readers
        // and this viewer can both open.
        await downloadBundleFile(
            Object.assign({}, state.bundle, {version: BUNDLE_VERSION, app_state: collectSessionState()}),
            exportBaseName() + '_session',
            'Saved session'
        );
    }
"""

def _selection_presets_js() -> str:
    """Ten keyboard-addressable selection presets with a hover preview.

    Plain (non-f) string like `_layout_worker_js`; see `_session_state_js`.
    """
    return r"""
    // ---------------------------------------------------------------------
    // Selection presets
    //
    // Ten slots above the hierarchy view. Press 0-9 to recall a slot (or click
    // it), Shift+0-9 to store the current selection into one, and hover a slot
    // to outline its nodes in gray without disturbing the live selection.
    //
    // Storing is bound to Shift rather than Ctrl/Cmd because browsers reserve
    // Ctrl/Cmd+1-9 for tab switching and Ctrl/Cmd+0 for zoom reset; a page
    // cannot intercept those. Slots are keyed off event.code (Digit0..Digit9)
    // so neither the keyboard layout nor Shift's symbol mapping matters.
    // ---------------------------------------------------------------------

    function presetPreviewNodeSet() {
        if (state.presetPreviewSlot === null) { return null; }
        const preset = state.selectionPresets.get(state.presetPreviewSlot);
        return preset && preset.nodeIndices.size > 0 ? preset.nodeIndices : null;
    }

    function renderPresetSlots() {
        const container = document.getElementById('preset-slots');
        if (!container) { return; }
        container.querySelectorAll('[data-preset-slot]').forEach(button => {
            const slot = Number(button.dataset.presetSlot);
            const preset = state.selectionPresets.get(slot);
            const count = preset ? preset.nodeIndices.size : 0;
            button.classList.toggle('preset-slot-filled', count > 0);
            button.disabled = !state.bundle;
            button.title = count > 0
                ? 'Preset ' + slot + ': ' + count.toLocaleString() + ' nodes.' +
                  ' Click (or press ' + slot + ') to replace the selection with it;' +
                  ' Shift-click to add it; Alt+Shift-click to subtract it;' +
                  ' hover to preview; Shift+' + slot + ' to overwrite it.'
                : 'Preset ' + slot + ' is empty. Select nodes, then press Shift+' + slot + ' to store them here.';
        });
    }

    function storeSelectionPreset(slot) {
        if (!state.bundle) { return; }
        if (state.selectedNodeIndices.size === 0) {
            // Storing nothing is how you clear a slot.
            state.selectionPresets.delete(slot);
            setStatus('Cleared preset ' + slot + ' (nothing was selected).');
        } else {
            state.selectionPresets.set(slot, {
                nodeIndices: new Set(state.selectedNodeIndices),
                savedAt: new Date().toISOString(),
            });
            setStatus('Stored ' + state.selectedNodeIndices.size.toLocaleString() +
                ' nodes in preset ' + slot + '.');
        }
        renderPresetSlots();
        scheduleClusterRender();
    }

    // mode: 'replace' (the default), 'add' to union the preset into the current
    // selection, or 'subtract' to remove it. Shift/Alt+Shift on a slot button pick
    // the latter two, matching the canvas box-select gesture where Shift selects
    // and adding Alt deselects.
    function recallSelectionPreset(slot, mode = 'replace') {
        if (!state.bundle) { return; }
        const preset = state.selectionPresets.get(slot);
        if (!preset) {
            setStatus('Preset ' + slot + ' is empty.');
            return;
        }
        const before = state.selectedNodeIndices.size;
        if (mode === 'add') {
            preset.nodeIndices.forEach(nodeIndex => state.selectedNodeIndices.add(nodeIndex));
        } else if (mode === 'subtract') {
            preset.nodeIndices.forEach(nodeIndex => state.selectedNodeIndices.delete(nodeIndex));
        } else {
            state.selectedNodeIndices = new Set(preset.nodeIndices);
        }
        const after = state.selectedNodeIndices.size;
        resetMetadataPage();
        clearMetadataRowSelection();
        if (mode === 'add') {
            setStatus('Added preset ' + slot + ' to the selection (+' +
                (after - before).toLocaleString() + ' nodes, ' + after.toLocaleString() + ' total).');
        } else if (mode === 'subtract') {
            setStatus('Subtracted preset ' + slot + ' from the selection (-' +
                (before - after).toLocaleString() + ' nodes, ' + after.toLocaleString() + ' remaining).');
        } else {
            setStatus('Recalled preset ' + slot + ' (' + after.toLocaleString() + ' nodes).');
        }
        renderClusterView();
        updateMetadataTable();
    }

    function presetClickMode(event) {
        if (!event.shiftKey) { return 'replace'; }
        return event.altKey ? 'subtract' : 'add';
    }

    function setPresetPreviewSlot(slot) {
        if (state.presetPreviewSlot === slot) { return; }
        state.presetPreviewSlot = slot;
        scheduleClusterRender();
    }

    // True when a keystroke belongs to a control the user is typing into, in
    // which case the bare-digit shortcuts must not fire.
    function keyboardTargetIsTextEntry(target) {
        if (!target || !target.tagName) { return false; }
        if (target.isContentEditable) { return true; }
        return ['INPUT', 'SELECT', 'TEXTAREA', 'OPTION'].includes(target.tagName);
    }

    function handlePresetKeydown(event) {
        if (!state.bundle) { return; }
        if (event.ctrlKey || event.metaKey || event.altKey) { return; }
        if (keyboardTargetIsTextEntry(event.target)) { return; }
        if (colorPickerIsOpen() || metadataPasteDialogIsOpen() || nameDialogIsOpen()
            || columnChartIsOpen()) { return; }
        const match = /^Digit([0-9])$/.exec(event.code || '');
        if (!match) { return; }
        const slot = Number(match[1]);
        event.preventDefault();
        if (event.shiftKey) {
            storeSelectionPreset(slot);
        } else {
            recallSelectionPreset(slot);
        }
    }

    function setupSelectionPresets() {
        const container = document.getElementById('preset-slots');
        container.addEventListener('click', event => {
            const button = event.target.closest('[data-preset-slot]');
            if (button) {
                recallSelectionPreset(Number(button.dataset.presetSlot), presetClickMode(event));
            }
        });
        // mouseover/mouseout (not mouseenter/mouseleave) so one delegated pair
        // covers all ten buttons.
        container.addEventListener('mouseover', event => {
            const button = event.target.closest('[data-preset-slot]');
            if (button) { setPresetPreviewSlot(Number(button.dataset.presetSlot)); }
        });
        container.addEventListener('mouseout', event => {
            const button = event.target.closest('[data-preset-slot]');
            if (button) { setPresetPreviewSlot(null); }
        });
        // Keyboard focus should preview too, so the slots are usable without a mouse.
        container.addEventListener('focusin', event => {
            const button = event.target.closest('[data-preset-slot]');
            if (button) { setPresetPreviewSlot(Number(button.dataset.presetSlot)); }
        });
        container.addEventListener('focusout', () => setPresetPreviewSlot(null));
        document.addEventListener('keydown', handlePresetKeydown);
        renderPresetSlots();
    }
"""

def _table_editing_js() -> str:
    """Editable metadata cells, user-added columns, and bulk fills.

    Plain (non-f) string like `_layout_worker_js`; see `_session_state_js`.
    """
    return r"""
    // ---------------------------------------------------------------------
    // Metadata editing
    //
    // Edits write straight into state.metadataByNodeIndex, which *is*
    // bundle.metadata.rows (installBundle assigns it by reference). A saved
    // session is therefore an ordinary bundle whose metadata simply reflects
    // the annotations -- no overlay to reconcile, and Python readers
    // (ssn_bundle, ssn_navigator) see the edits with no extra work.
    //
    // node_id is not editable: it is the graph key that node indices, the
    // hierarchy and every saved selection resolve through.
    // ---------------------------------------------------------------------

    function columnIsEditable(columnKey) {
        return Boolean(state.bundle) && columnKey !== 'node_id' && Boolean(metadataColumn(columnKey));
    }

    function columnIsViewerAdded(columnKey) {
        const column = metadataColumn(columnKey);
        return Boolean(column && column.origin === 'viewer');
    }

    // Parse user text into the column's declared type. Returns {ok, value, widenTo}
    // so callers can reject rather than silently coercing "abc" in a numeric
    // column to NaN (or quietly turning the whole column into text).
    //
    // Typing 12.5 into an int column widens the column to float rather than
    // truncating to 12: the type was only inferred from the values that happened
    // to be in the source TSV, so widening is what Python would have inferred had
    // the data looked like this to begin with, and it loses nothing.
    function parseMetadataInput(columnKey, rawText) {
        const text = String(rawText === null || rawText === undefined ? '' : rawText).trim();
        if (text === '' || text === '—') {
            return {ok: true, value: null};
        }
        const columnType = metadataColumnType(columnKey);
        if (columnType === 'int' || columnType === 'float') {
            // Reject "12abc": Number() is lenient about surrounding whitespace but
            // not about trailing garbage, which is the failure mode worth catching.
            const numeric = Number(text);
            if (!Number.isFinite(numeric)) {
                return {ok: false, reason: '"' + text + '" is not a number'};
            }
            if (columnType === 'int' && !Number.isInteger(numeric)) {
                return {ok: true, value: numeric, widenTo: 'float'};
            }
            return {ok: true, value: numeric};
        }
        return {ok: true, value: text};
    }

    // Applies a parse result's requested type widening, if any. Returns whether
    // the column list changed, so callers can refresh the column-driven menus.
    function applyColumnWidening(columnName, parsed) {
        if (!parsed.widenTo) { return false; }
        const column = metadataColumn(columnName);
        if (!column || column.type === parsed.widenTo) { return false; }
        column.type = parsed.widenTo;
        setStatus('Column "' + columnName + '" now holds fractional values; its type widened to float.');
        return true;
    }

    // Writes the value without touching caches; callers batch a single
    // refreshAfterMetadataEdit() so a bulk fill rebuilds caches once, not once
    // per row (rebuildMetadataCaches is O(rows x columns)).
    function writeMetadataValue(nodeIndex, columnName, value) {
        const columnIndex = state.metadataColumnIndexByName.get(columnName);
        if (columnIndex === undefined) { return false; }
        const row = state.metadataByNodeIndex[nodeIndex];
        if (!row) { return false; }
        row[columnIndex] = value;
        return true;
    }

    // Everything downstream of a metadata value: the search index and numeric
    // min/max live in rebuildMetadataCaches, node colors are cached separately.
    function refreshAfterMetadataEdit(options) {
        const settings = options || {};
        rebuildMetadataCaches();
        if (settings.columnsChanged) {
            populateMetadataControlsPreservingSelection();
            populateEditableColumnMenus();
        }
        rebuildNodeColorCache();
        updateMetadataTable();
        renderClusterView();
        if (colorPickerIsOpen()) {
            openColorPicker();
        }
    }

    // populateMetadataControls resets color-by/label-by to the bundle defaults,
    // which would undo the user's choice every time a column is added.
    function populateMetadataControlsPreservingSelection() {
        const colorBy = document.getElementById('color-by');
        const labelBy = document.getElementById('label-by');
        const previousColor = colorBy.value;
        const previousLabel = labelBy.value;
        populateMetadataControls();
        if (Array.from(colorBy.options).some(option => option.value === previousColor)) {
            colorBy.value = previousColor;
        }
        if (Array.from(labelBy.options).some(option => option.value === previousLabel)) {
            labelBy.value = previousLabel;
        }
    }

    function setMetadataValue(nodeIndex, columnName, rawText) {
        if (!columnIsEditable(columnName)) { return false; }
        const parsed = parseMetadataInput(columnName, rawText);
        if (!parsed.ok) {
            setStatus('Cannot set "' + columnName + '": ' + parsed.reason + '.');
            return false;
        }
        if (!writeMetadataValue(nodeIndex, columnName, parsed.value)) { return false; }
        const widened = applyColumnWidening(columnName, parsed);
        refreshAfterMetadataEdit({columnsChanged: widened});
        return true;
    }

    // ---- in-place cell editing ----

    function beginMetadataCellEdit(cell) {
        if (state.metadataEditCell) { return; }
        const row = cell.closest('tr[data-node-index]');
        const columnKey = cell.dataset.columnKey;
        if (!row || !columnIsEditable(columnKey)) { return; }
        const nodeIndex = Number(row.dataset.nodeIndex);
        const value = metadataValue(nodeIndex, columnKey);

        const input = document.createElement('input');
        input.type = 'text';
        input.className = 'metadata-cell-input';
        input.value = isMissingMetadataValue(value) ? '' : String(value);
        input.setAttribute('aria-label', 'Edit ' + columnKey);
        cell.textContent = '';
        cell.appendChild(input);
        state.metadataEditCell = {cell, nodeIndex, columnKey};
        input.focus();
        input.select();

        const commit = () => finishMetadataCellEdit(true);
        input.addEventListener('blur', commit);
        input.addEventListener('keydown', event => {
            if (event.key === 'Enter') {
                event.preventDefault();
                input.removeEventListener('blur', commit);
                finishMetadataCellEdit(true);
            } else if (event.key === 'Escape') {
                event.preventDefault();
                input.removeEventListener('blur', commit);
                finishMetadataCellEdit(false);
            }
        });
    }

    function finishMetadataCellEdit(commit) {
        const editing = state.metadataEditCell;
        if (!editing) { return; }
        state.metadataEditCell = null;
        const input = editing.cell.querySelector('input');
        const rawText = input ? input.value : '';
        // Clicking away from an untouched cell is the common case, and it should
        // not pay for a cache rebuild -- it only needs the <input> swapped out.
        const parsed = commit ? parseMetadataInput(editing.columnKey, rawText) : null;
        const unchanged = parsed && parsed.ok && !parsed.widenTo &&
            parsed.value === metadataValue(editing.nodeIndex, editing.columnKey);
        if (!commit || unchanged ||
                !setMetadataValue(editing.nodeIndex, editing.columnKey, rawText)) {
            // Canceled, unchanged, or rejected: repaint to restore the cell text.
            updateMetadataTable();
        }
    }

    // ---- adding and removing columns ----

    // Shared by add and rename. `node_id` is reserved: it is the table's synthetic
    // first column, so a metadata column of that name would be unreachable.
    function validateColumnName(name) {
        const columnName = String(name || '').trim();
        if (columnName === '') {
            return {ok: false, reason: 'Enter a column name.'};
        }
        if (columnName === 'node_id' || state.metadataColumnIndexByName.has(columnName)) {
            return {ok: false, reason: 'A column named "' + columnName + '" already exists.'};
        }
        return {ok: true, name: columnName};
    }

    function addMetadataColumn(name, columnType) {
        if (!state.bundle) { return false; }
        const validated = validateColumnName(name);
        if (!validated.ok) {
            setStatus(validated.reason);
            return false;
        }
        const columnName = validated.name;
        // `origin` marks columns created here so they can be deleted again. It is
        // additive: build_ssn_viewer writes only {name, type} and no reader
        // inspects extra keys, so a bundle carrying it stays valid.
        state.metadataColumns.push({
            name: columnName,
            type: columnType === 'number' ? 'float' : 'str',
            origin: 'viewer',
        });
        state.metadataByNodeIndex.forEach(row => { row.push(null); });
        refreshAfterMetadataEdit({columnsChanged: true});
        setStatus('Added column "' + columnName + '".');
        return true;
    }

    // Menus whose option values are column names. Which column each one points at
    // is held in the DOM rather than in `state`, so it is one more thing keyed by
    // column name that a rename has to carry across.
    const COLUMN_NAME_MENU_IDS = [
        'color-by',
        'label-by',
        'metadata-fill-column',
        'metadata-paste-column',
        'metadata-rename-column',
        'metadata-delete-column',
        'metadata-select-column',
    ];

    // Retarget every menu currently pointing at `oldName`. Rewriting the option
    // updates its select's value, so the rebuild that follows -- which restores
    // whatever each menu was showing -- carries the choice over on its own.
    function renameColumnInMenus(oldName, newName) {
        COLUMN_NAME_MENU_IDS.forEach(menuId => {
            Array.from(document.getElementById(menuId).options).forEach(option => {
                if (option.value === oldName) {
                    option.value = newName;
                }
            });
        });
    }

    // The viewer's counterpart to `build_ssn.py --lb <threshold> --cluster`: writes
    // the connected-component number at the current threshold into a column. The
    // default name matches that tool's CLUSTER_COLUMN so the two are interchangeable
    // downstream, and an existing column of the same name is overwritten, as there.
    function addClusterColumn(name) {
        if (!state.bundle) { return false; }
        const columnName = String(name || '').trim() || 'SSN_cluster';
        if (columnName === 'node_id') {
            setStatus('"node_id" is reserved; choose another column name.');
            return false;
        }
        const existing = metadataColumn(columnName);
        if (existing && existing.origin !== 'viewer' &&
                !window.confirm('Overwrite the column "' + columnName + '" with cluster numbers?\n\n' +
                    'It came from this bundle\'s metadata, not from the viewer, and its current values will be lost.')) {
            return false;
        }

        const clusterNumbers = clusterNumbersAtCurrentThreshold();
        if (!existing) {
            state.metadataColumns.push({name: columnName, type: 'int', origin: 'viewer'});
            state.metadataByNodeIndex.forEach(row => { row.push(null); });
        } else {
            existing.type = 'int';
        }
        rebuildMetadataCaches();  // so writeMetadataValue can find the new column index
        state.allNodeIndices.forEach(nodeIndex => {
            writeMetadataValue(nodeIndex, columnName, clusterNumbers.numberFor(nodeIndex));
        });
        // Cluster ids are labels, not magnitudes, so color them as discrete
        // categories -- the same thing build_ssn_viewer.py's --categorical asks for.
        state.categoricalColumns.add(columnName);

        refreshAfterMetadataEdit({columnsChanged: !existing});
        const stop = currentSliderStop();
        setStatus('Wrote ' + clusterNumbers.clusterCount.toLocaleString() + ' cluster numbers into "' +
            columnName + '" at threshold ' + (stop ? stop.threshold_label : '∞') + '.');
        return true;
    }

    // Moves any state that is keyed by column name along with the column. Values
    // themselves are positional, so the rows need no rewriting.
    function renameMetadataColumn(oldName, newName) {
        if (!state.bundle || !columnIsEditable(oldName)) { return false; }
        const trimmed = String(newName || '').trim();
        if (trimmed === oldName) { return true; }
        const validated = validateColumnName(trimmed);
        if (!validated.ok) {
            setStatus(validated.reason);
            return false;
        }
        const columnName = validated.name;

        // Both palette spellings, because paletteKey() suffixes a numeric column
        // that is currently colored categorically.
        ['', '\u0000categorical'].forEach(suffix => {
            const fromKey = oldName + suffix;
            if (Object.prototype.hasOwnProperty.call(state.customPalettes, fromKey)) {
                state.customPalettes[columnName + suffix] = state.customPalettes[fromKey];
                delete state.customPalettes[fromKey];
            }
        });
        if (state.categoricalColumns.delete(oldName)) {
            state.categoricalColumns.add(columnName);
        }
        if (state.metadataColumnWidths.has(oldName)) {
            state.metadataColumnWidths.set(columnName, state.metadataColumnWidths.get(oldName));
            state.metadataColumnWidths.delete(oldName);
        }
        if (state.metadataSort.columnKey === oldName) {
            state.metadataSort.columnKey = columnName;
        }
        renameColumnInMenus(oldName, columnName);
        metadataColumn(oldName).name = columnName;

        refreshAfterMetadataEdit({columnsChanged: true});
        setStatus('Renamed column "' + oldName + '" to "' + columnName + '".');
        return true;
    }

    function deleteMetadataColumn(columnName) {
        if (!columnIsEditable(columnName)) { return false; }
        // Dropping a column that came from the source data is the one destructive,
        // un-undoable edit in the viewer -- and it only really bites once the
        // session is saved over the original. Columns created here are cheap to
        // recreate, so they go without a prompt.
        if (!columnIsViewerAdded(columnName) &&
                !window.confirm('Delete the column "' + columnName + '"?\n\n' +
                    'It came from this bundle\'s metadata, not from the viewer. ' +
                    'Reload the bundle to get it back; a session saved after this will not contain it.')) {
            return false;
        }
        const columnIndex = state.metadataColumnIndexByName.get(columnName);
        state.metadataColumns.splice(columnIndex, 1);
        state.metadataByNodeIndex.forEach(row => { row.splice(columnIndex, 1); });
        state.metadataColumnWidths.delete(columnName);
        // Palettes are stored under paletteKey(), which suffixes a numeric column
        // that is being colored categorically; clear both spellings.
        delete state.customPalettes[columnName];
        delete state.customPalettes[columnName + '\u0000categorical'];
        state.categoricalColumns.delete(columnName);
        if (state.metadataSort.columnKey === columnName) {
            state.metadataSort = {columnKey: null, direction: 'asc'};
        }
        refreshAfterMetadataEdit({columnsChanged: true});
        setStatus('Deleted column "' + columnName + '".');
        return true;
    }

    // ---- bulk fill ----

    // The three fill targets. "selection" is the annotation workflow (select a
    // cluster on the canvas, label it); "rows" uses the table's staged rows;
    // "page" matches what the per-column Copy button copies.
    function bulkFillTargetNodeIndices(target) {
        if (target === 'rows') {
            return Array.from(state.selectedMetadataNodeIndices).sort((left, right) => left - right);
        }
        if (target === 'page') {
            return state.renderedNodeIndices.slice();
        }
        if (target === 'all') {
            return state.allNodeIndices.slice();
        }
        return Array.from(state.selectedNodeIndices).sort((left, right) => left - right);
    }

    function bulkFillTargetLabel(target) {
        if (target === 'rows') { return 'staged table rows'; }
        if (target === 'page') { return 'rows on this page'; }
        if (target === 'all') { return 'nodes in the network'; }
        return 'selected nodes';
    }

    function applyBulkFill() {
        if (!state.bundle) { return; }
        const columnName = document.getElementById('metadata-fill-column').value;
        const target = document.getElementById('metadata-fill-target').value;
        const rawText = document.getElementById('metadata-fill-value').value;
        if (!columnIsEditable(columnName)) {
            setStatus('Choose a column to fill.');
            return;
        }
        const nodeIndices = bulkFillTargetNodeIndices(target);
        if (nodeIndices.length === 0) {
            setStatus('No ' + bulkFillTargetLabel(target) + ' to fill.');
            return;
        }
        const parsed = parseMetadataInput(columnName, rawText);
        if (!parsed.ok) {
            setStatus('Cannot fill "' + columnName + '": ' + parsed.reason + '.');
            return;
        }
        nodeIndices.forEach(nodeIndex => writeMetadataValue(nodeIndex, columnName, parsed.value));
        const widened = applyColumnWidening(columnName, parsed);
        refreshAfterMetadataEdit({columnsChanged: widened});
        setStatus('Set "' + columnName + '" on ' + nodeIndices.length.toLocaleString() + ' ' +
            bulkFillTargetLabel(target) + '.');
    }

    // Paste applies to the rendered page in display order -- the exact inverse
    // of the per-column Copy button. Counts must match: a short or long paste
    // almost always means the page, filter or sort moved since the copy, and
    // silently filling a prefix would mislabel rows.
    function applyPastedColumn() {
        if (!state.bundle) { return; }
        const columnName = document.getElementById('metadata-paste-column').value;
        const text = document.getElementById('metadata-paste-values').value;
        if (!columnIsEditable(columnName)) {
            setStatus('Choose a column to paste into.');
            return;
        }
        const lines = text.split(/\r?\n/);
        while (lines.length > 0 && lines[lines.length - 1].trim() === '') {
            lines.pop();
        }
        const nodeIndices = state.renderedNodeIndices;
        if (lines.length !== nodeIndices.length) {
            setStatus('Paste has ' + lines.length.toLocaleString() + ' values but this page shows ' +
                nodeIndices.length.toLocaleString() + ' rows; they must match.');
            return;
        }
        const parsedValues = [];
        let widenTo = null;
        for (let i = 0; i < lines.length; i++) {
            // Tolerate a pasted multi-column block by taking the first field.
            const parsed = parseMetadataInput(columnName, lines[i].split('\t')[0]);
            if (!parsed.ok) {
                setStatus('Cannot paste into "' + columnName + '" at line ' + (i + 1) + ': ' + parsed.reason + '.');
                return;
            }
            widenTo = widenTo || parsed.widenTo || null;
            parsedValues.push(parsed.value);
        }
        nodeIndices.forEach((nodeIndex, position) => {
            writeMetadataValue(nodeIndex, columnName, parsedValues[position]);
        });
        const widened = applyColumnWidening(columnName, {widenTo});
        refreshAfterMetadataEdit({columnsChanged: widened});
        closeMetadataPasteDialog();
        setStatus('Pasted ' + parsedValues.length.toLocaleString() + ' values into "' + columnName + '".');
    }

    // ---- editing controls ----

    function populateEditableColumnMenus() {
        const editable = state.metadataColumns.map(column => column.name);
        [
            document.getElementById('metadata-fill-column'),
            document.getElementById('metadata-paste-column'),
            document.getElementById('metadata-rename-column'),
        ].forEach(select => {
            const previous = select.value;
            select.innerHTML = '';
            editable.forEach(name => {
                const option = document.createElement('option');
                option.value = name;
                option.textContent = name;
                select.appendChild(option);
            });
            if (editable.includes(previous)) { select.value = previous; }
            select.disabled = editable.length === 0;
        });

        const deleteSelect = document.getElementById('metadata-delete-column');
        const previousDelete = deleteSelect.value;
        deleteSelect.innerHTML = '';
        state.metadataColumns.forEach(column => {
            const option = document.createElement('option');
            option.value = column.name;
            // Say where each column came from, since deleting a source column
            // discards data the viewer cannot regenerate.
            option.textContent = column.name +
                (column.origin === 'viewer' ? ' (added here)' : ' (from bundle)');
            deleteSelect.appendChild(option);
        });
        if (editable.includes(previousDelete)) { deleteSelect.value = previousDelete; }
        deleteSelect.disabled = editable.length === 0;
        document.getElementById('metadata-delete-column-apply').disabled = editable.length === 0;
        document.getElementById('metadata-fill-apply').disabled = editable.length === 0;
        document.getElementById('metadata-paste-open').disabled = editable.length === 0;
        document.getElementById('metadata-rename-apply').disabled = editable.length === 0;
        updateDeleteColumnNote();
        populateMetadataSelectMenus();
        updateMetadataSelectOperands();
        updateMetadataSelectNote();
    }

    function updateDeleteColumnNote() {
        const note = document.getElementById('metadata-delete-note');
        const columnName = document.getElementById('metadata-delete-column').value;
        if (!columnIsEditable(columnName)) {
            note.textContent = 'There are no columns to delete.';
        } else if (columnIsViewerAdded(columnName)) {
            note.textContent = 'Created in the viewer.';
        } else {
            note.textContent = 'From the bundle — reload it to restore this column.';
        }
    }


    // ---------------------------------------------------------------------
    // Select by value
    //
    // A column/comparison/value query that edits the node selection. Almost
    // everything downstream -- the per-column charts, "Set column", Export table
    // TSV, Save extraction -- is driven by state.selectedNodeIndices, and until
    // now the only way to build that set was to find the nodes on the canvas by
    // eye or to sort the table and tick rows.
    //
    // Two conventions, both surfaced in the panel's own note:
    //   * Text comparisons run against the *stored* value, not the table's
    //     formatted rendering of it, so a thousands separator never has to be
    //     typed to match the cell that displays as "1,234".
    //   * Ordering comparisons go through compareMetadataValues(), the very
    //     comparator the table's column sort uses -- so "greater than" means what
    //     sorting that column already showed, on text columns as well as numeric
    //     ones. Cells with no value never match anything.
    // ---------------------------------------------------------------------

    // node_id is offered alongside the metadata columns: it is the one field every
    // node has, and an accession prefix is a common way into a network.
    const METADATA_SELECT_NODE_ID = '__node_id__';

    // `operands` drives which value fields the panel shows, so how many bounds a
    // comparison takes is stated once instead of re-derived at each use.
    const METADATA_SELECT_OPS = [
        {id: 'contains', label: 'contains', operands: 1},
        {id: 'exact', label: 'is exactly', operands: 1},
        {id: 'gt', label: 'greater than', operands: 1},
        {id: 'lt', label: 'less than', operands: 1},
        {id: 'between', label: 'between', operands: 2},
        {id: 'regex', label: 'matches regex', operands: 1},
    ];

    function metadataSelectOp(opId) {
        return METADATA_SELECT_OPS.find(op => op.id === opId) || METADATA_SELECT_OPS[0];
    }

    function metadataSelectFieldLabel(field) {
        return field === METADATA_SELECT_NODE_ID ? 'node_id' : field;
    }

    function metadataSelectFieldValue(nodeIndex, field) {
        return field === METADATA_SELECT_NODE_ID
            ? nodeId(nodeIndex)
            : metadataValue(nodeIndex, field);
    }

    // Returns {test} on success, or {error} carrying a message for the panel note.
    // `test` takes the stored value and its String() form, since the substring and
    // regex comparisons want text while the ordering ones want the value itself.
    function metadataSelectMatcher(opId, firstText, secondText) {
        const first = String(firstText ?? '').trim();
        const second = String(secondText ?? '').trim();
        if (first === '') {
            return {error: 'Type a value to compare against.'};
        }
        if (opId === 'regex') {
            let pattern;
            try {
                // Case-insensitive, and unanchored, so it reads like the search box
                // rather than like a fullmatch.
                pattern = new RegExp(first, 'i');
            } catch (error) {
                return {error: 'Not a valid regular expression: ' + error.message};
            }
            return {test: (raw, text) => pattern.test(text)};
        }
        if (opId === 'contains') {
            const needle = first.toLowerCase();
            return {test: (raw, text) => text.toLowerCase().includes(needle)};
        }
        if (opId === 'exact') {
            const wanted = first.toLowerCase();
            // Numeric equality as well as text, so "5.0" and "5" both match a float
            // column's 5 and the match does not depend on how the number was typed.
            const wantedNumber = Number(first);
            const numeric = Number.isFinite(wantedNumber);
            return {test: (raw, text) => text.toLowerCase() === wanted
                || (numeric && typeof raw === 'number' && raw === wantedNumber)};
        }
        if (opId === 'gt') {
            return {test: raw => compareMetadataValues(raw, first) > 0};
        }
        if (opId === 'lt') {
            return {test: raw => compareMetadataValues(raw, first) < 0};
        }
        if (opId === 'between') {
            if (second === '') {
                return {error: 'Type both bounds, or pick another comparison.'};
            }
            // Inclusive at both ends, and tolerant of bounds typed the wrong way
            // round -- ordering them by the same comparator that will test them.
            const reversed = compareMetadataValues(first, second) > 0;
            const low = reversed ? second : first;
            const high = reversed ? first : second;
            return {test: raw => compareMetadataValues(raw, low) >= 0
                && compareMetadataValues(raw, high) <= 0};
        }
        return {error: 'Unknown comparison.'};
    }

    function metadataSelectMatches(nodeIndices, field, matcher) {
        return nodeIndices.filter(nodeIndex => {
            const raw = metadataSelectFieldValue(nodeIndex, field);
            // A blank cell is an absence, not a value: it matches no comparison,
            // including "contains ''" (which the matcher rejects outright anyway).
            if (isMissingMetadataValue(raw)) {
                return false;
            }
            return matcher.test(raw, String(raw));
        });
    }

    function setMetadataSelectNote(text) {
        document.getElementById('metadata-select-note').textContent = text;
    }

    // Says how the chosen column will be compared, which is the one thing about
    // this panel that is not visible from its controls.
    function updateMetadataSelectNote() {
        const field = document.getElementById('metadata-select-column').value;
        if (!field) {
            setMetadataSelectNote('');
            return;
        }
        const columnType = field === METADATA_SELECT_NODE_ID ? 'string' : metadataColumnType(field);
        const kind = (columnType === 'int' || columnType === 'float') ? 'numbers' : 'text';
        setMetadataSelectNote('Compares ' + metadataSelectFieldLabel(field) + ' as ' + kind +
            '; blank cells never match.');
    }

    function updateMetadataSelectOperands() {
        const op = metadataSelectOp(document.getElementById('metadata-select-op').value);
        const twoOperands = op.operands === 2;
        document.getElementById('metadata-select-and').hidden = !twoOperands;
        document.getElementById('metadata-select-value2').hidden = !twoOperands;
        document.getElementById('metadata-select-value').placeholder = twoOperands ? 'lower bound' : 'value';
    }

    function populateMetadataSelectMenus() {
        const opSelect = document.getElementById('metadata-select-op');
        if (opSelect.options.length === 0) {
            METADATA_SELECT_OPS.forEach(op => {
                const option = document.createElement('option');
                option.value = op.id;
                option.textContent = op.label;
                opSelect.appendChild(option);
            });
        }
        const select = document.getElementById('metadata-select-column');
        const previous = select.value;
        select.innerHTML = '';
        const nodeIdOption = document.createElement('option');
        nodeIdOption.value = METADATA_SELECT_NODE_ID;
        nodeIdOption.textContent = 'node_id';
        select.appendChild(nodeIdOption);
        state.metadataColumns.forEach(column => {
            const option = document.createElement('option');
            option.value = column.name;
            option.textContent = column.name;
            select.appendChild(option);
        });
        const available = Array.from(select.options).map(option => option.value);
        select.value = available.includes(previous) ? previous : METADATA_SELECT_NODE_ID;
    }

    // `mode` is 'add' (search the whole network and union the matches in),
    // 'remove' (search the current selection and drop the matches) or 'subset'
    // (search the current selection and keep only the matches). Only 'add' can
    // grow the selection; that asymmetry is what makes the three buttons compose
    // into an intersection of queries.
    function applyMetadataSelectByValue(mode) {
        if (!state.bundle) {
            return false;
        }
        const field = document.getElementById('metadata-select-column').value;
        const opId = document.getElementById('metadata-select-op').value;
        const matcher = metadataSelectMatcher(
            opId,
            document.getElementById('metadata-select-value').value,
            document.getElementById('metadata-select-value2').value,
        );
        if (matcher.error) {
            setMetadataSelectNote(matcher.error);
            return false;
        }
        const scope = mode === 'add'
            ? state.allNodeIndices
            : Array.from(state.selectedNodeIndices).sort((left, right) => left - right);
        const matches = metadataSelectMatches(scope, field, matcher);
        const before = state.selectedNodeIndices.size;
        if (mode === 'add') {
            matches.forEach(nodeIndex => state.selectedNodeIndices.add(nodeIndex));
        } else if (mode === 'remove') {
            matches.forEach(nodeIndex => state.selectedNodeIndices.delete(nodeIndex));
        } else {
            state.selectedNodeIndices = new Set(matches);
        }
        const verb = {add: 'Matched', remove: 'Removed', subset: 'Kept'}[mode];
        const scopeLabel = mode === 'add' ? 'the network' : 'the selection';
        setMetadataSelectNote(verb + ' ' + matches.length.toLocaleString() + ' of ' +
            scope.length.toLocaleString() + ' node' + (scope.length === 1 ? '' : 's') +
            ' in ' + scopeLabel + '.');
        setStatus(metadataSelectFieldLabel(field) + ' ' + metadataSelectOp(opId).label +
            ': ' + verb.toLowerCase() + ' ' + matches.length.toLocaleString() + ' node' +
            (matches.length === 1 ? '' : 's') + '; ' +
            state.selectedNodeIndices.size.toLocaleString() + ' now selected (was ' +
            before.toLocaleString() + ').');
        resetMetadataPage();
        renderClusterView();
        updateMetadataTable();
        return true;
    }

    // One panel at a time: showing every editing control at once made the panel
    // busy, and these are occasional actions rather than per-row ones.
    const METADATA_EDIT_PANELS = ['add', 'set', 'rename', 'delete', 'select'];

    function toggleMetadataEditPanel(name) {
        const wasOpen = document.getElementById('metadata-panel-' + name)
            .getAttribute('aria-expanded') === 'true';
        METADATA_EDIT_PANELS.forEach(panelName => {
            const open = !wasOpen && panelName === name;
            document.getElementById('metadata-panel-' + panelName)
                .setAttribute('aria-expanded', open ? 'true' : 'false');
            document.getElementById('metadata-' + panelName + '-panel').hidden = !open;
        });
        if (wasOpen) { return; }
        populateEditableColumnMenus();
        const focusTarget = {
            add: 'metadata-new-column-name',
            set: 'metadata-fill-value',
            rename: 'metadata-rename-value',
            delete: 'metadata-delete-column',
            select: 'metadata-select-value',
        }[name];
        document.getElementById(focusTarget).focus();
    }

    function closeMetadataEditPanels() {
        METADATA_EDIT_PANELS.forEach(panelName => {
            document.getElementById('metadata-panel-' + panelName).setAttribute('aria-expanded', 'false');
            document.getElementById('metadata-' + panelName + '-panel').hidden = true;
        });
    }

    function applyColumnRename() {
        const input = document.getElementById('metadata-rename-value');
        if (renameMetadataColumn(document.getElementById('metadata-rename-column').value, input.value)) {
            input.value = '';
        }
    }

    function openMetadataPasteDialog() {
        populateEditableColumnMenus();
        document.getElementById('metadata-paste-count').textContent =
            state.renderedNodeIndices.length.toLocaleString();
        document.getElementById('metadata-paste-values').value = '';
        document.getElementById('metadata-paste-overlay').hidden = false;
        document.getElementById('metadata-paste-values').focus();
    }

    function closeMetadataPasteDialog() {
        document.getElementById('metadata-paste-overlay').hidden = true;
    }

    function metadataPasteDialogIsOpen() {
        return !document.getElementById('metadata-paste-overlay').hidden;
    }

    function setupMetadataEditing() {
        METADATA_EDIT_PANELS.forEach(panelName => {
            document.getElementById('metadata-panel-' + panelName)
                .addEventListener('click', () => toggleMetadataEditPanel(panelName));
        });
        document.getElementById('metadata-delete-column')
            .addEventListener('change', updateDeleteColumnNote);
        document.getElementById('metadata-select-column')
            .addEventListener('change', updateMetadataSelectNote);
        document.getElementById('metadata-select-op').addEventListener('change', () => {
            updateMetadataSelectOperands();
            updateMetadataSelectNote();
        });
        ['add', 'remove', 'subset'].forEach(mode => {
            document.getElementById('metadata-select-' + mode)
                .addEventListener('click', () => applyMetadataSelectByValue(mode));
        });
        ['metadata-select-value', 'metadata-select-value2'].forEach(inputId => {
            document.getElementById(inputId).addEventListener('keydown', event => {
                if (event.key === 'Enter') {
                    event.preventDefault();
                    // Enter runs the widening action; the two narrowing ones are
                    // destructive enough to be worth an explicit click.
                    applyMetadataSelectByValue('add');
                }
            });
        });
        document.getElementById('metadata-rename-apply').addEventListener('click', applyColumnRename);
        document.getElementById('metadata-rename-value').addEventListener('keydown', event => {
            if (event.key === 'Enter') {
                event.preventDefault();
                applyColumnRename();
            }
        });

        const tbody = document.querySelector('#metadata-table tbody');
        tbody.addEventListener('dblclick', event => {
            const cell = event.target.closest('td[data-column-key]');
            if (cell) {
                event.preventDefault();
                beginMetadataCellEdit(cell);
            }
        });

        document.getElementById('metadata-add-column').addEventListener('click', () => {
            const nameInput = document.getElementById('metadata-new-column-name');
            if (addMetadataColumn(nameInput.value, document.getElementById('metadata-new-column-type').value)) {
                nameInput.value = '';
            }
        });
        document.getElementById('metadata-add-cluster-column').addEventListener('click', () => {
            const nameInput = document.getElementById('metadata-new-column-name');
            if (addClusterColumn(nameInput.value)) {
                nameInput.value = '';
            }
        });
        document.getElementById('metadata-new-column-name').addEventListener('keydown', event => {
            if (event.key === 'Enter') {
                event.preventDefault();
                document.getElementById('metadata-add-column').click();
            }
        });
        document.getElementById('metadata-delete-column-apply').addEventListener('click', () => {
            deleteMetadataColumn(document.getElementById('metadata-delete-column').value);
        });
        document.getElementById('metadata-fill-apply').addEventListener('click', applyBulkFill);
        document.getElementById('metadata-fill-value').addEventListener('keydown', event => {
            if (event.key === 'Enter') {
                event.preventDefault();
                applyBulkFill();
            }
        });
        document.getElementById('metadata-paste-open').addEventListener('click', openMetadataPasteDialog);
        document.getElementById('metadata-paste-cancel').addEventListener('click', closeMetadataPasteDialog);
        document.getElementById('metadata-paste-apply').addEventListener('click', applyPastedColumn);
        document.getElementById('metadata-paste-overlay').addEventListener('click', event => {
            if (event.target === event.currentTarget) { closeMetadataPasteDialog(); }
        });
        document.addEventListener('keydown', event => {
            if (event.key === 'Escape' && metadataPasteDialogIsOpen()) {
                closeMetadataPasteDialog();
            }
        });
    }
"""


def _gradient_stops_js() -> str:
    """Editable gradient stops: two mandatory ends plus any number between.

    Plain (non-f) string like `_layout_worker_js`; see `_session_state_js`.
    """
    return r"""
    // ---------------------------------------------------------------------
    // Gradient stops
    //
    // A numeric palette is an ordered list of {value, color} stops: two
    // mandatory ends plus any number of intermediates. Coloring interpolates
    // between the bracketing pair and holds the end colors flat outside the
    // ends, so two stops give a plain ramp and three reproduce the old
    // low/mid/high midpoint exactly.
    //
    // `state.gradientStops` is the working list the dialog edits; committing it
    // writes a sorted copy into the column's stored palette. The histogram, the
    // gradient bar and the slider all read that one list, so they cannot drift
    // apart, and every edit -- typing, dragging a knob, adding or removing a
    // stop -- funnels through commitGradientStops().
    // ---------------------------------------------------------------------

    const MIN_GRADIENT_STOPS = 2;
    // A cap, not a limitation of the model: past a dozen the rows stop fitting
    // in the dialog and the knobs start overlapping on the track.
    const MAX_GRADIENT_STOPS = 12;
    const DEFAULT_NO_VALUE_COLOR = '#b3a89d';

    function sortedGradientStops(stops) {
        return stops.slice().sort((left, right) => left.value - right.value);
    }

    function copyGradientStops(stops) {
        return stops.map(stop => ({value: stop.value, color: stop.color}));
    }

    // Read a stored palette as stops. Sessions written before stops existed
    // carry {low, mid, high} with separate lowValue/midValue/highValue bounds;
    // converting here means nothing downstream has to know about that shape.
    function numericPaletteStops(palette, minValue, maxValue) {
        if (!palette) {
            return null;
        }
        if (Array.isArray(palette.stops) && palette.stops.length >= MIN_GRADIENT_STOPS) {
            return sortedGradientStops(palette.stops
                .filter(stop => stop && Number.isFinite(Number(stop.value)) && stop.color)
                .map(stop => ({value: Number(stop.value), color: stop.color})));
        }
        if (!palette.low || !palette.high) {
            return null;
        }
        const low = (palette.lowValue !== null && palette.lowValue !== undefined)
            ? palette.lowValue : minValue;
        const high = (palette.highValue !== null && palette.highValue !== undefined)
            ? palette.highValue : maxValue;
        const stops = [{value: low, color: palette.low}, {value: high, color: palette.high}];
        if (palette.mid) {
            const mid = (palette.midValue !== null && palette.midValue !== undefined)
                ? palette.midValue : ((low + high) / 2);
            stops.push({value: mid, color: palette.mid});
        }
        return sortedGradientStops(stops);
    }

    // Rewrite any legacy numeric palettes in place, so a session saved by an
    // older viewer keeps its colors instead of silently reverting to defaults.
    function normalizeCustomPalettes(palettes) {
        Object.keys(palettes || {}).forEach(key => {
            const palette = palettes[key];
            if (!palette || palette.type !== 'numeric' || Array.isArray(palette.stops)) {
                return;
            }
            const stops = numericPaletteStops(palette, 0, 1);
            palettes[key] = {
                type: 'numeric',
                stops: stops || [],
                nullColor: palette.nullColor || null,
            };
        });
        return palettes;
    }

    function currentGradientStops() {
        return Array.isArray(state.gradientStops) ? state.gradientStops : [];
    }

    // ---------------------------------------------------------------------
    // The shared axis
    //
    // The histogram's binned extent, not the outermost stops: narrowing the
    // ramp then moves the knobs inward instead of rescaling the chart under
    // them.
    // ---------------------------------------------------------------------
    function gradientRangeDomain() {
        const histogram = state.colorHistogram;
        if (histogram && histogram.highEdge > histogram.lowEdge) {
            return {low: histogram.lowEdge, high: histogram.highEdge};
        }
        const info = colorInfo(currentColorField());
        if (info && Number.isFinite(info.min) && Number.isFinite(info.max) && info.max > info.min) {
            return {low: info.min, high: info.max};
        }
        return null;
    }

    // Where a value sits on the axis, as a percentage. Deliberately unclamped:
    // CSS gradient stops accept out-of-range percentages and extend the ramp
    // past the box, which is exactly how a stop outside the data range should
    // render. Knob placement clamps separately.
    function gradientRangePercent(value, domain) {
        return ((value - domain.low) / (domain.high - domain.low)) * 100;
    }

    function gradientValueAtPercent(percent, domain) {
        return domain.low + ((percent / 100) * (domain.high - domain.low));
    }

    // Integer columns snap to whole numbers; float columns round to a fixed
    // number of decimals chosen from the span, so dragging writes 3.42 rather
    // than 3.4166666666666665 into the number input.
    function roundGradientValue(value, domain) {
        const histogram = state.colorHistogram;
        if (histogram && histogram.integer) {
            return Math.round(value);
        }
        const span = domain.high - domain.low;
        const decimals = span > 0
            ? Math.max(0, Math.min(6, 3 - Math.floor(Math.log10(span))))
            : 2;
        return Number(value.toFixed(decimals));
    }

    function gradientKnobStep(domain) {
        const histogram = state.colorHistogram;
        if (histogram && histogram.integer) {
            return 1;
        }
        return (domain.high - domain.low) / 100;
    }

    // A stop is held between its neighbors so the ramp cannot fold over during
    // a drag. Typing into the number inputs is left unconstrained -- values are
    // re-sorted on commit instead, so a stop can be typed past its neighbor.
    function clampGradientStopValue(index, value, domain) {
        const stops = currentGradientStops();
        let lowLimit = domain.low;
        let highLimit = domain.high;
        if (index > 0) {
            lowLimit = Math.max(lowLimit, stops[index - 1].value);
        }
        if (index < stops.length - 1) {
            highLimit = Math.min(highLimit, stops[index + 1].value);
        }
        if (highLimit < lowLimit) {
            highLimit = lowLimit;
        }
        return Math.max(lowLimit, Math.min(highLimit, value));
    }

    // ---------------------------------------------------------------------
    // Committing an edit
    // ---------------------------------------------------------------------
    function commitGradientStops() {
        const columnName = currentColorField();
        if (!columnName) {
            return;
        }
        const stored = ensureNumericPalette(columnName);
        // Stored sorted so numericColor -- which runs once per node -- can walk
        // the stops without sorting them itself.
        stored.stops = copyGradientStops(sortedGradientStops(currentGradientStops()));
        stored.nullColor = document.getElementById('color-null').value;
        updateGradientPreview();
        drawColorHistogram(columnName);
        rebuildNodeColorCache();
        scheduleClusterRender();
    }

    function gradientStopRole(index, total) {
        if (index === 0) {
            return 'Min';
        }
        if (index === total - 1) {
            return 'Max';
        }
        return 'Stop ' + (index + 1);
    }

    function renderGradientStopRows() {
        const stops = currentGradientStops();
        const list = document.getElementById('color-stop-list');
        list.innerHTML = stops.map((stop, index) => {
            const role = gradientStopRole(index, stops.length);
            const removable = index > 0 && index < stops.length - 1;
            return '<div class="cp-stop-row" data-stop-index="' + index + '">'
                + '<span class="cp-stop-role">' + role + '</span>'
                + '<input type="color" data-stop-field="color" value="' + stop.color
                + '" aria-label="' + role + ' color" />'
                + '<input type="text" class="cp-stop-hex" data-stop-field="hex"'
                + ' value="' + normalizeColorHex(stop.color) + '" spellcheck="false"'
                + ' autocapitalize="off" autocomplete="off"'
                + ' aria-label="' + role + ' hex code" />'
                + '<input type="number" step="any" data-stop-field="value" value="' + stop.value
                + '" aria-label="' + role + ' value" />'
                + (removable
                    ? '<button type="button" class="cp-stop-remove" data-stop-remove'
                        + ' aria-label="Remove ' + role + '" title="Remove this stop">×</button>'
                    : '<span class="cp-stop-remove-gap"></span>')
                + '</div>';
        }).join('');
        const addButton = document.getElementById('color-add-stop');
        addButton.disabled = stops.length >= MAX_GRADIENT_STOPS;
        addButton.title = addButton.disabled
            ? ('At the limit of ' + MAX_GRADIENT_STOPS + ' stops.')
            : 'Add a stop in the widest gap, in the color the ramp already has there';
    }

    function addGradientStop() {
        const stops = currentGradientStops();
        const domain = gradientRangeDomain();
        if (stops.length < MIN_GRADIENT_STOPS || !domain) {
            return;
        }
        if (stops.length >= MAX_GRADIENT_STOPS) {
            setStatus('A gradient holds at most ' + MAX_GRADIENT_STOPS + ' stops.');
            return;
        }
        // Split the widest gap: with no other information that is where an extra
        // stop buys the most control.
        let gapIndex = 1;
        let widest = -Infinity;
        for (let index = 1; index < stops.length; index++) {
            const gap = stops[index].value - stops[index - 1].value;
            if (gap > widest) {
                widest = gap;
                gapIndex = index;
            }
        }
        const lower = stops[gapIndex - 1];
        const upper = stops[gapIndex];
        const midpoint = (lower.value + upper.value) / 2;
        const rounded = roundGradientValue(midpoint, domain);
        // Rounding an integer column would land the new stop on top of a
        // neighbor when the gap is a single unit; keep the exact midpoint there.
        const value = (rounded > lower.value && rounded < upper.value) ? rounded : midpoint;
        stops.splice(gapIndex, 0, {
            // Taking the color the ramp already has here means adding a stop
            // changes nothing until it is moved: you shape a ramp, not reset it.
            value,
            color: lerpHexColor(lower.color, upper.color, 0.5),
        });
        renderGradientStopRows();
        commitGradientStops();
        setStatus('Added a gradient stop at ' + formatValue(value) + '.');
    }

    function removeGradientStop(index) {
        const stops = currentGradientStops();
        if (index <= 0 || index >= stops.length - 1) {
            return;   // the two ends define the ramp and are not removable
        }
        const removed = stops.splice(index, 1)[0];
        renderGradientStopRows();
        commitGradientStops();
        setStatus('Removed the gradient stop at ' + formatValue(removed.value) + '.');
    }

    // "Reset values to data range": refit the ramp onto the data while keeping
    // the relative spacing of any stops the user has placed.
    function resetGradientStopValues() {
        const info = colorInfo(currentColorField());
        const stops = currentGradientStops();
        const domain = gradientRangeDomain();
        if (!info || stops.length < MIN_GRADIENT_STOPS) {
            return;
        }
        const first = stops[0].value;
        const span = stops[stops.length - 1].value - first;
        stops.forEach((stop, index) => {
            // Even spacing is the only sensible reading when every stop
            // currently sits on the same value.
            const fraction = span > 0
                ? (stop.value - first) / span
                : index / (stops.length - 1);
            const value = info.min + (fraction * (info.max - info.min));
            stop.value = domain ? roundGradientValue(value, domain) : value;
        });
        renderGradientStopRows();
        commitGradientStops();
    }

    function setGradientStopHexField(index, color) {
        const field = document.querySelector(
            '.cp-stop-row[data-stop-index="' + index + '"] [data-stop-field="hex"]');
        if (field) {
            field.value = normalizeColorHex(color) || color;
        }
    }

    function setupGradientStopEditing() {
        const list = document.getElementById('color-stop-list');
        const stopIndexOf = target => {
            const row = target.closest('[data-stop-index]');
            return row ? Number(row.dataset.stopIndex) : -1;
        };
        list.addEventListener('input', event => {
            const field = event.target.dataset ? event.target.dataset.stopField : null;
            const index = stopIndexOf(event.target);
            const stop = currentGradientStops()[index];
            if (!field || !stop) {
                return;
            }
            if (field === 'color') {
                stop.color = event.target.value;
                setGradientStopHexField(index, stop.color);
            } else if (field === 'hex') {
                const hex = normalizeColorHex(event.target.value);
                if (hex === null) {
                    return;   // mid-edit ("#", "#1a2"); the change handler has the say
                }
                // Stored lower case, which is the one form the color input and
                // lerpHexColor both produce; the field itself shows the canonical
                // upper-case spelling.
                stop.color = hex.toLowerCase();
                const swatch = document.querySelector(
                    '.cp-stop-row[data-stop-index="' + index + '"] [data-stop-field="color"]');
                if (swatch) { swatch.value = stop.color; }
            } else {
                const value = parseNumberOrNull(event.target.value);
                if (value === null) {
                    return;   // mid-edit ("", "-", "1e"); wait for something numeric
                }
                stop.value = value;
            }
            commitGradientStops();
        });
        // Leaving the hex field rewrites whatever was typed into the canonical
        // spelling; text that is not a color at all reverts to the stop's color
        // rather than being silently dropped or left sitting there looking valid.
        list.addEventListener('change', event => {
            if (!event.target.dataset || event.target.dataset.stopField !== 'hex') {
                return;
            }
            const index = stopIndexOf(event.target);
            const stop = currentGradientStops()[index];
            if (!stop) {
                return;
            }
            const typed = event.target.value.trim();
            // Clearing the field and tabbing away reads as canceling the edit, so
            // that reverts quietly; anything else that is not a color is worth saying.
            if (typed !== '' && normalizeColorHex(typed) === null) {
                setStatus('"' + typed + '" is not a hex color (expected #RGB or #RRGGBB).');
            }
            setGradientStopHexField(index, stop.color);
        });
        // Re-sorting mid-keystroke would yank the row out from under the cursor,
        // so a stop typed past its neighbor is only reordered on commit.
        list.addEventListener('change', event => {
            if (!event.target.dataset || event.target.dataset.stopField !== 'value') {
                return;
            }
            const stops = currentGradientStops();
            const sorted = sortedGradientStops(stops);
            if (sorted.some((stop, index) => stop !== stops[index])) {
                state.gradientStops = sorted;
                renderGradientStopRows();
                commitGradientStops();
            }
        });
        list.addEventListener('click', event => {
            const button = event.target.closest('[data-stop-remove]');
            if (button) {
                removeGradientStop(stopIndexOf(button));
            }
        });
        document.getElementById('color-add-stop').addEventListener('click', addGradientStop);
    }

    // ---------------------------------------------------------------------
    // The slider: one knob per stop
    // ---------------------------------------------------------------------
    function updateGradientSlider() {
        const slider = document.getElementById('color-range-slider');
        const domain = gradientRangeDomain();
        const stops = currentGradientStops();
        if (!domain || stops.length < MIN_GRADIENT_STOPS) {
            slider.hidden = true;
            return;
        }
        slider.hidden = false;
        let knobs = slider.querySelectorAll('.cp-range-knob');
        if (knobs.length !== stops.length) {
            // Rebuilt only when the count changes, so a drag is never
            // interrupted by its own knob being replaced underneath it.
            knobs.forEach(knob => knob.remove());
            stops.forEach((stop, index) => {
                const knob = document.createElement('div');
                knob.className = 'cp-range-knob';
                knob.dataset.stopIndex = String(index);
                knob.setAttribute('role', 'slider');
                knob.setAttribute('tabindex', '0');
                slider.appendChild(knob);
            });
            knobs = slider.querySelectorAll('.cp-range-knob');
        }
        knobs.forEach((knob, index) => {
            const stop = stops[index];
            const role = gradientStopRole(index, stops.length);
            const isEnd = index === 0 || index === stops.length - 1;
            knob.classList.toggle('cp-range-knob-mid', !isEnd);
            knob.style.left = Math.max(0, Math.min(100, gradientRangePercent(stop.value, domain))) + '%';
            knob.style.background = stop.color;
            knob.setAttribute('aria-label', role + ' value');
            knob.setAttribute('aria-valuemin', String(domain.low));
            knob.setAttribute('aria-valuemax', String(domain.high));
            knob.setAttribute('aria-valuenow', String(stop.value));
            knob.setAttribute('aria-valuetext', formatValue(stop.value));
        });
        // Shade the part of the axis the ramp actually spans.
        const span = document.getElementById('color-range-span');
        const lowPercent = Math.max(0, Math.min(100, gradientRangePercent(stops[0].value, domain)));
        const highPercent = Math.max(0, Math.min(100,
            gradientRangePercent(stops[stops.length - 1].value, domain)));
        span.hidden = false;
        span.style.left = Math.min(lowPercent, highPercent) + '%';
        span.style.width = Math.abs(highPercent - lowPercent) + '%';
    }

    // Move a stop from the slider and push the result through the same commit
    // path as typing, so there is one way for a gradient to change.
    function setGradientStopFromSlider(index, rawValue, domain) {
        const stop = currentGradientStops()[index];
        if (!stop) {
            return;
        }
        const value = roundGradientValue(clampGradientStopValue(index, rawValue, domain), domain);
        if (value === stop.value) {
            return;
        }
        stop.value = value;
        const input = document.querySelector(
            '.cp-stop-row[data-stop-index="' + index + '"] [data-stop-field="value"]');
        if (input) {
            input.value = value;
        }
        commitGradientStops();
    }

    function setupGradientRangeSlider() {
        const slider = document.getElementById('color-range-slider');
        slider.addEventListener('pointerdown', event => {
            const knob = event.target.closest('.cp-range-knob');
            const domain = gradientRangeDomain();
            if (!knob || !domain) {
                return;
            }
            const index = Number(knob.dataset.stopIndex);
            event.preventDefault();
            knob.focus();
            // Pointer capture keeps the drag alive past the knob's 18px, and
            // past the dialog edge, without a document-level move listener.
            knob.setPointerCapture(event.pointerId);
            const drag = moveEvent => {
                const rect = slider.getBoundingClientRect();
                if (rect.width <= 0) {
                    return;
                }
                const percent = ((moveEvent.clientX - rect.left) / rect.width) * 100;
                setGradientStopFromSlider(index, gradientValueAtPercent(percent, domain), domain);
            };
            const stop = () => {
                knob.removeEventListener('pointermove', drag);
                knob.removeEventListener('pointerup', stop);
                knob.removeEventListener('pointercancel', stop);
            };
            knob.addEventListener('pointermove', drag);
            knob.addEventListener('pointerup', stop);
            knob.addEventListener('pointercancel', stop);
            drag(event);
        });
        slider.addEventListener('keydown', event => {
            const knob = event.target.closest('.cp-range-knob');
            const domain = gradientRangeDomain();
            if (!knob || !domain) {
                return;
            }
            const index = Number(knob.dataset.stopIndex);
            const stop = currentGradientStops()[index];
            if (!stop) {
                return;
            }
            const step = gradientKnobStep(domain);
            let next = null;
            if (event.key === 'ArrowLeft' || event.key === 'ArrowDown') { next = stop.value - step; }
            else if (event.key === 'ArrowRight' || event.key === 'ArrowUp') { next = stop.value + step; }
            else if (event.key === 'PageDown') { next = stop.value - (step * 10); }
            else if (event.key === 'PageUp') { next = stop.value + (step * 10); }
            else if (event.key === 'Home') { next = domain.low; }
            else if (event.key === 'End') { next = domain.high; }
            else if (event.key === 'Delete' || event.key === 'Backspace') {
                event.preventDefault();
                removeGradientStop(index);
                return;
            }
            if (next === null) {
                return;
            }
            event.preventDefault();
            setGradientStopFromSlider(index, next, domain);
        });
    }
"""


def _extraction_js() -> str:
    """"Save extraction": a new bundle over just the selected nodes.

    Plain (non-f) string like `_layout_worker_js`; see `_session_state_js`.
    The hierarchy and merge-series maths here are ports of `ssn_hierarchy.py`
    (`build_mst_component_hierarchy`, `component_size_summary_by_threshold`,
    `threshold_merge_event_rows`, `merge_event_moving_sum`,
    `filter_merge_event_rows`) and must stay in step with them.
    """
    return r"""
    // ---------------------------------------------------------------------
    // Save extraction
    //
    // Writes a self-contained bundle over just the selected nodes.
    //
    // Nodes that the MST joins must stay joined in the selection, wherever the
    // slider happens to sit. The MST is a lossy summary: it kept one path
    // between any two nodes and discarded every other edge. So if the selection
    // drops the nodes along that path, the bundle holds no evidence of how the
    // remaining pieces relate, and writing them out as separate components
    // would assert an absence of similarity that only the original all-vs-all
    // matrix could establish. Subset that matrix with
    // build_ssn_viewer.py --subset instead, which measures the relationship
    // rather than guessing at it.
    //
    // Nodes in *different* components of the original network are exempt: they
    // were already unrelated before the extraction, so keeping them apart
    // invents nothing. That is what lets a whole multi-component network be
    // extracted, and it is checked per original component rather than globally.
    // ---------------------------------------------------------------------
    const MOVING_SUM_WINDOW_FRACTION = 0.05;
    const MOVING_SUM_GRID_POINTS = 800;
    const DEFAULT_MAX_MERGE_EVENTS = 500;
    // The viewer's own default. It picks its events from the chart's current window, so
    // a much smaller number shows more: zoom in and the cap re-spends itself on the
    // stretch being looked at. See ssn_hierarchy.DEFAULT_WINDOWED_MAX_MERGE_EVENTS.
    const DEFAULT_WINDOWED_MAX_MERGE_EVENTS = 50;
    const MERGE_EVENT_DENSITY_BINS = 20;

    function makeUnionFind(size) {
        const parent = new Int32Array(size);
        for (let i = 0; i < size; i++) { parent[i] = i; }
        function find(index) {
            let root = index;
            while (parent[root] !== root) { root = parent[root]; }
            while (parent[index] !== index) {
                const next = parent[index];
                parent[index] = root;
                index = next;
            }
            return root;
        }
        return {find, union: (a, b) => { parent[find(a)] = find(b); }};
    }

    // Port of ssn_hierarchy.format_threshold_value.
    function formatThresholdValue(threshold) {
        return Number.isFinite(threshold) ? Number(threshold).toFixed(2) : '∞';
    }

    // The selected nodes plus the MST edges with both ends inside, re-indexed into
    // the extraction. `pieces` is how many MST-connected components the selection
    // falls into; anything above 1 is what makes an extraction unfaithful.
    function extractionInducedGraph() {
        const nodeIndices = Array.from(state.selectedNodeIndices).sort((left, right) => left - right);
        const newIndexByOld = new Map(nodeIndices.map((oldIndex, newIndex) => [oldIndex, newIndex]));
        const edges = [];
        // Original order is weight-descending, which every consumer below relies on.
        state.bundle.graph.mst_edges.forEach(edge => {
            const left = newIndexByOld.get(edge[0]);
            const right = newIndexByOld.get(edge[1]);
            if (left !== undefined && right !== undefined) {
                edges.push([left, right, edge[2]]);
            }
        });

        const unionFind = makeUnionFind(nodeIndices.length);
        edges.forEach(edge => unionFind.union(edge[0], edge[1]));

        // Split each original network component's share of the selection into the
        // pieces it fell into. A component contributing more than one piece is the
        // failure case: the MST joined those nodes and the selection broke them up.
        const originalRoot = originalComponentByNode();
        const piecesByComponent = new Map();
        nodeIndices.forEach((oldIndex, newIndex) => {
            const component = originalRoot[oldIndex];
            if (!piecesByComponent.has(component)) { piecesByComponent.set(component, new Set()); }
            piecesByComponent.get(component).add(unionFind.find(newIndex));
        });
        let brokenComponents = 0;
        let brokenPieces = 0;
        piecesByComponent.forEach(pieces => {
            if (pieces.size > 1) {
                brokenComponents += 1;
                brokenPieces += pieces.size;
            }
        });
        return {nodeIndices, newIndexByOld, edges, brokenComponents, brokenPieces};
    }

    // Which top-level component of the full MST forest each node belongs to.
    function originalComponentByNode() {
        const hierarchy = state.bundle.graph.hierarchy;
        const componentOf = new Int32Array(state.bundle.graph.nodes.length).fill(-1);
        hierarchy.roots.forEach(rootId => {
            const root = hierarchy.nodes[rootId];
            for (let i = root.leaf_start; i < root.leaf_start + root.leaf_count; i++) {
                componentOf[hierarchy.leaf_order[i]] = rootId;
            }
        });
        return componentOf;
    }

    // Port of ssn_hierarchy.build_mst_component_hierarchy.
    function buildExtractionHierarchy(nodeCount, edges) {
        const unionFind = makeUnionFind(nodeCount);
        const nodes = [];
        const componentIdByRoot = new Map();
        const sizeById = [];
        const minLeafById = [];
        for (let nodeIndex = 0; nodeIndex < nodeCount; nodeIndex++) {
            nodes.push({id: nodeIndex, kind: 'leaf', node_index: nodeIndex, size: 1, parent: null});
            componentIdByRoot.set(nodeIndex, nodeIndex);
            sizeById.push(1);
            minLeafById.push(nodeIndex);
        }

        edges.forEach(edge => {
            const leftRoot = unionFind.find(edge[0]);
            const rightRoot = unionFind.find(edge[1]);
            if (leftRoot === rightRoot) { return; }
            let leftId = componentIdByRoot.get(leftRoot);
            let rightId = componentIdByRoot.get(rightRoot);
            if (minLeafById[leftId] > minLeafById[rightId]) {
                const swap = leftId; leftId = rightId; rightId = swap;
            }
            const componentId = nodes.length;
            nodes[leftId].parent = componentId;
            nodes[rightId].parent = componentId;
            nodes.push({
                id: componentId,
                kind: 'cluster',
                left: leftId,
                right: rightId,
                threshold: edge[2],
                size: sizeById[leftId] + sizeById[rightId],
                parent: null,
            });
            sizeById.push(sizeById[leftId] + sizeById[rightId]);
            minLeafById.push(Math.min(minLeafById[leftId], minLeafById[rightId]));
            unionFind.union(leftRoot, rightRoot);
            componentIdByRoot.delete(leftRoot);
            componentIdByRoot.set(unionFind.find(rightRoot), componentId);
        });

        const roots = Array.from(componentIdByRoot.values())
            .sort((left, right) => nodes[right].size - nodes[left].size || left - right);

        const leafOrder = [];
        const stack = [];
        for (let i = roots.length - 1; i >= 0; i--) { stack.push([roots[i], false]); }
        while (stack.length > 0) {
            const [componentId, visited] = stack.pop();
            const node = nodes[componentId];
            if (node.kind === 'leaf') {
                node.leaf_start = leafOrder.length;
                node.leaf_count = 1;
                leafOrder.push(node.node_index);
                continue;
            }
            if (visited) {
                node.leaf_start = nodes[node.left].leaf_start;
                node.leaf_count = nodes[node.left].leaf_count + nodes[node.right].leaf_count;
                continue;
            }
            stack.push([componentId, true]);
            stack.push([node.right, false]);
            stack.push([node.left, false]);
        }
        return {nodes, roots, leaf_order: leafOrder};
    }

    // Port of ssn_hierarchy.component_size_summary_by_threshold followed by
    // threshold_merge_event_rows: one row per distinct MST edge weight.
    function mergeEventRows(nodeCount, edges, metric) {
        if (nodeCount === 0 || edges.length === 0) { return []; }
        const unionFind = makeUnionFind(nodeCount);
        const componentSizes = new Int32Array(nodeCount).fill(1);
        // Sizes only ever grow, so the running maximum is the largest component.
        let largest = nodeCount > 0 ? 1 : 0;
        let nonSingletonCount = 0;
        let nonSingletonSum = 0;
        const summary = [{threshold: Infinity, largest, avg: 0, impact: 0}];

        edges.forEach(edge => {
            const leftRoot = unionFind.find(edge[0]);
            const rightRoot = unionFind.find(edge[1]);
            let impact = 0;
            if (leftRoot !== rightRoot) {
                const leftSize = componentSizes[leftRoot];
                const rightSize = componentSizes[rightRoot];
                impact = metric === 'product' ? leftSize * rightSize : Math.min(leftSize, rightSize);
                if (leftSize > 1) { nonSingletonCount -= 1; nonSingletonSum -= leftSize; }
                if (rightSize > 1) { nonSingletonCount -= 1; nonSingletonSum -= rightSize; }
                const mergedSize = leftSize + rightSize;
                unionFind.union(leftRoot, rightRoot);
                componentSizes[unionFind.find(rightRoot)] = mergedSize;
                largest = Math.max(largest, mergedSize);
                nonSingletonCount += 1;
                nonSingletonSum += mergedSize;
            }
            summary.push({
                threshold: edge[2],
                largest,
                avg: nonSingletonCount > 0 ? nonSingletonSum / nonSingletonCount : 0,
                impact,
            });
        });

        const rows = [];
        let previous = summary[0];
        let rowIndex = 1;
        while (rowIndex < summary.length) {
            const firstRowIndex = rowIndex;
            const thresholdValue = summary[rowIndex].threshold;
            let mergeImpact = 0;
            let mergeCount = 0;
            let largestMerge = 0;
            const mergeSizeCounts = {};
            let lastRow = summary[rowIndex];
            while (rowIndex < summary.length && summary[rowIndex].threshold === thresholdValue) {
                const impact = summary[rowIndex].impact;
                mergeImpact += impact;
                if (impact > 0) {
                    // A zero impact means both ends were already in one component.
                    const key = Number.isInteger(impact) ? String(impact) : String(impact);
                    mergeSizeCounts[key] = (mergeSizeCounts[key] || 0) + 1;
                    mergeCount += 1;
                    largestMerge = Math.max(largestMerge, impact);
                }
                lastRow = summary[rowIndex];
                rowIndex += 1;
            }
            rows.push({
                edge_index: firstRowIndex - 2,
                summary_row_from: firstRowIndex,
                summary_row_to: rowIndex - 1,
                threshold_from_value: previous.threshold,
                threshold_from: formatThresholdValue(previous.threshold),
                threshold_to: formatThresholdValue(lastRow.threshold),
                threshold_value: lastRow.threshold,
                merge_impact: mergeImpact,
                merge_size_counts: mergeSizeCounts,
                largest_merge: largestMerge,
                merge_count: mergeCount,
                delta_largest: Math.abs(lastRow.largest - previous.largest),
                delta_avg_non_singleton: Math.abs(lastRow.avg - previous.avg),
            });
            previous = lastRow;
        }
        return rows;
    }

    // Port of ssn_hierarchy.merge_event_rank_key: strongest first.
    function compareMergeEventRank(left, right) {
        return right.merge_impact - left.merge_impact ||
            right.delta_largest - left.delta_largest ||
            right.delta_avg_non_singleton - left.delta_avg_non_singleton ||
            left.edge_index - right.edge_index;
    }

    // Port of ssn_hierarchy.merge_event_density_bin.
    function mergeEventDensityBin(thresholdValue, lo, hi, densityBins) {
        if (densityBins < 1 || !Number.isFinite(lo) || !Number.isFinite(hi) || hi <= lo) {
            return -1;
        }
        const position = Math.floor(((Number(thresholdValue) - lo) / (hi - lo)) * densityBins);
        return Math.max(0, Math.min(densityBins - 1, position));
    }

    // Port of ssn_hierarchy.filter_merge_event_rows: the strongest maxMergeEvents events
    // *within `window`*, plus the strongest event in each otherwise-empty band of the
    // whole axis, plus `pinnedThreshold`'s event. See that function for the reasoning;
    // in short, the window is what makes a small cap show more rather than less, the
    // bands keep every stretch of the slider reachable, and the pin keeps the cut
    // currently in effect from being filtered out from under the user.
    function filterMergeEventRows(rows, maxMergeEvents = DEFAULT_MAX_MERGE_EVENTS,
                                  densityBins = MERGE_EVENT_DENSITY_BINS,
                                  window = null, pinnedThreshold = null) {
        if (maxMergeEvents === 0 || rows.length <= maxMergeEvents) { return rows.slice(); }
        const ranked = rows.slice().sort(compareMergeEventRank);

        let pool = ranked;
        if (window) {
            const low = Math.min(window.min, window.max);
            const high = Math.max(window.min, window.max);
            pool = ranked.filter(row => {
                const value = Number(row.threshold_value);
                return value >= low && value <= high;
            });
        }
        const filtered = pool.slice(0, maxMergeEvents);
        const selected = new Set(filtered);

        // Band edges come from the whole series, because that is what the axis spans --
        // a window narrows what the cap ranks, never what the axis has to cover.
        // Scanned rather than spread through Math.min: `rows` is the UNFILTERED series,
        // which on a large network runs to six figures and would overflow the argument
        // list.
        let lo = Infinity;
        let hi = -Infinity;
        let finiteCount = 0;
        for (const row of rows) {
            const value = Number(row.threshold_value);
            if (!Number.isFinite(value)) { continue; }
            finiteCount += 1;
            if (value < lo) { lo = value; }
            if (value > hi) { hi = value; }
        }
        if (finiteCount > 0 && densityBins > 0) {
            const covered = new Set(filtered.map(row =>
                mergeEventDensityBin(row.threshold_value, lo, hi, densityBins)));
            // `ranked` is strongest-first, so the first row seen in an empty band is the
            // strongest one available to represent it. Candidates come from the whole
            // series, not the window: a band outside the window still needs its stop.
            for (const row of ranked) {
                if (selected.has(row)) { continue; }
                const bin = mergeEventDensityBin(row.threshold_value, lo, hi, densityBins);
                if (bin < 0 || covered.has(bin)) { continue; }
                covered.add(bin);
                filtered.push(row);
                selected.add(row);
            }
        }

        if (pinnedThreshold !== null && Number.isFinite(pinnedThreshold)
                && !filtered.some(row => Number(row.threshold_value) === pinnedThreshold)) {
            const pinned = ranked.find(row => Number(row.threshold_value) === pinnedThreshold);
            if (pinned) { filtered.push(pinned); }
        }

        return filtered.sort((left, right) => left.edge_index - right.edge_index);
    }

    // Port of ssn_hierarchy.merge_event_moving_sum. Runs over the UNFILTERED rows:
    // filtering keeps the strongest events, so summing afterwards would undercount
    // exactly the small ones this series exists to show.
    function mergeEventMovingSum(rows) {
        const empty = {window: 0, x: [], y: []};
        if (rows.length === 0) { return empty; }
        const points = rows
            .filter(row => Number.isFinite(row.threshold_value) && Number.isFinite(row.merge_impact))
            .map(row => ({threshold: row.threshold_value, impact: row.merge_impact}))
            .sort((left, right) => left.threshold - right.threshold);
        if (points.length === 0) { return empty; }
        const lo = points[0].threshold;
        const hi = points[points.length - 1].threshold;
        if (hi === lo) { return empty; }

        const window = MOVING_SUM_WINDOW_FRACTION * (hi - lo);
        const halfWindow = window / 2;
        const sortedThresholds = points.map(point => point.threshold);
        const cumulative = [0];
        points.forEach(point => cumulative.push(cumulative[cumulative.length - 1] + point.impact));

        // Inclusive at both ends, matching numpy searchsorted 'left'/'right'.
        const lowerBound = value => {
            let low = 0, high = sortedThresholds.length;
            while (low < high) {
                const mid = (low + high) >> 1;
                if (sortedThresholds[mid] < value) { low = mid + 1; } else { high = mid; }
            }
            return low;
        };
        const upperBound = value => {
            let low = 0, high = sortedThresholds.length;
            while (low < high) {
                const mid = (low + high) >> 1;
                if (sortedThresholds[mid] <= value) { low = mid + 1; } else { high = mid; }
            }
            return low;
        };

        const x = [];
        const y = [];
        for (let i = 0; i < MOVING_SUM_GRID_POINTS; i++) {
            const gridValue = lo + ((hi - lo) * i) / (MOVING_SUM_GRID_POINTS - 1);
            x.push(gridValue);
            y.push(cumulative[upperBound(gridValue + halfWindow)] - cumulative[lowerBound(gridValue - halfWindow)]);
        }
        return {window, x, y};
    }

    // Port of ssn_hierarchy.threshold_slider_stops, including its trailing floor stop:
    // every other stop excludes its own tie group, so without one strictly below the
    // weakest edge the fully merged network cannot be reached on the slider.
    function buildSliderStops(rows, edges) {
        const stops = [{edge_index: -1, threshold_label: '\u221e', threshold_value: null}];
        rows.forEach(row => {
            stops.push({
                edge_index: row.edge_index,
                threshold_label: row.threshold_to,
                threshold_value: row.threshold_value,
            });
        });
        if (edges.length > 0) {
            // edges are weight-descending.
            const floorValue = floorThresholdValue(edges[edges.length - 1][2], edges[0][2]);
            stops.push({
                edge_index: edges.length - 1,   // every MST edge is strictly above
                threshold_label: formatThresholdValue(floorValue),
                threshold_value: floorValue,
            });
        }
        return stops;
    }

    // The window-independent half of the series: the full event list and its moving sum,
    // derived from the merge order the bundle carries. From bundle v6 the file stores
    // none of it -- it is a pure function of (node count, mst_edges, metric), and
    // deriving it here is what lets both the cap and the chart's window decide what is
    // shown, neither of which is knowable when the file is written. Pre-v6 files still
    // carry a capped copy; it is ignored in favor of this, so an old bundle is not
    // stuck with whatever cap it happened to be built with.
    //
    // This is the expensive half -- one union-find replay over the MST edges -- so it
    // runs once per bundle. Zooming re-runs only selectMergeEvents below.
    function deriveMergeSeries(nodeCount, edges, metric) {
        const eventRows = mergeEventRows(nodeCount, edges, metric);
        return {
            edges,
            eventRows,
            // From the UNFILTERED rows, so the chart's axis and the slider's track span
            // every event whatever the cap and the window are.
            movingSum: mergeEventMovingSum(eventRows),
            total: eventRows.length,
            // Filled in by selectMergeEvents, which needs a window this does not have.
            selectedRows: [],
            sliderStops: [],
            cap: 0,
        };
    }

    // The window-dependent half: which of those events the chart plots and the slider
    // offers stops for. Cheap (a sort and a scan), so it re-runs on every zoom and pan.
    function selectMergeEvents(series, maxMergeEvents, window, pinnedThreshold) {
        series.selectedRows = filterMergeEventRows(
            series.eventRows, maxMergeEvents, MERGE_EVENT_DENSITY_BINS, window, pinnedThreshold);
        series.sliderStops = buildSliderStops(series.selectedRows, series.edges);
        series.cap = maxMergeEvents;
        return series;
    }

    // The full threshold range of the series, which is what the chart's axis spans and
    // what a zoom window is clamped against. Taken from the UNFILTERED rows (and the
    // moving sum, which runs over them too) rather than from the current selection --
    // the selection depends on the window, so deriving the range from it would be
    // circular, and would make the axis breathe as the user zoomed.
    function splitSeriesDataRange(series) {
        let min = Infinity;
        let max = -Infinity;
        const consider = value => {
            if (!Number.isFinite(value)) { return; }
            if (value < min) { min = value; }
            if (value > max) { max = value; }
        };
        for (const row of series.eventRows) { consider(Number(row.threshold_value)); }
        for (const value of (series.movingSum.x || [])) { consider(Number(value)); }
        return Number.isFinite(min) ? {min, max} : {min: 0, max: 0};
    }

    // Port of ssn_hierarchy.floor_threshold_value.
    const FLOOR_THRESHOLD_SPAN_FRACTION = 0.01;

    function floorThresholdValue(lowestMstWeight, highestMstWeight) {
        const span = highestMstWeight - lowestMstWeight;
        if (span > 0) {
            return lowestMstWeight - (FLOOR_THRESHOLD_SPAN_FRACTION * span);
        }
        const step = Math.abs(lowestMstWeight) * FLOOR_THRESHOLD_SPAN_FRACTION;
        return lowestMstWeight - (step > 0 ? step : 1);
    }

    // The session state, with every saved node id that did not survive dropped, so
    // the extraction opens without complaining about ids it does not contain.
    function extractionAppState(keptIds) {
        const appState = collectSessionState();
        const keep = ids => (ids || []).filter(id => keptIds.has(id));
        appState.selection.node_ids = keep(appState.selection.node_ids);
        const presets = appState.selection.presets || {};
        Object.keys(presets).forEach(slot => {
            const surviving = keep(presets[slot].node_ids);
            if (surviving.length === 0) {
                delete presets[slot];
            } else {
                presets[slot].node_ids = surviving;
            }
        });
        return appState;
    }

    function buildExtractionBundle() {
        const induced = extractionInducedGraph();
        if (induced.nodeIndices.length === 0) {
            return {error: 'Select some nodes or clusters first.'};
        }
        if (induced.brokenComponents > 0) {
            return {error: 'The selection splits ' + induced.brokenComponents.toLocaleString() +
                ' network component' + (induced.brokenComponents === 1 ? '' : 's') + ' into ' +
                induced.brokenPieces.toLocaleString() + ' pieces by leaving out the nodes that ' +
                'join them in the MST. The bundle records no similarity between those pieces, so ' +
                'extracting them would assert they are unrelated. Select the nodes that link them ' +
                '(raising the threshold shows where the links are), or subset the original matrix ' +
                'with build_ssn_viewer.py --subset.'};
        }

        const nodeCount = induced.nodeIndices.length;
        const metric = state.bundle.graph.merge_impact_metric;
        const keptIds = new Set(induced.nodeIndices.map(nodeIndex => nodeId(nodeIndex)));

        return {
            bundle: {
                format: BUNDLE_FORMAT,
                version: BUNDLE_VERSION,
                name: ((state.bundle.name || 'network') + '_extraction'),
                domainator_version: DOMAINATOR_VERSION,
                // A v6 graph: only what cannot be derived from the merge order. The
                // split-event series, its moving sum and the slider's stops are all a
                // function of these four keys, so an extraction no longer computes a
                // frozen copy of them -- whoever opens it derives them, at whatever cap
                // they are looking at it with.
                graph: {
                    nodes: induced.nodeIndices.map(nodeIndex => nodeId(nodeIndex)),
                    mst_edges: induced.edges,
                    merge_impact_metric: metric,
                    hierarchy: buildExtractionHierarchy(nodeCount, induced.edges),
                },
                metadata: {
                    columns: state.metadataColumns.map(column => ({...column})),
                    rows: induced.nodeIndices.map(nodeIndex => state.metadataByNodeIndex[nodeIndex].slice()),
                },
                defaults: {
                    color_by: currentColorField() || null,
                    label_by: currentLabelField() === '__node_id__' ? null : (currentLabelField() || null),
                    categorical_columns: Array.from(state.categoricalColumns),
                },
                app_state: extractionAppState(keptIds),
            },
            nodeCount,
        };
    }

    function saveExtractionFile() {
        if (!state.bundle) { return; }
        // Validate before asking for a name, so a selection that cannot be extracted
        // is refused without making the user name it first.
        const built = buildExtractionBundle();
        if (built.error) {
            setStatus(built.error);
            window.alert(built.error);
            return;
        }
        openNameDialog({
            title: 'Save extraction',
            label: 'Name for the extracted network',
            note: 'Names the new bundle and the file it is saved to. The extraction contains '
                + built.nodeCount.toLocaleString() + ' of '
                + state.bundle.graph.nodes.length.toLocaleString() + ' nodes.',
            value: built.bundle.name,
            confirmLabel: 'Save',
            onConfirm: name => {
                built.bundle.name = name;
                downloadBundleFile(
                    built.bundle,
                    name.replace(/\s+/g, '_'),
                    'Saved "' + name + '": an extraction of ' + built.nodeCount.toLocaleString() + ' nodes'
                );
            },
        });
    }
"""


def _split_chart_hover_js() -> str:
    """Hover readout and click-to-jump for the split chart.

    Plain (non-f) string like `_layout_worker_js`; see `_session_state_js`.
    The readout carries the same quantities matrix_report's Plotly hovertemplates
    do, so the canvas chart and the Plotly one tell the same story.
    """
    return r"""
    // ---------------------------------------------------------------------
    // Split chart hover + click
    //
    // The chart is a canvas, so there are no DOM nodes to attach handlers to and
    // no hit-testing for free. Rather than recomputing where each mark landed --
    // the divergence buildClusterViewSVG warns about -- drawSplitChart() records
    // its own geometry as it paints, and everything here reads that record. A
    // mark the hover can find is therefore a mark that was actually drawn.
    //
    // Clicking jumps the threshold to the event under the cursor, or to the stop
    // nearest the cursor when the pointer is out on the moving-sum line. That is
    // the same action the arrow buttons take, so it keeps the current pan/zoom.
    // ---------------------------------------------------------------------

    // Canvas pixels. The event radius is generous because a stem is 1.5px wide and
    // the beads are 4px: asking for pixel accuracy on a 1100px-wide chart holding
    // hundreds of events would make the readout unusable.
    const SPLIT_CHART_HOVER_RADIUS = 14;
    const SPLIT_CHART_BEAD_RADIUS = 10;

    function recordSplitChartGeometry(geometry) {
        state.splitChartHit = geometry;
        // The only two on-screen paint paths call this, so the zoom controls are
        // refreshed exactly when the window they describe can have changed.
        updateSplitChartZoomControls();
    }

    // ---------------------------------------------------------------------
    // Split chart zoom + pan
    //
    // The window is threshold state, not pixels: state.splitChartZoom is {min, max}
    // in threshold units, or null for "the whole series". Pixels would be wrong the
    // moment the chart is exported at another size, and a scale factor would drift
    // as the range it was taken against changed under it.
    //
    // Nothing here computes geometry of its own either. Every gesture reads the
    // window the paint pass recorded and hands a new one back through
    // setSplitChartWindow, so a gesture cannot act on a chart other than the one
    // the user is looking at.
    // ---------------------------------------------------------------------

    // 500x separates two adjacent lollipops even on a 160,000-event axis, and stops
    // the window shrinking to a width where the arithmetic stops resolving.
    const SPLIT_CHART_MAX_ZOOM = 500;
    // Canvas pixels of travel before a press counts as a pan rather than a click,
    // like the cluster canvas's own drag slop.
    const SPLIT_CHART_DRAG_SLOP = 3;

    // Read by splitChartLayout on every paint. Clamping on read rather than on write
    // is what keeps a stale window harmless: load another bundle, or restore a session
    // saved against a different --max_merge_events, and the window is re-fitted to the
    // range that actually exists instead of framing empty axis.
    function splitChartVisibleWindow(dataMin, dataMax) {
        const dataSpan = dataMax - dataMin;
        const zoom = state.splitChartZoom;
        if (!zoom || !(dataSpan > 0)) {
            return {min: dataMin, max: dataMax, zoomed: false};
        }
        const span = Math.min(dataSpan, Math.max(dataSpan / SPLIT_CHART_MAX_ZOOM, zoom.max - zoom.min));
        const min = Math.min(Math.max(zoom.min, dataMin), dataMax - span);
        return {min, max: min + span, zoomed: span < dataSpan};
    }

    // The one writer of the window, so "zoomed all the way out is not a window at
    // all" is decided once -- and the Reset button, which reads state.splitChartZoom,
    // grays itself out the moment a zoom-out gesture reaches the full range.
    function setSplitChartWindow(min, max, dataMin, dataMax) {
        const dataSpan = dataMax - dataMin;
        if (!(dataSpan > 0) || (max - min) >= dataSpan) {
            state.splitChartZoom = null;
        } else {
            const span = Math.max(dataSpan / SPLIT_CHART_MAX_ZOOM, max - min);
            const clampedMin = Math.min(Math.max(min, dataMin), dataMax - span);
            state.splitChartZoom = {min: clampedMin, max: clampedMin + span};
        }
        // The slider's positions are warped by this window, so it moves with it: zooming
        // the chart into a band magnifies that band on the track below. See
        // positionSliderStops.
        repositionSliderStops();
        drawSplitChart();
    }

    function resetSplitChartZoom() {
        if (!state.splitChartZoom) { return; }
        state.splitChartZoom = null;
        hideSplitChartTip();
        repositionSliderStops();
        drawSplitChart();
    }

    function updateSplitChartZoomControls() {
        const hit = state.splitChartHit;
        const zoomed = Boolean(state.splitChartZoom) && Boolean(hit);
        const button = document.getElementById('split-chart-reset-zoom');
        if (button) { button.disabled = !zoomed; }
        const hint = document.getElementById('split-chart-zoom-hint');
        if (!hint) { return; }
        // The hint doubles as the window's readout: before a gesture it teaches the
        // gestures, and after one it answers the question a zoomed axis raises, which
        // is what part of the range is still on screen.
        hint.textContent = zoomed
            ? 'Showing ' + splitChartWindowLabel(hit.minThreshold, hit.thresholdSpan) + ' – ' +
                splitChartWindowLabel(hit.minThreshold + hit.thresholdSpan, hit.thresholdSpan) +
                ' · double-click to reset'
            : 'Scroll to zoom · drag to pan';
    }

    // Roughly two significant digits of the window's own width, so the two ends read
    // as different numbers however far in the zoom has gone. formatValue's fixed two
    // decimals would print a tight window as "0.87 – 0.87", which is the same lie the
    // axis ticks used to tell before decimalsForTickStep.
    function splitChartWindowLabel(value, span) {
        const decimals = span > 0
            ? Math.min(8, Math.max(2, Math.ceil(-Math.log10(span / 100))))
            : 2;
        return value.toLocaleString(undefined, {
            minimumFractionDigits: decimals,
            maximumFractionDigits: decimals,
        });
    }

    // Zooms about `canvasX`, so the threshold under the pointer stays under the
    // pointer: the gesture reads as moving a lens over the axis rather than as
    // re-centering it somewhere the user did not ask for.
    function zoomSplitChartAt(canvasX, factor) {
        const hit = state.splitChartHit;
        if (!hit || !Number.isFinite(hit.dataMin)) { return; }
        const anchor = splitChartThresholdAt(canvasX);
        if (anchor === null || !(factor > 0)) { return; }
        const span = hit.thresholdSpan / factor;
        const fraction = (anchor - hit.minThreshold) / hit.thresholdSpan;
        setSplitChartWindow(anchor - (fraction * span), anchor + ((1 - fraction) * span),
            hit.dataMin, hit.dataMax);
    }

    // `dx` is pointer travel in canvas pixels: the window moves against it, so the
    // marks follow the pointer.
    function panSplitChartByPixels(dx) {
        const hit = state.splitChartHit;
        if (!hit || hit.plotWidth <= 0 || !Number.isFinite(hit.dataMin)) { return; }
        // Nothing to pan when the whole series is already showing.
        if (!state.splitChartZoom) { return; }
        const delta = -(dx / hit.plotWidth) * hit.thresholdSpan;
        setSplitChartWindow(hit.minThreshold + delta, hit.minThreshold + hit.thresholdSpan + delta,
            hit.dataMin, hit.dataMax);
    }

    function handleSplitChartWheel(event) {
        if (!state.bundle || !state.splitChartHit) { return; }
        event.preventDefault();
        hideSplitChartTip();
        // Shift-wheel, and a trackpad's horizontal wheel, scroll the window; a plain
        // wheel zooms it, with the cluster canvas's own rate so the two feel alike.
        const horizontal = Math.abs(event.deltaX) > Math.abs(event.deltaY);
        if (event.shiftKey || horizontal) {
            panSplitChartByPixels(-(horizontal ? event.deltaX : event.deltaY));
            return;
        }
        const point = splitCanvasCoordinatesFromEvent(event);
        zoomSplitChartAt(point.x, Math.exp(-event.deltaY * 0.0012));
    }

    function handleSplitChartPointerDown(event) {
        if (!state.bundle || event.button !== 0 || !state.splitChartHit) { return; }
        state.splitChartDrag = {
            pointerId: event.pointerId,
            startX: event.clientX,
            lastX: event.clientX,
            moved: false,
        };
        // Captured so a fast drag that leaves the canvas keeps panning, and so the
        // matching pointerup arrives here however far the pointer has traveled.
        if (splitCanvas.setPointerCapture) { splitCanvas.setPointerCapture(event.pointerId); }
    }

    function endSplitChartDrag(event) {
        const drag = state.splitChartDrag;
        if (!drag || drag.pointerId !== event.pointerId) { return; }
        state.splitChartDrag = null;
        splitCanvas.style.cursor = '';
        if (splitCanvas.releasePointerCapture && splitCanvas.hasPointerCapture(event.pointerId)) {
            splitCanvas.releasePointerCapture(event.pointerId);
        }
        // A drag is not a click. Without this, letting go after panning would also
        // jump the threshold to wherever the pointer happened to land.
        state.splitChartSuppressClick = drag.moved;
    }

    // A click sets the threshold and a double-click resets the zoom, so the two clicks
    // inside a double-click have already moved the threshold by the time this fires --
    // a click cannot know a second one is coming. Making every click wait out the
    // double-click interval would make the primary gesture feel broken, so the jump is
    // undone here instead: what the user asked for was the reset, and nothing else.
    function handleSplitChartDoubleClick() {
        if (state.stopBeforeSplitChartClick) {
            snapSliderToStop(state.stopBeforeSplitChartClick);
            scheduleThresholdUI(false);
        }
        resetSplitChartZoom();
    }

    function handleSplitChartPointerCancel(event) {
        endSplitChartDrag(event);
        state.splitChartSuppressClick = false;
        hideSplitChartTip();
    }

    function splitChartTipElement() {
        return document.getElementById('split-chart-tip');
    }

    function hideSplitChartTip() {
        const tip = splitChartTipElement();
        if (tip) { tip.hidden = true; }
        splitCanvas.style.cursor = '';
    }

    function splitCanvasCoordinatesFromEvent(event) {
        // The canvas is laid out at width:100% over a fixed backing store, so client
        // coordinates have to be scaled into canvas space.
        const rect = splitCanvas.getBoundingClientRect();
        return {
            x: (event.clientX - rect.left) * (splitCanvas.width / rect.width),
            y: (event.clientY - rect.top) * (splitCanvas.height / rect.height),
        };
    }

    function splitChartThresholdAt(x) {
        const hit = state.splitChartHit;
        if (!hit || hit.plotWidth <= 0) { return null; }
        const fraction = (x - hit.plotLeft) / hit.plotWidth;
        return hit.minThreshold + (fraction * hit.thresholdSpan);
    }

    // The trace is drawn as a step function, holding each sample's value until the
    // next sample's x. So the value at `threshold` is the last sample at or before it.
    function splitChartMovingSumAt(threshold) {
        const hit = state.splitChartHit;
        if (!hit || !hit.movingSumX || hit.movingSumX.length === 0 || threshold === null) {
            return null;
        }
        const xs = hit.movingSumX;
        if (threshold < xs[0]) { return null; }
        let low = 0;
        let high = xs.length - 1;
        while (low < high) {
            const middle = Math.ceil((low + high) / 2);
            if (xs[middle] <= threshold) { low = middle; } else { high = middle - 1; }
        }
        return hit.movingSumY[low];
    }

    // {kind: 'event', entry, bead} for a split event, {kind: 'movingSum', threshold,
    // value} out on the line, or null when the pointer is outside the plot area.
    function splitChartHitAt(x, y) {
        const hit = state.splitChartHit;
        if (!hit) { return null; }
        const slack = 4;
        if (x < hit.plotLeft - slack || x > hit.plotLeft + hit.plotWidth + slack) { return null; }
        if (y < hit.plotTop - slack || y > hit.plotTop + hit.plotHeight + slack) { return null; }

        let nearest = null;
        let nearestDistance = Infinity;
        for (const entry of hit.events) {
            const distance = Math.abs(entry.x - x);
            if (distance < nearestDistance) {
                nearestDistance = distance;
                nearest = entry;
            }
        }
        if (nearest && nearestDistance <= SPLIT_CHART_HOVER_RADIUS) {
            let bead = null;
            let beadDistance = Infinity;
            for (const candidate of nearest.beads) {
                const distance = Math.abs(candidate.y - y);
                if (distance < beadDistance) {
                    beadDistance = distance;
                    bead = candidate;
                }
            }
            return {
                kind: 'event',
                entry: nearest,
                bead: beadDistance <= SPLIT_CHART_BEAD_RADIUS ? bead : null,
            };
        }
        const threshold = splitChartThresholdAt(x);
        const value = splitChartMovingSumAt(threshold);
        if (threshold === null || value === null) { return null; }
        return {kind: 'movingSum', threshold, value};
    }

    // A min_child impact is a count of nodes; a product impact is not, and must not
    // be labeled as one. The phrasing comes from merge_impact_axis_labels, which
    // matrix_report's hovertemplates read too, so the two charts word it alike.
    function splitImpactAmount(value) {
        const amount = Number(value);
        const phrase = (splitAxisLabels().impactAmount || '{} nodes')
            .replace('{}', amount.toLocaleString());
        // "1 nodes" reads badly, and only the trailing-noun form can be singularized.
        return amount === 1 && phrase.endsWith(' nodes') ? phrase.slice(0, -1) : phrase;
    }

    function splitChartMovingSumLine(threshold) {
        const hit = state.splitChartHit;
        const value = splitChartMovingSumAt(threshold);
        if (value === null) { return null; }
        const halfWindow = (hit.movingSumWindow || 0) / 2;
        const label = splitAxisLabels().movingSumHover || 'Moving sum within';
        return label + ' \u00b1' + formatValue(halfWindow) + ': ' + Number(value).toLocaleString();
    }

    function splitChartTipLines(hit) {
        if (hit.kind === 'movingSum') {
            const lines = ['Threshold ' + formatValue(hit.threshold)];
            const sum = splitChartMovingSumLine(hit.threshold);
            if (sum) { lines.push(sum); }
            lines.push('Click to jump to the nearest split');
            return lines;
        }

        const event = hit.entry.event;
        const lines = [];
        // threshold_from is the cut just above this one, so the pair reads as the
        // interval the merge happened in -- as matrix_report's hover does.
        lines.push(event.threshold_from && event.threshold_from !== event.threshold_to
            ? 'Threshold ' + event.threshold_to + ' (from ' + event.threshold_from + ')'
            : 'Threshold ' + formatValue(event.threshold_value));
        if (hit.bead) {
            lines.push(hit.bead.count.toLocaleString() + ' split' +
                (hit.bead.count === 1 ? '' : 's') + ' of ' + splitImpactAmount(hit.bead.size));
        } else {
            lines.push('Largest single split: ' + splitImpactAmount(event.largest_merge));
        }
        lines.push(event.merge_count.toLocaleString() + ' split' +
            (event.merge_count === 1 ? '' : 's') + ' here, ' +
            splitImpactAmount(event.merge_impact) + ' total');
        const sum = splitChartMovingSumLine(event.threshold_value);
        if (sum) { lines.push(sum); }
        lines.push('Click to jump here');
        return lines;
    }

    function showSplitChartTip(hit, event) {
        const tip = splitChartTipElement();
        if (!tip) { return; }
        tip.textContent = '';
        const lines = splitChartTipLines(hit);
        lines.forEach((text, index) => {
            const row = document.createElement('div');
            // The last line is the affordance, not data, so it is set apart.
            row.className = index === lines.length - 1 ? 'split-tip-hint' : 'split-tip-row';
            row.textContent = text;
            tip.appendChild(row);
        });
        tip.hidden = false;

        // Fixed positioning, flipped away from the viewport edges, like the column
        // chart menu. pointer-events:none keeps it from ever stealing the click.
        const offset = 14;
        const rect = tip.getBoundingClientRect();
        let left = event.clientX + offset;
        let top = event.clientY + offset;
        if (left + rect.width > window.innerWidth - 8) {
            left = Math.max(8, event.clientX - offset - rect.width);
        }
        if (top + rect.height > window.innerHeight - 8) {
            top = Math.max(8, event.clientY - offset - rect.height);
        }
        tip.style.left = left + 'px';
        tip.style.top = top + 'px';
    }

    function handleSplitChartPointerMove(event) {
        if (!state.bundle) { hideSplitChartTip(); return; }
        const drag = state.splitChartDrag;
        if (drag && drag.pointerId === event.pointerId) {
            // Client pixels scaled into canvas space, the same conversion
            // splitCanvasCoordinatesFromEvent makes, because the canvas is laid out at
            // width:100% over a fixed backing store.
            const rect = splitCanvas.getBoundingClientRect();
            const dx = (event.clientX - drag.lastX) * (splitCanvas.width / rect.width);
            drag.lastX = event.clientX;
            if (Math.abs(event.clientX - drag.startX) > SPLIT_CHART_DRAG_SLOP) { drag.moved = true; }
            if (drag.moved) {
                hideSplitChartTip();
                splitCanvas.style.cursor = 'grabbing';
                panSplitChartByPixels(dx);
            }
            return;
        }
        const point = splitCanvasCoordinatesFromEvent(event);
        const hit = splitChartHitAt(point.x, point.y);
        if (!hit) { hideSplitChartTip(); return; }
        splitCanvas.style.cursor = 'pointer';
        showSplitChartTip(hit, event);
    }

    function handleSplitChartClick(event) {
        if (!state.bundle) { return; }
        if (state.splitChartSuppressClick) {
            state.splitChartSuppressClick = false;
            return;
        }
        // Recorded before the jump, and only on the opening click of a gesture, so a
        // double-click can undo the jumps its own clicks made. See
        // handleSplitChartDoubleClick.
        if (event.detail <= 1) {
            state.stopBeforeSplitChartClick = currentSliderStop();
        }
        const point = splitCanvasCoordinatesFromEvent(event);
        const hit = splitChartHitAt(point.x, point.y);
        if (!hit) { return; }
        const threshold = hit.kind === 'event' ? hit.entry.event.threshold_value : hit.threshold;
        const stop = nearestStopForThreshold(threshold);
        if (!stop) { return; }
        snapSliderToStop(stop);
        // Keep the user's pan/zoom, as the threshold arrows do: clicking along the
        // chart to watch one region break up is the point.
        scheduleThresholdUI(false);
    }

    function setupSplitChartHover() {
        splitCanvas.addEventListener('pointermove', handleSplitChartPointerMove);
        splitCanvas.addEventListener('pointerleave', hideSplitChartTip);
        splitCanvas.addEventListener('pointercancel', handleSplitChartPointerCancel);
        splitCanvas.addEventListener('click', handleSplitChartClick);
        // passive:false because the wheel gesture is the zoom, so the page must not
        // scroll out from under it.
        splitCanvas.addEventListener('wheel', handleSplitChartWheel, {passive: false});
        splitCanvas.addEventListener('pointerdown', handleSplitChartPointerDown);
        splitCanvas.addEventListener('pointerup', endSplitChartDrag);
        splitCanvas.addEventListener('dblclick', handleSplitChartDoubleClick);
        document.getElementById('split-chart-reset-zoom')
            .addEventListener('click', resetSplitChartZoom);
        // The readout is positioned from client coordinates, so it is stale the
        // moment the page moves under it.
        window.addEventListener('scroll', hideSplitChartTip, true);
        window.addEventListener('resize', hideSplitChartTip);
    }
"""


def _column_charts_js() -> str:
    """Per-column charts and frequency tables from the metadata table's headers.

    Plain (non-f) string like `_layout_worker_js`; see `_session_state_js`.
    """
    return r"""
    // ---------------------------------------------------------------------
    // Per-column charts and tables
    //
    // Each metadata column header carries a chart glyph that opens a menu of
    // chart kinds appropriate to that column, and picking one opens a dialog.
    //
    // Two invariants hold the feature together:
    //
    // 1. Every chart summarizes exactly the rows the header's own Copy button
    //    copies (columnChartNodeIndices), so a chart and a copied column can
    //    never disagree about what "these nodes" means.
    // 2. Every chart is one SVG string, consumed by three sinks: innerHTML for
    //    the preview, a Blob for the SVG export, and Image->canvas->toBlob for
    //    the PNG export. A canvas preview would need a second renderer for the
    //    SVG export, which is exactly the drift buildClusterViewSVG warns
    //    about. Charts are tens to hundreds of marks, so SVG DOM cost is
    //    irrelevant here (unlike renderClusterView's thousands of nodes) and
    //    axis labels come out as real selectable text.
    //
    // Colors are inherited from the column's own palette via customPalette(),
    // which routes through paletteKey() -- so a numeric column flipped to
    // discrete coloring gets its categorical palette, and a column that has
    // never been the color-by column still gets the same colors it would get if
    // it were (DEFAULT_CATEGORICAL_PALETTE, when nothing has been customized).
    // ---------------------------------------------------------------------

    const COLUMN_CHART_PIE_MAX_SLICES = 12;
    const COLUMN_CHART_BAR_MAX_BARS = 24;
    // The frequency table is a table, so it can afford many more rows than a
    // chart can afford marks -- but not unbounded ones, since it goes through
    // innerHTML. The TSV export always carries every value.
    const COLUMN_CHART_TABLE_MAX_ROWS = 500;
    // An ECDF over a large selection would emit one step per value. Sampling
    // caps the path at a size an SVG editor can still open; the curve is
    // visually identical because adjacent steps land on the same pixel.
    const COLUMN_CHART_ECDF_MAX_POINTS = 2000;

    // Single source of truth for the menu, the dialog's kind <select>, and for
    // validating a kind carried over when the dialog reopens on another column.
    const COLUMN_CHART_KINDS = [
        {id: 'frequency', label: 'Frequency table', numericOnly: false},
        {id: 'bar', label: 'Bar chart', numericOnly: false},
        {id: 'pie', label: 'Pie chart', numericOnly: false},
        {id: 'histogram', label: 'Histogram', numericOnly: true},
        {id: 'box', label: 'Box plot', numericOnly: true},
        {id: 'ecdf', label: 'Cumulative distribution', numericOnly: true},
        {id: 'summary', label: 'Summary statistics', numericOnly: false},
    ];

    // Histogram/box/ECDF are offered for any int/float column even while it is
    // colored as discrete categories: that toggle is a coloring choice, not a
    // claim that the numbers are not numbers.
    function columnChartKindsFor(columnName) {
        const numeric = columnIsNumericType(columnName);
        return COLUMN_CHART_KINDS.filter(kind => numeric || !kind.numericOnly);
    }

    function columnChartKindLabel(kindId) {
        const kind = COLUMN_CHART_KINDS.find(entry => entry.id === kindId);
        return kind ? kind.label : 'Chart';
    }

    // ---- Scope ----

    // The rows every chart summarizes: the table's current rows after the
    // search filter, across all pages. Identical to what copyMetadataColumn
    // exports, and metadataBaseNodeIndices already encodes "the selection, or
    // every node when nothing is selected". Pagination is a rendering artifact,
    // so charting only the visible page would be a trap.
    function columnChartNodeIndices() {
        return metadataDisplayNodeIndices(metadataBaseNodeIndices());
    }

    // Baked into every chart's subtitle so an exported SVG/PNG says what it
    // counted, not just how many.
    function columnChartScopeLabel(rowCount) {
        const total = state.bundle ? state.bundle.graph.nodes.length : 0;
        const parts = [rowCount.toLocaleString() + ' of ' + total.toLocaleString() + ' nodes'];
        parts.push(state.selectedNodeIndices.size > 0 ? 'selection' : 'all nodes');
        const filterText = metadataFilterText();
        if (filterText) {
            parts.push('filter "' + filterText + '"');
        }
        return parts.join(' · ');
    }

    // ---- Data model ----

    // Counts for the same column over every node in the bundle, not just the
    // scoped rows. A chart's own percentages read one way -- "a third of this
    // selection is alpha" -- and a selection immediately raises the other -- "and
    // that is a twelfth of every alpha in the network". Only these counts can
    // answer the second question, because the scope cannot see what it excluded.
    //
    // Recomputed per chart rather than cached: metadata is editable, and one pass
    // over a column costs what the scope's own pass beside it already costs.
    function columnGlobalCounts(columnName) {
        const counts = {byKey: new Map(), nullCount: 0, total: 0};
        const columnIndex = state.metadataColumnIndexByName.get(columnName);
        if (columnIndex === undefined || !state.bundle) {
            return counts;
        }
        counts.total = state.bundle.graph.nodes.length;
        for (let nodeIndex = 0; nodeIndex < counts.total; nodeIndex++) {
            const raw = state.metadataByNodeIndex[nodeIndex]?.[columnIndex] ?? null;
            if (isMissingMetadataValue(raw)) {
                counts.nullCount += 1;
                continue;
            }
            const key = String(raw);
            counts.byKey.set(key, (counts.byKey.get(key) || 0) + 1);
        }
        return counts;
    }

    // Scoped analog of distinctColumnValues (which walks every node): same
    // String(raw) keying, same formatValue labels, and the same count-desc then
    // key-asc ordering, so a chart's category order matches the color picker's
    // swatch order and the legend's rows.
    //
    // options.topN caps the entry count and rolls the tail into one "Other"
    // entry (the pie and bar charts cap; the frequency table lists everything).
    // options.includeNull appends a trailing no-value entry, which every chart
    // wants and only the summary table, which counts nulls in its own row,
    // turns off.
    function columnValueDistribution(columnName, nodeIndices, options = {}) {
        const topN = options.topN === undefined ? Infinity : options.topN;
        const includeNull = options.includeNull !== false;
        const palette = customPalette(columnName);
        const model = {
            entries: [],
            nullCount: 0,
            otherCount: 0,
            otherDistinct: 0,
            distinctTotal: 0,
            rowCount: nodeIndices.length,
            truncated: false,
            globalTotal: 0,
        };
        const columnIndex = state.metadataColumnIndexByName.get(columnName);
        if (columnIndex === undefined) {
            return model;
        }
        const byKey = new Map();
        nodeIndices.forEach(nodeIndex => {
            const raw = state.metadataByNodeIndex[nodeIndex]?.[columnIndex] ?? null;
            if (isMissingMetadataValue(raw)) {
                model.nullCount += 1;
                return;
            }
            const key = String(raw);
            const existing = byKey.get(key);
            if (existing) {
                existing.count += 1;
            } else {
                byKey.set(key, {key, label: formatValue(raw), raw, count: 1});
            }
        });
        model.distinctTotal = byKey.size;
        const globals = columnGlobalCounts(columnName);
        model.globalTotal = globals.total;
        const sorted = Array.from(byKey.values()).sort(
            (a, b) => (b.count - a.count) || a.key.localeCompare(b.key)
        );
        const kept = Number.isFinite(topN) ? sorted.slice(0, topN) : sorted;
        const tail = Number.isFinite(topN) ? sorted.slice(topN) : [];
        model.entries = kept.map(entry => Object.assign({}, entry, {
            color: categoricalColor(entry.raw, palette),
            globalCount: globals.byKey.get(entry.key) || 0,
            isNull: false,
            isOther: false,
        }));
        if (tail.length > 0) {
            model.otherDistinct = tail.length;
            model.otherCount = tail.reduce((sum, entry) => sum + entry.count, 0);
            model.truncated = true;
            // The two synthetic entries below get NUL-prefixed keys so they can
            // never collide with a real value's String(raw), the same trick
            // paletteKey uses to keep a column's two palettes apart.
            model.entries.push({
                key: '\u0000other',
                label: 'Other (' + tail.length.toLocaleString() + ' values)',
                raw: null,
                count: model.otherCount,
                // The rollup's own global count is the sum over the values rolled
                // up, not every value outside the top N: "Other" is defined by what
                // this scope saw, and a value absent from the scope is not in it.
                globalCount: tail.reduce((sum, entry) => sum + (globals.byKey.get(entry.key) || 0), 0),
                // The same neutral gray utils.get_palette gives values with no
                // color of their own, so a rolled-up slice reads as "not a
                // value" rather than as a category with that hue.
                color: PALETTE_NO_VALUE_COLOR,
                isNull: false,
                isOther: true,
            });
        }
        if (includeNull && model.nullCount > 0) {
            model.entries.push({
                key: '\u0000null',
                label: '—',
                raw: null,
                count: model.nullCount,
                globalCount: globals.nullCount,
                // categoricalColor's own null fallback, so an empty cell is the
                // same color in the chart as its node is on the canvas. Note
                // this is a different gray from the "Other" rollup above.
                color: (palette && palette.nullColor) || '#b3a89d',
                isNull: true,
                isOther: false,
            });
        }
        return model;
    }

    function columnNumericValues(columnName, nodeIndices) {
        const result = {values: [], missing: 0};
        const columnIndex = state.metadataColumnIndexByName.get(columnName);
        if (columnIndex === undefined) {
            return result;
        }
        nodeIndices.forEach(nodeIndex => {
            const value = state.metadataByNodeIndex[nodeIndex]?.[columnIndex] ?? null;
            if (typeof value === 'number' && Number.isFinite(value)) {
                result.values.push(value);
            } else {
                result.missing += 1;
            }
        });
        result.values.sort((left, right) => left - right);
        return result;
    }

    // Same shape as columnHistogram, and the same binning rules, but over the
    // scoped rows with the scoped min/max. The two are deliberately separate:
    // columnHistogram bins the *whole column* because the gradient dialog's
    // slider knobs ride that axis, so narrowing it to a selection would move
    // the knobs out from under the user.
    //
    // colorMin/colorMax carry the whole-column range through anyway, so bars
    // are colored the way the nodes are colored (see buildHistogramSVG).
    function scopedColumnHistogram(columnName, nodeIndices) {
        const info = colorInfo(columnName);
        if (!state.bundle || !info || info.baseType !== 'numeric') {
            return null;
        }
        const numeric = columnNumericValues(columnName, nodeIndices);
        const values = numeric.values;
        if (values.length === 0) {
            return null;
        }
        const column = state.metadataColumnByName.get(columnName);
        const min = values[0];
        const max = values[values.length - 1];
        const MAX_BINS = 48;
        const integer = !!(column && column.type === 'int'
            && Number.isInteger(min) && Number.isInteger(max));
        let binCount;
        let lowEdge;
        let highEdge;
        if (max <= min) {
            binCount = 1;
            lowEdge = min - 0.5;
            highEdge = min + 0.5;
        } else if (integer && (max - min) + 1 <= MAX_BINS) {
            // One bin per integer: evenly dividing a handful of distinct
            // integers leaves empty gaps and doubled-up bars that read as
            // structure when they are an artifact of the binning.
            binCount = (max - min) + 1;
            lowEdge = min;
            highEdge = max;
        } else {
            binCount = Math.min(MAX_BINS, Math.max(8, Math.ceil(Math.sqrt(values.length))));
            lowEdge = min;
            highEdge = max;
        }
        const span = highEdge - lowEdge;
        const counts = new Array(binCount).fill(0);
        values.forEach(value => {
            const slot = Math.floor(((value - lowEdge) / span) * binCount);
            counts[Math.max(0, Math.min(binCount - 1, slot))] += 1;
        });
        return {
            counts, lowEdge, highEdge, binCount, integer,
            total: values.length, missing: numeric.missing, min, max,
            colorMin: info.min, colorMax: info.max,
        };
    }

    // Quantiles by linear interpolation between order statistics (numpy's
    // default, and R's type 7), so the numbers match what a user would get
    // running the exported TSV through pandas. `values` must be sorted.
    function quantileOf(values, fraction) {
        if (values.length === 0) {
            return NaN;
        }
        if (values.length === 1) {
            return values[0];
        }
        const position = (values.length - 1) * fraction;
        const lower = Math.floor(position);
        const upper = Math.ceil(position);
        if (lower === upper) {
            return values[lower];
        }
        return values[lower] + ((values[upper] - values[lower]) * (position - lower));
    }

    // Five-number summary plus mean/SD and Tukey 1.5*IQR whiskers, shared by
    // the box plot and the summary table so the two can never disagree.
    // `values` must be sorted ascending.
    function numericSummary(values) {
        if (values.length === 0) {
            return null;
        }
        const count = values.length;
        const min = values[0];
        const max = values[count - 1];
        const q1 = quantileOf(values, 0.25);
        const median = quantileOf(values, 0.5);
        const q3 = quantileOf(values, 0.75);
        const iqr = q3 - q1;
        const mean = values.reduce((sum, value) => sum + value, 0) / count;
        // Sample standard deviation (n-1), undefined for a single value.
        const variance = count > 1
            ? values.reduce((sum, value) => sum + ((value - mean) ** 2), 0) / (count - 1)
            : 0;
        const stdev = count > 1 ? Math.sqrt(variance) : NaN;
        const lowFence = q1 - (1.5 * iqr);
        const highFence = q3 + (1.5 * iqr);
        // Whiskers reach the most extreme value still inside the fences, never
        // past the data -- so a distribution with no outliers draws whiskers at
        // min and max. Seeded crossed so the first value inside wins; at least
        // one always is, since q1 <= median <= q3 sits inside both fences.
        let lowWhisker = max;
        let highWhisker = min;
        const outliers = [];
        values.forEach(value => {
            if (value < lowFence || value > highFence) {
                outliers.push(value);
                return;
            }
            lowWhisker = Math.min(lowWhisker, value);
            highWhisker = Math.max(highWhisker, value);
        });
        return {
            count, min, q1, median, q3, max, iqr, mean, stdev,
            lowWhisker, highWhisker, outliers,
        };
    }

    // Mirrors formatMetadataDisplayValue's float branch, but is never routed
    // through a column's declared type: the quartiles, mean and SD of an int
    // column are not integers, and the int formatter would truncate a first
    // quartile of 2.25 to "2".
    function formatStatNumber(value) {
        if (!Number.isFinite(value)) {
            return '—';
        }
        if (Number.isInteger(value)) {
            return value.toLocaleString();
        }
        const magnitude = Math.abs(value);
        if (magnitude >= 10000 || magnitude < 0.001) {
            return value.toExponential(3);
        }
        return new Intl.NumberFormat(undefined, {maximumFractionDigits: 4}).format(value);
    }

    // [label, text] rows for the summary table.
    function columnSummaryRows(columnName, nodeIndices) {
        const model = columnValueDistribution(columnName, nodeIndices, {includeNull: false});
        const rows = [
            ['Rows', model.rowCount.toLocaleString()],
            ['With a value', (model.rowCount - model.nullCount).toLocaleString()],
            ['No value', model.nullCount.toLocaleString()],
            ['Distinct values', model.distinctTotal.toLocaleString()],
        ];
        if (model.entries.length > 0) {
            const top = model.entries[0];
            rows.push(['Most common', top.label + ' (' + top.count.toLocaleString() + ')']);
        }
        if (!columnIsNumericType(columnName)) {
            return rows;
        }
        const summary = numericSummary(columnNumericValues(columnName, nodeIndices).values);
        if (!summary) {
            return rows;
        }
        rows.push(['Minimum', formatStatNumber(summary.min)]);
        rows.push(['1st quartile', formatStatNumber(summary.q1)]);
        rows.push(['Median', formatStatNumber(summary.median)]);
        rows.push(['3rd quartile', formatStatNumber(summary.q3)]);
        rows.push(['Maximum', formatStatNumber(summary.max)]);
        rows.push(['Mean', formatStatNumber(summary.mean)]);
        rows.push(['Std. deviation', formatStatNumber(summary.stdev)]);
        return rows;
    }

    // {color, label, count, percent, globalCount, globalPercent} rows shared by
    // the HTML preview, the SVG export and the TSV, so all three show the same
    // numbers. The two percentages have different denominators on purpose:
    // `percent` is the value's share of the charted rows, `globalPercent` is the
    // share of that value's own network-wide population which those rows caught.
    // With nothing selected and no filter the scope is the whole bundle, so every
    // globalPercent is 100% -- which is the honest reading, not a bug.
    function frequencyTableRows(model) {
        const denominator = model.rowCount > 0 ? model.rowCount : 1;
        return model.entries.map(entry => {
            const globalCount = entry.globalCount || 0;
            return {
                color: entry.color,
                label: entry.label,
                count: entry.count,
                percent: (entry.count / denominator) * 100,
                globalCount,
                globalPercent: globalCount > 0 ? (entry.count / globalCount) * 100 : 0,
            };
        });
    }

    // ---- SVG scaffolding ----

    const CHART_TITLE_FONT = 18;
    const CHART_SUBTITLE_FONT = 12;
    const CHART_LABEL_FONT = 13;
    const CHART_PAD = 14;
    // Tall enough to clear the title and subtitle baselines with room to spare:
    // the histogram and ECDF label their topmost gridline, which sits here.
    const CHART_HEAD = 64;

    // The same width approximation buildLegendSVG uses. Measuring text properly
    // would mean laying it out in the DOM first, which the export path (a
    // detached string) cannot do.
    function estimateTextWidth(text, fontSize) {
        return (String(text).length * fontSize) / 2;
    }

    function truncateChartLabel(text, maxChars) {
        const value = String(text);
        return value.length <= maxChars ? value : (value.slice(0, Math.max(1, maxChars - 1)) + '…');
    }

    // Emit one color as an SVG-editor-safe fill attribute pair. svgColorParts
    // normalizes through the canvas parser (the default schemes emit hsl(),
    // which Illustrator and Inkscape choke on) and splits any alpha out into
    // its own attribute.
    function chartSvgFill(color) {
        const parts = svgColorParts(color);
        return 'fill="' + parts.color + '" fill-opacity="' + parts.opacity + '"';
    }

    // Deliberately no <?xml ...?> prolog: the same string is assigned to
    // innerHTML for the preview, which refuses a processing instruction. The
    // export path prepends it (see exportColumnChartSVG).
    function chartSvgOpen(parts, width, height, title, subtitle) {
        parts.push('<svg xmlns="http://www.w3.org/2000/svg" width="' + width + '" height="' + height +
            '" viewBox="0 0 ' + width + ' ' + height + '">');
        parts.push('<g font-family="Arial, sans-serif">');
        parts.push('<rect x="0" y="0" width="' + width + '" height="' + height + '" fill="white"/>');
        parts.push('<text x="' + CHART_PAD + '" y="' + (CHART_PAD + CHART_TITLE_FONT - 4) +
            '" font-size="' + CHART_TITLE_FONT + '" font-weight="bold">' + escapeXml(title) + '</text>');
        parts.push('<text x="' + CHART_PAD + '" y="' + (CHART_PAD + CHART_TITLE_FONT + CHART_SUBTITLE_FONT + 2) +
            '" font-size="' + CHART_SUBTITLE_FONT + '" fill="#5c6a70">' + escapeXml(subtitle) + '</text>');
    }

    function chartSvgClose(parts) {
        parts.push('</g>');
        parts.push('</svg>');
        return parts.join('\n');
    }

    // "Nice" 1/2/2.5/5 tick values spanning [min, max]. Used by the histogram,
    // box plot and ECDF, all of which need a readable numeric axis rather than
    // the legend's fixed end ticks.
    function chartAxisTicks(min, max, targetCount) {
        if (!Number.isFinite(min) || !Number.isFinite(max) || max <= min) {
            return [min];
        }
        const rawStep = (max - min) / Math.max(1, targetCount);
        const magnitude = Math.pow(10, Math.floor(Math.log10(rawStep)));
        const normalized = rawStep / magnitude;
        // 2.5 is on the ladder because a bare 1/2/5/10 rounds a raw step of
        // 2.2 all the way up to 5, leaving a wide axis with two ticks on it.
        const step = magnitude * (normalized <= 1 ? 1
            : (normalized <= 2 ? 2 : (normalized <= 2.5 ? 2.5 : (normalized <= 5 ? 5 : 10))));
        const ticks = [];
        const first = Math.ceil(min / step) * step;
        for (let value = first; value <= max + (step / 1000); value += step) {
            // Re-round to kill the float drift that accumulates over the loop.
            ticks.push(Number((Math.round(value / step) * step).toPrecision(12)));
        }
        if (ticks.length === 0) {
            return [min, max];
        }
        return ticks;
    }

    // SVG arc path for a pie slice. Angles are in radians, measured clockwise
    // from twelve o'clock so slices read in the order the legend lists them.
    function describeArcPath(cx, cy, radius, startAngle, endAngle) {
        const point = angle => {
            const x = cx + (radius * Math.sin(angle));
            const y = cy - (radius * Math.cos(angle));
            return x.toFixed(2) + ' ' + y.toFixed(2);
        };
        const largeArc = (endAngle - startAngle) > Math.PI ? 1 : 0;
        return 'M ' + cx.toFixed(2) + ' ' + cy.toFixed(2) +
            ' L ' + point(startAngle) +
            ' A ' + radius + ' ' + radius + ' 0 ' + largeArc + ' 1 ' + point(endAngle) + ' Z';
    }

    // ---- Chart builders ----
    // Each returns {svg, width, height}.

    function buildBarChartSVG(model, meta) {
        const rows = frequencyTableRows(model);
        const barHeight = 22;
        const barGap = 6;
        const labelChars = 28;
        const labelWidth = Math.min(
            260,
            Math.max(90, rows.reduce(
                (widest, row) => Math.max(widest, estimateTextWidth(truncateChartLabel(row.label, labelChars), CHART_LABEL_FONT)),
                0
            ) + 8)
        );
        const countWidth = 96;
        const plotWidth = 420;
        const width = CHART_PAD + labelWidth + plotWidth + countWidth + CHART_PAD;
        const height = CHART_HEAD + (rows.length * (barHeight + barGap)) + CHART_PAD;
        const parts = [];
        chartSvgOpen(parts, width, height, meta.title, meta.subtitle);
        const maxCount = rows.reduce((most, row) => Math.max(most, row.count), 0) || 1;
        const barX = CHART_PAD + labelWidth;
        rows.forEach((row, index) => {
            const y = CHART_HEAD + (index * (barHeight + barGap));
            parts.push('<text x="' + (barX - 8) + '" y="' + (y + (barHeight / 2)) +
                '" font-size="' + CHART_LABEL_FONT + '" text-anchor="end" dominant-baseline="central">' +
                escapeXml(truncateChartLabel(row.label, labelChars)) + '</text>');
            // A floor of 1px keeps a count of 1 visible next to a dominant bar.
            const barWidth = Math.max(1, (row.count / maxCount) * plotWidth);
            parts.push('<rect x="' + barX + '" y="' + y + '" width="' + barWidth.toFixed(2) +
                '" height="' + barHeight + '" ' + chartSvgFill(row.color) +
                ' stroke="#1e2a2f" stroke-width="0.75"/>');
            parts.push('<text x="' + (barX + barWidth + 8) + '" y="' + (y + (barHeight / 2)) +
                '" font-size="' + CHART_LABEL_FONT + '" fill="#5c6a70" dominant-baseline="central">' +
                escapeXml(row.count.toLocaleString() + ' (' + row.percent.toFixed(1) + '%)') + '</text>');
        });
        return {svg: chartSvgClose(parts), width, height};
    }

    function buildPieChartSVG(model, meta) {
        const rows = frequencyTableRows(model);
        const radius = 130;
        const diameter = radius * 2;
        const swatch = 14;
        const legendRowHeight = 22;
        const labelChars = 30;
        const legendWidth = Math.max(140, rows.reduce(
            (widest, row) => Math.max(widest, estimateTextWidth(
                truncateChartLabel(row.label, labelChars) + '  ' + row.count.toLocaleString() + ' (' + row.percent.toFixed(1) + '%)',
                CHART_LABEL_FONT
            )),
            0
        ) + swatch + 16);
        const width = CHART_PAD + diameter + 24 + legendWidth + CHART_PAD;
        const height = Math.max(
            CHART_HEAD + diameter + CHART_PAD,
            CHART_HEAD + (rows.length * legendRowHeight) + CHART_PAD
        );
        const parts = [];
        chartSvgOpen(parts, width, height, meta.title, meta.subtitle);
        const cx = CHART_PAD + radius;
        const cy = CHART_HEAD + radius;
        const total = rows.reduce((sum, row) => sum + row.count, 0);
        if (total > 0) {
            if (rows.length === 1) {
                // A single 100% slice: an arc whose start and end angles are
                // equal draws nothing, so the whole circle has to be a circle.
                parts.push('<circle cx="' + cx + '" cy="' + cy + '" r="' + radius + '" ' +
                    chartSvgFill(rows[0].color) + ' stroke="#1e2a2f" stroke-width="0.75"/>');
            } else {
                let angle = 0;
                rows.forEach(row => {
                    const sweep = (row.count / total) * Math.PI * 2;
                    parts.push('<path d="' + describeArcPath(cx, cy, radius, angle, angle + sweep) + '" ' +
                        chartSvgFill(row.color) + ' stroke="#1e2a2f" stroke-width="0.75"/>');
                    angle += sweep;
                });
            }
        }
        const legendX = CHART_PAD + diameter + 24;
        rows.forEach((row, index) => {
            const y = CHART_HEAD + (index * legendRowHeight);
            parts.push('<rect x="' + legendX + '" y="' + (y + ((legendRowHeight - swatch) / 2)) +
                '" width="' + swatch + '" height="' + swatch + '" ' + chartSvgFill(row.color) +
                ' stroke="#1e2a2f" stroke-width="0.75"/>');
            parts.push('<text x="' + (legendX + swatch + 8) + '" y="' + (y + (legendRowHeight / 2)) +
                '" font-size="' + CHART_LABEL_FONT + '" dominant-baseline="central">' +
                escapeXml(truncateChartLabel(row.label, labelChars) + '  ' +
                    row.count.toLocaleString() + ' (' + row.percent.toFixed(1) + '%)') + '</text>');
        });
        return {svg: chartSvgClose(parts), width, height};
    }

    // Bars are filled with the color their bin's values receive on the network
    // canvas -- so this doubles as a check on the column's gradient. The color
    // domain is the whole column (colorMin/colorMax), not the scoped range,
    // because that is what nodeColor uses.
    function buildHistogramSVG(histogram, meta) {
        const palette = customPalette(meta.columnName);
        const plotWidth = 620;
        const plotHeight = 240;
        const axisSpace = 44;
        const yLabelWidth = 56;
        const width = CHART_PAD + yLabelWidth + plotWidth + CHART_PAD;
        const height = CHART_HEAD + plotHeight + axisSpace;
        const parts = [];
        chartSvgOpen(parts, width, height, meta.title, meta.subtitle);
        const plotX = CHART_PAD + yLabelWidth;
        const baselineY = CHART_HEAD + plotHeight;
        const maxCount = Math.max(...histogram.counts) || 1;
        const slotWidth = plotWidth / histogram.binCount;
        const gap = slotWidth > 8 ? Math.min(3, slotWidth * 0.16) : 0;
        const binSpan = (histogram.highEdge - histogram.lowEdge) / histogram.binCount;

        // Y axis: zero and the peak count are the two numbers that give the
        // bars their scale.
        [0, maxCount].forEach(count => {
            const y = baselineY - ((count / maxCount) * plotHeight);
            parts.push('<line x1="' + plotX + '" y1="' + y.toFixed(2) + '" x2="' + (plotX + plotWidth) +
                '" y2="' + y.toFixed(2) + '" stroke="#d8dce2" stroke-width="1"/>');
            parts.push('<text x="' + (plotX - 8) + '" y="' + y.toFixed(2) +
                '" font-size="' + CHART_LABEL_FONT + '" fill="#5c6a70" text-anchor="end" dominant-baseline="central">' +
                escapeXml(count.toLocaleString()) + '</text>');
        });

        histogram.counts.forEach((count, binIndex) => {
            if (count === 0) {
                return;
            }
            // A floor of 2px keeps rare bins from vanishing beside a dominant one.
            const barHeight = Math.max(2, (count / maxCount) * plotHeight);
            const center = histogram.lowEdge + ((binIndex + 0.5) * binSpan);
            const color = numericColor(center, histogram.colorMin, histogram.colorMax, palette);
            parts.push('<rect x="' + (plotX + (binIndex * slotWidth) + (gap / 2)).toFixed(2) +
                '" y="' + (baselineY - barHeight).toFixed(2) +
                '" width="' + Math.max(1, slotWidth - gap).toFixed(2) +
                '" height="' + barHeight.toFixed(2) + '" ' + chartSvgFill(color) + '/>');
        });

        parts.push('<line x1="' + plotX + '" y1="' + baselineY + '" x2="' + (plotX + plotWidth) +
            '" y2="' + baselineY + '" stroke="#1e2a2f" stroke-width="1.25"/>');

        // X ticks, dropping any that would collide with the label before it --
        // the same pruning buildLegendSVG does on its gradient ticks.
        const tickFont = 12;
        let lastRight = -Infinity;
        chartAxisTicks(histogram.lowEdge, histogram.highEdge, 6).forEach(value => {
            const fraction = (value - histogram.lowEdge) / (histogram.highEdge - histogram.lowEdge);
            const x = plotX + (fraction * plotWidth);
            const text = formatValue(value);
            const textWidth = estimateTextWidth(text, tickFont);
            if (x - (textWidth / 2) < lastRight + 6) {
                return;
            }
            parts.push('<line x1="' + x.toFixed(2) + '" y1="' + baselineY + '" x2="' + x.toFixed(2) +
                '" y2="' + (baselineY + 5) + '" stroke="#1e2a2f" stroke-width="1"/>');
            parts.push('<text x="' + x.toFixed(2) + '" y="' + (baselineY + 5 + tickFont + 2) +
                '" font-size="' + tickFont + '" text-anchor="middle">' + escapeXml(text) + '</text>');
            lastRight = x + (textWidth / 2);
        });
        parts.push('<text x="' + (plotX + (plotWidth / 2)) + '" y="' + (height - 4) +
            '" font-size="' + tickFont + '" fill="#5c6a70" text-anchor="middle">' +
            escapeXml(meta.columnName) + '</text>');
        return {svg: chartSvgClose(parts), width, height};
    }

    function buildBoxPlotSVG(summary, meta) {
        const palette = customPalette(meta.columnName);
        const info = colorInfo(meta.columnName);
        const colorMin = info ? info.min : summary.min;
        const colorMax = info ? info.max : summary.max;
        const plotWidth = 620;
        const boxHeight = 74;
        const axisSpace = 44;
        const width = CHART_PAD + plotWidth + CHART_PAD;
        const height = CHART_HEAD + boxHeight + axisSpace;
        const parts = [];
        chartSvgOpen(parts, width, height, meta.title, meta.subtitle);
        const plotX = CHART_PAD;
        // Pad the domain so a whisker or outlier at the extreme is not clipped
        // by the plot edge.
        const dataLow = Math.min(summary.lowWhisker, summary.min);
        const dataHigh = Math.max(summary.highWhisker, summary.max);
        const domainSpan = dataHigh - dataLow;
        const pad = domainSpan > 0 ? domainSpan * 0.04 : 0.5;
        const low = dataLow - pad;
        const high = dataHigh + pad;
        const xOf = value => plotX + (((value - low) / (high - low)) * plotWidth);
        const midY = CHART_HEAD + (boxHeight / 2);

        parts.push('<line x1="' + xOf(summary.lowWhisker).toFixed(2) + '" y1="' + midY +
            '" x2="' + xOf(summary.highWhisker).toFixed(2) + '" y2="' + midY +
            '" stroke="#1e2a2f" stroke-width="1.25"/>');
        [summary.lowWhisker, summary.highWhisker].forEach(value => {
            parts.push('<line x1="' + xOf(value).toFixed(2) + '" y1="' + (midY - 14) +
                '" x2="' + xOf(value).toFixed(2) + '" y2="' + (midY + 14) +
                '" stroke="#1e2a2f" stroke-width="1.25"/>');
        });
        const boxLeft = xOf(summary.q1);
        const boxRight = xOf(summary.q3);
        parts.push('<rect x="' + boxLeft.toFixed(2) + '" y="' + (midY - 26) +
            '" width="' + Math.max(1, boxRight - boxLeft).toFixed(2) + '" height="52" ' +
            chartSvgFill(numericColor(summary.median, colorMin, colorMax, palette)) +
            ' stroke="#1e2a2f" stroke-width="1.25"/>');
        parts.push('<line x1="' + xOf(summary.median).toFixed(2) + '" y1="' + (midY - 26) +
            '" x2="' + xOf(summary.median).toFixed(2) + '" y2="' + (midY + 26) +
            '" stroke="#1e2a2f" stroke-width="2.5"/>');
        summary.outliers.forEach(value => {
            parts.push('<circle cx="' + xOf(value).toFixed(2) + '" cy="' + midY + '" r="3.5" ' +
                chartSvgFill(numericColor(value, colorMin, colorMax, palette)) +
                ' stroke="#1e2a2f" stroke-width="0.75"/>');
        });

        const baselineY = CHART_HEAD + boxHeight;
        parts.push('<line x1="' + plotX + '" y1="' + baselineY + '" x2="' + (plotX + plotWidth) +
            '" y2="' + baselineY + '" stroke="#1e2a2f" stroke-width="1.25"/>');
        const tickFont = 12;
        let lastRight = -Infinity;
        chartAxisTicks(low, high, 6).forEach(value => {
            const x = xOf(value);
            const text = formatValue(value);
            const textWidth = estimateTextWidth(text, tickFont);
            if (x - (textWidth / 2) < lastRight + 6) {
                return;
            }
            parts.push('<line x1="' + x.toFixed(2) + '" y1="' + baselineY + '" x2="' + x.toFixed(2) +
                '" y2="' + (baselineY + 5) + '" stroke="#1e2a2f" stroke-width="1"/>');
            parts.push('<text x="' + x.toFixed(2) + '" y="' + (baselineY + 5 + tickFont + 2) +
                '" font-size="' + tickFont + '" text-anchor="middle">' + escapeXml(text) + '</text>');
            lastRight = x + (textWidth / 2);
        });
        parts.push('<text x="' + (plotX + (plotWidth / 2)) + '" y="' + (height - 4) +
            '" font-size="' + tickFont + '" fill="#5c6a70" text-anchor="middle">' +
            escapeXml(meta.columnName) + '</text>');
        return {svg: chartSvgClose(parts), width, height};
    }

    // Empirical cumulative distribution: the fraction of values at or below
    // each x. The right chart for picking a cutoff, since you can read the
    // proportion kept straight off the curve.
    function buildEcdfSVG(values, meta) {
        const plotWidth = 620;
        const plotHeight = 240;
        const axisSpace = 44;
        const yLabelWidth = 56;
        const width = CHART_PAD + yLabelWidth + plotWidth + CHART_PAD;
        const height = CHART_HEAD + plotHeight + axisSpace;
        const parts = [];
        chartSvgOpen(parts, width, height, meta.title, meta.subtitle);
        const plotX = CHART_PAD + yLabelWidth;
        const baselineY = CHART_HEAD + plotHeight;
        const low = values[0];
        const high = values[values.length - 1];
        const xOf = value => (high > low
            ? plotX + (((value - low) / (high - low)) * plotWidth)
            : plotX + (plotWidth / 2));
        const yOf = fraction => baselineY - (fraction * plotHeight);

        [0, 0.25, 0.5, 0.75, 1].forEach(fraction => {
            const y = yOf(fraction);
            parts.push('<line x1="' + plotX + '" y1="' + y.toFixed(2) + '" x2="' + (plotX + plotWidth) +
                '" y2="' + y.toFixed(2) + '" stroke="#d8dce2" stroke-width="1"/>');
            parts.push('<text x="' + (plotX - 8) + '" y="' + y.toFixed(2) +
                '" font-size="' + CHART_LABEL_FONT + '" fill="#5c6a70" text-anchor="end" dominant-baseline="central">' +
                escapeXml((fraction * 100).toFixed(0) + '%') + '</text>');
        });

        // Sample when there are more values than the cap: adjacent steps would
        // land on the same pixel anyway, and the full list would emit a path an
        // SVG editor struggles to open.
        const stride = Math.max(1, Math.ceil(values.length / COLUMN_CHART_ECDF_MAX_POINTS));
        const path = ['M ' + xOf(low).toFixed(2) + ' ' + yOf(0).toFixed(2)];
        for (let index = 0; index < values.length; index += stride) {
            const fraction = (index + 1) / values.length;
            path.push('L ' + xOf(values[index]).toFixed(2) + ' ' + yOf(fraction - (1 / values.length)).toFixed(2));
            path.push('L ' + xOf(values[index]).toFixed(2) + ' ' + yOf(fraction).toFixed(2));
        }
        // Always land on the last value at 100%, whatever the stride skipped.
        path.push('L ' + xOf(high).toFixed(2) + ' ' + yOf(1).toFixed(2));
        parts.push('<path d="' + path.join(' ') + '" fill="none" stroke="#c8553d" stroke-width="2"/>');

        parts.push('<line x1="' + plotX + '" y1="' + baselineY + '" x2="' + (plotX + plotWidth) +
            '" y2="' + baselineY + '" stroke="#1e2a2f" stroke-width="1.25"/>');
        const tickFont = 12;
        let lastRight = -Infinity;
        chartAxisTicks(low, high, 6).forEach(value => {
            const x = xOf(value);
            const text = formatValue(value);
            const textWidth = estimateTextWidth(text, tickFont);
            if (x - (textWidth / 2) < lastRight + 6) {
                return;
            }
            parts.push('<line x1="' + x.toFixed(2) + '" y1="' + baselineY + '" x2="' + x.toFixed(2) +
                '" y2="' + (baselineY + 5) + '" stroke="#1e2a2f" stroke-width="1"/>');
            parts.push('<text x="' + x.toFixed(2) + '" y="' + (baselineY + 5 + tickFont + 2) +
                '" font-size="' + tickFont + '" text-anchor="middle">' + escapeXml(text) + '</text>');
            lastRight = x + (textWidth / 2);
        });
        parts.push('<text x="' + (plotX + (plotWidth / 2)) + '" y="' + (height - 4) +
            '" font-size="' + tickFont + '" fill="#5c6a70" text-anchor="middle">' +
            escapeXml(meta.columnName) + '</text>');
        return {svg: chartSvgClose(parts), width, height};
    }

    // Table renderer shared by the frequency and summary exports. `columns` is
    // [{key, label, align, swatch}]; a swatch column draws each row's color
    // instead of text.
    function buildTableSVG(rows, columns, meta) {
        const rowHeight = 24;
        const cellPad = 10;
        const swatchSize = 14;
        const widths = columns.map(column => {
            if (column.swatch) {
                return swatchSize + (cellPad * 2);
            }
            const widest = rows.reduce(
                (most, row) => Math.max(most, estimateTextWidth(row[column.key], CHART_LABEL_FONT)),
                estimateTextWidth(column.label, CHART_LABEL_FONT)
            );
            return Math.min(320, widest + (cellPad * 2));
        });
        const width = CHART_PAD + widths.reduce((sum, value) => sum + value, 0) + CHART_PAD;
        const height = CHART_HEAD + ((rows.length + 1) * rowHeight) + CHART_PAD;
        const parts = [];
        chartSvgOpen(parts, width, height, meta.title, meta.subtitle);
        const offsets = [];
        let cursor = CHART_PAD;
        widths.forEach(value => {
            offsets.push(cursor);
            cursor += value;
        });
        const cellText = (text, x, y, align, weight) => {
            const anchor = align === 'right' ? 'end' : 'start';
            const textX = align === 'right' ? (x - cellPad) : (x + cellPad);
            parts.push('<text x="' + textX + '" y="' + y + '" font-size="' + CHART_LABEL_FONT +
                '" text-anchor="' + anchor + '" dominant-baseline="central"' +
                (weight ? ' font-weight="' + weight + '"' : '') + '>' + escapeXml(text) + '</text>');
        };
        columns.forEach((column, columnIndex) => {
            if (!column.label) {
                // The swatch column has no heading, so emit no text node at all.
                return;
            }
            const x = column.align === 'right' ? offsets[columnIndex] + widths[columnIndex] : offsets[columnIndex];
            cellText(column.label, x, CHART_HEAD + (rowHeight / 2), column.align, 'bold');
        });
        parts.push('<line x1="' + CHART_PAD + '" y1="' + (CHART_HEAD + rowHeight) +
            '" x2="' + (width - CHART_PAD) + '" y2="' + (CHART_HEAD + rowHeight) +
            '" stroke="#1e2a2f" stroke-width="1.25"/>');
        rows.forEach((row, rowIndex) => {
            const y = CHART_HEAD + ((rowIndex + 1.5) * rowHeight);
            columns.forEach((column, columnIndex) => {
                if (column.swatch) {
                    parts.push('<rect x="' + (offsets[columnIndex] + cellPad) + '" y="' + (y - (swatchSize / 2)) +
                        '" width="' + swatchSize + '" height="' + swatchSize + '" ' + chartSvgFill(row.color) +
                        ' stroke="#1e2a2f" stroke-width="0.75"/>');
                    return;
                }
                const x = column.align === 'right' ? offsets[columnIndex] + widths[columnIndex] : offsets[columnIndex];
                cellText(row[column.key], x, y, column.align, null);
            });
            parts.push('<line x1="' + CHART_PAD + '" y1="' + (y + (rowHeight / 2)) +
                '" x2="' + (width - CHART_PAD) + '" y2="' + (y + (rowHeight / 2)) +
                '" stroke="#d8dce2" stroke-width="0.75"/>');
        });
        return {svg: chartSvgClose(parts), width, height};
    }

    // ---- HTML table previews ----
    // The two table kinds preview as real HTML rather than as their export SVG
    // so the numbers can be selected and copied out of the dialog.

    function frequencyTableHTML(rows) {
        // cssToHex, not just htmlEscape: this color lands inside a style
        // attribute, and a hex literal cannot carry extra CSS declarations with
        // it however the palette it came from was populated.
        const body = rows.map(row =>
            '<tr><td class="cc-swatch-cell"><span class="cc-swatch" style="background:' +
            cssToHex(row.color) + '"></span></td>' +
            '<td>' + htmlEscape(row.label) + '</td>' +
            '<td class="cc-num">' + htmlEscape(row.count.toLocaleString()) + '</td>' +
            '<td class="cc-num">' + htmlEscape(row.percent.toFixed(1) + '%') + '</td>' +
            '<td class="cc-num">' + htmlEscape(row.globalCount.toLocaleString()) + '</td>' +
            '<td class="cc-num">' + htmlEscape(row.globalPercent.toFixed(1) + '%') + '</td></tr>'
        ).join('');
        // "Percent" is read against the charted rows; "% of all" against that
        // value's whole population, whose size is the column beside it so the
        // denominator is never a mystery.
        return '<table class="cc-freq-table"><thead><tr><th class="cc-swatch-cell"></th>' +
            '<th>Value</th><th class="cc-num">Count</th><th class="cc-num">Percent</th>' +
            '<th class="cc-num">All nodes</th><th class="cc-num">% of all</th>' +
            '</tr></thead><tbody>' + body + '</tbody></table>';
    }

    function summaryTableHTML(rows) {
        const body = rows.map(row =>
            '<tr><td>' + htmlEscape(row[0]) + '</td>' +
            '<td class="cc-num">' + htmlEscape(row[1]) + '</td></tr>'
        ).join('');
        return '<table class="cc-freq-table"><thead><tr><th>Statistic</th>' +
            '<th class="cc-num">Value</th></tr></thead><tbody>' + body + '</tbody></table>';
    }

    function chartTSV(header, rows) {
        return [header.join('\t')].concat(rows.map(row => row.join('\t'))).join('\n') + '\n';
    }

    // ---- The one dispatcher ----

    // Returns {kind, columnName, svg, width, height, html, tsv, notes, empty}.
    // `html` is what the preview shows when set, otherwise `svg` is; `svg` is
    // always present so every kind can be exported as an image.
    function buildColumnChartArtifact(columnName, kind, nodeIndices) {
        const subtitle = columnChartScopeLabel(nodeIndices.length);
        const meta = {
            columnName,
            title: columnName + ' — ' + columnChartKindLabel(kind).toLowerCase(),
            subtitle,
        };
        const notes = [subtitle];
        const empty = {kind, columnName, svg: null, html: null, tsv: '', notes, empty: true};

        if (kind === 'histogram' || kind === 'box' || kind === 'ecdf') {
            const numeric = columnNumericValues(columnName, nodeIndices);
            if (numeric.values.length === 0) {
                notes.push('No numeric values in these rows.');
                return empty;
            }
            if (numeric.missing > 0) {
                notes.push(numeric.missing.toLocaleString() + ' with no value (not plotted)');
            }
            if (kind === 'histogram') {
                const histogram = scopedColumnHistogram(columnName, nodeIndices);
                if (!histogram) {
                    return empty;
                }
                notes.push(histogram.total.toLocaleString() +
                    ' value' + (histogram.total === 1 ? '' : 's') +
                    ' in ' + histogram.binCount + ' bin' + (histogram.binCount === 1 ? '' : 's'));
                const built = buildHistogramSVG(histogram, meta);
                // One bin per integer means the bins *are* values, and the
                // edges the binning uses to separate them (min + k*span/bins)
                // are an implementation detail nobody wants in a TSV.
                const perInteger = histogram.integer
                    && histogram.binCount === (histogram.max - histogram.min) + 1;
                const binSpan = (histogram.highEdge - histogram.lowEdge) / histogram.binCount;
                const tsv = perInteger
                    ? chartTSV(['value', 'count'], histogram.counts.map((count, index) => [
                        histogram.min + index, count,
                    ]))
                    : chartTSV(['bin_low', 'bin_high', 'count'], histogram.counts.map((count, index) => [
                        histogram.lowEdge + (index * binSpan),
                        histogram.lowEdge + ((index + 1) * binSpan),
                        count,
                    ]));
                return Object.assign({kind, columnName, html: null, notes, empty: false}, built, {tsv});
            }
            if (kind === 'box') {
                const summary = numericSummary(numeric.values);
                const built = buildBoxPlotSVG(summary, meta);
                if (summary.outliers.length > 0) {
                    notes.push(summary.outliers.length.toLocaleString() +
                        ' outlier' + (summary.outliers.length === 1 ? '' : 's') + ' beyond 1.5×IQR');
                }
                return Object.assign({kind, columnName, html: null, notes, empty: false}, built, {
                    tsv: chartTSV(['statistic', 'value'], [
                        ['count', summary.count],
                        ['minimum', summary.min],
                        ['low_whisker', summary.lowWhisker],
                        ['q1', summary.q1],
                        ['median', summary.median],
                        ['q3', summary.q3],
                        ['high_whisker', summary.highWhisker],
                        ['maximum', summary.max],
                        ['iqr', summary.iqr],
                        ['outliers', summary.outliers.join(';')],
                    ]),
                });
            }
            const built = buildEcdfSVG(numeric.values, meta);
            return Object.assign({kind, columnName, html: null, notes, empty: false}, built, {
                tsv: chartTSV(['value', 'cumulative_fraction'], numeric.values.map((value, index) => [
                    value,
                    (index + 1) / numeric.values.length,
                ])),
            });
        }

        if (kind === 'summary') {
            const rows = columnSummaryRows(columnName, nodeIndices);
            const built = buildTableSVG(
                rows.map(row => ({statistic: row[0], value: row[1]})),
                [{key: 'statistic', label: 'Statistic'}, {key: 'value', label: 'Value', align: 'right'}],
                meta
            );
            return Object.assign({kind, columnName, notes, empty: false}, built, {
                html: summaryTableHTML(rows),
                tsv: chartTSV(['statistic', 'value'], rows),
            });
        }

        const topN = kind === 'pie'
            ? COLUMN_CHART_PIE_MAX_SLICES
            : (kind === 'bar' ? COLUMN_CHART_BAR_MAX_BARS : Infinity);
        const model = columnValueDistribution(columnName, nodeIndices, {topN});
        if (model.entries.length === 0) {
            notes.push('No values in these rows.');
            return empty;
        }
        notes.push(model.distinctTotal.toLocaleString() + ' distinct value' +
            (model.distinctTotal === 1 ? '' : 's'));
        if (model.truncated) {
            notes.push('Showing the ' + topN + ' most common; ' + model.otherDistinct.toLocaleString() +
                ' more rolled into “Other”');
        }
        if (model.nullCount > 0) {
            notes.push(model.nullCount.toLocaleString() + ' with no value (shown as —)');
        }
        const rows = frequencyTableRows(model);
        const tsv = chartTSV(
            ['value', 'count', 'percent', 'count_all_nodes', 'percent_of_all', 'color'],
            rows.map(row => [
                row.label, row.count, row.percent.toFixed(4),
                row.globalCount, row.globalPercent.toFixed(4), cssToHex(row.color),
            ])
        );
        if (kind === 'bar') {
            return Object.assign({kind, columnName, html: null, notes, empty: false},
                buildBarChartSVG(model, meta), {tsv});
        }
        if (kind === 'pie') {
            return Object.assign({kind, columnName, html: null, notes, empty: false},
                buildPieChartSVG(model, meta), {tsv});
        }
        // Frequency table. The preview caps its row count (it goes through
        // innerHTML); the TSV above always carries every value.
        //
        // The table is the one kind that shows both denominators, so it is the one
        // that has to say what the second is. Only worth spelling out once the scope
        // is narrower than the bundle: when it is not, every "% of all" reads 100%.
        if (model.rowCount < model.globalTotal) {
            notes.push('“% of all” counts each value against its own network-wide total, ' +
                'not against these ' + model.rowCount.toLocaleString() + ' rows');
        }
        const previewRows = rows.slice(0, COLUMN_CHART_TABLE_MAX_ROWS);
        if (rows.length > previewRows.length) {
            notes.push('Showing the first ' + COLUMN_CHART_TABLE_MAX_ROWS.toLocaleString() +
                ' rows; the TSV has all ' + rows.length.toLocaleString() + '.');
        }
        const built = buildTableSVG(
            previewRows.map(row => ({
                color: row.color,
                label: row.label,
                count: row.count.toLocaleString(),
                percent: row.percent.toFixed(1) + '%',
                globalCount: row.globalCount.toLocaleString(),
                globalPercent: row.globalPercent.toFixed(1) + '%',
            })),
            [
                {key: 'color', label: '', swatch: true},
                {key: 'label', label: 'Value'},
                {key: 'count', label: 'Count', align: 'right'},
                {key: 'percent', label: 'Percent', align: 'right'},
                {key: 'globalCount', label: 'All nodes', align: 'right'},
                {key: 'globalPercent', label: '% of all', align: 'right'},
            ],
            meta
        );
        return Object.assign({kind, columnName, notes, empty: false}, built, {
            html: frequencyTableHTML(previewRows),
            tsv,
        });
    }

    // ---- Anchored menu ----

    function columnChartMenuElement() {
        return document.getElementById('column-chart-menu');
    }

    function columnChartMenuIsOpen() {
        const menu = columnChartMenuElement();
        return Boolean(menu) && !menu.hidden;
    }

    function closeColumnChartMenu() {
        const menu = columnChartMenuElement();
        if (!menu || menu.hidden) {
            return;
        }
        menu.hidden = true;
        menu.innerHTML = '';
        const anchor = state.columnChartMenu && state.columnChartMenu.anchor;
        if (anchor && anchor.isConnected) {
            anchor.setAttribute('aria-expanded', 'false');
        }
        state.columnChartMenu = null;
    }

    function openColumnChartMenu(columnName, anchorButton) {
        const menu = columnChartMenuElement();
        if (!menu || !state.bundle) {
            return;
        }
        const wasOpenFor = state.columnChartMenu && state.columnChartMenu.columnName;
        closeColumnChartMenu();
        if (wasOpenFor === columnName) {
            // A second click on the same glyph closes the menu.
            return;
        }
        menu.innerHTML = columnChartKindsFor(columnName).map(kind =>
            '<button type="button" role="menuitem" data-chart-kind="' + htmlEscape(kind.id) + '">' +
            htmlEscape(kind.label) + '</button>'
        ).join('');
        menu.hidden = false;
        anchorButton.setAttribute('aria-expanded', 'true');
        // Positioned after unhiding so offsetWidth/Height are real, and flipped
        // rather than clipped when the header sits near a viewport edge.
        const anchorRect = anchorButton.getBoundingClientRect();
        const menuWidth = menu.offsetWidth;
        const menuHeight = menu.offsetHeight;
        let left = anchorRect.left;
        if (left + menuWidth > window.innerWidth - 8) {
            left = Math.max(8, anchorRect.right - menuWidth);
        }
        let top = anchorRect.bottom + 4;
        if (top + menuHeight > window.innerHeight - 8) {
            top = Math.max(8, anchorRect.top - menuHeight - 4);
        }
        menu.style.left = left + 'px';
        menu.style.top = top + 'px';
        // anchorLeft/anchorTop record where the anchor sat when we positioned
        // against it, so a later scroll or resize can tell "the anchor moved"
        // from "an event arrived".
        state.columnChartMenu = {
            columnName,
            anchor: anchorButton,
            anchorLeft: anchorRect.left,
            anchorTop: anchorRect.top,
        };
        const first = menu.querySelector('button');
        if (first) {
            // preventScroll matters: the menu is already positioned inside the
            // viewport, and letting focus() scroll an ancestor to "reveal" it
            // would fire the scroll listener below and close the menu again.
            first.focus({preventScroll: true});
        }
    }

    // The menu is anchored in viewport coordinates, so it is only stale once its
    // anchor has actually moved. Closing on the scroll event itself is not
    // enough of a test: focusing the glyph can start a scroll whose event
    // arrives *after* the menu opens, and that would close it immediately.
    function closeColumnChartMenuIfAnchorMoved() {
        const open = state.columnChartMenu;
        if (!columnChartMenuIsOpen() || !open) {
            return;
        }
        if (!open.anchor || !open.anchor.isConnected) {
            closeColumnChartMenu();
            return;
        }
        const rect = open.anchor.getBoundingClientRect();
        if (Math.abs(rect.left - open.anchorLeft) > 1 || Math.abs(rect.top - open.anchorTop) > 1) {
            closeColumnChartMenu();
        }
    }

    function moveColumnChartMenuFocus(delta) {
        const menu = columnChartMenuElement();
        if (!menu || menu.hidden) {
            return;
        }
        const items = Array.from(menu.querySelectorAll('button'));
        if (items.length === 0) {
            return;
        }
        const current = items.indexOf(document.activeElement);
        let next;
        if (delta === 'first') {
            next = 0;
        } else if (delta === 'last') {
            next = items.length - 1;
        } else {
            next = ((current < 0 ? 0 : current) + delta + items.length) % items.length;
        }
        // preventScroll for the same reason openColumnChartMenu uses it: a
        // scroll here would move the anchor and close the menu mid-navigation.
        items[next].focus({preventScroll: true});
    }

    // ---- Dialog ----

    function columnChartIsOpen() {
        const overlay = document.getElementById('column-chart-overlay');
        return Boolean(overlay) && !overlay.hidden;
    }

    // `returnFocusTo` is the header glyph the dialog was opened from, so
    // closing hands focus back where it started. It may have been detached in
    // the meantime -- updateMetadataTable rebuilds the whole header row.
    function openColumnChart(columnName, kind, returnFocusTo) {
        if (!state.bundle) {
            return;
        }
        const kinds = columnChartKindsFor(columnName);
        const resolved = kinds.some(entry => entry.id === kind) ? kind : kinds[0].id;
        state.columnChart = {columnName, kind: resolved, returnFocusTo: returnFocusTo || null};
        const select = document.getElementById('column-chart-kind');
        select.innerHTML = kinds.map(entry =>
            '<option value="' + htmlEscape(entry.id) + '"' +
            (entry.id === resolved ? ' selected' : '') + '>' + htmlEscape(entry.label) + '</option>'
        ).join('');
        document.getElementById('column-chart-overlay').hidden = false;
        renderColumnChart();
        document.getElementById('column-chart-close').focus();
    }

    function closeColumnChart() {
        const overlay = document.getElementById('column-chart-overlay');
        if (!overlay || overlay.hidden) {
            return;
        }
        overlay.hidden = true;
        document.getElementById('column-chart-preview').innerHTML = '';
        const returnFocusTo = state.columnChart && state.columnChart.returnFocusTo;
        state.columnChart = null;
        state.columnChartArtifact = null;
        if (returnFocusTo && returnFocusTo.isConnected) {
            returnFocusTo.focus({preventScroll: true});
        }
    }

    function renderColumnChart() {
        if (!state.columnChart || !state.bundle) {
            return;
        }
        const {columnName, kind} = state.columnChart;
        // A column can be renamed or deleted from under an open dialog.
        if (!state.metadataColumnIndexByName.has(columnName)) {
            closeColumnChart();
            return;
        }
        const artifact = buildColumnChartArtifact(columnName, kind, columnChartNodeIndices());
        state.columnChartArtifact = artifact;
        document.getElementById('column-chart-title').textContent =
            columnName + ' — ' + columnChartKindLabel(kind);
        const preview = document.getElementById('column-chart-preview');
        // notes[0] is always the scope line; anything after it on an empty
        // artifact is the reason there is nothing to draw, and belongs in the
        // preview box rather than leaving it blank.
        preview.innerHTML = artifact.html || artifact.svg
            || ('<p class="note">' + htmlEscape(artifact.notes.slice(1).join(' ')) + '</p>');
        const notes = artifact.notes.slice();
        if (!artifact.empty) {
            notes.push('PNG export uses the resolution set in the view toolbar.');
        }
        document.getElementById('column-chart-note').textContent = notes.join(' · ');
        ['column-chart-copy-tsv', 'column-chart-download-tsv',
            'column-chart-export-svg', 'column-chart-export-png'].forEach(id => {
            document.getElementById(id).disabled = artifact.empty;
        });
    }

    function refreshColumnChartIfOpen() {
        if (columnChartIsOpen()) {
            renderColumnChart();
        }
    }

    // ---- Exports ----

    function columnChartFileBase() {
        const artifact = state.columnChartArtifact;
        if (!artifact) {
            return exportBaseName();
        }
        return exportBaseName() + '_' + artifact.columnName.replace(/\s+/g, '_') + '_' + artifact.kind;
    }

    function exportColumnChartSVG() {
        const artifact = state.columnChartArtifact;
        if (!artifact || !artifact.svg) {
            return;
        }
        // The prolog is added only here: chartSvgOpen leaves it out so the same
        // string can be assigned to innerHTML for the preview.
        const svg = '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n' + artifact.svg;
        const href = URL.createObjectURL(new Blob([svg], {type: 'image/svg+xml'}));
        triggerDownload(href, columnChartFileBase() + '.svg', true);
    }

    function exportColumnChartPNG() {
        const artifact = state.columnChartArtifact;
        if (!artifact || !artifact.svg) {
            return;
        }
        const scaleFactor = selectedPngScale();
        const svgUrl = URL.createObjectURL(new Blob([artifact.svg], {type: 'image/svg+xml'}));
        const image = new Image();
        image.onload = () => {
            const target = document.createElement('canvas');
            target.width = Math.max(1, Math.round(artifact.width * scaleFactor));
            target.height = Math.max(1, Math.round(artifact.height * scaleFactor));
            const targetContext = target.getContext('2d');
            if (!targetContext) {
                URL.revokeObjectURL(svgUrl);
                window.alert('Could not export the chart PNG at ' + scaleFactor + '×; try a lower resolution.');
                return;
            }
            targetContext.drawImage(image, 0, 0, target.width, target.height);
            target.toBlob(blob => {
                URL.revokeObjectURL(svgUrl);
                if (!blob) {
                    window.alert('The chart is too large to export as a ' + scaleFactor + '× PNG.');
                    return;
                }
                const href = URL.createObjectURL(blob);
                const suffix = scaleFactor > 1 ? '@' + scaleFactor + 'x.png' : '.png';
                triggerDownload(href, columnChartFileBase() + suffix, true);
            }, 'image/png');
        };
        image.onerror = () => {
            URL.revokeObjectURL(svgUrl);
            window.alert('Could not rasterize the chart for PNG export.');
        };
        image.src = svgUrl;
    }

    function columnChartTSV() {
        return state.columnChartArtifact ? state.columnChartArtifact.tsv : '';
    }

    async function copyColumnChartTSV() {
        const tsv = columnChartTSV();
        if (!tsv) {
            return;
        }
        try {
            await writeTextToClipboard(tsv);
            setStatus('Copied the ' + columnChartKindLabel(state.columnChartArtifact.kind).toLowerCase() +
                ' for "' + state.columnChartArtifact.columnName + '".');
        } catch (error) {
            console.error(error);
            setStatus('Failed to copy the chart data: ' + error.message);
        }
    }

    function downloadColumnChartTSV() {
        const tsv = columnChartTSV();
        if (!tsv) {
            return;
        }
        const href = URL.createObjectURL(new Blob([tsv], {type: 'text/tab-separated-values'}));
        triggerDownload(href, columnChartFileBase() + '.tsv', true);
    }

    // ---- Wiring ----

    function setupColumnCharts() {
        const menu = columnChartMenuElement();
        const overlay = document.getElementById('column-chart-overlay');

        menu.addEventListener('click', event => {
            const item = event.target.closest('[data-chart-kind]');
            if (!item || !state.columnChartMenu) {
                return;
            }
            const {columnName, anchor} = state.columnChartMenu;
            closeColumnChartMenu();
            openColumnChart(columnName, item.dataset.chartKind, anchor);
        });

        document.addEventListener('pointerdown', event => {
            if (!columnChartMenuIsOpen()) {
                return;
            }
            if (event.target.closest('#column-chart-menu') || event.target.closest('[data-chart-column]')) {
                return;
            }
            closeColumnChartMenu();
        });
        // The scroll listener has to capture: the metadata table scrolls inside
        // .table-wrap, whose scroll events never reach window by bubbling.
        window.addEventListener('resize', closeColumnChartMenuIfAnchorMoved);
        window.addEventListener('scroll', closeColumnChartMenuIfAnchorMoved, true);

        document.addEventListener('keydown', event => {
            if (columnChartMenuIsOpen()) {
                if (event.key === 'Escape') {
                    const anchor = state.columnChartMenu && state.columnChartMenu.anchor;
                    event.preventDefault();
                    closeColumnChartMenu();
                    if (anchor && anchor.isConnected) {
                        anchor.focus({preventScroll: true});
                    }
                    return;
                }
                if (event.key === 'ArrowDown' || event.key === 'ArrowUp'
                    || event.key === 'Home' || event.key === 'End') {
                    event.preventDefault();
                    moveColumnChartMenuFocus(
                        event.key === 'ArrowDown' ? 1
                            : (event.key === 'ArrowUp' ? -1 : (event.key === 'Home' ? 'first' : 'last'))
                    );
                }
                return;
            }
            if (event.key === 'Escape' && columnChartIsOpen()) {
                event.preventDefault();
                closeColumnChart();
            }
        });

        overlay.addEventListener('click', event => {
            if (event.target === event.currentTarget) {
                closeColumnChart();
            }
        });
        document.getElementById('column-chart-close').addEventListener('click', closeColumnChart);
        document.getElementById('column-chart-kind').addEventListener('change', event => {
            if (state.columnChart) {
                state.columnChart.kind = event.target.value;
                renderColumnChart();
            }
        });
        document.getElementById('column-chart-copy-tsv').addEventListener('click', copyColumnChartTSV);
        document.getElementById('column-chart-download-tsv').addEventListener('click', downloadColumnChartTSV);
        document.getElementById('column-chart-export-svg').addEventListener('click', exportColumnChartSVG);
        document.getElementById('column-chart-export-png').addEventListener('click', exportColumnChartPNG);
    }
"""


def ssn_viewer_html(
    title: str = VIEWER_APP_NAME,
    embedded_bundle_json: bytes | None = None,
) -> str:
    def escape_html(text):
        return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

    # The browser tab keeps the bare network name: it is the part that
    # distinguishes one open viewer from another once the tab is truncated.
    escaped_title = escape_html(title)
    escaped_heading = escape_html(viewer_heading(title))
    embedded_bundle_base64 = None
    if embedded_bundle_json is not None:
        embedded_bundle_base64 = base64.b64encode(embedded_bundle_json).decode("ascii")
    layout_worker_code_json = json.dumps(_layout_worker_js())
    layout_core_js = _layout_core_js()
    # Plain-string JS modules (single braces, no f-string escaping). They are
    # interpolated into the one <script> scope below, so they share `state` and
    # can call the functions defined around them; declarations hoist.
    session_state_js = _session_state_js()
    selection_presets_js = _selection_presets_js()
    table_editing_js = _table_editing_js()
    extraction_js = _extraction_js()
    gradient_stops_js = _gradient_stops_js()
    column_charts_js = _column_charts_js()
    split_chart_hover_js = _split_chart_hover_js()
    # Generated from the same helper matrix_report uses, keyed by metric so a bundle built
    # with --merge_impact_metric product gets titles that match what its numbers mean.
    split_axis_labels_js = json.dumps({
        metric: {
            "largest": merge_impact_axis_labels(metric)["largest"],
            "movingSum": merge_impact_axis_labels(metric)["moving_sum_short"],
            # One-placeholder templates for the hover readout, so it names an impact
            # the same way matrix_report's hovertemplates do.
            "impactAmount": merge_impact_axis_labels(metric)["impact_amount"],
            "movingSumHover": merge_impact_axis_labels(metric)["moving_sum_hover"],
        }
        for metric in MERGE_IMPACT_CHOICES
    })
    # Shared with get_palette (build_ssn.py and friends) so "Domainator distinct"
    # in the viewer means the same colors those tools write.
    categorical_palettes_js = json.dumps([
        {
            "name": palette["name"],
            "label": palette["label"],
            "note": palette["note"],
            "colors": list(palette["colors"]),
        }
        for palette in NAMED_CATEGORICAL_PALETTES
    ])
    other_color_js = json.dumps(OTHER_COLOR)
    # Kept in lockstep with ssn_bundle.py so the browser and the Python readers
    # agree on which bundle revisions they understand.
    bundle_format_js = json.dumps(SSN_VIEWER_BUNDLE_FORMAT)
    bundle_version_js = json.dumps(SSN_VIEWER_BUNDLE_VERSION)
    supported_bundle_versions_js = json.dumps(list(SUPPORTED_SSN_VIEWER_BUNDLE_VERSIONS))
    domainator_version_js = json.dumps(__version__)
    viewer_app_name_js = json.dumps(VIEWER_APP_NAME)
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width, initial-scale=1.0" />
<title>{escaped_title}</title>
<style>
    :root {{
        color-scheme: light;
        --bg: #f4f5f7;
        --panel: #ffffff;
        --panel-strong: #ffffff;
        --ink: #1e2a2f;
        --muted: #5c6a70;
        --accent: #c8553d;
        --accent-soft: #f3cdb9;
        --line: #d8dce2;
        --shadow: 0 16px 40px rgba(15, 23, 42, 0.08);
    }}
    * {{ box-sizing: border-box; }}
    body {{
        margin: 0;
        font-family: Georgia, "Iowan Old Style", "Palatino Linotype", serif;
        color: var(--ink);
        background:
            radial-gradient(circle at top left, rgba(200, 85, 61, 0.06), transparent 30%),
            linear-gradient(180deg, #ffffff 0%, var(--bg) 100%);
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
        background: rgba(255, 255, 255, 0.85);
        border: 1px solid rgba(216, 220, 226, 0.8);
        box-shadow: var(--shadow);
        border-radius: 18px;
        backdrop-filter: blur(8px);
    }}
    .hero h1 {{
        margin: 0;
        font-size: clamp(1.25rem, 2vw, 1.75rem);
        font-weight: 600;
        letter-spacing: -0.02em;
    }}
    /* The name lives in the heading, so the control that edits it sits there too. */
    .hero-title {{
        display: flex;
        align-items: center;
        gap: 10px;
        flex-wrap: wrap;
    }}
    /* Sits outside `.toolbar`, so it carries its own chrome rather than inheriting. */
    button.title-edit {{
        border: 1px solid var(--line);
        background: var(--panel-strong);
        border-radius: 9px;
        padding: 2px 10px;
        font: inherit;
        font-size: 1.2rem;
        line-height: 1.4;
        color: var(--muted);
        cursor: pointer;
        transition: background 140ms ease, color 140ms ease;
    }}
    button.title-edit:hover:not(:disabled) {{
        background: #fff;
        color: var(--ink);
    }}
    button.title-edit:disabled {{
        opacity: 0.45;
        cursor: not-allowed;
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
    /* `display: flex` above would otherwise defeat the `hidden` attribute. */
    .checkbox[hidden] {{ display: none; }}
    /* A sub-option: indented under the toggle that governs it. */
    .checkbox.sub-option {{ margin-left: 22px; }}
    /* :disabled dims the box but not the <label> beside it, which is the part
       that has to look unavailable. */
    .checkbox.sub-option input:disabled + label {{ opacity: 0.5; }}
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
    .toolbar button[aria-pressed="true"]:hover:not(:disabled) {{
        background: rgba(200, 85, 61, 0.2);
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
    .cp-overlay {{
        position: fixed;
        inset: 0;
        z-index: 1000;
        display: flex;
        align-items: center;
        justify-content: center;
        background: rgba(15, 23, 42, 0.22);
    }}
    .cp-overlay[hidden] {{
        display: none;
    }}
    .cp-dialog {{
        width: min(560px, 92vw);
        max-height: 86vh;
        display: flex;
        flex-direction: column;
        background: var(--panel);
        border: 1px solid var(--line);
        border-radius: 18px;
        box-shadow: var(--shadow);
        padding: 20px 22px;
    }}
    .cp-head {{
        display: flex;
        align-items: center;
        justify-content: space-between;
        margin-bottom: 12px;
    }}
    .cp-head h3 {{
        margin: 0;
    }}
    .cp-close {{
        width: auto;
        min-width: 0;
        border: none;
        background: transparent;
        font-size: 1.4rem;
        line-height: 1;
        cursor: pointer;
        color: var(--muted);
    }}
    .cp-body {{
        overflow-y: auto;
        flex: 1 1 auto;
    }}
    .cp-pager {{
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: 10px;
        margin-bottom: 10px;
    }}
    /* As with .checkbox: `display: flex` would otherwise defeat the `hidden`
       attribute, leaving a dead pager above every column of under 100 values. */
    .cp-pager[hidden] {{
        display: none;
    }}
    .cp-pager button {{
        width: auto;
        min-width: 0;
        flex: 0 0 auto;
    }}
    .cp-pager button:disabled {{
        opacity: 0.45;
        cursor: not-allowed;
    }}
    .cp-page-status {{
        color: var(--muted);
        font-size: 0.85rem;
        flex: 1 1 auto;
        text-align: center;
    }}
    .cp-swatches {{
        display: flex;
        flex-direction: column;
        gap: 6px;
        max-height: 52vh;
        overflow-y: auto;
    }}
    .cp-swatch-row {{
        display: flex;
        align-items: center;
        gap: 10px;
    }}
    .cp-swatch-row input[type="color"] {{
        width: 38px;
        height: 28px;
        padding: 0;
        flex: 0 0 auto;
        border-radius: 8px;
    }}
    .cp-swatch-label {{
        flex: 1 1 auto;
        overflow: hidden;
        text-overflow: ellipsis;
        white-space: nowrap;
    }}
    .cp-swatch-count {{
        color: var(--muted);
        font-size: 0.85rem;
        flex: 0 0 auto;
    }}
    .cp-dialog input[type="color"] {{
        padding: 2px;
        height: 34px;
        min-height: 34px;
        cursor: pointer;
    }}
    /* The histogram, the gradient bar and the range slider are one stacked axis.
       Side padding of half a knob leaves room for the end knobs to overhang, and
       every child is `width: 100%` with a 1px border under `box-sizing: border-box`,
       so all three content boxes line up to the pixel. */
    .cp-ramp {{
        padding: 0 9px;
    }}
    .cp-histogram {{
        display: block;
        width: 100%;
        height: 90px;
        border: 1px solid var(--line);
        border-radius: 8px 8px 0 0;
        border-bottom: none;
        background: #fff;
    }}
    .cp-histogram-note {{
        margin: 5px 0 12px;
        font-size: 0.85rem;
    }}
    .cp-range {{
        position: relative;
        height: 22px;
        /* Transparent border, purely to match the 1px borders above it so the
           percentage positions inside resolve against the same content box. */
        border: 1px solid transparent;
        width: 100%;
        margin-top: 2px;
    }}
    .cp-range[hidden] {{
        display: none;
    }}
    .cp-range-track {{
        position: absolute;
        left: 0;
        right: 0;
        top: 9px;
        height: 4px;
        border-radius: 2px;
        background: var(--line);
    }}
    .cp-range-span {{
        position: absolute;
        top: 9px;
        height: 4px;
        border-radius: 2px;
        background: rgba(200, 85, 61, 0.55);
    }}
    .cp-range-knob {{
        position: absolute;
        top: 2px;
        width: 18px;
        height: 18px;
        margin-left: -9px;
        border: 1px solid var(--muted);
        border-radius: 50%;
        background: #fff;
        box-shadow: 0 1px 3px rgba(15, 23, 42, 0.22);
        cursor: ew-resize;
        touch-action: none;
    }}
    .cp-range-knob:hover {{
        border-color: var(--accent);
    }}
    .cp-range-knob:focus-visible {{
        outline: 2px solid var(--accent);
        outline-offset: 2px;
    }}
    /* Intermediate stops read as secondary to the two that bound the ramp. */
    .cp-range-knob-mid {{
        width: 14px;
        height: 14px;
        margin-left: -7px;
        top: 4px;
    }}
    .cp-stop-list {{
        display: flex;
        flex-direction: column;
        gap: 6px;
        margin-bottom: 10px;
    }}
    .cp-stop-row {{
        display: flex;
        align-items: center;
        gap: 10px;
    }}
    .cp-stop-role {{
        flex: 0 0 74px;
        color: var(--muted);
        font-size: 0.9rem;
    }}
    .cp-stop-row input[type="color"] {{
        flex: 0 0 46px;
        width: 46px;
    }}
    .cp-stop-hex {{
        flex: 0 0 104px;
        width: 104px;
        font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
        font-size: 0.88rem;
        letter-spacing: 0.02em;
        text-transform: uppercase;
        border: 1px solid var(--line);
        background: var(--panel-strong);
        border-radius: 10px;
        padding: 7px 8px;
        color: var(--ink);
    }}
    .cp-stop-row input[type="number"] {{
        flex: 1 1 auto;
        min-width: 0;
        border: 1px solid var(--line);
        background: var(--panel-strong);
        border-radius: 10px;
        padding: 7px 9px;
        font: inherit;
        color: var(--ink);
    }}
    .cp-stop-remove, .cp-stop-remove-gap {{
        flex: 0 0 30px;
        width: 30px;
        height: 30px;
    }}
    .cp-stop-remove {{
        border: 1px solid var(--line);
        border-radius: 8px;
        background: var(--panel-strong);
        color: var(--muted);
        font-size: 1.1rem;
        line-height: 1;
        cursor: pointer;
        padding: 0;
    }}
    .cp-stop-remove:hover {{
        border-color: var(--accent);
        color: var(--accent);
    }}
    .cp-stop-actions {{
        display: flex;
        gap: 10px;
    }}
    .cp-stop-actions button {{
        width: auto;
        min-width: 0;
        border: 1px solid var(--line);
        border-radius: 10px;
        background: var(--panel-strong);
        padding: 8px 12px;
        font: inherit;
        color: var(--ink);
        cursor: pointer;
    }}
    .cp-stop-actions button:hover:not(:disabled) {{
        background: #fff;
        border-color: var(--accent);
    }}
    .cp-stop-actions button:disabled {{
        opacity: 0.45;
        cursor: not-allowed;
    }}
    .cp-gradient {{
        height: 22px;
        border-radius: 0 0 8px 8px;
        border: 1px solid var(--line);
    }}
    .cp-extent {{
        display: flex;
        justify-content: space-between;
        color: var(--muted);
        font-size: 0.85rem;
    }}
    .cp-foot {{
        margin-top: 14px;
        margin-bottom: 0;
    }}
    .cp-hidden-file {{
        position: absolute;
        width: 1px;
        height: 1px;
        opacity: 0;
        pointer-events: none;
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
        border: 1px solid rgba(216, 220, 226, 0.7);
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
    /* Follows the pointer over the split chart, so it is positioned from client
       coordinates rather than parented to the (overflow:hidden) canvas wrapper.
       pointer-events:none keeps it from ever swallowing the click it advertises. */
    .split-tip {{
        position: fixed;
        z-index: 900;
        pointer-events: none;
        max-width: 320px;
        padding: 8px 10px;
        border: 1px solid var(--line);
        border-radius: 10px;
        background: var(--panel-strong);
        box-shadow: 0 10px 26px rgba(30, 42, 47, 0.18);
        color: var(--ink);
        font-size: 0.85rem;
        line-height: 1.45;
    }}
    .split-tip[hidden] {{
        display: none;
    }}
    .split-tip-row {{
        white-space: nowrap;
    }}
    .split-tip-hint {{
        margin-top: 4px;
        color: var(--muted);
        font-size: 0.78rem;
        white-space: nowrap;
    }}
    .canvas-wrap {{
        border-radius: 16px;
        overflow: hidden;
        background: #ffffff;
        border: 1px solid rgba(216, 220, 226, 0.85);
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
    /* Lives in a flex toolbar rather than under one, so it drops .note's leading gap. */
    .split-zoom-hint {{
        margin-top: 0;
        align-self: center;
    }}
    .table-wrap {{
        max-height: 680px;
        overflow: auto;
        border-radius: 12px;
        border: 1px solid rgba(216, 220, 226, 0.8);
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
        border-bottom: 1px solid rgba(216, 220, 226, 0.65);
        text-align: left;
        vertical-align: top;
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
    }}
    thead th {{
        position: sticky;
        top: 0;
        background: #ffffff;
        z-index: 1;
    }}
    th.metadata-header {{
        position: sticky;
        top: 0;
        padding: 0;
        background: #ffffff;
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
        border: 1px solid rgba(216, 220, 226, 0.9);
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
        background: rgba(216, 220, 226, 0.95);
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
    .metadata-chart-button {{
        flex: 0 0 auto;
        display: inline-flex;
        align-items: center;
        justify-content: center;
        border: 1px solid rgba(216, 220, 226, 0.9);
        border-radius: 999px;
        background: rgba(255, 255, 255, 0.9);
        color: var(--muted);
        padding: 4px 6px;
        cursor: pointer;
        transition: background 140ms ease, color 140ms ease, border-color 140ms ease;
    }}
    .metadata-chart-button:hover,
    .metadata-chart-button[aria-expanded="true"] {{
        color: #7b2f22;
        border-color: rgba(200, 85, 61, 0.4);
        background: rgba(200, 85, 61, 0.10);
    }}
    .metadata-chart-glyph {{
        display: block;
        width: 11px;
        height: 11px;
        fill: currentColor;
    }}
    /* Anchored in viewport coordinates by openColumnChartMenu. Above the
       sticky thead (z-index 1) and the modal overlays (1000). */
    .cc-menu {{
        position: fixed;
        z-index: 1100;
        min-width: 210px;
        padding: 6px;
        background: var(--panel);
        border: 1px solid var(--line);
        border-radius: 12px;
        box-shadow: var(--shadow);
    }}
    .cc-menu[hidden] {{ display: none; }}
    .cc-menu button {{
        display: block;
        width: 100%;
        padding: 7px 10px;
        border: 0;
        border-radius: 8px;
        background: transparent;
        color: inherit;
        font: inherit;
        font-size: 0.9rem;
        text-align: left;
        cursor: pointer;
    }}
    .cc-menu button:hover,
    .cc-menu button:focus-visible {{
        background: rgba(200, 85, 61, 0.10);
        color: #7b2f22;
        outline: none;
    }}
    .cc-dialog {{
        width: min(900px, 94vw);
    }}
    .cc-preview {{
        overflow: auto;
        max-height: 52vh;
        background: #ffffff;
        border: 1px solid var(--line);
        border-radius: 12px;
        padding: 10px;
    }}
    .cc-preview svg {{
        max-width: 100%;
        height: auto;
    }}
    .cc-freq-table {{
        width: 100%;
        border-collapse: collapse;
        table-layout: auto;
        font-size: 0.9rem;
    }}
    .cc-freq-table th,
    .cc-freq-table td {{
        padding: 5px 8px;
        border-bottom: 1px solid rgba(216, 220, 226, 0.65);
        text-align: left;
        white-space: nowrap;
    }}
    .cc-freq-table tbody tr:hover {{
        background: rgba(92, 106, 112, 0.06);
    }}
    .cc-freq-table tbody tr {{
        cursor: default;
    }}
    .cc-freq-table .cc-num {{
        text-align: right;
        font-variant-numeric: tabular-nums;
    }}
    .cc-swatch-cell {{
        width: 34px;
    }}
    .cc-swatch {{
        display: inline-block;
        width: 14px;
        height: 14px;
        border: 1px solid rgba(0, 0, 0, 0.25);
        border-radius: 3px;
        vertical-align: middle;
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
    .split-event-note {{
        display: flex;
        align-items: baseline;
        justify-content: space-between;
        gap: 12px;
        flex-wrap: wrap;
    }}
    .split-event-note .note {{
        flex: 1 1 320px;
    }}
    .split-event-cap {{
        display: inline-flex;
        align-items: center;
        gap: 8px;
    }}
    .split-event-cap label {{
        color: var(--muted);
        font-size: 0.9rem;
    }}
    .split-event-cap input {{
        width: 84px;
        border: 1px solid var(--line);
        background: var(--panel-strong);
        border-radius: 10px;
        padding: 4px 8px;
        font: inherit;
        color: var(--ink);
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
    .threshold-jump .threshold-step {{
        width: 32px;
        border: 1px solid var(--line);
        background: var(--panel-strong);
        border-radius: 8px;
        padding: 5px 0;
        font: inherit;
        line-height: 1.2;
        color: var(--ink);
        cursor: pointer;
    }}
    .threshold-jump .threshold-step:hover:not(:disabled) {{
        background: #fff;
    }}
    .threshold-jump .threshold-step:disabled {{
        opacity: 0.45;
        cursor: not-allowed;
    }}
    .legend-line::before,
    .legend-sum::before,
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
    .legend-sum::before {{
        width: 18px;
        height: 2px;
        background: #2f6f8f;
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
        align-items: center;
    }}
    .toolbar button.metadata-edit-toggle[aria-expanded="true"] {{
        border-color: #1e2a2f;
        background: #1e2a2f;
        color: #ffffff;
    }}
    .toolbar button.metadata-edit-toggle[aria-expanded="true"]:hover:not(:disabled) {{
        background: #3a4a52;
    }}
    .metadata-edit-panel {{
        display: flex;
        flex-wrap: wrap;
        align-items: center;
        gap: 8px;
        margin: 10px 0 0 0;
        padding: 12px 14px;
        border: 1px solid #d0d5db;
        border-radius: 10px;
        background: #eef0f3;
    }}
    .metadata-edit-panel[hidden] {{
        display: none;
    }}
    .metadata-edit-label {{
        font-size: 0.9rem;
        color: #5c6a70;
    }}
    .metadata-cell-editable {{
        cursor: cell;
    }}
    .metadata-cell-input {{
        width: 100%;
        box-sizing: border-box;
        font: inherit;
        padding: 1px 3px;
        border: 1px solid #1e2a2f;
        border-radius: 3px;
    }}
    .metadata-paste-values {{
        width: 100%;
        box-sizing: border-box;
        font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
        font-size: 0.85rem;
    }}
    .preset-slots {{
        display: inline-flex;
        align-items: center;
        gap: 4px;
    }}
    .preset-slots-label {{
        font-size: 0.9rem;
        color: #5c6a70;
        margin-right: 4px;
    }}
    .preset-slot {{
        width: 28px;
        height: 28px;
        padding: 0;
        border: 1px solid #c4cad2;
        border-radius: 6px;
        background: #ffffff;
        color: #5c6a70;
        font: inherit;
        font-size: 0.85rem;
        cursor: pointer;
    }}
    .preset-slot:disabled {{
        opacity: 0.45;
        cursor: default;
    }}
    .preset-slot:not(:disabled):hover {{
        border-color: #9aa4ad;
        background: #f0f2f4;
    }}
    .preset-slot-filled {{
        border-color: #1e2a2f;
        background: #1e2a2f;
        color: #ffffff;
        font-weight: 600;
    }}
    .preset-slot-filled:not(:disabled):hover {{
        background: #3a4a52;
        border-color: #3a4a52;
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
        <div class="hero-title">
            <h1 id="viewer-title" title="Double-click to rename this network">{escaped_heading}</h1>
            <button id="rename-network" type="button" class="title-edit" disabled aria-label="Rename this network" title="Rename this network. The name appears in the heading and the browser tab, is stored in the bundle, and names every file saved from here. You can also double-click the heading.">&#9998;</button>
        </div>
        <div class="loader">
            <input id="bundle-file" type="file" accept=".dsnv,.ssnv,.gz,.json,.ssnview" />
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
                <!-- aria-hidden: a pointer affordance duplicating data the chart's own
                     caption, the metadata table and the TSV exports already carry. -->
                <div id="split-chart-tip" class="split-tip" aria-hidden="true" hidden></div>
                <div class="split-legend">
                    <div class="split-legend-items">
                        <span class="legend-line">Largest single split</span>
                        <span class="legend-sum">Moving sum (5% window)</span>
                        <span class="legend-dot legend-threshold-readout">Current threshold <span id="threshold-label" class="legend-threshold-value">∞</span></span>
                    </div>
                    <div class="threshold-jump">
                        <button id="threshold-step-down" type="button" class="threshold-step" title="Previous split threshold (lower)" aria-label="Previous split threshold" disabled>&larr;</button>
                        <button id="threshold-step-up" type="button" class="threshold-step" title="Next split threshold (higher)" aria-label="Next split threshold" disabled>&rarr;</button>
                        <label for="threshold-input">Jump to threshold</label>
                        <!-- Deliberately type="text": a number input's spinner steps by a
                             fixed amount, which is meaningless on an axis whose splits can
                             be 0.001 or 30 apart. The buttons above and the arrow keys step
                             stop to stop instead. -->
                        <input id="threshold-input" type="text" inputmode="decimal" placeholder="Nearest split" disabled />
                    </div>
                </div>
                <div class="split-event-note">
                    <div class="note" id="split-event-count"></div>
                    <div class="split-event-cap">
                        <label for="split-event-cap">Split events</label>
                        <!-- The cap used to be build_ssn_viewer.py --max_merge_events,
                             baked into the bundle. The series is derived on load now, so
                             it is a control here instead: raising it costs one replay of
                             the merge order, not a rebuild of the file. -->
                        <input id="split-event-cap" type="number" min="0" step="50"
                               title="How many split events the chart plots and the slider offers stops for: the strongest this many by impact, plus a few so no 5% band of the axis with an event to show is left empty. 0 plots every event."
                               disabled />
                    </div>
                </div>
                <div class="toolbar">
                    <button id="split-chart-reset-zoom" type="button" disabled title="Show the whole threshold range again (or double-click the chart)">Reset zoom</button>
                    <button id="export-split-png" type="button" disabled title="Download this chart as a PNG at the resolution chosen in View Settings">Export chart PNG</button>
                    <button id="export-split-svg" type="button" disabled title="Download this chart as an editable SVG">Export chart SVG</button>
                    <!-- Teaches the gestures before one is used, and reads out the window
                         after one is; see updateSplitChartZoomControls. -->
                    <span id="split-chart-zoom-hint" class="note split-zoom-hint">Scroll to zoom · drag to pan</span>
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
            <!-- <h2>Hierarchy View</h2> -->
            <div class="panel-body">
                <div class="network-summary">
                    <div id="preset-slots" class="preset-slots" role="group" aria-label="Selection presets">
                        <span class="preset-slots-label">Presets</span>
                        <button type="button" class="preset-slot" data-preset-slot="0" disabled>0</button>
                        <button type="button" class="preset-slot" data-preset-slot="1" disabled>1</button>
                        <button type="button" class="preset-slot" data-preset-slot="2" disabled>2</button>
                        <button type="button" class="preset-slot" data-preset-slot="3" disabled>3</button>
                        <button type="button" class="preset-slot" data-preset-slot="4" disabled>4</button>
                        <button type="button" class="preset-slot" data-preset-slot="5" disabled>5</button>
                        <button type="button" class="preset-slot" data-preset-slot="6" disabled>6</button>
                        <button type="button" class="preset-slot" data-preset-slot="7" disabled>7</button>
                        <button type="button" class="preset-slot" data-preset-slot="8" disabled>8</button>
                        <button type="button" class="preset-slot" data-preset-slot="9" disabled>9</button>
                    </div>
                    <div id="hidden-summary" class="pill">0 nodes hidden</div>
                    <div id="selection-summary" class="pill">0 nodes selected</div>
                </div>
                <div class="canvas-wrap"><canvas id="cluster-view" width="1100" height="700"></canvas></div>
                <div class="note">Wheel to zoom, drag the background to pan, Shift-drag a box to select clusters, and Ctrl-click a node to toggle it individually. Hold Ctrl while Shift-dragging a box to select individual nodes within it. Hold Alt while Shift-dragging to deselect instead (combine with Ctrl to deselect individual nodes). Click a cluster bubble to toggle every node inside it. Press 0-9 to recall a selection preset, Shift+0-9 to store the current selection into one, and hover a preset button to outline its nodes in gray. Shift-click a preset button to add it to the current selection, or Alt+Shift-click to subtract it.</div>
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
                            <option value="grid">Grid (no edges)</option>
                            <option value="packed" selected>Packed (no edges)</option>
                            <option value="treemap">Treemap (no edges)</option>
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
                    <div class="control checkbox sub-option" title="Contract a chain of small clusters that only passes through into a single dashed edge carrying the chain's weakest link">
                        <input id="collapse-long-paths" type="checkbox" disabled />
                        <label for="collapse-long-paths">Collapse long paths</label>
                    </div>
                </div>
                <div class="toolbar">
                    <button id="sort-components-by-size" type="button" aria-pressed="true" disabled>Sort clusters by size: On</button>
                    <button id="focus-largest-cluster" type="button" disabled>Focus largest cluster</button>
                    <button id="focus-selection" type="button" disabled title="Center the view on the selected nodes, zooming out only if they do not already fit">Focus selection</button>
                    <button id="reset-view" disabled>Reset view</button>
                    <button id="clear-selection" disabled>Clear selection</button>
                    <button id="export-png" type="button" disabled>Export view PNG</button>
                    <select id="export-png-scale" title="PNG resolution" disabled>
                        <option value="1">1× (screen)</option>
                        <option value="2" selected>2×</option>
                        <option value="4">4×</option>
                        <option value="8">8×</option>
                    </select>
                    <button id="export-svg" type="button" disabled>Export view SVG</button>
                    <button id="customize-colors" type="button" disabled>Customize colors…</button>
                    <button id="save-session" type="button" disabled title="Download this bundle plus the current viewer state as a .dsnv session file">Save session…</button>
                    <button id="save-extraction" type="button" disabled title="Download a new .dsnv bundle containing only the selected nodes. The selection must be connected in the MST.">Save extraction…</button>
                </div>
                <div class="note" id="selection-note">Load a bundle to begin exploring metadata.</div>
            </div>
        </section>
        <section class="panel full-width">
            <h2>Node Metadata</h2>
            <div class="panel-body">
                <div class="toolbar">
                    <button id="metadata-panel-add" type="button" class="metadata-edit-toggle" aria-expanded="false" aria-controls="metadata-add-panel" data-edit-panel="add" disabled>Add column</button>
                    <button id="metadata-panel-set" type="button" class="metadata-edit-toggle" aria-expanded="false" aria-controls="metadata-set-panel" data-edit-panel="set" disabled>Set column</button>
                    <button id="metadata-panel-rename" type="button" class="metadata-edit-toggle" aria-expanded="false" aria-controls="metadata-rename-panel" data-edit-panel="rename" disabled>Rename column</button>
                    <button id="metadata-panel-delete" type="button" class="metadata-edit-toggle" aria-expanded="false" aria-controls="metadata-delete-panel" data-edit-panel="delete" disabled>Delete column</button>
                    <button id="metadata-panel-select" type="button" class="metadata-edit-toggle" aria-expanded="false" aria-controls="metadata-select-panel" data-edit-panel="select" disabled>Select by value</button>
                </div>
                <div id="metadata-add-panel" class="metadata-edit-panel" hidden>
                    <label for="metadata-new-column-name" class="metadata-edit-label">Name</label>
                    <input id="metadata-new-column-name" type="text" placeholder="new column name" />
                    <select id="metadata-new-column-type">
                        <option value="text" selected>Text</option>
                        <option value="number">Number</option>
                    </select>
                    <button id="metadata-add-column" type="button">Add</button>
                    <span class="metadata-edit-label">or</span>
                    <button id="metadata-add-cluster-column" type="button" title="Number every node by the cluster it belongs to at the current threshold, largest cluster first — the same column build_ssn.py --cluster writes. Defaults to the name SSN_cluster.">Fill with cluster numbers</button>
                </div>
                <div id="metadata-set-panel" class="metadata-edit-panel" hidden>
                    <select id="metadata-fill-column"></select>
                    <span class="metadata-edit-label">=</span>
                    <input id="metadata-fill-value" type="text" placeholder="value (blank clears)" />
                    <label for="metadata-fill-target" class="metadata-edit-label">for</label>
                    <select id="metadata-fill-target">
                        <option value="selection" selected>Selected nodes</option>
                        <option value="rows">Staged table rows</option>
                        <option value="page">Rows on this page</option>
                        <option value="all">All nodes</option>
                    </select>
                    <button id="metadata-fill-apply" type="button">Apply</button>
                    <button id="metadata-paste-open" type="button">Paste column…</button>
                </div>
                <div id="metadata-rename-panel" class="metadata-edit-panel" hidden>
                    <select id="metadata-rename-column"></select>
                    <span class="metadata-edit-label">&rarr;</span>
                    <input id="metadata-rename-value" type="text" placeholder="new column name" />
                    <button id="metadata-rename-apply" type="button">Rename</button>
                </div>
                <div id="metadata-delete-panel" class="metadata-edit-panel" hidden>
                    <select id="metadata-delete-column"></select>
                    <button id="metadata-delete-column-apply" type="button">Delete</button>
                    <span class="metadata-edit-label" id="metadata-delete-note"></span>
                </div>
                <div id="metadata-select-panel" class="metadata-edit-panel" hidden>
                    <select id="metadata-select-column" aria-label="Column to compare"></select>
                    <select id="metadata-select-op" aria-label="Comparison"></select>
                    <input id="metadata-select-value" type="text" placeholder="value" aria-label="Value to compare against" />
                    <span class="metadata-edit-label" id="metadata-select-and" hidden>and</span>
                    <input id="metadata-select-value2" type="text" placeholder="upper bound" aria-label="Upper bound" hidden />
                    <button id="metadata-select-add" type="button" title="Search every node in the network and add the matches to the selection.">Add to selection</button>
                    <button id="metadata-select-remove" type="button" title="Search the selected nodes and drop the matches from the selection." disabled>Remove from selection</button>
                    <button id="metadata-select-subset" type="button" title="Search the selected nodes and keep only the matches." disabled>Subset selection</button>
                    <!-- The note is the only channel for "that regex is malformed", so it
                         announces itself rather than waiting to be read. -->
                    <span class="metadata-edit-label" id="metadata-select-note" role="status"></span>
                </div>
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
                <div class="note">Double-click any metadata cell to edit it. node_id is read-only.</div>
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
<div id="name-overlay" class="cp-overlay" hidden>
    <div class="cp-dialog" role="dialog" aria-modal="true" aria-labelledby="name-dialog-title">
        <div class="cp-head">
            <h3 id="name-dialog-title">Rename network</h3>
            <button id="name-cancel" type="button" class="cp-close" aria-label="Close">×</button>
        </div>
        <div class="cp-body">
            <div class="control">
                <label for="name-value" id="name-dialog-label">Network name</label>
                <input id="name-value" type="text" />
            </div>
            <p class="note" id="name-dialog-note"></p>
        </div>
        <div class="toolbar cp-foot">
            <button id="name-apply" type="button">Apply</button>
        </div>
    </div>
</div>
<div id="metadata-paste-overlay" class="cp-overlay" hidden>
    <div class="cp-dialog" role="dialog" aria-modal="true" aria-labelledby="metadata-paste-title">
        <div class="cp-head">
            <h3 id="metadata-paste-title">Paste a column of values</h3>
            <button id="metadata-paste-cancel" type="button" class="cp-close" aria-label="Close">×</button>
        </div>
        <div class="cp-body">
            <div class="control">
                <label for="metadata-paste-column">Paste into column</label>
                <select id="metadata-paste-column"></select>
            </div>
            <p class="note">One value per line, applied to the rows on this page in the order they are displayed — the inverse of a column's Copy button. This page shows <span id="metadata-paste-count">0</span> rows, and the paste must have exactly that many lines.</p>
            <textarea id="metadata-paste-values" rows="10" class="metadata-paste-values" placeholder="one value per line"></textarea>
        </div>
        <div class="toolbar cp-foot">
            <button id="metadata-paste-apply" type="button">Apply</button>
        </div>
    </div>
</div>
<div id="color-picker-overlay" class="cp-overlay" hidden>
    <div id="color-picker-dialog" class="cp-dialog" role="dialog" aria-modal="true" aria-labelledby="color-picker-title">
        <div class="cp-head">
            <h3 id="color-picker-title">Customize colors</h3>
            <button id="color-picker-close" type="button" class="cp-close" aria-label="Close">×</button>
        </div>
        <div class="cp-body">
            <p class="note" id="color-picker-empty-note" hidden>Choose a column in “Color by” to customize its colors.</p>
            <div class="control checkbox" id="color-categorical-control" hidden>
                <input id="color-as-categorical" type="checkbox" />
                <label for="color-as-categorical">Treat values as categories (discrete colors)</label>
            </div>
            <p class="note" id="color-categorical-note" hidden></p>
            <div id="color-picker-discrete" hidden>
                <div class="control">
                    <label for="color-palette">Color palette</label>
                    <select id="color-palette"></select>
                </div>
                <p class="note" id="color-palette-note"></p>
                <div class="cp-pager" id="color-picker-pager" hidden>
                    <button id="color-picker-prev" type="button">‹ Prev</button>
                    <span id="color-picker-page-status" class="cp-page-status"></span>
                    <button id="color-picker-next" type="button">Next ›</button>
                </div>
                <div id="color-picker-swatch-list" class="cp-swatches"></div>
            </div>
            <div id="color-picker-continuous" hidden>
                <div id="color-histogram-wrap" class="cp-ramp">
                    <canvas id="color-histogram" class="cp-histogram" width="1040" height="180"></canvas>
                    <div class="cp-gradient" id="color-gradient-preview"></div>
                    <div class="cp-range" id="color-range-slider">
                        <div class="cp-range-track"></div>
                        <div class="cp-range-span" id="color-range-span"></div>
                    </div>
                    <div class="cp-extent"><span id="color-min-label"></span><span id="color-max-label"></span></div>
                    <p class="note cp-histogram-note" id="color-histogram-note"></p>
                </div>
                <div id="color-stop-list" class="cp-stop-list"></div>
                <div class="cp-stop-actions">
                    <button id="color-add-stop" type="button">Add stop</button>
                    <button id="color-reset-values" type="button">Reset values to data range</button>
                </div>
            </div>
            <div class="control" id="color-null-control" hidden><label for="color-null">No-value color</label><input id="color-null" type="color" /></div>
        </div>
        <div class="toolbar cp-foot">
            <button id="color-picker-reset" type="button">Reset to defaults</button>
            <input id="color-table-file" type="file" accept=".tsv,.txt" class="cp-hidden-file" />
            <button id="load-color-table" type="button">Load color table…</button>
            <button id="save-color-table" type="button">Save color table</button>
            <button id="export-legend-svg" type="button">Export legend SVG</button>
            <button id="export-legend-png" type="button">Export legend PNG</button>
        </div>
    </div>
</div>
<div id="column-chart-menu" class="cc-menu" role="menu" aria-label="Charts and tables for this column" hidden></div>
<div id="column-chart-overlay" class="cp-overlay" hidden>
    <div id="column-chart-dialog" class="cp-dialog cc-dialog" role="dialog" aria-modal="true" aria-labelledby="column-chart-title">
        <div class="cp-head">
            <h3 id="column-chart-title">Column chart</h3>
            <button id="column-chart-close" type="button" class="cp-close" aria-label="Close">×</button>
        </div>
        <div class="cp-body">
            <div class="toolbar cc-controls">
                <label for="column-chart-kind">Chart</label>
                <select id="column-chart-kind"></select>
            </div>
            <div id="column-chart-preview" class="cc-preview"></div>
            <p class="note" id="column-chart-note"></p>
        </div>
        <div class="toolbar cp-foot">
            <button id="column-chart-copy-tsv" type="button">Copy TSV</button>
            <button id="column-chart-download-tsv" type="button">Download TSV</button>
            <button id="column-chart-export-svg" type="button">Export SVG</button>
            <button id="column-chart-export-png" type="button">Export PNG</button>
        </div>
    </div>
</div>
<script>
    const state = {{
        bundle: null,
        bundleVersionWarning: '',
        metadataColumns: [],
        metadataRows: [],
        metadataByNodeIndex: [],
        metadataColumnByName: new Map(),
        metadataColumnIndexByName: new Map(),
        metadataColorInfoByName: new Map(),
        customPalettes: {{}},
        // Derived cache of defaultCategoricalPalette() results, keyed like
        // customPalettes. Not user data, so it is never saved into a session;
        // rebuildMetadataCaches() drops it.
        defaultPalettes: new Map(),
        // Names of numeric metadata columns the viewer should color as discrete
        // categories instead of a gradient (e.g. integer cluster numbers).
        categoricalColumns: new Set(),
        colorPickerPage: 0,
        // Where drawSplitChart() last painted each mark, for hover hit-testing.
        splitChartHit: null,
        // The stop the slider was last moved to deliberately; see currentSliderStop().
        selectedStop: null,
        // The stop in effect before the current split-chart click gesture began.
        stopBeforeSplitChartClick: null,
        // {{min, max}} in threshold units while the split chart is zoomed, else null.
        splitChartZoom: null,
        // {{pointerId, startX, lastX, moved}} while the split chart is being panned.
        splitChartDrag: null,
        splitChartSuppressClick: false,
        // {{columnName, anchor}} while a column header's chart menu is open.
        columnChartMenu: null,
        // {{columnName, kind}} while the chart dialog is open, and the built
        // artifact the four export buttons hand out (so they never rebuild it).
        columnChart: null,
        columnChartArtifact: null,
        // Binned counts for the gradient dialog's histogram. Cached because binning walks
        // every node, while recoloring the bars happens on every keystroke in the dialog.
        colorHistogram: null,
        // Working list of {{value, color}} gradient stops while the picker is open.
        gradientStops: [],
        metadataColumnWidths: new Map(),
        metadataSearchTextByNodeIndex: [],
        metadataPage: 0,
        activeClusters: [],
        visibleClusters: [],
        selectedNodeIndices: new Set(),
        selectedMetadataNodeIndices: new Set(),
        selectionPresets: new Map(),
        presetPreviewSlot: null,
        metadataEditCell: null,
        pendingNameAction: null,
        nodeIndexById: null,
        pendingRestoreViewTransform: null,
        sessionMissingNodeIds: 0,
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
    const BUNDLE_FORMAT = {bundle_format_js};
    const BUNDLE_VERSION = {bundle_version_js};
    const SUPPORTED_BUNDLE_VERSIONS = {supported_bundle_versions_js};
    const DOMAINATOR_VERSION = {domainator_version_js};
    const VIEWER_APP_NAME = {viewer_app_name_js};

    // The shell viewer is built before any bundle exists, and the file picker can
    // swap bundles at any time, so the heading and tab follow the loaded bundle.
    function updateViewerTitle(networkName) {{
        const name = (networkName || '').trim();
        document.getElementById('viewer-title').textContent =
            name === '' ? VIEWER_APP_NAME : VIEWER_APP_NAME + ': ' + name;
        if (name !== '') {{
            document.title = name;
        }}
    }}

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
        if (bundle.format !== BUNDLE_FORMAT) {{
            throw new Error('Unsupported bundle format: ' + bundle.format);
        }}
        // A wrong `format` is fatal, but an unrecognized `version` is not: bundle
        // revisions so far have only added sections, so loading and reporting beats
        // refusing a file this viewer can very likely still render.
        state.bundleVersionWarning = SUPPORTED_BUNDLE_VERSIONS.includes(bundle.version)
            ? ''
            : ' (bundle version ' + bundle.version + ' is not one this viewer knows about ('
                + SUPPORTED_BUNDLE_VERSIONS.join(', ') + '); loading anyway, some features may misbehave)';
        state.bundle = bundle;
        updateViewerTitle(bundle.name);
        state.metadataColumns = bundle.metadata.columns || [];
        state.metadataRows = bundle.metadata.rows || [];
        state.metadataByNodeIndex = state.metadataRows;
        state.customPalettes = {{}};
        state.categoricalColumns = new Set(bundle.defaults.categorical_columns || []);
        state.metadataSort = {{columnKey: null, direction: 'asc'}};
        state.metadataColumnWidths = new Map();
        state.metadataPage = 0;
        state.metadataResize = null;
        state.selectedNodeIndices = new Set();
        state.selectedMetadataNodeIndices = new Set();
        state.selectionPresets = new Map();
        state.metadataRowSelectionAnchor = null;
        state.nodeIndexById = null;
        state.pendingRestoreViewTransform = null;
        state.dotLayoutCache = new Map();
        state.layoutCache = new Map();
        // A window onto the previous bundle's threshold range means nothing here; a
        // session restore sets its own afterwards. Cleared before the series is built,
        // because the selection is taken against whatever window is in force.
        state.splitChartZoom = null;
        // Derived here rather than read from the file: see deriveMergeSeries. Every
        // consumer below reads state.series, never the bundle's own copy (a pre-v6 file
        // has one, already capped by whoever built it).
        rebuildMergeSeries(DEFAULT_WINDOWED_MAX_MERGE_EVENTS);
        state.allNodeIndices = bundle.graph.nodes.map((_, i) => i);
        rebuildMetadataCaches();

        // Must precede rebuildNodeColorCache: the cache is keyed off the "Color by"
        // selection, which populateMetadataControls sets from bundle.defaults.
        populateMetadataControls();
        rebuildNodeColorCache();
        const slider = document.getElementById('threshold-slider');
        slider.max = String(state.sliderModel.maxPosition);
        slider.value = String(state.sliderModel.initialPosition);
        slider.disabled = state.sliderModel.stops.length === 0;
        document.getElementById('threshold-min-label').textContent = state.sliderModel.minLabel;
        document.getElementById('threshold-max-label').textContent = state.sliderModel.maxLabel;
        document.getElementById('sort-components-by-size').disabled = false;
        document.getElementById('focus-largest-cluster').disabled = false;
        document.getElementById('reset-view').disabled = false;
        document.getElementById('export-png').disabled = false;
        document.getElementById('export-png-scale').disabled = false;
        document.getElementById('export-svg').disabled = false;
        document.getElementById('export-split-png').disabled = false;
        document.getElementById('export-split-svg').disabled = false;
        document.getElementById('customize-colors').disabled = false;
        document.getElementById('save-session').disabled = false;
        document.getElementById('rename-network').disabled = false;
        state.presetPreviewSlot = null;
        if (colorPickerIsOpen()) {{
            openColorPicker();
        }}
        document.getElementById('threshold-input').disabled = state.sliderModel.stops.length === 0;
        updateThresholdStepButtons();
        document.getElementById('metadata-select-nodes').disabled = true;
        document.getElementById('metadata-deselect-rows').disabled = true;
        document.getElementById('metadata-reset-sort').disabled = true;
        document.getElementById('metadata-filter').disabled = false;
        document.getElementById('metadata-filter').value = '';
        document.getElementById('metadata-null-order').disabled = false;
        document.getElementById('metadata-null-order').value = 'last';
        document.getElementById('metadata-rows-per-page').disabled = false;
        document.getElementById('metadata-panel-add').disabled = false;
        document.getElementById('metadata-panel-set').disabled = false;
        document.getElementById('metadata-panel-select').disabled = false;
        document.getElementById('metadata-panel-rename').disabled = false;
        document.getElementById('metadata-panel-delete').disabled = false;
        closeMetadataEditPanels();
        populateEditableColumnMenus();
        document.getElementById('metadata-prev-page').disabled = true;
        document.getElementById('metadata-next-page').disabled = true;
        // Must run after the slider's `max` is set (a value above the old max
        // would be clamped) and after populateMetadataControls has built the
        // color-by/label-by options, but before the color cache is used.
        let sessionNote = '';
        if (bundle.app_state) {{
            sessionNote = applySessionState(bundle.app_state);
            rebuildNodeColorCache();
            if (colorPickerIsOpen()) {{
                openColorPicker();
            }}
        }}
        setStatus(
            'Loaded ' + (bundle.name || 'bundle') + ' with ' + bundle.graph.nodes.length.toLocaleString() + ' nodes.'
            + state.bundleVersionWarning + sessionNote
        );
        renderPresetSlots();
        updateSplitEventCount();
        updateThresholdUI();
        updateMetadataTable();
    }}

    function rebuildMetadataCaches() {{
        // Which default color a value gets depends on the whole sorted set of
        // values in its column, so any change to the metadata invalidates it.
        state.defaultPalettes = new Map();
        state.metadataColumnByName = new Map();
        state.metadataColumnIndexByName = new Map();
        state.metadataColorInfoByName = new Map();

        state.metadataColumns.forEach((column, columnIndex) => {{
            state.metadataColumnByName.set(column.name, column);
            state.metadataColumnIndexByName.set(column.name, columnIndex);
            if (column.type === 'int' || column.type === 'float') {{
                // `baseType` is the column's intrinsic kind; `type` is how it is colored
                // right now, which the categorical toggle can flip for numeric columns.
                state.metadataColorInfoByName.set(column.name, {{type: 'numeric', baseType: 'numeric', min: Infinity, max: -Infinity}});
                return;
            }}
            state.metadataColorInfoByName.set(column.name, {{type: 'categorical', baseType: 'categorical'}});
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

        state.metadataColorInfoByName.forEach((info, columnName) => {{
            if (info.baseType !== 'numeric') {{
                return;
            }}
            if (!Number.isFinite(info.min) || !Number.isFinite(info.max)) {{
                info.min = 0;
                info.max = 0;
            }}
            // min/max stay available so toggling back to the gradient is free.
            info.type = state.categoricalColumns.has(columnName) ? 'categorical' : 'numeric';
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
        // A palette edit that recolors the nodes should recolor an open chart
        // of the same column too.
        refreshColumnChartIfOpen();
    }}

    // ---- Slider positions ----
    //
    // Positions are a warp of the threshold axis rather than a straight scale of it,
    // because two things pull on them at once:
    //
    //   * Every stop should be reachable by dragging, and 1001 integer positions have
    //     to carry every stop there is. On a network whose merges bunch into a narrow
    //     band of scores, 401 stops once landed on 103 positions -- eleven sharing one
    //     -- which left 74% of them impossible to select by dragging at all.
    //   * The track should still read like the chart above it, so that dragging toward
    //     a spike on the chart moves toward that spike.
    //
    // Separating crowded stops has to borrow track from somewhere, so the two cannot
    // both hold exactly. What holds instead, in two steps:
    //
    //   1. Every stop is given one position of its own, which is what makes it
    //      draggable, and
    //   2. the track left over is spent saying where the stops are -- proportionally to
    //      threshold, except that the stretch the chart is zoomed into carries the
    //      chart's own magnification times the track density of the rest, so zooming
    //      the chart magnifies that stretch of the slider too.
    //
    // Where stops are spread out, step 1 costs almost nothing and the result is the
    // plain proportional map this started as. Where they crowd, step 1 is the whole
    // budget and they come out evenly spaced -- which is the most a thousand positions
    // can say about a thousand stops.

    // Finite stops share 0..920; the infinity stop sits alone at 1000, far enough above
    // them to read as the separate thing it is.
    const SLIDER_MAX_FINITE_POSITION = 920;
    const SLIDER_INFINITY_POSITION = 1000;

    // Threshold -> fraction of the finite track. Piecewise linear: one density inside
    // the chart's window and another outside it, integrated along the axis, so the map
    // stays monotone and stays proportional *within* each region.
    function sliderValueWarp(lowValue, highValue) {{
        const span = highValue - lowValue;
        if (!(span > 0)) {{
            return () => 0;
        }}
        const proportional = value => (value - lowValue) / span;
        const zoom = state.splitChartZoom;
        if (!zoom) {{
            return proportional;
        }}
        const windowLow = Math.min(Math.max(zoom.min, lowValue), highValue);
        const windowHigh = Math.max(Math.min(zoom.max, highValue), lowValue);
        const windowSpan = windowHigh - windowLow;
        if (!(windowSpan > 0)) {{
            return proportional;
        }}
        // The chart's own magnification. Because the outside keeps its own width, the
        // window's share of the track tends to half however far the zoom goes, rather
        // than swallowing the track and stranding everything else on one position.
        const weight = span / windowSpan;
        const total = (span - windowSpan) + (weight * windowSpan);
        return value => {{
            const below = Math.min(Math.max(value, lowValue), windowLow) - lowValue;
            const inside = Math.min(Math.max(value, windowLow), windowHigh) - windowLow;
            const above = Math.min(Math.max(value, windowHigh), highValue) - windowHigh;
            return (below + (weight * inside) + above) / total;
        }};
    }}

    // Assigns `sliderPosition` in place, so a stop object held elsewhere -- notably
    // state.selectedStop -- keeps pointing at the same stop across a remap.
    //
    // Every stop is given a slot of its own first, and only the track left over after
    // that is spent on saying where the stops are, through the warp. That ordering is
    // what makes the two pulls compatible: the reserved slot is what makes a stop
    // draggable, and the surplus is what keeps the sparse stretches proportional and
    // hands the chart's window its extra room. Positions come out strictly increasing
    // by construction -- consecutive stops differ by at least the one slot, and rounding
    // a non-decreasing sequence cannot close a gap of one -- so no separation pass is
    // needed, and none of the cascading that a pass has to do when a crowd overflows
    // the end of the track can undo the warp.
    function positionSliderStops(finiteStops) {{
        if (finiteStops.length === 0) {{
            return;
        }}
        const lastIndex = finiteStops.length - 1;
        const lowValue = finiteStops[0].threshold_value;
        const highValue = finiteStops[lastIndex].threshold_value;
        const surplus = SLIDER_MAX_FINITE_POSITION - lastIndex;
        if (surplus < 0 || !(highValue - lowValue > 0)) {{
            // More stops than the track has positions -- only reachable by raising
            // --max_merge_events past ~900 -- or every stop at one threshold. Neither
            // leaves anything to be proportional to, so fall back to rank.
            finiteStops.forEach((stop, stopIndex) => {{
                stop.sliderPosition = Math.round(
                    (stopIndex / Math.max(1, lastIndex)) * SLIDER_MAX_FINITE_POSITION);
            }});
            return;
        }}
        const warp = sliderValueWarp(lowValue, highValue);
        finiteStops.forEach((stop, stopIndex) => {{
            stop.sliderPosition = stopIndex + Math.round(surplus * warp(stop.threshold_value));
        }});
    }}

    // Both which stops exist and where they sit are functions of the chart's window, so
    // both are redone whenever that window moves -- and the thumb put back on the
    // threshold it was already on, which is a different position, and possibly a
    // different stop object, under the new map.
    //
    // The selection is re-run because the cap spends itself on the window: zoom in and
    // the same `maxMergeEvents` buys stops among the events actually on screen, which is
    // what makes small events reachable at all. The threshold in effect is captured
    // first and pinned through the re-selection, so the cut never changes under a zoom.
    function repositionSliderStops() {{
        if (!state.series || !state.sliderModel || state.sliderModel.stops.length === 0) {{
            return;
        }}
        const threshold = selectedThresholdValue();
        applyMergeSelection(threshold);
        restoreThreshold(threshold);
        updateThresholdStepButtons();
        updateSplitEventCount();
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

        const positionedStops = finiteStops.map(stop => ({{...stop}}));
        positionSliderStops(positionedStops);
        positionedStops.push({{...infinityStop, sliderPosition: SLIDER_INFINITY_POSITION}});

        return {{
            stops: positionedStops,
            maxPosition: SLIDER_INFINITY_POSITION,
            initialPosition: 0,
            minLabel: finiteStops[0].threshold_label,
            maxLabel: '∞',
        }};
    }}

    // The slider's position is the source of truth, but not a sufficient one, so a stop
    // chosen deliberately is remembered and honored while the slider still sits where
    // that choice put it. Two reasons it has to be:
    //
    //   * positionSliderStops re-lays-out the track whenever the chart's window moves,
    //     and the stop the thumb is on has to survive being given a new position; and
    //   * with more stops than the track has positions -- only reachable by raising
    //     --max_merge_events past ~900 -- stops share positions again, and resolving by
    //     position alone always hands back the first stop at that position. Snapping to
    //     any of the others and asking again would then return a different stop than the
    //     one just chosen, which is what made the step arrows look stuck partway up the
    //     axis before the track was warped.
    //
    // Dragging the slider moves it off that position and the nearest-position search
    // below takes over again, which is right: a drag picks a position, not a stop.
    function currentSliderStop() {{
        if (!state.sliderModel || state.sliderModel.stops.length === 0) {{
            return null;
        }}
        const sliderPosition = Number(document.getElementById('threshold-slider').value);
        if (state.selectedStop && state.selectedStop.sliderPosition === sliderPosition) {{
            return state.selectedStop;
        }}
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
        state.selectedStop = stop;
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

    // Static menu of the named palettes. The disabled first entry is what a
    // hand-edited column shows; "Reset to defaults" is the way back to
    // DEFAULT_CATEGORICAL_PALETTE, so the menu carries no pseudo-entry for it.
    function populateColorPaletteMenu() {{
        const select = document.getElementById('color-palette');
        select.innerHTML = '';
        const customOption = document.createElement('option');
        customOption.value = '';
        customOption.textContent = 'Custom colors';
        customOption.disabled = true;
        select.appendChild(customOption);
        CATEGORICAL_PALETTES.forEach(scheme => {{
            const option = document.createElement('option');
            option.value = scheme.name;
            option.textContent = scheme.label;
            option.title = scheme.note;
            select.appendChild(option);
        }});
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

        // "Label by" can also use the sequence's node_id (not a metadata column).
        const nodeIdOption = document.createElement('option');
        nodeIdOption.value = '__node_id__';
        nodeIdOption.textContent = 'node_id';
        labelBy.appendChild(nodeIdOption);

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
            // Strictly-above, matching `build_ssn --lb thresholdValue`: a component that
            // merged exactly at the threshold is split back apart here.
            if (component.threshold <= thresholdValue) {{
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

    function collapseLongPathsEnabled() {{
        return document.getElementById('collapse-long-paths').checked;
    }}

    function updateCollapseLongPathsControl() {{
        // Sub-option of leaf pruning, and only meaningful there: with leaf pruning off
        // every below-minimum cluster is dropped outright, so no pass-through chain is
        // left to contract.
        document.getElementById('collapse-long-paths').disabled = !leafPruningOnlyEnabled();
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
        if (layoutMode === 'treemap') {{
            // Treemap node positions are fixed by the lattice, so sorting has no effect.
            button.disabled = true;
            button.textContent = 'Sort clusters by size: n/a';
            button.title = 'Treemap nodes follow a fixed lattice; size sorting does not apply.';
            return;
        }}
        button.disabled = false;
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

    // Distinct values of a column, ordered by frequency (desc) then key (asc). Keyed by
    // String(value) -- the same key categoricalColor hashes and the TSV/legend use. Null/
    // empty cells are counted separately (nullCount) and never become a key. `cap` bounds
    // the returned list (the picker/legend cap at 200; the TSV save path passes Infinity so
    // the file round-trips); `overflow` flags that more distinct values exist than returned.
    function distinctColumnValues(columnName, cap = 200) {{
        const result = {{values: [], nullCount: 0, overflow: false, total: 0}};
        const columnIndex = state.metadataColumnIndexByName.get(columnName);
        if (columnIndex === undefined) {{
            return result;
        }}
        const byKey = new Map();
        state.metadataByNodeIndex.forEach(row => {{
            const raw = row ? (row[columnIndex] ?? null) : null;
            if (raw === null || raw === undefined || raw === '') {{
                result.nullCount += 1;
                return;
            }}
            const key = String(raw);
            const existing = byKey.get(key);
            if (existing) {{
                existing.count += 1;
            }} else {{
                byKey.set(key, {{key, label: formatValue(raw), raw, count: 1}});
            }}
        }});
        result.total = byKey.size;
        const sorted = Array.from(byKey.values()).sort(
            (a, b) => (b.count - a.count) || a.key.localeCompare(b.key)
        );
        if (sorted.length > cap) {{
            result.overflow = true;
            result.values = sorted.slice(0, cap);
        }} else {{
            result.values = sorted;
        }}
        return result;
    }}

    // Normalize a color string to '#RRGGBB' (uppercase), mirroring Python's
    // normalize_color_hex. Returns null for anything that isn't a 6-digit hex.
    // Mirrors color_genbank.normalize_color_hex: an optional '#' and six hex digits,
    // canonicalized to '#RRGGBB'. Loosened by one case for hand-typed input -- CSS-style
    // three-digit shorthand, where each digit doubles. Returns null for anything else.
    function normalizeColorHex(text) {{
        const trimmed = String(text).trim().toUpperCase();
        const short = /^#?([0-9A-F]{{3}})$/.exec(trimmed);
        if (short) {{
            const digits = short[1];
            return '#' + digits.charAt(0) + digits.charAt(0)
                + digits.charAt(1) + digits.charAt(1)
                + digits.charAt(2) + digits.charAt(2);
        }}
        const match = /^#?([0-9A-F]{{6}})$/.exec(trimmed);
        return match ? '#' + match[1] : null;
    }}

    // Interpolate two '#rrggbb' colors channel-wise in RGB. t is clamped to [0, 1].
    function lerpHexColor(a, b, t) {{
        const clampT = Math.max(0, Math.min(1, t));
        const parse = hex => {{
            const value = parseInt(hex.slice(1), 16);
            return [(value >> 16) & 255, (value >> 8) & 255, value & 255];
        }};
        const left = parse(a);
        const right = parse(b);
        const channel = index => Math.round(left[index] + ((right[index] - left[index]) * clampT));
        const toHex = n => n.toString(16).padStart(2, '0');
        return '#' + toHex(channel(0)) + toHex(channel(1)) + toHex(channel(2));
    }}

    // Named qualitative palettes (from domainator.utils.NAMED_CATEGORICAL_PALETTES) and
    // the neutral gray get_palette gives values with no color of their own.
    const CATEGORICAL_PALETTES = {categorical_palettes_js};
    const PALETTE_NO_VALUE_COLOR = {other_color_js};

    function paletteSchemeByName(schemeName) {{
        return CATEGORICAL_PALETTES.find(scheme => scheme.name === schemeName) || null;
    }}

    // Python's str comparison is by code point; localeCompare is not, and the assignment
    // order has to match domainator.utils.sort_palette_values to reproduce its colors.
    function comparePaletteText(left, right) {{
        if (left < right) {{
            return -1;
        }}
        return left > right ? 1 : 0;
    }}

    // Mirror of domainator.utils.sort_palette_values: numeric-looking keys sort
    // numerically (ties broken by text), everything else sorts by text. Nulls never
    // reach here -- distinctColumnValues counts them separately.
    function sortPaletteKeys(keys) {{
        const numeric = [];
        const other = [];
        keys.forEach(key => {{
            const number = Number(key);
            if (key.trim() !== '' && Number.isFinite(number)) {{
                numeric.push({{number, key}});
            }} else {{
                other.push(key);
            }}
        }});
        numeric.sort((a, b) => (a.number - b.number) || comparePaletteText(a.key, b.key));
        other.sort(comparePaletteText);
        return numeric.map(entry => entry.key).concat(other);
    }}

    // Fill a column's palette from a named scheme: colors are handed out in
    // sort_palette_values order and cycle once the scheme runs out, exactly like
    // get_palette. Values keep their assigned color afterwards, so individual swatches,
    // the color-table TSV and the legend all stay editable/exportable.
    function applyNamedPalette(columnName, schemeName) {{
        const scheme = paletteSchemeByName(schemeName);
        if (!scheme) {{
            return null;
        }}
        const distinct = distinctColumnValues(columnName, Infinity);
        const colors = {{}};
        sortPaletteKeys(distinct.values.map(entry => entry.key)).forEach((key, index) => {{
            colors[key] = scheme.colors[index % scheme.colors.length];
        }});
        const palette = ensureCategoricalPalette(columnName);
        palette.colors = colors;
        palette.nullColor = PALETTE_NO_VALUE_COLOR;
        palette.scheme = scheme.name;
        return {{scheme, assigned: distinct.values.length}};
    }}

    // The palette a column has before anyone customizes it. Values used to be
    // colored by a hue hashed from their text -- stable, but arbitrary, and free to
    // put two adjacent categories on nearly the same color. They now get
    // get_palette()'s own 64 distinct colors in sort_palette_values order, so an
    // untouched viewer agrees with what build_ssn.py would have drawn.
    const DEFAULT_CATEGORICAL_PALETTE = 'domainator';

    // Derived, not user data: never written into a session, and never reported as
    // "custom colors" by the picker. Keyed like state.customPalettes (by
    // paletteKey, so a numeric column flipped to discrete gets its own entry) and
    // dropped by rebuildMetadataCaches(), since the assignment depends on the whole
    // set of values in the column.
    function defaultCategoricalPalette(columnName) {{
        const key = paletteKey(columnName);
        const cached = state.defaultPalettes.get(key);
        if (cached) {{
            return cached;
        }}
        const scheme = paletteSchemeByName(DEFAULT_CATEGORICAL_PALETTE);
        const colors = {{}};
        if (scheme) {{
            sortPaletteKeys(distinctColumnValues(columnName, Infinity).values.map(entry => entry.key))
                .forEach((valueKey, index) => {{
                    colors[valueKey] = scheme.colors[index % scheme.colors.length];
                }});
        }}
        // nullColor is deliberately left unset: what color an empty cell gets is a
        // separate choice from the categorical assignment, and categoricalColor's
        // own fallback is already the color the rest of the viewer uses for one.
        const palette = {{
            type: 'categorical',
            colors,
            nullColor: null,
            scheme: scheme ? scheme.name : null,
        }};
        state.defaultPalettes.set(key, palette);
        return palette;
    }}

    // `palette` is a categorical palette ({{colors, nullColor}}) -- normally the
    // column's own, from customPalette(), which now always supplies one for a
    // discrete column. The hue hashed from the value's text is the last resort for
    // a value the palette says nothing about: one typed in after a color table was
    // loaded, say, or after the palette was hand-edited.
    function categoricalColor(value, palette) {{
        if (value === null || value === undefined || value === '') {{
            return (palette && palette.nullColor) || '#b3a89d';
        }}
        const text = String(value);
        if (palette && palette.colors && palette.colors[text]) {{
            return palette.colors[text];
        }}
        let hash = 0;
        for (let i = 0; i < text.length; i++) {{
            hash = ((hash << 5) - hash) + text.charCodeAt(i);
            hash |= 0;
        }}
        const hue = Math.abs(hash) % 360;
        return 'hsl(' + hue + ' 58% 54%)';
    }}

    // `palette` is an optional numeric custom palette ({{stops, nullColor}}).
    // With low+high set, interpolate in RGB (low->mid->high split at 0.5 when mid is
    // present); otherwise fall back to the default cyan->orange hsl gradient.
    function numericColor(value, minValue, maxValue, palette) {{
        if (value === null || value === undefined || Number.isNaN(value)) {{
            return (palette && palette.nullColor) || '#b3a89d';
        }}
        // A custom palette is a list of {{value, color}} stops, kept sorted by whoever
        // stored it. Interpolate between the pair that brackets the value and hold the
        // end colors flat outside the ends -- so two stops are a plain ramp, three
        // reproduce a low/mid/high midpoint, and more shape the ramp arbitrarily.
        const stops = palette && palette.stops;
        if (stops && stops.length >= 2) {{
            if (value <= stops[0].value) {{
                return stops[0].color;
            }}
            const last = stops[stops.length - 1];
            if (value >= last.value) {{
                return last.color;
            }}
            for (let index = 1; index < stops.length; index++) {{
                const upper = stops[index];
                if (value <= upper.value) {{
                    const lower = stops[index - 1];
                    const span = upper.value - lower.value;
                    return span > 0
                        ? lerpHexColor(lower.color, upper.color, (value - lower.value) / span)
                        : upper.color;
                }}
            }}
            return last.color;
        }}
        // No custom palette: the built-in hue ramp across the column's own range.
        const fraction = maxValue <= minValue ? 0.5 : (value - minValue) / (maxValue - minValue);
        const clamped = Math.max(0, Math.min(1, fraction));
        const hue = 200 - (160 * clamped);
        const light = 72 - (24 * clamped);
        return 'hsl(' + hue + ' 72% ' + light + '%)';
    }}

    function colorInfo(columnName) {{
        return state.metadataColorInfoByName.get(columnName) || null;
    }}

    // True when `columnName` holds numbers, and so can be colored either as a gradient
    // or as discrete categories.
    function columnIsNumericType(columnName) {{
        const info = colorInfo(columnName);
        return !!info && info.baseType === 'numeric';
    }}

    function columnIsCategoricalNumeric(columnName) {{
        return columnIsNumericType(columnName) && state.categoricalColumns.has(columnName);
    }}

    // Switch a numeric column between gradient and discrete coloring. No-op for columns
    // that are not numeric (those are always categorical).
    function setColumnCategorical(columnName, categorical) {{
        const info = colorInfo(columnName);
        if (!info || info.baseType !== 'numeric') {{
            return;
        }}
        if (categorical) {{
            state.categoricalColumns.add(columnName);
        }} else {{
            state.categoricalColumns.delete(columnName);
        }}
        info.type = categorical ? 'categorical' : 'numeric';
    }}

    // Custom palettes are stored per column AND per coloring mode, so flipping a numeric
    // column to categorical and back does not discard the gradient the user set up.
    function paletteKey(columnName) {{
        return columnIsCategoricalNumeric(columnName) ? columnName + '\\u0000categorical' : columnName;
    }}

    function customPalette(columnName) {{
        const stored = state.customPalettes[paletteKey(columnName)];
        if (stored) {{
            return stored;
        }}
        // A discrete column with nothing stored is not uncolored -- it is on the
        // default named palette. Only gradient columns fall through to null, where
        // numericColor's own built-in ramp takes over.
        const info = colorInfo(columnName);
        return info && info.type === 'categorical' ? defaultCategoricalPalette(columnName) : null;
    }}

    function nodeColor(nodeIndex) {{
        const columnName = currentColorField();
        const info = colorInfo(columnName);
        const value = metadataValue(nodeIndex, columnName);
        if (!info) {{
            return '#d88f3d';
        }}
        const palette = customPalette(columnName);
        if (info.type === 'numeric') {{
            return numericColor(value, info.min, info.max, palette);
        }}
        return categoricalColor(value, palette);
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

    // Contract chains of pass-through clusters into a single edge.
    //
    // Leaf pruning has already dropped every below-minimum cluster that merely hangs
    // off the graph, so the small clusters still standing are the ones that bridge. A
    // run of them joined end to end is a spindle: it takes up most of the canvas and
    // says nothing except "these two big clusters are related, weakly". Replacing the
    // run with a single edge carrying the run's *weakest* link keeps exactly that
    // statement -- the weakest link is all a path can support, the same
    // min-over-the-path rule the threshold slider itself applies -- and hides the
    // spindle, which is what lets the relationships between the big clusters read.
    //
    // The cluster graph is a contraction of the MST and therefore a forest, and that is
    // what makes this safe: a chain cannot close on itself, contracting a degree-2
    // cluster can neither create a parallel edge nor change any surviving cluster's
    // degree (so no new leaves appear and one pass is enough). The guards below say
    // what happens if a hand-built bundle breaks that assumption.
    //
    // `visibleSet` is narrowed in place: a contracted cluster is no longer drawn, so it
    // counts as hidden by the minimum cluster size like any other.
    function collapseLongPathLinks(links, visibleSet, minClusterSize) {{
        const hierarchyNodes = state.bundle.graph.hierarchy.nodes;
        const adjacency = new Map();
        links.forEach(link => {{
            if (!visibleSet.has(link.sourceId) || !visibleSet.has(link.targetId)) {{ return; }}
            if (!adjacency.has(link.sourceId)) {{ adjacency.set(link.sourceId, []); }}
            if (!adjacency.has(link.targetId)) {{ adjacency.set(link.targetId, []); }}
            adjacency.get(link.sourceId).push({{other: link.targetId, weight: link.weight}});
            adjacency.get(link.targetId).push({{other: link.sourceId, weight: link.weight}});
        }});

        const isPassThrough = componentId => visibleSet.has(componentId)
            && (adjacency.get(componentId) || []).length === 2
            && hierarchyNodes[componentId].size < minClusterSize;

        // Walked in ascending cluster id so the result does not depend on link order.
        const starts = Array.from(adjacency.keys()).filter(isPassThrough).sort((leftId, rightId) => leftId - rightId);
        const consumed = new Set();
        const collapsedLinks = new Map();
        let collapsedPaths = 0;

        starts.forEach(startId => {{
            if (consumed.has(startId)) {{ return; }}
            consumed.add(startId);
            const chain = [startId];
            const ends = [];
            let weakest = Infinity;
            // Two steps out of a degree-2 cluster: walk each way to the first cluster
            // that is not itself pass-through.
            (adjacency.get(startId) || []).forEach(step => {{
                let previousId = startId;
                let currentId = step.other;
                weakest = Math.min(weakest, step.weight);
                while (isPassThrough(currentId) && !consumed.has(currentId)) {{
                    consumed.add(currentId);
                    chain.push(currentId);
                    const next = (adjacency.get(currentId) || []).find(edge => edge.other !== previousId);
                    if (!next) {{ break; }}
                    weakest = Math.min(weakest, next.weight);
                    previousId = currentId;
                    currentId = next.other;
                }}
                ends.push(isPassThrough(currentId) ? null : currentId);
            }});

            chain.forEach(componentId => visibleSet.delete(componentId));
            const leftEnd = ends[0];
            const rightEnd = ends[1];
            if (leftEnd === null || rightEnd === null || leftEnd === rightEnd) {{
                // Only reachable if the forest assumption above does not hold. The chain
                // is gone either way; drawing a self-loop or a dangling edge is worse
                // than drawing nothing.
                return;
            }}
            collapsedPaths += 1;
            const sourceId = Math.min(leftEnd, rightEnd);
            const targetId = Math.max(leftEnd, rightEnd);
            const key = sourceId + ':' + targetId;
            const existing = collapsedLinks.get(key);
            if (existing) {{
                existing.weight = Math.min(existing.weight, weakest);
                existing.collapsed += chain.length;
                return;
            }}
            collapsedLinks.set(key, {{sourceId, targetId, weight: weakest, collapsed: chain.length}});
        }});

        const kept = links.filter(link => visibleSet.has(link.sourceId) && visibleSet.has(link.targetId));
        return {{
            links: kept.concat(Array.from(collapsedLinks.values())),
            collapsedPaths,
            collapsedClusters: consumed.size,
        }};
    }}

    function mstLinksForActiveClusters(activeClusterIds, minClusterSize, leafPruningOnly, collapseLongPaths) {{
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
                collapsedPaths: 0,
                collapsedClusters: 0,
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
                collapsedPaths: 0,
                collapsedClusters: 0,
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

        // Runs after pruning and narrows visibleSet further, so visibleIds and
        // hiddenNodes below already account for the contracted clusters.
        const collapse = collapseLongPaths
            ? collapseLongPathLinks(allLinks, visibleSet, minClusterSize)
            : null;

        const visibleIds = activeClusterIds.filter(componentId => visibleSet.has(componentId));
        const hiddenNodes = activeClusterIds.reduce((sum, componentId) => {{
            if (visibleSet.has(componentId)) {{
                return sum;
            }}
            return sum + state.bundle.graph.hierarchy.nodes[componentId].size;
        }}, 0);

        return {{
            visibleIds,
            links: collapse
                ? collapse.links
                : allLinks.filter(link => visibleSet.has(link.sourceId) && visibleSet.has(link.targetId)),
            hiddenNodes,
            collapsedPaths: collapse ? collapse.collapsedPaths : 0,
            collapsedClusters: collapse ? collapse.collapsedClusters : 0,
        }};
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
                ? packingRadius * 0.82
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

    // Every treemap node is drawn at this fixed world size on one global fixed-pitch lattice, so
    // a node looks identical everywhere and its position never changes with the threshold.
    // TREEMAP_CELL = node square + surrounding padding; cluster boundaries are drawn on the
    // padding between cells so they never overlap a node.
    const TREEMAP_NODE = 13;   // node square side (world units)
    const TREEMAP_CELL = 18;   // lattice pitch = node + padding (world units)
    const TREEMAP_ORIGIN = 58; // world offset of lattice cell (0,0)

    function sign2(value) {{
        return (value > 0) - (value < 0);
    }}

    // Generalized Hilbert ("gilbert") curve: enumerate every cell of a width x height grid in a
    // Hilbert-like order that works for non-power-of-two, non-square rectangles. Contiguous runs
    // of the returned order are spatially compact and connected, which is what makes a contiguous
    // dendrogram leaf range render as a compact blob with a single staircase boundary. Port of
    // Jakub Cerveny's gilbert2d. Recursion depth is O(log(width*height)).
    function gilbertCurve(width, height) {{
        const out = [];
        const generate = (x, y, ax, ay, bx, by) => {{
            const w = Math.abs(ax + ay);
            const h = Math.abs(bx + by);
            const dax = sign2(ax), day = sign2(ay);
            const dbx = sign2(bx), dby = sign2(by);
            if (h === 1) {{
                for (let i = 0; i < w; i++) {{ out.push([x, y]); x += dax; y += day; }}
                return;
            }}
            if (w === 1) {{
                for (let i = 0; i < h; i++) {{ out.push([x, y]); x += dbx; y += dby; }}
                return;
            }}
            let ax2 = Math.floor(ax / 2), ay2 = Math.floor(ay / 2);
            let bx2 = Math.floor(bx / 2), by2 = Math.floor(by / 2);
            const w2 = Math.abs(ax2 + ay2);
            const h2 = Math.abs(bx2 + by2);
            if (2 * w > 3 * h) {{
                if ((w2 % 2) && w > 2) {{ ax2 += dax; ay2 += day; }}
                generate(x, y, ax2, ay2, bx, by);
                generate(x + ax2, y + ay2, ax - ax2, ay - ay2, bx, by);
            }} else {{
                if ((h2 % 2) && h > 2) {{ bx2 += dbx; by2 += dby; }}
                generate(x, y, bx2, by2, ax2, ay2);
                generate(x + bx2, y + by2, ax, ay, bx - bx2, by - by2);
                generate(x + (ax - dax) + (bx2 - dbx), y + (ay - day) + (by2 - dby),
                    -bx2, -by2, -(ax - ax2), -(ay - ay2));
            }}
        }};
        if (width >= height) {{
            generate(0, 0, width, 0, 0, height);
        }} else {{
            generate(0, 0, 0, height, width, 0);
        }}
        return out;
    }}

    // Build (and cache) the global node lattice for the whole network. Every leaf (node) is
    // assigned a fixed cell along a generalized Hilbert (gilbert) curve, whose locality means any
    // contiguous leaf-order range -- i.e. any cluster at any threshold -- maps to a compact,
    // near-square region. The arrangement is threshold-independent: changing the threshold only
    // regroups these fixed cells. Cached by node count.
    function ensureLatticeGlobal() {{
        const leafOrder = state.bundle.graph.hierarchy.leaf_order;
        const nodeCount = leafOrder.length;
        if (state.latticeGlobal && state.latticeGlobal.nodeCount === nodeCount) {{
            return state.latticeGlobal;
        }}
        const aspect = clusterCanvas.width / Math.max(clusterCanvas.height, 1);
        const width = Math.max(1, Math.round(Math.sqrt(nodeCount * aspect)) || 1);
        const height = Math.max(1, Math.ceil(nodeCount / width));
        const curve = gilbertCurve(width, height);
        const colOf = new Int32Array(nodeCount);
        const rowOf = new Int32Array(nodeCount);
        const cellNode = new Int32Array(width * height).fill(-1);
        for (let position = 0; position < nodeCount; position++) {{
            const cell = curve[position];
            const nodeIndex = leafOrder[position];
            colOf[nodeIndex] = cell[0];
            rowOf[nodeIndex] = cell[1];
            cellNode[(cell[1] * width) + cell[0]] = nodeIndex;
        }}
        state.latticeGlobal = {{nodeCount, width, height, colOf, rowOf, cellNode}};
        return state.latticeGlobal;
    }}

    function latticeCellWorld(column, row) {{
        return {{
            x: TREEMAP_ORIGIN + ((column + 0.5) * TREEMAP_CELL),
            y: TREEMAP_ORIGIN + ((row + 0.5) * TREEMAP_CELL),
        }};
    }}

    // Node index under a world-space point, or -1. O(1) via the fixed lattice. The whole cell
    // (node square + its surrounding padding) counts as a hit, not just the drawn square: each
    // cell's padding splits evenly to its neighbors, so a click in the gap between two nodes
    // snaps to the nearer one. This keeps clicks between nodes selecting their cluster instead of
    // falling through. Empty cells (no node) still miss.
    function latticeNodeAtWorld(worldX, worldY) {{
        const lattice = state.latticeGlobal;
        if (!lattice) {{ return -1; }}
        const column = Math.floor((worldX - TREEMAP_ORIGIN) / TREEMAP_CELL);
        const row = Math.floor((worldY - TREEMAP_ORIGIN) / TREEMAP_CELL);
        if (column < 0 || row < 0 || column >= lattice.width || row >= lattice.height) {{ return -1; }}
        return lattice.cellNode[(row * lattice.width) + column];
    }}

    // Fixed-lattice "treemap". Node positions come from the global node lattice (built once and
    // never moving with the threshold); this only turns the active clusters into lightweight items
    // carrying a bounding box (for culling/fit) and the leaf range (for boundary tracing and
    // hit-testing). Every active cluster is emitted regardless of minimum cluster size -- min size
    // only governs which boundaries are drawn (see renderClusterView), never whether nodes shown.
    function treemapClusterLayout(activeClusterIds) {{
        if (activeClusterIds.length === 0) {{
            clusterCanvas.height = 760;
            return [];
        }}
        clusterCanvas.height = Math.max(760, Math.min(1180, Math.round(window.innerHeight * 0.8)));
        const nodes = state.bundle.graph.hierarchy.nodes;
        const leafOrder = state.bundle.graph.hierarchy.leaf_order;
        const lattice = ensureLatticeGlobal();

        const orderedIds = [...activeClusterIds].sort(
            (a, b) => nodes[a].leaf_start - nodes[b].leaf_start || a - b);
        return orderedIds.map(id => {{
            const node = nodes[id];
            const start = node.leaf_start;
            const end = start + node.leaf_count;
            let minCol = Infinity, minRow = Infinity, maxCol = -Infinity, maxRow = -Infinity;
            for (let position = start; position < end; position++) {{
                const nodeIndex = leafOrder[position];
                const column = lattice.colOf[nodeIndex];
                const row = lattice.rowOf[nodeIndex];
                if (column < minCol) {{ minCol = column; }}
                if (column > maxCol) {{ maxCol = column; }}
                if (row < minRow) {{ minRow = row; }}
                if (row > maxRow) {{ maxRow = row; }}
            }}
            const x0 = TREEMAP_ORIGIN + (minCol * TREEMAP_CELL);
            const y0 = TREEMAP_ORIGIN + (minRow * TREEMAP_CELL);
            const x1 = TREEMAP_ORIGIN + ((maxCol + 1) * TREEMAP_CELL);
            const y1 = TREEMAP_ORIGIN + ((maxRow + 1) * TREEMAP_CELL);
            return {{
                componentId: id,
                shape: 'lattice',
                leafStart: start,
                leafCount: node.leaf_count,
                x0, y0, x1, y1,
                x: (x0 + x1) / 2,
                y: (y0 + y1) / 2,
                radius: Math.hypot((x1 - x0) / 2, (y1 - y0) / 2),
            }};
        }});
    }}

    function forceDirectedForestLayout(visibleIds, links, hierarchyNodes, sortBySizeEnabled) {{
        if (visibleIds.length === 0) {{
            clusterCanvas.height = 760;
            return [];
        }}
        clusterCanvas.height = Math.max(760, Math.min(1180, Math.round(window.innerHeight * 0.8)));
        const {{adjacency, components}} = clusterGraphComponents(visibleIds, links, hierarchyNodes, sortBySizeEnabled);
        // Same tuning the Worker uses (forestForceOptions lives in the shared layout core),
        // so the fallback path cannot drift into drawing a different layout.
        const componentLayouts = components.map(componentIds => simulateComponentLayout(
            componentIds, adjacency, forestForceOptions(componentIds.length), hierarchyNodes));
        return packLayouts(componentLayouts, FOREST_PACK_OPTIONS);
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
            if (algorithm === 'treemap') {{
                return treemapClusterLayout(visibleIds);
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

    // Axis titles for the split chart. A `product` merge impact is a product of two
    // component sizes, not a count of nodes, so it must not be labeled as one.
    function splitAxisLabels() {{
        const labels = {split_axis_labels_js};
        const metric = state.bundle && state.bundle.graph
            ? state.bundle.graph.merge_impact_metric
            : null;
        return labels[metric] || labels['{MERGE_IMPACT_MIN_CHILD}'];
    }}

    // ---- Split chart ----
    //
    // One model, three consumers. splitChartLayout() computes the axis box, the scales,
    // every tick value and every mark's position; drawSplitChart() paints it to a canvas,
    // buildSplitChartSVG() emits the same thing as SVG, and the hover handlers hit-test
    // against the positions the paint pass recorded. buildClusterViewSVG's header warns
    // about what happens when two renderers each compute their own geometry -- this is
    // that warning taken seriously, so the exported chart cannot drift from the drawn one.

    const SPLIT_CHART_COLORS = {{
        axis: 'rgba(92,106,112,0.35)',
        tick: 'rgba(92,106,112,0.4)',
        text: '#5c6a70',
        stem: '#dd9687',
        bead: '#c8553d',
        movingSum: '#2f6f8f',
        marker: '#1e2a2f',
        markerDot: '#e29b4b',
    }};
    const SPLIT_CHART_TITLE_FONT = 12;
    const SPLIT_CHART_TICK_FONT = 11;
    const SPLIT_CHART_MESSAGE_FONT = 18;
    const SPLIT_CHART_TICK_LENGTH = 6;
    const SPLIT_CHART_BEAD_DRAW_RADIUS = 4;

    // Decimals needed to print `step` exactly, so a tick label names the value the tick
    // sits at. Without this the axis lies by rounding: a tick at 0.8736 printed as "0.87"
    // reads as misplaced next to a lollipop whose own readout also says 0.87.
    function decimalsForTickStep(step) {{
        const magnitude = Math.abs(step);
        if (!Number.isFinite(magnitude) || magnitude === 0) {{ return 2; }}
        for (let decimals = 0; decimals <= 8; decimals++) {{
            if (Math.abs(Number(magnitude.toFixed(decimals)) - magnitude) <= magnitude * 1e-9) {{
                return decimals;
            }}
        }}
        return 8;
    }}

    // Tick values on round numbers, with the precision their step requires. Ticks used to
    // sit at fixed fractions of the domain, which put them at arbitrary values that
    // formatValue then rounded to two decimals -- both halves of the misalignment.
    // `options.integer` marks an axis of counts, where a tick at 2.5 nodes names nothing.
    function splitAxisTicks(min, max, targetCount, options = {{}}) {{
        let values = chartAxisTicks(min, max, targetCount)
            .filter(value => value >= min - 1e-12 && value <= max + 1e-12);
        if (options.integer && (max - min) < targetCount) {{
            // Too narrow for a nice fractional step to mean anything: one tick per integer.
            values = [];
            for (let value = Math.ceil(min); value <= Math.floor(max); value++) {{
                values.push(value);
            }}
        }}
        if (values.length === 0) {{ values = [min]; }}
        const step = values.length > 1 ? values[1] - values[0] : (Math.abs(values[0]) || 1);
        const decimals = decimalsForTickStep(step);
        return values.map(value => ({{
            value,
            // toLocaleString rather than toFixed so a 50,000-node axis keeps its
            // separators, the way formatValue renders integers elsewhere.
            label: value.toLocaleString(undefined, {{
                minimumFractionDigits: decimals,
                maximumFractionDigits: decimals,
            }}),
        }}));
    }}

    function splitChartLayout(viewWidth, viewHeight) {{
        const margin = {{top: 26, right: 80, bottom: 62, left: 76}};
        const width = viewWidth - margin.left - margin.right;
        const height = viewHeight - margin.top - margin.bottom;
        const axisLabels = splitAxisLabels();
        const frame = {{
            viewWidth,
            viewHeight,
            margin,
            width,
            height,
            titles: {{x: 'Threshold', y: axisLabels.largest, y2: axisLabels.movingSum}},
        }};

        if (!state.bundle) {{
            return Object.assign(frame, {{
                message: 'Load a bundle to render split events.',
                axes: false,
            }});
        }}

        const events = state.series ? state.series.selectedRows : [];
        if (events.length === 0) {{
            return Object.assign(frame, {{
                message: 'No split events in this bundle.',
                axes: true,
            }});
        }}

        // Stems show the LARGEST single merge at each threshold, not the tie group's sum:
        // one big cluster splitting off and a swarm of singletons are opposite stories that
        // the sum renders identically.
        const maxImpact = Math.max(...events.map(event => event.largest_merge), 1);
        const movingSum = state.series.movingSum || {{x: [], y: []}};
        const movingSumX = movingSum.x || [];
        const movingSumY = movingSum.y || [];
        // The axis spans the whole series, never just what is selected -- the selection
        // is a function of the window, so an axis derived from it would chase its own
        // tail as the user zoomed.
        const dataRange = splitSeriesDataRange(state.series);
        const dataMin = dataRange.min;
        const dataMax = dataRange.max;
        // The slice of the series the chart is showing: all of it unless it has been
        // zoomed. Every scale, tick, mark and marker below is derived from this window,
        // so it is the single place the zoom can be got wrong -- and the single place
        // it is clamped back into the data.
        const view = splitChartVisibleWindow(dataMin, dataMax);
        const minThreshold = view.min;
        const maxThreshold = view.max;
        const thresholdSpan = Math.max(1e-9, maxThreshold - minThreshold || 1);
        const maxMovingSum = Math.max(...movingSumY, 1);

        const xFor = value => margin.left + ((value - minThreshold) / thresholdSpan) * width;
        const yFor = value => margin.top + height - (value / maxImpact) * height;
        const y2For = value => margin.top + height - (value / maxMovingSum) * height;

        // Marks outside the window are dropped rather than drawn and clipped. The
        // painters clip the plot box anyway, so this is about the hover: "the nearest
        // event" must not be allowed to mean one the chart is not showing.
        const markSlack = (thresholdSpan * SPLIT_CHART_BEAD_DRAW_RADIUS) / Math.max(1, width);
        const marks = events.filter(event =>
            event.threshold_value >= minThreshold - markSlack &&
            event.threshold_value <= maxThreshold + markSlack
        ).map(event => ({{
            event,
            x: xFor(event.threshold_value),
            stemY: yFor(event.largest_merge),
            // JSON object keys are always strings, so this Number() is required.
            beads: Object.keys(event.merge_size_counts || {{}}).map(size => ({{
                size: Number(size),
                count: event.merge_size_counts[size],
                y: yFor(Number(size)),
            }})),
        }}));

        // The trace is a step function -- a moving sum over discrete events genuinely is
        // one, and smoothing misrepresents where the events sit -- so the vertices are
        // built once here and both painters stroke the same polyline.
        const movingSumPoints = [];
        for (let index = 0; index < movingSumX.length; index++) {{
            const x = xFor(movingSumX[index]);
            if (index > 0) {{
                movingSumPoints.push({{x, y: y2For(movingSumY[index - 1])}});
            }}
            movingSumPoints.push({{x, y: y2For(movingSumY[index])}});
        }}

        const stop = currentSliderStop();
        // Clamped into the plotted range first. The infinity stop sits above every score
        // in the series and the floor stop below every one, so both belong at an end of
        // the axis rather than off it -- the floor stop's marker used to be drawn outside
        // the plot box, to the left of the y-axis. Clamping first also means a zoomed
        // window drops them like any other out-of-window threshold.
        const markerValue = stop === null ? null : Math.min(Math.max(
            stop.threshold_value === null ? dataMax : stop.threshold_value, dataMin), dataMax);

        return Object.assign(frame, {{
            message: null,
            axes: true,
            minThreshold,
            maxThreshold,
            thresholdSpan,
            // The full extent of the series, which the zoom gestures clamp against.
            dataMin,
            dataMax,
            zoomed: view.zoomed,
            maxImpact,
            maxMovingSum,
            xFor,
            yFor,
            y2For,
            // 8 rather than 6: chartAxisTicks' 1/2/2.5/5/10 ladder rounds a raw step that
            // lands just above a power of ten all the way up to the next rung (0.1005
            // becomes 0.2), halving the tick count on exactly the axes -- similarity scores
            // spanning a little over half a decade -- this chart usually draws.
            xTicks: splitAxisTicks(minThreshold, maxThreshold, 8)
                .map(tick => Object.assign({{x: xFor(tick.value)}}, tick)),
            // Both vertical axes count nodes (or, under --merge_impact_metric product, a
            // product of counts); either way a fractional tick is not a quantity.
            yTicks: splitAxisTicks(0, maxImpact, 6, {{integer: true}})
                .map(tick => Object.assign({{y: yFor(tick.value)}}, tick)),
            y2Ticks: splitAxisTicks(0, maxMovingSum, 6, {{integer: true}})
                .map(tick => Object.assign({{y: y2For(tick.value)}}, tick)),
            marks,
            movingSum: {{window: movingSum.window || 0, x: movingSumX, y: movingSumY, points: movingSumPoints}},
            // Dropped once the window no longer contains it: a marker pinned to the
            // edge of a zoomed axis would claim the threshold is sitting there.
            markerX: markerValue === null || markerValue < minThreshold || markerValue > maxThreshold
                ? null
                : xFor(markerValue),
        }});
    }}

    function drawSplitChart(ctx = splitContext, viewWidth = splitCanvas.width,
                            viewHeight = splitCanvas.height, options = {{}}) {{
        // Aliased so the body can paint into any context -- e.g. a scaled offscreen canvas
        // for the high-resolution PNG export -- without rewriting every draw call.
        const context = ctx;
        const onScreen = options.onScreen !== false;
        context.clearRect(0, 0, viewWidth, viewHeight);
        context.fillStyle = '#ffffff';
        context.fillRect(0, 0, viewWidth, viewHeight);

        const layout = splitChartLayout(viewWidth, viewHeight);
        const margin = layout.margin;
        const width = layout.width;
        const height = layout.height;

        if (layout.axes) {{
            context.strokeStyle = SPLIT_CHART_COLORS.axis;
            context.lineWidth = 1;
            context.beginPath();
            context.moveTo(margin.left, margin.top + height);
            context.lineTo(margin.left + width, margin.top + height);
            context.moveTo(margin.left, margin.top);
            context.lineTo(margin.left, margin.top + height);
            context.stroke();

            context.fillStyle = SPLIT_CHART_COLORS.text;
            context.font = SPLIT_CHART_TITLE_FONT + 'px Georgia';
            context.textAlign = 'center';
            context.textBaseline = 'alphabetic';
            context.fillText(layout.titles.x, margin.left + (width / 2), viewHeight - 14);
            context.save();
            context.translate(20, margin.top + (height / 2));
            context.rotate(-Math.PI / 2);
            context.fillText(layout.titles.y, 0, 0);
            context.restore();
            context.save();
            context.translate(viewWidth - 18, margin.top + (height / 2));
            context.rotate(-Math.PI / 2);
            context.fillStyle = SPLIT_CHART_COLORS.movingSum;
            context.fillText(layout.titles.y2, 0, 0);
            context.restore();
        }}

        if (layout.message) {{
            if (onScreen) {{
                recordSplitChartGeometry(null);
                hideSplitChartTip();
            }}
            context.fillStyle = SPLIT_CHART_COLORS.text;
            context.font = SPLIT_CHART_MESSAGE_FONT + 'px Georgia';
            // Set explicitly: the context is reused across calls, so inheriting whatever
            // alignment the last paint left behind once centered this message on x=30 and
            // ran half of it off the left edge.
            context.textAlign = 'left';
            context.textBaseline = 'alphabetic';
            context.fillText(layout.message, 30, 50);
            return;
        }}

        context.fillStyle = SPLIT_CHART_COLORS.text;
        context.font = SPLIT_CHART_TICK_FONT + 'px Georgia';
        context.textBaseline = 'top';
        context.textAlign = 'center';
        context.strokeStyle = SPLIT_CHART_COLORS.tick;
        context.lineWidth = 1;
        layout.xTicks.forEach(tick => {{
            context.beginPath();
            context.moveTo(tick.x, margin.top + height);
            context.lineTo(tick.x, margin.top + height + SPLIT_CHART_TICK_LENGTH);
            context.stroke();
            context.fillText(tick.label, tick.x, margin.top + height + SPLIT_CHART_TICK_LENGTH + 4);
        }});

        context.textAlign = 'right';
        context.textBaseline = 'middle';
        layout.yTicks.forEach(tick => {{
            context.beginPath();
            context.moveTo(margin.left - SPLIT_CHART_TICK_LENGTH, tick.y);
            context.lineTo(margin.left, tick.y);
            context.stroke();
            context.fillText(tick.label, margin.left - SPLIT_CHART_TICK_LENGTH - 6, tick.y);
        }});

        // Right axis (moving sum scale).
        const rightAxisX = margin.left + width;
        context.strokeStyle = SPLIT_CHART_COLORS.axis;
        context.beginPath();
        context.moveTo(rightAxisX, margin.top);
        context.lineTo(rightAxisX, margin.top + height);
        context.stroke();
        context.fillStyle = SPLIT_CHART_COLORS.movingSum;
        context.strokeStyle = SPLIT_CHART_COLORS.movingSum;
        context.textAlign = 'left';
        layout.y2Ticks.forEach(tick => {{
            context.beginPath();
            context.moveTo(rightAxisX, tick.y);
            context.lineTo(rightAxisX + SPLIT_CHART_TICK_LENGTH, tick.y);
            context.stroke();
            context.fillText(tick.label, rightAxisX + SPLIT_CHART_TICK_LENGTH + 4, tick.y);
        }});

        // Clipped to the plot box for as long as the data marks are being painted: a bead
        // is 4px wide, so one sitting against the edge of a zoomed window would otherwise
        // spill over the axis and into the tick labels.
        context.save();
        context.beginPath();
        context.rect(margin.left, margin.top, width, height);
        context.clip();

        // Stem to the largest single split, then one bead per distinct merge size. Beads are
        // drawn unoutlined: on a stem carrying several close beads the outlines would merge
        // into a band that erases the stem. Bead and stem differ by shade, not by outline.
        layout.marks.forEach(mark => {{
            context.strokeStyle = SPLIT_CHART_COLORS.stem;
            context.lineWidth = 1.5;
            context.beginPath();
            context.moveTo(mark.x, margin.top + height);
            context.lineTo(mark.x, mark.stemY);
            context.stroke();
            context.fillStyle = SPLIT_CHART_COLORS.bead;
            mark.beads.forEach(bead => {{
                context.beginPath();
                context.arc(mark.x, bead.y, SPLIT_CHART_BEAD_DRAW_RADIUS, 0, Math.PI * 2);
                context.fill();
            }});
        }});

        if (layout.movingSum.points.length > 0) {{
            context.strokeStyle = SPLIT_CHART_COLORS.movingSum;
            context.lineWidth = 1.75;
            context.beginPath();
            layout.movingSum.points.forEach((point, index) => {{
                if (index === 0) {{
                    context.moveTo(point.x, point.y);
                }} else {{
                    context.lineTo(point.x, point.y);
                }}
            }});
            context.stroke();
        }}
        context.restore();

        // The dashed line and its dot mark where the slider is sitting right now. That is
        // UI state, not something the chart measures, so neither export draws it --
        // buildSplitChartSVG omits it for the same reason.
        if (onScreen && layout.markerX !== null) {{
            context.strokeStyle = SPLIT_CHART_COLORS.marker;
            context.lineWidth = 1;
            context.setLineDash([6, 6]);
            context.beginPath();
            context.moveTo(layout.markerX, margin.top);
            context.lineTo(layout.markerX, margin.top + height);
            context.stroke();
            context.setLineDash([]);
            context.fillStyle = SPLIT_CHART_COLORS.markerDot;
            context.beginPath();
            context.arc(layout.markerX, margin.top + height + 10, 5, 0, Math.PI * 2);
            context.fill();
        }}

        // Recorded from the layout the paint pass just used, so a mark the hover finds is a
        // mark that was actually drawn. See setupSplitChartHover / splitChartHitAt.
        if (onScreen) {{
            recordSplitChartGeometry({{
                plotLeft: margin.left,
                plotTop: margin.top,
                plotWidth: width,
                plotHeight: height,
                minThreshold: layout.minThreshold,
                thresholdSpan: layout.thresholdSpan,
                // The gestures zoom and pan against the window the paint pass used and
                // the extent it was clamped into, so they cannot disagree with the
                // chart in front of the user.
                dataMin: layout.dataMin,
                dataMax: layout.dataMax,
                events: layout.marks.map(mark => ({{event: mark.event, x: mark.x, beads: mark.beads}})),
                movingSumWindow: layout.movingSum.window,
                movingSumX: layout.movingSum.x,
                movingSumY: layout.movingSum.y,
            }});
        }}
    }}

    function buildSplitChartSVG() {{
        const viewWidth = splitCanvas.width;
        const viewHeight = splitCanvas.height;
        const fmt = value => Math.round(value * 100) / 100;
        const layout = splitChartLayout(viewWidth, viewHeight);
        const margin = layout.margin;
        const width = layout.width;
        const height = layout.height;
        const parts = [];
        parts.push('<svg xmlns="http://www.w3.org/2000/svg" width="' + viewWidth + '" height="' + viewHeight +
            '" viewBox="0 0 ' + viewWidth + ' ' + viewHeight + '" font-family="Georgia, serif">');
        parts.push('<rect x="0" y="0" width="' + viewWidth + '" height="' + viewHeight + '" fill="#ffffff"/>');

        const rotatedTitle = (text, x, y, color) =>
            '<text transform="translate(' + fmt(x) + ' ' + fmt(y) + ') rotate(-90)" text-anchor="middle"' +
            ' font-size="' + SPLIT_CHART_TITLE_FONT + '" fill="' + color + '">' + escapeXml(text) + '</text>';

        if (layout.axes) {{
            const axis = svgColorParts(SPLIT_CHART_COLORS.axis);
            parts.push('<path d="M' + margin.left + ' ' + margin.top + ' L' + margin.left + ' ' +
                (margin.top + height) + ' L' + (margin.left + width) + ' ' + (margin.top + height) +
                '" fill="none" stroke="' + axis.color + '" stroke-opacity="' + axis.opacity + '" stroke-width="1"/>');
            parts.push('<text x="' + fmt(margin.left + (width / 2)) + '" y="' + (viewHeight - 14) +
                '" text-anchor="middle" font-size="' + SPLIT_CHART_TITLE_FONT + '" fill="' +
                SPLIT_CHART_COLORS.text + '">' + escapeXml(layout.titles.x) + '</text>');
            parts.push(rotatedTitle(layout.titles.y, 20, margin.top + (height / 2), SPLIT_CHART_COLORS.text));
            parts.push(rotatedTitle(layout.titles.y2, viewWidth - 18, margin.top + (height / 2),
                SPLIT_CHART_COLORS.movingSum));
        }}

        if (layout.message) {{
            parts.push('<text x="30" y="50" font-size="' + SPLIT_CHART_MESSAGE_FONT + '" fill="' +
                SPLIT_CHART_COLORS.text + '">' + escapeXml(layout.message) + '</text>');
            parts.push('</svg>');
            return parts.join('\\n');
        }}

        // The same clip the canvas painter applies while drawing the data marks.
        parts.push('<defs><clipPath id="split-plot-clip"><rect x="' + margin.left + '" y="' + margin.top +
            '" width="' + width + '" height="' + height + '"/></clipPath></defs>');

        // Ticks. Values, positions and label text all come from the layout, so the export
        // cannot round or place them differently from the canvas.
        const tickColor = svgColorParts(SPLIT_CHART_COLORS.tick);
        const tickMarks = [];
        const tickLabels = [];
        layout.xTicks.forEach(tick => {{
            tickMarks.push('<path d="M' + fmt(tick.x) + ' ' + (margin.top + height) + ' V' +
                (margin.top + height + SPLIT_CHART_TICK_LENGTH) + '"/>');
            tickLabels.push('<text x="' + fmt(tick.x) + '" y="' +
                (margin.top + height + SPLIT_CHART_TICK_LENGTH + 4) + '" text-anchor="middle"' +
                ' dominant-baseline="hanging">' + escapeXml(tick.label) + '</text>');
        }});
        layout.yTicks.forEach(tick => {{
            tickMarks.push('<path d="M' + (margin.left - SPLIT_CHART_TICK_LENGTH) + ' ' + fmt(tick.y) +
                ' H' + margin.left + '"/>');
            tickLabels.push('<text x="' + (margin.left - SPLIT_CHART_TICK_LENGTH - 6) + '" y="' + fmt(tick.y) +
                '" text-anchor="end" dominant-baseline="central">' + escapeXml(tick.label) + '</text>');
        }});
        parts.push('<g fill="none" stroke="' + tickColor.color + '" stroke-opacity="' + tickColor.opacity +
            '" stroke-width="1">' + tickMarks.join('') + '</g>');
        parts.push('<g font-size="' + SPLIT_CHART_TICK_FONT + '" fill="' + SPLIT_CHART_COLORS.text + '">' +
            tickLabels.join('') + '</g>');

        // Right axis (moving sum scale).
        const rightAxisX = margin.left + width;
        const axisColor = svgColorParts(SPLIT_CHART_COLORS.axis);
        parts.push('<path d="M' + rightAxisX + ' ' + margin.top + ' V' + (margin.top + height) +
            '" fill="none" stroke="' + axisColor.color + '" stroke-opacity="' + axisColor.opacity +
            '" stroke-width="1"/>');
        const rightTickMarks = [];
        const rightTickLabels = [];
        layout.y2Ticks.forEach(tick => {{
            rightTickMarks.push('<path d="M' + rightAxisX + ' ' + fmt(tick.y) + ' H' +
                (rightAxisX + SPLIT_CHART_TICK_LENGTH) + '"/>');
            rightTickLabels.push('<text x="' + (rightAxisX + SPLIT_CHART_TICK_LENGTH + 4) + '" y="' +
                fmt(tick.y) + '" dominant-baseline="central">' + escapeXml(tick.label) + '</text>');
        }});
        parts.push('<g fill="none" stroke="' + SPLIT_CHART_COLORS.movingSum + '" stroke-width="1">' +
            rightTickMarks.join('') + '</g>');
        parts.push('<g font-size="' + SPLIT_CHART_TICK_FONT + '" fill="' + SPLIT_CHART_COLORS.movingSum + '">' +
            rightTickLabels.join('') + '</g>');

        // Stems, then beads, in the canvas's own order so overlaps stack the same way.
        const stems = layout.marks.map(mark =>
            '<path d="M' + fmt(mark.x) + ' ' + (margin.top + height) + ' V' + fmt(mark.stemY) + '"/>');
        parts.push('<g clip-path="url(#split-plot-clip)" fill="none" stroke="' + SPLIT_CHART_COLORS.stem +
            '" stroke-width="1.5">' + stems.join('') + '</g>');
        const beads = [];
        layout.marks.forEach(mark => {{
            mark.beads.forEach(bead => {{
                beads.push('<circle cx="' + fmt(mark.x) + '" cy="' + fmt(bead.y) + '" r="' +
                    SPLIT_CHART_BEAD_DRAW_RADIUS + '"/>');
            }});
        }});
        parts.push('<g clip-path="url(#split-plot-clip)" fill="' + SPLIT_CHART_COLORS.bead + '">' +
            beads.join('') + '</g>');

        if (layout.movingSum.points.length > 0) {{
            const d = layout.movingSum.points
                .map((point, index) => (index === 0 ? 'M' : ' L') + fmt(point.x) + ' ' + fmt(point.y))
                .join('');
            parts.push('<path d="' + d + '" clip-path="url(#split-plot-clip)" fill="none" stroke="' +
                SPLIT_CHART_COLORS.movingSum + '" stroke-width="1.75"/>');
        }}

        // No threshold marker: this builder is only ever an export, and the dashed line
        // and orange dot say where the slider is, which is not a property of the data.
        // See the matching `onScreen` guard in drawSplitChart.

        parts.push('</svg>');
        return parts.join('\\n');
    }}

    function exportSplitChartSVG() {{
        if (!state.bundle) {{
            return;
        }}
        const blob = new Blob([buildSplitChartSVG()], {{type: 'image/svg+xml'}});
        triggerDownload(URL.createObjectURL(blob), exportBaseName() + '_split_events.svg', true);
    }}

    function exportSplitChartPNG() {{
        if (!state.bundle) {{
            return;
        }}
        // Re-rasterized at the chosen density rather than snapshotted and upscaled, and
        // offscreen so the hover geometry the on-screen chart recorded is left alone.
        const scaleFactor = selectedPngScale();
        const target = document.createElement('canvas');
        target.width = Math.round(splitCanvas.width * scaleFactor);
        target.height = Math.round(splitCanvas.height * scaleFactor);
        const targetContext = target.getContext('2d');
        if (!targetContext) {{
            window.alert('Could not export PNG at ' + scaleFactor + '×; try a lower resolution.');
            return;
        }}
        targetContext.scale(scaleFactor, scaleFactor);
        drawSplitChart(targetContext, splitCanvas.width, splitCanvas.height, {{onScreen: false}});
        target.toBlob(blob => {{
            if (!blob) {{
                window.alert('The chart is too large to export as a ' + scaleFactor + '× PNG (' +
                    target.width + '×' + target.height + ' px). Try a lower resolution.');
                return;
            }}
            const suffix = scaleFactor > 1 ? '_split_events@' + scaleFactor + 'x.png' : '_split_events.png';
            triggerDownload(URL.createObjectURL(blob), exportBaseName() + suffix, true);
        }}, 'image/png');
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
            if (item.shape === 'rect' || item.shape === 'lattice') {{
                minX = Math.min(minX, item.x0);
                minY = Math.min(minY, item.y0);
                maxX = Math.max(maxX, item.x1);
                maxY = Math.max(maxY, item.y1);
            }} else {{
                minX = Math.min(minX, item.x - item.radius);
                minY = Math.min(minY, item.y - item.radius);
                maxX = Math.max(maxX, item.x + item.radius);
                maxY = Math.max(maxY, item.y + item.radius);
            }}
        }});
        return {{minX, minY, maxX, maxY}};
    }}

    // Margin left between the fitted content and the canvas edge, shared by
    // "Reset view" and "Focus selection" so both leave the same breathing room.
    const VIEW_FIT_PADDING = 42;

    function fitClusterViewToLayout() {{
        const bounds = clusterLayoutBounds();
        if (!bounds) {{
            state.viewTransform.scale = 1;
            state.viewTransform.offsetX = 0;
            state.viewTransform.offsetY = 0;
            return;
        }}
        const padding = VIEW_FIT_PADDING;
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

    // World-space bounding box of the selected nodes, grown by each node's own radius so
    // that "fits on screen" means the drawn dot fits, not merely its center. Selected nodes
    // sitting in clusters the current view hides (minimum cluster size, threshold) have no
    // position at all and are skipped; returns null when nothing selected is on screen.
    function selectedNodeBounds() {{
        if (!state.bundle || state.selectedNodeIndices.size === 0) {{
            return null;
        }}
        let minX = Infinity;
        let minY = Infinity;
        let maxX = -Infinity;
        let maxY = -Infinity;
        let count = 0;
        state.visibleLayout.forEach(item => {{
            const members = componentMembers(item.componentId);
            // Skip the member layout entirely for clusters holding nothing selected.
            if (!componentSelectionState(members).anySelected) {{
                return;
            }}
            const component = state.bundle.graph.hierarchy.nodes[item.componentId];
            for (const dot of componentMemberLayout(component, item)) {{
                if (!state.selectedNodeIndices.has(dot.memberIndex)) {{
                    continue;
                }}
                minX = Math.min(minX, dot.x - dot.radius);
                minY = Math.min(minY, dot.y - dot.radius);
                maxX = Math.max(maxX, dot.x + dot.radius);
                maxY = Math.max(maxY, dot.y + dot.radius);
                count += 1;
            }}
        }});
        return count === 0 ? null : {{minX, minY, maxX, maxY, count}};
    }}

    // Center the view on the selection, zooming out only when it does not already fit.
    // Holding the zoom whenever it fits is deliberate: stepping between selection presets
    // is the common case, and a zoom that changes on every step is disorienting.
    // Centering uses the bounding box center rather than the centroid of the nodes --
    // with a lopsided selection the centroid can leave outliers off screen even at a
    // scale that fits, which would contradict the "all fit" guarantee.
    function focusSelection() {{
        const bounds = selectedNodeBounds();
        if (!bounds) {{
            return false;
        }}
        const padding = VIEW_FIT_PADDING;
        const width = Math.max(1, bounds.maxX - bounds.minX);
        const height = Math.max(1, bounds.maxY - bounds.minY);
        const usableWidth = Math.max(1, clusterCanvas.width - (padding * 2));
        const usableHeight = Math.max(1, clusterCanvas.height - (padding * 2));
        const fitScale = Math.min(usableWidth / width, usableHeight / height);
        if (fitScale < state.viewTransform.scale) {{
            // Clamped like every other zoom; a selection wider than minScale allows stays
            // partly off screen, which beats silently exceeding the zoom range.
            state.viewTransform.scale = Math.max(
                state.viewTransform.minScale,
                Math.min(state.viewTransform.maxScale, fitScale),
            );
        }}
        const scale = state.viewTransform.scale;
        const centerX = (bounds.minX + bounds.maxX) / 2;
        const centerY = (bounds.minY + bounds.maxY) / 2;
        state.viewTransform.offsetX = (clusterCanvas.width / 2) - (centerX * scale);
        state.viewTransform.offsetY = (clusterCanvas.height / 2) - (centerY * scale);
        return true;
    }}

    function drawBadge(ctx, text, x, y) {{
        const clusterContext = ctx;
        clusterContext.save();
        clusterContext.font = '11px Georgia';
        const textWidth = clusterContext.measureText(text).width;
        const badgeWidth = textWidth + 12;
        const badgeHeight = 18;
        clusterContext.fillStyle = 'rgba(255, 255, 255, 0.92)';
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
        // Trim each end to its own bubble radius so the edge meets each boundary exactly. (Using a
        // shared distance/2 cap would stop the edge at the midpoint — i.e. *inside* a much larger
        // bubble.) If the two bubbles overlap, shrink proportionally so the edge never penetrates.
        let leftOffset = link.left.radius;
        let rightOffset = link.right.radius;
        if (leftOffset + rightOffset > distance) {{
            const k = distance / (leftOffset + rightOffset);
            leftOffset *= k;
            rightOffset *= k;
        }}
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

    // Where an edge's score badge belongs: the halfway point *along the line as drawn*,
    // plus the unit normal of the segment it lands on.
    //
    // This used to be the midpoint of the two bubble centers, which is not on the drawn line.
    // The link is trimmed to each bubble's own radius, so that midpoint is off by
    // (rightRadius - leftRadius) / 2; and in the Tree layout the link is drawn as a
    // three-segment elbow the straight chord does not follow at all. Offsetting along the
    // normal rather than straight up also keeps the badge clear of a vertical elbow riser.
    //
    // `length` is the drawn length, which is what the fit test has to measure -- two large
    // bubbles nearly touching are far apart center to center but have almost no edge showing.
    function linkLabelAnchor(link) {{
        const segments = renderedLinkSegments(link);
        if (segments.length === 0) {{
            return null;
        }}
        const lengths = segments.map(segment =>
            Math.hypot(segment.endX - segment.startX, segment.endY - segment.startY));
        const total = lengths.reduce((sum, value) => sum + value, 0);
        let remaining = total / 2;
        for (let index = 0; index < segments.length; index++) {{
            const length = lengths[index];
            if (remaining > length && index < segments.length - 1) {{
                remaining -= length;
                continue;
            }}
            const segment = segments[index];
            const fraction = length > 1e-9 ? remaining / length : 0;
            const unitX = length > 1e-9 ? (segment.endX - segment.startX) / length : 1;
            const unitY = length > 1e-9 ? (segment.endY - segment.startY) / length : 0;
            return {{
                x: segment.startX + ((segment.endX - segment.startX) * fraction),
                y: segment.startY + ((segment.endY - segment.startY) * fraction),
                // Normal pointing "up" for a left-to-right segment, matching the -8px nudge
                // the badge used to get unconditionally.
                normalX: unitY,
                normalY: -unitX,
                length: total,
            }};
        }}
        return null;
    }}

    // ---- Label level-of-detail ----
    //
    // A label is worth drawing when it fits the mark it labels, which is a property
    // of that mark's size *on screen*. These used to be one global zoom threshold
    // each, which had it backwards: a network of two enormous clusters fits the
    // viewport at a small zoom, so the labels were suppressed while every bubble had
    // hundreds of pixels of room going spare. Keying off the mark makes the rule
    // scale itself, and applying it per mark means a crowded layout still drops only
    // the labels that genuinely do not fit.
    //
    // Shared by renderClusterView and buildClusterViewSVG, which have to agree about
    // which labels exist, and using estimateTextWidth rather than the canvas
    // measurer for the same reason -- the SVG path has no canvas to measure with.
    const CLUSTER_COUNT_LABEL_FONT = 12;
    // Dash pattern for a collapsed path, in screen pixels. Shared by renderClusterView
    // and buildClusterViewSVG so the export draws the dashes the screen drew.
    const LINK_DASH_ON = 6;
    const LINK_DASH_OFF = 4;
    const EDGE_SCORE_LABEL_FONT = 11;
    // drawBadge pads the text by 12px; the extra 8 keeps the badge off the link's
    // endpoints, where it would sit on top of the cluster bubbles it joins.
    const EDGE_SCORE_BADGE_PADDING = 12;
    const EDGE_SCORE_LINK_MARGIN = 8;
    // How far the badge sits off the line it labels, along that segment's normal.
    const EDGE_SCORE_LABEL_OFFSET = 8;
    // Below this a dot is not a mark any more, just a tinted pixel, and labeling it
    // points at nothing the user can see.
    const MIN_LABELED_DOT_SCREEN_RADIUS = 1.5;

    // On-screen extent of a layout item's footprint. Lattice and rect items are boxes,
    // so their own width is the room available -- `radius` is a half-diagonal there and
    // would overstate it for a tall, narrow cluster.
    function itemScreenExtent(item) {{
        const scale = state.viewTransform.scale;
        if (item.shape === 'lattice' || item.shape === 'rect') {{
            return {{width: (item.x1 - item.x0) * scale, height: (item.y1 - item.y0) * scale}};
        }}
        const diameter = item.radius * 2 * scale;
        return {{width: diameter, height: diameter}};
    }}

    function clusterCountLabelFits(text, item) {{
        const extent = itemScreenExtent(item);
        if (extent.height < CLUSTER_COUNT_LABEL_FONT + 4) {{
            return false;
        }}
        return estimateTextWidth(text, CLUSTER_COUNT_LABEL_FONT) + 6 <= extent.width;
    }}

    function edgeScoreLabelFits(text, linkScreenLength) {{
        const badgeWidth = estimateTextWidth(text, EDGE_SCORE_LABEL_FONT) + EDGE_SCORE_BADGE_PADDING;
        return linkScreenLength >= badgeWidth + (EDGE_SCORE_LINK_MARGIN * 2);
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

    // Treemap member layout: every member sits at its fixed cell on the global gilbert lattice,
    // so node positions are identical in every cluster and never move when the threshold changes.
    function componentSquareLayout(component, item) {{
        const members = componentMembers(item.componentId);
        const results = [];
        const lattice = state.latticeGlobal;
        if (members.length === 0 || !lattice) {{ return results; }}
        const dotRadius = TREEMAP_NODE / 2;
        for (let i = 0; i < members.length; i++) {{
            const nodeIndex = members[i];
            const center = latticeCellWorld(lattice.colOf[nodeIndex], lattice.rowOf[nodeIndex]);
            results.push({{
                memberIndex: nodeIndex,
                x: center.x,
                y: center.y,
                radius: dotRadius,
                square: true,
            }});
        }}
        return results;
    }}

    // Shape-aware member placement: lattice tiles read from the global node lattice, everything
    // else keeps the circular phyllotaxis packing. Rendering, SVG export, and hit-testing all
    // route through here so treemap and bubble modes stay in sync.
    function componentMemberLayout(component, item) {{
        if (item.shape === 'lattice') {{
            return componentSquareLayout(component, item);
        }}
        return componentDotLayout(component, item);
    }}

    function componentSelectionState(members, nodeSet = state.selectedNodeIndices) {{
        let selectedCount = 0;
        members.forEach(nodeIndex => {{
            if (nodeSet.has(nodeIndex)) {{
                selectedCount += 1;
            }}
        }});
        return {{
            selectedCount,
            anySelected: selectedCount > 0,
            allSelected: members.length > 0 && selectedCount === members.length,
        }};
    }}

    function renderClusterView(ctx = clusterContext, viewWidth = clusterCanvas.width, viewHeight = clusterCanvas.height, options = {{}}) {{
        // Aliased so the body below can target any context (e.g. a scaled
        // offscreen canvas for high-resolution PNG export) without rewriting
        // every draw call. Defaults render to the on-screen canvas.
        const clusterContext = ctx;
        const clusterCanvas = {{width: viewWidth, height: viewHeight}};
        clusterContext.clearRect(0, 0, clusterCanvas.width, clusterCanvas.height);
        clusterContext.fillStyle = '#ffffff';
        clusterContext.fillRect(0, 0, clusterCanvas.width, clusterCanvas.height);

        if (!state.bundle) {{
            clusterContext.fillStyle = '#5c6a70';
            clusterContext.font = '18px Georgia';
            clusterContext.fillText('Load a bundle to render cluster bubbles.', 30, 50);
            return;
        }}

        const selectedNodeOutlines = [];
        // Hovering a preset slot outlines that preset's nodes/clusters in gray.
        // It is a transient hover cue, so exports opt out via options.preview.
        const previewNodeSet = options.preview === false ? null : presetPreviewNodeSet();
        const previewNodeOutlines = [];
        const previewBoundsPaths = [];
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
            // A collapsed path is not an MST edge between these two clusters -- it is the
            // weakest link along a chain that was contracted away -- so it is dashed
            // rather than drawn as if it were a measured similarity between them.
            clusterContext.setLineDash(link.collapsed
                ? [LINK_DASH_ON / drawScale, LINK_DASH_OFF / drawScale]
                : []);
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
        clusterContext.setLineDash([]);

        // Collect dots batched by color for batch drawing (#3). Circular member dots and
        // square treemap members are batched separately so each can use its own draw call.
        const dotColorBuckets = new Map();
        const squareColorBuckets = new Map();

        const showClusterBounds = renderClusterBoundsEnabled();
        const showNodes = renderNodesEnabled();
        // In the lattice ("treemap") layout nodes are never hidden; the minimum cluster size only
        // decides which cluster boundaries get drawn (small clusters just go un-outlined).
        const minOutlineSize = Math.max(1, Number(document.getElementById('min-cluster-size').value) || 1);
        const lattice = state.latticeGlobal;

        // Per-node ("Label by") labels attach to individual member dots. Collect candidates that are
        // in the viewport during the dot pass; we draw them only when few enough nodes are visible
        // (so labels appear progressively as you zoom in). Overlap is acceptable.
        const dotLabelField = currentLabelField();
        const MAX_LABELED_DOTS = 250;
        // Whether a dot is drawn large enough to be worth pointing at is decided per
        // dot below; MAX_LABELED_DOTS is what keeps a crowded view readable.
        const collectDotLabels = dotLabelField !== '' && showNodes;
        const dotLabelCandidates = [];
        let dotLabelsOverflow = false;

        state.visibleLayout.forEach(item => {{
            // Viewport culling (#4)
            if (item.x + item.radius < worldMinX || item.x - item.radius > worldMaxX ||
                item.y + item.radius < worldMinY || item.y - item.radius > worldMaxY) {{
                return;
            }}

            const component = state.bundle.graph.hierarchy.nodes[item.componentId];
            const members = componentMembers(item.componentId);
            const selectionState = componentSelectionState(members);
            const previewState = previewNodeSet ? componentSelectionState(members, previewNodeSet) : null;

            if (item.shape === 'lattice') {{
                // Lattice clusters are outlined with a staircase along the padding between cells.
                // Nodes are never hidden; a cluster is only outlined when it is big enough (or
                // fully selected). Node squares themselves are drawn in the shared dot pass below.
                const outline = showClusterBounds && component.size >= minOutlineSize;
                if ((outline || selectionState.allSelected) && lattice) {{
                    const width = lattice.width;
                    const cellSet = new Set();
                    for (const nodeIndex of members) {{
                        cellSet.add((lattice.rowOf[nodeIndex] * width) + lattice.colOf[nodeIndex]);
                    }}
                    clusterContext.strokeStyle = selectionState.allSelected ? '#1e2a2f'
                        : (outline ? '#9aa4ad' : 'rgba(30, 42, 47, 0.45)');
                    clusterContext.lineWidth = (selectionState.allSelected ? 2.4 : 1.2) / drawScale;
                    clusterContext.beginPath();
                    for (const nodeIndex of members) {{
                        const column = lattice.colOf[nodeIndex];
                        const row = lattice.rowOf[nodeIndex];
                        const ex0 = TREEMAP_ORIGIN + (column * TREEMAP_CELL);
                        const ey0 = TREEMAP_ORIGIN + (row * TREEMAP_CELL);
                        const ex1 = ex0 + TREEMAP_CELL;
                        const ey1 = ey0 + TREEMAP_CELL;
                        if (!cellSet.has((row * width) + (column - 1))) {{ clusterContext.moveTo(ex0, ey0); clusterContext.lineTo(ex0, ey1); }}
                        if (!cellSet.has((row * width) + (column + 1))) {{ clusterContext.moveTo(ex1, ey0); clusterContext.lineTo(ex1, ey1); }}
                        if (!cellSet.has(((row - 1) * width) + column)) {{ clusterContext.moveTo(ex0, ey0); clusterContext.lineTo(ex1, ey0); }}
                        if (!cellSet.has(((row + 1) * width) + column)) {{ clusterContext.moveTo(ex0, ey1); clusterContext.lineTo(ex1, ey1); }}
                    }}
                    clusterContext.stroke();
                }}
            }} else {{
                if (previewState && previewState.allSelected) {{
                    previewBoundsPaths.push(item);
                }}
                const isRect = item.shape === 'rect';
                const traceBounds = () => {{
                    clusterContext.beginPath();
                    if (isRect) {{
                        clusterContext.rect(item.x0, item.y0, item.x1 - item.x0, item.y1 - item.y0);
                    }} else {{
                        clusterContext.arc(item.x, item.y, item.radius, 0, TAU);
                    }}
                }};
                if (showClusterBounds) {{
                    clusterContext.fillStyle = 'rgba(245, 246, 248, 0.95)';
                    clusterContext.strokeStyle = selectionState.allSelected ? '#1e2a2f' : '#c4cad2';
                    clusterContext.lineWidth = (selectionState.allSelected ? 3 : 1.5) / drawScale;
                    traceBounds();
                    clusterContext.fill();
                    clusterContext.stroke();
                }} else if (selectionState.allSelected) {{
                    // Cluster bounds are hidden, but keep a faint outline so the
                    // selection is still discernible.
                    clusterContext.strokeStyle = 'rgba(30, 42, 47, 0.45)';
                    clusterContext.lineWidth = 1.5 / drawScale;
                    traceBounds();
                    clusterContext.stroke();
                }}
            }}

            // When nodes are hidden we still need to know where the selected
            // ones are so we can hint at them, otherwise skip the dot layout.
            const needsSelectionHint = !selectionState.allSelected && selectionState.anySelected;
            if (!showNodes && !needsSelectionHint) {{
                return;
            }}

            const dotLayout = componentMemberLayout(component, item);
            for (const dot of dotLayout) {{
                if (showNodes) {{
                    // Use pre-computed color cache (#2)
                    const color = state.nodeColorCache[dot.memberIndex] ?? nodeColor(dot.memberIndex);
                    const buckets = dot.square ? squareColorBuckets : dotColorBuckets;
                    let bucket = buckets.get(color);
                    if (bucket === undefined) {{
                        bucket = [];
                        buckets.set(color, bucket);
                    }}
                    bucket.push(dot);
                }}
                if (!selectionState.allSelected && state.selectedNodeIndices.has(dot.memberIndex)) {{
                    selectedNodeOutlines.push({{x: dot.x, y: dot.y, radius: dot.radius, faint: !showNodes}});
                }}
                if (previewNodeSet && previewNodeSet.has(dot.memberIndex)) {{
                    previewNodeOutlines.push({{x: dot.x, y: dot.y, radius: dot.radius}});
                }}
                if (collectDotLabels && !dotLabelsOverflow &&
                    dot.x >= worldMinX && dot.x <= worldMaxX && dot.y >= worldMinY && dot.y <= worldMaxY) {{
                    if (dot.radius * state.viewTransform.scale < MIN_LABELED_DOT_SCREEN_RADIUS) {{
                        // Sub-pixel dot: a label beside it would point at nothing visible.
                    }} else if (dotLabelCandidates.length >= MAX_LABELED_DOTS) {{
                        dotLabelsOverflow = true;
                    }} else {{
                        dotLabelCandidates.push({{x: dot.x, y: dot.y, r: dot.radius, memberIndex: dot.memberIndex}});
                    }}
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

        // Batch draw treemap square members grouped by color.
        squareColorBuckets.forEach((squares, color) => {{
            clusterContext.fillStyle = color;
            clusterContext.beginPath();
            for (const square of squares) {{
                const side = square.radius * 2;
                clusterContext.rect(square.x - square.radius, square.y - square.radius, side, side);
            }}
            clusterContext.fill();
        }});

        clusterContext.restore();

        // Drawn before the selection outlines and at a wider radius, so a node that
        // is both previewed and selected shows a gray ring outside a black one.
        if (previewNodeOutlines.length > 0 || previewBoundsPaths.length > 0) {{
            clusterContext.save();
            clusterContext.strokeStyle = '#9aa4ad';
            clusterContext.lineWidth = 2;
            previewBoundsPaths.forEach(item => {{
                clusterContext.beginPath();
                if (item.shape === 'rect') {{
                    const topLeft = worldToScreenPoint(item.x0, item.y0);
                    const bottomRight = worldToScreenPoint(item.x1, item.y1);
                    clusterContext.rect(topLeft.x, topLeft.y, bottomRight.x - topLeft.x, bottomRight.y - topLeft.y);
                }} else {{
                    const center = worldToScreenPoint(item.x, item.y);
                    clusterContext.arc(center.x, center.y, (item.radius * state.viewTransform.scale) + 3, 0, TAU);
                }}
                clusterContext.stroke();
            }});
            previewNodeOutlines.forEach(outline => {{
                const screenPoint = worldToScreenPoint(outline.x, outline.y);
                const screenRadius = Math.max(4, (outline.radius * state.viewTransform.scale) + 3.4);
                clusterContext.beginPath();
                clusterContext.arc(screenPoint.x, screenPoint.y, screenRadius, 0, TAU);
                clusterContext.stroke();
            }});
            clusterContext.restore();
        }}

        if (selectedNodeOutlines.length > 0) {{
            clusterContext.save();
            selectedNodeOutlines.forEach(outline => {{
                const screenPoint = worldToScreenPoint(outline.x, outline.y);
                const screenRadius = Math.max(3, (outline.radius * state.viewTransform.scale) + 1.6);
                clusterContext.beginPath();
                clusterContext.arc(screenPoint.x, screenPoint.y, screenRadius, 0, TAU);
                if (outline.faint) {{
                    // Nodes are hidden, so leave a faint mark where the selected
                    // node sits instead of outlining a dot that isn't drawn.
                    clusterContext.fillStyle = 'rgba(30, 42, 47, 0.30)';
                    clusterContext.fill();
                }} else {{
                    clusterContext.strokeStyle = '#1e2a2f';
                    clusterContext.lineWidth = 1.8;
                    clusterContext.stroke();
                }}
            }});
            clusterContext.restore();
        }}

        if (showEdgeScoresEnabled()) {{
            state.splitLinks.forEach(link => {{
                const anchor = linkLabelAnchor(link);
                if (anchor === null) {{
                    return;
                }}
                const text = formatValue(link.threshold);
                if (!edgeScoreLabelFits(text, anchor.length * state.viewTransform.scale)) {{
                    return;
                }}
                const point = worldToScreenPoint(anchor.x, anchor.y);
                drawBadge(clusterContext, text,
                    point.x + (anchor.normalX * EDGE_SCORE_LABEL_OFFSET),
                    point.y + (anchor.normalY * EDGE_SCORE_LABEL_OFFSET));
            }});
        }}

        // Node-count labels: one per cluster bubble, centered inside, wherever the
        // bubble is big enough on screen to hold the text.
        if (showNodeCountsEnabled()) {{
            clusterContext.fillStyle = '#5c6a70';
            clusterContext.font = '600 ' + CLUSTER_COUNT_LABEL_FONT + 'px Georgia';
            clusterContext.textAlign = 'center';
            clusterContext.textBaseline = 'middle';
            state.visibleLayout.forEach(item => {{
                const screenPoint = worldToScreenPoint(item.x, item.y);
                const screenRadius = item.radius * state.viewTransform.scale;
                if (screenPoint.x + screenRadius < 0 || screenPoint.x - screenRadius > clusterCanvas.width ||
                    screenPoint.y + screenRadius < 0 || screenPoint.y - screenRadius > clusterCanvas.height) {{
                    return;
                }}
                const component = state.bundle.graph.hierarchy.nodes[item.componentId];
                const text = component.size.toLocaleString();
                if (!clusterCountLabelFits(text, item)) {{
                    return;
                }}
                clusterContext.fillText(text, screenPoint.x, screenPoint.y + 4);
            }});
        }}

        // "Label by" metadata labels: one per member node, drawn next to its dot. Shown only when
        // few enough nodes are in view (collected above with a cap), so they appear as you zoom in.
        if (collectDotLabels && !dotLabelsOverflow && dotLabelCandidates.length > 0) {{
            clusterContext.fillStyle = '#1e2a2f';
            clusterContext.font = '12px Georgia';
            clusterContext.textAlign = 'center';
            clusterContext.textBaseline = 'middle';
            dotLabelCandidates.forEach(candidate => {{
                let text;
                if (dotLabelField === '__node_id__') {{
                    text = nodeId(candidate.memberIndex);
                }} else {{
                    const value = metadataValue(candidate.memberIndex, dotLabelField);
                    // Null/empty values are simply not labeled (no node_id fallback).
                    if (value === null || value === undefined || value === '') {{
                        return;
                    }}
                    text = formatValue(value);
                }}
                if (!text) {{
                    return;
                }}
                const screenPoint = worldToScreenPoint(candidate.x, candidate.y);
                clusterContext.fillText(String(text).slice(0, 24), screenPoint.x, screenPoint.y);
            }});
        }}

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
            clusterContext.fillStyle = 'rgba(255, 255, 255, 0.72)';
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
        const links = (layoutAlgorithm === 'grid' || layoutAlgorithm === 'packed' || layoutAlgorithm === 'treemap')
            ? []
            : visibleGraph.links
                .map(link => {{
                    const left = layoutById.get(link.sourceId);
                    const right = layoutById.get(link.targetId);
                    if (!left || !right) {{ return null; }}
                    return {{left, right, threshold: link.weight, collapsed: link.collapsed || 0}};
                }})
                .filter(Boolean);
        state.visibleLayout = layout;
        state.splitLinks = links;
        // For the lattice layout, map every node to the visible cluster item that owns it so a
        // click can resolve cell -> node -> cluster in O(1) (bounding boxes overlap, so we can't
        // hit-test clusters by bbox). Rebuilt whenever the visible cluster set changes.
        if (layout.length > 0 && layout[0].shape === 'lattice') {{
            const leafOrder = state.bundle.graph.hierarchy.leaf_order;
            const nodeItem = new Int32Array(leafOrder.length).fill(-1);
            layout.forEach((item, itemIndex) => {{
                const end = item.leafStart + item.leafCount;
                for (let position = item.leafStart; position < end; position++) {{
                    nodeItem[leafOrder[position]] = itemIndex;
                }}
            }});
            state.latticeNodeItem = nodeItem;
        }} else {{
            state.latticeNodeItem = null;
        }}
        document.getElementById('stat-clusters').textContent = layout.length.toLocaleString();
        document.getElementById('stat-links').textContent = links.length.toLocaleString();
        if (state.pendingRestoreViewTransform) {{
            // A session restored a pan/zoom; adopt it instead of auto-fitting,
            // which would otherwise immediately discard it.
            const restored = state.pendingRestoreViewTransform;
            state.pendingRestoreViewTransform = null;
            state.viewTransform.scale = restored.scale;
            state.viewTransform.offsetX = restored.offsetX;
            state.viewTransform.offsetY = restored.offsetY;
        }} else if (resetView) {{
            fitClusterViewToLayout();
        }}
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
        const hierarchyNodes = state.bundle.graph.hierarchy.nodes;
        const sortBySizeEnabled = sortComponentsBySizeEnabled();

        // The lattice ("treemap") layout never hides nodes -- every active cluster is laid out and
        // the minimum cluster size only suppresses small clusters' outlines. Every other layout
        // filters nodes by minimum cluster size as usual.
        const latticeMode = layoutAlgorithm === 'treemap';
        const visibleGraph = latticeMode
            ? {{visibleIds: activeClusterIds, links: [], hiddenNodes: 0, collapsedPaths: 0, collapsedClusters: 0}}
            : mstLinksForActiveClusters(activeClusterIds, minClusterSize,
                leafPruningOnlyEnabled(), collapseLongPathsEnabled());

        state.activeClusters = activeClusterIds;
        state.visibleClusters = visibleGraph.visibleIds;

        const hidden = visibleGraph.hiddenNodes;
        const shown = state.bundle.graph.nodes.length - hidden;
        // Contracted paths are part of the hidden count, but they are hidden for a second
        // reason and are worth naming: the dashed edges they leave behind are the only
        // places the drawing shows a path weight rather than one MST edge.
        const collapsedPaths = visibleGraph.collapsedPaths || 0;
        const collapsedNote = collapsedPaths > 0
            ? ', ' + collapsedPaths.toLocaleString() + ' path' + (collapsedPaths === 1 ? '' : 's') + ' collapsed'
            : '';
        document.getElementById('hidden-summary').textContent = latticeMode
            ? 'All nodes shown; minimum cluster size hides only small-cluster outlines'
            : (hidden.toLocaleString() + ' nodes hidden by minimum cluster size' + collapsedNote);
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

{layout_core_js}
{session_state_js}
{selection_presets_js}
{table_editing_js}
{extraction_js}
{gradient_stops_js}
{column_charts_js}
{split_chart_hover_js}

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
            // Checked before the sort button: the chart glyph must never be
            // given a data-column-key, or clicking it would also sort.
            const chartBtn = event.target.closest('[data-chart-column]');
            if (chartBtn) {{
                event.preventDefault();
                event.stopPropagation();
                openColumnChartMenu(chartBtn.dataset.chartColumn, chartBtn);
                return;
            }}
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
            // The second click of a double-click opens the cell editor; letting it
            // through would toggle the row selection on the way in.
            if (event.detail > 1) {{ return; }}
            if (!event.shiftKey && metadataTableHasActiveTextSelection()) {{ return; }}
            const row = event.target.closest('tr[data-node-index]');
            if (!row) {{ return; }}
            toggleMetadataRowSelection(Number(row.dataset.nodeIndex), state.renderedNodeIndices, {{range: event.shiftKey}});
        }});

        tbody.addEventListener('keydown', event => {{
            if (event.key !== 'Enter' && event.key !== ' ') {{ return; }}
            if (event.target.closest('button, a, input, select, textarea')) {{ return; }}
            const row = event.target.closest('tr[data-node-index]');
            if (!row) {{ return; }}
            event.preventDefault();
            toggleMetadataRowSelection(Number(row.dataset.nodeIndex), state.renderedNodeIndices, {{range: event.shiftKey}});
        }});
    }}

    function updateMetadataTable() {{
        // The table body is re-rendered wholesale, so any in-progress cell edit
        // is destroyed along with its <input>. The header goes the same way, so
        // an open chart menu would be left anchored to a detached button.
        state.metadataEditCell = null;
        closeColumnChartMenu();
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
                // No chart glyph on node_id: it is the graph key, unique by
                // construction, so every chart of it is degenerate (N rows of
                // count 1, N equal slices). A disabled button would only invite
                // a "why?" with no answer.
                + (label === 'node_id' ? ''
                    : '<button type="button" class="metadata-chart-button" data-chart-column="' + escapedLabel
                        + '" aria-haspopup="menu" aria-expanded="false" title="Charts and tables for ' + escapedLabel + '">'
                        + '<svg class="metadata-chart-glyph" viewBox="0 0 12 12" aria-hidden="true" focusable="false">'
                        + '<rect x="1" y="6" width="2.6" height="5"></rect>'
                        + '<rect x="4.7" y="3" width="2.6" height="8"></rect>'
                        + '<rect x="8.4" y="1" width="2.6" height="10"></rect></svg></button>')
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
                    const escapedName = htmlEscape(column.name);
                    const displayText = htmlEscape(formatMetadataDisplayValue(column.name, value));
                    cellsHtml += '<td data-column-key="' + escapedName + '"'
                        + ' class="metadata-cell-editable' + (cellClass ? ' ' + cellClass : '') + '"'
                        + ' title="' + displayText + ' (double-click to edit)">' + displayText + '</td>';
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
        document.getElementById('focus-selection').disabled = selected.length === 0;
        document.getElementById('save-extraction').disabled = !state.bundle || selected.length === 0;
        applyMetadataTableRowHighlights();
        // Charts summarize the rows this function just resolved, so a changed
        // selection, filter or edit has to reach an open chart.
        refreshColumnChartIfOpen();
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
        // "Add to selection" searches the whole network, so it works from a standing
        // start; the two narrowing actions have nothing to narrow without a selection.
        document.getElementById('metadata-select-add').disabled = !state.bundle;
        document.getElementById('metadata-select-remove').disabled = !state.bundle || state.selectedNodeIndices.size === 0;
        document.getElementById('metadata-select-subset').disabled = !state.bundle || state.selectedNodeIndices.size === 0;
    }}

    // Cluster numbers for the current threshold, keyed by node index: connected
    // components ranked by size with the largest numbered 1. This mirrors
    // build_ssn.py's --cluster (CLUSTER_COLUMN, rename_labels_by_frequency + 1),
    // so a column made here lines up with what that tool writes. Equal-sized
    // clusters are ordered by hierarchy position rather than by lowest member
    // index, so their numbers can differ from build_ssn.py's on ties.
    function clusterNumbersAtCurrentThreshold() {{
        const activeClusterIds = activeClustersAtThreshold(selectedThresholdValue());
        const componentAssignments = activeClusterAssignments(activeClusterIds);
        const rankedClusterIds = [...activeClusterIds].sort((leftId, rightId) => {{
            const leftNode = state.bundle.graph.hierarchy.nodes[leftId];
            const rightNode = state.bundle.graph.hierarchy.nodes[rightId];
            return rightNode.size - leftNode.size || leftNode.leaf_start - rightNode.leaf_start || leftId - rightId;
        }});
        const clusterByComponentId = new Map(rankedClusterIds.map((componentId, clusterIndex) => [componentId, clusterIndex + 1]));
        return {{
            clusterCount: activeClusterIds.length,
            numberFor: nodeIndex => clusterByComponentId.get(componentAssignments[nodeIndex]) ?? null,
        }};
    }}

    function exportSelection() {{
        const selected = Array.from(state.selectedNodeIndices).sort((left, right) => left - right);
        if (!state.bundle) {{
            return;
        }}
        const clusterNumbers = clusterNumbersAtCurrentThreshold();
        const metadataColumns = state.metadataColumns.filter(column => column.name !== 'SSN_cluster');
        const exportedNodeIndices = metadataDisplayNodeIndices(selected.length > 0
            ? selected
            : state.bundle.graph.nodes.map((_, nodeIndex) => nodeIndex));
        const header = ['node_id', 'SSN_cluster', ...metadataColumns.map(column => column.name)];
        const lines = [header.join('\\t')];
        exportedNodeIndices.forEach(nodeIndex => {{
            const row = [
                nodeId(nodeIndex),
                clusterNumbers.numberFor(nodeIndex) ?? '',
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

    function exportBaseName() {{
        return ((state.bundle && state.bundle.name) || 'network').replace(/\\s+/g, '_');
    }}

    function triggerDownload(href, filename, revoke) {{
        const link = document.createElement('a');
        link.href = href;
        link.download = filename;
        link.click();
        if (revoke) {{
            URL.revokeObjectURL(href);
        }}
    }}

    function selectedPngScale() {{
        const select = document.getElementById('export-png-scale');
        const value = select ? Number(select.value) : 1;
        return Number.isFinite(value) && value >= 1 ? value : 1;
    }}

    function exportClusterPNG() {{
        if (!state.bundle) {{
            return;
        }}
        const scaleFactor = selectedPngScale();
        // Render into an offscreen canvas scaled by the chosen factor so the PNG
        // is a true high-resolution re-rasterization (not an upscaled snapshot).
        // The view transform and logical dimensions are unchanged, so the image
        // shows exactly the current view at higher pixel density. Rendering
        // offscreen also keeps any transient selection box out of the export.
        const target = document.createElement('canvas');
        target.width = Math.round(clusterCanvas.width * scaleFactor);
        target.height = Math.round(clusterCanvas.height * scaleFactor);
        const targetContext = target.getContext('2d');
        if (!targetContext) {{
            window.alert('Could not export PNG at ' + scaleFactor + '×; try a lower resolution.');
            return;
        }}
        targetContext.scale(scaleFactor, scaleFactor);
        renderClusterView(targetContext, clusterCanvas.width, clusterCanvas.height, {{preview: false}});
        target.toBlob(blob => {{
            if (!blob) {{
                // Browsers (notably Safari, ~16.7M px) refuse to encode an
                // oversized canvas and hand back null. Tell the user instead of
                // failing silently.
                window.alert('The view is too large to export as a ' + scaleFactor +
                    '× PNG (' + target.width + '×' + target.height + ' px). Try a lower resolution or zoom in.');
                return;
            }}
            const href = URL.createObjectURL(blob);
            const suffix = scaleFactor > 1 ? '_view@' + scaleFactor + 'x.png' : '_view.png';
            triggerDownload(href, exportBaseName() + suffix, true);
        }}, 'image/png');
    }}

    function escapeXml(text) {{
        return String(text)
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;');
    }}

    // Normalize any CSS color (hex, rgb, rgba, hsl) into an SVG-friendly
    // {{color, opacity}} pair using the canvas's own color parser. SVG editors
    // such as Illustrator/Inkscape choke on rgba()/space-separated hsl(), so we
    // split out the alpha into a separate opacity attribute.
    function svgColorParts(css) {{
        clusterContext.fillStyle = '#000000';
        clusterContext.fillStyle = css;
        const normalized = clusterContext.fillStyle;
        const match = /^rgba\\(\\s*([\\d.]+)\\s*,\\s*([\\d.]+)\\s*,\\s*([\\d.]+)\\s*,\\s*([\\d.]+)\\s*\\)$/.exec(normalized);
        if (match) {{
            return {{color: 'rgb(' + match[1] + ',' + match[2] + ',' + match[3] + ')', opacity: Number(match[4])}};
        }}
        return {{color: normalized, opacity: 1}};
    }}

    // Normalize any CSS color (hex, rgb, hsl) to '#rrggbb' using the canvas parser.
    // <input type="color"> only accepts '#rrggbb', but the default schemes emit
    // hsl(...), so color inputs must be seeded through this. Alpha is dropped.
    function cssToHex(css) {{
        clusterContext.fillStyle = '#000000';
        clusterContext.fillStyle = css;
        const normalized = clusterContext.fillStyle;
        if (normalized.charAt(0) === '#') {{
            return normalized.length === 7 ? normalized : normalized.slice(0, 7);
        }}
        const match = /rgba?\\(\\s*([\\d.]+)\\s*,\\s*([\\d.]+)\\s*,\\s*([\\d.]+)/.exec(normalized);
        if (!match) {{
            return '#000000';
        }}
        const toHex = n => Math.round(Number(n)).toString(16).padStart(2, '0');
        return '#' + toHex(match[1]) + toHex(match[2]) + toHex(match[3]);
    }}

    // Render the current cluster view to a standalone SVG document. Geometry is
    // baked into screen space (every point passed through worldToScreenPoint)
    // so the markup needs no transforms and matches the canvas pixel-for-pixel.
    // This mirrors renderClusterView(); keep the two in sync.
    function buildClusterViewSVG() {{
        const width = clusterCanvas.width;
        const height = clusterCanvas.height;
        const fmt = value => Math.round(value * 100) / 100;
        const parts = [];
        parts.push('<svg xmlns="http://www.w3.org/2000/svg" width="' + width + '" height="' + height +
            '" viewBox="0 0 ' + width + ' ' + height + '" font-family="Georgia, serif">');
        parts.push('<rect x="0" y="0" width="' + width + '" height="' + height + '" fill="#ffffff"/>');

        if (!state.bundle) {{
            parts.push('<text x="30" y="50" font-size="18" fill="#5c6a70">Load a bundle to render cluster bubbles.</text>');
            parts.push('</svg>');
            return parts.join('\\n');
        }}

        const TAU = Math.PI * 2;
        const drawScale = Math.max(state.viewTransform.scale, 1e-9);
        const worldMinX = (0 - state.viewTransform.offsetX) / drawScale;
        const worldMaxX = (width - state.viewTransform.offsetX) / drawScale;
        const worldMinY = (0 - state.viewTransform.offsetY) / drawScale;
        const worldMaxY = (height - state.viewTransform.offsetY) / drawScale;

        // Split links (edges). Collapsed paths go in their own dashed group, matching
        // renderClusterView.
        const linkColor = svgColorParts('rgba(92, 106, 112, 0.42)');
        const linkPaths = [];
        const collapsedLinkPaths = [];
        state.splitLinks.forEach(link => {{
            const segments = renderedLinkSegments(link);
            if (segments.length === 0) {{
                return;
            }}
            let d = '';
            segments.forEach((segment, index) => {{
                const start = worldToScreenPoint(segment.startX, segment.startY);
                const end = worldToScreenPoint(segment.endX, segment.endY);
                if (index === 0) {{
                    d += 'M' + fmt(start.x) + ' ' + fmt(start.y);
                }} else {{
                    d += ' L' + fmt(start.x) + ' ' + fmt(start.y);
                }}
                d += ' L' + fmt(end.x) + ' ' + fmt(end.y);
            }});
            (link.collapsed ? collapsedLinkPaths : linkPaths).push('<path d="' + d + '"/>');
        }});
        if (linkPaths.length > 0) {{
            parts.push('<g fill="none" stroke="' + linkColor.color + '" stroke-opacity="' + linkColor.opacity +
                '" stroke-width="1.6">');
            parts.push(linkPaths.join(''));
            parts.push('</g>');
        }}
        if (collapsedLinkPaths.length > 0) {{
            parts.push('<g fill="none" stroke="' + linkColor.color + '" stroke-opacity="' + linkColor.opacity +
                '" stroke-width="1.6" stroke-dasharray="' + LINK_DASH_ON + ' ' + LINK_DASH_OFF + '">');
            parts.push(collapsedLinkPaths.join(''));
            parts.push('</g>');
        }}

        const showClusterBounds = renderClusterBoundsEnabled();
        const showNodes = renderNodesEnabled();
        const minOutlineSize = Math.max(1, Number(document.getElementById('min-cluster-size').value) || 1);
        const lattice = state.latticeGlobal;
        const dotLabelField = currentLabelField();
        const MAX_LABELED_DOTS = 250;
        // Whether a dot is drawn large enough to be worth pointing at is decided per
        // dot below; MAX_LABELED_DOTS is what keeps a crowded view readable.
        const collectDotLabels = dotLabelField !== '' && showNodes;
        const dotLabelCandidates = [];
        let dotLabelsOverflow = false;

        const boundsFill = svgColorParts('rgba(245, 246, 248, 0.95)');
        const boundsParts = [];
        const dotColorBuckets = new Map();
        const squareColorBuckets = new Map();
        const selectedNodeOutlines = [];

        state.visibleLayout.forEach(item => {{
            if (item.x + item.radius < worldMinX || item.x - item.radius > worldMaxX ||
                item.y + item.radius < worldMinY || item.y - item.radius > worldMaxY) {{
                return;
            }}

            const component = state.bundle.graph.hierarchy.nodes[item.componentId];
            const members = componentMembers(item.componentId);
            const selectionState = componentSelectionState(members);
            const center = worldToScreenPoint(item.x, item.y);
            const screenRadius = item.radius * drawScale;
            const isLattice = item.shape === 'lattice';
            const isRect = item.shape === 'rect';
            const rectTopLeft = isRect ? worldToScreenPoint(item.x0, item.y0) : null;
            const rectWidth = isRect ? (item.x1 - item.x0) * drawScale : 0;
            const rectHeight = isRect ? (item.y1 - item.y0) * drawScale : 0;
            const boundsShape = (fillAttrs, strokeAttrs) => (isRect
                ? '<rect x="' + fmt(rectTopLeft.x) + '" y="' + fmt(rectTopLeft.y) + '" width="' + fmt(rectWidth) +
                    '" height="' + fmt(rectHeight) + '" ' + fillAttrs + ' ' + strokeAttrs + '/>'
                : '<circle cx="' + fmt(center.x) + '" cy="' + fmt(center.y) + '" r="' + fmt(screenRadius) +
                    '" ' + fillAttrs + ' ' + strokeAttrs + '/>');

            if (isLattice) {{
                // Staircase outline traced along the padding between cells (screen space).
                const outline = showClusterBounds && component.size >= minOutlineSize;
                if ((outline || selectionState.allSelected) && lattice) {{
                    const width = lattice.width;
                    const cellSet = new Set();
                    for (const nodeIndex of members) {{
                        cellSet.add((lattice.rowOf[nodeIndex] * width) + lattice.colOf[nodeIndex]);
                    }}
                    let d = '';
                    const edge = (wx0, wy0, wx1, wy1) => {{
                        const a = worldToScreenPoint(wx0, wy0);
                        const b = worldToScreenPoint(wx1, wy1);
                        d += 'M' + fmt(a.x) + ' ' + fmt(a.y) + 'L' + fmt(b.x) + ' ' + fmt(b.y);
                    }};
                    for (const nodeIndex of members) {{
                        const column = lattice.colOf[nodeIndex];
                        const row = lattice.rowOf[nodeIndex];
                        const ex0 = TREEMAP_ORIGIN + (column * TREEMAP_CELL);
                        const ey0 = TREEMAP_ORIGIN + (row * TREEMAP_CELL);
                        const ex1 = ex0 + TREEMAP_CELL;
                        const ey1 = ey0 + TREEMAP_CELL;
                        if (!cellSet.has((row * width) + (column - 1))) {{ edge(ex0, ey0, ex0, ey1); }}
                        if (!cellSet.has((row * width) + (column + 1))) {{ edge(ex1, ey0, ex1, ey1); }}
                        if (!cellSet.has(((row - 1) * width) + column)) {{ edge(ex0, ey0, ex1, ey0); }}
                        if (!cellSet.has(((row + 1) * width) + column)) {{ edge(ex0, ey1, ex1, ey1); }}
                    }}
                    if (d !== '') {{
                        const stroke = svgColorParts(selectionState.allSelected ? '#1e2a2f' : (outline ? '#9aa4ad' : 'rgba(30, 42, 47, 0.45)'));
                        boundsParts.push('<path d="' + d + '" fill="none" stroke="' + stroke.color +
                            (stroke.opacity !== 1 ? '" stroke-opacity="' + stroke.opacity : '') +
                            '" stroke-width="' + (selectionState.allSelected ? 2.4 : 1.2) + '"/>');
                    }}
                }}
            }} else if (showClusterBounds) {{
                const stroke = svgColorParts(selectionState.allSelected ? '#1e2a2f' : '#c4cad2');
                boundsParts.push(boundsShape(
                    'fill="' + boundsFill.color + '" fill-opacity="' + boundsFill.opacity + '"',
                    'stroke="' + stroke.color + '" stroke-width="' + (selectionState.allSelected ? 3 : 1.5) + '"'));
            }} else if (selectionState.allSelected) {{
                const stroke = svgColorParts('rgba(30, 42, 47, 0.45)');
                boundsParts.push(boundsShape(
                    'fill="none"',
                    'stroke="' + stroke.color + '" stroke-opacity="' + stroke.opacity + '" stroke-width="1.5"'));
            }}

            const needsSelectionHint = !selectionState.allSelected && selectionState.anySelected;
            if (!showNodes && !needsSelectionHint) {{
                return;
            }}

            const dotLayout = componentMemberLayout(component, item);
            for (const dot of dotLayout) {{
                if (showNodes) {{
                    const color = state.nodeColorCache[dot.memberIndex] ?? nodeColor(dot.memberIndex);
                    const buckets = dot.square ? squareColorBuckets : dotColorBuckets;
                    let bucket = buckets.get(color);
                    if (bucket === undefined) {{
                        bucket = [];
                        buckets.set(color, bucket);
                    }}
                    const dotCenter = worldToScreenPoint(dot.x, dot.y);
                    bucket.push({{x: dotCenter.x, y: dotCenter.y, r: dot.radius * drawScale}});
                }}
                if (!selectionState.allSelected && state.selectedNodeIndices.has(dot.memberIndex)) {{
                    selectedNodeOutlines.push({{x: dot.x, y: dot.y, radius: dot.radius, faint: !showNodes}});
                }}
                if (collectDotLabels && !dotLabelsOverflow &&
                    dot.x >= worldMinX && dot.x <= worldMaxX && dot.y >= worldMinY && dot.y <= worldMaxY) {{
                    if (dot.radius * state.viewTransform.scale < MIN_LABELED_DOT_SCREEN_RADIUS) {{
                        // Sub-pixel dot: a label beside it would point at nothing visible.
                    }} else if (dotLabelCandidates.length >= MAX_LABELED_DOTS) {{
                        dotLabelsOverflow = true;
                    }} else {{
                        dotLabelCandidates.push({{x: dot.x, y: dot.y, memberIndex: dot.memberIndex}});
                    }}
                }}
            }}
        }});

        if (boundsParts.length > 0) {{
            parts.push(boundsParts.join(''));
        }}

        dotColorBuckets.forEach((dots, color) => {{
            const fill = svgColorParts(color);
            const circles = dots.map(dot => '<circle cx="' + fmt(dot.x) + '" cy="' + fmt(dot.y) + '" r="' + fmt(dot.r) + '"/>').join('');
            parts.push('<g fill="' + fill.color + '"' + (fill.opacity !== 1 ? ' fill-opacity="' + fill.opacity + '"' : '') + '>' + circles + '</g>');
        }});

        squareColorBuckets.forEach((squares, color) => {{
            const fill = svgColorParts(color);
            const rects = squares.map(square => '<rect x="' + fmt(square.x - square.r) + '" y="' + fmt(square.y - square.r) +
                '" width="' + fmt(square.r * 2) + '" height="' + fmt(square.r * 2) + '"/>').join('');
            parts.push('<g fill="' + fill.color + '"' + (fill.opacity !== 1 ? ' fill-opacity="' + fill.opacity + '"' : '') + '>' + rects + '</g>');
        }});

        if (selectedNodeOutlines.length > 0) {{
            const outlineParts = selectedNodeOutlines.map(outline => {{
                const screenPoint = worldToScreenPoint(outline.x, outline.y);
                const screenRadius = Math.max(3, (outline.radius * state.viewTransform.scale) + 1.6);
                if (outline.faint) {{
                    const fill = svgColorParts('rgba(30, 42, 47, 0.30)');
                    return '<circle cx="' + fmt(screenPoint.x) + '" cy="' + fmt(screenPoint.y) + '" r="' + fmt(screenRadius) +
                        '" fill="' + fill.color + '" fill-opacity="' + fill.opacity + '"/>';
                }}
                return '<circle cx="' + fmt(screenPoint.x) + '" cy="' + fmt(screenPoint.y) + '" r="' + fmt(screenRadius) +
                    '" fill="none" stroke="#1e2a2f" stroke-width="1.8"/>';
            }});
            parts.push(outlineParts.join(''));
        }}

        // Edge score badges. The fit test is edgeScoreLabelFits, shared with
        // renderClusterView so the export shows the labels the screen showed.
        if (showEdgeScoresEnabled()) {{
            clusterContext.save();
            clusterContext.font = EDGE_SCORE_LABEL_FONT + 'px Georgia';
            state.splitLinks.forEach(link => {{
                const anchor = linkLabelAnchor(link);
                if (anchor === null) {{
                    return;
                }}
                const text = formatValue(link.threshold);
                if (!edgeScoreLabelFits(text, anchor.length * state.viewTransform.scale)) {{
                    return;
                }}
                const point = worldToScreenPoint(anchor.x, anchor.y);
                const x = point.x + (anchor.normalX * EDGE_SCORE_LABEL_OFFSET);
                const y = point.y + (anchor.normalY * EDGE_SCORE_LABEL_OFFSET);
                // Measured here rather than estimated: the badge is being drawn, and
                // the rectangle has to actually enclose the glyphs.
                const textWidth = clusterContext.measureText(text).width;
                const badgeWidth = textWidth + EDGE_SCORE_BADGE_PADDING;
                const badgeHeight = 18;
                parts.push('<rect x="' + fmt(x - (badgeWidth / 2)) + '" y="' + fmt(y - (badgeHeight / 2)) +
                    '" width="' + fmt(badgeWidth) + '" height="' + badgeHeight + '" rx="8" ry="8"' +
                    ' fill="#ffffff" fill-opacity="0.92" stroke="rgb(92,106,112)" stroke-opacity="0.32" stroke-width="1"/>');
                parts.push('<text x="' + fmt(x) + '" y="' + fmt(y + 0.5) + '" font-size="11" fill="#334147"' +
                    ' text-anchor="middle" dominant-baseline="central">' + escapeXml(text) + '</text>');
            }});
            clusterContext.restore();
        }}

        // Node-count labels (one per bubble), under the same fit test as the canvas.
        if (showNodeCountsEnabled()) {{
            state.visibleLayout.forEach(item => {{
                const screenPoint = worldToScreenPoint(item.x, item.y);
                const screenRadius = item.radius * state.viewTransform.scale;
                if (screenPoint.x + screenRadius < 0 || screenPoint.x - screenRadius > width ||
                    screenPoint.y + screenRadius < 0 || screenPoint.y - screenRadius > height) {{
                    return;
                }}
                const component = state.bundle.graph.hierarchy.nodes[item.componentId];
                const text = component.size.toLocaleString();
                if (!clusterCountLabelFits(text, item)) {{
                    return;
                }}
                parts.push('<text x="' + fmt(screenPoint.x) + '" y="' + fmt(screenPoint.y + 4) + '" font-size="' +
                    CLUSTER_COUNT_LABEL_FONT + '" font-weight="600"' +
                    ' fill="#5c6a70" text-anchor="middle" dominant-baseline="central">' + escapeXml(text) + '</text>');
            }});
        }}

        // Per-node "Label by" metadata labels.
        if (collectDotLabels && !dotLabelsOverflow && dotLabelCandidates.length > 0) {{
            dotLabelCandidates.forEach(candidate => {{
                let text;
                if (dotLabelField === '__node_id__') {{
                    text = nodeId(candidate.memberIndex);
                }} else {{
                    const value = metadataValue(candidate.memberIndex, dotLabelField);
                    if (value === null || value === undefined || value === '') {{
                        return;
                    }}
                    text = formatValue(value);
                }}
                if (!text) {{
                    return;
                }}
                const screenPoint = worldToScreenPoint(candidate.x, candidate.y);
                parts.push('<text x="' + fmt(screenPoint.x) + '" y="' + fmt(screenPoint.y) + '" font-size="12" fill="#1e2a2f"' +
                    ' text-anchor="middle" dominant-baseline="central">' + escapeXml(String(text).slice(0, 24)) + '</text>');
            }});
        }}

        parts.push('</svg>');
        return parts.join('\\n');
    }}

    function exportClusterSVG() {{
        if (!state.bundle) {{
            return;
        }}
        const svg = buildClusterViewSVG();
        const blob = new Blob([svg], {{type: 'image/svg+xml'}});
        const href = URL.createObjectURL(blob);
        triggerDownload(href, exportBaseName() + '_view.svg', true);
    }}

    // ---- Custom color palette: load/save color table TSV (categorical only) ----

    // Guard shared by the TSV load/save paths: a 2-column value->color table only
    // describes discrete mappings, so it applies to categorical columns. Returns the
    // column name, or null after setting an explanatory status message.
    function categoricalColumnForColorTable(action) {{
        const columnName = currentColorField();
        const info = colorInfo(columnName);
        if (!columnName || !info) {{
            setStatus(action + ' color table: choose a "Color by" column first.');
            return null;
        }}
        if (info.type !== 'categorical') {{
            setStatus(action + ' color table: only categorical columns use a value/color table (use the picker for "' + columnName + '").');
            return null;
        }}
        return columnName;
    }}

    function loadColorTableText(text) {{
        const columnName = categoricalColumnForColorTable('Load');
        if (!columnName) {{
            return;
        }}
        const colors = {{}};
        let nullColor = null;
        let loaded = 0;
        let skipped = 0;
        String(text).split('\\n').forEach(rawLine => {{
            const line = rawLine.trim();
            if (line === '') {{
                return;
            }}
            const parts = line.split('\\t');
            const hex = parts.length === 2 ? normalizeColorHex(parts[1]) : null;
            if (parts.length !== 2 || !hex) {{
                skipped += 1;
                return;
            }}
            if (parts[0] === '\\u2014') {{
                nullColor = hex;
            }} else {{
                colors[parts[0]] = hex;
            }}
            loaded += 1;
        }});
        state.customPalettes[paletteKey(columnName)] = {{type: 'categorical', colors, nullColor}};
        rebuildNodeColorCache();
        renderClusterView();
        if (colorPickerIsOpen()) {{
            openColorPicker();
        }}
        setStatus('Loaded ' + loaded + ' colors for "' + columnName + '"' +
            (skipped > 0 ? ', skipped ' + skipped + ' malformed line(s).' : '.'));
    }}

    function saveColorTableTSV() {{
        const columnName = categoricalColumnForColorTable('Save');
        if (!columnName) {{
            return;
        }}
        const palette = customPalette(columnName);
        const distinct = distinctColumnValues(columnName, Infinity);
        // Write the effective on-screen color for every distinct value (custom override
        // or the current default), so the file matches the view and round-trips on load.
        const rows = distinct.values.map(entry => ({{
            key: entry.key,
            color: (palette && palette.colors && palette.colors[entry.key]) || cssToHex(categoricalColor(entry.raw, palette)),
        }}));
        rows.sort((a, b) => a.key.localeCompare(b.key));
        if (distinct.nullCount > 0) {{
            // Em-dash row carries the no-value color; written last (null sorts last).
            rows.push({{key: '\\u2014', color: (palette && palette.nullColor) || cssToHex('#b3a89d')}});
        }}
        const tsv = rows.map(row => row.key + '\\t' + row.color).join('\\n') + '\\n';
        const blob = new Blob([tsv], {{type: 'text/tab-separated-values'}});
        const href = URL.createObjectURL(blob);
        triggerDownload(href, exportBaseName() + '_' + columnName.replace(/\\s+/g, '_') + '_colors.tsv', true);
    }}

    // ---- Custom color palette: legend export (SVG + PNG) ----

    // Build a standalone legend SVG for the current color column. Returns
    // {{svg, width, height}}. Discrete columns get a box+label per value (ported from
    // color_table_to_legend.py); numeric columns get a sampled gradient bar with
    // min/mid/max ticks. All fills pass through svgColorParts so default hsl() colors
    // become SVG-editor-safe rgb()+opacity. No external refs (keeps PNG untainted).
    function buildLegendSVG() {{
        const columnName = currentColorField();
        const info = colorInfo(columnName);
        const title = columnName || 'Legend';
        const titleFontSize = 24;
        const itemFontSize = 20;
        const strokeWidth = 2;
        const padding = 10;
        const titleSpace = 40;
        const parts = [];
        const open = (width, height) => {{
            parts.push('<?xml version="1.0" encoding="UTF-8" standalone="no"?>');
            parts.push('<svg xmlns="http://www.w3.org/2000/svg" width="' + width + '" height="' + height +
                '" viewBox="0 0 ' + width + ' ' + height + '">');
            parts.push('<g font-family="Arial, sans-serif">');
            parts.push('<rect x="0" y="0" width="' + width + '" height="' + height +
                '" fill="white" stroke="black" stroke-width="' + strokeWidth + '"/>');
            parts.push('<text x="' + padding + '" y="' + titleFontSize + '" font-size="' + titleFontSize +
                '" font-weight="bold">' + escapeXml(title) + '</text>');
        }};
        const swatchRow = (color, label, y, boxHeight, boxWidth) => {{
            const fill = svgColorParts(color);
            parts.push('<rect x="' + padding + '" y="' + y + '" width="' + boxWidth + '" height="' + (boxHeight - 10) +
                '" fill="' + fill.color + '" fill-opacity="' + fill.opacity +
                '" stroke="black" stroke-width="' + strokeWidth + '"/>');
            parts.push('<text x="' + (boxWidth + padding * 2) + '" y="' + (y + boxHeight / 2) +
                '" font-size="' + itemFontSize + '" dominant-baseline="central">' + escapeXml(label) + '</text>');
        }};

        if (!info || info.type !== 'numeric') {{
            // Discrete legend.
            const palette = (info && customPalette(columnName)) || null;
            // Enumerate every distinct value (no cap) so the saved legend is complete.
            const distinct = distinctColumnValues(columnName, Infinity);
            const labels = distinct.values.map(entry => ({{color: categoricalColor(entry.raw, palette), label: entry.label}}));
            if (distinct.nullCount > 0) {{
                labels.push({{color: (palette && palette.nullColor) || '#b3a89d', label: '\\u2014'}});
            }}
            const boxHeight = itemFontSize * 2;
            const boxWidth = boxHeight;
            const longest = labels.reduce((max, item) => Math.max(max, item.label.length), 1);
            const height = titleSpace + (labels.length * boxHeight) + (2 * padding);
            const width = Math.max(
                boxWidth + (padding * 3) + ((longest * itemFontSize) / 2),
                (padding * 3) + ((title.length * titleFontSize) / 2)
            );
            open(width, height);
            let y = titleSpace;
            labels.forEach(item => {{
                swatchRow(item.color, item.label, y, boxHeight, boxWidth);
                y += boxHeight;
            }});
            parts.push('</g>');
            parts.push('</svg>');
            return {{svg: parts.join('\\n'), width, height}};
        }}

        // Continuous legend: sampled gradient bar with low/mid/high ticks. Custom
        // bounds (lb/ub) and mid value reshape the domain shown on the bar.
        const palette = customPalette(columnName);
        const distinct = distinctColumnValues(columnName, 1);
        const paletteStops = numericPaletteStops(palette, info.min, info.max) || [
            {{value: info.min, color: cssToHex(numericColor(info.min, info.min, info.max))}},
            {{value: info.max, color: cssToHex(numericColor(info.max, info.min, info.max))}},
        ];
        const lb = paletteStops[0].value;
        const ub = paletteStops[paletteStops.length - 1].value;
        const barX = padding;
        const barY = titleSpace;
        const barWidth = 320;
        const barHeight = 30;
        const tickFontSize = 16;
        // One SVG stop per palette stop. The ramp is piecewise linear between them
        // and SVG interpolates the same way, so this is exact -- sampling the ramp at
        // a fixed number of points would round off the corners of a many-stop ramp.
        const stops = paletteStops.map((stop, index) => {{
            const fraction = ub > lb ? (stop.value - lb) / (ub - lb)
                : (index / Math.max(1, paletteStops.length - 1));
            const fill = svgColorParts(stop.color);
            return '<stop offset="' + (fraction * 100) + '%" stop-color="' + fill.color +
                '" stop-opacity="' + fill.opacity + '"/>';
        }});
        const nullRow = distinct.nullCount > 0;
        const height = titleSpace + barHeight + tickFontSize + (padding * 2) + (nullRow ? itemFontSize * 2 : 0);
        const width = Math.max(barWidth + (padding * 2), (padding * 3) + ((title.length * titleFontSize) / 2));
        open(width, height);
        parts.push('<defs><linearGradient id="legend-gradient" x1="0%" y1="0%" x2="100%" y2="0%">' +
            stops.join('') + '</linearGradient></defs>');
        parts.push('<rect x="' + barX + '" y="' + barY + '" width="' + barWidth + '" height="' + barHeight +
            '" fill="url(#legend-gradient)" stroke="black" stroke-width="' + strokeWidth + '"/>');
        const tickY = barY + barHeight + tickFontSize;
        const tick = (x, text, anchor) => parts.push('<text x="' + x + '" y="' + tickY + '" font-size="' + tickFontSize +
            '" text-anchor="' + anchor + '">' + escapeXml(text) + '</text>');
        // One tick per stop. Intermediates that would collide with the label before
        // them are dropped, so a many-stop ramp still reads at any stop count; the two
        // ends always draw, since they are what give the bar its scale.
        let lastRight = -Infinity;
        paletteStops.forEach((stop, index) => {{
            const isEnd = index === 0 || index === paletteStops.length - 1;
            const fraction = ub > lb ? (stop.value - lb) / (ub - lb) : 0.5;
            const x = barX + (fraction * barWidth);
            const text = formatValue(stop.value);
            const textWidth = (text.length * tickFontSize) / 2;
            const anchor = index === 0 ? 'start' : (isEnd ? 'end' : 'middle');
            const left = anchor === 'start' ? x : (anchor === 'end' ? x - textWidth : x - (textWidth / 2));
            if (!isEnd && left < lastRight + 4) {{
                return;
            }}
            tick(x, text, anchor);
            lastRight = left + textWidth;
        }});
        if (nullRow) {{
            swatchRow((palette && palette.nullColor) || '#b3a89d', '\\u2014', tickY + padding, itemFontSize * 2, itemFontSize * 2);
        }}
        parts.push('</g>');
        parts.push('</svg>');
        return {{svg: parts.join('\\n'), width, height}};
    }}

    function exportLegendSVG() {{
        if (!state.bundle) {{
            return;
        }}
        const legend = buildLegendSVG();
        const blob = new Blob([legend.svg], {{type: 'image/svg+xml'}});
        const href = URL.createObjectURL(blob);
        triggerDownload(href, exportBaseName() + '_legend.svg', true);
    }}

    function exportLegendPNG() {{
        if (!state.bundle) {{
            return;
        }}
        const legend = buildLegendSVG();
        const scaleFactor = selectedPngScale();
        const svgUrl = URL.createObjectURL(new Blob([legend.svg], {{type: 'image/svg+xml'}}));
        const image = new Image();
        image.onload = () => {{
            const target = document.createElement('canvas');
            target.width = Math.max(1, Math.round(legend.width * scaleFactor));
            target.height = Math.max(1, Math.round(legend.height * scaleFactor));
            const targetContext = target.getContext('2d');
            if (!targetContext) {{
                URL.revokeObjectURL(svgUrl);
                window.alert('Could not export legend PNG at ' + scaleFactor + 'x; try a lower resolution.');
                return;
            }}
            targetContext.drawImage(image, 0, 0, target.width, target.height);
            target.toBlob(blob => {{
                URL.revokeObjectURL(svgUrl);
                if (!blob) {{
                    window.alert('The legend is too large to export as a ' + scaleFactor + 'x PNG. Try a lower resolution.');
                    return;
                }}
                triggerDownload(URL.createObjectURL(blob), exportBaseName() + '_legend.png', true);
            }}, 'image/png');
        }};
        image.onerror = () => {{
            URL.revokeObjectURL(svgUrl);
            setStatus('Legend PNG export failed.');
        }};
        image.src = svgUrl;
    }}

    // ---- Custom color palette: interactive picker (modal dialog) ----

    function colorPickerIsOpen() {{
        const overlay = document.getElementById('color-picker-overlay');
        return overlay ? !overlay.hidden : false;
    }}

    function closeColorPicker() {{
        const overlay = document.getElementById('color-picker-overlay');
        if (overlay) {{
            overlay.hidden = true;
        }}
    }}

    function ensureCategoricalPalette(columnName) {{
        let palette = state.customPalettes[paletteKey(columnName)];
        if (!palette || palette.type !== 'categorical') {{
            // Seeded from the default assignment. Without that, storing one
            // hand-picked swatch would shadow the default palette and drop every
            // other value in the column back to an unassigned color.
            const base = defaultCategoricalPalette(columnName);
            palette = {{
                type: 'categorical',
                colors: {{...base.colors}},
                nullColor: base.nullColor,
                scheme: base.scheme,
            }};
            state.customPalettes[paletteKey(columnName)] = palette;
        }}
        return palette;
    }}

    function ensureNumericPalette(columnName) {{
        let palette = state.customPalettes[paletteKey(columnName)];
        if (!palette || palette.type !== 'numeric') {{
            palette = {{type: 'numeric', stops: [], nullColor: null}};
            state.customPalettes[paletteKey(columnName)] = palette;
        }}
        return palette;
    }}

    const COLOR_PICKER_PAGE_SIZE = 100;

    function renderDiscretePicker(columnName) {{
        const palette = customPalette(columnName);
        // Enumerate every distinct value and page through them, so columns with
        // thousands of categories stay editable (not just the first 200).
        const distinct = distinctColumnValues(columnName, Infinity);
        updateColorPaletteControl(columnName, palette, distinct.values.length);
        const total = distinct.values.length;
        const pageCount = Math.max(1, Math.ceil(total / COLOR_PICKER_PAGE_SIZE));
        const page = Math.max(0, Math.min(state.colorPickerPage, pageCount - 1));
        state.colorPickerPage = page;
        const start = page * COLOR_PICKER_PAGE_SIZE;
        const pageValues = distinct.values.slice(start, start + COLOR_PICKER_PAGE_SIZE);

        const pager = document.getElementById('color-picker-pager');
        if (total > COLOR_PICKER_PAGE_SIZE) {{
            document.getElementById('color-picker-page-status').textContent =
                'Values ' + (start + 1).toLocaleString() + '\\u2013' + (start + pageValues.length).toLocaleString() +
                ' of ' + total.toLocaleString() + ' (page ' + (page + 1) + '/' + pageCount + ')';
            document.getElementById('color-picker-prev').disabled = page === 0;
            document.getElementById('color-picker-next').disabled = page >= pageCount - 1;
            pager.hidden = false;
        }} else {{
            pager.hidden = true;
        }}

        const list = document.getElementById('color-picker-swatch-list');
        list.textContent = '';
        pageValues.forEach(entry => {{
            const row = document.createElement('div');
            row.className = 'cp-swatch-row';
            const input = document.createElement('input');
            input.type = 'color';
            const override = palette && palette.colors ? palette.colors[entry.key] : undefined;
            input.value = override || cssToHex(categoricalColor(entry.raw, palette));
            input.addEventListener('input', () => {{
                const edited = ensureCategoricalPalette(columnName);
                edited.colors[entry.key] = input.value;
                // No longer exactly a named palette; the menu says "Custom colors" next
                // time the picker renders.
                edited.scheme = null;
                rebuildNodeColorCache();
                scheduleClusterRender();
            }});
            const label = document.createElement('span');
            label.className = 'cp-swatch-label';
            label.textContent = entry.label;
            label.title = entry.label;
            const count = document.createElement('span');
            count.className = 'cp-swatch-count';
            count.textContent = entry.count.toLocaleString();
            row.appendChild(input);
            row.appendChild(label);
            row.appendChild(count);
            list.appendChild(row);
        }});
        // Shared no-value input.
        const nullInput = document.getElementById('color-null');
        nullInput.value = (palette && palette.nullColor) || cssToHex('#b3a89d');
        nullInput.oninput = () => {{
            ensureCategoricalPalette(columnName).nullColor = nullInput.value;
            rebuildNodeColorCache();
            scheduleClusterRender();
        }};
    }}

    function parseNumberOrNull(text) {{
        if (text === '' || text === null || text === undefined) {{
            return null;
        }}
        const value = Number(text);
        return Number.isNaN(value) ? null : value;
    }}

    // Bin a numeric column's values for the gradient dialog's histogram. Bins span the
    // data range, not the low/high gradient bounds, so that values the ramp clamps stay
    // visible -- seeing them pile up flat against one end is the point of the chart.
    function columnHistogram(columnName) {{
        const info = colorInfo(columnName);
        const columnIndex = state.metadataColumnIndexByName.get(columnName);
        if (!state.bundle || !info || info.baseType !== 'numeric' || columnIndex === undefined) {{
            return null;
        }}
        const values = [];
        let missing = 0;
        for (let nodeIndex = 0; nodeIndex < state.bundle.graph.nodes.length; nodeIndex++) {{
            const value = state.metadataByNodeIndex[nodeIndex]?.[columnIndex] ?? null;
            if (typeof value === 'number' && Number.isFinite(value)) {{
                values.push(value);
            }} else {{
                missing += 1;
            }}
        }}
        if (values.length === 0) {{
            return null;
        }}
        const column = state.metadataColumnByName.get(columnName);
        const min = info.min;
        const max = info.max;
        const MAX_BINS = 48;
        // Whether the column holds whole numbers, which decides both the binning
        // below and whether the range slider snaps to integers.
        const integer = !!(column && column.type === 'int'
            && Number.isInteger(min) && Number.isInteger(max));
        let binCount;
        let lowEdge;
        let highEdge;
        if (max <= min) {{
            // Every value identical: one bin, centered, rather than a zero-width domain.
            binCount = 1;
            lowEdge = min - 0.5;
            highEdge = min + 0.5;
        }} else if (integer && (max - min) + 1 <= MAX_BINS) {{
            // One bin per integer. Spreading a handful of distinct integers over evenly
            // divided bins leaves empty gaps and doubled-up bars that read as structure
            // in the data when they are only an artifact of the binning.
            //
            // span + 1 bins across [min, max] separates every integer -- value min + k
            // lands in bin floor(k + k/span), which is k for every k below span, and the
            // top value clamps into the last bin. Widening the edges to min - 0.5 and
            // max + 0.5 would separate them too, but it would push the axis half a bin
            // past the data at both ends, so the default bounds would no longer sit at
            // the ends of the bar the knobs ride on.
            binCount = (max - min) + 1;
            lowEdge = min;
            highEdge = max;
        }} else {{
            binCount = Math.min(MAX_BINS, Math.max(8, Math.ceil(Math.sqrt(values.length))));
            lowEdge = min;
            highEdge = max;
        }}
        const span = highEdge - lowEdge;
        const counts = new Array(binCount).fill(0);
        values.forEach(value => {{
            const slot = Math.floor(((value - lowEdge) / span) * binCount);
            counts[Math.max(0, Math.min(binCount - 1, slot))] += 1;
        }});
        return {{counts, lowEdge, highEdge, binCount, integer, total: values.length, missing, min, max}};
    }}

    // Paint the cached bins, filling each bar with the color that bin's values currently
    // receive. That makes the chart double as gradient feedback: a run of flat-colored
    // bars at either end is the ramp clamping, and a pale stretch with no bars under it
    // is ramp spent on a part of the range where there is no data.
    function drawColorHistogram(columnName) {{
        const wrap = document.getElementById('color-histogram-wrap');
        const histogram = state.colorHistogram;
        if (!histogram) {{
            wrap.hidden = true;
            return;
        }}
        wrap.hidden = false;
        const canvas = document.getElementById('color-histogram');
        const ctx = canvas.getContext('2d');
        ctx.clearRect(0, 0, canvas.width, canvas.height);

        const palette = customPalette(columnName);
        // No side padding: the drawing surface *is* the shared axis, so bar x
        // positions match the gradient bar and the slider knobs below.
        const paddingX = 0;
        const paddingTop = 8;
        const baselineY = canvas.height - 7;
        const plotWidth = canvas.width - (paddingX * 2);
        const plotHeight = baselineY - paddingTop;
        const slotWidth = plotWidth / histogram.binCount;
        // Hairline gaps only while the bars are wide enough to keep one.
        const gap = slotWidth > 8 ? Math.min(3, slotWidth * 0.16) : 0;
        const maxCount = Math.max(...histogram.counts);
        const binSpan = (histogram.highEdge - histogram.lowEdge) / histogram.binCount;

        histogram.counts.forEach((count, binIndex) => {{
            if (count === 0) {{
                return;
            }}
            // A floor of 2px keeps rare bins from disappearing next to a dominant one.
            const barHeight = Math.max(2, (count / maxCount) * plotHeight);
            const center = histogram.lowEdge + ((binIndex + 0.5) * binSpan);
            ctx.fillStyle = numericColor(center, histogram.min, histogram.max, palette);
            ctx.fillRect(
                paddingX + (binIndex * slotWidth) + (gap / 2),
                baselineY - barHeight,
                Math.max(1, slotWidth - gap),
                barHeight,
            );
        }});

        ctx.strokeStyle = '#d8dce2';
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.moveTo(paddingX, baselineY + 1);
        ctx.lineTo(canvas.width - paddingX, baselineY + 1);
        ctx.stroke();

        const parts = [
            histogram.total.toLocaleString() + ' value' + (histogram.total === 1 ? '' : 's')
                + ' in ' + histogram.binCount + ' bin' + (histogram.binCount === 1 ? '' : 's'),
            'peak bin ' + maxCount.toLocaleString(),
        ];
        if (histogram.missing > 0) {{
            parts.push(histogram.missing.toLocaleString() + ' with no value');
        }}
        document.getElementById('color-histogram-note').textContent = parts.join(' · ');
    }}

    function updateGradientPreview() {{
        const preview = document.getElementById('color-gradient-preview');
        // Sorted for display only: a stop typed past its neighbor is not
        // reordered until the field commits, and the bar should still read left
        // to right in the meantime.
        const stops = sortedGradientStops(currentGradientStops());
        if (stops.length < MIN_GRADIENT_STOPS) {{
            return;
        }}
        // The bar shares the histogram's axis (the column's data range) rather than
        // spanning the outermost stops, so a stop sits above the values it colors.
        // Percentages are deliberately left unclamped: CSS extends a gradient past
        // 0%/100% and holds the end colors flat outside the stops, which is exactly
        // how numericColor treats values beyond the ends.
        const domain = gradientRangeDomain();
        const parts = stops.map((stop, index) => {{
            const percent = domain
                ? gradientRangePercent(stop.value, domain)
                : ((index / (stops.length - 1)) * 100);
            return stop.color + ' ' + percent + '%';
        }});
        preview.style.background = 'linear-gradient(90deg, ' + parts.join(', ') + ')';
        updateGradientSlider();
    }}

    function renderContinuousPicker(columnName) {{
        const palette = customPalette(columnName);
        const info = colorInfo(columnName);
        const min = info.min;
        const max = info.max;
        const nullInput = document.getElementById('color-null');
        nullInput.value = (palette && palette.nullColor) || cssToHex(DEFAULT_NO_VALUE_COLOR);
        // Binning first: the axis every other part of the panel is drawn against
        // comes from the histogram.
        state.colorHistogram = columnHistogram(columnName);
        // With no stored palette the ramp starts as the two ends of the built-in
        // gradient, which is what the column already looks like on screen.
        state.gradientStops = numericPaletteStops(palette, min, max) || [
            {{value: min, color: cssToHex(numericColor(min, min, max))}},
            {{value: max, color: cssToHex(numericColor(max, min, max))}},
        ];
        // The labels sit under the shared axis, so they report the data extent. The
        // ramp's own ends are the Min/Max rows and the outermost knobs.
        document.getElementById('color-min-label').textContent = formatValue(min);
        document.getElementById('color-max-label').textContent = formatValue(max);
        renderGradientStopRows();
        updateGradientPreview();
        drawColorHistogram(columnName);

        nullInput.oninput = () => {{
            ensureNumericPalette(columnName).nullColor = nullInput.value;
            rebuildNodeColorCache();
            scheduleClusterRender();
            drawColorHistogram(columnName);
        }};
        document.getElementById('color-reset-values').onclick = resetGradientStopValues;
    }}

    function updateColorPaletteControl(columnName, palette, valueCount) {{
        const select = document.getElementById('color-palette');
        const note = document.getElementById('color-palette-note');
        const scheme = palette && palette.scheme ? paletteSchemeByName(palette.scheme) : null;
        const hasCustomColors = !!(palette && palette.colors && Object.keys(palette.colors).length > 0);
        select.value = scheme ? scheme.name : '';
        if (!scheme) {{
            const fallback = paletteSchemeByName(DEFAULT_CATEGORICAL_PALETTE);
            note.textContent = (hasCustomColors
                ? 'Colors were set by hand or loaded from a color table.'
                : 'No palette assigned.') +
                (fallback ? ' Reset to defaults returns this column to ' + fallback.label + '.' : '');
            return;
        }}
        const cycles = valueCount > scheme.colors.length;
        note.textContent = scheme.note + ' ' + scheme.colors.length.toLocaleString() + ' colors for ' +
            valueCount.toLocaleString() + ' value' + (valueCount === 1 ? '' : 's') +
            (cycles ? ', so colors repeat.' : '.');
    }}

    function changeColorPalette() {{
        const columnName = currentColorField();
        if (!columnName) {{
            return;
        }}
        // Every entry in the menu is a named palette now, so there is no
        // pseudo-scheme to special-case; "Reset to defaults" handles going back.
        if (!applyNamedPalette(columnName, document.getElementById('color-palette').value)) {{
            return;
        }}
        rebuildNodeColorCache();
        renderClusterView();
        renderDiscretePicker(columnName);
    }}

    // Step the discrete swatch pager (no-op for numeric/empty columns).
    function changeColorPickerPage(delta) {{
        const columnName = currentColorField();
        const info = colorInfo(columnName);
        if (!columnName || !info || info.type === 'numeric') {{
            return;
        }}
        state.colorPickerPage += delta;
        renderDiscretePicker(columnName);
    }}

    function openColorPicker() {{
        if (!state.bundle) {{
            return;
        }}
        // A fresh open (column switch, reset, TSV load) starts at the first page.
        state.colorPickerPage = 0;
        const overlay = document.getElementById('color-picker-overlay');
        const columnName = currentColorField();
        const info = colorInfo(columnName);
        const discretePanel = document.getElementById('color-picker-discrete');
        const continuousPanel = document.getElementById('color-picker-continuous');
        const nullControl = document.getElementById('color-null-control');
        const emptyNote = document.getElementById('color-picker-empty-note');
        document.getElementById('color-picker-title').textContent =
            columnName ? ('Customize colors: ' + columnName) : 'Customize colors';
        updateCategoricalToggle(columnName);
        if (!columnName || !info) {{
            discretePanel.hidden = true;
            continuousPanel.hidden = true;
            nullControl.hidden = true;
            emptyNote.hidden = false;
        }} else if (info.type === 'numeric') {{
            emptyNote.hidden = true;
            discretePanel.hidden = true;
            continuousPanel.hidden = false;
            nullControl.hidden = false;
            renderContinuousPicker(columnName);
        }} else {{
            emptyNote.hidden = true;
            continuousPanel.hidden = true;
            discretePanel.hidden = false;
            nullControl.hidden = false;
            renderDiscretePicker(columnName);
        }}
        overlay.hidden = false;
    }}

    // The gradient/discrete switch is only meaningful for numeric columns, so it is
    // shown (and its distinct-value hint filled in) only for those.
    function updateCategoricalToggle(columnName) {{
        const control = document.getElementById('color-categorical-control');
        const toggle = document.getElementById('color-as-categorical');
        const note = document.getElementById('color-categorical-note');
        if (!columnName || !columnIsNumericType(columnName)) {{
            control.hidden = true;
            note.hidden = true;
            toggle.checked = false;
            return;
        }}
        control.hidden = false;
        toggle.checked = columnIsCategoricalNumeric(columnName);
        if (toggle.checked) {{
            const total = distinctColumnValues(columnName, 0).total;
            note.textContent = 'Coloring ' + total.toLocaleString() + ' distinct value' +
                (total === 1 ? '' : 's') + ' as separate categories.';
            note.hidden = false;
        }} else {{
            note.hidden = true;
        }}
    }}

    function toggleColorCategorical() {{
        const columnName = currentColorField();
        if (!columnName || !columnIsNumericType(columnName)) {{
            return;
        }}
        setColumnCategorical(columnName, document.getElementById('color-as-categorical').checked);
        rebuildNodeColorCache();
        renderClusterView();
        openColorPicker();
    }}

    function resetColorPicker() {{
        const columnName = currentColorField();
        if (columnName) {{
            delete state.customPalettes[paletteKey(columnName)];
            rebuildNodeColorCache();
            scheduleClusterRender();
            openColorPicker();
        }}
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
        updateThresholdStepButtons();
        drawSplitChart();
        drawClusterView(resetView);
    }}

    // The bundle carries the merge order; everything the threshold slider and the split
    // chart show is derived from it, here, on load. See deriveMergeSeries.
    function normalizeMaxMergeEvents(value) {{
        const cap = Number(value);
        if (!Number.isFinite(cap) || cap < 0) {{ return DEFAULT_WINDOWED_MAX_MERGE_EVENTS; }}
        return Math.floor(cap);
    }}

    // The part of the axis the split chart is currently showing, in threshold units, or
    // null when the whole series is on screen. Clamped by splitChartVisibleWindow against
    // the series' real range, so a stale zoom left over from another bundle cannot narrow
    // the selection to a stretch that no longer exists.
    function currentSplitWindow() {{
        if (!state.series || !state.splitChartZoom) {{ return null; }}
        const range = splitSeriesDataRange(state.series);
        const view = splitChartVisibleWindow(range.min, range.max);
        return view.zoomed ? {{min: view.min, max: view.max}} : null;
    }}

    // Re-select which events the chart plots and the slider stops at, and re-lay-out the
    // slider, without redrawing. Every caller that changes the cap or the window goes
    // through here.
    //
    // `pinnedThreshold` is the cut currently in effect. It has to survive the
    // re-selection: if zooming could filter out the stop the thumb is on, the slider
    // would land on a neighboring stop and the network would silently re-cluster.
    function applyMergeSelection(pinnedThreshold) {{
        selectMergeEvents(state.series, state.maxMergeEvents, currentSplitWindow(),
                          Number.isFinite(pinnedThreshold) ? pinnedThreshold : null);
        state.sliderModel = buildSliderModel(state.series.sliderStops);
        // Point into the model just replaced, so they cannot outlive it.
        state.selectedStop = null;
        state.stopBeforeSplitChartClick = null;

        const slider = document.getElementById('threshold-slider');
        slider.max = String(state.sliderModel.maxPosition);
        slider.disabled = state.sliderModel.stops.length === 0;
        document.getElementById('threshold-min-label').textContent = state.sliderModel.minLabel;
        document.getElementById('threshold-max-label').textContent = state.sliderModel.maxLabel;
        document.getElementById('threshold-input').disabled = state.sliderModel.stops.length === 0;
        const capInput = document.getElementById('split-event-cap');
        capInput.disabled = false;
        capInput.value = String(state.maxMergeEvents);
    }}

    // Put the thumb back on `threshold` under whatever stop list now exists.
    function restoreThreshold(threshold) {{
        snapSliderToStop(Number.isFinite(threshold)
            ? nearestStopForThreshold(threshold)
            : state.sliderModel.stops.find(stop => stop.threshold_value === null));
    }}

    // A new bundle, or a session restore setting the cap. Replays the merge order, which
    // is the expensive half, then selects against the current window.
    //
    // A session restore uses this directly: the threshold_value field that follows it in
    // the registry does the snapping, so snapping here would be undone a moment later.
    function rebuildMergeSeries(maxMergeEvents) {{
        state.maxMergeEvents = normalizeMaxMergeEvents(maxMergeEvents);
        const graph = state.bundle.graph;
        state.series = deriveMergeSeries(
            graph.nodes.length, graph.mst_edges || [], graph.merge_impact_metric);
        applyMergeSelection(null);
    }}

    // The user moved the cap. The stops are re-laid-out, so the old selection's slider
    // position means nothing; the threshold *value* is what survives -- exactly as it
    // does across a session save, and for the same reason. No replay here: the cap only
    // changes which of the already-derived events are selected.
    function applyMaxMergeEvents(maxMergeEvents) {{
        if (!state.bundle || !state.series) {{ return; }}
        const previousThreshold = selectedThresholdValue();
        state.maxMergeEvents = normalizeMaxMergeEvents(maxMergeEvents);
        applyMergeSelection(previousThreshold);
        restoreThreshold(previousThreshold);
        updateSplitEventCount();
        updateThresholdUI(false);
    }}

    // How much of the network's split-event series the chart is actually showing.
    // Without this the blank stretches of a large network's axis are unreadable: the
    // axis and the moving sum span every merge, but the stems are a capped selection,
    // so "nothing plotted here" and "nothing happens here" look identical.
    //
    // The series is derived on load from the bundle's merge order, so the denominator
    // is always known -- on a bundle of any version, including the pre-v6 files whose
    // stored series was already capped when it was written.
    function updateSplitEventCount() {{
        const note = document.getElementById('split-event-count');
        if (!state.bundle || !state.series) {{
            note.textContent = '';
            note.removeAttribute('title');
            return;
        }}
        const plotted = state.series.selectedRows.length;
        const total = state.series.total;
        const cap = state.series.cap;
        const bandPercent = 100 / MERGE_EVENT_DENSITY_BINS;
        const plural = count => (count === 1 ? '' : 's');

        if (plotted >= total) {{
            note.textContent = 'All ' + total.toLocaleString() + ' merge event' + plural(total) + ' plotted.';
        }} else {{
            const backfilled = cap === 0 ? 0 : Math.max(0, plotted - cap);
            note.textContent = plotted.toLocaleString() + ' of ' + total.toLocaleString() +
                ' merge events plotted' +
                (backfilled > 0
                    ? ' \u2014 the strongest ' + cap.toLocaleString() + ' by impact, plus ' +
                      backfilled.toLocaleString() + ' so that every ' + bandPercent + '% of the axis'
                      + ' with an event to show has one.'
                    : ' \u2014 the strongest by impact.');
        }}
        note.title = plotted >= total
            ? 'Every split event in this network is drawn.'
            : 'The split chart draws a capped selection of the network\u2019s ' +
              total.toLocaleString() + ' split events, set by the "Split events" box ' +
              '(0 plots them all). The axis and the moving-sum line always span every ' +
              'event, so a stretch with no stems is a stretch whose events were too ' +
              'small to make the cut, not necessarily a quiet one.';
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

    // The next stop in `delta`'s direction that actually names a different threshold,
    // or -1 when there is none. What an arrow press promises is a threshold the view
    // has not just been at; two stops naming one threshold describe the same cut, so
    // landing on the second would move the slider and change nothing on screen. The
    // event rows the stops are built from are grouped by exact threshold today, so this
    // skips nothing -- it is here because the promise belongs to the arrows rather than
    // to an invariant held one module away in ssn_hierarchy.
    //
    // (The arrows' own stuck-looking behavior came from somewhere else: many stops
    // share one slider position. See currentSliderStop.)
    function nextDistinctStopIndex(stops, index, delta) {{
        const current = stops[index].threshold_value;
        for (let probe = index + delta; probe >= 0 && probe < stops.length; probe += delta) {{
            if (stops[probe].threshold_value !== current) {{
                return probe;
            }}
        }}
        return -1;
    }}

    // Walk the slider's stop list. Those stops are the only thresholds the view can
    // actually take -- every one is an edge weight at which the graph splits -- so
    // stepping to the next distinct one is the only step size that always lands
    // somewhere new. `delta` is +1 toward higher thresholds (rightward on the split
    // plot, ending at the infinity stop) and -1 toward lower ones.
    function stepThreshold(delta) {{
        if (!state.bundle || !state.sliderModel) {{
            return;
        }}
        const stops = state.sliderModel.stops;
        if (stops.length === 0) {{
            return;
        }}
        // currentSliderStop() returns an element of `stops`, so identity search is safe.
        const index = stops.indexOf(currentSliderStop());
        const target = index < 0 ? 0 : nextDistinctStopIndex(stops, index, delta);
        if (target < 0) {{
            return;
        }}
        snapSliderToStop(stops[target]);
        // Keep the user's pan/zoom, like releasing the slider does; stepping through
        // stops to watch one cluster break up is the point of these buttons.
        scheduleThresholdUI(false);
    }}

    // The ends of the stop list are dead ends, so say so rather than no-op silently.
    // Asked the same question stepThreshold answers, so an arrow is enabled exactly
    // when pressing it would move the threshold somewhere new.
    function updateThresholdStepButtons() {{
        const stops = state.sliderModel ? state.sliderModel.stops : [];
        const index = state.bundle && stops.length > 0 ? stops.indexOf(currentSliderStop()) : -1;
        document.getElementById('threshold-step-down').disabled =
            index < 0 || nextDistinctStopIndex(stops, index, -1) < 0;
        document.getElementById('threshold-step-up').disabled =
            index < 0 || nextDistinctStopIndex(stops, index, 1) < 0;
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
        if (state.latticeNodeItem) {{
            const nodeIndex = latticeNodeAtWorld(worldPoint.x, worldPoint.y);
            return nodeIndex >= 0 ? nodeIndex : null;
        }}
        for (let index = state.visibleLayout.length - 1; index >= 0; index--) {{
            const item = state.visibleLayout[index];
            if (item.shape === 'rect') {{
                if (worldPoint.x < item.x0 || worldPoint.x > item.x1 ||
                    worldPoint.y < item.y0 || worldPoint.y > item.y1) {{
                    continue;
                }}
            }} else {{
                const dx = worldPoint.x - item.x;
                const dy = worldPoint.y - item.y;
                if ((dx * dx) + (dy * dy) > item.radius * item.radius) {{
                    continue;
                }}
            }}
            const component = state.bundle.graph.hierarchy.nodes[item.componentId];
            const dotLayout = componentMemberLayout(component, item);
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
        if (state.latticeNodeItem) {{
            const nodeIndex = latticeNodeAtWorld(worldPoint.x, worldPoint.y);
            if (nodeIndex < 0) {{ return null; }}
            const itemIndex = state.latticeNodeItem[nodeIndex];
            return itemIndex >= 0 ? state.visibleLayout[itemIndex] : null;
        }}
        for (let index = state.visibleLayout.length - 1; index >= 0; index--) {{
            const item = state.visibleLayout[index];
            if (item.shape === 'rect') {{
                if (worldPoint.x >= item.x0 && worldPoint.x <= item.x1 &&
                    worldPoint.y >= item.y0 && worldPoint.y <= item.y1) {{
                    return item;
                }}
            }} else {{
                const dx = worldPoint.x - item.x;
                const dy = worldPoint.y - item.y;
                if ((dx * dx) + (dy * dy) <= item.radius * item.radius) {{
                    return item;
                }}
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

    // Axis-aligned overlap between a screen-space rectangle [ax0,ay0]-[ax1,ay1] (corners in
    // any order) and a selection box {{left, top, width, height}}.
    function rectIntersectsRect(ax0, ay0, ax1, ay1, rect) {{
        const left = Math.min(ax0, ax1);
        const right = Math.max(ax0, ax1);
        const top = Math.min(ay0, ay1);
        const bottom = Math.max(ay0, ay1);
        return left <= (rect.left + rect.width) && right >= rect.left &&
            top <= (rect.top + rect.height) && bottom >= rect.top;
    }}

    function itemIntersectsSelectionBox(item, rect) {{
        if (item.shape === 'rect' || item.shape === 'lattice') {{
            const p0 = worldToScreenPoint(item.x0, item.y0);
            const p1 = worldToScreenPoint(item.x1, item.y1);
            return rectIntersectsRect(p0.x, p0.y, p1.x, p1.y, rect);
        }}
        const screenPoint = worldToScreenPoint(item.x, item.y);
        const screenRadius = item.radius * state.viewTransform.scale;
        return circleIntersectsRect(screenPoint.x, screenPoint.y, screenRadius, rect);
    }}

    function selectComponentsInBox(rect, deselect) {{
        let changed = false;
        state.visibleLayout.forEach(item => {{
            if (!itemIntersectsSelectionBox(item, rect)) {{
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
            if (!itemIntersectsSelectionBox(item, rect)) {{
                return;
            }}
            const component = state.bundle.graph.hierarchy.nodes[item.componentId];
            const dotLayout = componentMemberLayout(component, item);
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
    document.getElementById('split-event-cap').addEventListener('change', event => {{
        // Re-deriving the whole series is one pass over the MST edges, so this is a
        // plain 'change' handler rather than anything debounced or deferred.
        applyMaxMergeEvents(event.target.value);
    }});
    document.getElementById('threshold-input').addEventListener('change', jumpToThresholdValue);
    document.getElementById('threshold-input').addEventListener('keydown', event => {{
        // Up/Down step stops -- the behavior a number input's spinner would have had
        // if its step were the distance to the next split. Left/Right are left alone
        // so they still move the caret inside the field.
        if (event.key === 'ArrowUp' || event.key === 'ArrowDown') {{
            event.preventDefault();
            stepThreshold(event.key === 'ArrowUp' ? 1 : -1);
            return;
        }}
        if (event.key !== 'Enter') {{
            return;
        }}
        event.preventDefault();
        jumpToThresholdValue();
    }});
    document.getElementById('threshold-step-down').addEventListener('click', () => stepThreshold(-1));
    document.getElementById('threshold-step-up').addEventListener('click', () => stepThreshold(1));
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
    document.getElementById('render-cluster-bounds').addEventListener('change', () => renderClusterView());
    document.getElementById('render-nodes').addEventListener('change', () => renderClusterView());
    document.getElementById('leaf-pruning-only').addEventListener('change', () => {{
        updateCollapseLongPathsControl();
        scheduleThresholdUI(true);
    }});
    document.getElementById('collapse-long-paths').addEventListener('change', () => {{
        scheduleThresholdUI(true);
    }});
    document.getElementById('color-by').addEventListener('change', () => {{
        rebuildNodeColorCache();
        renderClusterView();
        if (colorPickerIsOpen()) {{
            openColorPicker();
        }}
    }});
    document.getElementById('label-by').addEventListener('change', () => renderClusterView());
    document.getElementById('show-node-counts').addEventListener('change', () => renderClusterView());
    document.getElementById('show-edge-scores').addEventListener('change', () => renderClusterView());
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
    document.getElementById('export-png').addEventListener('click', exportClusterPNG);
    document.getElementById('export-split-png').addEventListener('click', exportSplitChartPNG);
    document.getElementById('export-split-svg').addEventListener('click', exportSplitChartSVG);
    document.getElementById('export-svg').addEventListener('click', exportClusterSVG);
    document.getElementById('customize-colors').addEventListener('click', openColorPicker);
    document.getElementById('save-session').addEventListener('click', saveSessionFile);
    document.getElementById('rename-network').addEventListener('click', openRenameNetworkDialog);
    document.getElementById('save-extraction').addEventListener('click', saveExtractionFile);
    document.getElementById('color-picker-close').addEventListener('click', closeColorPicker);
    document.getElementById('color-picker-overlay').addEventListener('click', event => {{
        if (event.target === event.currentTarget) {{
            closeColorPicker();
        }}
    }});
    document.addEventListener('keydown', event => {{
        if (event.key === 'Escape' && colorPickerIsOpen()) {{
            closeColorPicker();
        }}
    }});
    document.getElementById('color-palette').addEventListener('change', changeColorPalette);
    document.getElementById('color-as-categorical').addEventListener('change', toggleColorCategorical);
    document.getElementById('color-picker-reset').addEventListener('click', resetColorPicker);
    document.getElementById('color-picker-prev').addEventListener('click', () => changeColorPickerPage(-1));
    document.getElementById('color-picker-next').addEventListener('click', () => changeColorPickerPage(1));
    document.getElementById('load-color-table').addEventListener('click', () => document.getElementById('color-table-file').click());
    document.getElementById('color-table-file').addEventListener('change', async event => {{
        const file = event.target.files && event.target.files[0];
        if (!file) {{
            return;
        }}
        try {{
            const text = await file.text();
            loadColorTableText(text);
        }} catch (error) {{
            console.error(error);
            setStatus('Failed to load color table: ' + error.message);
        }} finally {{
            event.target.value = '';
        }}
    }});
    document.getElementById('save-color-table').addEventListener('click', saveColorTableTSV);
    document.getElementById('export-legend-svg').addEventListener('click', exportLegendSVG);
    document.getElementById('export-legend-png').addEventListener('click', exportLegendPNG);
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
    document.getElementById('focus-selection').addEventListener('click', () => {{
        if (!focusSelection()) {{
            setStatus('No selected nodes are visible at the current threshold and minimum cluster size.');
            return;
        }}
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
    setupSelectionPresets();
    setupMetadataEditing();
    setupColumnCharts();
    setupSplitChartHover();
    setupNameDialog();
    setupGradientRangeSlider();
    setupGradientStopEditing();
    populateColorPaletteMenu();
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
    title: str = VIEWER_APP_NAME,
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