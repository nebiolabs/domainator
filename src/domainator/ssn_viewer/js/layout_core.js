/*
 * Layout math shared by the layout Web Worker and the main-thread fallback.
 *
 * Pure functions only -- no `state`, no DOM -- because this same source is
 * compiled into the Worker, where neither exists. It goes both into the worker
 * source (with layout_worker.js) and into the page's one `<script>` scope (its
 * `// @include` in viewer.js). Function declarations hoist, so the order below
 * is for reading, not for dependency.
 *
 * This used to exist twice -- once in minified form for the worker and once
 * inline in the page -- and the two copies had already drifted. Keep it single-copy.
 */

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
// Arc each bubble is allotted along its ring, over and above its own diameter. Matches the
// bubblePadding that refineLayoutGeometry separates overlapping bubbles to.
const RADIAL_SEED_ARC_PAD = 17;

function radialTreeSeed(componentIds, adjacency, hierarchyNodes, ringGap) {
    const rootId = treeCenter(componentIds, adjacency, hierarchyNodes);
    const tree = rootedTree(rootId, adjacency, hierarchyNodes);
    const radii = new Map(componentIds.map(id => [id, componentRadiusForSize(hierarchyNodes[id].size)]));
    // Concentric rings, one radius per depth, with the step from each ring to the next sized
    // by the biggest bubbles actually on those two depths.
    //
    // Both halves of that matter. Spacing every ring by the largest bubble anywhere in the
    // component (what this used to do) let one big cluster push apart rings that only join
    // singletons, and anchorStrength then held the simulation in that spread-out arrangement.
    // But giving each parent->child pair its own step is worse: nodes at the same depth end
    // up at different radii, the wedges stop being disjoint annular sectors, and subtrees
    // interleave -- on a network with one 14k-member cluster that produced 216 edge crossings
    // and 11 overlapping bubbles. Per depth keeps the rings concentric, so sibling subtrees
    // stay in their own sectors, while still being sized by what is actually there.
    const gap = Math.max(ringGap || 0, 60);
    const depthByNode = new Map([[rootId, 0]]);
    tree.order.forEach(nodeId => {
        (tree.children.get(nodeId) || []).forEach(c => depthByNode.set(c, (depthByNode.get(nodeId) || 0) + 1));
    });
    let maxDepth = 0;
    componentIds.forEach(id => { maxDepth = Math.max(maxDepth, depthByNode.get(id) || 0); });
    const ringMaxRadius = new Array(maxDepth + 1).fill(0);
    componentIds.forEach(id => {
        const depth = depthByNode.get(id) || 0;
        ringMaxRadius[depth] = Math.max(ringMaxRadius[depth], radii.get(id) || 10);
    });
    // A ring also has to be long enough to hold everything standing on it. Sizing rings by
    // radius alone let a cluster with hundreds of neighbours seed them shoulder to shoulder
    // on a circle with no room for them, and the collision forces then blew the ring apart
    // into a one-sided spike -- on a network whose main cluster had 351 neighbours, all 351
    // ended up in a single beam off one side of it.
    const ringArcNeed = new Array(maxDepth + 1).fill(0);
    componentIds.forEach(id => {
        const depth = depthByNode.get(id) || 0;
        ringArcNeed[depth] += (2 * (radii.get(id) || 10)) + RADIAL_SEED_ARC_PAD;
    });
    const depthRadius = new Array(maxDepth + 1).fill(0);
    for (let depth = 1; depth <= maxDepth; depth++) {
        depthRadius[depth] = Math.max(
            depthRadius[depth - 1]
                + Math.max(gap, ringMaxRadius[depth - 1] + ringMaxRadius[depth] + RADIAL_SEED_RING_PAD),
            ringArcNeed[depth] / (2 * Math.PI));
    }
    const ringRadius = new Map(componentIds.map(id => [id, depthRadius[depthByNode.get(id) || 0]]));
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
