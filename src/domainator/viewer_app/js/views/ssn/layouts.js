    function activeClusterAssignments(activeClusterIds) {
        const assignments = new Array(state.bundle.nodes.length).fill(-1);
        activeClusterIds.forEach(componentId => {
            const members = componentMembers(componentId);
            members.forEach(nodeIndex => {
                assignments[nodeIndex] = componentId;
            });
        });
        return assignments;
    }

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
    function collapseLongPathLinks(links, visibleSet, minClusterSize) {
        const hierarchyNodes = state.bundle.views.ssn.hierarchy.nodes;
        const adjacency = new Map();
        links.forEach(link => {
            if (!visibleSet.has(link.sourceId) || !visibleSet.has(link.targetId)) { return; }
            if (!adjacency.has(link.sourceId)) { adjacency.set(link.sourceId, []); }
            if (!adjacency.has(link.targetId)) { adjacency.set(link.targetId, []); }
            adjacency.get(link.sourceId).push({other: link.targetId, weight: link.weight});
            adjacency.get(link.targetId).push({other: link.sourceId, weight: link.weight});
        });

        const isPassThrough = componentId => visibleSet.has(componentId)
            && (adjacency.get(componentId) || []).length === 2
            && hierarchyNodes[componentId].size < minClusterSize;

        // Walked in ascending cluster id so the result does not depend on link order.
        const starts = Array.from(adjacency.keys()).filter(isPassThrough).sort((leftId, rightId) => leftId - rightId);
        const consumed = new Set();
        const collapsedLinks = new Map();
        let collapsedPaths = 0;

        starts.forEach(startId => {
            if (consumed.has(startId)) { return; }
            consumed.add(startId);
            const chain = [startId];
            const ends = [];
            let weakest = Infinity;
            // Two steps out of a degree-2 cluster: walk each way to the first cluster
            // that is not itself pass-through.
            (adjacency.get(startId) || []).forEach(step => {
                let previousId = startId;
                let currentId = step.other;
                weakest = Math.min(weakest, step.weight);
                while (isPassThrough(currentId) && !consumed.has(currentId)) {
                    consumed.add(currentId);
                    chain.push(currentId);
                    const next = (adjacency.get(currentId) || []).find(edge => edge.other !== previousId);
                    if (!next) { break; }
                    weakest = Math.min(weakest, next.weight);
                    previousId = currentId;
                    currentId = next.other;
                }
                ends.push(isPassThrough(currentId) ? null : currentId);
            });

            chain.forEach(componentId => visibleSet.delete(componentId));
            const leftEnd = ends[0];
            const rightEnd = ends[1];
            if (leftEnd === null || rightEnd === null || leftEnd === rightEnd) {
                // Only reachable if the forest assumption above does not hold. The chain
                // is gone either way; drawing a self-loop or a dangling edge is worse
                // than drawing nothing.
                return;
            }
            collapsedPaths += 1;
            const sourceId = Math.min(leftEnd, rightEnd);
            const targetId = Math.max(leftEnd, rightEnd);
            const key = sourceId + ':' + targetId;
            const existing = collapsedLinks.get(key);
            if (existing) {
                existing.weight = Math.min(existing.weight, weakest);
                existing.collapsed += chain.length;
                return;
            }
            collapsedLinks.set(key, {sourceId, targetId, weight: weakest, collapsed: chain.length});
        });

        const kept = links.filter(link => visibleSet.has(link.sourceId) && visibleSet.has(link.targetId));
        return {
            links: kept.concat(Array.from(collapsedLinks.values())),
            collapsedPaths,
            collapsedClusters: consumed.size,
        };
    }

    function mstLinksForActiveClusters(activeClusterIds, minClusterSize, leafPruningOnly, collapseLongPaths) {
        const assignments = activeClusterAssignments(activeClusterIds);
        const linkMap = new Map();

        state.bundle.views.ssn.mst_edges.forEach(edge => {
            const sourceComponentId = assignments[edge[0]];
            const targetComponentId = assignments[edge[1]];
            if (sourceComponentId < 0 || targetComponentId < 0 || sourceComponentId === targetComponentId) {
                return;
            }
            const leftId = sourceComponentId < targetComponentId ? sourceComponentId : targetComponentId;
            const rightId = sourceComponentId < targetComponentId ? targetComponentId : sourceComponentId;
            const key = leftId + ':' + rightId;
            if (!linkMap.has(key)) {
                linkMap.set(key, {sourceId: leftId, targetId: rightId, weight: edge[2]});
            }
        });

        if (minClusterSize <= 1) {
            return {
                visibleIds: [...activeClusterIds],
                links: Array.from(linkMap.values()),
                hiddenNodes: 0,
                collapsedPaths: 0,
                collapsedClusters: 0,
            };
        }

        const allLinks = Array.from(linkMap.values());
        if (!leafPruningOnly) {
            const visibleIds = activeClusterIds.filter(componentId => state.bundle.views.ssn.hierarchy.nodes[componentId].size >= minClusterSize);
            const visibleSet = new Set(visibleIds);
            const hiddenNodes = activeClusterIds.reduce((sum, componentId) => {
                if (visibleSet.has(componentId)) {
                    return sum;
                }
                return sum + state.bundle.views.ssn.hierarchy.nodes[componentId].size;
            }, 0);
            return {
                visibleIds,
                links: allLinks.filter(link => visibleSet.has(link.sourceId) && visibleSet.has(link.targetId)),
                hiddenNodes,
                collapsedPaths: 0,
                collapsedClusters: 0,
            };
        }

        const visibleSet = new Set(activeClusterIds);
        const adjacency = new Map();
        activeClusterIds.forEach(componentId => adjacency.set(componentId, []));
        allLinks.forEach(link => {
            adjacency.get(link.sourceId)?.push(link.targetId);
            adjacency.get(link.targetId)?.push(link.sourceId);
        });

        const degree = new Map();
        activeClusterIds.forEach(componentId => degree.set(componentId, (adjacency.get(componentId) || []).length));
        let removedAny = true;
        while (removedAny) {
            removedAny = false;
            const pruneIds = [];
            visibleSet.forEach(componentId => {
                const size = state.bundle.views.ssn.hierarchy.nodes[componentId].size;
                const currentDegree = degree.get(componentId) || 0;
                if (currentDegree <= 1 && size < minClusterSize) {
                    pruneIds.push(componentId);
                }
            });
            if (pruneIds.length === 0) {
                break;
            }
            removedAny = true;
            pruneIds.forEach(componentId => {
                visibleSet.delete(componentId);
            });
            pruneIds.forEach(componentId => {
                (adjacency.get(componentId) || []).forEach(neighborId => {
                    if (!visibleSet.has(neighborId)) {
                        return;
                    }
                    degree.set(neighborId, Math.max(0, (degree.get(neighborId) || 0) - 1));
                });
                degree.set(componentId, 0);
            });
        }

        // Runs after pruning and narrows visibleSet further, so visibleIds and
        // hiddenNodes below already account for the contracted clusters.
        const collapse = collapseLongPaths
            ? collapseLongPathLinks(allLinks, visibleSet, minClusterSize)
            : null;

        const visibleIds = activeClusterIds.filter(componentId => visibleSet.has(componentId));
        const hiddenNodes = activeClusterIds.reduce((sum, componentId) => {
            if (visibleSet.has(componentId)) {
                return sum;
            }
            return sum + state.bundle.views.ssn.hierarchy.nodes[componentId].size;
        }, 0);

        return {
            visibleIds,
            links: collapse
                ? collapse.links
                : allLinks.filter(link => visibleSet.has(link.sourceId) && visibleSet.has(link.targetId)),
            hiddenNodes,
            collapsedPaths: collapse ? collapse.collapsedPaths : 0,
            collapsedClusters: collapse ? collapse.collapsedClusters : 0,
        };
    }

    function componentDotCount(componentSize) {
        // Every node in the component is drawn as a dot.
        return componentSize;
    }

    function componentDotRadius(componentSize, bubbleRadius) {
        return 1.95;
    }

    function componentDotOffset(sampleIndex, sampleCount, packingRadius) {
        if (sampleCount <= 1 || packingRadius <= 0) {
            return {x: 0, y: 0};
        }

        if (sampleCount <= 6) {
            const ringRadius = sampleCount === 2
                ? packingRadius * 0.82
                : sampleCount <= 4
                    ? packingRadius * 0.72
                    : packingRadius * 0.78;
            const angleOffset = sampleCount === 4 ? Math.PI / 4 : -Math.PI / 2;
            const angle = angleOffset + (sampleIndex * ((Math.PI * 2) / sampleCount));
            return {
                x: Math.cos(angle) * ringRadius,
                y: Math.sin(angle) * ringRadius,
            };
        }

        const angle = sampleIndex * 2.399963229728653;
        const distance = Math.sqrt((sampleIndex + 0.5) / sampleCount) * packingRadius;
        return {
            x: Math.cos(angle) * distance,
            y: Math.sin(angle) * distance,
        };
    }


    function sampledMembersForComponent(componentId, sampleCount) {
        const hierarchy = state.bundle.views.ssn.hierarchy;
        const component = hierarchy.nodes[componentId];
        const sampledMembers = [];
        for (let sampleIndex = 0; sampleIndex < sampleCount; sampleIndex++) {
            const leafPosition = component.leaf_start + Math.floor(sampleIndex * component.leaf_count / sampleCount);
            sampledMembers.push({
                nodeIndex: hierarchy.leaf_order[leafPosition],
                leafPosition,
            });
        }
        return sampledMembers;
    }

    function dotLayoutCacheKey(componentId, sampleCount, arrangement, spacingKey = 'base') {
        return arrangement + ':' + componentId + ':' + sampleCount + ':' + spacingKey;
    }

    function radialDotPositions(sampleCount) {
        const positions = [];
        for (let sampleIndex = 0; sampleIndex < sampleCount; sampleIndex++) {
            const offset = componentDotOffset(sampleIndex, sampleCount, 1);
            positions.push({x: offset.x, y: offset.y, sampleIndex});
        }
        return positions;
    }

    function radialDotLayout(sampledMembers) {
        const positions = radialDotPositions(sampledMembers.length);
        return sampledMembers.map((member, sampleIndex) => {
            const position = positions[sampleIndex];
            return {memberIndex: member.nodeIndex, x: position.x, y: position.y};
        });
    }

    const phyllotaxisCoordCache = new Map();

    function phyllotaxisFor(n) {
        // Flat coordinate arrays for the phyllotaxis (sunflower) positions of size n.
        // Positions for a given n are identical across components, so cache by n. The
        // grouped layout only reads these, never mutates them, so sharing is safe.
        let coords = phyllotaxisCoordCache.get(n);
        if (coords === undefined) {
            const points = radialDotPositions(n);
            const posX = new Float64Array(n);
            const posY = new Float64Array(n);
            for (let i = 0; i < n; i++) {
                posX[i] = points[i].x;
                posY[i] = points[i].y;
            }
            coords = {posX: posX, posY: posY};
            phyllotaxisCoordCache.set(n, coords);
        }
        return coords;
    }

    function selectByCoord(idx, start, end, leftCount, posX, posY, useY) {
        // In-place Hoare quickselect: rearrange idx[start, end) so the leftCount entries
        // with the smallest coordinate on the chosen axis occupy [start, start+leftCount).
        // Coordinates are read indirectly via the shared position arrays. Deterministic
        // (median-of-three pivot); no allocation, no closures.
        const target = start + leftCount;
        let lo = start;
        let hi = end;
        while (hi - lo > 1) {
            const a = lo;
            const b = (lo + hi) >> 1;
            const c = hi - 1;
            const va = useY ? posY[idx[a]] : posX[idx[a]];
            const vb = useY ? posY[idx[b]] : posX[idx[b]];
            const vc = useY ? posY[idx[c]] : posX[idx[c]];
            let pivot;
            if (va < vb) {
                pivot = vb < vc ? vb : (va < vc ? vc : va);
            } else {
                pivot = va < vc ? va : (vb < vc ? vc : vb);
            }
            let i = lo;
            let j = hi - 1;
            while (i <= j) {
                while ((useY ? posY[idx[i]] : posX[idx[i]]) < pivot) { i++; }
                while ((useY ? posY[idx[j]] : posX[idx[j]]) > pivot) { j--; }
                if (i <= j) {
                    const t = idx[i];
                    idx[i] = idx[j];
                    idx[j] = t;
                    i++;
                    j--;
                }
            }
            if (target <= j) {
                hi = j + 1;
            } else if (target >= i) {
                lo = i;
            } else {
                break;
            }
        }
    }

    function regionLongerAxisIsY(idx, pStart, pEnd, posX, posY) {
        let minX = Infinity;
        let maxX = -Infinity;
        let minY = Infinity;
        let maxY = -Infinity;
        for (let k = pStart; k < pEnd; k++) {
            const pi = idx[k];
            const xx = posX[pi];
            const yy = posY[pi];
            if (xx < minX) { minX = xx; }
            if (xx > maxX) { maxX = xx; }
            if (yy < minY) { minY = yy; }
            if (yy > maxY) { maxY = yy; }
        }
        return (maxY - minY) > (maxX - minX);
    }

    function regionPrincipalAxis(idx, pStart, pEnd, posX, posY) {
        // Principal (largest-variance) axis of the points in [pStart, pEnd), via the
        // eigenvector of the 2x2 covariance matrix. Cutting perpendicular to this axis
        // keeps each child region rounder than the bounding-box axis when the region is
        // diagonally elongated. Used only when "reduce elongation" is enabled.
        const m = pEnd - pStart;
        let mx = 0;
        let my = 0;
        for (let k = pStart; k < pEnd; k++) {
            const pi = idx[k];
            mx += posX[pi];
            my += posY[pi];
        }
        mx /= m;
        my /= m;
        let sxx = 0;
        let sxy = 0;
        let syy = 0;
        for (let k = pStart; k < pEnd; k++) {
            const pi = idx[k];
            const dx = posX[pi] - mx;
            const dy = posY[pi] - my;
            sxx += dx * dx;
            sxy += dx * dy;
            syy += dy * dy;
        }
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
        if ((ax1 * ax1 + ay1 * ay1) >= (ax2 * ax2 + ay2 * ay2)) { ax = ax1; ay = ay1; } else { ax = ax2; ay = ay2; }
        const norm = Math.hypot(ax, ay);
        if (norm < 1e-12) { return {ax: 1, ay: 0}; }   // isotropic; any axis is fine
        return {ax: ax / norm, ay: ay / norm};
    }

    function selectByProjection(idx, start, end, leftCount, posX, posY, ax, ay) {
        // Like selectByCoord but partitions by projection onto an arbitrary axis (ax, ay)
        // instead of a single coordinate. In-place Hoare quickselect, median-of-three.
        const target = start + leftCount;
        let lo = start;
        let hi = end;
        while (hi - lo > 1) {
            const a = lo;
            const b = (lo + hi) >> 1;
            const c = hi - 1;
            const va = posX[idx[a]] * ax + posY[idx[a]] * ay;
            const vb = posX[idx[b]] * ax + posY[idx[b]] * ay;
            const vc = posX[idx[c]] * ax + posY[idx[c]] * ay;
            let pivot;
            if (va < vb) {
                pivot = vb < vc ? vb : (va < vc ? vc : va);
            } else {
                pivot = va < vc ? va : (vb < vc ? vc : vb);
            }
            let i = lo;
            let j = hi - 1;
            while (i <= j) {
                while ((posX[idx[i]] * ax + posY[idx[i]] * ay) < pivot) { i++; }
                while ((posX[idx[j]] * ax + posY[idx[j]] * ay) > pivot) { j--; }
                if (i <= j) {
                    const t = idx[i];
                    idx[i] = idx[j];
                    idx[j] = t;
                    i++;
                    j--;
                }
            }
            if (target <= j) {
                hi = j + 1;
            } else if (target >= i) {
                lo = i;
            } else {
                break;
            }
        }
    }

    function collectMajorChildren(nodes, sampledMembers, nodeId, lo, hi, frac) {
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
        while (true) {
            const node = cur >= 0 ? nodes[cur] : null;
            const m = curHi - curLo;
            if (!node || node.kind === 'leaf' || m <= 1) {
                head.push({node: cur, lo: curLo, hi: curHi});
                break;
            }
            const leftNode = nodes[node.left];
            const leftBoundary = leftNode.leaf_start + leftNode.leaf_count;
            // lowerBound in [curLo, curHi): first member whose leafPosition >= leftBoundary.
            let sLo = curLo;
            let sHi = curHi;
            while (sLo < sHi) {
                const mid = (sLo + sHi) >> 1;
                if (sampledMembers[mid].leafPosition >= leftBoundary) { sHi = mid; } else { sLo = mid + 1; }
            }
            const splitMember = sLo;
            const leftCount = splitMember - curLo;
            const rightCount = curHi - splitMember;
            if (leftCount <= 0 || rightCount <= 0) {
                // Sampling can't resolve this boundary; keep cur as one atomic child.
                head.push({node: cur, lo: curLo, hi: curHi});
                break;
            }
            if (Math.min(leftCount, rightCount) >= frac * m) {
                // Balanced split: emit both children (left precedes right) and stop.
                head.push({node: node.left, lo: curLo, hi: splitMember});
                head.push({node: node.right, lo: splitMember, hi: curHi});
                break;
            }
            // Lopsided: peel the smaller sibling, descend into the bigger.
            if (leftCount <= rightCount) {
                head.push({node: node.left, lo: curLo, hi: splitMember});   // left precedes -> head
                cur = node.right; curLo = splitMember;
            } else {
                tail.push({node: node.right, lo: splitMember, hi: curHi});   // right follows -> tail
                cur = node.left; curHi = splitMember;
            }
        }
        // tail holds right-peels in decreasing leaf_start; reverse to increasing and append.
        for (let k = tail.length - 1; k >= 0; k--) { head.push(tail[k]); }
        return head;
    }

    function partitionAmongChildren(children, idx, posX, posY, pStart0, pEnd0, superStack, usePca) {
        // Partition the point slice [pStart0, pEnd0) among the ordered child list by
        // recursively bisecting the LIST near its member-count median. Every list cut
        // falls between two children (a real sibling boundary), so no subtree is split.
        // Points are split by the longer bounding-box axis, or (usePca) by the region's
        // principal axis to reduce elongation. Emits a super-node frame per child. Own
        // explicit stack.
        const partStack = [{ci: 0, cj: children.length, pStart: pStart0, pEnd: pEnd0}];
        while (partStack.length > 0) {
            const f = partStack.pop();
            const ci = f.ci;
            const cj = f.cj;
            const ps = f.pStart;
            const pe = f.pEnd;
            if (cj - ci === 1) {
                const ch = children[ci];
                superStack.push({lo: ch.lo, hi: ch.hi, pStart: ps, pEnd: pe, nodeId: ch.node});
                continue;
            }
            const mlo = children[ci].lo;
            const mhi = children[cj - 1].hi;
            const target = mlo + ((mhi - mlo) >> 1);
            // First list index in (ci, cj) whose child.lo >= target (balanced member split).
            let aLo = ci + 1;
            let aHi = cj;
            while (aLo < aHi) {
                const mid = (aLo + aHi) >> 1;
                if (children[mid].lo >= target) { aHi = mid; } else { aLo = mid + 1; }
            }
            let s = aLo;
            if (s <= ci) { s = ci + 1; }
            if (s >= cj) { s = cj - 1; }
            const leftMembers = children[s].lo - mlo;
            const rightMembers = mhi - children[s].lo;
            if (usePca) {
                const axis = regionPrincipalAxis(idx, ps, pe, posX, posY);
                selectByProjection(idx, ps, pe, leftMembers, posX, posY, axis.ax, axis.ay);
            } else {
                const useY = regionLongerAxisIsY(idx, ps, pe, posX, posY);
                selectByCoord(idx, ps, pe, leftMembers, posX, posY, useY);
            }
            const pSplit = ps + leftMembers;
            // Push the larger group first so the smaller is processed first (bounded stack).
            if (leftMembers >= rightMembers) {
                partStack.push({ci: ci, cj: s, pStart: ps, pEnd: pSplit});
                partStack.push({ci: s, cj: cj, pStart: pSplit, pEnd: pe});
            } else {
                partStack.push({ci: s, cj: cj, pStart: pSplit, pEnd: pe});
                partStack.push({ci: ci, cj: s, pStart: ps, pEnd: pSplit});
            }
        }
    }

    function groupedDotLayout(componentId, sampledMembers, usePca) {
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
        if (n <= 6) {
            // Keep the tuned ring placement for tiny clusters.
            return radialDotLayout(sampledMembers);
        }
        const nodes = state.bundle.views.ssn.hierarchy.nodes;
        const phy = phyllotaxisFor(n);
        const posX = phy.posX;
        const posY = phy.posY;
        const FRAC = 0.18;

        const idx = new Int32Array(n);
        for (let i = 0; i < n; i++) { idx[i] = i; }
        const out = new Array(n);

        // Super-node frame invariant: pEnd - pStart === hi - lo (members paired one-to-one
        // with the points allocated to this subtree). idx is partitioned in place.
        const superStack = [{lo: 0, hi: n, pStart: 0, pEnd: n, nodeId: componentId}];
        while (superStack.length > 0) {
            const fr = superStack.pop();
            const lo = fr.lo;
            const hi = fr.hi;
            const pStart = fr.pStart;
            const pEnd = fr.pEnd;
            const m = hi - lo;
            if (m <= 0) { continue; }
            if (m === 1) {
                const p = idx[pStart];
                out[lo] = {memberIndex: sampledMembers[lo].nodeIndex, x: posX[p], y: posY[p]};
                continue;
            }
            const children = collectMajorChildren(nodes, sampledMembers, fr.nodeId, lo, hi, FRAC);
            if (children.length <= 1) {
                // Atomic block (leaf with m>1, or unresolved sampling): assign in order.
                for (let k = 0; k < m; k++) {
                    const p = idx[pStart + k];
                    out[lo + k] = {memberIndex: sampledMembers[lo + k].nodeIndex, x: posX[p], y: posY[p]};
                }
                continue;
            }
            partitionAmongChildren(children, idx, posX, posY, pStart, pEnd, superStack, usePca);
        }
        return out;
    }

    function normalizedComponentDotLayout(componentId, sampleCount, minimumDistance = 0) {
        // Dots are always laid out with the grouped (chain-contraction) algorithm; the
        // arrangement dropdown was removed. The optional PCA axis reduces elongation.
        const usePca = reduceElongationEnabled();
        const spacingKey = usePca ? 'grouped-pca' : 'grouped-phyllotaxis';
        const cacheKey = dotLayoutCacheKey(componentId, sampleCount, 'grouped', spacingKey);
        const cached = state.dotLayoutCache.get(cacheKey);
        if (cached) {
            return cached;
        }

        const sampledMembers = sampledMembersForComponent(componentId, sampleCount);
        const layout = groupedDotLayout(componentId, sampledMembers, usePca);
        state.dotLayoutCache.set(cacheKey, layout);
        return layout;
    }

    function gridClusterLayout(visibleIds) {
        if (visibleIds.length === 0) {
            clusterCanvas.height = 760;
            return [];
        }

        clusterCanvas.height = Math.max(760, Math.min(1180, Math.round(window.innerHeight * 0.8)));
        const orderedIds = [...visibleIds].sort(compareVisibleClusterIds);
        const items = orderedIds.map(componentId => {
            const node = state.bundle.views.ssn.hierarchy.nodes[componentId];
            return {
                componentId,
                radius: componentRadiusForSize(node.size),
            };
        });
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

        return items.map((item, index) => {
            const column = index % columnCount;
            const row = Math.floor(index / columnCount);
            return {
                componentId: item.componentId,
                radius: item.radius,
                x: outerPadding + (column * cellSize) + (cellSize / 2),
                y: outerPadding + (row * cellSize) + (cellSize / 2),
            };
        });
    }

    // Place circle c externally tangent to already-placed circles a and b (Wang front-chain).
    // Parameter order (b, a, c) mirrors d3's place() so the front-chain insertion stays overlap-free.
    function packPlace(b, a, c) {
        const dx = b.x - a.x;
        const dy = b.y - a.y;
        const d2 = (dx * dx) + (dy * dy);
        if (d2 > 1e-9) {
            const a2 = (a.radius + c.radius) * (a.radius + c.radius);
            const b2 = (b.radius + c.radius) * (b.radius + c.radius);
            if (a2 > b2) {
                const x = (d2 + b2 - a2) / (2 * d2);
                const y = Math.sqrt(Math.max(0, (b2 / d2) - (x * x)));
                c.x = b.x - (x * dx) - (y * dy);
                c.y = b.y - (x * dy) + (y * dx);
            } else {
                const x = (d2 + a2 - b2) / (2 * d2);
                const y = Math.sqrt(Math.max(0, (a2 / d2) - (x * x)));
                c.x = a.x + (x * dx) - (y * dy);
                c.y = a.y + (x * dy) + (y * dx);
            }
        } else {
            c.x = a.x + c.radius;
            c.y = a.y;
        }
    }

    function packIntersects(a, b) {
        const dr = a.radius + b.radius - 1e-6;
        const dx = b.x - a.x;
        const dy = b.y - a.y;
        return dr > 0 && (dr * dr) > ((dx * dx) + (dy * dy));
    }

    function packScore(node) {
        const a = node._;
        const b = node.next._;
        const ab = a.radius + b.radius;
        const dx = ((a.x * b.radius) + (b.x * a.radius)) / ab;
        const dy = ((a.y * b.radius) + (b.y * a.radius)) / ab;
        return (dx * dx) + (dy * dy);
    }

    // Front-chain circle packing (the algorithm behind d3's packSiblings). Packs circles in the
    // order given, compactly from the center outward, with no overlaps. Assigns x,y in place.
    function packSiblingsTight(circles) {
        const n = circles.length;
        if (n === 0) { return; }
        let a = circles[0];
        a.x = 0; a.y = 0;
        if (n === 1) { return; }
        let b = circles[1];
        a.x = -b.radius; b.x = a.radius; b.y = 0;
        if (n === 2) { return; }
        let c = circles[2];
        packPlace(b, a, c);

        let na = {_: a, next: null, previous: null};
        let nb = {_: b, next: null, previous: null};
        let nc = {_: c, next: null, previous: null};
        na.next = nc.previous = nb;
        nb.next = na.previous = nc;
        nc.next = nb.previous = na;

        pack: for (let i = 3; i < n; i++) {
            c = circles[i];
            packPlace(na._, nb._, c);
            nc = {_: c, next: null, previous: null};

            let j = nb.next;
            let k = na.previous;
            let sj = nb._.radius;
            let sk = na._.radius;
            do {
                if (sj <= sk) {
                    if (packIntersects(j._, c)) {
                        nb = j; na.next = nb; nb.previous = na; i--;
                        continue pack;
                    }
                    sj += j._.radius; j = j.next;
                } else {
                    if (packIntersects(k._, c)) {
                        na = k; na.next = nb; nb.previous = na; i--;
                        continue pack;
                    }
                    sk += k._.radius; k = k.previous;
                }
            } while (j !== k.next);

            nc.previous = na; nc.next = nb; na.next = nb.previous = nb = nc;

            let aa = packScore(na);
            let nn = nc;
            while ((nn = nn.next) !== nb) {
                const ca = packScore(nn);
                if (ca < aa) { na = nn; aa = ca; }
            }
            nb = na.next;
        }
    }

    function packedClusterLayout(visibleIds) {
        if (visibleIds.length === 0) {
            clusterCanvas.height = 760;
            return [];
        }
        clusterCanvas.height = Math.max(760, Math.min(1180, Math.round(window.innerHeight * 0.8)));
        const orderedIds = [...visibleIds].sort(compareVisibleClusterIds);
        const circles = orderedIds.map(componentId => ({
            componentId,
            radius: componentRadiusForSize(state.bundle.views.ssn.hierarchy.nodes[componentId].size),
            x: 0,
            y: 0,
        }));
        packSiblingsTight(circles);
        return normalizeComponentLayout(circles, 40).items;
    }

    // Every treemap node is drawn at this fixed world size on one global fixed-pitch lattice, so
    // a node looks identical everywhere and its position never changes with the threshold.
    // TREEMAP_CELL = node square + surrounding padding; cluster boundaries are drawn on the
    // padding between cells so they never overlap a node.
    const TREEMAP_NODE = 13;   // node square side (world units)
    const TREEMAP_CELL = 18;   // lattice pitch = node + padding (world units)
    const TREEMAP_ORIGIN = 58; // world offset of lattice cell (0,0)

    function sign2(value) {
        return (value > 0) - (value < 0);
    }

    // Generalized Hilbert ("gilbert") curve: enumerate every cell of a width x height grid in a
    // Hilbert-like order that works for non-power-of-two, non-square rectangles. Contiguous runs
    // of the returned order are spatially compact and connected, which is what makes a contiguous
    // dendrogram leaf range render as a compact blob with a single staircase boundary. Port of
    // Jakub Cerveny's gilbert2d. Recursion depth is O(log(width*height)).
    function gilbertCurve(width, height) {
        const out = [];
        const generate = (x, y, ax, ay, bx, by) => {
            const w = Math.abs(ax + ay);
            const h = Math.abs(bx + by);
            const dax = sign2(ax), day = sign2(ay);
            const dbx = sign2(bx), dby = sign2(by);
            if (h === 1) {
                for (let i = 0; i < w; i++) { out.push([x, y]); x += dax; y += day; }
                return;
            }
            if (w === 1) {
                for (let i = 0; i < h; i++) { out.push([x, y]); x += dbx; y += dby; }
                return;
            }
            let ax2 = Math.floor(ax / 2), ay2 = Math.floor(ay / 2);
            let bx2 = Math.floor(bx / 2), by2 = Math.floor(by / 2);
            const w2 = Math.abs(ax2 + ay2);
            const h2 = Math.abs(bx2 + by2);
            if (2 * w > 3 * h) {
                if ((w2 % 2) && w > 2) { ax2 += dax; ay2 += day; }
                generate(x, y, ax2, ay2, bx, by);
                generate(x + ax2, y + ay2, ax - ax2, ay - ay2, bx, by);
            } else {
                if ((h2 % 2) && h > 2) { bx2 += dbx; by2 += dby; }
                generate(x, y, bx2, by2, ax2, ay2);
                generate(x + bx2, y + by2, ax, ay, bx - bx2, by - by2);
                generate(x + (ax - dax) + (bx2 - dbx), y + (ay - day) + (by2 - dby),
                    -bx2, -by2, -(ax - ax2), -(ay - ay2));
            }
        };
        if (width >= height) {
            generate(0, 0, width, 0, 0, height);
        } else {
            generate(0, 0, 0, height, width, 0);
        }
        return out;
    }

    // Build (and cache) the global node lattice for the whole network. Every leaf (node) is
    // assigned a fixed cell along a generalized Hilbert (gilbert) curve, whose locality means any
    // contiguous leaf-order range -- i.e. any cluster at any threshold -- maps to a compact,
    // near-square region. The arrangement is threshold-independent: changing the threshold only
    // regroups these fixed cells. Cached by node count.
    function ensureLatticeGlobal() {
        const leafOrder = state.bundle.views.ssn.hierarchy.leaf_order;
        const nodeCount = leafOrder.length;
        if (state.latticeGlobal && state.latticeGlobal.nodeCount === nodeCount) {
            return state.latticeGlobal;
        }
        const aspect = clusterCanvas.width / Math.max(clusterCanvas.height, 1);
        const width = Math.max(1, Math.round(Math.sqrt(nodeCount * aspect)) || 1);
        const height = Math.max(1, Math.ceil(nodeCount / width));
        const curve = gilbertCurve(width, height);
        const colOf = new Int32Array(nodeCount);
        const rowOf = new Int32Array(nodeCount);
        const cellNode = new Int32Array(width * height).fill(-1);
        for (let position = 0; position < nodeCount; position++) {
            const cell = curve[position];
            const nodeIndex = leafOrder[position];
            colOf[nodeIndex] = cell[0];
            rowOf[nodeIndex] = cell[1];
            cellNode[(cell[1] * width) + cell[0]] = nodeIndex;
        }
        state.latticeGlobal = {nodeCount, width, height, colOf, rowOf, cellNode};
        return state.latticeGlobal;
    }

    function latticeCellWorld(column, row) {
        return {
            x: TREEMAP_ORIGIN + ((column + 0.5) * TREEMAP_CELL),
            y: TREEMAP_ORIGIN + ((row + 0.5) * TREEMAP_CELL),
        };
    }

    // Node index under a world-space point, or -1. O(1) via the fixed lattice. The whole cell
    // (node square + its surrounding padding) counts as a hit, not just the drawn square: each
    // cell's padding splits evenly to its neighbors, so a click in the gap between two nodes
    // snaps to the nearer one. This keeps clicks between nodes selecting their cluster instead of
    // falling through. Empty cells (no node) still miss.
    function latticeNodeAtWorld(worldX, worldY) {
        const lattice = state.latticeGlobal;
        if (!lattice) { return -1; }
        const column = Math.floor((worldX - TREEMAP_ORIGIN) / TREEMAP_CELL);
        const row = Math.floor((worldY - TREEMAP_ORIGIN) / TREEMAP_CELL);
        if (column < 0 || row < 0 || column >= lattice.width || row >= lattice.height) { return -1; }
        return lattice.cellNode[(row * lattice.width) + column];
    }

    // Fixed-lattice "treemap". Node positions come from the global node lattice (built once and
    // never moving with the threshold); this only turns the active clusters into lightweight items
    // carrying a bounding box (for culling/fit) and the leaf range (for boundary tracing and
    // hit-testing). Every active cluster is emitted regardless of minimum cluster size -- min size
    // only governs which boundaries are drawn (see renderClusterView), never whether nodes shown.
    function treemapClusterLayout(activeClusterIds) {
        if (activeClusterIds.length === 0) {
            clusterCanvas.height = 760;
            return [];
        }
        clusterCanvas.height = Math.max(760, Math.min(1180, Math.round(window.innerHeight * 0.8)));
        const nodes = state.bundle.views.ssn.hierarchy.nodes;
        const leafOrder = state.bundle.views.ssn.hierarchy.leaf_order;
        const lattice = ensureLatticeGlobal();

        const orderedIds = [...activeClusterIds].sort(
            (a, b) => nodes[a].leaf_start - nodes[b].leaf_start || a - b);
        return orderedIds.map(id => {
            const node = nodes[id];
            const start = node.leaf_start;
            const end = start + node.leaf_count;
            let minCol = Infinity, minRow = Infinity, maxCol = -Infinity, maxRow = -Infinity;
            for (let position = start; position < end; position++) {
                const nodeIndex = leafOrder[position];
                const column = lattice.colOf[nodeIndex];
                const row = lattice.rowOf[nodeIndex];
                if (column < minCol) { minCol = column; }
                if (column > maxCol) { maxCol = column; }
                if (row < minRow) { minRow = row; }
                if (row > maxRow) { maxRow = row; }
            }
            const x0 = TREEMAP_ORIGIN + (minCol * TREEMAP_CELL);
            const y0 = TREEMAP_ORIGIN + (minRow * TREEMAP_CELL);
            const x1 = TREEMAP_ORIGIN + ((maxCol + 1) * TREEMAP_CELL);
            const y1 = TREEMAP_ORIGIN + ((maxRow + 1) * TREEMAP_CELL);
            return {
                componentId: id,
                shape: 'lattice',
                leafStart: start,
                leafCount: node.leaf_count,
                x0, y0, x1, y1,
                x: (x0 + x1) / 2,
                y: (y0 + y1) / 2,
                radius: Math.hypot((x1 - x0) / 2, (y1 - y0) / 2),
            };
        });
    }

    function forceDirectedForestLayout(visibleIds, links, hierarchyNodes, sortBySizeEnabled) {
        if (visibleIds.length === 0) {
            clusterCanvas.height = 760;
            return [];
        }
        clusterCanvas.height = Math.max(760, Math.min(1180, Math.round(window.innerHeight * 0.8)));
        const {adjacency, components} = clusterGraphComponents(visibleIds, links, hierarchyNodes, sortBySizeEnabled);
        // Same tuning the Worker uses (forestForceOptions lives in the shared layout core),
        // so the fallback path cannot drift into drawing a different layout.
        const componentLayouts = components.map(componentIds => simulateComponentLayout(
            componentIds, adjacency, forestForceOptions(componentIds.length), hierarchyNodes));
        return packLayouts(componentLayouts, FOREST_PACK_OPTIONS);
    }

    function tidyForestLayout(visibleIds, links, hierarchyNodes, sortBySizeEnabled) {
        if (visibleIds.length === 0) {
            clusterCanvas.height = 760;
            return [];
        }

        const {adjacency, components} = clusterGraphComponents(visibleIds, links, hierarchyNodes, sortBySizeEnabled);

        clusterCanvas.height = Math.max(760, Math.min(1180, Math.round(window.innerHeight * 0.8)));
        const outerPadding = 58;
        const componentGap = 110;
        const layout = [];
        let currentX = outerPadding;

        components.forEach(componentIds => {
            const tidy = tidyComponentLayout(componentIds, adjacency, hierarchyNodes, currentX);
            tidy.order.forEach(nodeId => layout.push(tidy.positionById.get(nodeId)));
            currentX += tidy.width + outerPadding + componentGap;
        });

        return layout;
    }

    function layoutCacheKey(visibleIds, links, algorithm) {
        const flags = sortComponentsBySizeEnabled() ? 1 : 0;
        let h = flags;
        for (let i = 0; i < visibleIds.length; i++) {
            h = (Math.imul(h, 1664525) + visibleIds[i] + 1013904223) | 0;
        }
        let lh = 0;
        for (let i = 0; i < links.length; i++) {
            const link = links[i];
            lh ^= (Math.imul((link.sourceId * 31 + link.targetId) | 0, 2654435761) | 0) + Math.round((link.weight || 0) * 10000);
        }
        h = (Math.imul(h, 1664525) + lh + 1013904223) | 0;
        return algorithm + ':' + flags + ':' + visibleIds.length + ':' + links.length + ':' + h;
    }

    function computeVisibleLayout(visibleIds, links, algorithm, hierarchyNodes, sortBySizeEnabled) {
        const key = layoutCacheKey(visibleIds, links, algorithm);
        const cached = state.layoutCache.get(key);
        if (cached !== undefined) {
            return {layout: cached.map(item => ({...item})), key, async: false};
        }

        if (algorithm === 'force' && state.layoutWorker) {
            const workerHierarchyNodes = {};
            visibleIds.forEach(id => {
                const node = hierarchyNodes[id];
                workerHierarchyNodes[id] = {size: node.size, leaf_start: node.leaf_start};
            });
            const requestId = ++state.pendingLayoutRequestId;
            state.layoutWorker.postMessage({
                type: 'computeLayout',
                requestId,
                key,
                algorithm,
                visibleIds,
                links: links.map(link => ({sourceId: link.sourceId, targetId: link.targetId, weight: link.weight})),
                hierarchyNodes: workerHierarchyNodes,
                sortBySizeEnabled,
            });
            return {layout: null, key, async: true};
        }

        const layout = (() => {
            if (algorithm === 'grid') {
                return gridClusterLayout(visibleIds);
            }
            if (algorithm === 'packed') {
                return packedClusterLayout(visibleIds);
            }
            if (algorithm === 'treemap') {
                return treemapClusterLayout(visibleIds);
            }
            if (algorithm === 'force') {
                return forceDirectedForestLayout(visibleIds, links, hierarchyNodes, sortBySizeEnabled);
            }
            return tidyForestLayout(visibleIds, links, hierarchyNodes, sortBySizeEnabled);
        })();
        state.layoutCache.set(key, layout.map(item => ({...item})));
        if (state.layoutCache.size > 24) {
            const oldestKey = state.layoutCache.keys().next().value;
            state.layoutCache.delete(oldestKey);
        }
        return {layout, key, async: false};
    }

