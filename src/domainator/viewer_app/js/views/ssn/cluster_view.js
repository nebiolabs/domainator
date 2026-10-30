    function clusterLayoutBounds() {
        if (state.visibleLayout.length === 0) {
            return null;
        }
        let minX = Infinity;
        let minY = Infinity;
        let maxX = -Infinity;
        let maxY = -Infinity;
        state.visibleLayout.forEach(item => {
            if (item.shape === 'rect' || item.shape === 'lattice') {
                minX = Math.min(minX, item.x0);
                minY = Math.min(minY, item.y0);
                maxX = Math.max(maxX, item.x1);
                maxY = Math.max(maxY, item.y1);
            } else {
                minX = Math.min(minX, item.x - item.radius);
                minY = Math.min(minY, item.y - item.radius);
                maxX = Math.max(maxX, item.x + item.radius);
                maxY = Math.max(maxY, item.y + item.radius);
            }
        });
        return {minX, minY, maxX, maxY};
    }

    // Margin left between the fitted content and the canvas edge, shared by
    // "Reset view" and "Focus selection" so both leave the same breathing room.
    const VIEW_FIT_PADDING = 42;

    function fitClusterViewToLayout() {
        const bounds = clusterLayoutBounds();
        if (!bounds) {
            state.viewTransform.scale = 1;
            state.viewTransform.offsetX = 0;
            state.viewTransform.offsetY = 0;
            return;
        }
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
    }

    // World-space bounding box of the selected nodes, grown by each node's own radius so
    // that "fits on screen" means the drawn dot fits, not merely its center. Selected nodes
    // sitting in clusters the current view hides (minimum cluster size, threshold) have no
    // position at all and are skipped; returns null when nothing selected is on screen.
    function selectedNodeBounds() {
        if (!state.bundle || state.selectedNodeIndices.size === 0) {
            return null;
        }
        let minX = Infinity;
        let minY = Infinity;
        let maxX = -Infinity;
        let maxY = -Infinity;
        let count = 0;
        state.visibleLayout.forEach(item => {
            const members = componentMembers(item.componentId);
            // Skip the member layout entirely for clusters holding nothing selected.
            if (!componentSelectionState(members).anySelected) {
                return;
            }
            const component = state.bundle.views.ssn.hierarchy.nodes[item.componentId];
            for (const dot of componentMemberLayout(component, item)) {
                if (!state.selectedNodeIndices.has(dot.memberIndex)) {
                    continue;
                }
                minX = Math.min(minX, dot.x - dot.radius);
                minY = Math.min(minY, dot.y - dot.radius);
                maxX = Math.max(maxX, dot.x + dot.radius);
                maxY = Math.max(maxY, dot.y + dot.radius);
                count += 1;
            }
        });
        return count === 0 ? null : {minX, minY, maxX, maxY, count};
    }

    // Center the view on the selection, zooming out only when it does not already fit.
    // Holding the zoom whenever it fits is deliberate: stepping between selection presets
    // is the common case, and a zoom that changes on every step is disorienting.
    // Centering uses the bounding box center rather than the centroid of the nodes --
    // with a lopsided selection the centroid can leave outliers off screen even at a
    // scale that fits, which would contradict the "all fit" guarantee.
    function focusSelection() {
        const bounds = selectedNodeBounds();
        if (!bounds) {
            return false;
        }
        const padding = VIEW_FIT_PADDING;
        const width = Math.max(1, bounds.maxX - bounds.minX);
        const height = Math.max(1, bounds.maxY - bounds.minY);
        const usableWidth = Math.max(1, clusterCanvas.width - (padding * 2));
        const usableHeight = Math.max(1, clusterCanvas.height - (padding * 2));
        const fitScale = Math.min(usableWidth / width, usableHeight / height);
        if (fitScale < state.viewTransform.scale) {
            // Clamped like every other zoom; a selection wider than minScale allows stays
            // partly off screen, which beats silently exceeding the zoom range.
            state.viewTransform.scale = Math.max(
                state.viewTransform.minScale,
                Math.min(state.viewTransform.maxScale, fitScale),
            );
        }
        const scale = state.viewTransform.scale;
        const centerX = (bounds.minX + bounds.maxX) / 2;
        const centerY = (bounds.minY + bounds.maxY) / 2;
        state.viewTransform.offsetX = (clusterCanvas.width / 2) - (centerX * scale);
        state.viewTransform.offsetY = (clusterCanvas.height / 2) - (centerY * scale);
        return true;
    }

    function drawBadge(ctx, text, x, y) {
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
    }

    function trimmedLinkEndpoints(link) {
        const dx = link.right.x - link.left.x;
        const dy = link.right.y - link.left.y;
        const distance = Math.hypot(dx, dy);
        if (distance < 1e-6) {
            return {startX: link.left.x, startY: link.left.y, endX: link.right.x, endY: link.right.y};
        }
        const unitX = dx / distance;
        const unitY = dy / distance;
        // Trim each end to its own bubble radius so the edge meets each boundary exactly. (Using a
        // shared distance/2 cap would stop the edge at the midpoint — i.e. *inside* a much larger
        // bubble.) If the two bubbles overlap, shrink proportionally so the edge never penetrates.
        let leftOffset = link.left.radius;
        let rightOffset = link.right.radius;
        if (leftOffset + rightOffset > distance) {
            const k = distance / (leftOffset + rightOffset);
            leftOffset *= k;
            rightOffset *= k;
        }
        return {
            startX: link.left.x + (unitX * leftOffset),
            startY: link.left.y + (unitY * leftOffset),
            endX: link.right.x - (unitX * rightOffset),
            endY: link.right.y - (unitY * rightOffset),
        };
    }

    function renderedLinkSegments(link) {
        const trimmed = trimmedLinkEndpoints(link);
        if (currentLayoutAlgorithm() !== 'tree') {
            return [trimmed];
        }
        const middleX = (trimmed.startX + trimmed.endX) / 2;
        return [
            {startX: trimmed.startX, startY: trimmed.startY, endX: middleX, endY: trimmed.startY},
            {startX: middleX, startY: trimmed.startY, endX: middleX, endY: trimmed.endY},
            {startX: middleX, startY: trimmed.endY, endX: trimmed.endX, endY: trimmed.endY},
        ].filter(segment => Math.hypot(segment.endX - segment.startX, segment.endY - segment.startY) > 1e-6);
    }

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
    function linkLabelAnchor(link) {
        const segments = renderedLinkSegments(link);
        if (segments.length === 0) {
            return null;
        }
        const lengths = segments.map(segment =>
            Math.hypot(segment.endX - segment.startX, segment.endY - segment.startY));
        const total = lengths.reduce((sum, value) => sum + value, 0);
        let remaining = total / 2;
        for (let index = 0; index < segments.length; index++) {
            const length = lengths[index];
            if (remaining > length && index < segments.length - 1) {
                remaining -= length;
                continue;
            }
            const segment = segments[index];
            const fraction = length > 1e-9 ? remaining / length : 0;
            const unitX = length > 1e-9 ? (segment.endX - segment.startX) / length : 1;
            const unitY = length > 1e-9 ? (segment.endY - segment.startY) / length : 0;
            return {
                x: segment.startX + ((segment.endX - segment.startX) * fraction),
                y: segment.startY + ((segment.endY - segment.startY) * fraction),
                // Normal pointing "up" for a left-to-right segment, matching the -8px nudge
                // the badge used to get unconditionally.
                normalX: unitY,
                normalY: -unitX,
                length: total,
            };
        }
        return null;
    }

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
    function itemScreenExtent(item, scale = state.viewTransform.scale) {
        if (item.shape === 'lattice' || item.shape === 'rect') {
            return {width: (item.x1 - item.x0) * scale, height: (item.y1 - item.y0) * scale};
        }
        const diameter = item.radius * 2 * scale;
        return {width: diameter, height: diameter};
    }

    function clusterCountLabelFits(text, item, scale = state.viewTransform.scale) {
        const extent = itemScreenExtent(item, scale);
        if (extent.height < CLUSTER_COUNT_LABEL_FONT + 4) {
            return false;
        }
        return estimateTextWidth(text, CLUSTER_COUNT_LABEL_FONT) + 6 <= extent.width;
    }

    function edgeScoreLabelFits(text, linkScreenLength) {
        const badgeWidth = estimateTextWidth(text, EDGE_SCORE_LABEL_FONT) + EDGE_SCORE_BADGE_PADDING;
        return linkScreenLength >= badgeWidth + (EDGE_SCORE_LINK_MARGIN * 2);
    }

    function componentDotGeometry(component, item) {
        const sampleCount = componentDotCount(component.size);
        const dotRadius = componentDotRadius(component.size, item.radius);
        const packingRadius = Math.max(0, item.radius - dotRadius - 0.6);
        return {sampleCount, dotRadius, packingRadius};
    }

    function componentDotLayout(component, item) {
        const {sampleCount, dotRadius, packingRadius} = componentDotGeometry(component, item);
        const minimumDistance = packingRadius > 0 ? Math.min(0.22, (dotRadius * 2.12) / packingRadius) : 0;
        const normalizedLayout = normalizedComponentDotLayout(component.id, sampleCount, minimumDistance);
        return normalizedLayout.map(layout => ({
            memberIndex: layout.memberIndex,
            x: item.x + (layout.x * packingRadius),
            y: item.y + (layout.y * packingRadius),
            radius: dotRadius,
        }));
    }

    // Treemap member layout: every member sits at its fixed cell on the global gilbert lattice,
    // so node positions are identical in every cluster and never move when the threshold changes.
    function componentSquareLayout(component, item) {
        const members = componentMembers(item.componentId);
        const results = [];
        const lattice = state.latticeGlobal;
        if (members.length === 0 || !lattice) { return results; }
        const dotRadius = TREEMAP_NODE / 2;
        for (let i = 0; i < members.length; i++) {
            const nodeIndex = members[i];
            const center = latticeCellWorld(lattice.colOf[nodeIndex], lattice.rowOf[nodeIndex]);
            results.push({
                memberIndex: nodeIndex,
                x: center.x,
                y: center.y,
                radius: dotRadius,
                square: true,
            });
        }
        return results;
    }

    // Shape-aware member placement: lattice tiles read from the global node lattice, everything
    // else keeps the circular phyllotaxis packing. Rendering, SVG export, and hit-testing all
    // route through here so treemap and bubble modes stay in sync.
    function componentMemberLayout(component, item) {
        if (item.shape === 'lattice') {
            return componentSquareLayout(component, item);
        }
        return componentDotLayout(component, item);
    }

    function componentSelectionState(members, nodeSet = state.selectedNodeIndices) {
        let selectedCount = 0;
        members.forEach(nodeIndex => {
            if (nodeSet.has(nodeIndex)) {
                selectedCount += 1;
            }
        });
        return {
            selectedCount,
            anySelected: selectedCount > 0,
            allSelected: members.length > 0 && selectedCount === members.length,
        };
    }

    function renderClusterView(ctx = clusterContext, viewWidth = clusterCanvas.width, viewHeight = clusterCanvas.height, options = {}) {
        // Aliased so the body below can target any context (e.g. a scaled
        // offscreen canvas for high-resolution PNG export) without rewriting
        // every draw call. Defaults render to the on-screen canvas.
        const clusterContext = ctx;
        const clusterCanvas = {width: viewWidth, height: viewHeight};
        const transform = options.transform || state.viewTransform;
        const toScreen = (x, y) => worldToScreenPoint(x, y, transform);
        // Every pan and zoom ends in a repaint. An export repaints the same transform,
        // so refreshing the box then too is harmless.
        if (options.interactive !== false) { updateZoomLevelControl(); }
        clusterContext.clearRect(0, 0, clusterCanvas.width, clusterCanvas.height);
        clusterContext.fillStyle = '#ffffff';
        clusterContext.fillRect(0, 0, clusterCanvas.width, clusterCanvas.height);

        if (!state.bundle) {
            clusterContext.fillStyle = '#5c6a70';
            clusterContext.font = '18px Georgia';
            clusterContext.fillText('Load a bundle to render cluster bubbles.', 30, 50);
            return;
        }

        const selectedNodeOutlines = [];
        // Hovering a preset slot outlines that preset's nodes/clusters in gray.
        // It is a transient hover cue, so exports opt out via options.preview.
        const previewNodeSet = options.preview === false ? null : presetPreviewNodeSet();
        const previewNodeOutlines = [];
        const previewBoundsPaths = [];
        const TAU = Math.PI * 2;
        const drawScale = Math.max(transform.scale, 1e-9);

        // World-space viewport bounds for culling (#4)
        const worldMinX = (0 - transform.offsetX) / drawScale;
        const worldMaxX = (clusterCanvas.width - transform.offsetX) / drawScale;
        const worldMinY = (0 - transform.offsetY) / drawScale;
        const worldMaxY = (clusterCanvas.height - transform.offsetY) / drawScale;

        clusterContext.save();
        clusterContext.translate(transform.offsetX, transform.offsetY);
        clusterContext.scale(drawScale, drawScale);

        state.splitLinks.forEach(link => {
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
            segments.forEach((segment, index) => {
                if (index === 0) {
                    clusterContext.moveTo(segment.startX, segment.startY);
                } else {
                    clusterContext.lineTo(segment.startX, segment.startY);
                }
                clusterContext.lineTo(segment.endX, segment.endY);
            });
            clusterContext.stroke();
        });
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

        state.visibleLayout.forEach(item => {
            // Viewport culling (#4)
            if (item.x + item.radius < worldMinX || item.x - item.radius > worldMaxX ||
                item.y + item.radius < worldMinY || item.y - item.radius > worldMaxY) {
                return;
            }

            const component = state.bundle.views.ssn.hierarchy.nodes[item.componentId];
            const members = componentMembers(item.componentId);
            const selectionState = componentSelectionState(members);
            const previewState = previewNodeSet ? componentSelectionState(members, previewNodeSet) : null;

            if (item.shape === 'lattice') {
                // Lattice clusters are outlined with a staircase along the padding between cells.
                // Nodes are never hidden; a cluster is only outlined when it is big enough (or
                // fully selected). Node squares themselves are drawn in the shared dot pass below.
                const outline = showClusterBounds && component.size >= minOutlineSize;
                if ((outline || selectionState.allSelected) && lattice) {
                    const width = lattice.width;
                    const cellSet = new Set();
                    for (const nodeIndex of members) {
                        cellSet.add((lattice.rowOf[nodeIndex] * width) + lattice.colOf[nodeIndex]);
                    }
                    clusterContext.strokeStyle = selectionState.allSelected ? '#1e2a2f'
                        : (outline ? '#9aa4ad' : 'rgba(30, 42, 47, 0.45)');
                    clusterContext.lineWidth = (selectionState.allSelected ? 2.4 : 1.2) / drawScale;
                    clusterContext.beginPath();
                    for (const nodeIndex of members) {
                        const column = lattice.colOf[nodeIndex];
                        const row = lattice.rowOf[nodeIndex];
                        const ex0 = TREEMAP_ORIGIN + (column * TREEMAP_CELL);
                        const ey0 = TREEMAP_ORIGIN + (row * TREEMAP_CELL);
                        const ex1 = ex0 + TREEMAP_CELL;
                        const ey1 = ey0 + TREEMAP_CELL;
                        if (!cellSet.has((row * width) + (column - 1))) { clusterContext.moveTo(ex0, ey0); clusterContext.lineTo(ex0, ey1); }
                        if (!cellSet.has((row * width) + (column + 1))) { clusterContext.moveTo(ex1, ey0); clusterContext.lineTo(ex1, ey1); }
                        if (!cellSet.has(((row - 1) * width) + column)) { clusterContext.moveTo(ex0, ey0); clusterContext.lineTo(ex1, ey0); }
                        if (!cellSet.has(((row + 1) * width) + column)) { clusterContext.moveTo(ex0, ey1); clusterContext.lineTo(ex1, ey1); }
                    }
                    clusterContext.stroke();
                }
            } else {
                if (previewState && previewState.allSelected) {
                    previewBoundsPaths.push(item);
                }
                const isRect = item.shape === 'rect';
                const traceBounds = () => {
                    clusterContext.beginPath();
                    if (isRect) {
                        clusterContext.rect(item.x0, item.y0, item.x1 - item.x0, item.y1 - item.y0);
                    } else {
                        clusterContext.arc(item.x, item.y, item.radius, 0, TAU);
                    }
                };
                if (showClusterBounds) {
                    clusterContext.fillStyle = 'rgba(245, 246, 248, 0.95)';
                    clusterContext.strokeStyle = selectionState.allSelected ? '#1e2a2f' : '#c4cad2';
                    clusterContext.lineWidth = (selectionState.allSelected ? 3 : 1.5) / drawScale;
                    traceBounds();
                    clusterContext.fill();
                    clusterContext.stroke();
                } else if (selectionState.allSelected) {
                    // Cluster bounds are hidden, but keep a faint outline so the
                    // selection is still discernible.
                    clusterContext.strokeStyle = 'rgba(30, 42, 47, 0.45)';
                    clusterContext.lineWidth = 1.5 / drawScale;
                    traceBounds();
                    clusterContext.stroke();
                }
            }

            // When nodes are hidden we still need to know where the selected
            // ones are so we can hint at them, otherwise skip the dot layout.
            const needsSelectionHint = !selectionState.allSelected && selectionState.anySelected;
            if (!showNodes && !needsSelectionHint) {
                return;
            }

            const dotLayout = componentMemberLayout(component, item);
            for (const dot of dotLayout) {
                if (showNodes) {
                    // Use pre-computed color cache (#2)
                    const color = state.nodeColorCache[dot.memberIndex] ?? nodeColor(dot.memberIndex);
                    const shape = shapeForNode(state.nodeStyles[0], dot.memberIndex, dot.square ? 'cell' : 'circle');
                    const buckets = shape === 'cell' ? squareColorBuckets : dotColorBuckets;
                    const key = shape === 'cell' ? color : color + '|' + shape;
                    let bucket = buckets.get(key);
                    if (bucket === undefined) {
                        bucket = [];
                        buckets.set(key, bucket);
                    }
                    bucket.push(dot);
                }
                if (!selectionState.allSelected && state.selectedNodeIndices.has(dot.memberIndex)) {
                    selectedNodeOutlines.push({x: dot.x, y: dot.y, radius: dot.radius, faint: !showNodes});
                }
                if (previewNodeSet && previewNodeSet.has(dot.memberIndex)) {
                    previewNodeOutlines.push({x: dot.x, y: dot.y, radius: dot.radius});
                }
                if (collectDotLabels && !dotLabelsOverflow &&
                    dot.x >= worldMinX && dot.x <= worldMaxX && dot.y >= worldMinY && dot.y <= worldMaxY) {
                    if (dot.radius * drawScale < MIN_LABELED_DOT_SCREEN_RADIUS) {
                        // Sub-pixel dot: a label beside it would point at nothing visible.
                    } else if (dotLabelCandidates.length >= MAX_LABELED_DOTS) {
                        dotLabelsOverflow = true;
                    } else {
                        dotLabelCandidates.push({x: dot.x, y: dot.y, r: dot.radius, memberIndex: dot.memberIndex});
                    }
                }
            }
        });

        // Batch draw all dots grouped by color (#3)
        dotColorBuckets.forEach((dots, key) => {
            const [color, shape] = key.split('|');
            clusterContext.fillStyle = color;
            clusterContext.beginPath();
            for (const dot of dots) {
                appendNodeGlyph(clusterContext, shape, dot.x, dot.y, dot.radius);
            }
            clusterContext.fill();
        });

        // Batch draw treemap square members grouped by color.
        squareColorBuckets.forEach((squares, color) => {
            clusterContext.fillStyle = color;
            clusterContext.beginPath();
            for (const square of squares) {
                const side = square.radius * 2;
                clusterContext.rect(square.x - square.radius, square.y - square.radius, side, side);
            }
            clusterContext.fill();
        });

        clusterContext.restore();

        // Drawn before the selection outlines and at a wider radius, so a node that
        // is both previewed and selected shows a gray ring outside a black one.
        if (previewNodeOutlines.length > 0 || previewBoundsPaths.length > 0) {
            clusterContext.save();
            clusterContext.strokeStyle = '#9aa4ad';
            clusterContext.lineWidth = 2;
            previewBoundsPaths.forEach(item => {
                clusterContext.beginPath();
                if (item.shape === 'rect') {
                    const topLeft = toScreen(item.x0, item.y0);
                    const bottomRight = toScreen(item.x1, item.y1);
                    clusterContext.rect(topLeft.x, topLeft.y, bottomRight.x - topLeft.x, bottomRight.y - topLeft.y);
                } else {
                    const center = toScreen(item.x, item.y);
                    clusterContext.arc(center.x, center.y, (item.radius * drawScale) + 3, 0, TAU);
                }
                clusterContext.stroke();
            });
            previewNodeOutlines.forEach(outline => {
                const screenPoint = toScreen(outline.x, outline.y);
                const screenRadius = Math.max(4, (outline.radius * drawScale) + 3.4);
                clusterContext.beginPath();
                clusterContext.arc(screenPoint.x, screenPoint.y, screenRadius, 0, TAU);
                clusterContext.stroke();
            });
            clusterContext.restore();
        }

        if (selectedNodeOutlines.length > 0) {
            clusterContext.save();
            selectedNodeOutlines.forEach(outline => {
                const screenPoint = toScreen(outline.x, outline.y);
                const screenRadius = Math.max(3, (outline.radius * drawScale) + 1.6);
                clusterContext.beginPath();
                clusterContext.arc(screenPoint.x, screenPoint.y, screenRadius, 0, TAU);
                if (outline.faint) {
                    // Nodes are hidden, so leave a faint mark where the selected
                    // node sits instead of outlining a dot that isn't drawn.
                    clusterContext.fillStyle = 'rgba(30, 42, 47, 0.30)';
                    clusterContext.fill();
                } else {
                    clusterContext.strokeStyle = '#1e2a2f';
                    clusterContext.lineWidth = 1.8;
                    clusterContext.stroke();
                }
            });
            clusterContext.restore();
        }

        if (showEdgeScoresEnabled()) {
            state.splitLinks.forEach(link => {
                const anchor = linkLabelAnchor(link);
                if (anchor === null) {
                    return;
                }
                const text = formatValue(link.threshold);
                if (!edgeScoreLabelFits(text, anchor.length * drawScale)) {
                    return;
                }
                const point = toScreen(anchor.x, anchor.y);
                drawBadge(clusterContext, text,
                    point.x + (anchor.normalX * EDGE_SCORE_LABEL_OFFSET),
                    point.y + (anchor.normalY * EDGE_SCORE_LABEL_OFFSET));
            });
        }

        // Node-count labels: one per cluster bubble, centered inside, wherever the
        // bubble is big enough on screen to hold the text.
        if (showNodeCountsEnabled()) {
            clusterContext.fillStyle = '#5c6a70';
            clusterContext.font = '600 ' + CLUSTER_COUNT_LABEL_FONT + 'px Georgia';
            clusterContext.textAlign = 'center';
            clusterContext.textBaseline = 'middle';
            state.visibleLayout.forEach(item => {
                const screenPoint = toScreen(item.x, item.y);
                const screenRadius = item.radius * drawScale;
                if (screenPoint.x + screenRadius < 0 || screenPoint.x - screenRadius > clusterCanvas.width ||
                    screenPoint.y + screenRadius < 0 || screenPoint.y - screenRadius > clusterCanvas.height) {
                    return;
                }
                const component = state.bundle.views.ssn.hierarchy.nodes[item.componentId];
                const text = component.size.toLocaleString();
                if (!clusterCountLabelFits(text, item, drawScale)) {
                    return;
                }
                clusterContext.fillText(text, screenPoint.x, screenPoint.y + 4);
            });
        }

        // "Label by" metadata labels: one per member node, drawn next to its dot. Shown only when
        // few enough nodes are in view (collected above with a cap), so they appear as you zoom in.
        if (collectDotLabels && !dotLabelsOverflow && dotLabelCandidates.length > 0) {
            clusterContext.fillStyle = '#1e2a2f';
            clusterContext.font = '12px Georgia';
            clusterContext.textAlign = 'center';
            clusterContext.textBaseline = 'middle';
            dotLabelCandidates.forEach(candidate => {
                let text;
                if (dotLabelField === '__node_id__') {
                    text = nodeId(candidate.memberIndex);
                } else {
                    const value = metadataValue(candidate.memberIndex, dotLabelField);
                    // Null/empty values are simply not labeled (no node_id fallback).
                    if (value === null || value === undefined || value === '') {
                        return;
                    }
                    text = formatValue(value);
                }
                if (!text) {
                    return;
                }
                const screenPoint = toScreen(candidate.x, candidate.y);
                clusterContext.fillText(String(text).slice(0, 24), screenPoint.x, screenPoint.y);
            });
        }

        if (state.selectionBox && options.interactive !== false) {
            const box = state.selectionBox;
            clusterContext.save();
            clusterContext.fillStyle = 'rgba(200, 85, 61, 0.10)';
            clusterContext.strokeStyle = 'rgba(200, 85, 61, 0.62)';
            clusterContext.setLineDash([8, 6]);
            clusterContext.lineWidth = 1.5;
            clusterContext.fillRect(box.left, box.top, box.width, box.height);
            clusterContext.strokeRect(box.left, box.top, box.width, box.height);
            clusterContext.restore();
        }

        if (state.layoutComputing && options.interactive !== false) {
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
        }
    }

    function scheduleClusterRender() {
        if (state.pendingClusterRenderFrame !== null) {
            return;
        }
        state.pendingClusterRenderFrame = window.requestAnimationFrame(() => {
            state.pendingClusterRenderFrame = null;
            renderClusterView();
        });
    }

    function applyComputedLayout(layout, visibleGraph, layoutAlgorithm, resetView) {
        const layoutById = new Map(layout.map(item => [item.componentId, item]));
        const links = (layoutAlgorithm === 'grid' || layoutAlgorithm === 'packed' || layoutAlgorithm === 'treemap')
            ? []
            : visibleGraph.links
                .map(link => {
                    const left = layoutById.get(link.sourceId);
                    const right = layoutById.get(link.targetId);
                    if (!left || !right) { return null; }
                    return {left, right, threshold: link.weight, collapsed: link.collapsed || 0};
                })
                .filter(Boolean);
        state.visibleLayout = layout;
        state.splitLinks = links;
        // For the lattice layout, map every node to the visible cluster item that owns it so a
        // click can resolve cell -> node -> cluster in O(1) (bounding boxes overlap, so we can't
        // hit-test clusters by bbox). Rebuilt whenever the visible cluster set changes.
        if (layout.length > 0 && layout[0].shape === 'lattice') {
            const leafOrder = state.bundle.views.ssn.hierarchy.leaf_order;
            const nodeItem = new Int32Array(leafOrder.length).fill(-1);
            layout.forEach((item, itemIndex) => {
                const end = item.leafStart + item.leafCount;
                for (let position = item.leafStart; position < end; position++) {
                    nodeItem[leafOrder[position]] = itemIndex;
                }
            });
            state.latticeNodeItem = nodeItem;
        } else {
            state.latticeNodeItem = null;
        }
        document.getElementById('stat-clusters').textContent = layout.length.toLocaleString();
        document.getElementById('stat-links').textContent = links.length.toLocaleString();
        if (state.pendingRestoreViewTransform) {
            // A session restored a pan/zoom; adopt it instead of auto-fitting,
            // which would otherwise immediately discard it.
            const restored = state.pendingRestoreViewTransform;
            state.pendingRestoreViewTransform = null;
            state.viewTransform.scale = restored.scale;
            state.viewTransform.offsetX = restored.offsetX;
            state.viewTransform.offsetY = restored.offsetY;
        } else if (resetView) {
            fitClusterViewToLayout();
        }
        renderClusterView();
    }

    function drawClusterView(resetView = true) {
        if (!state.bundle) {
            renderClusterView();
            return;
        }

        const minClusterSize = Math.max(1, Number(document.getElementById('min-cluster-size').value) || 1);
        const layoutAlgorithm = currentLayoutAlgorithm();
        const thresholdValue = selectedThresholdValue();
        const activeClusterIds = activeClustersAtThreshold(thresholdValue);
        const hierarchyNodes = state.bundle.views.ssn.hierarchy.nodes;
        const sortBySizeEnabled = sortComponentsBySizeEnabled();

        // The lattice ("treemap") layout never hides nodes -- every active cluster is laid out and
        // the minimum cluster size only suppresses small clusters' outlines. Every other layout
        // filters nodes by minimum cluster size as usual.
        const latticeMode = layoutAlgorithm === 'treemap';
        const visibleGraph = latticeMode
            ? {visibleIds: activeClusterIds, links: [], hiddenNodes: 0, collapsedPaths: 0, collapsedClusters: 0}
            : mstLinksForActiveClusters(activeClusterIds, minClusterSize,
                leafPruningOnlyEnabled(), collapseLongPathsEnabled());

        state.activeClusters = activeClusterIds;
        state.visibleClusters = visibleGraph.visibleIds;

        const hidden = visibleGraph.hiddenNodes;
        const shown = state.bundle.nodes.length - hidden;
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
        if (result.async) {
            state._pendingVisibleGraph = visibleGraph;
            state._pendingLayoutAlgorithm = layoutAlgorithm;
            state._pendingLayoutResetView = resetView;
            state.layoutComputing = true;
            renderClusterView();
            return;
        }
        state.layoutComputing = false;
        applyComputedLayout(result.layout, visibleGraph, layoutAlgorithm, resetView);
    }

    function setupLayoutWorker() {
        try {
            const blob = new Blob([LAYOUT_WORKER_CODE], {type: 'application/javascript'});
            const workerUrl = URL.createObjectURL(blob);
            const worker = new Worker(workerUrl);
            worker.onmessage = function(event) {
                const {requestId, key, layout} = event.data;
                if (requestId !== state.pendingLayoutRequestId) { return; }
                state.layoutCache.set(key, layout.map(item => ({...item})));
                if (state.layoutCache.size > 24) {
                    state.layoutCache.delete(state.layoutCache.keys().next().value);
                }
                state.layoutComputing = false;
                applyComputedLayout(layout, state._pendingVisibleGraph, state._pendingLayoutAlgorithm, state._pendingLayoutResetView);
            };
            worker.onerror = function(error) {
                console.warn('Layout worker error, falling back to main thread:', error);
                state.layoutWorker = null;
                state.layoutComputing = false;
                renderClusterView();
            };
            state.layoutWorker = worker;
        } catch (error) {
            console.warn('Layout worker unavailable, using main thread:', error);
        }
    }

    // Cluster numbers for the current threshold, keyed by node index: connected
    // components ranked by size with the largest numbered 1. This mirrors
    // domainator build_ssn's --cluster (CLUSTER_COLUMN, rename_labels_by_frequency + 1),
    // so a column made here lines up with what that tool writes. Equal-sized
    // clusters are ordered by hierarchy position rather than by lowest member
    // index, so their numbers can differ from domainator build_ssn's on ties.
    function clusterNumbersAtCurrentThreshold() {
        const activeClusterIds = activeClustersAtThreshold(selectedThresholdValue());
        const componentAssignments = activeClusterAssignments(activeClusterIds);
        const rankedClusterIds = [...activeClusterIds].sort((leftId, rightId) => {
            const leftNode = state.bundle.views.ssn.hierarchy.nodes[leftId];
            const rightNode = state.bundle.views.ssn.hierarchy.nodes[rightId];
            return rightNode.size - leftNode.size || leftNode.leaf_start - rightNode.leaf_start || leftId - rightId;
        });
        const clusterByComponentId = new Map(rankedClusterIds.map((componentId, clusterIndex) => [componentId, clusterIndex + 1]));
        return {
            clusterCount: activeClusterIds.length,
            numberFor: nodeIndex => clusterByComponentId.get(componentAssignments[nodeIndex]) ?? null,
        };
    }

    function clusterExportSize() {
        const size = {width: clusterCanvas.width, height: clusterCanvas.height, transform: {...state.viewTransform}};
        if (entireContentExport()) {
            if (state.layoutComputing) { throw new Error('Wait for the network layout to finish before exporting entire content.'); }
            const bounds = clusterLayoutBounds();
            if (bounds) {
                const scale = state.viewTransform.scale;
                const measurer = document.createElement('canvas').getContext('2d');
                if (!measurer) { throw new Error('Could not measure network labels.'); }
                let padding = 12;
                const field = currentLabelField();
                measurer.font = '12px Georgia';
                for (const item of state.visibleLayout) {
                    const component = state.bundle.views.ssn.hierarchy.nodes[item.componentId];
                    for (const dot of componentMemberLayout(component, item)) {
                        if (!renderNodesEnabled() && !state.selectedNodeIndices.has(dot.memberIndex)) { continue; }
                        bounds.minX = Math.min(bounds.minX, dot.x - dot.radius);
                        bounds.minY = Math.min(bounds.minY, dot.y - dot.radius);
                        bounds.maxX = Math.max(bounds.maxX, dot.x + dot.radius);
                        bounds.maxY = Math.max(bounds.maxY, dot.y + dot.radius);
                        if (field && renderNodesEnabled()) {
                            const value = field === '__node_id__' ? nodeId(dot.memberIndex) : metadataValue(dot.memberIndex, field);
                            padding = Math.max(padding, measurer.measureText(String(value ?? '').slice(0, 24)).width / 2 + 4);
                        }
                    }
                }
                measurer.font = EDGE_SCORE_LABEL_FONT + 'px Georgia';
                for (const link of state.splitLinks) {
                    for (const segment of renderedLinkSegments(link)) {
                        bounds.minX = Math.min(bounds.minX, segment.startX, segment.endX);
                        bounds.maxX = Math.max(bounds.maxX, segment.startX, segment.endX);
                        bounds.minY = Math.min(bounds.minY, segment.startY, segment.endY);
                        bounds.maxY = Math.max(bounds.maxY, segment.startY, segment.endY);
                    }
                    if (showEdgeScoresEnabled()) {
                        padding = Math.max(padding, measurer.measureText(formatValue(link.threshold)).width / 2 +
                            EDGE_SCORE_BADGE_PADDING / 2 + EDGE_SCORE_LABEL_OFFSET + 2, 20);
                    }
                }
                size.width = Math.ceil((bounds.maxX - bounds.minX) * scale + 2 * padding);
                size.height = Math.ceil((bounds.maxY - bounds.minY) * scale + 2 * padding);
                size.transform = {scale, offsetX: padding - bounds.minX * scale, offsetY: padding - bounds.minY * scale};
            }
        }
        validateImageExportSize(size.width, size.height);
        return size;
    }

    function exportClusterPNG() {
        if (!state.bundle) { return; }
        const size = clusterExportSize();
        return exportCanvasPNG(size, ctx => renderClusterView(ctx, size.width, size.height,
            {transform: size.transform, preview: false, interactive: false}), imageExportSuffix('png', true));
    }

    // Render the current cluster view to a standalone SVG document. Geometry is
    // baked into screen space (every point passed through worldToScreenPoint)
    // so the markup needs no transforms and matches the canvas pixel-for-pixel.
    // This mirrors renderClusterView(); keep the two in sync.
    function buildClusterViewSVG(size = clusterExportSize()) {
        const {width, height, transform} = size;
        const toScreen = (x, y) => worldToScreenPoint(x, y, transform);
        const fmt = value => Math.round(value * 100) / 100;
        const parts = [];
        parts.push('<svg xmlns="http://www.w3.org/2000/svg" width="' + width + '" height="' + height +
            '" viewBox="0 0 ' + width + ' ' + height + '" font-family="Georgia, serif">');
        parts.push('<rect x="0" y="0" width="' + width + '" height="' + height + '" fill="#ffffff"/>');

        if (!state.bundle) {
            parts.push('<text x="30" y="50" font-size="18" fill="#5c6a70">Load a bundle to render cluster bubbles.</text>');
            parts.push('</svg>');
            return parts.join('\n');
        }

        const TAU = Math.PI * 2;
        const drawScale = Math.max(transform.scale, 1e-9);
        const worldMinX = (0 - transform.offsetX) / drawScale;
        const worldMaxX = (width - transform.offsetX) / drawScale;
        const worldMinY = (0 - transform.offsetY) / drawScale;
        const worldMaxY = (height - transform.offsetY) / drawScale;

        // Split links (edges). Collapsed paths go in their own dashed group, matching
        // renderClusterView.
        const linkColor = svgColorParts('rgba(92, 106, 112, 0.42)');
        const linkPaths = [];
        const collapsedLinkPaths = [];
        state.splitLinks.forEach(link => {
            const segments = renderedLinkSegments(link);
            if (segments.length === 0) {
                return;
            }
            let d = '';
            segments.forEach((segment, index) => {
                const start = toScreen(segment.startX, segment.startY);
                const end = toScreen(segment.endX, segment.endY);
                if (index === 0) {
                    d += 'M' + fmt(start.x) + ' ' + fmt(start.y);
                } else {
                    d += ' L' + fmt(start.x) + ' ' + fmt(start.y);
                }
                d += ' L' + fmt(end.x) + ' ' + fmt(end.y);
            });
            (link.collapsed ? collapsedLinkPaths : linkPaths).push('<path d="' + d + '"/>');
        });
        if (linkPaths.length > 0) {
            parts.push('<g fill="none" stroke="' + linkColor.color + '" stroke-opacity="' + linkColor.opacity +
                '" stroke-width="1.6">');
            parts.push(linkPaths.join(''));
            parts.push('</g>');
        }
        if (collapsedLinkPaths.length > 0) {
            parts.push('<g fill="none" stroke="' + linkColor.color + '" stroke-opacity="' + linkColor.opacity +
                '" stroke-width="1.6" stroke-dasharray="' + LINK_DASH_ON + ' ' + LINK_DASH_OFF + '">');
            parts.push(collapsedLinkPaths.join(''));
            parts.push('</g>');
        }

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

        state.visibleLayout.forEach(item => {
            if (item.x + item.radius < worldMinX || item.x - item.radius > worldMaxX ||
                item.y + item.radius < worldMinY || item.y - item.radius > worldMaxY) {
                return;
            }

            const component = state.bundle.views.ssn.hierarchy.nodes[item.componentId];
            const members = componentMembers(item.componentId);
            const selectionState = componentSelectionState(members);
            const center = toScreen(item.x, item.y);
            const screenRadius = item.radius * drawScale;
            const isLattice = item.shape === 'lattice';
            const isRect = item.shape === 'rect';
            const rectTopLeft = isRect ? toScreen(item.x0, item.y0) : null;
            const rectWidth = isRect ? (item.x1 - item.x0) * drawScale : 0;
            const rectHeight = isRect ? (item.y1 - item.y0) * drawScale : 0;
            const boundsShape = (fillAttrs, strokeAttrs) => (isRect
                ? '<rect x="' + fmt(rectTopLeft.x) + '" y="' + fmt(rectTopLeft.y) + '" width="' + fmt(rectWidth) +
                    '" height="' + fmt(rectHeight) + '" ' + fillAttrs + ' ' + strokeAttrs + '/>'
                : '<circle cx="' + fmt(center.x) + '" cy="' + fmt(center.y) + '" r="' + fmt(screenRadius) +
                    '" ' + fillAttrs + ' ' + strokeAttrs + '/>');

            if (isLattice) {
                // Staircase outline traced along the padding between cells (screen space).
                const outline = showClusterBounds && component.size >= minOutlineSize;
                if ((outline || selectionState.allSelected) && lattice) {
                    const width = lattice.width;
                    const cellSet = new Set();
                    for (const nodeIndex of members) {
                        cellSet.add((lattice.rowOf[nodeIndex] * width) + lattice.colOf[nodeIndex]);
                    }
                    let d = '';
                    const edge = (wx0, wy0, wx1, wy1) => {
                        const a = toScreen(wx0, wy0);
                        const b = toScreen(wx1, wy1);
                        d += 'M' + fmt(a.x) + ' ' + fmt(a.y) + 'L' + fmt(b.x) + ' ' + fmt(b.y);
                    };
                    for (const nodeIndex of members) {
                        const column = lattice.colOf[nodeIndex];
                        const row = lattice.rowOf[nodeIndex];
                        const ex0 = TREEMAP_ORIGIN + (column * TREEMAP_CELL);
                        const ey0 = TREEMAP_ORIGIN + (row * TREEMAP_CELL);
                        const ex1 = ex0 + TREEMAP_CELL;
                        const ey1 = ey0 + TREEMAP_CELL;
                        if (!cellSet.has((row * width) + (column - 1))) { edge(ex0, ey0, ex0, ey1); }
                        if (!cellSet.has((row * width) + (column + 1))) { edge(ex1, ey0, ex1, ey1); }
                        if (!cellSet.has(((row - 1) * width) + column)) { edge(ex0, ey0, ex1, ey0); }
                        if (!cellSet.has(((row + 1) * width) + column)) { edge(ex0, ey1, ex1, ey1); }
                    }
                    if (d !== '') {
                        const stroke = svgColorParts(selectionState.allSelected ? '#1e2a2f' : (outline ? '#9aa4ad' : 'rgba(30, 42, 47, 0.45)'));
                        boundsParts.push('<path d="' + d + '" fill="none" stroke="' + stroke.color +
                            (stroke.opacity !== 1 ? '" stroke-opacity="' + stroke.opacity : '') +
                            '" stroke-width="' + (selectionState.allSelected ? 2.4 : 1.2) + '"/>');
                    }
                }
            } else if (showClusterBounds) {
                const stroke = svgColorParts(selectionState.allSelected ? '#1e2a2f' : '#c4cad2');
                boundsParts.push(boundsShape(
                    'fill="' + boundsFill.color + '" fill-opacity="' + boundsFill.opacity + '"',
                    'stroke="' + stroke.color + '" stroke-width="' + (selectionState.allSelected ? 3 : 1.5) + '"'));
            } else if (selectionState.allSelected) {
                const stroke = svgColorParts('rgba(30, 42, 47, 0.45)');
                boundsParts.push(boundsShape(
                    'fill="none"',
                    'stroke="' + stroke.color + '" stroke-opacity="' + stroke.opacity + '" stroke-width="1.5"'));
            }

            const needsSelectionHint = !selectionState.allSelected && selectionState.anySelected;
            if (!showNodes && !needsSelectionHint) {
                return;
            }

            const dotLayout = componentMemberLayout(component, item);
            for (const dot of dotLayout) {
                if (showNodes) {
                    const color = state.nodeColorCache[dot.memberIndex] ?? nodeColor(dot.memberIndex);
                    const shape = shapeForNode(state.nodeStyles[0], dot.memberIndex, dot.square ? 'cell' : 'circle');
                    const buckets = shape === 'cell' ? squareColorBuckets : dotColorBuckets;
                    const key = shape === 'cell' ? color : color + '|' + shape;
                    let bucket = buckets.get(key);
                    if (bucket === undefined) {
                        bucket = [];
                        buckets.set(key, bucket);
                    }
                    const dotCenter = toScreen(dot.x, dot.y);
                    bucket.push({x: dotCenter.x, y: dotCenter.y, r: dot.radius * drawScale});
                }
                if (!selectionState.allSelected && state.selectedNodeIndices.has(dot.memberIndex)) {
                    selectedNodeOutlines.push({x: dot.x, y: dot.y, radius: dot.radius, faint: !showNodes});
                }
                if (collectDotLabels && !dotLabelsOverflow &&
                    dot.x >= worldMinX && dot.x <= worldMaxX && dot.y >= worldMinY && dot.y <= worldMaxY) {
                    if (dot.radius * drawScale < MIN_LABELED_DOT_SCREEN_RADIUS) {
                        // Sub-pixel dot: a label beside it would point at nothing visible.
                    } else if (dotLabelCandidates.length >= MAX_LABELED_DOTS) {
                        dotLabelsOverflow = true;
                    } else {
                        dotLabelCandidates.push({x: dot.x, y: dot.y, memberIndex: dot.memberIndex});
                    }
                }
            }
        });

        if (boundsParts.length > 0) {
            parts.push(boundsParts.join(''));
        }

        dotColorBuckets.forEach((dots, key) => {
            const [color, shape] = key.split('|');
            const fill = svgColorParts(color);
            const circles = dots.map(dot => nodeGlyphSVG(shape, dot.x, dot.y, dot.r, fmt)).join('');
            parts.push('<g fill="' + fill.color + '"' + (fill.opacity !== 1 ? ' fill-opacity="' + fill.opacity + '"' : '') + '>' + circles + '</g>');
        });

        squareColorBuckets.forEach((squares, color) => {
            const fill = svgColorParts(color);
            const rects = squares.map(square => '<rect x="' + fmt(square.x - square.r) + '" y="' + fmt(square.y - square.r) +
                '" width="' + fmt(square.r * 2) + '" height="' + fmt(square.r * 2) + '"/>').join('');
            parts.push('<g fill="' + fill.color + '"' + (fill.opacity !== 1 ? ' fill-opacity="' + fill.opacity + '"' : '') + '>' + rects + '</g>');
        });

        if (selectedNodeOutlines.length > 0) {
            const outlineParts = selectedNodeOutlines.map(outline => {
                const screenPoint = toScreen(outline.x, outline.y);
                const screenRadius = Math.max(3, (outline.radius * drawScale) + 1.6);
                if (outline.faint) {
                    const fill = svgColorParts('rgba(30, 42, 47, 0.30)');
                    return '<circle cx="' + fmt(screenPoint.x) + '" cy="' + fmt(screenPoint.y) + '" r="' + fmt(screenRadius) +
                        '" fill="' + fill.color + '" fill-opacity="' + fill.opacity + '"/>';
                }
                return '<circle cx="' + fmt(screenPoint.x) + '" cy="' + fmt(screenPoint.y) + '" r="' + fmt(screenRadius) +
                    '" fill="none" stroke="#1e2a2f" stroke-width="1.8"/>';
            });
            parts.push(outlineParts.join(''));
        }

        // Edge score badges. The fit test is edgeScoreLabelFits, shared with
        // renderClusterView so the export shows the labels the screen showed.
        if (showEdgeScoresEnabled()) {
            clusterContext.save();
            clusterContext.font = EDGE_SCORE_LABEL_FONT + 'px Georgia';
            state.splitLinks.forEach(link => {
                const anchor = linkLabelAnchor(link);
                if (anchor === null) {
                    return;
                }
                const text = formatValue(link.threshold);
                if (!edgeScoreLabelFits(text, anchor.length * drawScale)) {
                    return;
                }
                const point = toScreen(anchor.x, anchor.y);
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
            });
            clusterContext.restore();
        }

        // Node-count labels (one per bubble), under the same fit test as the canvas.
        if (showNodeCountsEnabled()) {
            state.visibleLayout.forEach(item => {
                const screenPoint = toScreen(item.x, item.y);
                const screenRadius = item.radius * drawScale;
                if (screenPoint.x + screenRadius < 0 || screenPoint.x - screenRadius > width ||
                    screenPoint.y + screenRadius < 0 || screenPoint.y - screenRadius > height) {
                    return;
                }
                const component = state.bundle.views.ssn.hierarchy.nodes[item.componentId];
                const text = component.size.toLocaleString();
                if (!clusterCountLabelFits(text, item, drawScale)) {
                    return;
                }
                parts.push('<text x="' + fmt(screenPoint.x) + '" y="' + fmt(screenPoint.y + 4) + '" font-size="' +
                    CLUSTER_COUNT_LABEL_FONT + '" font-weight="600"' +
                    ' fill="#5c6a70" text-anchor="middle" dominant-baseline="central">' + escapeXml(text) + '</text>');
            });
        }

        // Per-node "Label by" metadata labels.
        if (collectDotLabels && !dotLabelsOverflow && dotLabelCandidates.length > 0) {
            dotLabelCandidates.forEach(candidate => {
                let text;
                if (dotLabelField === '__node_id__') {
                    text = nodeId(candidate.memberIndex);
                } else {
                    const value = metadataValue(candidate.memberIndex, dotLabelField);
                    if (value === null || value === undefined || value === '') {
                        return;
                    }
                    text = formatValue(value);
                }
                if (!text) {
                    return;
                }
                const screenPoint = toScreen(candidate.x, candidate.y);
                parts.push('<text x="' + fmt(screenPoint.x) + '" y="' + fmt(screenPoint.y) + '" font-size="12" fill="#1e2a2f"' +
                    ' text-anchor="middle" dominant-baseline="central">' + escapeXml(String(text).slice(0, 24)) + '</text>');
            });
        }

        parts.push('</svg>');
        return parts.join('\n');
    }

    function exportClusterSVG() {
        if (!state.bundle) {
            return;
        }
        const svg = buildClusterViewSVG();
        const blob = new Blob([svg], {type: 'image/svg+xml'});
        const href = URL.createObjectURL(blob);
        triggerDownload(href, exportBaseName() + imageExportSuffix('svg'), true);
    }

    function toggleSelectionForComponent(componentId) {
        const members = componentMembers(componentId);
        const allSelected = members.every(nodeIndex => state.selectedNodeIndices.has(nodeIndex));
        members.forEach(nodeIndex => {
            if (allSelected) {
                state.selectedNodeIndices.delete(nodeIndex);
            } else {
                state.selectedNodeIndices.add(nodeIndex);
            }
        });
        renderClusterView();
        updateMetadataTable();
    }

    function toggleSelectionForNode(nodeIndex) {
        if (state.selectedNodeIndices.has(nodeIndex)) {
            state.selectedNodeIndices.delete(nodeIndex);
        } else {
            state.selectedNodeIndices.add(nodeIndex);
        }
        renderClusterView();
        updateMetadataTable();
    }

    function hitTestNodeAt(screenX, screenY) {
        const worldPoint = screenToWorldPoint(screenX, screenY);
        if (state.latticeNodeItem) {
            const nodeIndex = latticeNodeAtWorld(worldPoint.x, worldPoint.y);
            return nodeIndex >= 0 ? nodeIndex : null;
        }
        for (let index = state.visibleLayout.length - 1; index >= 0; index--) {
            const item = state.visibleLayout[index];
            if (item.shape === 'rect') {
                if (worldPoint.x < item.x0 || worldPoint.x > item.x1 ||
                    worldPoint.y < item.y0 || worldPoint.y > item.y1) {
                    continue;
                }
            } else {
                const dx = worldPoint.x - item.x;
                const dy = worldPoint.y - item.y;
                if ((dx * dx) + (dy * dy) > item.radius * item.radius) {
                    continue;
                }
            }
            const component = state.bundle.views.ssn.hierarchy.nodes[item.componentId];
            const dotLayout = componentMemberLayout(component, item);
            const dotRadius = dotLayout.length > 0 ? dotLayout[0].radius : componentDotGeometry(component, item).dotRadius;
            const hitRadius = Math.max(dotRadius, 5 / Math.max(state.viewTransform.scale, 0.1));
            let nearestNodeIndex = null;
            let nearestDistanceSq = Infinity;
            for (const dot of dotLayout) {
                const nodeDx = worldPoint.x - dot.x;
                const nodeDy = worldPoint.y - dot.y;
                const distanceSq = (nodeDx * nodeDx) + (nodeDy * nodeDy);
                if (distanceSq <= hitRadius * hitRadius && distanceSq < nearestDistanceSq) {
                    nearestDistanceSq = distanceSq;
                    nearestNodeIndex = dot.memberIndex;
                }
            }
            if (nearestNodeIndex !== null) {
                return nearestNodeIndex;
            }
        }
        return null;
    }

    function hitTestComponentAt(screenX, screenY) {
        const worldPoint = screenToWorldPoint(screenX, screenY);
        if (state.latticeNodeItem) {
            const nodeIndex = latticeNodeAtWorld(worldPoint.x, worldPoint.y);
            if (nodeIndex < 0) { return null; }
            const itemIndex = state.latticeNodeItem[nodeIndex];
            return itemIndex >= 0 ? state.visibleLayout[itemIndex] : null;
        }
        for (let index = state.visibleLayout.length - 1; index >= 0; index--) {
            const item = state.visibleLayout[index];
            if (item.shape === 'rect') {
                if (worldPoint.x >= item.x0 && worldPoint.x <= item.x1 &&
                    worldPoint.y >= item.y0 && worldPoint.y <= item.y1) {
                    return item;
                }
            } else {
                const dx = worldPoint.x - item.x;
                const dy = worldPoint.y - item.y;
                if ((dx * dx) + (dy * dy) <= item.radius * item.radius) {
                    return item;
                }
            }
        }
        return null;
    }

    function itemIntersectsSelectionBox(item, rect) {
        if (item.shape === 'rect' || item.shape === 'lattice') {
            const p0 = worldToScreenPoint(item.x0, item.y0);
            const p1 = worldToScreenPoint(item.x1, item.y1);
            return rectIntersectsRect(p0.x, p0.y, p1.x, p1.y, rect);
        }
        const screenPoint = worldToScreenPoint(item.x, item.y);
        const screenRadius = item.radius * state.viewTransform.scale;
        return circleIntersectsRect(screenPoint.x, screenPoint.y, screenRadius, rect);
    }

    function selectComponentsInBox(rect, deselect) {
        let changed = false;
        state.visibleLayout.forEach(item => {
            if (!itemIntersectsSelectionBox(item, rect)) {
                return;
            }
            componentMembers(item.componentId).forEach(nodeIndex => {
                if (deselect) {
                    if (state.selectedNodeIndices.delete(nodeIndex)) {
                        changed = true;
                    }
                } else if (!state.selectedNodeIndices.has(nodeIndex)) {
                    state.selectedNodeIndices.add(nodeIndex);
                    changed = true;
                }
            });
        });
        if (changed) {
            renderClusterView();
            updateMetadataTable();
        }
    }

    function selectNodesInBox(rect, deselect) {
        let changed = false;
        state.visibleLayout.forEach(item => {
            if (!itemIntersectsSelectionBox(item, rect)) {
                return;
            }
            const component = state.bundle.views.ssn.hierarchy.nodes[item.componentId];
            const dotLayout = componentMemberLayout(component, item);
            dotLayout.forEach(dot => {
                const dotScreen = worldToScreenPoint(dot.x, dot.y);
                const dotRadius = dot.radius * state.viewTransform.scale;
                if (!circleIntersectsRect(dotScreen.x, dotScreen.y, dotRadius, rect)) {
                    return;
                }
                if (deselect) {
                    if (state.selectedNodeIndices.delete(dot.memberIndex)) {
                        changed = true;
                    }
                } else if (!state.selectedNodeIndices.has(dot.memberIndex)) {
                    state.selectedNodeIndices.add(dot.memberIndex);
                    changed = true;
                }
            });
        });
        if (changed) {
            renderClusterView();
            updateMetadataTable();
        }
    }
