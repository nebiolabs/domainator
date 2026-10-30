    // ---- Tree view: drawing ----
    //
    // drawTree and drawMsa draw onto any canvas-like context (see createSvgContext), so the
    // screen, the PNG export and the SVG export are one drawing. Both take a geometry from
    // treeGeometry, which places the rows: the alignment's rows are the tree's rows.

    const TREE_PADDING_LEFT = 14;
    const TREE_PADDING_RIGHT = 8;
    const TREE_PADDING_BOTTOM = 30;
    const TREE_HEADER_PLAIN = 10;
    // The alignment's ruler and conservation track, above the rows of both panes.
    const MSA_HEADER = 40;
    const MSA_RULER_HEIGHT = 16;
    const TREE_LABEL_MIN_ROW = 7;
    const TREE_DOT_MIN_ROW = 2.5;
    const TREE_SUPPORT_MIN_ROW = 10;
    const MSA_LETTER_MIN_SIZE = 7;
    const TREE_LABEL_FONT_MAX = 12;
    const TREE_COLOR_STRIP = 8;
    const TREE_SELECTION_BAND = 'rgba(66, 133, 244, 0.16)';
    const TREE_SELECTION_BRANCH = 'rgba(66, 133, 244, 0.45)';

    function treeHeaderHeight() {
        return state.tree.msa ? MSA_HEADER : TREE_HEADER_PLAIN;
    }

    function treeLeafLabel(bundleIndex) {
        const field = currentLabelField();
        if (field && field !== '__node_id__') {
            const value = metadataValue(bundleIndex, field);
            return value === null || value === undefined ? '' : String(value);
        }
        return nodeId(bundleIndex);
    }

    // What a row is called: a leaf's label, or a collapsed clade's label and size.
    function treeRowLabel(node) {
        const tree = state.tree;
        if (tree.children[node].length === 0) {
            return treeLeafLabel(tree.treeToBundle[node]);
        }
        const label = tree.view.label[node];
        return (label ? label + ' ' : '') + '(' + tree.leafCount[node].toLocaleString() + ' leaves)';
    }

    function treeShowLabels() {
        return document.getElementById('tree-show-labels').checked;
    }

    function treeShowSupport() {
        return document.getElementById('tree-show-support').checked;
    }

    function treeViewHeight(height) {
        return Math.max(1, height - treeHeaderHeight() - TREE_PADDING_BOTTOM);
    }

    // The row height that shows every row at once.
    function treeFitRowHeight(height) {
        return Math.max(TREE_ROW_HEIGHT_MIN, treeViewHeight(height) / Math.max(1, state.tree.layout.rowNode.length));
    }

    // Rows can shrink until the whole tree fits (or to the default height, for a tree that fits
    // with room to spare) and grow to TREE_ROW_HEIGHT_MAX; the scroll stays within the rows.
    function treeClampView(height) {
        const tree = state.tree;
        const smallest = Math.min(treeFitRowHeight(height), TREE_ROW_HEIGHT_DEFAULT);
        tree.rowHeight = Math.max(smallest, Math.min(TREE_ROW_HEIGHT_MAX, tree.rowHeight));
        const maxScroll = Math.max(0, tree.layout.rowNode.length * tree.rowHeight - treeViewHeight(height));
        tree.scrollY = Math.max(0, Math.min(maxScroll, tree.scrollY));
    }

    function treeGeometry(width, height, options = {}) {
        const tree = state.tree;
        const layout = tree.layout;
        const top = treeHeaderHeight();
        const viewHeight = treeViewHeight(height);
        const rowHeight = tree.rowHeight;
        const labels = treeShowLabels() && rowHeight >= TREE_LABEL_MIN_ROW;
        const fontSize = Math.min(TREE_LABEL_FONT_MAX, rowHeight * 0.8);
        let labelWidth = TREE_COLOR_STRIP + 10;
        if (labels) {
            let longest = 0;
            layout.rowNode.forEach(node => { longest = Math.max(longest, treeRowLabel(node).length); });
            labelWidth = Math.min(width * 0.45, longest * fontSize * 0.6 + 16);
        }
        const treeWidth = options.treeWidth ?? Math.max(40, width - TREE_PADDING_LEFT - TREE_PADDING_RIGHT - labelWidth);
        const left = options.left ?? TREE_PADDING_LEFT;
        const scrollY = options.scrollY ?? tree.scrollY;
        const rowCount = layout.rowNode.length;
        return {
            width, height, top, viewHeight, rowHeight, labels, fontSize, treeWidth,
            xOf: value => left + value / layout.maxX * treeWidth,
            xPx: node => left + layout.x[node] / layout.maxX * treeWidth,
            yPx: row => top + (row + 0.5) * rowHeight - scrollY,
            rowStart: Math.max(0, Math.floor(scrollY / rowHeight) - 1),
            rowEnd: Math.min(rowCount - 1, Math.ceil((scrollY + viewHeight) / rowHeight) + 1),
            msaScrollX: options.msaScrollX ?? tree.msaScrollX,
            preview: options.preview !== false,
        };
    }

    // Per node: the color every leaf below it shares (null when they differ or nothing is
    // colored), how many of its leaves are selected, and how many are in the previewed preset.
    function treeCladeSummaries(showPreview = true) {
        const tree = state.tree;
        const coloring = Boolean(colorInfo(currentColorField()));
        const preview = showPreview ? presetPreviewNodeSet() : null;
        const color = new Array(tree.count).fill(null);
        const selected = new Int32Array(tree.count);
        const previewed = new Int32Array(tree.count);
        for (let node = tree.count - 1; node >= 0; node--) {
            const children = tree.children[node];
            if (children.length === 0) {
                const bundleIndex = tree.treeToBundle[node];
                color[node] = coloring ? state.nodeColorCache[bundleIndex] : null;
                selected[node] = state.selectedNodeIndices.has(bundleIndex) ? 1 : 0;
                previewed[node] = preview && preview.has(bundleIndex) ? 1 : 0;
                continue;
            }
            let shared = color[children[0]];
            children.forEach(child => {
                if (color[child] !== shared) {
                    shared = null;
                }
                selected[node] += selected[child];
                previewed[node] += previewed[child];
            });
            color[node] = shared;
        }
        return {coloring, color, selected, previewed};
    }

    // Draws the row backgrounds of `summary`'s selected and previewed rows across `width`.
    function drawTreeRowBands(ctx, geometry, summary, width) {
        const tree = state.tree;
        const {rowNode} = tree.layout;
        for (let row = geometry.rowStart; row <= geometry.rowEnd; row++) {
            const node = rowNode[row];
            const y = geometry.yPx(row) - geometry.rowHeight / 2;
            if (summary.selected[node] > 0 && summary.selected[node] === tree.leafCount[node]) {
                ctx.fillStyle = TREE_SELECTION_BAND;
                ctx.fillRect(0, y, width, geometry.rowHeight);
            }
            if (summary.previewed[node] > 0) {
                ctx.strokeStyle = VIEW_PREVIEW_OUTLINE;
                ctx.lineWidth = 1;
                ctx.strokeRect(0.5, y + 0.5, width - 1, Math.max(1, geometry.rowHeight - 1));
            }
        }
    }

    function drawTree(ctx, geometry, summary) {
        const tree = state.tree;
        const layout = tree.layout;
        const view = tree.view;
        const {width, height, rowHeight, rowStart, rowEnd, xPx, yPx} = geometry;
        ctx.save();
        ctx.fillStyle = '#ffffff';
        ctx.fillRect(0, 0, width, height);
        drawTreeRowBands(ctx, geometry, summary, width);

        // Branches, batched by color: a clade's branches take the color all its leaves share.
        const segments = new Map();
        const selectedSegments = [];
        const addSegment = (color, x1, y1, x2, y2, selected) => {
            if (!segments.has(color)) {
                segments.set(color, []);
            }
            segments.get(color).push(x1, y1, x2, y2);
            if (selected) {
                selectedSegments.push(x1, y1, x2, y2);
            }
        };
        layout.visibleOrder.forEach(node => {
            const row = layout.nodeRow[node];
            const color = summary.color[node] || VIEW_INK;
            const selected = summary.selected[node] > 0 && summary.selected[node] === tree.leafCount[node];
            const parent = view.parent[node];
            if (parent >= 0 && row >= rowStart - 1 && row <= rowEnd + 1) {
                addSegment(color, xPx(parent), yPx(row), xPx(node), yPx(row), selected);
            }
            const children = layout.displayChildren[node];
            if (children.length > 0 && !tree.collapsed.has(node)) {
                const firstRow = layout.nodeRow[children[0]];
                const lastRow = layout.nodeRow[children[children.length - 1]];
                if (lastRow >= rowStart - 1 && firstRow <= rowEnd + 1) {
                    addSegment(color, xPx(node), yPx(firstRow), xPx(node), yPx(lastRow), selected);
                }
            }
        });
        ctx.lineCap = 'round';
        if (selectedSegments.length > 0) {
            ctx.strokeStyle = TREE_SELECTION_BRANCH;
            ctx.lineWidth = 5;
            ctx.beginPath();
            for (let i = 0; i < selectedSegments.length; i += 4) {
                ctx.moveTo(selectedSegments[i], selectedSegments[i + 1]);
                ctx.lineTo(selectedSegments[i + 2], selectedSegments[i + 3]);
            }
            ctx.stroke();
        }
        ctx.lineWidth = rowHeight < 2 ? 0.75 : 1.25;
        segments.forEach((points, color) => {
            ctx.strokeStyle = color;
            ctx.beginPath();
            for (let i = 0; i < points.length; i += 4) {
                ctx.moveTo(points[i], points[i + 1]);
                ctx.lineTo(points[i + 2], points[i + 3]);
            }
            ctx.stroke();
        });

        // Rows: a collapsed clade's wedge, a leaf's dot, the colors and the labels.
        const radius = Math.min(4, Math.max(1.2, rowHeight * 0.3));
        const stripX = geometry.xOf(layout.maxX) + 6;
        ctx.font = geometry.fontSize + 'px ' + VIEW_FONT;
        ctx.textBaseline = 'middle';
        ctx.textAlign = 'left';
        for (let row = rowStart; row <= rowEnd; row++) {
            const node = layout.rowNode[row];
            const y = yPx(row);
            const x = xPx(node);
            const leaf = tree.children[node].length === 0;
            const color = summary.color[node];
            let labelX = x + radius + 4;
            if (!leaf) {
                const end = geometry.xOf(layout.cladeMaxX[node]);
                const half = Math.max(0.5, rowHeight * 0.4);
                ctx.fillStyle = color || '#9aa3a8';
                ctx.globalAlpha = 0.7;
                ctx.beginPath();
                ctx.moveTo(x, y);
                ctx.lineTo(Math.max(end, x + 8), y - half);
                ctx.lineTo(Math.max(end, x + 8), y + half);
                ctx.closePath();
                ctx.fill();
                ctx.globalAlpha = 1;
                labelX = Math.max(end, x + 8) + 4;
            } else if (rowHeight >= TREE_DOT_MIN_ROW) {
                ctx.fillStyle = summary.coloring ? color : VIEW_MUTED;
                ctx.beginPath();
                ctx.arc(x, y, radius, 0, 2 * Math.PI);
                ctx.fill();
            }
            if (geometry.labels) {
                ctx.fillStyle = VIEW_INK;
                ctx.fillText(treeRowLabel(node), labelX, y);
            } else if (summary.coloring) {
                ctx.fillStyle = color || '#9aa3a8';
                ctx.fillRect(stripX, y - rowHeight / 2, TREE_COLOR_STRIP, Math.max(rowHeight, 1));
            }
        }

        // Support values and internal labels.
        if (treeShowSupport() && tree.hasSupport && rowHeight >= TREE_SUPPORT_MIN_ROW) {
            ctx.font = '9px ' + VIEW_FONT;
            ctx.fillStyle = VIEW_MUTED;
            ctx.textAlign = 'right';
            ctx.textBaseline = 'bottom';
            layout.visibleOrder.forEach(node => {
                if (tree.children[node].length === 0 || tree.collapsed.has(node)) {
                    return;
                }
                const row = layout.nodeRow[node];
                if (row < rowStart || row > rowEnd) {
                    return;
                }
                const text = view.label[node] || (view.support[node] === null || view.support[node] === undefined
                    ? '' : formatTreeLength(view.support[node]));
                if (text) {
                    ctx.fillText(text, xPx(node) - 3, yPx(row) - 2);
                }
            });
        }

        // The bottom margin: a scale bar for branch lengths.
        ctx.fillStyle = '#ffffff';
        ctx.fillRect(0, height - TREE_PADDING_BOTTOM, width, TREE_PADDING_BOTTOM);
        ctx.fillRect(0, 0, width, geometry.top);
        if (treeLayoutMode() === 'phylogram' && layout.maxX > 0) {
            const length = niceLength(layout.maxX * 0.2);
            const barY = height - 12;
            ctx.strokeStyle = VIEW_INK;
            ctx.lineWidth = 1.5;
            ctx.beginPath();
            ctx.moveTo(geometry.xOf(0), barY);
            ctx.lineTo(geometry.xOf(length), barY);
            ctx.stroke();
            ctx.font = '11px ' + VIEW_FONT;
            ctx.fillStyle = VIEW_INK;
            ctx.textAlign = 'left';
            ctx.textBaseline = 'middle';
            ctx.fillText(formatTreeLength(length), geometry.xOf(length) + 6, barY);
        }

        const cutX = treeCutX();
        if (Number.isFinite(cutX) && cutX > 0) {
            const x = geometry.xOf(cutX);
            ctx.strokeStyle = VIEW_ACCENT;
            ctx.lineWidth = 1.5;
            ctx.setLineDash([5, 4]);
            ctx.beginPath();
            ctx.moveTo(x, geometry.top);
            ctx.lineTo(x, height - TREE_PADDING_BOTTOM);
            ctx.stroke();
            ctx.setLineDash([]);
        }
        ctx.restore();
    }

    // ---- The alignment pane ----

    function msaColumnWidth() {
        const width = Number(document.getElementById('msa-column-width').value);
        return Number.isFinite(width) && width > 0 ? Math.min(40, Math.max(1, width)) : 10;
    }

    function msaClampScroll(width) {
        const tree = state.tree;
        const maxScroll = Math.max(0, tree.msa.length * msaColumnWidth() - width);
        tree.msaScrollX = Math.max(0, Math.min(maxScroll, tree.msaScrollX));
    }

    function msaColumnRange(width, scrollX = state.tree.msaScrollX) {
        const columnWidth = msaColumnWidth();
        const length = state.tree.msa.length;
        return {
            columnWidth,
            start: Math.max(0, Math.floor(scrollX / columnWidth)),
            end: Math.min(length - 1, Math.ceil((scrollX + width) / columnWidth)),
            xOf: column => column * columnWidth - scrollX,
        };
    }

    function drawMsa(ctx, geometry, width, summary) {
        const tree = state.tree;
        const {rowHeight, rowStart, rowEnd, top, height} = geometry;
        const columns = msaColumnRange(width, geometry.msaScrollX);
        const {columnWidth, start, end, xOf} = columns;
        const scheme = msaScheme();
        ctx.save();
        ctx.fillStyle = '#ffffff';
        ctx.fillRect(0, 0, width, height);

        // Rows narrower than a pixel are drawn a few at a time, from the first of each group.
        const stride = rowHeight < 1 ? Math.ceil(1 / rowHeight) : 1;
        const letters = columnWidth >= MSA_LETTER_MIN_SIZE && rowHeight >= 8;
        for (let row = rowStart - (rowStart % stride); row <= rowEnd; row += stride) {
            const node = tree.layout.rowNode[row];
            const sequence = msaRowSequence(node);
            const y = geometry.yPx(row) - rowHeight / 2;
            const rowSpan = rowHeight * stride;
            if (sequence === null) {
                ctx.fillStyle = '#fafafa';
                ctx.fillRect(0, y, width, rowSpan);
                continue;
            }
            // Runs of one color are one rectangle.
            let runColor = null;
            let runStart = start;
            for (let column = start; column <= end + 1; column++) {
                const color = column <= end ? msaResidueColor(sequence[column], column, scheme) : null;
                if (color !== runColor) {
                    if (runColor !== null) {
                        ctx.fillStyle = runColor;
                        ctx.fillRect(xOf(runStart), y, (column - runStart) * columnWidth, rowSpan);
                    }
                    runColor = color;
                    runStart = column;
                }
            }
            if (letters) {
                ctx.fillStyle = VIEW_INK;
                ctx.font = Math.min(columnWidth, rowHeight) * 0.75 + 'px ui-monospace, Menlo, Consolas, monospace';
                ctx.textAlign = 'center';
                ctx.textBaseline = 'middle';
                for (let column = start; column <= end; column++) {
                    if (!MSA_GAPS.has(sequence[column])) {
                        ctx.fillText(sequence[column], xOf(column) + columnWidth / 2, y + rowHeight / 2);
                    }
                }
            }
        }
        drawTreeRowBands(ctx, geometry, summary, width);

        // The header: column numbers, and how conserved each column is.
        ctx.fillStyle = '#ffffff';
        ctx.fillRect(0, 0, width, top);
        ctx.fillRect(0, height - TREE_PADDING_BOTTOM, width, TREE_PADDING_BOTTOM);
        const consensus = msaConsensus();
        const trackHeight = top - MSA_RULER_HEIGHT - 4;
        ctx.fillStyle = MSA_IDENTICAL_COLOR;
        for (let column = start; column <= end; column++) {
            const barHeight = consensus.identity[column] * trackHeight;
            ctx.fillRect(xOf(column), top - 2 - barHeight, Math.max(columnWidth - (columnWidth > 3 ? 1 : 0), 0.5), barHeight);
        }
        const step = [1, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000].find(value => value * columnWidth >= 40) || 10000;
        ctx.fillStyle = VIEW_MUTED;
        ctx.strokeStyle = VIEW_MUTED;
        ctx.lineWidth = 1;
        ctx.font = '10px ' + VIEW_FONT;
        ctx.textAlign = 'center';
        ctx.textBaseline = 'top';
        for (let position = Math.ceil((start + 1) / step) * step; position <= end + 1; position += step) {
            const x = xOf(position - 1) + columnWidth / 2;
            ctx.fillText(String(position), x, 1);
            ctx.beginPath();
            ctx.moveTo(x, MSA_RULER_HEIGHT - 3);
            ctx.lineTo(x, MSA_RULER_HEIGHT);
            ctx.stroke();
        }

        // Where the view is along the alignment.
        const total = tree.msa.length * columnWidth;
        if (total > width) {
            const barY = height - 14;
            ctx.fillStyle = '#e8ebef';
            ctx.fillRect(4, barY, width - 8, 6);
            ctx.fillStyle = '#9aa3a8';
            ctx.fillRect(4 + (width - 8) * tree.msaScrollX / total, barY, Math.max(6, (width - 8) * width / total), 6);
        }
        ctx.restore();
    }

    // ---- Painting the panes ----

    function renderTreeView() {
        const tree = state.tree;
        if (!tree || !tree.layout || state.activeViewName !== 'tree') {
            return;
        }
        if (tree.pendingFrame !== null) {
            cancelAnimationFrame(tree.pendingFrame);
            tree.pendingFrame = null;
        }
        const treePane = viewCanvasContext(document.getElementById('tree-canvas'));
        if (!treePane) {
            return;
        }
        if (tree.rowHeight === null) {
            tree.rowHeight = Math.min(TREE_ROW_HEIGHT_DEFAULT, treeFitRowHeight(treePane.height));
        }
        treeClampView(treePane.height);
        const summary = treeCladeSummaries();
        drawTree(treePane.ctx, treeGeometry(treePane.width, treePane.height), summary);
        if (tree.msa) {
            const msaPane = viewCanvasContext(document.getElementById('msa-canvas'));
            if (msaPane) {
                msaClampScroll(msaPane.width);
                drawMsa(msaPane.ctx, treeGeometry(treePane.width, msaPane.height), msaPane.width, summary);
            }
        }
        const rowHeightInput = document.getElementById('tree-row-height');
        if (document.activeElement !== rowHeightInput) {
            rowHeightInput.value = String(Number(tree.rowHeight.toPrecision(3)));
        }
        drawTreeSelectionBox();
    }

    function scheduleTreeRender() {
        const tree = state.tree;
        if (tree && tree.pendingFrame === null) {
            tree.pendingFrame = requestAnimationFrame(() => {
                tree.pendingFrame = null;
                renderTreeView();
            });
        }
    }

    // ---- Exports ----

    // Both panes as they are shown, drawn onto `ctx` side by side.
    function drawTreeViewExport(ctx, treeWidth, msaWidth, height, options = {}) {
        const summary = treeCladeSummaries(options.preview !== false);
        const geometry = treeGeometry(treeWidth, height, options);
        drawTree(ctx, geometry, summary);
        if (msaWidth > 0) {
            ctx.save();
            ctx.translate(treeWidth, 0);
            drawMsa(ctx, geometry, msaWidth, summary);
            ctx.restore();
        }
    }

    function treeExportSize() {
        const treeRect = document.getElementById('tree-canvas').getBoundingClientRect();
        const msaRect = state.tree.msa ? document.getElementById('msa-canvas').getBoundingClientRect() : {width: 0};
        const size = {treeWidth: treeRect.width, msaWidth: msaRect.width, height: treeRect.height, options: {}};
        if (entireContentExport()) {
            const tree = state.tree;
            const geometry = treeGeometry(treeRect.width, treeRect.height);
            const measurer = document.createElement('canvas').getContext('2d');
            if (!measurer) { throw new Error('Could not measure tree labels.'); }
            let labelWidth = TREE_COLOR_STRIP + 10;
            if (geometry.labels) {
                measurer.font = geometry.fontSize + 'px ' + VIEW_FONT;
                for (const node of tree.layout.rowNode) {
                    labelWidth = Math.max(labelWidth, measurer.measureText(treeRowLabel(node)).width + 24);
                }
            }
            let left = TREE_PADDING_LEFT;
            if (treeShowSupport() && tree.rowHeight >= TREE_SUPPORT_MIN_ROW) {
                measurer.font = '9px ' + VIEW_FONT;
                for (const node of tree.layout.visibleOrder) {
                    const text = tree.view.label[node] || (tree.view.support[node] === null ? '' : formatTreeLength(tree.view.support[node]));
                    left = Math.max(left, measurer.measureText(text).width + 8);
                }
            }
            size.treeWidth = Math.ceil(left + geometry.treeWidth + labelWidth + TREE_PADDING_RIGHT);
            measurer.font = '10px ' + VIEW_FONT;
            size.msaWidth = tree.msa ? Math.ceil(tree.msa.length * msaColumnWidth() +
                measurer.measureText(String(tree.msa.length)).width / 2 + 4) : 0;
            size.height = Math.ceil(treeHeaderHeight() + tree.layout.rowNode.length * tree.rowHeight + TREE_PADDING_BOTTOM);
            size.options = {treeWidth: geometry.treeWidth, left, scrollY: 0, msaScrollX: 0, preview: false};
        }
        size.width = size.treeWidth + size.msaWidth;
        validateImageExportSize(Math.round(size.width), Math.round(size.height));
        return size;
    }

    function exportTreePNG() {
        if (!state.tree) { return; }
        const size = treeExportSize();
        return exportCanvasPNG(size, ctx => drawTreeViewExport(ctx, size.treeWidth, size.msaWidth, size.height, size.options));
    }

    function buildTreeViewSVG() {
        const {treeWidth, msaWidth, height, options} = treeExportSize();
        const ctx = createSvgContext(Math.round(treeWidth + msaWidth), Math.round(height));
        drawTreeViewExport(ctx, treeWidth, msaWidth, height, options);
        return ctx.toSVG();
    }

    function exportTreeSVG() {
        if (!state.tree) { return; }
        const blob = new Blob([buildTreeViewSVG()], {type: 'image/svg+xml'});
        triggerDownload(URL.createObjectURL(blob), exportBaseName() + imageExportSuffix('svg'), true);
    }
