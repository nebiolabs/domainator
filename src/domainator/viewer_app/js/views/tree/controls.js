    // ---- Tree view: pointer and controls ----

    // A press that moves less than this many pixels is a click, not a drag.
    const TREE_CLICK_SLOP = 3;
    // How far from a branch or a node, in pixels, a click still lands on it.
    const TREE_HIT_SLOP = 5;
    const TREE_LINE_PIXELS = 16;

    function treePanePoint(canvas, event) {
        const rect = canvas.getBoundingClientRect();
        return {x: event.clientX - rect.left, y: event.clientY - rect.top, width: rect.width, height: rect.height};
    }

    // The row under a point, or -1 above or below the rows.
    function treeRowAt(point, geometry) {
        if (point.y < geometry.top || point.y > point.height - TREE_PADDING_BOTTOM) {
            return -1;
        }
        const row = Math.floor((point.y - geometry.top + state.tree.scrollY) / geometry.rowHeight);
        return row >= 0 && row < state.tree.layout.rowNode.length ? row : -1;
    }

    // The node a point of the tree pane is on: a row's leaf or collapsed clade when the point is
    // right of its tip, an internal node at the point, the branch under the point, or the clade
    // whose vertical line the point is beside; null when none is.
    function treeHitTest(point) {
        const tree = state.tree;
        const layout = tree.layout;
        const geometry = treeGeometry(point.width, point.height);
        const row = treeRowAt(point, geometry);
        if (row >= 0 && point.x >= geometry.xPx(layout.rowNode[row]) - TREE_HIT_SLOP) {
            return layout.rowNode[row];
        }
        const internalNodes = layout.visibleOrder.filter(node => layout.displayChildren[node].length > 0 && !tree.collapsed.has(node));
        let best = null;
        let bestDistance = TREE_HIT_SLOP;
        internalNodes.forEach(node => {
            const distance = Math.hypot(geometry.xPx(node) - point.x, geometry.yPx(layout.nodeRow[node]) - point.y);
            if (distance <= bestDistance) {
                best = node;
                bestDistance = distance;
            }
        });
        if (best !== null) {
            return best;
        }
        if (row >= 0) {
            for (let ancestor = layout.rowNode[row]; ancestor > 0; ancestor = tree.view.parent[ancestor]) {
                const parentX = geometry.xPx(tree.view.parent[ancestor]);
                if (Math.abs(geometry.yPx(layout.nodeRow[ancestor]) - point.y) <= TREE_HIT_SLOP &&
                        point.x >= parentX && point.x <= geometry.xPx(ancestor) + TREE_HIT_SLOP) {
                    return ancestor;
                }
            }
        }
        bestDistance = TREE_HIT_SLOP;
        internalNodes.forEach(node => {
            const children = layout.displayChildren[node];
            const distance = Math.abs(geometry.xPx(node) - point.x);
            const top = geometry.yPx(layout.nodeRow[children[0]]);
            const bottom = geometry.yPx(layout.nodeRow[children[children.length - 1]]);
            if (distance <= bestDistance && point.y >= top - TREE_HIT_SLOP && point.y <= bottom + TREE_HIT_SLOP) {
                best = node;
                bestDistance = distance;
            }
        });
        return best;
    }

    function msaHitTest(point) {
        const geometry = treeGeometry(point.width, point.height);
        const row = treeRowAt(point, geometry);
        return row >= 0 ? state.tree.layout.rowNode[row] : null;
    }

    // Selects every leaf of a clade, or deselects them all when they all are selected.
    function treeToggleClade(node) {
        const leaves = treeCladeLeaves(node);
        const allSelected = leaves.every(index => state.selectedNodeIndices.has(index));
        leaves.forEach(index => {
            if (allSelected) {
                state.selectedNodeIndices.delete(index);
            } else {
                state.selectedNodeIndices.add(index);
            }
        });
        renderTreeView();
        updateMetadataTable();
    }

    function treeToggleCollapsed(node) {
        const tree = state.tree;
        if (node === 0 || tree.children[node].length === 0) {
            return;
        }
        if (tree.collapsed.has(node)) {
            tree.collapsed.delete(node);
        } else {
            tree.collapsed.add(node);
        }
        treeComputeLayout();
        renderTreeView();
    }

    // Selects (or with `deselect`, deselects) the leaves of the rows a box spans. In the tree
    // pane the box must also reach a row's tip.
    function treeSelectInBox(box, inTreePane) {
        const tree = state.tree;
        const geometry = treeGeometry(box.width, box.height);
        const y0 = Math.min(box.y0, box.y1);
        const y1 = Math.max(box.y0, box.y1);
        const x1 = Math.max(box.x0, box.x1);
        let changed = false;
        for (let row = geometry.rowStart; row <= geometry.rowEnd; row++) {
            const y = geometry.yPx(row);
            const node = tree.layout.rowNode[row];
            if (y < y0 || y > y1 || (inTreePane && x1 < geometry.xPx(node) - TREE_HIT_SLOP)) {
                continue;
            }
            treeCladeLeaves(node).forEach(index => {
                if (box.deselect) {
                    changed = state.selectedNodeIndices.delete(index) || changed;
                } else if (!state.selectedNodeIndices.has(index)) {
                    state.selectedNodeIndices.add(index);
                    changed = true;
                }
            });
        }
        if (changed) {
            updateMetadataTable();
        }
    }

    function treeTooltipText(node, column) {
        const tree = state.tree;
        const lines = [];
        if (tree.children[node].length === 0) {
            const bundleIndex = tree.treeToBundle[node];
            lines.push(nodeId(bundleIndex));
            [currentLabelField(), currentColorField()].forEach(field => {
                if (field && field !== '__node_id__') {
                    lines.push(field + ': ' + formatValue(metadataValue(bundleIndex, field)));
                }
            });
        } else {
            lines.push((tree.view.label[node] ? tree.view.label[node] + ': ' : '') + 'clade of ' +
                tree.leafCount[node].toLocaleString() + ' leaves' + (tree.collapsed.has(node) ? ' (collapsed)' : ''));
            if (tree.view.support[node] !== null && tree.view.support[node] !== undefined) {
                lines.push('support: ' + formatTreeLength(tree.view.support[node]));
            }
        }
        if (column !== null) {
            const sequence = msaRowSequence(node);
            lines.push('column ' + (column + 1).toLocaleString() + ': ' + (sequence === null ? 'no sequence' : sequence[column]));
        }
        return lines.join('\n');
    }

    // Zooms the rows around a point of either pane, keeping the row under it in place.
    function treeZoomRows(factor, point) {
        const tree = state.tree;
        const top = treeHeaderHeight();
        const anchor = (tree.scrollY + point.y - top) / tree.rowHeight;
        tree.rowHeight = tree.rowHeight * factor;
        treeClampView(point.height);
        tree.scrollY = anchor * tree.rowHeight - (point.y - top);
    }

    function msaZoomColumns(factor, point) {
        const tree = state.tree;
        const anchor = (tree.msaScrollX + point.x) / msaColumnWidth();
        const width = Math.min(40, Math.max(1, msaColumnWidth() * factor));
        document.getElementById('msa-column-width').value = String(Number(width.toPrecision(3)));
        tree.msaScrollX = anchor * width - point.x;
    }

    function attachTreePointer(canvas, inTreePane) {
        let press = null;
        canvas.addEventListener('pointerdown', event => {
            if (event.button !== 0 || !state.tree || !state.tree.layout) {
                return;
            }
            const point = treePanePoint(canvas, event);
            press = {
                id: event.pointerId, point, moved: false,
                box: event.shiftKey, deselect: event.altKey,
                scrollY: state.tree.scrollY, scrollX: state.tree.msaScrollX,
            };
        });
        canvas.addEventListener('pointermove', event => {
            if (!state.tree || !state.tree.layout) {
                return;
            }
            const point = treePanePoint(canvas, event);
            if (press && event.pointerId === press.id) {
                const dx = point.x - press.point.x;
                const dy = point.y - press.point.y;
                if (!press.moved && Math.hypot(dx, dy) >= TREE_CLICK_SLOP) {
                    press.moved = true;
                    canvas.setPointerCapture(event.pointerId);
                    if (!press.box) {
                        canvas.classList.add('panning');
                    }
                    hideViewTip();
                }
                if (press.moved) {
                    if (press.box) {
                        state.tree.box = {canvas, x0: press.point.x, y0: press.point.y, x1: point.x, y1: point.y};
                    } else {
                        state.tree.scrollY = press.scrollY - dy;
                        if (!inTreePane) {
                            state.tree.msaScrollX = press.scrollX - dx;
                        }
                    }
                    scheduleTreeRender();
                    return;
                }
            }
            const node = inTreePane ? treeHitTest(point) : msaHitTest(point);
            if (node === null) {
                hideViewTip();
                return;
            }
            let column = null;
            if (!inTreePane) {
                column = Math.floor((state.tree.msaScrollX + point.x) / msaColumnWidth());
                if (column < 0 || column >= state.tree.msa.length) {
                    column = null;
                }
            }
            showViewTip(event, treeTooltipText(node, column));
        });
        const release = event => {
            if (!press || event.pointerId !== press.id) {
                return;
            }
            const point = treePanePoint(canvas, event);
            if (press.moved && press.box) {
                treeSelectInBox({
                    x0: press.point.x, y0: press.point.y, x1: point.x, y1: point.y,
                    width: point.width, height: point.height, deselect: press.deselect,
                }, inTreePane);
                state.tree.box = null;
                renderTreeView();
            } else if (!press.moved && event.type === 'pointerup') {
                const node = inTreePane ? treeHitTest(point) : msaHitTest(point);
                if (node !== null) {
                    treeToggleClade(node);
                }
            }
            press = null;
            canvas.classList.remove('panning');
        };
        canvas.addEventListener('pointerup', release);
        canvas.addEventListener('pointercancel', release);
        canvas.addEventListener('pointerleave', hideViewTip);
        canvas.addEventListener('dblclick', event => {
            if (inTreePane && state.tree && state.tree.layout) {
                const node = treeHitTest(treePanePoint(canvas, event));
                if (node !== null) {
                    treeToggleCollapsed(node);
                }
            }
        });
        canvas.addEventListener('wheel', event => {
            if (!state.tree || !state.tree.layout) {
                return;
            }
            event.preventDefault();
            const point = treePanePoint(canvas, event);
            const unit = event.deltaMode === 1 ? TREE_LINE_PIXELS : event.deltaMode === 2 ? point.height : 1;
            const deltaX = event.deltaX * unit;
            const deltaY = event.deltaY * unit;
            if (event.ctrlKey || event.metaKey) {
                const factor = Math.pow(2, -deltaY * 0.004);
                if (inTreePane) {
                    treeZoomRows(factor, point);
                } else {
                    msaZoomColumns(factor, point);
                }
            } else if (!inTreePane && (event.shiftKey || Math.abs(deltaX) > Math.abs(deltaY))) {
                state.tree.msaScrollX += event.shiftKey ? deltaY : deltaX;
            } else {
                state.tree.scrollY += deltaY;
            }
            hideViewTip();
            scheduleTreeRender();
        }, {passive: false});
    }

    // Draws the box being dragged, over the pane it is dragged in.
    function drawTreeSelectionBox() {
        const box = state.tree && state.tree.box;
        if (!box) {
            return;
        }
        const pane = viewCanvasContext(box.canvas);
        if (!pane) {
            return;
        }
        pane.ctx.save();
        pane.ctx.strokeStyle = VIEW_ACCENT;
        pane.ctx.fillStyle = 'rgba(200, 85, 61, 0.08)';
        pane.ctx.setLineDash([4, 3]);
        const x = Math.min(box.x0, box.x1);
        const y = Math.min(box.y0, box.y1);
        pane.ctx.fillRect(x, y, Math.abs(box.x1 - box.x0), Math.abs(box.y1 - box.y0));
        pane.ctx.strokeRect(x, y, Math.abs(box.x1 - box.x0), Math.abs(box.y1 - box.y0));
        pane.ctx.restore();
    }

    // Scrolls the first selected row into the middle of the view.
    function treeFocusSelection() {
        const tree = state.tree;
        if (!tree) { return; }
        const summary = treeCladeSummaries();
        const row = tree.layout.rowNode.findIndex(node => summary.selected[node] > 0);
        if (row < 0) {
            setStatus('No leaves are selected.');
            return;
        }
        const height = document.getElementById('tree-canvas').getBoundingClientRect().height;
        tree.scrollY = (row + 0.5) * tree.rowHeight - treeViewHeight(height) / 2;
        renderTreeView();
    }

    function treeResetView() {
        const tree = state.tree;
        if (!tree) { return; }
        tree.rowHeight = null;
        tree.scrollY = 0;
        tree.msaScrollX = 0;
        renderTreeView();
    }

    function downloadTreeText(text, suffix, type) {
        triggerDownload(URL.createObjectURL(new Blob([text], {type})), exportBaseName() + suffix, true);
    }

    function setupTreeViewControls() {
        attachTreePointer(document.getElementById('tree-canvas'), true);
        attachTreePointer(document.getElementById('msa-canvas'), false);
        const relayout = () => {
            if (!state.tree) { return; }
            treeComputeLayout();
            updateTreeCutLabel();
            renderTreeView();
        };
        document.getElementById('tree-layout').addEventListener('change', relayout);
        document.getElementById('tree-ladderize').addEventListener('change', relayout);
        ['tree-show-labels', 'tree-show-support', 'msa-scheme', 'msa-column-width'].forEach(id => {
            document.getElementById(id).addEventListener('change', () => renderTreeView());
        });
        document.getElementById('tree-row-height').addEventListener('change', event => {
            const value = Number(event.target.value);
            if (state.tree && Number.isFinite(value) && value > 0) {
                state.tree.rowHeight = value;
            }
            renderTreeView();
        });
        document.getElementById('tree-cut').addEventListener('input', () => {
            if (!state.tree) { return; }
            updateTreeCutLabel();
            renderTreeView();
        });
        document.getElementById('tree-expand-all').addEventListener('click', () => {
            if (!state.tree) { return; }
            state.tree.collapsed.clear();
            relayout();
        });
        document.getElementById('export-newick').addEventListener('click', () => {
            if (!state.tree) { return; }
            const selection = state.selectedNodeIndices.size > 0 ? state.selectedNodeIndices : null;
            downloadTreeText(treeNewick(selection) + '\n', selection ? '_selected.nwk' : '.nwk', 'text/plain');
        });
        document.getElementById('export-fasta').addEventListener('click', () => {
            if (!state.tree || !state.tree.msa) { return; }
            const ungapped = document.getElementById('export-fasta-gaps').value === 'ungapped';
            const suffix = (state.selectedNodeIndices.size > 0 ? '_selected' : '') + (ungapped ? '.fasta' : '.afa');
            downloadTreeText(msaFasta(ungapped), suffix, 'text/plain');
        });
        window.addEventListener('resize', () => {
            if (state.bundle && state.activeViewName === 'tree') {
                scheduleTreeRender();
            }
        });
    }
