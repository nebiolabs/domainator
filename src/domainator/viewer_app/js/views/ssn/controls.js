    // The listeners of the SSN view's own controls: its canvas, the threshold slider and the
    // split chart, and the layout and drawing options.
    function setupSsnViewControls() {
        clusterCanvas.addEventListener('click', event => {
            if (!state.bundle) {
                return;
            }
            if (state.suppressClick) {
                state.suppressClick = false;
                return;
            }
            const point = canvasCoordinatesFromEvent(event);
            if (event.ctrlKey || event.metaKey) {
                const nodeIndex = hitTestNodeAt(point.x, point.y);
                if (nodeIndex !== null) {
                    toggleSelectionForNode(nodeIndex);
                }
                return;
            }
            const hit = hitTestComponentAt(point.x, point.y);
            if (hit) {
                toggleSelectionForComponent(hit.componentId);
            }
        });

        clusterCanvas.addEventListener('wheel', event => {
            if (!state.bundle) {
                return;
            }
            event.preventDefault();
            const point = canvasCoordinatesFromEvent(event);
            const anchor = screenToWorldPoint(point.x, point.y);
            const zoomFactor = Math.exp(-event.deltaY * 0.0012);
            const nextScale = Math.max(state.viewTransform.minScale, Math.min(state.viewTransform.maxScale, state.viewTransform.scale * zoomFactor));
            state.viewTransform.scale = nextScale;
            state.viewTransform.offsetX = point.x - (anchor.x * nextScale);
            state.viewTransform.offsetY = point.y - (anchor.y * nextScale);
            scheduleClusterRender();
        }, {passive: false});

        clusterCanvas.addEventListener('mousedown', event => {
            if (!state.bundle) {
                return;
            }
            const point = canvasCoordinatesFromEvent(event);
            if (event.shiftKey) {
                state.dragState = {mode: 'select', nodeMode: event.ctrlKey || event.metaKey, deselect: event.altKey, startX: point.x, startY: point.y, endX: point.x, endY: point.y, moved: false};
                state.selectionBox = normalizedSelectionBox(state.dragState);
            } else {
                state.dragState = {mode: 'pan', startX: point.x, startY: point.y, originOffsetX: state.viewTransform.offsetX, originOffsetY: state.viewTransform.offsetY, moved: false};
            }
            clusterCanvas.style.cursor = state.dragState.mode === 'select' ? 'crosshair' : 'grabbing';
        });

        window.addEventListener('mousemove', event => {
            if (!state.dragState) {
                return;
            }
            const point = canvasCoordinatesFromEvent(event);
            if (state.dragState.mode === 'pan') {
                const dx = point.x - state.dragState.startX;
                const dy = point.y - state.dragState.startY;
                state.dragState.moved = state.dragState.moved || Math.abs(dx) > 3 || Math.abs(dy) > 3;
                state.viewTransform.offsetX = state.dragState.originOffsetX + dx;
                state.viewTransform.offsetY = state.dragState.originOffsetY + dy;
                scheduleClusterRender();
                return;
            }
            state.dragState.endX = point.x;
            state.dragState.endY = point.y;
            state.dragState.moved = state.dragState.moved || Math.abs(point.x - state.dragState.startX) > 3 || Math.abs(point.y - state.dragState.startY) > 3;
            state.selectionBox = normalizedSelectionBox(state.dragState);
            scheduleClusterRender();
        });

        window.addEventListener('mouseup', () => {
            if (!state.dragState) {
                return;
            }
            state.suppressClick = Boolean(state.dragState.moved);
            if (state.dragState.mode === 'select' && state.selectionBox && (state.selectionBox.width > 4 || state.selectionBox.height > 4)) {
                if (state.dragState.nodeMode) {
                    selectNodesInBox(state.selectionBox, state.dragState.deselect);
                } else {
                    selectComponentsInBox(state.selectionBox, state.dragState.deselect);
                }
            }
            state.dragState = null;
            state.selectionBox = null;
            clusterCanvas.style.cursor = 'grab';
            scheduleClusterRender();
        });

        clusterCanvas.addEventListener('mousemove', event => {
            if (state.dragState) {
                return;
            }
            clusterCanvas.style.cursor = event.shiftKey ? 'crosshair' : 'grab';
        });

        clusterCanvas.addEventListener('mouseleave', () => {
            if (!state.dragState) {
                clusterCanvas.style.cursor = 'grab';
            }
        });

        document.getElementById('threshold-slider').addEventListener('input', () => {
            scheduleThresholdUI(false);
        });
        document.getElementById('threshold-slider').addEventListener('change', () => {
            const stop = currentSliderStop();
            snapSliderToStop(stop);
            // Keep the user's current pan/zoom when the slider is released; only snap to the
            // nearest stop. (Use the Reset view button to re-fit.)
            scheduleThresholdUI(false);
        });
        document.getElementById('split-event-cap').addEventListener('change', event => {
            // Re-deriving the whole series is one pass over the MST edges, so this is a
            // plain 'change' handler rather than anything debounced or deferred.
            applyMaxMergeEvents(event.target.value);
        });
        document.getElementById('threshold-input').addEventListener('change', jumpToThresholdValue);
        document.getElementById('threshold-input').addEventListener('keydown', event => {
            // Up/Down step stops -- the behavior a number input's spinner would have had
            // if its step were the distance to the next split. Left/Right are left alone
            // so they still move the caret inside the field.
            if (event.key === 'ArrowUp' || event.key === 'ArrowDown') {
                event.preventDefault();
                stepThreshold(event.key === 'ArrowUp' ? 1 : -1);
                return;
            }
            if (event.key !== 'Enter') {
                return;
            }
            event.preventDefault();
            jumpToThresholdValue();
        });
        document.getElementById('threshold-step-down').addEventListener('click', () => stepThreshold(-1));
        document.getElementById('threshold-step-up').addEventListener('click', () => stepThreshold(1));
        document.getElementById('min-cluster-size').addEventListener('input', () => {
            scheduleThresholdUI(true);
        });
        document.getElementById('layout-algorithm').addEventListener('change', () => {
            updateComponentSortButton();
            if (!state.bundle) {
                return;
            }
            drawClusterView(true);
        });
        document.getElementById('render-cluster-bounds').addEventListener('change', () => renderClusterView());
        document.getElementById('render-nodes').addEventListener('change', () => renderClusterView());
        document.getElementById('leaf-pruning-only').addEventListener('change', () => {
            updateCollapseLongPathsControl();
            scheduleThresholdUI(true);
        });
        document.getElementById('collapse-long-paths').addEventListener('change', () => {
            scheduleThresholdUI(true);
        });
        document.getElementById('show-node-counts').addEventListener('change', () => renderClusterView());
        document.getElementById('show-edge-scores').addEventListener('change', () => renderClusterView());
        document.getElementById('reduce-elongation').addEventListener('change', () => {
            if (!state.bundle) {
                return;
            }
            renderClusterView();
        });
        document.getElementById('sort-components-by-size').addEventListener('click', () => {
            const button = document.getElementById('sort-components-by-size');
            const enabled = button.getAttribute('aria-pressed') === 'true';
            button.setAttribute('aria-pressed', enabled ? 'false' : 'true');
            updateComponentSortButton();
            if (!state.bundle) {
                return;
            }
            drawClusterView(true);
        });
        document.getElementById('export-split-png').addEventListener('click', exportSplitChartPNG);
        document.getElementById('export-split-svg').addEventListener('click', exportSplitChartSVG);
        document.getElementById('save-extraction').addEventListener('click', saveExtractionFile);
        document.getElementById('focus-largest-cluster').addEventListener('click', () => {
            if (!state.bundle || state.visibleLayout.length === 0) { return; }
            const hierarchyNodes = state.bundle.views.ssn.hierarchy.nodes;
            let largestItem = state.visibleLayout[0];
            let largestSize = hierarchyNodes[largestItem.componentId].size;
            for (const item of state.visibleLayout) {
                const size = hierarchyNodes[item.componentId].size;
                if (size > largestSize) { largestSize = size; largestItem = item; }
            }
            state.viewTransform.offsetX = clusterCanvas.width / 2 - largestItem.x * state.viewTransform.scale;
            state.viewTransform.offsetY = clusterCanvas.height / 2 - largestItem.y * state.viewTransform.scale;
            scheduleClusterRender();
        });
        window.addEventListener('resize', () => {
            if (!state.bundle || state.activeViewName !== 'ssn') {
                return;
            }
            drawClusterView(true);
        });
    }
