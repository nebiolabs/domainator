    // ---- Projection view ----
    //
    // A bundle's views.projection is {algorithm, axes, coords}: coords[a][i] is the bundle's
    // node i on axes[a] (two axes, or three for a 3D projection, shown two at a time). The view
    // keeps the data point at the plot's center and the data units per pixel, the same on both
    // axes, so the projection is never stretched.

    const PROJECTION_MARGIN = {left: 64, right: 16, top: 16, bottom: 46};
    const PROJECTION_POINT_OPACITY = 0.75;
    const PROJECTION_CLICK_SLOP = 3;
    const PROJECTION_HOVER_SLOP = 3;
    // More points than this in view go unlabeled.
    const PROJECTION_MAX_LABELS = 2000;

    function defaultProjectionPointRadius(count) {
        if (count <= 500) {
            return 4;
        }
        if (count <= 5000) {
            return 3;
        }
        return count <= 50000 ? 2 : 1.5;
    }

    function createProjectionState(bundle) {
        return {
            view: bundle.views.projection,
            // The data at the plot's center and the data units per pixel; null until the first
            // paint fits the points to the canvas.
            center: null,
            unitsPerPixel: null,
            box: null,
            pendingFrame: null,
        };
    }

    function projectionPair() {
        const pair = document.getElementById('projection-axes').value.split(',').map(Number);
        return state.projection.view.axes.length === 3 ? pair : [0, 1];
    }

    function projectionPointRadius() {
        const radius = Number(document.getElementById('projection-point-size').value);
        return Number.isFinite(radius) && radius > 0 ? radius : 4;
    }

    function projectionPlotArea(width, height, margin = PROJECTION_MARGIN) {
        return {
            left: margin.left, top: margin.top,
            width: Math.max(10, width - margin.left - margin.right),
            height: Math.max(10, height - margin.top - margin.bottom),
        };
    }

    // Centers every point in the plot area, with equal units on both axes.
    function fitProjectionView(width, height) {
        const projection = state.projection;
        const [xs, ys] = projectionPair().map(axis => projection.view.coords[axis]);
        let [x0, x1, y0, y1] = [Infinity, -Infinity, Infinity, -Infinity];
        xs.forEach((x, i) => {
            x0 = Math.min(x0, x);
            x1 = Math.max(x1, x);
            y0 = Math.min(y0, ys[i]);
            y1 = Math.max(y1, ys[i]);
        });
        const area = projectionPlotArea(width, height);
        if (!(x1 >= x0)) {
            projection.center = {x: 0, y: 0};
            projection.unitsPerPixel = 1;
            return;
        }
        projection.center = {x: (x0 + x1) / 2, y: (y0 + y1) / 2};
        projection.unitsPerPixel = Math.max((x1 - x0) / area.width, (y1 - y0) / area.height) * 1.1 || 1;
    }

    // Functions from data to canvas pixels and back, for the current view.
    function projectionScales(width, height, projection = state.projection, margin = PROJECTION_MARGIN) {
        const area = projectionPlotArea(width, height, margin);
        const u = projection.unitsPerPixel;
        const x0 = projection.center.x - u * area.width / 2;
        const y1 = projection.center.y + u * area.height / 2;
        return {
            area, x0, x1: x0 + u * area.width, y0: y1 - u * area.height, y1,
            sx: x => area.left + (x - x0) / u,
            sy: y => area.top + (y1 - y) / u,
            dataX: px => x0 + (px - area.left) * u,
            dataY: py => y1 - (py - area.top) * u,
        };
    }

    // Calls visit(i, px, py) for each point inside the plot area.
    function forEachProjectionPoint(scales, radius, visit) {
        const [xs, ys] = projectionPair().map(axis => state.projection.view.coords[axis]);
        const {area} = scales;
        for (let i = 0; i < xs.length; i++) {
            const px = scales.sx(xs[i]);
            const py = scales.sy(ys[i]);
            if (px >= area.left - radius && px <= area.left + area.width + radius &&
                    py >= area.top - radius && py <= area.top + area.height + radius) {
                visit(i, px, py);
            }
        }
    }

    function niceProjectionTicks(low, high, count) {
        const span = high - low;
        if (!(span > 0)) {
            return [low];
        }
        const step = niceLength(span / count);
        const ticks = [];
        for (let value = Math.ceil(low / step) * step; value <= high + step * 1e-9; value += step) {
            ticks.push(Math.abs(value) < step * 1e-9 ? 0 : value);
        }
        return ticks;
    }

    function formatProjectionTick(value) {
        return Number(value.toPrecision(6)).toString();
    }

    function drawProjection(ctx, width, height, options = {}) {
        const projection = state.projection;
        const scales = projectionScales(width, height, options.view, options.margin);
        const {area} = scales;
        const radius = projectionPointRadius();
        const coloring = Boolean(colorInfo(currentColorField()));
        const preview = options.preview === false ? null : presetPreviewNodeSet();
        ctx.save();
        ctx.fillStyle = '#ffffff';
        ctx.fillRect(0, 0, width, height);

        // Points, batched by color.
        const byColor = new Map();
        const selected = [];
        const previewed = [];
        const drawn = [];
        forEachProjectionPoint(scales, radius, (i, px, py) => {
            drawn.push([i, px, py]);
            const color = coloring ? state.nodeColorCache[i] : '#1f77b4';
            const shape = shapeForNode(state.nodeStyles[0], i);
            const key = color + '|' + shape;
            if (!byColor.has(key)) {
                byColor.set(key, {color, shape, points: []});
            }
            byColor.get(key).points.push(px, py);
            if (state.selectedNodeIndices.has(i)) {
                selected.push(px, py);
            }
            if (preview && preview.has(i)) {
                previewed.push(px, py);
            }
        });
        ctx.globalAlpha = PROJECTION_POINT_OPACITY;
        byColor.forEach(({points, color, shape}) => {
            ctx.fillStyle = color;
            ctx.beginPath();
            for (let k = 0; k < points.length; k += 2) {
                appendNodeGlyph(ctx, shape, points[k], points[k + 1], radius);
            }
            ctx.fill();
        });
        ctx.globalAlpha = 1;
        [[selected, VIEW_INK, 1.5], [previewed, VIEW_PREVIEW_OUTLINE, 1]].forEach(([points, color, lineWidth]) => {
            if (points.length === 0) { return; }
            ctx.strokeStyle = color;
            ctx.lineWidth = lineWidth;
            ctx.beginPath();
            for (let k = 0; k < points.length; k += 2) {
                ctx.moveTo(points[k] + radius + 2, points[k + 1]);
                ctx.arc(points[k], points[k + 1], radius + 2, 0, 2 * Math.PI);
            }
            ctx.stroke();
        });

        // Labels, when there are few enough points in view to read them.
        const labelField = currentLabelField();
        if (labelField && drawn.length <= PROJECTION_MAX_LABELS) {
            ctx.fillStyle = VIEW_INK;
            ctx.font = '11px ' + VIEW_FONT;
            ctx.textAlign = 'left';
            ctx.textBaseline = 'middle';
            drawn.forEach(([i, px, py]) => {
                const value = labelField === '__node_id__' ? nodeId(i) : metadataValue(i, labelField);
                if (value !== null && value !== undefined && value !== '') {
                    ctx.fillText(String(value), px + radius + 3, py);
                }
            });
        }

        // The frame: margins over any point that strays into them, then axes and ticks.
        ctx.fillStyle = '#ffffff';
        ctx.fillRect(0, 0, width, area.top);
        ctx.fillRect(0, area.top + area.height, width, height - area.top - area.height);
        ctx.fillRect(0, 0, area.left, height);
        ctx.fillRect(area.left + area.width, 0, width - area.left - area.width, height);
        ctx.strokeStyle = VIEW_INK;
        ctx.lineWidth = 1;
        ctx.strokeRect(area.left, area.top, area.width, area.height);
        ctx.fillStyle = VIEW_INK;
        ctx.font = '11px ' + VIEW_FONT;
        ctx.textAlign = 'center';
        ctx.textBaseline = 'top';
        niceProjectionTicks(scales.x0, scales.x1, 8).forEach(value => {
            const x = scales.sx(value);
            ctx.beginPath();
            ctx.moveTo(x, area.top + area.height);
            ctx.lineTo(x, area.top + area.height + 5);
            ctx.stroke();
            ctx.fillText(formatProjectionTick(value), x, area.top + area.height + 7);
        });
        ctx.textAlign = 'right';
        ctx.textBaseline = 'middle';
        niceProjectionTicks(scales.y0, scales.y1, 8).forEach(value => {
            const y = scales.sy(value);
            ctx.beginPath();
            ctx.moveTo(area.left - 5, y);
            ctx.lineTo(area.left, y);
            ctx.stroke();
            ctx.fillText(formatProjectionTick(value), area.left - 7, y);
        });
        const [xAxis, yAxis] = projectionPair().map(axis => projection.view.axes[axis]);
        ctx.textAlign = 'center';
        ctx.textBaseline = 'bottom';
        ctx.fillText(xAxis, area.left + area.width / 2, height - 4);
        // Rotated text needs ctx.rotate, which an SVG export does not record, so the y axis is
        // labeled above the plot instead.
        ctx.textAlign = 'left';
        ctx.textBaseline = 'top';
        ctx.fillText(yAxis, 4, 2);
        ctx.restore();
    }

    function renderProjectionView() {
        const projection = state.projection;
        if (!projection || state.activeViewName !== 'projection') {
            return;
        }
        if (projection.pendingFrame !== null) {
            cancelAnimationFrame(projection.pendingFrame);
            projection.pendingFrame = null;
        }
        const canvas = document.getElementById('projection-canvas');
        const pane = viewCanvasContext(canvas);
        if (!pane) {
            return;
        }
        if (projection.center === null) {
            fitProjectionView(pane.width, pane.height);
        }
        drawProjection(pane.ctx, pane.width, pane.height);
        if (projection.box) {
            const box = projection.box;
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
    }

    function scheduleProjectionRender() {
        const projection = state.projection;
        if (projection && projection.pendingFrame === null) {
            projection.pendingFrame = requestAnimationFrame(() => {
                projection.pendingFrame = null;
                renderProjectionView();
            });
        }
    }

    function projectionCanvasPoint(event) {
        const rect = document.getElementById('projection-canvas').getBoundingClientRect();
        return {x: event.clientX - rect.left, y: event.clientY - rect.top, width: rect.width, height: rect.height};
    }

    function nearestProjectionPoint(point) {
        const scales = projectionScales(point.width, point.height);
        const radius = projectionPointRadius();
        const reach = radius + PROJECTION_HOVER_SLOP;
        let best = null;
        let bestDistance = reach * reach;
        forEachProjectionPoint(scales, radius, (i, px, py) => {
            const distance = (px - point.x) ** 2 + (py - point.y) ** 2;
            if (distance <= bestDistance) {
                best = i;
                bestDistance = distance;
            }
        });
        return best;
    }

    // Zooms by `factor` keeping the data under the canvas point (x, y) in place.
    function zoomProjectionAt(factor, point) {
        const projection = state.projection;
        const scales = projectionScales(point.width, point.height);
        const dataX = scales.dataX(point.x);
        const dataY = scales.dataY(point.y);
        const u = projection.unitsPerPixel / factor;
        projection.center = {
            x: dataX - (dataX - projection.center.x) * u / projection.unitsPerPixel,
            y: dataY - (dataY - projection.center.y) * u / projection.unitsPerPixel,
        };
        projection.unitsPerPixel = u;
        scheduleProjectionRender();
    }

    function projectionTooltipText(index) {
        const projection = state.projection;
        const lines = [nodeId(index)];
        projection.view.axes.forEach((axis, a) => lines.push(axis + ': ' + formatProjectionTick(projection.view.coords[a][index])));
        [currentLabelField(), currentColorField()].forEach(field => {
            if (field && field !== '__node_id__') {
                lines.push(field + ': ' + formatValue(metadataValue(index, field)));
            }
        });
        return lines.join('\n');
    }

    function projectionSelectInBox(box) {
        const scales = projectionScales(box.width, box.height);
        const [x0, x1] = [Math.min(box.x0, box.x1), Math.max(box.x0, box.x1)];
        const [y0, y1] = [Math.min(box.y0, box.y1), Math.max(box.y0, box.y1)];
        let changed = false;
        forEachProjectionPoint(scales, 0, (i, px, py) => {
            if (px < x0 || px > x1 || py < y0 || py > y1) {
                return;
            }
            if (box.deselect) {
                changed = state.selectedNodeIndices.delete(i) || changed;
            } else if (!state.selectedNodeIndices.has(i)) {
                state.selectedNodeIndices.add(i);
                changed = true;
            }
        });
        if (changed) {
            updateMetadataTable();
        }
    }

    // Centers the view on the selected points, zooming out only if they do not fit.
    function projectionFocusSelection() {
        const projection = state.projection;
        if (!projection || state.selectedNodeIndices.size === 0) {
            setStatus('No points are selected.');
            return;
        }
        const [xs, ys] = projectionPair().map(axis => projection.view.coords[axis]);
        let [x0, x1, y0, y1] = [Infinity, -Infinity, Infinity, -Infinity];
        state.selectedNodeIndices.forEach(i => {
            x0 = Math.min(x0, xs[i]);
            x1 = Math.max(x1, xs[i]);
            y0 = Math.min(y0, ys[i]);
            y1 = Math.max(y1, ys[i]);
        });
        const rect = document.getElementById('projection-canvas').getBoundingClientRect();
        const area = projectionPlotArea(rect.width, rect.height);
        projection.center = {x: (x0 + x1) / 2, y: (y0 + y1) / 2};
        projection.unitsPerPixel = Math.max(projection.unitsPerPixel, (x1 - x0) / area.width * 1.1, (y1 - y0) / area.height * 1.1);
        renderProjectionView();
    }

    function projectionResetView() {
        if (!state.projection) { return; }
        state.projection.center = null;
        renderProjectionView();
    }

    function projectionExportSize() {
        const rect = document.getElementById('projection-canvas').getBoundingClientRect();
        const size = {width: rect.width, height: rect.height, options: {}};
        if (entireContentExport()) {
            const projection = state.projection;
            const [xs, ys] = projectionPair().map(axis => projection.view.coords[axis]);
            let [x0, x1, y0, y1] = [Infinity, -Infinity, Infinity, -Infinity];
            xs.forEach((x, i) => {
                x0 = Math.min(x0, x); x1 = Math.max(x1, x);
                y0 = Math.min(y0, ys[i]); y1 = Math.max(y1, ys[i]);
            });
            if (!xs.length) { x0 = x1 = y0 = y1 = 0; }
            const measurer = document.createElement('canvas').getContext('2d');
            if (!measurer) { throw new Error('Could not measure projection labels.'); }
            measurer.font = '11px ' + VIEW_FONT;
            let padX = Math.max(16, projectionPointRadius() + 4);
            const padY = Math.max(16, projectionPointRadius() + 4);
            const field = currentLabelField();
            if (field && xs.length <= PROJECTION_MAX_LABELS) {
                for (let i = 0; i < xs.length; i++) {
                    const value = field === '__node_id__' ? nodeId(i) : metadataValue(i, field);
                    padX = Math.max(padX, measurer.measureText(String(value ?? '')).width + projectionPointRadius() + 4);
                }
            }
            const u = projection.unitsPerPixel;
            const margin = {...PROJECTION_MARGIN};
            for (const tick of niceProjectionTicks(y0 - padY * u, y1 + padY * u, 8)) {
                margin.left = Math.max(margin.left, Math.ceil(measurer.measureText(formatProjectionTick(tick)).width + 12));
            }
            const axisWidth = Math.max(...projectionPair().map(axis => measurer.measureText(projection.view.axes[axis]).width));
            size.width = Math.ceil(Math.max((x1 - x0) / u + padX * 2 + margin.left + margin.right, axisWidth + 8));
            size.height = Math.ceil((y1 - y0) / u + padY * 2 + margin.top + margin.bottom);
            size.options = {view: {center: {x: x0 / 2 + x1 / 2, y: y0 / 2 + y1 / 2}, unitsPerPixel: u},
                margin, preview: false};
        }
        validateImageExportSize(Math.round(size.width), Math.round(size.height));
        return size;
    }

    function exportProjectionPNG() {
        if (!state.projection) { return; }
        const size = projectionExportSize();
        return exportCanvasPNG(size, ctx => drawProjection(ctx, size.width, size.height, size.options));
    }

    function buildProjectionSVG() {
        const size = projectionExportSize();
        const ctx = createSvgContext(Math.round(size.width), Math.round(size.height));
        drawProjection(ctx, size.width, size.height, size.options);
        return ctx.toSVG();
    }

    function exportProjectionSVG() {
        if (!state.projection) { return; }
        triggerDownload(URL.createObjectURL(new Blob([buildProjectionSVG()], {type: 'image/svg+xml'})), exportBaseName() + imageExportSuffix('svg'), true);
    }

    function setupProjectionViewControls() {
        const canvas = document.getElementById('projection-canvas');
        let press = null;
        canvas.addEventListener('pointerdown', event => {
            if (event.button !== 0 || !state.projection) {
                return;
            }
            press = {
                id: event.pointerId, point: projectionCanvasPoint(event), moved: false,
                box: event.shiftKey, deselect: event.altKey, center: state.projection.center,
            };
        });
        canvas.addEventListener('pointermove', event => {
            if (!state.projection) {
                return;
            }
            const point = projectionCanvasPoint(event);
            if (press && event.pointerId === press.id) {
                const dx = point.x - press.point.x;
                const dy = point.y - press.point.y;
                if (!press.moved && Math.hypot(dx, dy) >= PROJECTION_CLICK_SLOP) {
                    press.moved = true;
                    canvas.setPointerCapture(event.pointerId);
                    if (!press.box) {
                        canvas.classList.add('panning');
                    }
                    hideViewTip();
                }
                if (press.moved) {
                    if (press.box) {
                        state.projection.box = {x0: press.point.x, y0: press.point.y, x1: point.x, y1: point.y};
                    } else {
                        const u = state.projection.unitsPerPixel;
                        state.projection.center = {x: press.center.x - dx * u, y: press.center.y + dy * u};
                    }
                    scheduleProjectionRender();
                    return;
                }
            }
            const index = nearestProjectionPoint(point);
            if (index === null) {
                hideViewTip();
            } else {
                showViewTip(event, projectionTooltipText(index));
            }
        });
        const release = event => {
            if (!press || event.pointerId !== press.id) {
                return;
            }
            const point = projectionCanvasPoint(event);
            if (press.moved && press.box) {
                projectionSelectInBox({
                    x0: press.point.x, y0: press.point.y, x1: point.x, y1: point.y,
                    width: point.width, height: point.height, deselect: press.deselect,
                });
                state.projection.box = null;
                renderProjectionView();
            } else if (!press.moved && event.type === 'pointerup') {
                const index = nearestProjectionPoint(point);
                if (index !== null) {
                    if (!state.selectedNodeIndices.delete(index)) {
                        state.selectedNodeIndices.add(index);
                    }
                    renderProjectionView();
                    updateMetadataTable();
                }
            }
            press = null;
            canvas.classList.remove('panning');
        };
        canvas.addEventListener('pointerup', release);
        canvas.addEventListener('pointercancel', release);
        canvas.addEventListener('pointerleave', hideViewTip);
        canvas.addEventListener('wheel', event => {
            if (!state.projection) {
                return;
            }
            event.preventDefault();
            const delta = -event.deltaY * (event.deltaMode === 1 ? 0.05 : event.deltaMode ? 1 : 0.002);
            zoomProjectionAt(2 ** delta, projectionCanvasPoint(event));
            hideViewTip();
        }, {passive: false});
        canvas.addEventListener('dblclick', event => {
            if (state.projection) {
                zoomProjectionAt(event.shiftKey ? 0.5 : 2, projectionCanvasPoint(event));
            }
        });
        document.getElementById('projection-point-size').addEventListener('input', () => renderProjectionView());
        document.getElementById('projection-axes').addEventListener('change', () => {
            if (!state.projection) { return; }
            state.projection.center = null;
            renderProjectionView();
        });
        window.addEventListener('resize', () => {
            if (state.bundle && state.activeViewName === 'projection') {
                scheduleProjectionRender();
            }
        });
    }
