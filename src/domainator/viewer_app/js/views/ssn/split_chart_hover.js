/*
 * Hover readout and click-to-jump for the split chart.
 *
 * Included into viewer.js's script scope; see session_state.js.
 * The readout carries the same quantities matrix_report's Plotly hovertemplates
 * do, so the canvas chart and the Plotly one tell the same story.
 */

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
