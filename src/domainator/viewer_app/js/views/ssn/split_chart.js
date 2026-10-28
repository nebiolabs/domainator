    // Axis titles for the split chart. A `product` merge impact is a product of two
    // component sizes, not a count of nodes, so it must not be labeled as one.
    function splitAxisLabels() {
        const labels = VIEWER_CONFIG.splitAxisLabels;
        const metric = state.bundle && state.bundle.views.ssn
            ? state.bundle.views.ssn.merge_impact_metric
            : null;
        return labels[metric] || labels[VIEWER_CONFIG.defaultMergeImpactMetric];
    }

    // ---- Split chart ----
    //
    // One model, three consumers. splitChartLayout() computes the axis box, the scales,
    // every tick value and every mark's position; drawSplitChart() paints it to a canvas,
    // buildSplitChartSVG() emits the same thing as SVG, and the hover handlers hit-test
    // against the positions the paint pass recorded. buildClusterViewSVG's header warns
    // about what happens when two renderers each compute their own geometry -- this is
    // that warning taken seriously, so the exported chart cannot drift from the drawn one.

    const SPLIT_CHART_COLORS = {
        axis: 'rgba(92,106,112,0.35)',
        tick: 'rgba(92,106,112,0.4)',
        text: '#5c6a70',
        stem: '#dd9687',
        bead: '#c8553d',
        movingSum: '#2f6f8f',
        marker: '#1e2a2f',
        markerDot: '#e29b4b',
    };
    const SPLIT_CHART_TITLE_FONT = 12;
    const SPLIT_CHART_TICK_FONT = 11;
    const SPLIT_CHART_MESSAGE_FONT = 18;
    const SPLIT_CHART_TICK_LENGTH = 6;
    const SPLIT_CHART_BEAD_DRAW_RADIUS = 4;

    // Decimals needed to print `step` exactly, so a tick label names the value the tick
    // sits at. Without this the axis lies by rounding: a tick at 0.8736 printed as "0.87"
    // reads as misplaced next to a lollipop whose own readout also says 0.87.
    function decimalsForTickStep(step) {
        const magnitude = Math.abs(step);
        if (!Number.isFinite(magnitude) || magnitude === 0) { return 2; }
        for (let decimals = 0; decimals <= 8; decimals++) {
            if (Math.abs(Number(magnitude.toFixed(decimals)) - magnitude) <= magnitude * 1e-9) {
                return decimals;
            }
        }
        return 8;
    }

    // Tick values on round numbers, with the precision their step requires. Ticks used to
    // sit at fixed fractions of the domain, which put them at arbitrary values that
    // formatValue then rounded to two decimals -- both halves of the misalignment.
    // `options.integer` marks an axis of counts, where a tick at 2.5 nodes names nothing.
    function splitAxisTicks(min, max, targetCount, options = {}) {
        let values = chartAxisTicks(min, max, targetCount)
            .filter(value => value >= min - 1e-12 && value <= max + 1e-12);
        if (options.integer && (max - min) < targetCount) {
            // Too narrow for a nice fractional step to mean anything: one tick per integer.
            values = [];
            for (let value = Math.ceil(min); value <= Math.floor(max); value++) {
                values.push(value);
            }
        }
        if (values.length === 0) { values = [min]; }
        const step = values.length > 1 ? values[1] - values[0] : (Math.abs(values[0]) || 1);
        const decimals = decimalsForTickStep(step);
        return values.map(value => ({
            value,
            // toLocaleString rather than toFixed so a 50,000-node axis keeps its
            // separators, the way formatValue renders integers elsewhere.
            label: value.toLocaleString(undefined, {
                minimumFractionDigits: decimals,
                maximumFractionDigits: decimals,
            }),
        }));
    }

    function splitChartLayout(viewWidth, viewHeight) {
        const margin = {top: 26, right: 80, bottom: 62, left: 76};
        const width = viewWidth - margin.left - margin.right;
        const height = viewHeight - margin.top - margin.bottom;
        const axisLabels = splitAxisLabels();
        const frame = {
            viewWidth,
            viewHeight,
            margin,
            width,
            height,
            titles: {x: 'Threshold', y: axisLabels.largest, y2: axisLabels.movingSum},
        };

        if (!state.bundle) {
            return Object.assign(frame, {
                message: 'Load a bundle to render split events.',
                axes: false,
            });
        }

        const events = state.series ? state.series.selectedRows : [];
        if (events.length === 0) {
            return Object.assign(frame, {
                message: 'No split events in this bundle.',
                axes: true,
            });
        }

        // Stems show the LARGEST single merge at each threshold, not the tie group's sum:
        // one big cluster splitting off and a swarm of singletons are opposite stories that
        // the sum renders identically.
        const maxImpact = Math.max(...events.map(event => event.largest_merge), 1);
        const movingSum = state.series.movingSum || {x: [], y: []};
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
        ).map(event => ({
            event,
            x: xFor(event.threshold_value),
            stemY: yFor(event.largest_merge),
            // JSON object keys are always strings, so this Number() is required.
            beads: Object.keys(event.merge_size_counts || {}).map(size => ({
                size: Number(size),
                count: event.merge_size_counts[size],
                y: yFor(Number(size)),
            })),
        }));

        // The trace is a step function -- a moving sum over discrete events genuinely is
        // one, and smoothing misrepresents where the events sit -- so the vertices are
        // built once here and both painters stroke the same polyline.
        const movingSumPoints = [];
        for (let index = 0; index < movingSumX.length; index++) {
            const x = xFor(movingSumX[index]);
            if (index > 0) {
                movingSumPoints.push({x, y: y2For(movingSumY[index - 1])});
            }
            movingSumPoints.push({x, y: y2For(movingSumY[index])});
        }

        const stop = currentSliderStop();
        // Clamped into the plotted range first. The infinity stop sits above every score
        // in the series and the floor stop below every one, so both belong at an end of
        // the axis rather than off it -- the floor stop's marker used to be drawn outside
        // the plot box, to the left of the y-axis. Clamping first also means a zoomed
        // window drops them like any other out-of-window threshold.
        const markerValue = stop === null ? null : Math.min(Math.max(
            stop.threshold_value === null ? dataMax : stop.threshold_value, dataMin), dataMax);

        return Object.assign(frame, {
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
                .map(tick => Object.assign({x: xFor(tick.value)}, tick)),
            // Both vertical axes count nodes (or, under --merge_impact_metric product, a
            // product of counts); either way a fractional tick is not a quantity.
            yTicks: splitAxisTicks(0, maxImpact, 6, {integer: true})
                .map(tick => Object.assign({y: yFor(tick.value)}, tick)),
            y2Ticks: splitAxisTicks(0, maxMovingSum, 6, {integer: true})
                .map(tick => Object.assign({y: y2For(tick.value)}, tick)),
            marks,
            movingSum: {window: movingSum.window || 0, x: movingSumX, y: movingSumY, points: movingSumPoints},
            // Dropped once the window no longer contains it: a marker pinned to the
            // edge of a zoomed axis would claim the threshold is sitting there.
            markerX: markerValue === null || markerValue < minThreshold || markerValue > maxThreshold
                ? null
                : xFor(markerValue),
        });
    }

    function drawSplitChart(ctx = splitContext, viewWidth = splitCanvas.width,
                            viewHeight = splitCanvas.height, options = {}) {
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

        if (layout.axes) {
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
        }

        if (layout.message) {
            if (onScreen) {
                recordSplitChartGeometry(null);
                hideSplitChartTip();
            }
            context.fillStyle = SPLIT_CHART_COLORS.text;
            context.font = SPLIT_CHART_MESSAGE_FONT + 'px Georgia';
            // Set explicitly: the context is reused across calls, so inheriting whatever
            // alignment the last paint left behind once centered this message on x=30 and
            // ran half of it off the left edge.
            context.textAlign = 'left';
            context.textBaseline = 'alphabetic';
            context.fillText(layout.message, 30, 50);
            return;
        }

        context.fillStyle = SPLIT_CHART_COLORS.text;
        context.font = SPLIT_CHART_TICK_FONT + 'px Georgia';
        context.textBaseline = 'top';
        context.textAlign = 'center';
        context.strokeStyle = SPLIT_CHART_COLORS.tick;
        context.lineWidth = 1;
        layout.xTicks.forEach(tick => {
            context.beginPath();
            context.moveTo(tick.x, margin.top + height);
            context.lineTo(tick.x, margin.top + height + SPLIT_CHART_TICK_LENGTH);
            context.stroke();
            context.fillText(tick.label, tick.x, margin.top + height + SPLIT_CHART_TICK_LENGTH + 4);
        });

        context.textAlign = 'right';
        context.textBaseline = 'middle';
        layout.yTicks.forEach(tick => {
            context.beginPath();
            context.moveTo(margin.left - SPLIT_CHART_TICK_LENGTH, tick.y);
            context.lineTo(margin.left, tick.y);
            context.stroke();
            context.fillText(tick.label, margin.left - SPLIT_CHART_TICK_LENGTH - 6, tick.y);
        });

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
        layout.y2Ticks.forEach(tick => {
            context.beginPath();
            context.moveTo(rightAxisX, tick.y);
            context.lineTo(rightAxisX + SPLIT_CHART_TICK_LENGTH, tick.y);
            context.stroke();
            context.fillText(tick.label, rightAxisX + SPLIT_CHART_TICK_LENGTH + 4, tick.y);
        });

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
        layout.marks.forEach(mark => {
            context.strokeStyle = SPLIT_CHART_COLORS.stem;
            context.lineWidth = 1.5;
            context.beginPath();
            context.moveTo(mark.x, margin.top + height);
            context.lineTo(mark.x, mark.stemY);
            context.stroke();
            context.fillStyle = SPLIT_CHART_COLORS.bead;
            mark.beads.forEach(bead => {
                context.beginPath();
                context.arc(mark.x, bead.y, SPLIT_CHART_BEAD_DRAW_RADIUS, 0, Math.PI * 2);
                context.fill();
            });
        });

        if (layout.movingSum.points.length > 0) {
            context.strokeStyle = SPLIT_CHART_COLORS.movingSum;
            context.lineWidth = 1.75;
            context.beginPath();
            layout.movingSum.points.forEach((point, index) => {
                if (index === 0) {
                    context.moveTo(point.x, point.y);
                } else {
                    context.lineTo(point.x, point.y);
                }
            });
            context.stroke();
        }
        context.restore();

        // The dashed line and its dot mark where the slider is sitting right now. That is
        // UI state, not something the chart measures, so neither export draws it --
        // buildSplitChartSVG omits it for the same reason.
        if (onScreen && layout.markerX !== null) {
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
        }

        // Recorded from the layout the paint pass just used, so a mark the hover finds is a
        // mark that was actually drawn. See setupSplitChartHover / splitChartHitAt.
        if (onScreen) {
            recordSplitChartGeometry({
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
                events: layout.marks.map(mark => ({event: mark.event, x: mark.x, beads: mark.beads})),
                movingSumWindow: layout.movingSum.window,
                movingSumX: layout.movingSum.x,
                movingSumY: layout.movingSum.y,
            });
        }
    }

    function buildSplitChartSVG() {
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

        if (layout.axes) {
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
        }

        if (layout.message) {
            parts.push('<text x="30" y="50" font-size="' + SPLIT_CHART_MESSAGE_FONT + '" fill="' +
                SPLIT_CHART_COLORS.text + '">' + escapeXml(layout.message) + '</text>');
            parts.push('</svg>');
            return parts.join('\n');
        }

        // The same clip the canvas painter applies while drawing the data marks.
        parts.push('<defs><clipPath id="split-plot-clip"><rect x="' + margin.left + '" y="' + margin.top +
            '" width="' + width + '" height="' + height + '"/></clipPath></defs>');

        // Ticks. Values, positions and label text all come from the layout, so the export
        // cannot round or place them differently from the canvas.
        const tickColor = svgColorParts(SPLIT_CHART_COLORS.tick);
        const tickMarks = [];
        const tickLabels = [];
        layout.xTicks.forEach(tick => {
            tickMarks.push('<path d="M' + fmt(tick.x) + ' ' + (margin.top + height) + ' V' +
                (margin.top + height + SPLIT_CHART_TICK_LENGTH) + '"/>');
            tickLabels.push('<text x="' + fmt(tick.x) + '" y="' +
                (margin.top + height + SPLIT_CHART_TICK_LENGTH + 4) + '" text-anchor="middle"' +
                ' dominant-baseline="hanging">' + escapeXml(tick.label) + '</text>');
        });
        layout.yTicks.forEach(tick => {
            tickMarks.push('<path d="M' + (margin.left - SPLIT_CHART_TICK_LENGTH) + ' ' + fmt(tick.y) +
                ' H' + margin.left + '"/>');
            tickLabels.push('<text x="' + (margin.left - SPLIT_CHART_TICK_LENGTH - 6) + '" y="' + fmt(tick.y) +
                '" text-anchor="end" dominant-baseline="central">' + escapeXml(tick.label) + '</text>');
        });
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
        layout.y2Ticks.forEach(tick => {
            rightTickMarks.push('<path d="M' + rightAxisX + ' ' + fmt(tick.y) + ' H' +
                (rightAxisX + SPLIT_CHART_TICK_LENGTH) + '"/>');
            rightTickLabels.push('<text x="' + (rightAxisX + SPLIT_CHART_TICK_LENGTH + 4) + '" y="' +
                fmt(tick.y) + '" dominant-baseline="central">' + escapeXml(tick.label) + '</text>');
        });
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
        layout.marks.forEach(mark => {
            mark.beads.forEach(bead => {
                beads.push('<circle cx="' + fmt(mark.x) + '" cy="' + fmt(bead.y) + '" r="' +
                    SPLIT_CHART_BEAD_DRAW_RADIUS + '"/>');
            });
        });
        parts.push('<g clip-path="url(#split-plot-clip)" fill="' + SPLIT_CHART_COLORS.bead + '">' +
            beads.join('') + '</g>');

        if (layout.movingSum.points.length > 0) {
            const d = layout.movingSum.points
                .map((point, index) => (index === 0 ? 'M' : ' L') + fmt(point.x) + ' ' + fmt(point.y))
                .join('');
            parts.push('<path d="' + d + '" clip-path="url(#split-plot-clip)" fill="none" stroke="' +
                SPLIT_CHART_COLORS.movingSum + '" stroke-width="1.75"/>');
        }

        // No threshold marker: this builder is only ever an export, and the dashed line
        // and orange dot say where the slider is, which is not a property of the data.
        // See the matching `onScreen` guard in drawSplitChart.

        parts.push('</svg>');
        return parts.join('\n');
    }

    function exportSplitChartSVG() {
        if (!state.bundle) {
            return;
        }
        const blob = new Blob([buildSplitChartSVG()], {type: 'image/svg+xml'});
        triggerDownload(URL.createObjectURL(blob), exportBaseName() + '_split_events.svg', true);
    }

    function exportSplitChartPNG() {
        if (!state.bundle) {
            return;
        }
        // Re-rasterized at the chosen density rather than snapshotted and upscaled, and
        // offscreen so the hover geometry the on-screen chart recorded is left alone.
        const scaleFactor = selectedPngScale();
        const target = document.createElement('canvas');
        target.width = Math.round(splitCanvas.width * scaleFactor);
        target.height = Math.round(splitCanvas.height * scaleFactor);
        const targetContext = target.getContext('2d');
        if (!targetContext) {
            window.alert('Could not export PNG at ' + scaleFactor + '×; try a lower resolution.');
            return;
        }
        targetContext.scale(scaleFactor, scaleFactor);
        drawSplitChart(targetContext, splitCanvas.width, splitCanvas.height, {onScreen: false});
        target.toBlob(blob => {
            if (!blob) {
                window.alert('The chart is too large to export as a ' + scaleFactor + '× PNG (' +
                    target.width + '×' + target.height + ' px). Try a lower resolution.');
                return;
            }
            const suffix = scaleFactor > 1 ? '_split_events@' + scaleFactor + 'x.png' : '_split_events.png';
            triggerDownload(URL.createObjectURL(blob), exportBaseName() + suffix, true);
        }, 'image/png');
    }

