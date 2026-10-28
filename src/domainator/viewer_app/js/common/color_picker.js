    // ---- Custom color palette: load/save color table TSV (categorical only) ----

    // Guard shared by the TSV load/save paths: a 2-column value->color table only
    // describes discrete mappings, so it applies to categorical columns. Returns the
    // column name, or null after setting an explanatory status message.
    function categoricalColumnForColorTable(action) {
        const columnName = currentColorField();
        const info = colorInfo(columnName);
        if (!columnName || !info) {
            setStatus(action + ' color table: choose a "Color by" column first.');
            return null;
        }
        if (info.type !== 'categorical') {
            setStatus(action + ' color table: only categorical columns use a value/color table (use the picker for "' + columnName + '").');
            return null;
        }
        return columnName;
    }

    function loadColorTableText(text) {
        const columnName = categoricalColumnForColorTable('Load');
        if (!columnName) {
            return;
        }
        const colors = {};
        let nullColor = null;
        let loaded = 0;
        let skipped = 0;
        String(text).split('\n').forEach(rawLine => {
            const line = rawLine.trim();
            if (line === '') {
                return;
            }
            const parts = line.split('\t');
            const hex = parts.length === 2 ? normalizeColorHex(parts[1]) : null;
            if (parts.length !== 2 || !hex) {
                skipped += 1;
                return;
            }
            if (parts[0] === '\u2014') {
                nullColor = hex;
            } else {
                colors[parts[0]] = hex;
            }
            loaded += 1;
        });
        state.customPalettes[paletteKey(columnName)] = {type: 'categorical', colors, nullColor};
        rebuildNodeColorCache();
        repaintNow();
        if (colorPickerIsOpen()) {
            openColorPicker();
        }
        setStatus('Loaded ' + loaded + ' colors for "' + columnName + '"' +
            (skipped > 0 ? ', skipped ' + skipped + ' malformed line(s).' : '.'));
    }

    function saveColorTableTSV() {
        const columnName = categoricalColumnForColorTable('Save');
        if (!columnName) {
            return;
        }
        const palette = customPalette(columnName);
        const distinct = distinctColumnValues(columnName, Infinity);
        // Write the effective on-screen color for every distinct value (custom override
        // or the current default), so the file matches the view and round-trips on load.
        const rows = distinct.values.map(entry => ({
            key: entry.key,
            color: (palette && palette.colors && palette.colors[entry.key]) || cssToHex(categoricalColor(entry.raw, palette)),
        }));
        rows.sort((a, b) => a.key.localeCompare(b.key));
        if (distinct.nullCount > 0) {
            // Em-dash row carries the no-value color; written last (null sorts last).
            rows.push({key: '\u2014', color: (palette && palette.nullColor) || cssToHex('#b3a89d')});
        }
        const tsv = rows.map(row => row.key + '\t' + row.color).join('\n') + '\n';
        const blob = new Blob([tsv], {type: 'text/tab-separated-values'});
        const href = URL.createObjectURL(blob);
        triggerDownload(href, exportBaseName() + '_' + columnName.replace(/\s+/g, '_') + '_colors.tsv', true);
    }

    // ---- Custom color palette: legend export (SVG + PNG) ----

    // Build a standalone legend SVG for the current color column. Returns
    // {svg, width, height}. Discrete columns get a box+label per value (ported from
    // domainator color_table_to_legend); numeric columns get a sampled gradient bar with
    // min/mid/max ticks. All fills pass through svgColorParts so default hsl() colors
    // become SVG-editor-safe rgb()+opacity. No external refs (keeps PNG untainted).
    function buildLegendSVG() {
        const columnName = currentColorField();
        const info = colorInfo(columnName);
        const title = columnName || 'Legend';
        const titleFontSize = 24;
        const itemFontSize = 20;
        const strokeWidth = 2;
        const padding = 10;
        const titleSpace = 40;
        const parts = [];
        const open = (width, height) => {
            parts.push('<?xml version="1.0" encoding="UTF-8" standalone="no"?>');
            parts.push('<svg xmlns="http://www.w3.org/2000/svg" width="' + width + '" height="' + height +
                '" viewBox="0 0 ' + width + ' ' + height + '">');
            parts.push('<g font-family="Arial, sans-serif">');
            parts.push('<rect x="0" y="0" width="' + width + '" height="' + height +
                '" fill="white" stroke="black" stroke-width="' + strokeWidth + '"/>');
            parts.push('<text x="' + padding + '" y="' + titleFontSize + '" font-size="' + titleFontSize +
                '" font-weight="bold">' + escapeXml(title) + '</text>');
        };
        const swatchRow = (color, label, y, boxHeight, boxWidth) => {
            const fill = svgColorParts(color);
            parts.push('<rect x="' + padding + '" y="' + y + '" width="' + boxWidth + '" height="' + (boxHeight - 10) +
                '" fill="' + fill.color + '" fill-opacity="' + fill.opacity +
                '" stroke="black" stroke-width="' + strokeWidth + '"/>');
            parts.push('<text x="' + (boxWidth + padding * 2) + '" y="' + (y + boxHeight / 2) +
                '" font-size="' + itemFontSize + '" dominant-baseline="central">' + escapeXml(label) + '</text>');
        };

        if (!info || info.type !== 'numeric') {
            // Discrete legend.
            const palette = (info && customPalette(columnName)) || null;
            // Enumerate every distinct value (no cap) so the saved legend is complete.
            const distinct = distinctColumnValues(columnName, Infinity);
            const labels = distinct.values.map(entry => ({color: categoricalColor(entry.raw, palette), label: entry.label}));
            if (distinct.nullCount > 0) {
                labels.push({color: (palette && palette.nullColor) || '#b3a89d', label: '\u2014'});
            }
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
            labels.forEach(item => {
                swatchRow(item.color, item.label, y, boxHeight, boxWidth);
                y += boxHeight;
            });
            parts.push('</g>');
            parts.push('</svg>');
            return {svg: parts.join('\n'), width, height};
        }

        // Continuous legend: sampled gradient bar with low/mid/high ticks. Custom
        // bounds (lb/ub) and mid value reshape the domain shown on the bar.
        const palette = customPalette(columnName);
        const distinct = distinctColumnValues(columnName, 1);
        const paletteStops = numericPaletteStops(palette, info.min, info.max) || [
            {value: info.min, color: cssToHex(numericColor(info.min, info.min, info.max))},
            {value: info.max, color: cssToHex(numericColor(info.max, info.min, info.max))},
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
        const stops = paletteStops.map((stop, index) => {
            const fraction = ub > lb ? (stop.value - lb) / (ub - lb)
                : (index / Math.max(1, paletteStops.length - 1));
            const fill = svgColorParts(stop.color);
            return '<stop offset="' + (fraction * 100) + '%" stop-color="' + fill.color +
                '" stop-opacity="' + fill.opacity + '"/>';
        });
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
        paletteStops.forEach((stop, index) => {
            const isEnd = index === 0 || index === paletteStops.length - 1;
            const fraction = ub > lb ? (stop.value - lb) / (ub - lb) : 0.5;
            const x = barX + (fraction * barWidth);
            const text = formatValue(stop.value);
            const textWidth = (text.length * tickFontSize) / 2;
            const anchor = index === 0 ? 'start' : (isEnd ? 'end' : 'middle');
            const left = anchor === 'start' ? x : (anchor === 'end' ? x - textWidth : x - (textWidth / 2));
            if (!isEnd && left < lastRight + 4) {
                return;
            }
            tick(x, text, anchor);
            lastRight = left + textWidth;
        });
        if (nullRow) {
            swatchRow((palette && palette.nullColor) || '#b3a89d', '\u2014', tickY + padding, itemFontSize * 2, itemFontSize * 2);
        }
        parts.push('</g>');
        parts.push('</svg>');
        return {svg: parts.join('\n'), width, height};
    }

    function exportLegendSVG() {
        if (!state.bundle) {
            return;
        }
        const legend = buildLegendSVG();
        const blob = new Blob([legend.svg], {type: 'image/svg+xml'});
        const href = URL.createObjectURL(blob);
        triggerDownload(href, exportBaseName() + '_legend.svg', true);
    }

    function exportLegendPNG() {
        if (!state.bundle) {
            return;
        }
        const legend = buildLegendSVG();
        const scaleFactor = selectedPngScale();
        const svgUrl = URL.createObjectURL(new Blob([legend.svg], {type: 'image/svg+xml'}));
        const image = new Image();
        image.onload = () => {
            const target = document.createElement('canvas');
            target.width = Math.max(1, Math.round(legend.width * scaleFactor));
            target.height = Math.max(1, Math.round(legend.height * scaleFactor));
            const targetContext = target.getContext('2d');
            if (!targetContext) {
                URL.revokeObjectURL(svgUrl);
                window.alert('Could not export legend PNG at ' + scaleFactor + 'x; try a lower resolution.');
                return;
            }
            targetContext.drawImage(image, 0, 0, target.width, target.height);
            target.toBlob(blob => {
                URL.revokeObjectURL(svgUrl);
                if (!blob) {
                    window.alert('The legend is too large to export as a ' + scaleFactor + 'x PNG. Try a lower resolution.');
                    return;
                }
                triggerDownload(URL.createObjectURL(blob), exportBaseName() + '_legend.png', true);
            }, 'image/png');
        };
        image.onerror = () => {
            URL.revokeObjectURL(svgUrl);
            setStatus('Legend PNG export failed.');
        };
        image.src = svgUrl;
    }

    // ---- Custom color palette: interactive picker (modal dialog) ----

    function colorPickerIsOpen() {
        const overlay = document.getElementById('color-picker-overlay');
        return overlay ? !overlay.hidden : false;
    }

    function closeColorPicker() {
        const overlay = document.getElementById('color-picker-overlay');
        if (overlay) {
            overlay.hidden = true;
        }
    }

    function ensureCategoricalPalette(columnName) {
        let palette = state.customPalettes[paletteKey(columnName)];
        if (!palette || palette.type !== 'categorical') {
            // Seeded from the default assignment. Without that, storing one
            // hand-picked swatch would shadow the default palette and drop every
            // other value in the column back to an unassigned color.
            const base = defaultCategoricalPalette(columnName);
            palette = {
                type: 'categorical',
                colors: {...base.colors},
                nullColor: base.nullColor,
                scheme: base.scheme,
                offset: 0,
            };
            state.customPalettes[paletteKey(columnName)] = palette;
        }
        return palette;
    }

    function ensureNumericPalette(columnName) {
        let palette = state.customPalettes[paletteKey(columnName)];
        if (!palette || palette.type !== 'numeric') {
            palette = {type: 'numeric', stops: [], nullColor: null};
            state.customPalettes[paletteKey(columnName)] = palette;
        }
        return palette;
    }

    const COLOR_PICKER_PAGE_SIZE = 100;

    function renderDiscretePicker(columnName) {
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
        if (total > COLOR_PICKER_PAGE_SIZE) {
            document.getElementById('color-picker-page-status').textContent =
                'Values ' + (start + 1).toLocaleString() + '\u2013' + (start + pageValues.length).toLocaleString() +
                ' of ' + total.toLocaleString() + ' (page ' + (page + 1) + '/' + pageCount + ')';
            document.getElementById('color-picker-prev').disabled = page === 0;
            document.getElementById('color-picker-next').disabled = page >= pageCount - 1;
            pager.hidden = false;
        } else {
            pager.hidden = true;
        }

        const list = document.getElementById('color-picker-swatch-list');
        list.textContent = '';
        pageValues.forEach(entry => {
            const row = document.createElement('div');
            row.className = 'cp-swatch-row';
            const input = document.createElement('input');
            input.type = 'color';
            const override = palette && palette.colors ? palette.colors[entry.key] : undefined;
            input.value = override || cssToHex(categoricalColor(entry.raw, palette));
            input.addEventListener('input', () => {
                const edited = ensureCategoricalPalette(columnName);
                edited.colors[entry.key] = input.value;
                // No longer exactly a named palette; the menu says "Custom colors" next
                // time the picker renders.
                edited.scheme = null;
                rebuildNodeColorCache();
                scheduleRepaint();
            });
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
        });
        // Shared no-value input.
        const nullInput = document.getElementById('color-null');
        nullInput.value = (palette && palette.nullColor) || cssToHex('#b3a89d');
        nullInput.oninput = () => {
            ensureCategoricalPalette(columnName).nullColor = nullInput.value;
            rebuildNodeColorCache();
            scheduleRepaint();
        };
    }

    function parseNumberOrNull(text) {
        if (text === '' || text === null || text === undefined) {
            return null;
        }
        const value = Number(text);
        return Number.isNaN(value) ? null : value;
    }

    // Bin a numeric column's values for the gradient dialog's histogram. Bins span the
    // data range, not the low/high gradient bounds, so that values the ramp clamps stay
    // visible -- seeing them pile up flat against one end is the point of the chart.
    function columnHistogram(columnName) {
        const info = colorInfo(columnName);
        const columnIndex = state.metadataColumnIndexByName.get(columnName);
        if (!state.bundle || !info || info.baseType !== 'numeric' || columnIndex === undefined) {
            return null;
        }
        const values = [];
        let missing = 0;
        for (let nodeIndex = 0; nodeIndex < viewerNodeCount(); nodeIndex++) {
            const value = state.metadataByNodeIndex[nodeIndex]?.[columnIndex] ?? null;
            if (typeof value === 'number' && Number.isFinite(value)) {
                values.push(value);
            } else {
                missing += 1;
            }
        }
        if (values.length === 0) {
            return null;
        }
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
        if (max <= min) {
            // Every value identical: one bin, centered, rather than a zero-width domain.
            binCount = 1;
            lowEdge = min - 0.5;
            highEdge = min + 0.5;
        } else if (integer && (max - min) + 1 <= MAX_BINS) {
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
        return {counts, lowEdge, highEdge, binCount, integer, total: values.length, missing, min, max};
    }

    // Paint the cached bins, filling each bar with the color that bin's values currently
    // receive. That makes the chart double as gradient feedback: a run of flat-colored
    // bars at either end is the ramp clamping, and a pale stretch with no bars under it
    // is ramp spent on a part of the range where there is no data.
    function drawColorHistogram(columnName) {
        const wrap = document.getElementById('color-histogram-wrap');
        const histogram = state.colorHistogram;
        if (!histogram) {
            wrap.hidden = true;
            return;
        }
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

        histogram.counts.forEach((count, binIndex) => {
            if (count === 0) {
                return;
            }
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
        });

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
        if (histogram.missing > 0) {
            parts.push(histogram.missing.toLocaleString() + ' with no value');
        }
        document.getElementById('color-histogram-note').textContent = parts.join(' · ');
    }

    function updateGradientPreview() {
        const preview = document.getElementById('color-gradient-preview');
        // Sorted for display only: a stop typed past its neighbor is not
        // reordered until the field commits, and the bar should still read left
        // to right in the meantime.
        const stops = sortedGradientStops(currentGradientStops());
        if (stops.length < MIN_GRADIENT_STOPS) {
            return;
        }
        // The bar shares the histogram's axis (the column's data range) rather than
        // spanning the outermost stops, so a stop sits above the values it colors.
        // Percentages are deliberately left unclamped: CSS extends a gradient past
        // 0%/100% and holds the end colors flat outside the stops, which is exactly
        // how numericColor treats values beyond the ends.
        const domain = gradientRangeDomain();
        const parts = stops.map((stop, index) => {
            const percent = domain
                ? gradientRangePercent(stop.value, domain)
                : ((index / (stops.length - 1)) * 100);
            return stop.color + ' ' + percent + '%';
        });
        preview.style.background = 'linear-gradient(90deg, ' + parts.join(', ') + ')';
        updateGradientSlider();
    }

    function renderContinuousPicker(columnName) {
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
            {value: min, color: cssToHex(numericColor(min, min, max))},
            {value: max, color: cssToHex(numericColor(max, min, max))},
        ];
        // The labels sit under the shared axis, so they report the data extent. The
        // ramp's own ends are the Min/Max rows and the outermost knobs.
        document.getElementById('color-min-label').textContent = formatValue(min);
        document.getElementById('color-max-label').textContent = formatValue(max);
        renderGradientStopRows();
        updateGradientPreview();
        drawColorHistogram(columnName);

        nullInput.oninput = () => {
            ensureNumericPalette(columnName).nullColor = nullInput.value;
            rebuildNodeColorCache();
            scheduleRepaint();
            drawColorHistogram(columnName);
        };
        document.getElementById('color-reset-values').onclick = resetGradientStopValues;
    }

    function updateColorPaletteControl(columnName, palette, valueCount) {
        const select = document.getElementById('color-palette');
        const note = document.getElementById('color-palette-note');
        const scheme = palette && palette.scheme ? paletteSchemeByName(palette.scheme) : null;
        const hasCustomColors = !!(palette && palette.colors && Object.keys(palette.colors).length > 0);
        select.value = scheme ? scheme.name : '';
        // An offset only means something for a named palette; hand-picked colors
        // have no position in one to start from.
        const offsetInput = document.getElementById('color-palette-offset');
        const offset = scheme ? paletteOffsetValue(palette.offset) : 0;
        offsetInput.value = String(offset);
        offsetInput.disabled = !scheme;
        if (!scheme) {
            const fallback = paletteSchemeByName(DEFAULT_CATEGORICAL_PALETTE);
            const builtWithTable = !!state.bundle.defaults.color_table
                && columnName === state.bundle.defaults.color_by;
            note.textContent = (hasCustomColors
                ? 'Colors were set by hand or loaded from a color table.'
                : 'No palette assigned.') +
                (builtWithTable
                    ? ' Reset to defaults returns this column to the color table it was built with.'
                    : (fallback ? ' Reset to defaults returns this column to ' + fallback.label + '.' : ''));
            return;
        }
        const cycles = valueCount > scheme.colors.length;
        const startNote = offset > 0
            ? ', starting at color ' + ((offset % scheme.colors.length) + 1).toLocaleString()
            : '';
        note.textContent = scheme.note + ' ' + scheme.colors.length.toLocaleString() + ' colors for ' +
            valueCount.toLocaleString() + ' value' + (valueCount === 1 ? '' : 's') + startNote +
            (cycles ? ', so colors repeat.' : '.');
    }

    // Serves both the palette menu and the start offset: either one re-applies the
    // menu's palette from the offset.
    function changeColorPalette() {
        const columnName = currentColorField();
        if (!columnName) {
            return;
        }
        // Every entry in the menu is a named palette now, so there is no
        // pseudo-scheme to special-case; "Reset to defaults" handles going back.
        if (!applyNamedPalette(columnName, document.getElementById('color-palette').value,
            document.getElementById('color-palette-offset').value)) {
            return;
        }
        rebuildNodeColorCache();
        repaintNow();
        renderDiscretePicker(columnName);
    }

    // Step the discrete swatch pager (no-op for numeric/empty columns).
    function changeColorPickerPage(delta) {
        const columnName = currentColorField();
        const info = colorInfo(columnName);
        if (!columnName || !info || info.type === 'numeric') {
            return;
        }
        state.colorPickerPage += delta;
        renderDiscretePicker(columnName);
    }

    function openColorPicker() {
        if (!state.bundle) {
            return;
        }
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
        if (!columnName || !info) {
            discretePanel.hidden = true;
            continuousPanel.hidden = true;
            nullControl.hidden = true;
            emptyNote.hidden = false;
        } else if (info.type === 'numeric') {
            emptyNote.hidden = true;
            discretePanel.hidden = true;
            continuousPanel.hidden = false;
            nullControl.hidden = false;
            renderContinuousPicker(columnName);
        } else {
            emptyNote.hidden = true;
            continuousPanel.hidden = true;
            discretePanel.hidden = false;
            nullControl.hidden = false;
            renderDiscretePicker(columnName);
        }
        overlay.hidden = false;
    }

    // The gradient/discrete switch is only meaningful for numeric columns, so it is
    // shown (and its distinct-value hint filled in) only for those.
    function updateCategoricalToggle(columnName) {
        const control = document.getElementById('color-categorical-control');
        const toggle = document.getElementById('color-as-categorical');
        const note = document.getElementById('color-categorical-note');
        if (!columnName || !columnIsNumericType(columnName)) {
            control.hidden = true;
            note.hidden = true;
            toggle.checked = false;
            return;
        }
        control.hidden = false;
        toggle.checked = columnIsCategoricalNumeric(columnName);
        if (toggle.checked) {
            const total = distinctColumnValues(columnName, 0).total;
            note.textContent = 'Coloring ' + total.toLocaleString() + ' distinct value' +
                (total === 1 ? '' : 's') + ' as separate categories.';
            note.hidden = false;
        } else {
            note.hidden = true;
        }
    }

    function toggleColorCategorical() {
        const columnName = currentColorField();
        if (!columnName || !columnIsNumericType(columnName)) {
            return;
        }
        setColumnCategorical(columnName, document.getElementById('color-as-categorical').checked);
        rebuildNodeColorCache();
        repaintNow();
        openColorPicker();
    }

    function resetColorPicker() {
        const columnName = currentColorField();
        if (columnName) {
            delete state.customPalettes[paletteKey(columnName)];
            rebuildNodeColorCache();
            scheduleRepaint();
            openColorPicker();
        }
    }

