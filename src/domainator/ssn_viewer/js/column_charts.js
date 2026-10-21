/*
 * Per-column charts and frequency tables from the metadata table's headers.
 *
 * Included into viewer.js's script scope; see session_state.js.
 */

    // ---------------------------------------------------------------------
    // Per-column charts and tables
    //
    // Each metadata column header carries a chart glyph that opens a menu of
    // chart kinds appropriate to that column, and picking one opens a dialog.
    //
    // Two invariants hold the feature together:
    //
    // 1. Every chart summarizes exactly the rows the header's own Copy button
    //    copies (columnChartNodeIndices), so a chart and a copied column can
    //    never disagree about what "these nodes" means.
    // 2. Every chart is one SVG string, consumed by three sinks: innerHTML for
    //    the preview, a Blob for the SVG export, and Image->canvas->toBlob for
    //    the PNG export. A canvas preview would need a second renderer for the
    //    SVG export, which is exactly the drift buildClusterViewSVG warns
    //    about. Charts are tens to hundreds of marks, so SVG DOM cost is
    //    irrelevant here (unlike renderClusterView's thousands of nodes) and
    //    axis labels come out as real selectable text.
    //
    // Colors are inherited from the column's own palette via customPalette(),
    // which routes through paletteKey() -- so a numeric column flipped to
    // discrete coloring gets its categorical palette, and a column that has
    // never been the color-by column still gets the same colors it would get if
    // it were (DEFAULT_CATEGORICAL_PALETTE, when nothing has been customized).
    // ---------------------------------------------------------------------

    const COLUMN_CHART_PIE_MAX_SLICES = 12;
    const COLUMN_CHART_BAR_MAX_BARS = 24;
    // The frequency table is a table, so it can afford many more rows than a
    // chart can afford marks -- but not unbounded ones, since it goes through
    // innerHTML. The TSV export always carries every value.
    const COLUMN_CHART_TABLE_MAX_ROWS = 500;
    // An ECDF over a large selection would emit one step per value. Sampling
    // caps the path at a size an SVG editor can still open; the curve is
    // visually identical because adjacent steps land on the same pixel.
    const COLUMN_CHART_ECDF_MAX_POINTS = 2000;

    // Single source of truth for the menu, the dialog's kind <select>, and for
    // validating a kind carried over when the dialog reopens on another column.
    const COLUMN_CHART_KINDS = [
        {id: 'frequency', label: 'Frequency table', numericOnly: false},
        {id: 'bar', label: 'Bar chart', numericOnly: false},
        {id: 'pie', label: 'Pie chart', numericOnly: false},
        {id: 'histogram', label: 'Histogram', numericOnly: true},
        {id: 'box', label: 'Box plot', numericOnly: true},
        {id: 'ecdf', label: 'Cumulative distribution', numericOnly: true},
        {id: 'summary', label: 'Summary statistics', numericOnly: false},
    ];

    // Histogram/box/ECDF are offered for any int/float column even while it is
    // colored as discrete categories: that toggle is a coloring choice, not a
    // claim that the numbers are not numbers.
    function columnChartKindsFor(columnName) {
        const numeric = columnIsNumericType(columnName);
        return COLUMN_CHART_KINDS.filter(kind => numeric || !kind.numericOnly);
    }

    function columnChartKindLabel(kindId) {
        const kind = COLUMN_CHART_KINDS.find(entry => entry.id === kindId);
        return kind ? kind.label : 'Chart';
    }

    // ---- Scope ----

    // The rows every chart summarizes: the table's current rows after the
    // search filter, across all pages. Identical to what copyMetadataColumn
    // exports, and metadataBaseNodeIndices already encodes "the selection, or
    // every node when nothing is selected". Pagination is a rendering artifact,
    // so charting only the visible page would be a trap.
    function columnChartNodeIndices() {
        return metadataDisplayNodeIndices(metadataBaseNodeIndices());
    }

    // Baked into every chart's subtitle so an exported SVG/PNG says what it
    // counted, not just how many.
    function columnChartScopeLabel(rowCount) {
        const total = state.bundle ? state.bundle.graph.nodes.length : 0;
        const parts = [rowCount.toLocaleString() + ' of ' + total.toLocaleString() + ' nodes'];
        parts.push(state.selectedNodeIndices.size > 0 ? 'selection' : 'all nodes');
        const filterText = metadataFilterText();
        if (filterText) {
            parts.push('filter "' + filterText + '"');
        }
        return parts.join(' · ');
    }

    // ---- Data model ----

    // Counts for the same column over every node in the bundle, not just the
    // scoped rows. A chart's own percentages read one way -- "a third of this
    // selection is alpha" -- and a selection immediately raises the other -- "and
    // that is a twelfth of every alpha in the network". Only these counts can
    // answer the second question, because the scope cannot see what it excluded.
    //
    // Recomputed per chart rather than cached: metadata is editable, and one pass
    // over a column costs what the scope's own pass beside it already costs.
    function columnGlobalCounts(columnName) {
        const counts = {byKey: new Map(), nullCount: 0, total: 0};
        const columnIndex = state.metadataColumnIndexByName.get(columnName);
        if (columnIndex === undefined || !state.bundle) {
            return counts;
        }
        counts.total = state.bundle.graph.nodes.length;
        for (let nodeIndex = 0; nodeIndex < counts.total; nodeIndex++) {
            const raw = state.metadataByNodeIndex[nodeIndex]?.[columnIndex] ?? null;
            if (isMissingMetadataValue(raw)) {
                counts.nullCount += 1;
                continue;
            }
            const key = String(raw);
            counts.byKey.set(key, (counts.byKey.get(key) || 0) + 1);
        }
        return counts;
    }

    // Scoped analog of distinctColumnValues (which walks every node): same
    // String(raw) keying, same formatValue labels, and the same count-desc then
    // key-asc ordering, so a chart's category order matches the color picker's
    // swatch order and the legend's rows.
    //
    // options.topN caps the entry count and rolls the tail into one "Other"
    // entry (the pie and bar charts cap; the frequency table lists everything).
    // options.includeNull appends a trailing no-value entry, which every chart
    // wants and only the summary table, which counts nulls in its own row,
    // turns off.
    function columnValueDistribution(columnName, nodeIndices, options = {}) {
        const topN = options.topN === undefined ? Infinity : options.topN;
        const includeNull = options.includeNull !== false;
        const palette = customPalette(columnName);
        const model = {
            entries: [],
            nullCount: 0,
            otherCount: 0,
            otherDistinct: 0,
            distinctTotal: 0,
            rowCount: nodeIndices.length,
            truncated: false,
            globalTotal: 0,
        };
        const columnIndex = state.metadataColumnIndexByName.get(columnName);
        if (columnIndex === undefined) {
            return model;
        }
        const byKey = new Map();
        nodeIndices.forEach(nodeIndex => {
            const raw = state.metadataByNodeIndex[nodeIndex]?.[columnIndex] ?? null;
            if (isMissingMetadataValue(raw)) {
                model.nullCount += 1;
                return;
            }
            const key = String(raw);
            const existing = byKey.get(key);
            if (existing) {
                existing.count += 1;
            } else {
                byKey.set(key, {key, label: formatValue(raw), raw, count: 1});
            }
        });
        model.distinctTotal = byKey.size;
        const globals = columnGlobalCounts(columnName);
        model.globalTotal = globals.total;
        const sorted = Array.from(byKey.values()).sort(
            (a, b) => (b.count - a.count) || a.key.localeCompare(b.key)
        );
        const kept = Number.isFinite(topN) ? sorted.slice(0, topN) : sorted;
        const tail = Number.isFinite(topN) ? sorted.slice(topN) : [];
        model.entries = kept.map(entry => Object.assign({}, entry, {
            color: categoricalColor(entry.raw, palette),
            globalCount: globals.byKey.get(entry.key) || 0,
            isNull: false,
            isOther: false,
        }));
        if (tail.length > 0) {
            model.otherDistinct = tail.length;
            model.otherCount = tail.reduce((sum, entry) => sum + entry.count, 0);
            model.truncated = true;
            // The two synthetic entries below get NUL-prefixed keys so they can
            // never collide with a real value's String(raw), the same trick
            // paletteKey uses to keep a column's two palettes apart.
            model.entries.push({
                key: '\u0000other',
                label: 'Other (' + tail.length.toLocaleString() + ' values)',
                raw: null,
                count: model.otherCount,
                // The rollup's own global count is the sum over the values rolled
                // up, not every value outside the top N: "Other" is defined by what
                // this scope saw, and a value absent from the scope is not in it.
                globalCount: tail.reduce((sum, entry) => sum + (globals.byKey.get(entry.key) || 0), 0),
                // The same neutral gray utils.get_palette gives values with no
                // color of their own, so a rolled-up slice reads as "not a
                // value" rather than as a category with that hue.
                color: PALETTE_NO_VALUE_COLOR,
                isNull: false,
                isOther: true,
            });
        }
        if (includeNull && model.nullCount > 0) {
            model.entries.push({
                key: '\u0000null',
                label: '—',
                raw: null,
                count: model.nullCount,
                globalCount: globals.nullCount,
                // categoricalColor's own null fallback, so an empty cell is the
                // same color in the chart as its node is on the canvas. Note
                // this is a different gray from the "Other" rollup above.
                color: (palette && palette.nullColor) || '#b3a89d',
                isNull: true,
                isOther: false,
            });
        }
        return model;
    }

    function columnNumericValues(columnName, nodeIndices) {
        const result = {values: [], missing: 0};
        const columnIndex = state.metadataColumnIndexByName.get(columnName);
        if (columnIndex === undefined) {
            return result;
        }
        nodeIndices.forEach(nodeIndex => {
            const value = state.metadataByNodeIndex[nodeIndex]?.[columnIndex] ?? null;
            if (typeof value === 'number' && Number.isFinite(value)) {
                result.values.push(value);
            } else {
                result.missing += 1;
            }
        });
        result.values.sort((left, right) => left - right);
        return result;
    }

    // Same shape as columnHistogram, and the same binning rules, but over the
    // scoped rows with the scoped min/max. The two are deliberately separate:
    // columnHistogram bins the *whole column* because the gradient dialog's
    // slider knobs ride that axis, so narrowing it to a selection would move
    // the knobs out from under the user.
    //
    // colorMin/colorMax carry the whole-column range through anyway, so bars
    // are colored the way the nodes are colored (see buildHistogramSVG).
    function scopedColumnHistogram(columnName, nodeIndices) {
        const info = colorInfo(columnName);
        if (!state.bundle || !info || info.baseType !== 'numeric') {
            return null;
        }
        const numeric = columnNumericValues(columnName, nodeIndices);
        const values = numeric.values;
        if (values.length === 0) {
            return null;
        }
        const column = state.metadataColumnByName.get(columnName);
        const min = values[0];
        const max = values[values.length - 1];
        const MAX_BINS = 48;
        const integer = !!(column && column.type === 'int'
            && Number.isInteger(min) && Number.isInteger(max));
        let binCount;
        let lowEdge;
        let highEdge;
        if (max <= min) {
            binCount = 1;
            lowEdge = min - 0.5;
            highEdge = min + 0.5;
        } else if (integer && (max - min) + 1 <= MAX_BINS) {
            // One bin per integer: evenly dividing a handful of distinct
            // integers leaves empty gaps and doubled-up bars that read as
            // structure when they are an artifact of the binning.
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
        return {
            counts, lowEdge, highEdge, binCount, integer,
            total: values.length, missing: numeric.missing, min, max,
            colorMin: info.min, colorMax: info.max,
        };
    }

    // Quantiles by linear interpolation between order statistics (numpy's
    // default, and R's type 7), so the numbers match what a user would get
    // running the exported TSV through pandas. `values` must be sorted.
    function quantileOf(values, fraction) {
        if (values.length === 0) {
            return NaN;
        }
        if (values.length === 1) {
            return values[0];
        }
        const position = (values.length - 1) * fraction;
        const lower = Math.floor(position);
        const upper = Math.ceil(position);
        if (lower === upper) {
            return values[lower];
        }
        return values[lower] + ((values[upper] - values[lower]) * (position - lower));
    }

    // Five-number summary plus mean/SD and Tukey 1.5*IQR whiskers, shared by
    // the box plot and the summary table so the two can never disagree.
    // `values` must be sorted ascending.
    function numericSummary(values) {
        if (values.length === 0) {
            return null;
        }
        const count = values.length;
        const min = values[0];
        const max = values[count - 1];
        const q1 = quantileOf(values, 0.25);
        const median = quantileOf(values, 0.5);
        const q3 = quantileOf(values, 0.75);
        const iqr = q3 - q1;
        const mean = values.reduce((sum, value) => sum + value, 0) / count;
        // Sample standard deviation (n-1), undefined for a single value.
        const variance = count > 1
            ? values.reduce((sum, value) => sum + ((value - mean) ** 2), 0) / (count - 1)
            : 0;
        const stdev = count > 1 ? Math.sqrt(variance) : NaN;
        const lowFence = q1 - (1.5 * iqr);
        const highFence = q3 + (1.5 * iqr);
        // Whiskers reach the most extreme value still inside the fences, never
        // past the data -- so a distribution with no outliers draws whiskers at
        // min and max. Seeded crossed so the first value inside wins; at least
        // one always is, since q1 <= median <= q3 sits inside both fences.
        let lowWhisker = max;
        let highWhisker = min;
        const outliers = [];
        values.forEach(value => {
            if (value < lowFence || value > highFence) {
                outliers.push(value);
                return;
            }
            lowWhisker = Math.min(lowWhisker, value);
            highWhisker = Math.max(highWhisker, value);
        });
        return {
            count, min, q1, median, q3, max, iqr, mean, stdev,
            lowWhisker, highWhisker, outliers,
        };
    }

    // Mirrors formatMetadataDisplayValue's float branch, but is never routed
    // through a column's declared type: the quartiles, mean and SD of an int
    // column are not integers, and the int formatter would truncate a first
    // quartile of 2.25 to "2".
    function formatStatNumber(value) {
        if (!Number.isFinite(value)) {
            return '—';
        }
        if (Number.isInteger(value)) {
            return value.toLocaleString();
        }
        const magnitude = Math.abs(value);
        if (magnitude >= 10000 || magnitude < 0.001) {
            return value.toExponential(3);
        }
        return new Intl.NumberFormat(undefined, {maximumFractionDigits: 4}).format(value);
    }

    // [label, text] rows for the summary table.
    function columnSummaryRows(columnName, nodeIndices) {
        const model = columnValueDistribution(columnName, nodeIndices, {includeNull: false});
        const rows = [
            ['Rows', model.rowCount.toLocaleString()],
            ['With a value', (model.rowCount - model.nullCount).toLocaleString()],
            ['No value', model.nullCount.toLocaleString()],
            ['Distinct values', model.distinctTotal.toLocaleString()],
        ];
        if (model.entries.length > 0) {
            const top = model.entries[0];
            rows.push(['Most common', top.label + ' (' + top.count.toLocaleString() + ')']);
        }
        if (!columnIsNumericType(columnName)) {
            return rows;
        }
        const summary = numericSummary(columnNumericValues(columnName, nodeIndices).values);
        if (!summary) {
            return rows;
        }
        rows.push(['Minimum', formatStatNumber(summary.min)]);
        rows.push(['1st quartile', formatStatNumber(summary.q1)]);
        rows.push(['Median', formatStatNumber(summary.median)]);
        rows.push(['3rd quartile', formatStatNumber(summary.q3)]);
        rows.push(['Maximum', formatStatNumber(summary.max)]);
        rows.push(['Mean', formatStatNumber(summary.mean)]);
        rows.push(['Std. deviation', formatStatNumber(summary.stdev)]);
        return rows;
    }

    // {color, label, count, percent, globalCount, globalPercent} rows shared by
    // the HTML preview, the SVG export and the TSV, so all three show the same
    // numbers. The two percentages have different denominators on purpose:
    // `percent` is the value's share of the charted rows, `globalPercent` is the
    // share of that value's own network-wide population which those rows caught.
    // With nothing selected and no filter the scope is the whole bundle, so every
    // globalPercent is 100% -- which is the honest reading, not a bug.
    function frequencyTableRows(model) {
        const denominator = model.rowCount > 0 ? model.rowCount : 1;
        return model.entries.map(entry => {
            const globalCount = entry.globalCount || 0;
            return {
                color: entry.color,
                label: entry.label,
                count: entry.count,
                percent: (entry.count / denominator) * 100,
                globalCount,
                globalPercent: globalCount > 0 ? (entry.count / globalCount) * 100 : 0,
            };
        });
    }

    // ---- SVG scaffolding ----

    const CHART_TITLE_FONT = 18;
    const CHART_SUBTITLE_FONT = 12;
    const CHART_LABEL_FONT = 13;
    const CHART_PAD = 14;
    // Tall enough to clear the title and subtitle baselines with room to spare:
    // the histogram and ECDF label their topmost gridline, which sits here.
    const CHART_HEAD = 64;

    // The same width approximation buildLegendSVG uses. Measuring text properly
    // would mean laying it out in the DOM first, which the export path (a
    // detached string) cannot do.
    function estimateTextWidth(text, fontSize) {
        return (String(text).length * fontSize) / 2;
    }

    function truncateChartLabel(text, maxChars) {
        const value = String(text);
        return value.length <= maxChars ? value : (value.slice(0, Math.max(1, maxChars - 1)) + '…');
    }

    // Emit one color as an SVG-editor-safe fill attribute pair. svgColorParts
    // normalizes through the canvas parser (the default schemes emit hsl(),
    // which Illustrator and Inkscape choke on) and splits any alpha out into
    // its own attribute.
    function chartSvgFill(color) {
        const parts = svgColorParts(color);
        return 'fill="' + parts.color + '" fill-opacity="' + parts.opacity + '"';
    }

    // Deliberately no <?xml ...?> prolog: the same string is assigned to
    // innerHTML for the preview, which refuses a processing instruction. The
    // export path prepends it (see exportColumnChartSVG).
    function chartSvgOpen(parts, width, height, title, subtitle) {
        parts.push('<svg xmlns="http://www.w3.org/2000/svg" width="' + width + '" height="' + height +
            '" viewBox="0 0 ' + width + ' ' + height + '">');
        parts.push('<g font-family="Arial, sans-serif">');
        parts.push('<rect x="0" y="0" width="' + width + '" height="' + height + '" fill="white"/>');
        parts.push('<text x="' + CHART_PAD + '" y="' + (CHART_PAD + CHART_TITLE_FONT - 4) +
            '" font-size="' + CHART_TITLE_FONT + '" font-weight="bold">' + escapeXml(title) + '</text>');
        parts.push('<text x="' + CHART_PAD + '" y="' + (CHART_PAD + CHART_TITLE_FONT + CHART_SUBTITLE_FONT + 2) +
            '" font-size="' + CHART_SUBTITLE_FONT + '" fill="#5c6a70">' + escapeXml(subtitle) + '</text>');
    }

    function chartSvgClose(parts) {
        parts.push('</g>');
        parts.push('</svg>');
        return parts.join('\n');
    }

    // "Nice" 1/2/2.5/5 tick values spanning [min, max]. Used by the histogram,
    // box plot and ECDF, all of which need a readable numeric axis rather than
    // the legend's fixed end ticks.
    function chartAxisTicks(min, max, targetCount) {
        if (!Number.isFinite(min) || !Number.isFinite(max) || max <= min) {
            return [min];
        }
        const rawStep = (max - min) / Math.max(1, targetCount);
        const magnitude = Math.pow(10, Math.floor(Math.log10(rawStep)));
        const normalized = rawStep / magnitude;
        // 2.5 is on the ladder because a bare 1/2/5/10 rounds a raw step of
        // 2.2 all the way up to 5, leaving a wide axis with two ticks on it.
        const step = magnitude * (normalized <= 1 ? 1
            : (normalized <= 2 ? 2 : (normalized <= 2.5 ? 2.5 : (normalized <= 5 ? 5 : 10))));
        const ticks = [];
        const first = Math.ceil(min / step) * step;
        for (let value = first; value <= max + (step / 1000); value += step) {
            // Re-round to kill the float drift that accumulates over the loop.
            ticks.push(Number((Math.round(value / step) * step).toPrecision(12)));
        }
        if (ticks.length === 0) {
            return [min, max];
        }
        return ticks;
    }

    // SVG arc path for a pie slice. Angles are in radians, measured clockwise
    // from twelve o'clock so slices read in the order the legend lists them.
    function describeArcPath(cx, cy, radius, startAngle, endAngle) {
        const point = angle => {
            const x = cx + (radius * Math.sin(angle));
            const y = cy - (radius * Math.cos(angle));
            return x.toFixed(2) + ' ' + y.toFixed(2);
        };
        const largeArc = (endAngle - startAngle) > Math.PI ? 1 : 0;
        return 'M ' + cx.toFixed(2) + ' ' + cy.toFixed(2) +
            ' L ' + point(startAngle) +
            ' A ' + radius + ' ' + radius + ' 0 ' + largeArc + ' 1 ' + point(endAngle) + ' Z';
    }

    // ---- Chart builders ----
    // Each returns {svg, width, height}.

    function buildBarChartSVG(model, meta) {
        const rows = frequencyTableRows(model);
        const barHeight = 22;
        const barGap = 6;
        const labelChars = 28;
        const labelWidth = Math.min(
            260,
            Math.max(90, rows.reduce(
                (widest, row) => Math.max(widest, estimateTextWidth(truncateChartLabel(row.label, labelChars), CHART_LABEL_FONT)),
                0
            ) + 8)
        );
        const countWidth = 96;
        const plotWidth = 420;
        const width = CHART_PAD + labelWidth + plotWidth + countWidth + CHART_PAD;
        const height = CHART_HEAD + (rows.length * (barHeight + barGap)) + CHART_PAD;
        const parts = [];
        chartSvgOpen(parts, width, height, meta.title, meta.subtitle);
        const maxCount = rows.reduce((most, row) => Math.max(most, row.count), 0) || 1;
        const barX = CHART_PAD + labelWidth;
        rows.forEach((row, index) => {
            const y = CHART_HEAD + (index * (barHeight + barGap));
            parts.push('<text x="' + (barX - 8) + '" y="' + (y + (barHeight / 2)) +
                '" font-size="' + CHART_LABEL_FONT + '" text-anchor="end" dominant-baseline="central">' +
                escapeXml(truncateChartLabel(row.label, labelChars)) + '</text>');
            // A floor of 1px keeps a count of 1 visible next to a dominant bar.
            const barWidth = Math.max(1, (row.count / maxCount) * plotWidth);
            parts.push('<rect x="' + barX + '" y="' + y + '" width="' + barWidth.toFixed(2) +
                '" height="' + barHeight + '" ' + chartSvgFill(row.color) +
                ' stroke="#1e2a2f" stroke-width="0.75"/>');
            parts.push('<text x="' + (barX + barWidth + 8) + '" y="' + (y + (barHeight / 2)) +
                '" font-size="' + CHART_LABEL_FONT + '" fill="#5c6a70" dominant-baseline="central">' +
                escapeXml(row.count.toLocaleString() + ' (' + row.percent.toFixed(1) + '%)') + '</text>');
        });
        return {svg: chartSvgClose(parts), width, height};
    }

    function buildPieChartSVG(model, meta) {
        const rows = frequencyTableRows(model);
        const radius = 130;
        const diameter = radius * 2;
        const swatch = 14;
        const legendRowHeight = 22;
        const labelChars = 30;
        const legendWidth = Math.max(140, rows.reduce(
            (widest, row) => Math.max(widest, estimateTextWidth(
                truncateChartLabel(row.label, labelChars) + '  ' + row.count.toLocaleString() + ' (' + row.percent.toFixed(1) + '%)',
                CHART_LABEL_FONT
            )),
            0
        ) + swatch + 16);
        const width = CHART_PAD + diameter + 24 + legendWidth + CHART_PAD;
        const height = Math.max(
            CHART_HEAD + diameter + CHART_PAD,
            CHART_HEAD + (rows.length * legendRowHeight) + CHART_PAD
        );
        const parts = [];
        chartSvgOpen(parts, width, height, meta.title, meta.subtitle);
        const cx = CHART_PAD + radius;
        const cy = CHART_HEAD + radius;
        const total = rows.reduce((sum, row) => sum + row.count, 0);
        if (total > 0) {
            if (rows.length === 1) {
                // A single 100% slice: an arc whose start and end angles are
                // equal draws nothing, so the whole circle has to be a circle.
                parts.push('<circle cx="' + cx + '" cy="' + cy + '" r="' + radius + '" ' +
                    chartSvgFill(rows[0].color) + ' stroke="#1e2a2f" stroke-width="0.75"/>');
            } else {
                let angle = 0;
                rows.forEach(row => {
                    const sweep = (row.count / total) * Math.PI * 2;
                    parts.push('<path d="' + describeArcPath(cx, cy, radius, angle, angle + sweep) + '" ' +
                        chartSvgFill(row.color) + ' stroke="#1e2a2f" stroke-width="0.75"/>');
                    angle += sweep;
                });
            }
        }
        const legendX = CHART_PAD + diameter + 24;
        rows.forEach((row, index) => {
            const y = CHART_HEAD + (index * legendRowHeight);
            parts.push('<rect x="' + legendX + '" y="' + (y + ((legendRowHeight - swatch) / 2)) +
                '" width="' + swatch + '" height="' + swatch + '" ' + chartSvgFill(row.color) +
                ' stroke="#1e2a2f" stroke-width="0.75"/>');
            parts.push('<text x="' + (legendX + swatch + 8) + '" y="' + (y + (legendRowHeight / 2)) +
                '" font-size="' + CHART_LABEL_FONT + '" dominant-baseline="central">' +
                escapeXml(truncateChartLabel(row.label, labelChars) + '  ' +
                    row.count.toLocaleString() + ' (' + row.percent.toFixed(1) + '%)') + '</text>');
        });
        return {svg: chartSvgClose(parts), width, height};
    }

    // Bars are filled with the color their bin's values receive on the network
    // canvas -- so this doubles as a check on the column's gradient. The color
    // domain is the whole column (colorMin/colorMax), not the scoped range,
    // because that is what nodeColor uses.
    function buildHistogramSVG(histogram, meta) {
        const palette = customPalette(meta.columnName);
        const plotWidth = 620;
        const plotHeight = 240;
        const axisSpace = 44;
        const yLabelWidth = 56;
        const width = CHART_PAD + yLabelWidth + plotWidth + CHART_PAD;
        const height = CHART_HEAD + plotHeight + axisSpace;
        const parts = [];
        chartSvgOpen(parts, width, height, meta.title, meta.subtitle);
        const plotX = CHART_PAD + yLabelWidth;
        const baselineY = CHART_HEAD + plotHeight;
        const maxCount = Math.max(...histogram.counts) || 1;
        const slotWidth = plotWidth / histogram.binCount;
        const gap = slotWidth > 8 ? Math.min(3, slotWidth * 0.16) : 0;
        const binSpan = (histogram.highEdge - histogram.lowEdge) / histogram.binCount;

        // Y axis: zero and the peak count are the two numbers that give the
        // bars their scale.
        [0, maxCount].forEach(count => {
            const y = baselineY - ((count / maxCount) * plotHeight);
            parts.push('<line x1="' + plotX + '" y1="' + y.toFixed(2) + '" x2="' + (plotX + plotWidth) +
                '" y2="' + y.toFixed(2) + '" stroke="#d8dce2" stroke-width="1"/>');
            parts.push('<text x="' + (plotX - 8) + '" y="' + y.toFixed(2) +
                '" font-size="' + CHART_LABEL_FONT + '" fill="#5c6a70" text-anchor="end" dominant-baseline="central">' +
                escapeXml(count.toLocaleString()) + '</text>');
        });

        histogram.counts.forEach((count, binIndex) => {
            if (count === 0) {
                return;
            }
            // A floor of 2px keeps rare bins from vanishing beside a dominant one.
            const barHeight = Math.max(2, (count / maxCount) * plotHeight);
            const center = histogram.lowEdge + ((binIndex + 0.5) * binSpan);
            const color = numericColor(center, histogram.colorMin, histogram.colorMax, palette);
            parts.push('<rect x="' + (plotX + (binIndex * slotWidth) + (gap / 2)).toFixed(2) +
                '" y="' + (baselineY - barHeight).toFixed(2) +
                '" width="' + Math.max(1, slotWidth - gap).toFixed(2) +
                '" height="' + barHeight.toFixed(2) + '" ' + chartSvgFill(color) + '/>');
        });

        parts.push('<line x1="' + plotX + '" y1="' + baselineY + '" x2="' + (plotX + plotWidth) +
            '" y2="' + baselineY + '" stroke="#1e2a2f" stroke-width="1.25"/>');

        // X ticks, dropping any that would collide with the label before it --
        // the same pruning buildLegendSVG does on its gradient ticks.
        const tickFont = 12;
        let lastRight = -Infinity;
        chartAxisTicks(histogram.lowEdge, histogram.highEdge, 6).forEach(value => {
            const fraction = (value - histogram.lowEdge) / (histogram.highEdge - histogram.lowEdge);
            const x = plotX + (fraction * plotWidth);
            const text = formatValue(value);
            const textWidth = estimateTextWidth(text, tickFont);
            if (x - (textWidth / 2) < lastRight + 6) {
                return;
            }
            parts.push('<line x1="' + x.toFixed(2) + '" y1="' + baselineY + '" x2="' + x.toFixed(2) +
                '" y2="' + (baselineY + 5) + '" stroke="#1e2a2f" stroke-width="1"/>');
            parts.push('<text x="' + x.toFixed(2) + '" y="' + (baselineY + 5 + tickFont + 2) +
                '" font-size="' + tickFont + '" text-anchor="middle">' + escapeXml(text) + '</text>');
            lastRight = x + (textWidth / 2);
        });
        parts.push('<text x="' + (plotX + (plotWidth / 2)) + '" y="' + (height - 4) +
            '" font-size="' + tickFont + '" fill="#5c6a70" text-anchor="middle">' +
            escapeXml(meta.columnName) + '</text>');
        return {svg: chartSvgClose(parts), width, height};
    }

    function buildBoxPlotSVG(summary, meta) {
        const palette = customPalette(meta.columnName);
        const info = colorInfo(meta.columnName);
        const colorMin = info ? info.min : summary.min;
        const colorMax = info ? info.max : summary.max;
        const plotWidth = 620;
        const boxHeight = 74;
        const axisSpace = 44;
        const width = CHART_PAD + plotWidth + CHART_PAD;
        const height = CHART_HEAD + boxHeight + axisSpace;
        const parts = [];
        chartSvgOpen(parts, width, height, meta.title, meta.subtitle);
        const plotX = CHART_PAD;
        // Pad the domain so a whisker or outlier at the extreme is not clipped
        // by the plot edge.
        const dataLow = Math.min(summary.lowWhisker, summary.min);
        const dataHigh = Math.max(summary.highWhisker, summary.max);
        const domainSpan = dataHigh - dataLow;
        const pad = domainSpan > 0 ? domainSpan * 0.04 : 0.5;
        const low = dataLow - pad;
        const high = dataHigh + pad;
        const xOf = value => plotX + (((value - low) / (high - low)) * plotWidth);
        const midY = CHART_HEAD + (boxHeight / 2);

        parts.push('<line x1="' + xOf(summary.lowWhisker).toFixed(2) + '" y1="' + midY +
            '" x2="' + xOf(summary.highWhisker).toFixed(2) + '" y2="' + midY +
            '" stroke="#1e2a2f" stroke-width="1.25"/>');
        [summary.lowWhisker, summary.highWhisker].forEach(value => {
            parts.push('<line x1="' + xOf(value).toFixed(2) + '" y1="' + (midY - 14) +
                '" x2="' + xOf(value).toFixed(2) + '" y2="' + (midY + 14) +
                '" stroke="#1e2a2f" stroke-width="1.25"/>');
        });
        const boxLeft = xOf(summary.q1);
        const boxRight = xOf(summary.q3);
        parts.push('<rect x="' + boxLeft.toFixed(2) + '" y="' + (midY - 26) +
            '" width="' + Math.max(1, boxRight - boxLeft).toFixed(2) + '" height="52" ' +
            chartSvgFill(numericColor(summary.median, colorMin, colorMax, palette)) +
            ' stroke="#1e2a2f" stroke-width="1.25"/>');
        parts.push('<line x1="' + xOf(summary.median).toFixed(2) + '" y1="' + (midY - 26) +
            '" x2="' + xOf(summary.median).toFixed(2) + '" y2="' + (midY + 26) +
            '" stroke="#1e2a2f" stroke-width="2.5"/>');
        summary.outliers.forEach(value => {
            parts.push('<circle cx="' + xOf(value).toFixed(2) + '" cy="' + midY + '" r="3.5" ' +
                chartSvgFill(numericColor(value, colorMin, colorMax, palette)) +
                ' stroke="#1e2a2f" stroke-width="0.75"/>');
        });

        const baselineY = CHART_HEAD + boxHeight;
        parts.push('<line x1="' + plotX + '" y1="' + baselineY + '" x2="' + (plotX + plotWidth) +
            '" y2="' + baselineY + '" stroke="#1e2a2f" stroke-width="1.25"/>');
        const tickFont = 12;
        let lastRight = -Infinity;
        chartAxisTicks(low, high, 6).forEach(value => {
            const x = xOf(value);
            const text = formatValue(value);
            const textWidth = estimateTextWidth(text, tickFont);
            if (x - (textWidth / 2) < lastRight + 6) {
                return;
            }
            parts.push('<line x1="' + x.toFixed(2) + '" y1="' + baselineY + '" x2="' + x.toFixed(2) +
                '" y2="' + (baselineY + 5) + '" stroke="#1e2a2f" stroke-width="1"/>');
            parts.push('<text x="' + x.toFixed(2) + '" y="' + (baselineY + 5 + tickFont + 2) +
                '" font-size="' + tickFont + '" text-anchor="middle">' + escapeXml(text) + '</text>');
            lastRight = x + (textWidth / 2);
        });
        parts.push('<text x="' + (plotX + (plotWidth / 2)) + '" y="' + (height - 4) +
            '" font-size="' + tickFont + '" fill="#5c6a70" text-anchor="middle">' +
            escapeXml(meta.columnName) + '</text>');
        return {svg: chartSvgClose(parts), width, height};
    }

    // Empirical cumulative distribution: the fraction of values at or below
    // each x. The right chart for picking a cutoff, since you can read the
    // proportion kept straight off the curve.
    function buildEcdfSVG(values, meta) {
        const plotWidth = 620;
        const plotHeight = 240;
        const axisSpace = 44;
        const yLabelWidth = 56;
        const width = CHART_PAD + yLabelWidth + plotWidth + CHART_PAD;
        const height = CHART_HEAD + plotHeight + axisSpace;
        const parts = [];
        chartSvgOpen(parts, width, height, meta.title, meta.subtitle);
        const plotX = CHART_PAD + yLabelWidth;
        const baselineY = CHART_HEAD + plotHeight;
        const low = values[0];
        const high = values[values.length - 1];
        const xOf = value => (high > low
            ? plotX + (((value - low) / (high - low)) * plotWidth)
            : plotX + (plotWidth / 2));
        const yOf = fraction => baselineY - (fraction * plotHeight);

        [0, 0.25, 0.5, 0.75, 1].forEach(fraction => {
            const y = yOf(fraction);
            parts.push('<line x1="' + plotX + '" y1="' + y.toFixed(2) + '" x2="' + (plotX + plotWidth) +
                '" y2="' + y.toFixed(2) + '" stroke="#d8dce2" stroke-width="1"/>');
            parts.push('<text x="' + (plotX - 8) + '" y="' + y.toFixed(2) +
                '" font-size="' + CHART_LABEL_FONT + '" fill="#5c6a70" text-anchor="end" dominant-baseline="central">' +
                escapeXml((fraction * 100).toFixed(0) + '%') + '</text>');
        });

        // Sample when there are more values than the cap: adjacent steps would
        // land on the same pixel anyway, and the full list would emit a path an
        // SVG editor struggles to open.
        const stride = Math.max(1, Math.ceil(values.length / COLUMN_CHART_ECDF_MAX_POINTS));
        const path = ['M ' + xOf(low).toFixed(2) + ' ' + yOf(0).toFixed(2)];
        for (let index = 0; index < values.length; index += stride) {
            const fraction = (index + 1) / values.length;
            path.push('L ' + xOf(values[index]).toFixed(2) + ' ' + yOf(fraction - (1 / values.length)).toFixed(2));
            path.push('L ' + xOf(values[index]).toFixed(2) + ' ' + yOf(fraction).toFixed(2));
        }
        // Always land on the last value at 100%, whatever the stride skipped.
        path.push('L ' + xOf(high).toFixed(2) + ' ' + yOf(1).toFixed(2));
        parts.push('<path d="' + path.join(' ') + '" fill="none" stroke="#c8553d" stroke-width="2"/>');

        parts.push('<line x1="' + plotX + '" y1="' + baselineY + '" x2="' + (plotX + plotWidth) +
            '" y2="' + baselineY + '" stroke="#1e2a2f" stroke-width="1.25"/>');
        const tickFont = 12;
        let lastRight = -Infinity;
        chartAxisTicks(low, high, 6).forEach(value => {
            const x = xOf(value);
            const text = formatValue(value);
            const textWidth = estimateTextWidth(text, tickFont);
            if (x - (textWidth / 2) < lastRight + 6) {
                return;
            }
            parts.push('<line x1="' + x.toFixed(2) + '" y1="' + baselineY + '" x2="' + x.toFixed(2) +
                '" y2="' + (baselineY + 5) + '" stroke="#1e2a2f" stroke-width="1"/>');
            parts.push('<text x="' + x.toFixed(2) + '" y="' + (baselineY + 5 + tickFont + 2) +
                '" font-size="' + tickFont + '" text-anchor="middle">' + escapeXml(text) + '</text>');
            lastRight = x + (textWidth / 2);
        });
        parts.push('<text x="' + (plotX + (plotWidth / 2)) + '" y="' + (height - 4) +
            '" font-size="' + tickFont + '" fill="#5c6a70" text-anchor="middle">' +
            escapeXml(meta.columnName) + '</text>');
        return {svg: chartSvgClose(parts), width, height};
    }

    // Table renderer shared by the frequency and summary exports. `columns` is
    // [{key, label, align, swatch}]; a swatch column draws each row's color
    // instead of text.
    function buildTableSVG(rows, columns, meta) {
        const rowHeight = 24;
        const cellPad = 10;
        const swatchSize = 14;
        const widths = columns.map(column => {
            if (column.swatch) {
                return swatchSize + (cellPad * 2);
            }
            const widest = rows.reduce(
                (most, row) => Math.max(most, estimateTextWidth(row[column.key], CHART_LABEL_FONT)),
                estimateTextWidth(column.label, CHART_LABEL_FONT)
            );
            return Math.min(320, widest + (cellPad * 2));
        });
        const width = CHART_PAD + widths.reduce((sum, value) => sum + value, 0) + CHART_PAD;
        const height = CHART_HEAD + ((rows.length + 1) * rowHeight) + CHART_PAD;
        const parts = [];
        chartSvgOpen(parts, width, height, meta.title, meta.subtitle);
        const offsets = [];
        let cursor = CHART_PAD;
        widths.forEach(value => {
            offsets.push(cursor);
            cursor += value;
        });
        const cellText = (text, x, y, align, weight) => {
            const anchor = align === 'right' ? 'end' : 'start';
            const textX = align === 'right' ? (x - cellPad) : (x + cellPad);
            parts.push('<text x="' + textX + '" y="' + y + '" font-size="' + CHART_LABEL_FONT +
                '" text-anchor="' + anchor + '" dominant-baseline="central"' +
                (weight ? ' font-weight="' + weight + '"' : '') + '>' + escapeXml(text) + '</text>');
        };
        columns.forEach((column, columnIndex) => {
            if (!column.label) {
                // The swatch column has no heading, so emit no text node at all.
                return;
            }
            const x = column.align === 'right' ? offsets[columnIndex] + widths[columnIndex] : offsets[columnIndex];
            cellText(column.label, x, CHART_HEAD + (rowHeight / 2), column.align, 'bold');
        });
        parts.push('<line x1="' + CHART_PAD + '" y1="' + (CHART_HEAD + rowHeight) +
            '" x2="' + (width - CHART_PAD) + '" y2="' + (CHART_HEAD + rowHeight) +
            '" stroke="#1e2a2f" stroke-width="1.25"/>');
        rows.forEach((row, rowIndex) => {
            const y = CHART_HEAD + ((rowIndex + 1.5) * rowHeight);
            columns.forEach((column, columnIndex) => {
                if (column.swatch) {
                    parts.push('<rect x="' + (offsets[columnIndex] + cellPad) + '" y="' + (y - (swatchSize / 2)) +
                        '" width="' + swatchSize + '" height="' + swatchSize + '" ' + chartSvgFill(row.color) +
                        ' stroke="#1e2a2f" stroke-width="0.75"/>');
                    return;
                }
                const x = column.align === 'right' ? offsets[columnIndex] + widths[columnIndex] : offsets[columnIndex];
                cellText(row[column.key], x, y, column.align, null);
            });
            parts.push('<line x1="' + CHART_PAD + '" y1="' + (y + (rowHeight / 2)) +
                '" x2="' + (width - CHART_PAD) + '" y2="' + (y + (rowHeight / 2)) +
                '" stroke="#d8dce2" stroke-width="0.75"/>');
        });
        return {svg: chartSvgClose(parts), width, height};
    }

    // ---- HTML table previews ----
    // The two table kinds preview as real HTML rather than as their export SVG
    // so the numbers can be selected and copied out of the dialog.

    function frequencyTableHTML(rows) {
        // cssToHex, not just htmlEscape: this color lands inside a style
        // attribute, and a hex literal cannot carry extra CSS declarations with
        // it however the palette it came from was populated.
        const body = rows.map(row =>
            '<tr><td class="cc-swatch-cell"><span class="cc-swatch" style="background:' +
            cssToHex(row.color) + '"></span></td>' +
            '<td>' + htmlEscape(row.label) + '</td>' +
            '<td class="cc-num">' + htmlEscape(row.count.toLocaleString()) + '</td>' +
            '<td class="cc-num">' + htmlEscape(row.percent.toFixed(1) + '%') + '</td>' +
            '<td class="cc-num">' + htmlEscape(row.globalCount.toLocaleString()) + '</td>' +
            '<td class="cc-num">' + htmlEscape(row.globalPercent.toFixed(1) + '%') + '</td></tr>'
        ).join('');
        // "Percent" is read against the charted rows; "% of all" against that
        // value's whole population, whose size is the column beside it so the
        // denominator is never a mystery.
        return '<table class="cc-freq-table"><thead><tr><th class="cc-swatch-cell"></th>' +
            '<th>Value</th><th class="cc-num">Count</th><th class="cc-num">Percent</th>' +
            '<th class="cc-num">All nodes</th><th class="cc-num">% of all</th>' +
            '</tr></thead><tbody>' + body + '</tbody></table>';
    }

    function summaryTableHTML(rows) {
        const body = rows.map(row =>
            '<tr><td>' + htmlEscape(row[0]) + '</td>' +
            '<td class="cc-num">' + htmlEscape(row[1]) + '</td></tr>'
        ).join('');
        return '<table class="cc-freq-table"><thead><tr><th>Statistic</th>' +
            '<th class="cc-num">Value</th></tr></thead><tbody>' + body + '</tbody></table>';
    }

    function chartTSV(header, rows) {
        return [header.join('\t')].concat(rows.map(row => row.join('\t'))).join('\n') + '\n';
    }

    // ---- The one dispatcher ----

    // Returns {kind, columnName, svg, width, height, html, tsv, notes, empty}.
    // `html` is what the preview shows when set, otherwise `svg` is; `svg` is
    // always present so every kind can be exported as an image.
    function buildColumnChartArtifact(columnName, kind, nodeIndices) {
        const subtitle = columnChartScopeLabel(nodeIndices.length);
        const meta = {
            columnName,
            title: columnName + ' — ' + columnChartKindLabel(kind).toLowerCase(),
            subtitle,
        };
        const notes = [subtitle];
        const empty = {kind, columnName, svg: null, html: null, tsv: '', notes, empty: true};

        if (kind === 'histogram' || kind === 'box' || kind === 'ecdf') {
            const numeric = columnNumericValues(columnName, nodeIndices);
            if (numeric.values.length === 0) {
                notes.push('No numeric values in these rows.');
                return empty;
            }
            if (numeric.missing > 0) {
                notes.push(numeric.missing.toLocaleString() + ' with no value (not plotted)');
            }
            if (kind === 'histogram') {
                const histogram = scopedColumnHistogram(columnName, nodeIndices);
                if (!histogram) {
                    return empty;
                }
                notes.push(histogram.total.toLocaleString() +
                    ' value' + (histogram.total === 1 ? '' : 's') +
                    ' in ' + histogram.binCount + ' bin' + (histogram.binCount === 1 ? '' : 's'));
                const built = buildHistogramSVG(histogram, meta);
                // One bin per integer means the bins *are* values, and the
                // edges the binning uses to separate them (min + k*span/bins)
                // are an implementation detail nobody wants in a TSV.
                const perInteger = histogram.integer
                    && histogram.binCount === (histogram.max - histogram.min) + 1;
                const binSpan = (histogram.highEdge - histogram.lowEdge) / histogram.binCount;
                const tsv = perInteger
                    ? chartTSV(['value', 'count'], histogram.counts.map((count, index) => [
                        histogram.min + index, count,
                    ]))
                    : chartTSV(['bin_low', 'bin_high', 'count'], histogram.counts.map((count, index) => [
                        histogram.lowEdge + (index * binSpan),
                        histogram.lowEdge + ((index + 1) * binSpan),
                        count,
                    ]));
                return Object.assign({kind, columnName, html: null, notes, empty: false}, built, {tsv});
            }
            if (kind === 'box') {
                const summary = numericSummary(numeric.values);
                const built = buildBoxPlotSVG(summary, meta);
                if (summary.outliers.length > 0) {
                    notes.push(summary.outliers.length.toLocaleString() +
                        ' outlier' + (summary.outliers.length === 1 ? '' : 's') + ' beyond 1.5×IQR');
                }
                return Object.assign({kind, columnName, html: null, notes, empty: false}, built, {
                    tsv: chartTSV(['statistic', 'value'], [
                        ['count', summary.count],
                        ['minimum', summary.min],
                        ['low_whisker', summary.lowWhisker],
                        ['q1', summary.q1],
                        ['median', summary.median],
                        ['q3', summary.q3],
                        ['high_whisker', summary.highWhisker],
                        ['maximum', summary.max],
                        ['iqr', summary.iqr],
                        ['outliers', summary.outliers.join(';')],
                    ]),
                });
            }
            const built = buildEcdfSVG(numeric.values, meta);
            return Object.assign({kind, columnName, html: null, notes, empty: false}, built, {
                tsv: chartTSV(['value', 'cumulative_fraction'], numeric.values.map((value, index) => [
                    value,
                    (index + 1) / numeric.values.length,
                ])),
            });
        }

        if (kind === 'summary') {
            const rows = columnSummaryRows(columnName, nodeIndices);
            const built = buildTableSVG(
                rows.map(row => ({statistic: row[0], value: row[1]})),
                [{key: 'statistic', label: 'Statistic'}, {key: 'value', label: 'Value', align: 'right'}],
                meta
            );
            return Object.assign({kind, columnName, notes, empty: false}, built, {
                html: summaryTableHTML(rows),
                tsv: chartTSV(['statistic', 'value'], rows),
            });
        }

        const topN = kind === 'pie'
            ? COLUMN_CHART_PIE_MAX_SLICES
            : (kind === 'bar' ? COLUMN_CHART_BAR_MAX_BARS : Infinity);
        const model = columnValueDistribution(columnName, nodeIndices, {topN});
        if (model.entries.length === 0) {
            notes.push('No values in these rows.');
            return empty;
        }
        notes.push(model.distinctTotal.toLocaleString() + ' distinct value' +
            (model.distinctTotal === 1 ? '' : 's'));
        if (model.truncated) {
            notes.push('Showing the ' + topN + ' most common; ' + model.otherDistinct.toLocaleString() +
                ' more rolled into “Other”');
        }
        if (model.nullCount > 0) {
            notes.push(model.nullCount.toLocaleString() + ' with no value (shown as —)');
        }
        const rows = frequencyTableRows(model);
        const tsv = chartTSV(
            ['value', 'count', 'percent', 'count_all_nodes', 'percent_of_all', 'color'],
            rows.map(row => [
                row.label, row.count, row.percent.toFixed(4),
                row.globalCount, row.globalPercent.toFixed(4), cssToHex(row.color),
            ])
        );
        if (kind === 'bar') {
            return Object.assign({kind, columnName, html: null, notes, empty: false},
                buildBarChartSVG(model, meta), {tsv});
        }
        if (kind === 'pie') {
            return Object.assign({kind, columnName, html: null, notes, empty: false},
                buildPieChartSVG(model, meta), {tsv});
        }
        // Frequency table. The preview caps its row count (it goes through
        // innerHTML); the TSV above always carries every value.
        //
        // The table is the one kind that shows both denominators, so it is the one
        // that has to say what the second is. Only worth spelling out once the scope
        // is narrower than the bundle: when it is not, every "% of all" reads 100%.
        if (model.rowCount < model.globalTotal) {
            notes.push('“% of all” counts each value against its own network-wide total, ' +
                'not against these ' + model.rowCount.toLocaleString() + ' rows');
        }
        const previewRows = rows.slice(0, COLUMN_CHART_TABLE_MAX_ROWS);
        if (rows.length > previewRows.length) {
            notes.push('Showing the first ' + COLUMN_CHART_TABLE_MAX_ROWS.toLocaleString() +
                ' rows; the TSV has all ' + rows.length.toLocaleString() + '.');
        }
        const built = buildTableSVG(
            previewRows.map(row => ({
                color: row.color,
                label: row.label,
                count: row.count.toLocaleString(),
                percent: row.percent.toFixed(1) + '%',
                globalCount: row.globalCount.toLocaleString(),
                globalPercent: row.globalPercent.toFixed(1) + '%',
            })),
            [
                {key: 'color', label: '', swatch: true},
                {key: 'label', label: 'Value'},
                {key: 'count', label: 'Count', align: 'right'},
                {key: 'percent', label: 'Percent', align: 'right'},
                {key: 'globalCount', label: 'All nodes', align: 'right'},
                {key: 'globalPercent', label: '% of all', align: 'right'},
            ],
            meta
        );
        return Object.assign({kind, columnName, notes, empty: false}, built, {
            html: frequencyTableHTML(previewRows),
            tsv,
        });
    }

    // ---- Anchored menu ----

    function columnChartMenuElement() {
        return document.getElementById('column-chart-menu');
    }

    function columnChartMenuIsOpen() {
        const menu = columnChartMenuElement();
        return Boolean(menu) && !menu.hidden;
    }

    function closeColumnChartMenu() {
        const menu = columnChartMenuElement();
        if (!menu || menu.hidden) {
            return;
        }
        menu.hidden = true;
        menu.innerHTML = '';
        const anchor = state.columnChartMenu && state.columnChartMenu.anchor;
        if (anchor && anchor.isConnected) {
            anchor.setAttribute('aria-expanded', 'false');
        }
        state.columnChartMenu = null;
    }

    function openColumnChartMenu(columnName, anchorButton) {
        const menu = columnChartMenuElement();
        if (!menu || !state.bundle) {
            return;
        }
        const wasOpenFor = state.columnChartMenu && state.columnChartMenu.columnName;
        closeColumnChartMenu();
        if (wasOpenFor === columnName) {
            // A second click on the same glyph closes the menu.
            return;
        }
        menu.innerHTML = columnChartKindsFor(columnName).map(kind =>
            '<button type="button" role="menuitem" data-chart-kind="' + htmlEscape(kind.id) + '">' +
            htmlEscape(kind.label) + '</button>'
        ).join('');
        menu.hidden = false;
        anchorButton.setAttribute('aria-expanded', 'true');
        // Positioned after unhiding so offsetWidth/Height are real, and flipped
        // rather than clipped when the header sits near a viewport edge.
        const anchorRect = anchorButton.getBoundingClientRect();
        const menuWidth = menu.offsetWidth;
        const menuHeight = menu.offsetHeight;
        let left = anchorRect.left;
        if (left + menuWidth > window.innerWidth - 8) {
            left = Math.max(8, anchorRect.right - menuWidth);
        }
        let top = anchorRect.bottom + 4;
        if (top + menuHeight > window.innerHeight - 8) {
            top = Math.max(8, anchorRect.top - menuHeight - 4);
        }
        menu.style.left = left + 'px';
        menu.style.top = top + 'px';
        // anchorLeft/anchorTop record where the anchor sat when we positioned
        // against it, so a later scroll or resize can tell "the anchor moved"
        // from "an event arrived".
        state.columnChartMenu = {
            columnName,
            anchor: anchorButton,
            anchorLeft: anchorRect.left,
            anchorTop: anchorRect.top,
        };
        const first = menu.querySelector('button');
        if (first) {
            // preventScroll matters: the menu is already positioned inside the
            // viewport, and letting focus() scroll an ancestor to "reveal" it
            // would fire the scroll listener below and close the menu again.
            first.focus({preventScroll: true});
        }
    }

    // The menu is anchored in viewport coordinates, so it is only stale once its
    // anchor has actually moved. Closing on the scroll event itself is not
    // enough of a test: focusing the glyph can start a scroll whose event
    // arrives *after* the menu opens, and that would close it immediately.
    function closeColumnChartMenuIfAnchorMoved() {
        const open = state.columnChartMenu;
        if (!columnChartMenuIsOpen() || !open) {
            return;
        }
        if (!open.anchor || !open.anchor.isConnected) {
            closeColumnChartMenu();
            return;
        }
        const rect = open.anchor.getBoundingClientRect();
        if (Math.abs(rect.left - open.anchorLeft) > 1 || Math.abs(rect.top - open.anchorTop) > 1) {
            closeColumnChartMenu();
        }
    }

    function moveColumnChartMenuFocus(delta) {
        const menu = columnChartMenuElement();
        if (!menu || menu.hidden) {
            return;
        }
        const items = Array.from(menu.querySelectorAll('button'));
        if (items.length === 0) {
            return;
        }
        const current = items.indexOf(document.activeElement);
        let next;
        if (delta === 'first') {
            next = 0;
        } else if (delta === 'last') {
            next = items.length - 1;
        } else {
            next = ((current < 0 ? 0 : current) + delta + items.length) % items.length;
        }
        // preventScroll for the same reason openColumnChartMenu uses it: a
        // scroll here would move the anchor and close the menu mid-navigation.
        items[next].focus({preventScroll: true});
    }

    // ---- Dialog ----

    function columnChartIsOpen() {
        const overlay = document.getElementById('column-chart-overlay');
        return Boolean(overlay) && !overlay.hidden;
    }

    // `returnFocusTo` is the header glyph the dialog was opened from, so
    // closing hands focus back where it started. It may have been detached in
    // the meantime -- updateMetadataTable rebuilds the whole header row.
    function openColumnChart(columnName, kind, returnFocusTo) {
        if (!state.bundle) {
            return;
        }
        const kinds = columnChartKindsFor(columnName);
        const resolved = kinds.some(entry => entry.id === kind) ? kind : kinds[0].id;
        state.columnChart = {columnName, kind: resolved, returnFocusTo: returnFocusTo || null};
        const select = document.getElementById('column-chart-kind');
        select.innerHTML = kinds.map(entry =>
            '<option value="' + htmlEscape(entry.id) + '"' +
            (entry.id === resolved ? ' selected' : '') + '>' + htmlEscape(entry.label) + '</option>'
        ).join('');
        document.getElementById('column-chart-overlay').hidden = false;
        renderColumnChart();
        document.getElementById('column-chart-close').focus();
    }

    function closeColumnChart() {
        const overlay = document.getElementById('column-chart-overlay');
        if (!overlay || overlay.hidden) {
            return;
        }
        overlay.hidden = true;
        document.getElementById('column-chart-preview').innerHTML = '';
        const returnFocusTo = state.columnChart && state.columnChart.returnFocusTo;
        state.columnChart = null;
        state.columnChartArtifact = null;
        if (returnFocusTo && returnFocusTo.isConnected) {
            returnFocusTo.focus({preventScroll: true});
        }
    }

    function renderColumnChart() {
        if (!state.columnChart || !state.bundle) {
            return;
        }
        const {columnName, kind} = state.columnChart;
        // A column can be renamed or deleted from under an open dialog.
        if (!state.metadataColumnIndexByName.has(columnName)) {
            closeColumnChart();
            return;
        }
        const artifact = buildColumnChartArtifact(columnName, kind, columnChartNodeIndices());
        state.columnChartArtifact = artifact;
        document.getElementById('column-chart-title').textContent =
            columnName + ' — ' + columnChartKindLabel(kind);
        const preview = document.getElementById('column-chart-preview');
        // notes[0] is always the scope line; anything after it on an empty
        // artifact is the reason there is nothing to draw, and belongs in the
        // preview box rather than leaving it blank.
        preview.innerHTML = artifact.html || artifact.svg
            || ('<p class="note">' + htmlEscape(artifact.notes.slice(1).join(' ')) + '</p>');
        const notes = artifact.notes.slice();
        if (!artifact.empty) {
            notes.push('PNG export uses the resolution set in the view toolbar.');
        }
        document.getElementById('column-chart-note').textContent = notes.join(' · ');
        ['column-chart-copy-tsv', 'column-chart-download-tsv',
            'column-chart-export-svg', 'column-chart-export-png'].forEach(id => {
            document.getElementById(id).disabled = artifact.empty;
        });
    }

    function refreshColumnChartIfOpen() {
        if (columnChartIsOpen()) {
            renderColumnChart();
        }
    }

    // ---- Exports ----

    function columnChartFileBase() {
        const artifact = state.columnChartArtifact;
        if (!artifact) {
            return exportBaseName();
        }
        return exportBaseName() + '_' + artifact.columnName.replace(/\s+/g, '_') + '_' + artifact.kind;
    }

    function exportColumnChartSVG() {
        const artifact = state.columnChartArtifact;
        if (!artifact || !artifact.svg) {
            return;
        }
        // The prolog is added only here: chartSvgOpen leaves it out so the same
        // string can be assigned to innerHTML for the preview.
        const svg = '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n' + artifact.svg;
        const href = URL.createObjectURL(new Blob([svg], {type: 'image/svg+xml'}));
        triggerDownload(href, columnChartFileBase() + '.svg', true);
    }

    function exportColumnChartPNG() {
        const artifact = state.columnChartArtifact;
        if (!artifact || !artifact.svg) {
            return;
        }
        const scaleFactor = selectedPngScale();
        const svgUrl = URL.createObjectURL(new Blob([artifact.svg], {type: 'image/svg+xml'}));
        const image = new Image();
        image.onload = () => {
            const target = document.createElement('canvas');
            target.width = Math.max(1, Math.round(artifact.width * scaleFactor));
            target.height = Math.max(1, Math.round(artifact.height * scaleFactor));
            const targetContext = target.getContext('2d');
            if (!targetContext) {
                URL.revokeObjectURL(svgUrl);
                window.alert('Could not export the chart PNG at ' + scaleFactor + '×; try a lower resolution.');
                return;
            }
            targetContext.drawImage(image, 0, 0, target.width, target.height);
            target.toBlob(blob => {
                URL.revokeObjectURL(svgUrl);
                if (!blob) {
                    window.alert('The chart is too large to export as a ' + scaleFactor + '× PNG.');
                    return;
                }
                const href = URL.createObjectURL(blob);
                const suffix = scaleFactor > 1 ? '@' + scaleFactor + 'x.png' : '.png';
                triggerDownload(href, columnChartFileBase() + suffix, true);
            }, 'image/png');
        };
        image.onerror = () => {
            URL.revokeObjectURL(svgUrl);
            window.alert('Could not rasterize the chart for PNG export.');
        };
        image.src = svgUrl;
    }

    function columnChartTSV() {
        return state.columnChartArtifact ? state.columnChartArtifact.tsv : '';
    }

    async function copyColumnChartTSV() {
        const tsv = columnChartTSV();
        if (!tsv) {
            return;
        }
        try {
            await writeTextToClipboard(tsv);
            setStatus('Copied the ' + columnChartKindLabel(state.columnChartArtifact.kind).toLowerCase() +
                ' for "' + state.columnChartArtifact.columnName + '".');
        } catch (error) {
            console.error(error);
            setStatus('Failed to copy the chart data: ' + error.message);
        }
    }

    function downloadColumnChartTSV() {
        const tsv = columnChartTSV();
        if (!tsv) {
            return;
        }
        const href = URL.createObjectURL(new Blob([tsv], {type: 'text/tab-separated-values'}));
        triggerDownload(href, columnChartFileBase() + '.tsv', true);
    }

    // ---- Wiring ----

    function setupColumnCharts() {
        const menu = columnChartMenuElement();
        const overlay = document.getElementById('column-chart-overlay');

        menu.addEventListener('click', event => {
            const item = event.target.closest('[data-chart-kind]');
            if (!item || !state.columnChartMenu) {
                return;
            }
            const {columnName, anchor} = state.columnChartMenu;
            closeColumnChartMenu();
            openColumnChart(columnName, item.dataset.chartKind, anchor);
        });

        document.addEventListener('pointerdown', event => {
            if (!columnChartMenuIsOpen()) {
                return;
            }
            if (event.target.closest('#column-chart-menu') || event.target.closest('[data-chart-column]')) {
                return;
            }
            closeColumnChartMenu();
        });
        // The scroll listener has to capture: the metadata table scrolls inside
        // .table-wrap, whose scroll events never reach window by bubbling.
        window.addEventListener('resize', closeColumnChartMenuIfAnchorMoved);
        window.addEventListener('scroll', closeColumnChartMenuIfAnchorMoved, true);

        document.addEventListener('keydown', event => {
            if (columnChartMenuIsOpen()) {
                if (event.key === 'Escape') {
                    const anchor = state.columnChartMenu && state.columnChartMenu.anchor;
                    event.preventDefault();
                    closeColumnChartMenu();
                    if (anchor && anchor.isConnected) {
                        anchor.focus({preventScroll: true});
                    }
                    return;
                }
                if (event.key === 'ArrowDown' || event.key === 'ArrowUp'
                    || event.key === 'Home' || event.key === 'End') {
                    event.preventDefault();
                    moveColumnChartMenuFocus(
                        event.key === 'ArrowDown' ? 1
                            : (event.key === 'ArrowUp' ? -1 : (event.key === 'Home' ? 'first' : 'last'))
                    );
                }
                return;
            }
            if (event.key === 'Escape' && columnChartIsOpen()) {
                event.preventDefault();
                closeColumnChart();
            }
        });

        overlay.addEventListener('click', event => {
            if (event.target === event.currentTarget) {
                closeColumnChart();
            }
        });
        document.getElementById('column-chart-close').addEventListener('click', closeColumnChart);
        document.getElementById('column-chart-kind').addEventListener('change', event => {
            if (state.columnChart) {
                state.columnChart.kind = event.target.value;
                renderColumnChart();
            }
        });
        document.getElementById('column-chart-copy-tsv').addEventListener('click', copyColumnChartTSV);
        document.getElementById('column-chart-download-tsv').addEventListener('click', downloadColumnChartTSV);
        document.getElementById('column-chart-export-svg').addEventListener('click', exportColumnChartSVG);
        document.getElementById('column-chart-export-png').addEventListener('click', exportColumnChartPNG);
    }
