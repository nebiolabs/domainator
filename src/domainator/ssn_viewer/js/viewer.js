    // Values the page is built with, from the dsnv-config block (see ssn_viewer_html.py).
    const DSNV_CONFIG = JSON.parse(document.getElementById('dsnv-config').textContent);
    const state = {
        bundle: null,
        bundleVersionWarning: '',
        metadataColumns: [],
        metadataRows: [],
        metadataByNodeIndex: [],
        metadataColumnByName: new Map(),
        metadataColumnIndexByName: new Map(),
        metadataColorInfoByName: new Map(),
        customPalettes: {},
        // Derived cache of defaultCategoricalPalette() results, keyed like
        // customPalettes. Not user data, so it is never saved into a session;
        // rebuildMetadataCaches() drops it.
        defaultPalettes: new Map(),
        // Names of numeric metadata columns the viewer should color as discrete
        // categories instead of a gradient (e.g. integer cluster numbers).
        categoricalColumns: new Set(),
        colorPickerPage: 0,
        // Where drawSplitChart() last painted each mark, for hover hit-testing.
        splitChartHit: null,
        // The stop the slider was last moved to deliberately; see currentSliderStop().
        selectedStop: null,
        // The stop in effect before the current split-chart click gesture began.
        stopBeforeSplitChartClick: null,
        // {min, max} in threshold units while the split chart is zoomed, else null.
        splitChartZoom: null,
        // {pointerId, startX, lastX, moved} while the split chart is being panned.
        splitChartDrag: null,
        splitChartSuppressClick: false,
        // {columnName, anchor} while a column header's chart menu is open.
        columnChartMenu: null,
        // {columnName, kind} while the chart dialog is open, and the built
        // artifact the four export buttons hand out (so they never rebuild it).
        columnChart: null,
        columnChartArtifact: null,
        // Binned counts for the gradient dialog's histogram. Cached because binning walks
        // every node, while recoloring the bars happens on every keystroke in the dialog.
        colorHistogram: null,
        // Working list of {value, color} gradient stops while the picker is open.
        gradientStops: [],
        metadataColumnWidths: new Map(),
        metadataSearchTextByNodeIndex: [],
        metadataPage: 0,
        activeClusters: [],
        visibleClusters: [],
        selectedNodeIndices: new Set(),
        selectedMetadataNodeIndices: new Set(),
        selectionPresets: new Map(),
        presetPreviewSlot: null,
        metadataEditCell: null,
        pendingNameAction: null,
        nodeIndexById: null,
        pendingRestoreViewTransform: null,
        sessionMissingNodeIds: 0,
        metadataRowSelectionAnchor: null,
        metadataResize: null,
        visibleLayout: [],
        splitLinks: [],
        sliderModel: null,
        metadataSort: {columnKey: null, direction: 'asc'},
        viewTransform: {scale: 1, offsetX: 0, offsetY: 0, minScale: 0.02, maxScale: 10},
        selectionBox: null,
        dragState: null,
        suppressClick: false,
        layoutCache: new Map(),
        dotLayoutCache: new Map(),
        pendingThresholdUIFrame: null,
        pendingThresholdUIResetView: false,
        pendingMetadataFilterTimer: null,
        pendingClusterRenderFrame: null,
        allNodeIndices: [],
        nodeColorCache: [],
        metadataBaseNoteText: '',
        renderedNodeIndices: [],
        layoutWorker: null,
        pendingLayoutRequestId: 0,
        _pendingVisibleGraph: null,
        _pendingLayoutAlgorithm: null,
        _pendingLayoutResetView: true,
        layoutComputing: false,
    };

    const splitCanvas = document.getElementById('split-chart');
    const splitContext = splitCanvas.getContext('2d');
    const clusterCanvas = document.getElementById('cluster-view');
    const clusterContext = clusterCanvas.getContext('2d');
    const DEFAULT_METADATA_PAGE_SIZE = 250;
    const EMBEDDED_BUNDLE_BASE64 = DSNV_CONFIG.embeddedBundleBase64;
    const LAYOUT_WORKER_CODE = DSNV_CONFIG.layoutWorkerCode;
    const BUNDLE_FORMAT = DSNV_CONFIG.bundleFormat;
    const BUNDLE_VERSION = DSNV_CONFIG.bundleVersion;
    const SUPPORTED_BUNDLE_VERSIONS = DSNV_CONFIG.supportedBundleVersions;
    const DOMAINATOR_VERSION = DSNV_CONFIG.domainatorVersion;
    const VIEWER_APP_NAME = DSNV_CONFIG.viewerAppName;

    // The shell viewer is built before any bundle exists, and the file picker can
    // swap bundles at any time, so the heading and tab follow the loaded bundle.
    function updateViewerTitle(networkName) {
        const name = (networkName || '').trim();
        document.getElementById('viewer-title').textContent =
            name === '' ? VIEWER_APP_NAME : VIEWER_APP_NAME + ': ' + name;
        if (name !== '') {
            document.title = name;
        }
    }

    function setStatus(message) {
        document.getElementById('bundle-status').textContent = message;
    }

    async function writeTextToClipboard(text) {
        if (navigator.clipboard?.writeText) {
            await navigator.clipboard.writeText(text);
            return;
        }
        const helper = document.createElement('textarea');
        helper.value = text;
        helper.setAttribute('readonly', 'readonly');
        helper.style.position = 'fixed';
        helper.style.opacity = '0';
        document.body.appendChild(helper);
        helper.select();
        document.execCommand('copy');
        document.body.removeChild(helper);
    }


    function browserSupportsBundleLoading() {
        return 'DecompressionStream' in window;
    }

    function warnIfUnsupported() {
        if (!browserSupportsBundleLoading()) {
            document.getElementById('browser-warning').textContent = EMBEDDED_BUNDLE_BASE64
                ? 'This browser lacks DecompressionStream support; bundled data will still open, but loading external gzipped bundles is unavailable.'
                : 'This browser lacks DecompressionStream support; gzipped bundles cannot be loaded.';
        }
    }

    async function decodeBundleFile(file) {
        if (browserSupportsBundleLoading()) {
            try {
                const stream = file.stream().pipeThrough(new DecompressionStream('gzip'));
                const text = await new Response(stream).text();
                return JSON.parse(text);
            } catch (error) {
                console.warn('gzip decode failed, trying plain JSON', error);
            }
        }
        const text = await file.text();
        return JSON.parse(text);
    }

    function base64ToBytes(base64Value) {
        return Uint8Array.from(atob(base64Value), c => c.charCodeAt(0));
    }

    async function autoloadEmbeddedBundle() {
        if (!EMBEDDED_BUNDLE_BASE64) {
            return;
        }
        setStatus('Loading bundled data...');
        try {
            const bytes = base64ToBytes(EMBEDDED_BUNDLE_BASE64);
            let text;
            if (browserSupportsBundleLoading()) {
                const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));
                text = await new Response(stream).text();
            } else {
                text = new TextDecoder().decode(bytes);
            }
            installBundle(JSON.parse(text));
        } catch (error) {
            console.error(error);
            setStatus('Failed to load bundled data: ' + error.message);
        }
    }

    function formatValue(value) {
        if (value === null || value === undefined || value === '') {
            return '—';
        }
        if (typeof value === 'number') {
            return Number.isInteger(value) ? value.toLocaleString() : value.toFixed(2);
        }
        return String(value);
    }

    function nodeId(nodeIndex) {
        return state.bundle.graph.nodes[nodeIndex];
    }

    function installBundle(bundle) {
        if (bundle.format !== BUNDLE_FORMAT) {
            throw new Error('Unsupported bundle format: ' + bundle.format);
        }
        // A wrong `format` is fatal, but an unrecognized `version` is not: bundle
        // revisions so far have only added sections, so loading and reporting beats
        // refusing a file this viewer can very likely still render.
        state.bundleVersionWarning = SUPPORTED_BUNDLE_VERSIONS.includes(bundle.version)
            ? ''
            : ' (bundle version ' + bundle.version + ' is not one this viewer knows about ('
                + SUPPORTED_BUNDLE_VERSIONS.join(', ') + '); loading anyway, some features may misbehave)';
        state.bundle = bundle;
        updateViewerTitle(bundle.name);
        state.metadataColumns = bundle.metadata.columns || [];
        state.metadataRows = bundle.metadata.rows || [];
        state.metadataByNodeIndex = state.metadataRows;
        state.customPalettes = {};
        state.categoricalColumns = new Set(bundle.defaults.categorical_columns || []);
        state.metadataSort = {columnKey: null, direction: 'asc'};
        state.metadataColumnWidths = new Map();
        state.metadataPage = 0;
        state.metadataResize = null;
        state.selectedNodeIndices = new Set();
        state.selectedMetadataNodeIndices = new Set();
        state.selectionPresets = new Map();
        state.metadataRowSelectionAnchor = null;
        state.nodeIndexById = null;
        state.pendingRestoreViewTransform = null;
        state.dotLayoutCache = new Map();
        state.layoutCache = new Map();
        state.latticeGlobal = null;
        // The previous bundle's drawing indexes its own hierarchy, and a force layout
        // is computed off the main thread, so the canvas repaints at least once before
        // this bundle has a layout: drop the old one rather than draw its cluster ids
        // against the new hierarchy. Bumping the request id discards a worker reply
        // still in flight for the old bundle.
        state.visibleLayout = [];
        state.splitLinks = [];
        state.latticeNodeItem = null;
        state.activeClusters = [];
        state.visibleClusters = [];
        state.pendingLayoutRequestId++;
        // A window onto the previous bundle's threshold range means nothing here; a
        // session restore sets its own afterwards. Cleared before the series is built,
        // because the selection is taken against whatever window is in force.
        state.splitChartZoom = null;
        // Derived here rather than read from the file: see deriveMergeSeries. Every
        // consumer below reads state.series, never the bundle's own copy (a pre-v6 file
        // has one, already capped by whoever built it).
        rebuildMergeSeries(DEFAULT_WINDOWED_MAX_MERGE_EVENTS);
        state.allNodeIndices = bundle.graph.nodes.map((_, i) => i);
        rebuildMetadataCaches();

        // Must precede rebuildNodeColorCache: the cache is keyed off the "Color by"
        // selection, which populateMetadataControls sets from bundle.defaults.
        populateMetadataControls();
        rebuildNodeColorCache();
        const slider = document.getElementById('threshold-slider');
        slider.max = String(state.sliderModel.maxPosition);
        slider.value = String(state.sliderModel.initialPosition);
        slider.disabled = state.sliderModel.stops.length === 0;
        document.getElementById('threshold-min-label').textContent = state.sliderModel.minLabel;
        document.getElementById('threshold-max-label').textContent = state.sliderModel.maxLabel;
        document.getElementById('sort-components-by-size').disabled = false;
        document.getElementById('focus-largest-cluster').disabled = false;
        document.getElementById('reset-view').disabled = false;
        document.getElementById('zoom-level').disabled = false;
        document.getElementById('export-png').disabled = false;
        document.getElementById('export-png-scale').disabled = false;
        document.getElementById('export-svg').disabled = false;
        document.getElementById('export-split-png').disabled = false;
        document.getElementById('export-split-svg').disabled = false;
        document.getElementById('customize-colors').disabled = false;
        document.getElementById('save-session').disabled = false;
        document.getElementById('rename-network').disabled = false;
        state.presetPreviewSlot = null;
        if (colorPickerIsOpen()) {
            openColorPicker();
        }
        document.getElementById('threshold-input').disabled = state.sliderModel.stops.length === 0;
        updateThresholdStepButtons();
        document.getElementById('metadata-select-nodes').disabled = true;
        document.getElementById('metadata-deselect-rows').disabled = true;
        document.getElementById('metadata-reset-sort').disabled = true;
        document.getElementById('metadata-filter').disabled = false;
        document.getElementById('metadata-filter').value = '';
        document.getElementById('metadata-null-order').disabled = false;
        document.getElementById('metadata-null-order').value = 'last';
        document.getElementById('metadata-rows-per-page').disabled = false;
        document.getElementById('metadata-panel-add').disabled = false;
        document.getElementById('metadata-panel-set').disabled = false;
        document.getElementById('metadata-panel-select').disabled = false;
        document.getElementById('metadata-panel-rename').disabled = false;
        document.getElementById('metadata-panel-delete').disabled = false;
        closeMetadataEditPanels();
        populateEditableColumnMenus();
        document.getElementById('metadata-prev-page').disabled = true;
        document.getElementById('metadata-next-page').disabled = true;
        // Must run after the slider's `max` is set (a value above the old max
        // would be clamped) and after populateMetadataControls has built the
        // color-by/label-by options, but before the color cache is used.
        let sessionNote = '';
        if (bundle.app_state) {
            sessionNote = applySessionState(bundle.app_state);
            rebuildNodeColorCache();
            if (colorPickerIsOpen()) {
                openColorPicker();
            }
        }
        setStatus(
            'Loaded ' + (bundle.name || 'bundle') + ' with ' + bundle.graph.nodes.length.toLocaleString() + ' nodes.'
            + state.bundleVersionWarning + sessionNote
        );
        renderPresetSlots();
        updateSplitEventCount();
        updateThresholdUI();
        updateMetadataTable();
    }

    function rebuildMetadataCaches() {
        // Which default color a value gets depends on the whole sorted set of
        // values in its column, so any change to the metadata invalidates it.
        state.defaultPalettes = new Map();
        state.metadataColumnByName = new Map();
        state.metadataColumnIndexByName = new Map();
        state.metadataColorInfoByName = new Map();

        state.metadataColumns.forEach((column, columnIndex) => {
            state.metadataColumnByName.set(column.name, column);
            state.metadataColumnIndexByName.set(column.name, columnIndex);
            if (column.type === 'int' || column.type === 'float') {
                // `baseType` is the column's intrinsic kind; `type` is how it is colored
                // right now, which the categorical toggle can flip for numeric columns.
                state.metadataColorInfoByName.set(column.name, {type: 'numeric', baseType: 'numeric', min: Infinity, max: -Infinity});
                return;
            }
            state.metadataColorInfoByName.set(column.name, {type: 'categorical', baseType: 'categorical'});
        });

        state.metadataByNodeIndex.forEach(row => {
            if (!row) {
                return;
            }
            state.metadataColumns.forEach((column, columnIndex) => {
                if (column.type !== 'int' && column.type !== 'float') {
                    return;
                }
                const value = row[columnIndex];
                if (typeof value !== 'number' || Number.isNaN(value)) {
                    return;
                }
                const info = state.metadataColorInfoByName.get(column.name);
                if (!info) {
                    return;
                }
                info.min = Math.min(info.min, value);
                info.max = Math.max(info.max, value);
            });
        });

        state.metadataColorInfoByName.forEach((info, columnName) => {
            if (info.baseType !== 'numeric') {
                return;
            }
            if (!Number.isFinite(info.min) || !Number.isFinite(info.max)) {
                info.min = 0;
                info.max = 0;
            }
            // min/max stay available so toggling back to the gradient is free.
            info.type = state.categoricalColumns.has(columnName) ? 'categorical' : 'numeric';
        });

        state.metadataSearchTextByNodeIndex = state.bundle.graph.nodes.map((nodeName, nodeIndex) => {
            const row = state.metadataByNodeIndex[nodeIndex] || [];
            const parts = [String(nodeName)];
            row.forEach(value => {
                if (value === null || value === undefined) {
                    return;
                }
                if (Array.isArray(value)) {
                    parts.push(value.join(' '));
                    return;
                }
                parts.push(String(value));
            });
            return parts.join(' ').toLowerCase();
        });
    }

    function rebuildNodeColorCache() {
        if (!state.bundle) {
            state.nodeColorCache = [];
            return;
        }
        const nodeCount = state.bundle.graph.nodes.length;
        const cache = new Array(nodeCount);
        for (let i = 0; i < nodeCount; i++) {
            cache[i] = nodeColor(i);
        }
        state.nodeColorCache = cache;
        // A palette edit that recolors the nodes should recolor an open chart
        // of the same column too.
        refreshColumnChartIfOpen();
    }

    // ---- Slider positions ----
    //
    // Positions are a warp of the threshold axis rather than a straight scale of it,
    // because two things pull on them at once:
    //
    //   * Every stop should be reachable by dragging, and 1001 integer positions have
    //     to carry every stop there is. On a network whose merges bunch into a narrow
    //     band of scores, 401 stops once landed on 103 positions -- eleven sharing one
    //     -- which left 74% of them impossible to select by dragging at all.
    //   * The track should still read like the chart above it, so that dragging toward
    //     a spike on the chart moves toward that spike.
    //
    // Separating crowded stops has to borrow track from somewhere, so the two cannot
    // both hold exactly. What holds instead, in two steps:
    //
    //   1. Every stop is given one position of its own, which is what makes it
    //      draggable, and
    //   2. the track left over is spent saying where the stops are -- proportionally to
    //      threshold, except that the stretch the chart is zoomed into carries the
    //      chart's own magnification times the track density of the rest, so zooming
    //      the chart magnifies that stretch of the slider too.
    //
    // Where stops are spread out, step 1 costs almost nothing and the result is the
    // plain proportional map this started as. Where they crowd, step 1 is the whole
    // budget and they come out evenly spaced -- which is the most a thousand positions
    // can say about a thousand stops.

    // Finite stops share 0..920; the infinity stop sits alone at 1000, far enough above
    // them to read as the separate thing it is.
    const SLIDER_MAX_FINITE_POSITION = 920;
    const SLIDER_INFINITY_POSITION = 1000;

    // Threshold -> fraction of the finite track. Piecewise linear: one density inside
    // the chart's window and another outside it, integrated along the axis, so the map
    // stays monotone and stays proportional *within* each region.
    function sliderValueWarp(lowValue, highValue) {
        const span = highValue - lowValue;
        if (!(span > 0)) {
            return () => 0;
        }
        const proportional = value => (value - lowValue) / span;
        const zoom = state.splitChartZoom;
        if (!zoom) {
            return proportional;
        }
        const windowLow = Math.min(Math.max(zoom.min, lowValue), highValue);
        const windowHigh = Math.max(Math.min(zoom.max, highValue), lowValue);
        const windowSpan = windowHigh - windowLow;
        if (!(windowSpan > 0)) {
            return proportional;
        }
        // The chart's own magnification. Because the outside keeps its own width, the
        // window's share of the track tends to half however far the zoom goes, rather
        // than swallowing the track and stranding everything else on one position.
        const weight = span / windowSpan;
        const total = (span - windowSpan) + (weight * windowSpan);
        return value => {
            const below = Math.min(Math.max(value, lowValue), windowLow) - lowValue;
            const inside = Math.min(Math.max(value, windowLow), windowHigh) - windowLow;
            const above = Math.min(Math.max(value, windowHigh), highValue) - windowHigh;
            return (below + (weight * inside) + above) / total;
        };
    }

    // Assigns `sliderPosition` in place, so a stop object held elsewhere -- notably
    // state.selectedStop -- keeps pointing at the same stop across a remap.
    //
    // Every stop is given a slot of its own first, and only the track left over after
    // that is spent on saying where the stops are, through the warp. That ordering is
    // what makes the two pulls compatible: the reserved slot is what makes a stop
    // draggable, and the surplus is what keeps the sparse stretches proportional and
    // hands the chart's window its extra room. Positions come out strictly increasing
    // by construction -- consecutive stops differ by at least the one slot, and rounding
    // a non-decreasing sequence cannot close a gap of one -- so no separation pass is
    // needed, and none of the cascading that a pass has to do when a crowd overflows
    // the end of the track can undo the warp.
    function positionSliderStops(finiteStops) {
        if (finiteStops.length === 0) {
            return;
        }
        const lastIndex = finiteStops.length - 1;
        const lowValue = finiteStops[0].threshold_value;
        const highValue = finiteStops[lastIndex].threshold_value;
        const surplus = SLIDER_MAX_FINITE_POSITION - lastIndex;
        if (surplus < 0 || !(highValue - lowValue > 0)) {
            // More stops than the track has positions -- only reachable by raising
            // --max_merge_events past ~900 -- or every stop at one threshold. Neither
            // leaves anything to be proportional to, so fall back to rank.
            finiteStops.forEach((stop, stopIndex) => {
                stop.sliderPosition = Math.round(
                    (stopIndex / Math.max(1, lastIndex)) * SLIDER_MAX_FINITE_POSITION);
            });
            return;
        }
        const warp = sliderValueWarp(lowValue, highValue);
        finiteStops.forEach((stop, stopIndex) => {
            stop.sliderPosition = stopIndex + Math.round(surplus * warp(stop.threshold_value));
        });
    }

    // Both which stops exist and where they sit are functions of the chart's window, so
    // both are redone whenever that window moves -- and the thumb put back on the
    // threshold it was already on, which is a different position, and possibly a
    // different stop object, under the new map.
    //
    // The selection is re-run because the cap spends itself on the window: zoom in and
    // the same `maxMergeEvents` buys stops among the events actually on screen, which is
    // what makes small events reachable at all. The threshold in effect is captured
    // first and pinned through the re-selection, so the cut never changes under a zoom.
    function repositionSliderStops() {
        if (!state.series || !state.sliderModel || state.sliderModel.stops.length === 0) {
            return;
        }
        const threshold = selectedThresholdValue();
        applyMergeSelection(threshold);
        restoreThreshold(threshold);
        updateThresholdStepButtons();
        updateSplitEventCount();
    }

    function buildSliderModel(sourceStops) {
        const stops = (sourceStops || []).map(stop => ({...stop}));
        if (stops.length === 0) {
            return {stops: [], maxPosition: 0, initialPosition: 0, minLabel: 'min', maxLabel: '∞'};
        }

        const infinityStop = stops.find(stop => stop.threshold_value === null) || {edge_index: -1, threshold_label: '∞', threshold_value: null};
        const finiteStops = stops
            .filter(stop => stop.threshold_value !== null)
            .sort((left, right) => left.threshold_value - right.threshold_value);

        if (finiteStops.length === 0) {
            return {
                stops: [{...infinityStop, sliderPosition: 0}],
                maxPosition: 0,
                initialPosition: 0,
                minLabel: '∞',
                maxLabel: '∞',
            };
        }

        const positionedStops = finiteStops.map(stop => ({...stop}));
        positionSliderStops(positionedStops);
        positionedStops.push({...infinityStop, sliderPosition: SLIDER_INFINITY_POSITION});

        return {
            stops: positionedStops,
            maxPosition: SLIDER_INFINITY_POSITION,
            initialPosition: 0,
            minLabel: finiteStops[0].threshold_label,
            maxLabel: '∞',
        };
    }

    // The slider's position is the source of truth, but not a sufficient one, so a stop
    // chosen deliberately is remembered and honored while the slider still sits where
    // that choice put it. Two reasons it has to be:
    //
    //   * positionSliderStops re-lays-out the track whenever the chart's window moves,
    //     and the stop the thumb is on has to survive being given a new position; and
    //   * with more stops than the track has positions -- only reachable by raising
    //     --max_merge_events past ~900 -- stops share positions again, and resolving by
    //     position alone always hands back the first stop at that position. Snapping to
    //     any of the others and asking again would then return a different stop than the
    //     one just chosen, which is what made the step arrows look stuck partway up the
    //     axis before the track was warped.
    //
    // Dragging the slider moves it off that position and the nearest-position search
    // below takes over again, which is right: a drag picks a position, not a stop.
    function currentSliderStop() {
        if (!state.sliderModel || state.sliderModel.stops.length === 0) {
            return null;
        }
        const sliderPosition = Number(document.getElementById('threshold-slider').value);
        if (state.selectedStop && state.selectedStop.sliderPosition === sliderPosition) {
            return state.selectedStop;
        }
        let nearestStop = state.sliderModel.stops[0];
        let nearestDistance = Math.abs(sliderPosition - nearestStop.sliderPosition);
        for (const stop of state.sliderModel.stops) {
            const distance = Math.abs(sliderPosition - stop.sliderPosition);
            if (distance < nearestDistance) {
                nearestStop = stop;
                nearestDistance = distance;
            }
        }
        return nearestStop;
    }

    function snapSliderToStop(stop) {
        if (!stop) {
            return;
        }
        document.getElementById('threshold-slider').value = String(stop.sliderPosition);
        state.selectedStop = stop;
    }

    function nearestStopForThreshold(targetValue) {
        if (!state.sliderModel || state.sliderModel.stops.length === 0 || !Number.isFinite(targetValue)) {
            return null;
        }
        const finiteStops = state.sliderModel.stops.filter(stop => stop.threshold_value !== null);
        if (finiteStops.length === 0) {
            return currentSliderStop();
        }
        let nearestStop = finiteStops[0];
        let nearestDistance = Math.abs(finiteStops[0].threshold_value - targetValue);
        finiteStops.forEach(stop => {
            const distance = Math.abs(stop.threshold_value - targetValue);
            if (distance < nearestDistance) {
                nearestStop = stop;
                nearestDistance = distance;
            }
        });
        return nearestStop;
    }

    // Static menu of the named palettes. The disabled first entry is what a
    // hand-edited column shows; "Reset to defaults" is the way back to
    // DEFAULT_CATEGORICAL_PALETTE, so the menu carries no pseudo-entry for it.
    function populateColorPaletteMenu() {
        const select = document.getElementById('color-palette');
        select.innerHTML = '';
        const customOption = document.createElement('option');
        customOption.value = '';
        customOption.textContent = 'Custom colors';
        customOption.disabled = true;
        select.appendChild(customOption);
        CATEGORICAL_PALETTES.forEach(scheme => {
            const option = document.createElement('option');
            option.value = scheme.name;
            option.textContent = scheme.label;
            option.title = scheme.note;
            select.appendChild(option);
        });
    }

    function populateMetadataControls() {
        const colorBy = document.getElementById('color-by');
        const labelBy = document.getElementById('label-by');
        colorBy.innerHTML = '';
        labelBy.innerHTML = '';

        const emptyOption = document.createElement('option');
        emptyOption.value = '';
        emptyOption.textContent = 'None';
        colorBy.appendChild(emptyOption.cloneNode(true));
        labelBy.appendChild(emptyOption.cloneNode(true));

        // "Label by" can also use the sequence's node_id (not a metadata column).
        const nodeIdOption = document.createElement('option');
        nodeIdOption.value = '__node_id__';
        nodeIdOption.textContent = 'node_id';
        labelBy.appendChild(nodeIdOption);

        state.metadataColumns.forEach(column => {
            const colorOption = document.createElement('option');
            colorOption.value = column.name;
            colorOption.textContent = column.name;
            colorBy.appendChild(colorOption);

            const labelOption = document.createElement('option');
            labelOption.value = column.name;
            labelOption.textContent = column.name;
            labelBy.appendChild(labelOption);
        });

        colorBy.value = state.bundle.defaults.color_by || '';
        labelBy.value = state.bundle.defaults.label_by || '';
    }

    function selectedThresholdValue() {
        const stop = currentSliderStop();
        if (!stop) {
            return Infinity;
        }
        return stop.threshold_value === null ? Infinity : stop.threshold_value;
    }

    function activeClustersAtThreshold(thresholdValue) {
        const hierarchy = state.bundle.graph.hierarchy;
        const active = [];
        const stack = [...hierarchy.roots].reverse();
        while (stack.length > 0) {
            const componentId = stack.pop();
            const component = hierarchy.nodes[componentId];
            if (component.kind === 'leaf') {
                active.push(componentId);
                continue;
            }
            // Strictly-above, matching `build_ssn --threshold thresholdValue`: a component that
            // merged exactly at the threshold is split back apart here.
            if (component.threshold <= thresholdValue) {
                stack.push(component.right);
                stack.push(component.left);
                continue;
            }
            active.push(componentId);
        }
        return active;
    }

    function componentMembers(componentId) {
        const hierarchy = state.bundle.graph.hierarchy;
        const component = hierarchy.nodes[componentId];
        const start = component.leaf_start;
        const count = component.leaf_count;
        return hierarchy.leaf_order.slice(start, start + count);
    }

    function leafPruningOnlyEnabled() {
        return document.getElementById('leaf-pruning-only').checked;
    }

    function collapseLongPathsEnabled() {
        return document.getElementById('collapse-long-paths').checked;
    }

    function updateCollapseLongPathsControl() {
        // Sub-option of leaf pruning, and only meaningful there: with leaf pruning off
        // every below-minimum cluster is dropped outright, so no pass-through chain is
        // left to contract.
        document.getElementById('collapse-long-paths').disabled = !leafPruningOnlyEnabled();
    }

    function currentLayoutAlgorithm() {
        return document.getElementById('layout-algorithm').value || 'grid';
    }

    function compareVisibleClusterIds(leftId, rightId) {
        const leftNode = state.bundle.graph.hierarchy.nodes[leftId];
        const rightNode = state.bundle.graph.hierarchy.nodes[rightId];
        if (sortComponentsBySizeEnabled()) {
            return rightNode.size - leftNode.size || leftNode.leaf_start - rightNode.leaf_start || leftId - rightId;
        }
        return leftNode.leaf_start - rightNode.leaf_start || rightNode.size - leftNode.size || leftId - rightId;
    }

    function sortComponentsBySizeEnabled() {
        return document.getElementById('sort-components-by-size').getAttribute('aria-pressed') === 'true';
    }

    function updateComponentSortButton() {
        const button = document.getElementById('sort-components-by-size');
        const enabled = sortComponentsBySizeEnabled();
        const layoutMode = currentLayoutAlgorithm();
        const gridMode = layoutMode === 'grid' || layoutMode === 'packed';
        if (layoutMode === 'treemap') {
            // Treemap node positions are fixed by the lattice, so sorting has no effect.
            button.disabled = true;
            button.textContent = 'Sort clusters by size: n/a';
            button.title = 'Treemap nodes follow a fixed lattice; size sorting does not apply.';
            return;
        }
        button.disabled = false;
        button.textContent = gridMode
            ? (enabled ? 'Sort clusters by size: On' : 'Sort clusters by size: Off')
            : (enabled ? 'Sort components by size: On' : 'Sort components by size: Off');
        button.title = gridMode
            ? (enabled
                ? 'Visible clusters are ordered by cluster size before grid placement.'
                : 'Visible clusters keep their hierarchy-derived order before grid placement.')
            : (enabled
                ? 'Connected components are ordered by total visible node count.'
                : 'Connected components keep their hierarchy-derived order.');
    }

    function showNodeCountsEnabled() {
        return document.getElementById('show-node-counts').checked;
    }

    function showEdgeScoresEnabled() {
        return document.getElementById('show-edge-scores').checked;
    }

    function renderClusterBoundsEnabled() {
        return document.getElementById('render-cluster-bounds').checked;
    }

    function renderNodesEnabled() {
        return document.getElementById('render-nodes').checked;
    }

    function reduceElongationEnabled() {
        return document.getElementById('reduce-elongation').checked;
    }

    function currentColorField() {
        return document.getElementById('color-by').value;
    }

    function currentLabelField() {
        return document.getElementById('label-by').value;
    }

    function metadataColumn(name) {
        return state.metadataColumnByName.get(name) || null;
    }

    function metadataColumnKeys() {
        return ['node_id', ...state.metadataColumns.map(column => column.name)];
    }

    function defaultMetadataColumnWidth(columnKey) {
        if (columnKey === 'node_id') {
            return 230;
        }
        const column = metadataColumn(columnKey);
        if (!column) {
            return 180;
        }
        if (column.type === 'int' || column.type === 'float') {
            return 140;
        }
        if (column.type === 'bool' || column.type === 'boolean') {
            return 120;
        }
        return 190;
    }

    function metadataColumnWidth(columnKey) {
        return state.metadataColumnWidths.get(columnKey) || defaultMetadataColumnWidth(columnKey);
    }

    function applyMetadataColumnWidths() {
        const colgroup = document.querySelector('#metadata-table colgroup');
        if (!colgroup) {
            return;
        }
        colgroup.innerHTML = '';
        metadataColumnKeys().forEach(columnKey => {
            const col = document.createElement('col');
            col.style.width = metadataColumnWidth(columnKey) + 'px';
            colgroup.appendChild(col);
        });
    }

    function updateMetadataColumnWidth(columnKey, nextWidth) {
        const clampedWidth = Math.max(96, Math.min(640, Math.round(nextWidth)));
        state.metadataColumnWidths.set(columnKey, clampedWidth);
        applyMetadataColumnWidths();
    }

    function startMetadataColumnResize(columnKey, event) {
        event.preventDefault();
        event.stopPropagation();
        state.metadataResize = {
            columnKey,
            startX: event.clientX,
            startWidth: metadataColumnWidth(columnKey),
        };
        document.body.style.cursor = 'col-resize';
        document.body.style.userSelect = 'none';
    }

    function metadataValue(nodeIndex, columnName) {
        if (!columnName) {
            return null;
        }
        const columnIndex = state.metadataColumnIndexByName.get(columnName);
        if (columnIndex === undefined) {
            return null;
        }
        return state.metadataByNodeIndex[nodeIndex]?.[columnIndex] ?? null;
    }

    // Distinct values of a column, ordered by frequency (desc) then key (asc). Keyed by
    // String(value) -- the same key categoricalColor hashes and the TSV/legend use. Null/
    // empty cells are counted separately (nullCount) and never become a key. `cap` bounds
    // the returned list (the picker/legend cap at 200; the TSV save path passes Infinity so
    // the file round-trips); `overflow` flags that more distinct values exist than returned.
    function distinctColumnValues(columnName, cap = 200) {
        const result = {values: [], nullCount: 0, overflow: false, total: 0};
        const columnIndex = state.metadataColumnIndexByName.get(columnName);
        if (columnIndex === undefined) {
            return result;
        }
        const byKey = new Map();
        state.metadataByNodeIndex.forEach(row => {
            const raw = row ? (row[columnIndex] ?? null) : null;
            if (raw === null || raw === undefined || raw === '') {
                result.nullCount += 1;
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
        result.total = byKey.size;
        const sorted = Array.from(byKey.values()).sort(
            (a, b) => (b.count - a.count) || a.key.localeCompare(b.key)
        );
        if (sorted.length > cap) {
            result.overflow = true;
            result.values = sorted.slice(0, cap);
        } else {
            result.values = sorted;
        }
        return result;
    }

    // Normalize a color string to '#RRGGBB' (uppercase), mirroring Python's
    // normalize_color_hex. Returns null for anything that isn't a 6-digit hex.
    // Mirrors color_genbank.normalize_color_hex: an optional '#' and six hex digits,
    // canonicalized to '#RRGGBB'. Loosened by one case for hand-typed input -- CSS-style
    // three-digit shorthand, where each digit doubles. Returns null for anything else.
    function normalizeColorHex(text) {
        const trimmed = String(text).trim().toUpperCase();
        const short = /^#?([0-9A-F]{3})$/.exec(trimmed);
        if (short) {
            const digits = short[1];
            return '#' + digits.charAt(0) + digits.charAt(0)
                + digits.charAt(1) + digits.charAt(1)
                + digits.charAt(2) + digits.charAt(2);
        }
        const match = /^#?([0-9A-F]{6})$/.exec(trimmed);
        return match ? '#' + match[1] : null;
    }

    // Interpolate two '#rrggbb' colors channel-wise in RGB. t is clamped to [0, 1].
    function lerpHexColor(a, b, t) {
        const clampT = Math.max(0, Math.min(1, t));
        const parse = hex => {
            const value = parseInt(hex.slice(1), 16);
            return [(value >> 16) & 255, (value >> 8) & 255, value & 255];
        };
        const left = parse(a);
        const right = parse(b);
        const channel = index => Math.round(left[index] + ((right[index] - left[index]) * clampT));
        const toHex = n => n.toString(16).padStart(2, '0');
        return '#' + toHex(channel(0)) + toHex(channel(1)) + toHex(channel(2));
    }

    // Named qualitative palettes (from domainator.utils.NAMED_CATEGORICAL_PALETTES) and
    // the neutral gray get_palette gives values with no color of their own.
    const CATEGORICAL_PALETTES = DSNV_CONFIG.categoricalPalettes;
    const PALETTE_NO_VALUE_COLOR = DSNV_CONFIG.otherColor;

    function paletteSchemeByName(schemeName) {
        return CATEGORICAL_PALETTES.find(scheme => scheme.name === schemeName) || null;
    }

    // Python's str comparison is by code point; localeCompare is not, and the assignment
    // order has to match domainator.utils.sort_palette_values to reproduce its colors.
    function comparePaletteText(left, right) {
        if (left < right) {
            return -1;
        }
        return left > right ? 1 : 0;
    }

    // Mirror of domainator.utils.sort_palette_values: numeric-looking keys sort
    // numerically (ties broken by text), everything else sorts by text. Nulls never
    // reach here -- distinctColumnValues counts them separately.
    function sortPaletteKeys(keys) {
        const numeric = [];
        const other = [];
        keys.forEach(key => {
            const number = Number(key);
            if (key.trim() !== '' && Number.isFinite(number)) {
                numeric.push({number, key});
            } else {
                other.push(key);
            }
        });
        numeric.sort((a, b) => (a.number - b.number) || comparePaletteText(a.key, b.key));
        other.sort(comparePaletteText);
        return numeric.map(entry => entry.key).concat(other);
    }

    // Fill a column's palette from a named scheme: colors are handed out in
    // sort_palette_values order and cycle once the scheme runs out, exactly like
    // get_palette. Values keep their assigned color afterwards, so individual swatches,
    // the color-table TSV and the legend all stay editable/exportable. `offset` skips
    // that many of the scheme's colors first, like get_palette's offset, so a second
    // column can avoid the colors a first one with that many values used.
    function applyNamedPalette(columnName, schemeName, offset = 0) {
        const scheme = paletteSchemeByName(schemeName);
        if (!scheme) {
            return null;
        }
        const start = paletteOffsetValue(offset);
        const distinct = distinctColumnValues(columnName, Infinity);
        const colors = {};
        sortPaletteKeys(distinct.values.map(entry => entry.key)).forEach((key, index) => {
            colors[key] = scheme.colors[(start + index) % scheme.colors.length];
        });
        const palette = ensureCategoricalPalette(columnName);
        palette.colors = colors;
        palette.nullColor = PALETTE_NO_VALUE_COLOR;
        palette.scheme = scheme.name;
        palette.offset = start;
        return {scheme, assigned: distinct.values.length};
    }

    // A whole number of colors to skip; anything else (blank, negative, fractional,
    // or missing from a palette saved before offsets existed) is clamped into one.
    function paletteOffsetValue(value) {
        const number = Math.floor(Number(value));
        return Number.isFinite(number) && number > 0 ? number : 0;
    }

    // The palette a column has before anyone customizes it. Values used to be
    // colored by a hue hashed from their text -- stable, but arbitrary, and free to
    // put two adjacent categories on nearly the same color. They now get
    // get_palette()'s own 64 distinct colors in sort_palette_values order, so an
    // untouched viewer agrees with what domainator build_ssn would have drawn.
    const DEFAULT_CATEGORICAL_PALETTE = 'domainator';

    // Derived, not user data: never written into a session, and never reported as
    // "custom colors" by the picker. Keyed like state.customPalettes (by
    // paletteKey, so a numeric column flipped to discrete gets its own entry) and
    // dropped by rebuildMetadataCaches(), since the assignment depends on the whole
    // set of values in the column.
    function defaultCategoricalPalette(columnName) {
        const key = paletteKey(columnName);
        const cached = state.defaultPalettes.get(key);
        if (cached) {
            return cached;
        }
        const bundleTable = state.bundle.defaults.color_table;
        if (bundleTable && columnName === state.bundle.defaults.color_by) {
            const palette = bundleColorTablePalette(columnName, bundleTable);
            state.defaultPalettes.set(key, palette);
            return palette;
        }
        const scheme = paletteSchemeByName(DEFAULT_CATEGORICAL_PALETTE);
        const colors = {};
        if (scheme) {
            sortPaletteKeys(distinctColumnValues(columnName, Infinity).values.map(entry => entry.key))
                .forEach((valueKey, index) => {
                    colors[valueKey] = scheme.colors[index % scheme.colors.length];
                });
        }
        // nullColor is deliberately left unset: what color an empty cell gets is a
        // separate choice from the categorical assignment, and categoricalColor's
        // own fallback is already the color the rest of the viewer uses for one.
        const palette = {
            type: 'categorical',
            colors,
            nullColor: null,
            scheme: scheme ? scheme.name : null,
        };
        state.defaultPalettes.set(key, palette);
        return palette;
    }

    // The palette from build_ssn_viewer --color_table. Like build_ssn, values the table
    // leaves out are gray; the em-dash row is the no-value color, as in loadColorTableText.
    function bundleColorTablePalette(columnName, table) {
        const colors = {};
        distinctColumnValues(columnName, Infinity).values.forEach(entry => {
            colors[entry.key] = PALETTE_NO_VALUE_COLOR;
        });
        let nullColor = null;
        Object.keys(table).forEach(value => {
            if (value === '—') {
                nullColor = table[value];
            } else {
                colors[value] = table[value];
            }
        });
        return {type: 'categorical', colors, nullColor, scheme: null};
    }

    // `palette` is a categorical palette ({colors, nullColor}) -- normally the
    // column's own, from customPalette(), which now always supplies one for a
    // discrete column. The hue hashed from the value's text is the last resort for
    // a value the palette says nothing about: one typed in after a color table was
    // loaded, say, or after the palette was hand-edited.
    function categoricalColor(value, palette) {
        if (value === null || value === undefined || value === '') {
            return (palette && palette.nullColor) || '#b3a89d';
        }
        const text = String(value);
        if (palette && palette.colors && palette.colors[text]) {
            return palette.colors[text];
        }
        let hash = 0;
        for (let i = 0; i < text.length; i++) {
            hash = ((hash << 5) - hash) + text.charCodeAt(i);
            hash |= 0;
        }
        const hue = Math.abs(hash) % 360;
        return 'hsl(' + hue + ' 58% 54%)';
    }

    // `palette` is an optional numeric custom palette ({stops, nullColor}).
    // With low+high set, interpolate in RGB (low->mid->high split at 0.5 when mid is
    // present); otherwise fall back to the default cyan->orange hsl gradient.
    function numericColor(value, minValue, maxValue, palette) {
        if (value === null || value === undefined || Number.isNaN(value)) {
            return (palette && palette.nullColor) || '#b3a89d';
        }
        // A custom palette is a list of {value, color} stops, kept sorted by whoever
        // stored it. Interpolate between the pair that brackets the value and hold the
        // end colors flat outside the ends -- so two stops are a plain ramp, three
        // reproduce a low/mid/high midpoint, and more shape the ramp arbitrarily.
        const stops = palette && palette.stops;
        if (stops && stops.length >= 2) {
            if (value <= stops[0].value) {
                return stops[0].color;
            }
            const last = stops[stops.length - 1];
            if (value >= last.value) {
                return last.color;
            }
            for (let index = 1; index < stops.length; index++) {
                const upper = stops[index];
                if (value <= upper.value) {
                    const lower = stops[index - 1];
                    const span = upper.value - lower.value;
                    return span > 0
                        ? lerpHexColor(lower.color, upper.color, (value - lower.value) / span)
                        : upper.color;
                }
            }
            return last.color;
        }
        // No custom palette: the built-in hue ramp across the column's own range.
        const fraction = maxValue <= minValue ? 0.5 : (value - minValue) / (maxValue - minValue);
        const clamped = Math.max(0, Math.min(1, fraction));
        const hue = 200 - (160 * clamped);
        const light = 72 - (24 * clamped);
        return 'hsl(' + hue + ' 72% ' + light + '%)';
    }

    function colorInfo(columnName) {
        return state.metadataColorInfoByName.get(columnName) || null;
    }

    // True when `columnName` holds numbers, and so can be colored either as a gradient
    // or as discrete categories.
    function columnIsNumericType(columnName) {
        const info = colorInfo(columnName);
        return !!info && info.baseType === 'numeric';
    }

    function columnIsCategoricalNumeric(columnName) {
        return columnIsNumericType(columnName) && state.categoricalColumns.has(columnName);
    }

    // Switch a numeric column between gradient and discrete coloring. No-op for columns
    // that are not numeric (those are always categorical).
    function setColumnCategorical(columnName, categorical) {
        const info = colorInfo(columnName);
        if (!info || info.baseType !== 'numeric') {
            return;
        }
        if (categorical) {
            state.categoricalColumns.add(columnName);
        } else {
            state.categoricalColumns.delete(columnName);
        }
        info.type = categorical ? 'categorical' : 'numeric';
    }

    // Custom palettes are stored per column AND per coloring mode, so flipping a numeric
    // column to categorical and back does not discard the gradient the user set up.
    function paletteKey(columnName) {
        return columnIsCategoricalNumeric(columnName) ? columnName + '\u0000categorical' : columnName;
    }

    function customPalette(columnName) {
        const stored = state.customPalettes[paletteKey(columnName)];
        if (stored) {
            return stored;
        }
        // A discrete column with nothing stored is not uncolored -- it is on the
        // default named palette. Only gradient columns fall through to null, where
        // numericColor's own built-in ramp takes over.
        const info = colorInfo(columnName);
        return info && info.type === 'categorical' ? defaultCategoricalPalette(columnName) : null;
    }

    function nodeColor(nodeIndex) {
        const columnName = currentColorField();
        const info = colorInfo(columnName);
        const value = metadataValue(nodeIndex, columnName);
        if (!info) {
            return '#d88f3d';
        }
        const palette = customPalette(columnName);
        if (info.type === 'numeric') {
            return numericColor(value, info.min, info.max, palette);
        }
        return categoricalColor(value, palette);
    }

    function activeClusterAssignments(activeClusterIds) {
        const assignments = new Array(state.bundle.graph.nodes.length).fill(-1);
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
        const hierarchyNodes = state.bundle.graph.hierarchy.nodes;
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

        state.bundle.graph.mst_edges.forEach(edge => {
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
            const visibleIds = activeClusterIds.filter(componentId => state.bundle.graph.hierarchy.nodes[componentId].size >= minClusterSize);
            const visibleSet = new Set(visibleIds);
            const hiddenNodes = activeClusterIds.reduce((sum, componentId) => {
                if (visibleSet.has(componentId)) {
                    return sum;
                }
                return sum + state.bundle.graph.hierarchy.nodes[componentId].size;
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
                const size = state.bundle.graph.hierarchy.nodes[componentId].size;
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
            return sum + state.bundle.graph.hierarchy.nodes[componentId].size;
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
        const hierarchy = state.bundle.graph.hierarchy;
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
        const nodes = state.bundle.graph.hierarchy.nodes;
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
            const node = state.bundle.graph.hierarchy.nodes[componentId];
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
            radius: componentRadiusForSize(state.bundle.graph.hierarchy.nodes[componentId].size),
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
        const leafOrder = state.bundle.graph.hierarchy.leaf_order;
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
        const nodes = state.bundle.graph.hierarchy.nodes;
        const leafOrder = state.bundle.graph.hierarchy.leaf_order;
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

    // Axis titles for the split chart. A `product` merge impact is a product of two
    // component sizes, not a count of nodes, so it must not be labeled as one.
    function splitAxisLabels() {
        const labels = DSNV_CONFIG.splitAxisLabels;
        const metric = state.bundle && state.bundle.graph
            ? state.bundle.graph.merge_impact_metric
            : null;
        return labels[metric] || labels[DSNV_CONFIG.defaultMergeImpactMetric];
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

    function worldToScreenPoint(x, y) {
        return {
            x: (x * state.viewTransform.scale) + state.viewTransform.offsetX,
            y: (y * state.viewTransform.scale) + state.viewTransform.offsetY,
        };
    }

    function screenToWorldPoint(x, y) {
        return {
            x: (x - state.viewTransform.offsetX) / Math.max(state.viewTransform.scale, 1e-9),
            y: (y - state.viewTransform.offsetY) / Math.max(state.viewTransform.scale, 1e-9),
        };
    }

    function canvasCoordinatesFromEvent(event) {
        const rect = clusterCanvas.getBoundingClientRect();
        const scaleX = clusterCanvas.width / rect.width;
        const scaleY = clusterCanvas.height / rect.height;
        return {
            x: (event.clientX - rect.left) * scaleX,
            y: (event.clientY - rect.top) * scaleY,
        };
    }

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
            const component = state.bundle.graph.hierarchy.nodes[item.componentId];
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

    // The zoom box shows the scale (canvas pixels per layout unit) and sets it, so a
    // zoom can be copied from one viewer into another. A typed zoom keeps the center of
    // the canvas in place, as the wheel keeps the point under the pointer.
    function formatZoomLevel(scale) {
        return String(Number(scale.toPrecision(4)));
    }

    function updateZoomLevelControl() {
        const input = document.getElementById('zoom-level');
        // Not while the user is typing in it: the force layout repaints on its own.
        if (document.activeElement === input) {
            return;
        }
        input.value = formatZoomLevel(state.viewTransform.scale);
    }

    function changeZoomLevel() {
        const input = document.getElementById('zoom-level');
        const scale = Number(input.value);
        // Number('') is 0, so a cleared box is refused here too.
        if (state.bundle && Number.isFinite(scale) && scale > 0) {
            const transform = state.viewTransform;
            const nextScale = Math.max(transform.minScale, Math.min(transform.maxScale, scale));
            const centerX = clusterCanvas.width / 2;
            const centerY = clusterCanvas.height / 2;
            const anchor = screenToWorldPoint(centerX, centerY);
            transform.scale = nextScale;
            transform.offsetX = centerX - (anchor.x * nextScale);
            transform.offsetY = centerY - (anchor.y * nextScale);
            renderClusterView();
        }
        // Shows the clamped scale, or puts back the current one after bad input.
        input.value = formatZoomLevel(state.viewTransform.scale);
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
    function itemScreenExtent(item) {
        const scale = state.viewTransform.scale;
        if (item.shape === 'lattice' || item.shape === 'rect') {
            return {width: (item.x1 - item.x0) * scale, height: (item.y1 - item.y0) * scale};
        }
        const diameter = item.radius * 2 * scale;
        return {width: diameter, height: diameter};
    }

    function clusterCountLabelFits(text, item) {
        const extent = itemScreenExtent(item);
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
        // Every pan and zoom ends in a repaint. An export repaints the same transform,
        // so refreshing the box then too is harmless.
        updateZoomLevelControl();
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
        const drawScale = Math.max(state.viewTransform.scale, 1e-9);

        // World-space viewport bounds for culling (#4)
        const worldMinX = (0 - state.viewTransform.offsetX) / drawScale;
        const worldMaxX = (clusterCanvas.width - state.viewTransform.offsetX) / drawScale;
        const worldMinY = (0 - state.viewTransform.offsetY) / drawScale;
        const worldMaxY = (clusterCanvas.height - state.viewTransform.offsetY) / drawScale;

        clusterContext.save();
        clusterContext.translate(state.viewTransform.offsetX, state.viewTransform.offsetY);
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

            const component = state.bundle.graph.hierarchy.nodes[item.componentId];
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
                    const buckets = dot.square ? squareColorBuckets : dotColorBuckets;
                    let bucket = buckets.get(color);
                    if (bucket === undefined) {
                        bucket = [];
                        buckets.set(color, bucket);
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
                    if (dot.radius * state.viewTransform.scale < MIN_LABELED_DOT_SCREEN_RADIUS) {
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
        dotColorBuckets.forEach((dots, color) => {
            clusterContext.fillStyle = color;
            clusterContext.beginPath();
            for (const dot of dots) {
                clusterContext.moveTo(dot.x + dot.radius, dot.y);
                clusterContext.arc(dot.x, dot.y, dot.radius, 0, TAU);
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
                    const topLeft = worldToScreenPoint(item.x0, item.y0);
                    const bottomRight = worldToScreenPoint(item.x1, item.y1);
                    clusterContext.rect(topLeft.x, topLeft.y, bottomRight.x - topLeft.x, bottomRight.y - topLeft.y);
                } else {
                    const center = worldToScreenPoint(item.x, item.y);
                    clusterContext.arc(center.x, center.y, (item.radius * state.viewTransform.scale) + 3, 0, TAU);
                }
                clusterContext.stroke();
            });
            previewNodeOutlines.forEach(outline => {
                const screenPoint = worldToScreenPoint(outline.x, outline.y);
                const screenRadius = Math.max(4, (outline.radius * state.viewTransform.scale) + 3.4);
                clusterContext.beginPath();
                clusterContext.arc(screenPoint.x, screenPoint.y, screenRadius, 0, TAU);
                clusterContext.stroke();
            });
            clusterContext.restore();
        }

        if (selectedNodeOutlines.length > 0) {
            clusterContext.save();
            selectedNodeOutlines.forEach(outline => {
                const screenPoint = worldToScreenPoint(outline.x, outline.y);
                const screenRadius = Math.max(3, (outline.radius * state.viewTransform.scale) + 1.6);
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
                if (!edgeScoreLabelFits(text, anchor.length * state.viewTransform.scale)) {
                    return;
                }
                const point = worldToScreenPoint(anchor.x, anchor.y);
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
                const screenPoint = worldToScreenPoint(item.x, item.y);
                const screenRadius = item.radius * state.viewTransform.scale;
                if (screenPoint.x + screenRadius < 0 || screenPoint.x - screenRadius > clusterCanvas.width ||
                    screenPoint.y + screenRadius < 0 || screenPoint.y - screenRadius > clusterCanvas.height) {
                    return;
                }
                const component = state.bundle.graph.hierarchy.nodes[item.componentId];
                const text = component.size.toLocaleString();
                if (!clusterCountLabelFits(text, item)) {
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
                const screenPoint = worldToScreenPoint(candidate.x, candidate.y);
                clusterContext.fillText(String(text).slice(0, 24), screenPoint.x, screenPoint.y);
            });
        }

        if (state.selectionBox) {
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

        if (state.layoutComputing) {
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
            const leafOrder = state.bundle.graph.hierarchy.leaf_order;
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
        const hierarchyNodes = state.bundle.graph.hierarchy.nodes;
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
        const shown = state.bundle.graph.nodes.length - hidden;
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

// @include layout_core.js
// @include session_state.js
// @include selection_presets.js
// @include table_editing.js
// @include extraction.js
// @include gradient_stops.js
// @include column_charts.js
// @include split_chart_hover.js

    function htmlEscape(value) {
        return String(value).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
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

    function setupMetadataTableDelegation() {
        const thead = document.querySelector('#metadata-table thead');
        const tbody = document.querySelector('#metadata-table tbody');

        thead.addEventListener('click', event => {
            // Checked before the sort button: the chart glyph must never be
            // given a data-column-key, or clicking it would also sort.
            const chartBtn = event.target.closest('[data-chart-column]');
            if (chartBtn) {
                event.preventDefault();
                event.stopPropagation();
                openColumnChartMenu(chartBtn.dataset.chartColumn, chartBtn);
                return;
            }
            const sortBtn = event.target.closest('[data-column-key]');
            if (sortBtn) { toggleMetadataSort(sortBtn.dataset.columnKey); return; }
            const copyBtn = event.target.closest('[data-copy-column]');
            if (copyBtn) { event.preventDefault(); event.stopPropagation(); copyMetadataColumn(copyBtn.dataset.copyColumn); }
        });

        thead.addEventListener('pointerdown', event => {
            const resizeHandle = event.target.closest('[data-resize-column]');
            if (resizeHandle) { startMetadataColumnResize(resizeHandle.dataset.resizeColumn, event); }
        });

        tbody.addEventListener('mousedown', event => {
            if (event.shiftKey && !event.target.closest('button, a, input, select, textarea')) {
                event.preventDefault();
            }
        });

        tbody.addEventListener('click', event => {
            if (event.target.closest('button, a, input, select, textarea')) { return; }
            // The second click of a double-click opens the cell editor; letting it
            // through would toggle the row selection on the way in.
            if (event.detail > 1) { return; }
            if (!event.shiftKey && metadataTableHasActiveTextSelection()) { return; }
            const row = event.target.closest('tr[data-node-index]');
            if (!row) { return; }
            toggleMetadataRowSelection(Number(row.dataset.nodeIndex), state.renderedNodeIndices, {range: event.shiftKey});
        });

        tbody.addEventListener('keydown', event => {
            if (event.key !== 'Enter' && event.key !== ' ') { return; }
            if (event.target.closest('button, a, input, select, textarea')) { return; }
            const row = event.target.closest('tr[data-node-index]');
            if (!row) { return; }
            event.preventDefault();
            toggleMetadataRowSelection(Number(row.dataset.nodeIndex), state.renderedNodeIndices, {range: event.shiftKey});
        });
    }

    function updateMetadataTable() {
        // The table body is re-rendered wholesale, so any in-progress cell edit
        // is destroyed along with its <input>. The header goes the same way, so
        // an open chart menu would be left anchored to a detached button.
        state.metadataEditCell = null;
        closeColumnChartMenu();
        const selected = Array.from(state.selectedNodeIndices).sort((left, right) => left - right);
        const baseNodeIndices = metadataBaseNodeIndices();
        const filteredNodeIndices = filteredMetadataNodeIndices(baseNodeIndices);
        const sortedNodeIndices = sortedMetadataNodeIndices(filteredNodeIndices);
        const pagination = metadataPagination(sortedNodeIndices.length);
        state.metadataPage = pagination.pageIndex;
        const renderedNodeIndices = sortedNodeIndices.slice(pagination.start, pagination.end);
        state.renderedNodeIndices = renderedNodeIndices;
        pruneMetadataRowSelection(renderedNodeIndices);
        applyMetadataColumnWidths();
        const thead = document.querySelector('#metadata-table thead');
        const tbody = document.querySelector('#metadata-table tbody');

        const columnKeys = metadataColumnKeys();
        const headerCells = columnKeys.map(label => {
            const escapedLabel = htmlEscape(label);
            const ariasort = state.metadataSort.columnKey === label
                ? (state.metadataSort.direction === 'asc' ? 'ascending' : 'descending')
                : 'none';
            const indicator = htmlEscape(metadataSortIndicator(label));
            return '<th class="metadata-header" aria-sort="' + ariasort + '">'
                + '<div class="metadata-header-cell">'
                + '<button type="button" class="metadata-sort-button" data-column-key="' + escapedLabel + '" title="Sort by ' + escapedLabel + '">'
                + '<span>' + escapedLabel + '</span><span class="metadata-sort-indicator">' + indicator + '</span></button>'
                // No chart glyph on node_id: it is the graph key, unique by
                // construction, so every chart of it is degenerate (N rows of
                // count 1, N equal slices). A disabled button would only invite
                // a "why?" with no answer.
                + (label === 'node_id' ? ''
                    : '<button type="button" class="metadata-chart-button" data-chart-column="' + escapedLabel
                        + '" aria-haspopup="menu" aria-expanded="false" title="Charts and tables for ' + escapedLabel + '">'
                        + '<svg class="metadata-chart-glyph" viewBox="0 0 12 12" aria-hidden="true" focusable="false">'
                        + '<rect x="1" y="6" width="2.6" height="5"></rect>'
                        + '<rect x="4.7" y="3" width="2.6" height="8"></rect>'
                        + '<rect x="8.4" y="1" width="2.6" height="10"></rect></svg></button>')
                + '<button type="button" class="metadata-copy-button" data-copy-column="' + escapedLabel + '" title="Copy the currently displayed values from ' + escapedLabel + '">Copy</button>'
                + '</div>'
                + '<div class="metadata-resize-handle" role="separator" aria-orientation="vertical" title="Drag to resize column" data-resize-column="' + escapedLabel + '"></div>'
                + '</th>';
        });
        thead.innerHTML = '<tr>' + headerCells.join('') + '</tr>';

        let bodyHtml;
        if (sortedNodeIndices.length === 0) {
            const emptyMsg = htmlEscape(metadataFilterText()
                ? 'No metadata rows match the current filter.'
                : 'No metadata rows to display.');
            bodyHtml = '<tr class="metadata-empty-row"><td colspan="' + (state.metadataColumns.length + 1) + '">' + emptyMsg + '</td></tr>';
        } else {
            bodyHtml = renderedNodeIndices.map(nodeIndex => {
                const isSelected = state.selectedMetadataNodeIndices.has(nodeIndex);
                const rowClass = isSelected ? ' class="metadata-row-selected"' : '';
                const ariaSelected = isSelected ? 'true' : 'false';
                const idText = htmlEscape(nodeId(nodeIndex));
                let cellsHtml = '<td title="' + idText + '">' + idText + '</td>';
                for (const column of state.metadataColumns) {
                    const value = metadataValue(nodeIndex, column.name);
                    const cellClass = metadataCellClass(column.name, value);
                    const escapedName = htmlEscape(column.name);
                    const displayText = htmlEscape(formatMetadataDisplayValue(column.name, value));
                    cellsHtml += '<td data-column-key="' + escapedName + '"'
                        + ' class="metadata-cell-editable' + (cellClass ? ' ' + cellClass : '') + '"'
                        + ' title="' + displayText + ' (double-click to edit)">' + displayText + '</td>';
                }
                return '<tr' + rowClass + ' tabindex="0" aria-selected="' + ariaSelected + '" data-node-index="' + nodeIndex + '">' + cellsHtml + '</tr>';
            }).join('');
        }
        tbody.innerHTML = bodyHtml;

        document.getElementById('selection-summary').textContent = selected.length.toLocaleString() + ' nodes selected';
        const pageStatus = document.getElementById('metadata-page-status');
        if (sortedNodeIndices.length === 0) {
            pageStatus.textContent = metadataFilterText()
                ? 'No rows match the current filter.'
                : 'No metadata rows available.';
        } else if (pagination.showAll) {
            pageStatus.textContent = 'Showing all ' + sortedNodeIndices.length.toLocaleString() + ' rows on one page.';
        } else {
            pageStatus.textContent = 'Showing rows ' + (pagination.start + 1).toLocaleString() + ' to ' + pagination.end.toLocaleString() +
                ' of ' + sortedNodeIndices.length.toLocaleString() + ' (page ' + (pagination.pageIndex + 1).toLocaleString() +
                ' of ' + pagination.pageCount.toLocaleString() + ').';
        }
        document.getElementById('metadata-prev-page').disabled = !state.bundle || sortedNodeIndices.length === 0 || pagination.pageIndex === 0;
        document.getElementById('metadata-next-page').disabled = !state.bundle || sortedNodeIndices.length === 0 || pagination.pageIndex >= pagination.pageCount - 1;
        const sortDescription = metadataSortDescription();
        const filterText = metadataFilterText();
        const filterDescription = filterText ? ' matching filter "' + filterText + '"' : '';
        const pagerActive = pagination.pageCount > 1;
        state.metadataBaseNoteText = (selected.length === 0
            ? (
                pagerActive
                    ? 'No clusters selected. Browse the full network table with the pager' + filterDescription + (sortDescription ? ', sorted by ' + sortDescription : '') + '.'
                    : 'No clusters selected. Showing the full network table' + filterDescription + (sortDescription ? ', sorted by ' + sortDescription : '') + '.'
            )
            : (
                pagerActive
                    ? 'Showing one page of the selected rows in the table' + filterDescription + (sortDescription ? ', sorted by ' + sortDescription : '') + '. Use the pager to browse more rows. Export includes the full filtered selection.'
                    : 'Ctrl-click a node to toggle it individually, click a cluster to toggle it, Shift-drag a box to add multiple clusters, click table rows to stage them, shift-click to select or deselect row ranges, use the search box to filter rows, and click a column header to sort' + (sortDescription ? ' by ' + sortDescription : '') + '.'
            ));
        document.getElementById('export-selected').disabled = !state.bundle;
        document.getElementById('metadata-reset-sort').disabled = !state.bundle || !state.metadataSort.columnKey;
        document.getElementById('clear-selection').disabled = selected.length === 0;
        document.getElementById('focus-selection').disabled = selected.length === 0;
        document.getElementById('save-extraction').disabled = !state.bundle || selected.length === 0;
        applyMetadataTableRowHighlights();
        // Charts summarize the rows this function just resolved, so a changed
        // selection, filter or edit has to reach an open chart.
        refreshColumnChartIfOpen();
    }

    function applyMetadataTableRowHighlights() {
        const tbody = document.querySelector('#metadata-table tbody');
        if (tbody) {
            tbody.querySelectorAll('tr[data-node-index]').forEach(row => {
                const idx = Number(row.dataset.nodeIndex);
                const sel = state.selectedMetadataNodeIndices.has(idx);
                row.className = sel ? 'metadata-row-selected' : '';
                row.setAttribute('aria-selected', sel ? 'true' : 'false');
            });
        }
        const metadataSelectionCount = state.selectedMetadataNodeIndices.size;
        const metadataSelectionDescription = metadataSelectionCount > 0
            ? metadataSelectionCount.toLocaleString() + ' table rows selected. Click Select nodes to promote them into the graph selection. '
            : '';
        document.getElementById('selection-note').textContent = metadataSelectionDescription + state.metadataBaseNoteText;
        document.getElementById('metadata-select-nodes').disabled = !state.bundle || metadataSelectionCount === 0;
        document.getElementById('metadata-deselect-rows').disabled = !state.bundle || metadataSelectionCount === 0;
        // "Add to selection" searches the whole network, so it works from a standing
        // start; the two narrowing actions have nothing to narrow without a selection.
        document.getElementById('metadata-select-add').disabled = !state.bundle;
        document.getElementById('metadata-select-remove').disabled = !state.bundle || state.selectedNodeIndices.size === 0;
        document.getElementById('metadata-select-subset').disabled = !state.bundle || state.selectedNodeIndices.size === 0;
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
            const leftNode = state.bundle.graph.hierarchy.nodes[leftId];
            const rightNode = state.bundle.graph.hierarchy.nodes[rightId];
            return rightNode.size - leftNode.size || leftNode.leaf_start - rightNode.leaf_start || leftId - rightId;
        });
        const clusterByComponentId = new Map(rankedClusterIds.map((componentId, clusterIndex) => [componentId, clusterIndex + 1]));
        return {
            clusterCount: activeClusterIds.length,
            numberFor: nodeIndex => clusterByComponentId.get(componentAssignments[nodeIndex]) ?? null,
        };
    }

    function exportSelection() {
        const selected = Array.from(state.selectedNodeIndices).sort((left, right) => left - right);
        if (!state.bundle) {
            return;
        }
        const clusterNumbers = clusterNumbersAtCurrentThreshold();
        const metadataColumns = state.metadataColumns.filter(column => column.name !== 'SSN_cluster');
        const exportedNodeIndices = metadataDisplayNodeIndices(selected.length > 0
            ? selected
            : state.bundle.graph.nodes.map((_, nodeIndex) => nodeIndex));
        const header = ['node_id', 'SSN_cluster', ...metadataColumns.map(column => column.name)];
        const lines = [header.join('\t')];
        exportedNodeIndices.forEach(nodeIndex => {
            const row = [
                nodeId(nodeIndex),
                clusterNumbers.numberFor(nodeIndex) ?? '',
                ...metadataColumns.map(column => {
                    const value = metadataValue(nodeIndex, column.name);
                    return value === null || value === undefined ? '' : String(value);
                }),
            ];
            lines.push(row.join('\t'));
        });
        const blob = new Blob([lines.join('\n')], {type: 'text/tab-separated-values'});
        const link = document.createElement('a');
        link.href = URL.createObjectURL(blob);
        link.download = (
            selected.length > 0
                ? ((state.bundle?.name || 'selection') + '_selected.tsv')
                : ((state.bundle?.name || 'network') + '_table.tsv')
        ).replace(/\s+/g, '_');
        link.click();
        URL.revokeObjectURL(link.href);
    }

    function exportBaseName() {
        return ((state.bundle && state.bundle.name) || 'network').replace(/\s+/g, '_');
    }

    function triggerDownload(href, filename, revoke) {
        const link = document.createElement('a');
        link.href = href;
        link.download = filename;
        link.click();
        if (revoke) {
            URL.revokeObjectURL(href);
        }
    }

    function selectedPngScale() {
        const select = document.getElementById('export-png-scale');
        const value = select ? Number(select.value) : 1;
        return Number.isFinite(value) && value >= 1 ? value : 1;
    }

    function exportClusterPNG() {
        if (!state.bundle) {
            return;
        }
        const scaleFactor = selectedPngScale();
        // Render into an offscreen canvas scaled by the chosen factor so the PNG
        // is a true high-resolution re-rasterization (not an upscaled snapshot).
        // The view transform and logical dimensions are unchanged, so the image
        // shows exactly the current view at higher pixel density. Rendering
        // offscreen also keeps any transient selection box out of the export.
        const target = document.createElement('canvas');
        target.width = Math.round(clusterCanvas.width * scaleFactor);
        target.height = Math.round(clusterCanvas.height * scaleFactor);
        const targetContext = target.getContext('2d');
        if (!targetContext) {
            window.alert('Could not export PNG at ' + scaleFactor + '×; try a lower resolution.');
            return;
        }
        targetContext.scale(scaleFactor, scaleFactor);
        renderClusterView(targetContext, clusterCanvas.width, clusterCanvas.height, {preview: false});
        target.toBlob(blob => {
            if (!blob) {
                // Browsers (notably Safari, ~16.7M px) refuse to encode an
                // oversized canvas and hand back null. Tell the user instead of
                // failing silently.
                window.alert('The view is too large to export as a ' + scaleFactor +
                    '× PNG (' + target.width + '×' + target.height + ' px). Try a lower resolution or zoom in.');
                return;
            }
            const href = URL.createObjectURL(blob);
            const suffix = scaleFactor > 1 ? '_view@' + scaleFactor + 'x.png' : '_view.png';
            triggerDownload(href, exportBaseName() + suffix, true);
        }, 'image/png');
    }

    function escapeXml(text) {
        return String(text)
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;');
    }

    // Normalize any CSS color (hex, rgb, rgba, hsl) into an SVG-friendly
    // {color, opacity} pair using the canvas's own color parser. SVG editors
    // such as Illustrator/Inkscape choke on rgba()/space-separated hsl(), so we
    // split out the alpha into a separate opacity attribute.
    function svgColorParts(css) {
        clusterContext.fillStyle = '#000000';
        clusterContext.fillStyle = css;
        const normalized = clusterContext.fillStyle;
        const match = /^rgba\(\s*([\d.]+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)\s*\)$/.exec(normalized);
        if (match) {
            return {color: 'rgb(' + match[1] + ',' + match[2] + ',' + match[3] + ')', opacity: Number(match[4])};
        }
        return {color: normalized, opacity: 1};
    }

    // Normalize any CSS color (hex, rgb, hsl) to '#rrggbb' using the canvas parser.
    // <input type="color"> only accepts '#rrggbb', but the default schemes emit
    // hsl(...), so color inputs must be seeded through this. Alpha is dropped.
    function cssToHex(css) {
        clusterContext.fillStyle = '#000000';
        clusterContext.fillStyle = css;
        const normalized = clusterContext.fillStyle;
        if (normalized.charAt(0) === '#') {
            return normalized.length === 7 ? normalized : normalized.slice(0, 7);
        }
        const match = /rgba?\(\s*([\d.]+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)/.exec(normalized);
        if (!match) {
            return '#000000';
        }
        const toHex = n => Math.round(Number(n)).toString(16).padStart(2, '0');
        return '#' + toHex(match[1]) + toHex(match[2]) + toHex(match[3]);
    }

    // Render the current cluster view to a standalone SVG document. Geometry is
    // baked into screen space (every point passed through worldToScreenPoint)
    // so the markup needs no transforms and matches the canvas pixel-for-pixel.
    // This mirrors renderClusterView(); keep the two in sync.
    function buildClusterViewSVG() {
        const width = clusterCanvas.width;
        const height = clusterCanvas.height;
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
        const drawScale = Math.max(state.viewTransform.scale, 1e-9);
        const worldMinX = (0 - state.viewTransform.offsetX) / drawScale;
        const worldMaxX = (width - state.viewTransform.offsetX) / drawScale;
        const worldMinY = (0 - state.viewTransform.offsetY) / drawScale;
        const worldMaxY = (height - state.viewTransform.offsetY) / drawScale;

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
                const start = worldToScreenPoint(segment.startX, segment.startY);
                const end = worldToScreenPoint(segment.endX, segment.endY);
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

            const component = state.bundle.graph.hierarchy.nodes[item.componentId];
            const members = componentMembers(item.componentId);
            const selectionState = componentSelectionState(members);
            const center = worldToScreenPoint(item.x, item.y);
            const screenRadius = item.radius * drawScale;
            const isLattice = item.shape === 'lattice';
            const isRect = item.shape === 'rect';
            const rectTopLeft = isRect ? worldToScreenPoint(item.x0, item.y0) : null;
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
                        const a = worldToScreenPoint(wx0, wy0);
                        const b = worldToScreenPoint(wx1, wy1);
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
                    const buckets = dot.square ? squareColorBuckets : dotColorBuckets;
                    let bucket = buckets.get(color);
                    if (bucket === undefined) {
                        bucket = [];
                        buckets.set(color, bucket);
                    }
                    const dotCenter = worldToScreenPoint(dot.x, dot.y);
                    bucket.push({x: dotCenter.x, y: dotCenter.y, r: dot.radius * drawScale});
                }
                if (!selectionState.allSelected && state.selectedNodeIndices.has(dot.memberIndex)) {
                    selectedNodeOutlines.push({x: dot.x, y: dot.y, radius: dot.radius, faint: !showNodes});
                }
                if (collectDotLabels && !dotLabelsOverflow &&
                    dot.x >= worldMinX && dot.x <= worldMaxX && dot.y >= worldMinY && dot.y <= worldMaxY) {
                    if (dot.radius * state.viewTransform.scale < MIN_LABELED_DOT_SCREEN_RADIUS) {
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

        dotColorBuckets.forEach((dots, color) => {
            const fill = svgColorParts(color);
            const circles = dots.map(dot => '<circle cx="' + fmt(dot.x) + '" cy="' + fmt(dot.y) + '" r="' + fmt(dot.r) + '"/>').join('');
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
                const screenPoint = worldToScreenPoint(outline.x, outline.y);
                const screenRadius = Math.max(3, (outline.radius * state.viewTransform.scale) + 1.6);
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
                if (!edgeScoreLabelFits(text, anchor.length * state.viewTransform.scale)) {
                    return;
                }
                const point = worldToScreenPoint(anchor.x, anchor.y);
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
                const screenPoint = worldToScreenPoint(item.x, item.y);
                const screenRadius = item.radius * state.viewTransform.scale;
                if (screenPoint.x + screenRadius < 0 || screenPoint.x - screenRadius > width ||
                    screenPoint.y + screenRadius < 0 || screenPoint.y - screenRadius > height) {
                    return;
                }
                const component = state.bundle.graph.hierarchy.nodes[item.componentId];
                const text = component.size.toLocaleString();
                if (!clusterCountLabelFits(text, item)) {
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
                const screenPoint = worldToScreenPoint(candidate.x, candidate.y);
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
        triggerDownload(href, exportBaseName() + '_view.svg', true);
    }

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
        renderClusterView();
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
                scheduleClusterRender();
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
            scheduleClusterRender();
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
        for (let nodeIndex = 0; nodeIndex < state.bundle.graph.nodes.length; nodeIndex++) {
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
            scheduleClusterRender();
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
        renderClusterView();
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
        renderClusterView();
        openColorPicker();
    }

    function resetColorPicker() {
        const columnName = currentColorField();
        if (columnName) {
            delete state.customPalettes[paletteKey(columnName)];
            rebuildNodeColorCache();
            scheduleClusterRender();
            openColorPicker();
        }
    }

    function updateThresholdUI(resetView = true) {
        if (!state.bundle) {
            return;
        }
        const stop = currentSliderStop();
        if (!stop) {
            return;
        }
        document.getElementById('threshold-label').textContent = stop.threshold_label;
        document.getElementById('threshold-input').value = stop.threshold_value === null ? '' : String(stop.threshold_value);
        updateThresholdStepButtons();
        drawSplitChart();
        drawClusterView(resetView);
    }

    // The bundle carries the merge order; everything the threshold slider and the split
    // chart show is derived from it, here, on load. See deriveMergeSeries.
    function normalizeMaxMergeEvents(value) {
        const cap = Number(value);
        if (!Number.isFinite(cap) || cap < 0) { return DEFAULT_WINDOWED_MAX_MERGE_EVENTS; }
        return Math.floor(cap);
    }

    // The part of the axis the split chart is currently showing, in threshold units, or
    // null when the whole series is on screen. Clamped by splitChartVisibleWindow against
    // the series' real range, so a stale zoom left over from another bundle cannot narrow
    // the selection to a stretch that no longer exists.
    function currentSplitWindow() {
        if (!state.series || !state.splitChartZoom) { return null; }
        const range = splitSeriesDataRange(state.series);
        const view = splitChartVisibleWindow(range.min, range.max);
        return view.zoomed ? {min: view.min, max: view.max} : null;
    }

    // Re-select which events the chart plots and the slider stops at, and re-lay-out the
    // slider, without redrawing. Every caller that changes the cap or the window goes
    // through here.
    //
    // `pinnedThreshold` is the cut currently in effect. It has to survive the
    // re-selection: if zooming could filter out the stop the thumb is on, the slider
    // would land on a neighboring stop and the network would silently re-cluster.
    function applyMergeSelection(pinnedThreshold) {
        selectMergeEvents(state.series, state.maxMergeEvents, currentSplitWindow(),
                          Number.isFinite(pinnedThreshold) ? pinnedThreshold : null);
        state.sliderModel = buildSliderModel(state.series.sliderStops);
        // Point into the model just replaced, so they cannot outlive it.
        state.selectedStop = null;
        state.stopBeforeSplitChartClick = null;

        const slider = document.getElementById('threshold-slider');
        slider.max = String(state.sliderModel.maxPosition);
        slider.disabled = state.sliderModel.stops.length === 0;
        document.getElementById('threshold-min-label').textContent = state.sliderModel.minLabel;
        document.getElementById('threshold-max-label').textContent = state.sliderModel.maxLabel;
        document.getElementById('threshold-input').disabled = state.sliderModel.stops.length === 0;
        const capInput = document.getElementById('split-event-cap');
        capInput.disabled = false;
        capInput.value = String(state.maxMergeEvents);
    }

    // Put the thumb back on `threshold` under whatever stop list now exists.
    function restoreThreshold(threshold) {
        snapSliderToStop(Number.isFinite(threshold)
            ? nearestStopForThreshold(threshold)
            : state.sliderModel.stops.find(stop => stop.threshold_value === null));
    }

    // A new bundle, or a session restore setting the cap. Replays the merge order, which
    // is the expensive half, then selects against the current window.
    //
    // A session restore uses this directly: the threshold_value field that follows it in
    // the registry does the snapping, so snapping here would be undone a moment later.
    function rebuildMergeSeries(maxMergeEvents) {
        state.maxMergeEvents = normalizeMaxMergeEvents(maxMergeEvents);
        const graph = state.bundle.graph;
        state.series = deriveMergeSeries(
            graph.nodes.length, graph.mst_edges || [], graph.merge_impact_metric);
        applyMergeSelection(null);
    }

    // The user moved the cap. The stops are re-laid-out, so the old selection's slider
    // position means nothing; the threshold *value* is what survives -- exactly as it
    // does across a session save, and for the same reason. No replay here: the cap only
    // changes which of the already-derived events are selected.
    function applyMaxMergeEvents(maxMergeEvents) {
        if (!state.bundle || !state.series) { return; }
        const previousThreshold = selectedThresholdValue();
        state.maxMergeEvents = normalizeMaxMergeEvents(maxMergeEvents);
        applyMergeSelection(previousThreshold);
        restoreThreshold(previousThreshold);
        updateSplitEventCount();
        updateThresholdUI(false);
    }

    // How much of the network's split-event series the chart is actually showing.
    // Without this the blank stretches of a large network's axis are unreadable: the
    // axis and the moving sum span every merge, but the stems are a capped selection,
    // so "nothing plotted here" and "nothing happens here" look identical.
    //
    // The series is derived on load from the bundle's merge order, so the denominator
    // is always known -- on a bundle of any version, including the pre-v6 files whose
    // stored series was already capped when it was written.
    function updateSplitEventCount() {
        const note = document.getElementById('split-event-count');
        if (!state.bundle || !state.series) {
            note.textContent = '';
            note.removeAttribute('title');
            return;
        }
        const plotted = state.series.selectedRows.length;
        const total = state.series.total;
        const cap = state.series.cap;
        const bandPercent = 100 / MERGE_EVENT_DENSITY_BINS;
        const plural = count => (count === 1 ? '' : 's');

        if (plotted >= total) {
            note.textContent = 'All ' + total.toLocaleString() + ' merge event' + plural(total) + ' plotted.';
        } else {
            const backfilled = cap === 0 ? 0 : Math.max(0, plotted - cap);
            note.textContent = plotted.toLocaleString() + ' of ' + total.toLocaleString() +
                ' merge events plotted' +
                (backfilled > 0
                    ? ' — the strongest ' + cap.toLocaleString() + ' by impact, plus ' +
                      backfilled.toLocaleString() + ' so that every ' + bandPercent + '% of the axis'
                      + ' with an event to show has one.'
                    : ' — the strongest by impact.');
        }
        note.title = plotted >= total
            ? 'Every split event in this network is drawn.'
            : 'The split chart draws a capped selection of the network’s ' +
              total.toLocaleString() + ' split events, set by the "Split events" box ' +
              '(0 plots them all). The axis and the moving-sum line always span every ' +
              'event, so a stretch with no stems is a stretch whose events were too ' +
              'small to make the cut, not necessarily a quiet one.';
    }

    function scheduleThresholdUI(resetView = true) {
        if (!state.bundle) {
            return;
        }
        state.pendingThresholdUIResetView = state.pendingThresholdUIResetView || resetView;
        if (state.pendingThresholdUIFrame !== null) {
            return;
        }
        state.pendingThresholdUIFrame = window.requestAnimationFrame(() => {
            const nextResetView = state.pendingThresholdUIResetView;
            state.pendingThresholdUIFrame = null;
            state.pendingThresholdUIResetView = false;
            updateThresholdUI(nextResetView);
        });
    }

    function metadataSortIndicator(columnKey) {
        if (state.metadataSort.columnKey !== columnKey) {
            return '↕';
        }
        return state.metadataSort.direction === 'asc' ? '↑' : '↓';
    }

    function metadataSortDescription() {
        if (!state.metadataSort.columnKey) {
            return null;
        }
        return state.metadataSort.columnKey + ' (' + (state.metadataSort.direction === 'asc' ? 'ascending' : 'descending') + ')';
    }

    function metadataFilterText() {
        const input = document.getElementById('metadata-filter');
        return input ? input.value.trim().toLowerCase() : '';
    }

    function metadataNullPlacement() {
        const select = document.getElementById('metadata-null-order');
        return select ? (select.value || 'last') : 'last';
    }

    function isMissingMetadataValue(value) {
        return value === null || value === undefined || value === '';
    }

    function metadataSortValue(nodeIndex, columnKey) {
        if (columnKey === 'node_id') {
            return nodeId(nodeIndex);
        }
        const columnIndex = state.metadataColumnIndexByName.get(columnKey);
        if (columnIndex === undefined) {
            return null;
        }
        return state.metadataByNodeIndex[nodeIndex]?.[columnIndex] ?? null;
    }

    function metadataMatchesFilter(nodeIndex, filterText) {
        if (!filterText) {
            return true;
        }
        return (state.metadataSearchTextByNodeIndex[nodeIndex] || '').includes(filterText);
    }

    function scheduleMetadataFilterUpdate(delayMs = 120) {
        if (state.pendingMetadataFilterTimer !== null) {
            window.clearTimeout(state.pendingMetadataFilterTimer);
        }
        state.pendingMetadataFilterTimer = window.setTimeout(() => {
            state.pendingMetadataFilterTimer = null;
            resetMetadataPage();
            updateMetadataTable();
        }, delayMs);
    }

    function filteredMetadataNodeIndices(nodeIndices) {
        const filterText = metadataFilterText();
        if (!filterText) {
            return [...nodeIndices];
        }
        return nodeIndices.filter(nodeIndex => metadataMatchesFilter(nodeIndex, filterText));
    }

    function metadataDisplayNodeIndices(nodeIndices) {
        return sortedMetadataNodeIndices(filteredMetadataNodeIndices(nodeIndices));
    }

    function metadataRowsPerPageSetting() {
        const select = document.getElementById('metadata-rows-per-page');
        const value = select ? select.value : String(DEFAULT_METADATA_PAGE_SIZE);
        if (value === 'all') {
            return {showAll: true, pageSize: Number.POSITIVE_INFINITY};
        }
        const pageSize = Number.parseInt(value, 10);
        if (!Number.isFinite(pageSize) || pageSize <= 0) {
            return {showAll: false, pageSize: DEFAULT_METADATA_PAGE_SIZE};
        }
        return {showAll: false, pageSize};
    }

    function metadataPagination(totalRowCount) {
        const rowsPerPage = metadataRowsPerPageSetting();
        if (rowsPerPage.showAll) {
            return {pageCount: 1, pageIndex: 0, start: 0, end: totalRowCount, showAll: true};
        }
        const pageCount = Math.max(1, Math.ceil(totalRowCount / rowsPerPage.pageSize));
        const pageIndex = Math.max(0, Math.min(state.metadataPage, pageCount - 1));
        const start = totalRowCount === 0 ? 0 : pageIndex * rowsPerPage.pageSize;
        const end = Math.min(totalRowCount, start + rowsPerPage.pageSize);
        return {pageCount, pageIndex, start, end, showAll: false};
    }

    function resetMetadataPage() {
        state.metadataPage = 0;
    }

    function stepMetadataPage(delta) {
        if (!state.bundle) {
            return;
        }
        const totalRowCount = metadataDisplayNodeIndices(metadataBaseNodeIndices()).length;
        const pagination = metadataPagination(totalRowCount);
        const nextPage = Math.max(0, Math.min(pagination.pageCount - 1, pagination.pageIndex + delta));
        if (nextPage === pagination.pageIndex) {
            return;
        }
        state.metadataPage = nextPage;
        updateMetadataTable();
    }

    function metadataBaseNodeIndices() {
        const selected = Array.from(state.selectedNodeIndices).sort((left, right) => left - right);
        if (selected.length > 0) {
            return selected;
        }
        return state.allNodeIndices;
    }

    function metadataColumnValues(nodeIndices, columnKey) {
        if (columnKey === 'node_id') {
            return nodeIndices.map(nodeIndex => nodeId(nodeIndex));
        }
        return nodeIndices.map(nodeIndex => formatMetadataDisplayValue(columnKey, metadataSortValue(nodeIndex, columnKey)));
    }

    async function copyMetadataColumn(columnKey) {
        if (!state.bundle) {
            return;
        }
        const nodeIndices = metadataDisplayNodeIndices(metadataBaseNodeIndices());
        const values = metadataColumnValues(nodeIndices, columnKey);
        try {
            await writeTextToClipboard(values.join('\n'));
            setStatus('Copied column "' + columnKey + '" for ' + values.length.toLocaleString() + ' rows.');
        } catch (error) {
            console.error(error);
            setStatus('Failed to copy column "' + columnKey + '": ' + error.message);
        }
    }

    function pruneMetadataRowSelection(visibleNodeIndices) {
        const visibleSet = new Set(visibleNodeIndices);
        if (state.metadataRowSelectionAnchor !== null && !visibleSet.has(state.metadataRowSelectionAnchor)) {
            state.metadataRowSelectionAnchor = null;
        }
        if (state.selectedMetadataNodeIndices.size === 0) {
            return;
        }
        state.selectedMetadataNodeIndices = new Set(
            Array.from(state.selectedMetadataNodeIndices).filter(nodeIndex => visibleSet.has(nodeIndex))
        );
    }

    function clearMetadataRowSelection() {
        state.selectedMetadataNodeIndices = new Set();
        state.metadataRowSelectionAnchor = null;
    }

    function metadataTableHasActiveTextSelection() {
        const selection = window.getSelection ? window.getSelection() : null;
        if (!selection || selection.rangeCount === 0 || selection.isCollapsed) {
            return false;
        }
        const metadataTable = document.getElementById('metadata-table');
        if (!metadataTable) {
            return false;
        }
        for (let rangeIndex = 0; rangeIndex < selection.rangeCount; rangeIndex += 1) {
            const range = selection.getRangeAt(rangeIndex);
            if (metadataTable.contains(range.commonAncestorContainer)) {
                return true;
            }
        }
        return false;
    }

    function toggleMetadataRowSelection(nodeIndex, renderedNodeIndices, options = {}) {
        const nextSelection = new Set(state.selectedMetadataNodeIndices);
        if (options.range && state.metadataRowSelectionAnchor !== null) {
            const anchorPosition = renderedNodeIndices.indexOf(state.metadataRowSelectionAnchor);
            const currentPosition = renderedNodeIndices.indexOf(nodeIndex);
            if (anchorPosition >= 0 && currentPosition >= 0) {
                const start = Math.min(anchorPosition, currentPosition);
                const end = Math.max(anchorPosition, currentPosition);
                const shouldSelect = !nextSelection.has(nodeIndex);
                for (let position = start; position <= end; position += 1) {
                    const rangeNodeIndex = renderedNodeIndices[position];
                    if (shouldSelect) {
                        nextSelection.add(rangeNodeIndex);
                    } else {
                        nextSelection.delete(rangeNodeIndex);
                    }
                }
            } else if (nextSelection.has(nodeIndex)) {
                nextSelection.delete(nodeIndex);
            } else {
                nextSelection.add(nodeIndex);
            }
        } else if (nextSelection.has(nodeIndex)) {
            nextSelection.delete(nodeIndex);
        } else {
            nextSelection.add(nodeIndex);
        }
        state.selectedMetadataNodeIndices = nextSelection;
        state.metadataRowSelectionAnchor = nodeIndex;
        applyMetadataTableRowHighlights();
    }

    function selectNodesFromMetadataRows() {
        if (state.selectedMetadataNodeIndices.size === 0) {
            return;
        }
        state.selectedNodeIndices = new Set(state.selectedMetadataNodeIndices);
        resetMetadataPage();
        clearMetadataRowSelection();
        renderClusterView();
        updateMetadataTable();
    }

    function metadataColumnType(columnKey) {
        if (columnKey === 'node_id') {
            return 'string';
        }
        return metadataColumn(columnKey)?.type || 'string';
    }

    function formatMetadataDisplayValue(columnKey, value) {
        if (isMissingMetadataValue(value)) {
            return '—';
        }
        const columnType = metadataColumnType(columnKey);
        if (columnType === 'int') {
            const numeric = typeof value === 'number' ? value : Number(value);
            return Number.isFinite(numeric) ? Math.trunc(numeric).toLocaleString() : String(value);
        }
        if (columnType === 'float') {
            const numeric = typeof value === 'number' ? value : Number(value);
            if (!Number.isFinite(numeric)) {
                return String(value);
            }
            const absValue = Math.abs(numeric);
            if (absValue >= 10000 || (absValue > 0 && absValue < 0.001)) {
                return numeric.toExponential(3);
            }
            return new Intl.NumberFormat(undefined, {maximumFractionDigits: 4}).format(numeric);
        }
        if (columnType === 'bool' || columnType === 'boolean') {
            return value ? 'True' : 'False';
        }
        if (Array.isArray(value)) {
            return value.join(', ');
        }
        return String(value);
    }

    function metadataCellClass(columnKey, value) {
        if (isMissingMetadataValue(value)) {
            return 'metadata-cell-null';
        }
        const columnType = metadataColumnType(columnKey);
        return columnType === 'int' || columnType === 'float' ? 'metadata-cell-number' : '';
    }

    function compareMetadataValues(leftValue, rightValue) {
        if (leftValue === rightValue) {
            return 0;
        }
        if (typeof leftValue === 'number' && typeof rightValue === 'number') {
            return leftValue - rightValue;
        }
        const leftText = String(leftValue);
        const rightText = String(rightValue);
        const leftNumeric = Number(leftText);
        const rightNumeric = Number(rightText);
        if (Number.isFinite(leftNumeric) && Number.isFinite(rightNumeric) && leftText.trim() !== '' && rightText.trim() !== '') {
            return leftNumeric - rightNumeric;
        }
        return leftText.localeCompare(rightText, undefined, {numeric: true, sensitivity: 'base'});
    }

    function sortedMetadataNodeIndices(nodeIndices) {
        if (!state.metadataSort.columnKey) {
            return [...nodeIndices];
        }
        const direction = state.metadataSort.direction === 'desc' ? -1 : 1;
        const nullPlacement = metadataNullPlacement();
        return [...nodeIndices].sort((leftIndex, rightIndex) => {
            const leftValue = metadataSortValue(leftIndex, state.metadataSort.columnKey);
            const rightValue = metadataSortValue(rightIndex, state.metadataSort.columnKey);
            const leftMissing = isMissingMetadataValue(leftValue);
            const rightMissing = isMissingMetadataValue(rightValue);
            if (leftMissing || rightMissing) {
                if (leftMissing && rightMissing) {
                    return nodeId(leftIndex).localeCompare(nodeId(rightIndex), undefined, {numeric: true, sensitivity: 'base'});
                }
                return leftMissing
                    ? (nullPlacement === 'first' ? -1 : 1)
                    : (nullPlacement === 'first' ? 1 : -1);
            }
            const comparison = compareMetadataValues(
                leftValue,
                rightValue,
            );
            if (comparison !== 0) {
                return comparison * direction;
            }
            return nodeId(leftIndex).localeCompare(nodeId(rightIndex), undefined, {numeric: true, sensitivity: 'base'});
        });
    }

    window.addEventListener('pointermove', event => {
        if (!state.metadataResize) {
            return;
        }
        updateMetadataColumnWidth(
            state.metadataResize.columnKey,
            state.metadataResize.startWidth + (event.clientX - state.metadataResize.startX),
        );
    });

    window.addEventListener('pointerup', () => {
        if (!state.metadataResize) {
            return;
        }
        state.metadataResize = null;
        document.body.style.cursor = '';
        document.body.style.userSelect = '';
    });

    window.addEventListener('pointercancel', () => {
        if (!state.metadataResize) {
            return;
        }
        state.metadataResize = null;
        document.body.style.cursor = '';
        document.body.style.userSelect = '';
    });

        document.getElementById('metadata-select-nodes').addEventListener('click', selectNodesFromMetadataRows);

    function toggleMetadataSort(columnKey) {
        if (state.metadataSort.columnKey === columnKey) {
            state.metadataSort.direction = state.metadataSort.direction === 'asc' ? 'desc' : 'asc';
        } else {
            state.metadataSort = {columnKey, direction: 'asc'};
        }
        resetMetadataPage();
        updateMetadataTable();
    }

    function resetMetadataSort() {
        state.metadataSort = {columnKey: null, direction: 'asc'};
        resetMetadataPage();
        updateMetadataTable();
    }

    function jumpToThresholdValue() {
        if (!state.bundle) {
            return;
        }
        const input = document.getElementById('threshold-input');
        const targetValue = Number(input.value);
        if (!Number.isFinite(targetValue)) {
            const stop = currentSliderStop();
            input.value = stop && stop.threshold_value !== null ? String(stop.threshold_value) : '';
            return;
        }
        const stop = nearestStopForThreshold(targetValue);
        if (!stop) {
            return;
        }
        snapSliderToStop(stop);
        updateThresholdUI();
    }

    // The next stop in `delta`'s direction that actually names a different threshold,
    // or -1 when there is none. What an arrow press promises is a threshold the view
    // has not just been at; two stops naming one threshold describe the same cut, so
    // landing on the second would move the slider and change nothing on screen. The
    // event rows the stops are built from are grouped by exact threshold today, so this
    // skips nothing -- it is here because the promise belongs to the arrows rather than
    // to an invariant held one module away in ssn_hierarchy.
    //
    // (The arrows' own stuck-looking behavior came from somewhere else: many stops
    // share one slider position. See currentSliderStop.)
    function nextDistinctStopIndex(stops, index, delta) {
        const current = stops[index].threshold_value;
        for (let probe = index + delta; probe >= 0 && probe < stops.length; probe += delta) {
            if (stops[probe].threshold_value !== current) {
                return probe;
            }
        }
        return -1;
    }

    // Walk the slider's stop list. Those stops are the only thresholds the view can
    // actually take -- every one is an edge weight at which the graph splits -- so
    // stepping to the next distinct one is the only step size that always lands
    // somewhere new. `delta` is +1 toward higher thresholds (rightward on the split
    // plot, ending at the infinity stop) and -1 toward lower ones.
    function stepThreshold(delta) {
        if (!state.bundle || !state.sliderModel) {
            return;
        }
        const stops = state.sliderModel.stops;
        if (stops.length === 0) {
            return;
        }
        // currentSliderStop() returns an element of `stops`, so identity search is safe.
        const index = stops.indexOf(currentSliderStop());
        const target = index < 0 ? 0 : nextDistinctStopIndex(stops, index, delta);
        if (target < 0) {
            return;
        }
        snapSliderToStop(stops[target]);
        // Keep the user's pan/zoom, like releasing the slider does; stepping through
        // stops to watch one cluster break up is the point of these buttons.
        scheduleThresholdUI(false);
    }

    // The ends of the stop list are dead ends, so say so rather than no-op silently.
    // Asked the same question stepThreshold answers, so an arrow is enabled exactly
    // when pressing it would move the threshold somewhere new.
    function updateThresholdStepButtons() {
        const stops = state.sliderModel ? state.sliderModel.stops : [];
        const index = state.bundle && stops.length > 0 ? stops.indexOf(currentSliderStop()) : -1;
        document.getElementById('threshold-step-down').disabled =
            index < 0 || nextDistinctStopIndex(stops, index, -1) < 0;
        document.getElementById('threshold-step-up').disabled =
            index < 0 || nextDistinctStopIndex(stops, index, 1) < 0;
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
            const component = state.bundle.graph.hierarchy.nodes[item.componentId];
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

    function normalizedSelectionBox(box) {
        return {
            left: Math.min(box.startX, box.endX),
            top: Math.min(box.startY, box.endY),
            width: Math.abs(box.endX - box.startX),
            height: Math.abs(box.endY - box.startY),
        };
    }

    function circleIntersectsRect(cx, cy, radius, rect) {
        const nearestX = Math.max(rect.left, Math.min(cx, rect.left + rect.width));
        const nearestY = Math.max(rect.top, Math.min(cy, rect.top + rect.height));
        const dx = cx - nearestX;
        const dy = cy - nearestY;
        return (dx * dx) + (dy * dy) <= radius * radius;
    }

    // Axis-aligned overlap between a screen-space rectangle [ax0,ay0]-[ax1,ay1] (corners in
    // any order) and a selection box {left, top, width, height}.
    function rectIntersectsRect(ax0, ay0, ax1, ay1, rect) {
        const left = Math.min(ax0, ax1);
        const right = Math.max(ax0, ax1);
        const top = Math.min(ay0, ay1);
        const bottom = Math.max(ay0, ay1);
        return left <= (rect.left + rect.width) && right >= rect.left &&
            top <= (rect.top + rect.height) && bottom >= rect.top;
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
            const component = state.bundle.graph.hierarchy.nodes[item.componentId];
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

    document.getElementById('bundle-file').addEventListener('change', async event => {
        const file = event.target.files?.[0];
        if (!file) {
            return;
        }
        try {
            const bundle = await decodeBundleFile(file);
            installBundle(bundle);
        } catch (error) {
            console.error(error);
            setStatus('Failed to load bundle: ' + error.message);
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
    document.getElementById('color-by').addEventListener('change', () => {
        rebuildNodeColorCache();
        renderClusterView();
        if (colorPickerIsOpen()) {
            openColorPicker();
        }
    });
    document.getElementById('label-by').addEventListener('change', () => renderClusterView());
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
    document.getElementById('export-selected').addEventListener('click', exportSelection);
    document.getElementById('export-png').addEventListener('click', exportClusterPNG);
    document.getElementById('export-split-png').addEventListener('click', exportSplitChartPNG);
    document.getElementById('export-split-svg').addEventListener('click', exportSplitChartSVG);
    document.getElementById('export-svg').addEventListener('click', exportClusterSVG);
    document.getElementById('customize-colors').addEventListener('click', openColorPicker);
    document.getElementById('save-session').addEventListener('click', saveSessionFile);
    document.getElementById('rename-network').addEventListener('click', openRenameNetworkDialog);
    document.getElementById('save-extraction').addEventListener('click', saveExtractionFile);
    document.getElementById('color-picker-close').addEventListener('click', closeColorPicker);
    document.getElementById('color-picker-overlay').addEventListener('click', event => {
        if (event.target === event.currentTarget) {
            closeColorPicker();
        }
    });
    document.addEventListener('keydown', event => {
        if (event.key === 'Escape' && colorPickerIsOpen()) {
            closeColorPicker();
        }
    });
    document.getElementById('color-palette').addEventListener('change', changeColorPalette);
    document.getElementById('color-palette-offset').addEventListener('change', changeColorPalette);
    document.getElementById('color-as-categorical').addEventListener('change', toggleColorCategorical);
    document.getElementById('color-picker-reset').addEventListener('click', resetColorPicker);
    document.getElementById('color-picker-prev').addEventListener('click', () => changeColorPickerPage(-1));
    document.getElementById('color-picker-next').addEventListener('click', () => changeColorPickerPage(1));
    document.getElementById('load-color-table').addEventListener('click', () => document.getElementById('color-table-file').click());
    document.getElementById('color-table-file').addEventListener('change', async event => {
        const file = event.target.files && event.target.files[0];
        if (!file) {
            return;
        }
        try {
            const text = await file.text();
            loadColorTableText(text);
        } catch (error) {
            console.error(error);
            setStatus('Failed to load color table: ' + error.message);
        } finally {
            event.target.value = '';
        }
    });
    document.getElementById('save-color-table').addEventListener('click', saveColorTableTSV);
    document.getElementById('export-legend-svg').addEventListener('click', exportLegendSVG);
    document.getElementById('export-legend-png').addEventListener('click', exportLegendPNG);
    document.getElementById('metadata-deselect-rows').addEventListener('click', () => {
        clearMetadataRowSelection();
        applyMetadataTableRowHighlights();
    });
    document.getElementById('metadata-reset-sort').addEventListener('click', resetMetadataSort);
    document.getElementById('metadata-filter').addEventListener('input', () => {
        scheduleMetadataFilterUpdate();
    });
    document.getElementById('metadata-null-order').addEventListener('change', () => {
        resetMetadataPage();
        updateMetadataTable();
    });
    document.getElementById('metadata-rows-per-page').addEventListener('change', () => {
        resetMetadataPage();
        updateMetadataTable();
    });
    document.getElementById('metadata-prev-page').addEventListener('click', () => {
        stepMetadataPage(-1);
    });
    document.getElementById('metadata-next-page').addEventListener('click', () => {
        stepMetadataPage(1);
    });
    document.getElementById('focus-largest-cluster').addEventListener('click', () => {
        if (!state.bundle || state.visibleLayout.length === 0) { return; }
        const hierarchyNodes = state.bundle.graph.hierarchy.nodes;
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
    document.getElementById('focus-selection').addEventListener('click', () => {
        if (!focusSelection()) {
            setStatus('No selected nodes are visible at the current threshold and minimum cluster size.');
            return;
        }
        scheduleClusterRender();
    });
    document.getElementById('zoom-level').addEventListener('change', changeZoomLevel);
    document.getElementById('reset-view').addEventListener('click', () => {
        fitClusterViewToLayout();
        renderClusterView();
    });
    document.getElementById('clear-selection').addEventListener('click', () => {
        state.selectedNodeIndices = new Set();
        resetMetadataPage();
        clearMetadataRowSelection();
        renderClusterView();
        updateMetadataTable();
    });
    window.addEventListener('resize', () => {
        if (!state.bundle) {
            return;
        }
        drawClusterView(true);
    });

    setupLayoutWorker();
    setupMetadataTableDelegation();
    setupSelectionPresets();
    setupMetadataEditing();
    setupColumnCharts();
    setupSplitChartHover();
    setupNameDialog();
    setupGradientRangeSlider();
    setupGradientStopEditing();
    populateColorPaletteMenu();
    warnIfUnsupported();
    updateComponentSortButton();
    drawSplitChart();
    drawClusterView();
    autoloadEmbeddedBundle();
