    // The SSN view: the hooks the shared code (js/common/) calls while a .dvb bundle holding a
    // network is open. See VIEWS in viewer.js.
    VIEWS.ssn = {
        opens: bundle => Object.keys(bundle.views)[0] === 'ssn',
        // Used in file names, as in "network_table.tsv".
        bundleNoun: 'network',
        hasClusters: true,
        canvas: clusterCanvas,
        nodeIds: () => state.bundle.nodes,
        repaintNow: renderClusterView,
        scheduleRepaint: scheduleClusterRender,
        // The connected-component number at the current threshold, the column domainator
        // build_ssn --threshold <threshold> --cluster writes.
        exportColumns: () => {
            const clusterNumbers = clusterNumbersAtCurrentThreshold();
            return [{name: 'SSN_cluster', valueFor: nodeIndex => clusterNumbers.numberFor(nodeIndex)}];
        },
        // The default name matches build_ssn's CLUSTER_COLUMN so the two are interchangeable
        // downstream.
        clusterAssignments: () => {
            const clusterNumbers = clusterNumbersAtCurrentThreshold();
            return {
                defaultName: 'SSN_cluster',
                clusterCount: clusterNumbers.clusterCount,
                numberFor: clusterNumbers.numberFor,
                description: () => {
                    const stop = currentSliderStop();
                    return 'at threshold ' + (stop ? stop.threshold_label : '∞');
                },
            };
        },
        metadataNoteText: ({selectedCount, pagerActive, filterDescription, sortDescription}) => {
            return (selectedCount === 0
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
        },
        // Clears what the previous bundle left behind, before the shared state is rebuilt.
        reset: () => {
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
        },
        // Called once the metadata controls are built, before a saved session is applied.
        enableControls: () => {
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
            document.getElementById('threshold-input').disabled = state.sliderModel.stops.length === 0;
            document.getElementById('metadata-add-cluster-column').title =
                'Number every node by the cluster it belongs to at the current threshold, largest cluster first — the same column domainator build_ssn --cluster writes. Defaults to the name SSN_cluster.';
            updateThresholdStepButtons();
        },
        // Controls whose appearance is derived from the values a session just restored.
        afterSessionApplied: () => {
            updateComponentSortButton();
            updateCollapseLongPathsControl();
            // Re-select and re-place under the window that was just restored. The registry
            // order already does this for a session carrying view.max_merge_events, so this
            // is idempotent there; it is load-bearing for an older session that carries a
            // window but no cap, whose selection was otherwise made before the window was.
            repositionSliderStops();
        },
        afterInstalled: () => {
            updateSplitEventCount();
            updateThresholdUI();
        },
        sessionFields: shared => [
            domValueField('view', 'layout_algorithm', 'layout-algorithm'),
            domValueField('view', 'min_cluster_size', 'min-cluster-size'),
            shared.colorBy,
            shared.labelBy,
            shared.exportPngScale,
            domCheckboxField('view', 'show_node_counts', 'show-node-counts'),
            domCheckboxField('view', 'show_edge_scores', 'show-edge-scores'),
            domCheckboxField('view', 'render_cluster_bounds', 'render-cluster-bounds'),
            domCheckboxField('view', 'render_nodes', 'render-nodes'),
            domCheckboxField('view', 'reduce_elongation', 'reduce-elongation'),
            domCheckboxField('view', 'leaf_pruning_only', 'leaf-pruning-only'),
            domCheckboxField('view', 'collapse_long_paths', 'collapse-long-paths'),
            {
                section: 'view',
                key: 'sort_components_by_size',
                get: () => sortComponentsBySizeEnabled(),
                set: value => {
                    document.getElementById('sort-components-by-size')
                        .setAttribute('aria-pressed', value ? 'true' : 'false');
                },
            },
            {
                // The split chart's x window, saved for the same reason the canvas
                // transform below is: reopening a session onto the whole axis would throw
                // away the part of it the session was about. Stored in threshold units, so
                // it survives any change of cap; splitChartLayout clamps it into whatever
                // range the bundle actually has.
                //
                // Restored before max_merge_events, because the window is what that cap is
                // spent on: selecting first and setting the window afterwards would leave a
                // zoomed session showing the unzoomed selection.
                section: 'view',
                key: 'split_chart_zoom',
                get: () => state.splitChartZoom
                    ? {min: state.splitChartZoom.min, max: state.splitChartZoom.max}
                    : null,
                set: value => {
                    if (value === null || value === undefined) {
                        state.splitChartZoom = null;
                        return;
                    }
                    if (!Number.isFinite(value.min) || !Number.isFinite(value.max) || !(value.max > value.min)) {
                        throw new Error('malformed split chart zoom');
                    }
                    state.splitChartZoom = {min: Number(value.min), max: Number(value.max)};
                },
            },
            {
                // Restored after split_chart_zoom (the window it is spent on) and before
                // threshold_value, because it decides which stops exist for that field to
                // snap to. Uses rebuildMergeSeries rather than applyMaxMergeEvents for the
                // same reason: snapping here would be undone.
                section: 'view',
                key: 'max_merge_events',
                get: () => state.maxMergeEvents,
                set: value => {
                    const cap = Number(value);
                    if (!Number.isFinite(cap) || cap < 0) { throw new Error('not a count: ' + value); }
                    rebuildMergeSeries(cap);
                },
            },
            {
                // The threshold is stored as its value, not as the slider position:
                // positions are derived from the slider's stops and shift whenever the
                // split-event cap changes.
                section: 'view',
                key: 'threshold_value',
                get: () => {
                    const value = selectedThresholdValue();
                    return Number.isFinite(value) ? value : null;
                },
                set: value => {
                    if (value === null || value === undefined) {
                        const infinityStop = state.sliderModel.stops.find(stop => stop.threshold_value === null);
                        snapSliderToStop(infinityStop);
                        return;
                    }
                    const target = Number(value);
                    if (!Number.isFinite(target)) { throw new Error('not a number: ' + value); }
                    snapSliderToStop(nearestStopForThreshold(target));
                },
            },
            shared.viewTransform,
            shared.filter,
            shared.nullOrder,
            shared.rowsPerPage,
            shared.sort,
            shared.columnWidths,
            shared.customPalettes,
            shared.categoricalColumns,
            shared.nodeIds,
            shared.presets,
        ],
        exportPng: exportClusterPNG,
        exportSvg: exportClusterSVG,
        focusSelection: () => {
            if (!focusSelection()) {
                setStatus('No selected nodes are visible at the current threshold and minimum cluster size.');
                return;
            }
            scheduleClusterRender();
        },
        resetView: () => {
            fitClusterViewToLayout();
            renderClusterView();
        },
        updateSelectionControls: selectedCount => {
            document.getElementById('focus-selection').disabled = selectedCount === 0;
            document.getElementById('save-extraction').disabled = !state.bundle || selectedCount === 0;
        },
    };
