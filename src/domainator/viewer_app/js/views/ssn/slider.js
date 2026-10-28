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

    function selectedThresholdValue() {
        const stop = currentSliderStop();
        if (!stop) {
            return Infinity;
        }
        return stop.threshold_value === null ? Infinity : stop.threshold_value;
    }

    function activeClustersAtThreshold(thresholdValue) {
        const hierarchy = state.bundle.views.ssn.hierarchy;
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
        const hierarchy = state.bundle.views.ssn.hierarchy;
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
        const leftNode = state.bundle.views.ssn.hierarchy.nodes[leftId];
        const rightNode = state.bundle.views.ssn.hierarchy.nodes[rightId];
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
        const network = state.bundle.views.ssn;
        state.series = deriveMergeSeries(
            state.bundle.nodes.length, network.mst_edges || [], network.merge_impact_metric);
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

