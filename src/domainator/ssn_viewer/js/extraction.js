/*
 * "Save extraction": a new bundle over just the selected nodes.
 *
 * Included into viewer.js's script scope; see session_state.js.
 * The hierarchy and merge-series maths here are ports of `ssn_hierarchy.py`
 * (`build_mst_component_hierarchy`, `component_size_summary_by_threshold`,
 * `threshold_merge_event_rows`, `merge_event_moving_sum`,
 * `filter_merge_event_rows`) and must stay in step with them.
 */

    // ---------------------------------------------------------------------
    // Save extraction
    //
    // Writes a self-contained bundle over just the selected nodes.
    //
    // Nodes that the MST joins must stay joined in the selection, wherever the
    // slider happens to sit. The MST is a lossy summary: it kept one path
    // between any two nodes and discarded every other edge. So if the selection
    // drops the nodes along that path, the bundle holds no evidence of how the
    // remaining pieces relate, and writing them out as separate components
    // would assert an absence of similarity that only the original all-vs-all
    // matrix could establish. Subset that matrix with
    // domainator build_ssn_viewer --subset instead, which measures the relationship
    // rather than guessing at it.
    //
    // Nodes in *different* components of the original network are exempt: they
    // were already unrelated before the extraction, so keeping them apart
    // invents nothing. That is what lets a whole multi-component network be
    // extracted, and it is checked per original component rather than globally.
    // ---------------------------------------------------------------------
    const MOVING_SUM_WINDOW_FRACTION = 0.05;
    const MOVING_SUM_GRID_POINTS = 800;
    const DEFAULT_MAX_MERGE_EVENTS = 500;
    // The viewer's own default. It picks its events from the chart's current window, so
    // a much smaller number shows more: zoom in and the cap re-spends itself on the
    // stretch being looked at. See ssn_hierarchy.DEFAULT_WINDOWED_MAX_MERGE_EVENTS.
    const DEFAULT_WINDOWED_MAX_MERGE_EVENTS = 50;
    const MERGE_EVENT_DENSITY_BINS = 20;

    function makeUnionFind(size) {
        const parent = new Int32Array(size);
        for (let i = 0; i < size; i++) { parent[i] = i; }
        function find(index) {
            let root = index;
            while (parent[root] !== root) { root = parent[root]; }
            while (parent[index] !== index) {
                const next = parent[index];
                parent[index] = root;
                index = next;
            }
            return root;
        }
        return {find, union: (a, b) => { parent[find(a)] = find(b); }};
    }

    // Port of ssn_hierarchy.format_threshold_value.
    function formatThresholdValue(threshold) {
        return Number.isFinite(threshold) ? Number(threshold).toFixed(2) : '∞';
    }

    // The selected nodes plus the MST edges with both ends inside, re-indexed into
    // the extraction. `pieces` is how many MST-connected components the selection
    // falls into; anything above 1 is what makes an extraction unfaithful.
    function extractionInducedGraph() {
        const nodeIndices = Array.from(state.selectedNodeIndices).sort((left, right) => left - right);
        const newIndexByOld = new Map(nodeIndices.map((oldIndex, newIndex) => [oldIndex, newIndex]));
        const edges = [];
        // Original order is weight-descending, which every consumer below relies on.
        state.bundle.graph.mst_edges.forEach(edge => {
            const left = newIndexByOld.get(edge[0]);
            const right = newIndexByOld.get(edge[1]);
            if (left !== undefined && right !== undefined) {
                edges.push([left, right, edge[2]]);
            }
        });

        const unionFind = makeUnionFind(nodeIndices.length);
        edges.forEach(edge => unionFind.union(edge[0], edge[1]));

        // Split each original network component's share of the selection into the
        // pieces it fell into. A component contributing more than one piece is the
        // failure case: the MST joined those nodes and the selection broke them up.
        const originalRoot = originalComponentByNode();
        const piecesByComponent = new Map();
        nodeIndices.forEach((oldIndex, newIndex) => {
            const component = originalRoot[oldIndex];
            if (!piecesByComponent.has(component)) { piecesByComponent.set(component, new Set()); }
            piecesByComponent.get(component).add(unionFind.find(newIndex));
        });
        let brokenComponents = 0;
        let brokenPieces = 0;
        piecesByComponent.forEach(pieces => {
            if (pieces.size > 1) {
                brokenComponents += 1;
                brokenPieces += pieces.size;
            }
        });
        return {nodeIndices, newIndexByOld, edges, brokenComponents, brokenPieces};
    }

    // Which top-level component of the full MST forest each node belongs to.
    function originalComponentByNode() {
        const hierarchy = state.bundle.graph.hierarchy;
        const componentOf = new Int32Array(state.bundle.graph.nodes.length).fill(-1);
        hierarchy.roots.forEach(rootId => {
            const root = hierarchy.nodes[rootId];
            for (let i = root.leaf_start; i < root.leaf_start + root.leaf_count; i++) {
                componentOf[hierarchy.leaf_order[i]] = rootId;
            }
        });
        return componentOf;
    }

    // Port of ssn_hierarchy.build_mst_component_hierarchy.
    function buildExtractionHierarchy(nodeCount, edges) {
        const unionFind = makeUnionFind(nodeCount);
        const nodes = [];
        const componentIdByRoot = new Map();
        const sizeById = [];
        const minLeafById = [];
        for (let nodeIndex = 0; nodeIndex < nodeCount; nodeIndex++) {
            nodes.push({id: nodeIndex, kind: 'leaf', node_index: nodeIndex, size: 1, parent: null});
            componentIdByRoot.set(nodeIndex, nodeIndex);
            sizeById.push(1);
            minLeafById.push(nodeIndex);
        }

        edges.forEach(edge => {
            const leftRoot = unionFind.find(edge[0]);
            const rightRoot = unionFind.find(edge[1]);
            if (leftRoot === rightRoot) { return; }
            let leftId = componentIdByRoot.get(leftRoot);
            let rightId = componentIdByRoot.get(rightRoot);
            if (minLeafById[leftId] > minLeafById[rightId]) {
                const swap = leftId; leftId = rightId; rightId = swap;
            }
            const componentId = nodes.length;
            nodes[leftId].parent = componentId;
            nodes[rightId].parent = componentId;
            nodes.push({
                id: componentId,
                kind: 'cluster',
                left: leftId,
                right: rightId,
                threshold: edge[2],
                size: sizeById[leftId] + sizeById[rightId],
                parent: null,
            });
            sizeById.push(sizeById[leftId] + sizeById[rightId]);
            minLeafById.push(Math.min(minLeafById[leftId], minLeafById[rightId]));
            unionFind.union(leftRoot, rightRoot);
            componentIdByRoot.delete(leftRoot);
            componentIdByRoot.set(unionFind.find(rightRoot), componentId);
        });

        const roots = Array.from(componentIdByRoot.values())
            .sort((left, right) => nodes[right].size - nodes[left].size || left - right);

        const leafOrder = [];
        const stack = [];
        for (let i = roots.length - 1; i >= 0; i--) { stack.push([roots[i], false]); }
        while (stack.length > 0) {
            const [componentId, visited] = stack.pop();
            const node = nodes[componentId];
            if (node.kind === 'leaf') {
                node.leaf_start = leafOrder.length;
                node.leaf_count = 1;
                leafOrder.push(node.node_index);
                continue;
            }
            if (visited) {
                node.leaf_start = nodes[node.left].leaf_start;
                node.leaf_count = nodes[node.left].leaf_count + nodes[node.right].leaf_count;
                continue;
            }
            stack.push([componentId, true]);
            stack.push([node.right, false]);
            stack.push([node.left, false]);
        }
        return {nodes, roots, leaf_order: leafOrder};
    }

    // Port of ssn_hierarchy.component_size_summary_by_threshold followed by
    // threshold_merge_event_rows: one row per distinct MST edge weight.
    function mergeEventRows(nodeCount, edges, metric) {
        if (nodeCount === 0 || edges.length === 0) { return []; }
        const unionFind = makeUnionFind(nodeCount);
        const componentSizes = new Int32Array(nodeCount).fill(1);
        // Sizes only ever grow, so the running maximum is the largest component.
        let largest = nodeCount > 0 ? 1 : 0;
        let nonSingletonCount = 0;
        let nonSingletonSum = 0;
        const summary = [{threshold: Infinity, largest, avg: 0, impact: 0}];

        edges.forEach(edge => {
            const leftRoot = unionFind.find(edge[0]);
            const rightRoot = unionFind.find(edge[1]);
            let impact = 0;
            if (leftRoot !== rightRoot) {
                const leftSize = componentSizes[leftRoot];
                const rightSize = componentSizes[rightRoot];
                impact = metric === 'product' ? leftSize * rightSize : Math.min(leftSize, rightSize);
                if (leftSize > 1) { nonSingletonCount -= 1; nonSingletonSum -= leftSize; }
                if (rightSize > 1) { nonSingletonCount -= 1; nonSingletonSum -= rightSize; }
                const mergedSize = leftSize + rightSize;
                unionFind.union(leftRoot, rightRoot);
                componentSizes[unionFind.find(rightRoot)] = mergedSize;
                largest = Math.max(largest, mergedSize);
                nonSingletonCount += 1;
                nonSingletonSum += mergedSize;
            }
            summary.push({
                threshold: edge[2],
                largest,
                avg: nonSingletonCount > 0 ? nonSingletonSum / nonSingletonCount : 0,
                impact,
            });
        });

        const rows = [];
        let previous = summary[0];
        let rowIndex = 1;
        while (rowIndex < summary.length) {
            const firstRowIndex = rowIndex;
            const thresholdValue = summary[rowIndex].threshold;
            let mergeImpact = 0;
            let mergeCount = 0;
            let largestMerge = 0;
            const mergeSizeCounts = {};
            let lastRow = summary[rowIndex];
            while (rowIndex < summary.length && summary[rowIndex].threshold === thresholdValue) {
                const impact = summary[rowIndex].impact;
                mergeImpact += impact;
                if (impact > 0) {
                    // A zero impact means both ends were already in one component.
                    const key = Number.isInteger(impact) ? String(impact) : String(impact);
                    mergeSizeCounts[key] = (mergeSizeCounts[key] || 0) + 1;
                    mergeCount += 1;
                    largestMerge = Math.max(largestMerge, impact);
                }
                lastRow = summary[rowIndex];
                rowIndex += 1;
            }
            rows.push({
                edge_index: firstRowIndex - 2,
                summary_row_from: firstRowIndex,
                summary_row_to: rowIndex - 1,
                threshold_from_value: previous.threshold,
                threshold_from: formatThresholdValue(previous.threshold),
                threshold_to: formatThresholdValue(lastRow.threshold),
                threshold_value: lastRow.threshold,
                merge_impact: mergeImpact,
                merge_size_counts: mergeSizeCounts,
                largest_merge: largestMerge,
                merge_count: mergeCount,
                delta_largest: Math.abs(lastRow.largest - previous.largest),
                delta_avg_non_singleton: Math.abs(lastRow.avg - previous.avg),
            });
            previous = lastRow;
        }
        return rows;
    }

    // Port of ssn_hierarchy.merge_event_rank_key: strongest first.
    function compareMergeEventRank(left, right) {
        return right.merge_impact - left.merge_impact ||
            right.delta_largest - left.delta_largest ||
            right.delta_avg_non_singleton - left.delta_avg_non_singleton ||
            left.edge_index - right.edge_index;
    }

    // Port of ssn_hierarchy.merge_event_density_bin.
    function mergeEventDensityBin(thresholdValue, lo, hi, densityBins) {
        if (densityBins < 1 || !Number.isFinite(lo) || !Number.isFinite(hi) || hi <= lo) {
            return -1;
        }
        const position = Math.floor(((Number(thresholdValue) - lo) / (hi - lo)) * densityBins);
        return Math.max(0, Math.min(densityBins - 1, position));
    }

    // Port of ssn_hierarchy.filter_merge_event_rows: the strongest maxMergeEvents events
    // *within `window`*, plus the strongest event in each otherwise-empty band of the
    // whole axis, plus `pinnedThreshold`'s event. See that function for the reasoning;
    // in short, the window is what makes a small cap show more rather than less, the
    // bands keep every stretch of the slider reachable, and the pin keeps the cut
    // currently in effect from being filtered out from under the user.
    function filterMergeEventRows(rows, maxMergeEvents = DEFAULT_MAX_MERGE_EVENTS,
                                  densityBins = MERGE_EVENT_DENSITY_BINS,
                                  window = null, pinnedThreshold = null) {
        if (maxMergeEvents === 0 || rows.length <= maxMergeEvents) { return rows.slice(); }
        const ranked = rows.slice().sort(compareMergeEventRank);

        let pool = ranked;
        if (window) {
            const low = Math.min(window.min, window.max);
            const high = Math.max(window.min, window.max);
            pool = ranked.filter(row => {
                const value = Number(row.threshold_value);
                return value >= low && value <= high;
            });
        }
        const filtered = pool.slice(0, maxMergeEvents);
        const selected = new Set(filtered);

        // Band edges come from the whole series, because that is what the axis spans --
        // a window narrows what the cap ranks, never what the axis has to cover.
        // Scanned rather than spread through Math.min: `rows` is the UNFILTERED series,
        // which on a large network runs to six figures and would overflow the argument
        // list.
        let lo = Infinity;
        let hi = -Infinity;
        let finiteCount = 0;
        for (const row of rows) {
            const value = Number(row.threshold_value);
            if (!Number.isFinite(value)) { continue; }
            finiteCount += 1;
            if (value < lo) { lo = value; }
            if (value > hi) { hi = value; }
        }
        if (finiteCount > 0 && densityBins > 0) {
            const covered = new Set(filtered.map(row =>
                mergeEventDensityBin(row.threshold_value, lo, hi, densityBins)));
            // `ranked` is strongest-first, so the first row seen in an empty band is the
            // strongest one available to represent it. Candidates come from the whole
            // series, not the window: a band outside the window still needs its stop.
            for (const row of ranked) {
                if (selected.has(row)) { continue; }
                const bin = mergeEventDensityBin(row.threshold_value, lo, hi, densityBins);
                if (bin < 0 || covered.has(bin)) { continue; }
                covered.add(bin);
                filtered.push(row);
                selected.add(row);
            }
        }

        if (pinnedThreshold !== null && Number.isFinite(pinnedThreshold)
                && !filtered.some(row => Number(row.threshold_value) === pinnedThreshold)) {
            const pinned = ranked.find(row => Number(row.threshold_value) === pinnedThreshold);
            if (pinned) { filtered.push(pinned); }
        }

        return filtered.sort((left, right) => left.edge_index - right.edge_index);
    }

    // Port of ssn_hierarchy.merge_event_moving_sum. Runs over the UNFILTERED rows:
    // filtering keeps the strongest events, so summing afterwards would undercount
    // exactly the small ones this series exists to show.
    function mergeEventMovingSum(rows) {
        const empty = {window: 0, x: [], y: []};
        if (rows.length === 0) { return empty; }
        const points = rows
            .filter(row => Number.isFinite(row.threshold_value) && Number.isFinite(row.merge_impact))
            .map(row => ({threshold: row.threshold_value, impact: row.merge_impact}))
            .sort((left, right) => left.threshold - right.threshold);
        if (points.length === 0) { return empty; }
        const lo = points[0].threshold;
        const hi = points[points.length - 1].threshold;
        if (hi === lo) { return empty; }

        const window = MOVING_SUM_WINDOW_FRACTION * (hi - lo);
        const halfWindow = window / 2;
        const sortedThresholds = points.map(point => point.threshold);
        const cumulative = [0];
        points.forEach(point => cumulative.push(cumulative[cumulative.length - 1] + point.impact));

        // Inclusive at both ends, matching numpy searchsorted 'left'/'right'.
        const lowerBound = value => {
            let low = 0, high = sortedThresholds.length;
            while (low < high) {
                const mid = (low + high) >> 1;
                if (sortedThresholds[mid] < value) { low = mid + 1; } else { high = mid; }
            }
            return low;
        };
        const upperBound = value => {
            let low = 0, high = sortedThresholds.length;
            while (low < high) {
                const mid = (low + high) >> 1;
                if (sortedThresholds[mid] <= value) { low = mid + 1; } else { high = mid; }
            }
            return low;
        };

        const x = [];
        const y = [];
        for (let i = 0; i < MOVING_SUM_GRID_POINTS; i++) {
            const gridValue = lo + ((hi - lo) * i) / (MOVING_SUM_GRID_POINTS - 1);
            x.push(gridValue);
            y.push(cumulative[upperBound(gridValue + halfWindow)] - cumulative[lowerBound(gridValue - halfWindow)]);
        }
        return {window, x, y};
    }

    // Port of ssn_hierarchy.threshold_slider_stops, including its trailing floor stop:
    // every other stop excludes its own tie group, so without one strictly below the
    // weakest edge the fully merged network cannot be reached on the slider.
    function buildSliderStops(rows, edges) {
        const stops = [{edge_index: -1, threshold_label: '\u221e', threshold_value: null}];
        rows.forEach(row => {
            stops.push({
                edge_index: row.edge_index,
                threshold_label: row.threshold_to,
                threshold_value: row.threshold_value,
            });
        });
        if (edges.length > 0) {
            // edges are weight-descending.
            const floorValue = floorThresholdValue(edges[edges.length - 1][2], edges[0][2]);
            stops.push({
                edge_index: edges.length - 1,   // every MST edge is strictly above
                threshold_label: formatThresholdValue(floorValue),
                threshold_value: floorValue,
            });
        }
        return stops;
    }

    // The window-independent half of the series: the full event list and its moving sum,
    // derived from the merge order the bundle carries. From bundle v6 the file stores
    // none of it -- it is a pure function of (node count, mst_edges, metric), and
    // deriving it here is what lets both the cap and the chart's window decide what is
    // shown, neither of which is knowable when the file is written. Pre-v6 files still
    // carry a capped copy; it is ignored in favor of this, so an old bundle is not
    // stuck with whatever cap it happened to be built with.
    //
    // This is the expensive half -- one union-find replay over the MST edges -- so it
    // runs once per bundle. Zooming re-runs only selectMergeEvents below.
    function deriveMergeSeries(nodeCount, edges, metric) {
        const eventRows = mergeEventRows(nodeCount, edges, metric);
        return {
            edges,
            eventRows,
            // From the UNFILTERED rows, so the chart's axis and the slider's track span
            // every event whatever the cap and the window are.
            movingSum: mergeEventMovingSum(eventRows),
            total: eventRows.length,
            // Filled in by selectMergeEvents, which needs a window this does not have.
            selectedRows: [],
            sliderStops: [],
            cap: 0,
        };
    }

    // The window-dependent half: which of those events the chart plots and the slider
    // offers stops for. Cheap (a sort and a scan), so it re-runs on every zoom and pan.
    function selectMergeEvents(series, maxMergeEvents, window, pinnedThreshold) {
        series.selectedRows = filterMergeEventRows(
            series.eventRows, maxMergeEvents, MERGE_EVENT_DENSITY_BINS, window, pinnedThreshold);
        series.sliderStops = buildSliderStops(series.selectedRows, series.edges);
        series.cap = maxMergeEvents;
        return series;
    }

    // The full threshold range of the series, which is what the chart's axis spans and
    // what a zoom window is clamped against. Taken from the UNFILTERED rows (and the
    // moving sum, which runs over them too) rather than from the current selection --
    // the selection depends on the window, so deriving the range from it would be
    // circular, and would make the axis breathe as the user zoomed.
    function splitSeriesDataRange(series) {
        let min = Infinity;
        let max = -Infinity;
        const consider = value => {
            if (!Number.isFinite(value)) { return; }
            if (value < min) { min = value; }
            if (value > max) { max = value; }
        };
        for (const row of series.eventRows) { consider(Number(row.threshold_value)); }
        for (const value of (series.movingSum.x || [])) { consider(Number(value)); }
        return Number.isFinite(min) ? {min, max} : {min: 0, max: 0};
    }

    // Port of ssn_hierarchy.floor_threshold_value.
    const FLOOR_THRESHOLD_SPAN_FRACTION = 0.01;

    function floorThresholdValue(lowestMstWeight, highestMstWeight) {
        const span = highestMstWeight - lowestMstWeight;
        if (span > 0) {
            return lowestMstWeight - (FLOOR_THRESHOLD_SPAN_FRACTION * span);
        }
        const step = Math.abs(lowestMstWeight) * FLOOR_THRESHOLD_SPAN_FRACTION;
        return lowestMstWeight - (step > 0 ? step : 1);
    }

    // The session state, with every saved node id that did not survive dropped, so
    // the extraction opens without complaining about ids it does not contain.
    function extractionAppState(keptIds) {
        const appState = collectSessionState();
        const keep = ids => (ids || []).filter(id => keptIds.has(id));
        appState.selection.node_ids = keep(appState.selection.node_ids);
        const presets = appState.selection.presets || {};
        Object.keys(presets).forEach(slot => {
            const surviving = keep(presets[slot].node_ids);
            if (surviving.length === 0) {
                delete presets[slot];
            } else {
                presets[slot].node_ids = surviving;
            }
        });
        return appState;
    }

    function buildExtractionBundle() {
        const induced = extractionInducedGraph();
        if (induced.nodeIndices.length === 0) {
            return {error: 'Select some nodes or clusters first.'};
        }
        if (induced.brokenComponents > 0) {
            return {error: 'The selection splits ' + induced.brokenComponents.toLocaleString() +
                ' network component' + (induced.brokenComponents === 1 ? '' : 's') + ' into ' +
                induced.brokenPieces.toLocaleString() + ' pieces by leaving out the nodes that ' +
                'join them in the MST. The bundle records no similarity between those pieces, so ' +
                'extracting them would assert they are unrelated. Select the nodes that link them ' +
                '(raising the threshold shows where the links are), or subset the original matrix ' +
                'with domainator build_ssn_viewer --subset.'};
        }

        const nodeCount = induced.nodeIndices.length;
        const metric = state.bundle.graph.merge_impact_metric;
        const keptIds = new Set(induced.nodeIndices.map(nodeIndex => nodeId(nodeIndex)));

        return {
            bundle: {
                format: BUNDLE_FORMAT,
                version: BUNDLE_VERSION,
                name: ((state.bundle.name || 'network') + '_extraction'),
                domainator_version: DOMAINATOR_VERSION,
                // A v6 graph: only what cannot be derived from the merge order. The
                // split-event series, its moving sum and the slider's stops are all a
                // function of these four keys, so an extraction no longer computes a
                // frozen copy of them -- whoever opens it derives them, at whatever cap
                // they are looking at it with.
                graph: {
                    nodes: induced.nodeIndices.map(nodeIndex => nodeId(nodeIndex)),
                    mst_edges: induced.edges,
                    merge_impact_metric: metric,
                    hierarchy: buildExtractionHierarchy(nodeCount, induced.edges),
                },
                metadata: {
                    columns: state.metadataColumns.map(column => ({...column})),
                    rows: induced.nodeIndices.map(nodeIndex => state.metadataByNodeIndex[nodeIndex].slice()),
                },
                defaults: {
                    color_by: currentColorField() || null,
                    label_by: currentLabelField() === '__node_id__' ? null : (currentLabelField() || null),
                    categorical_columns: Array.from(state.categoricalColumns),
                },
                app_state: extractionAppState(keptIds),
            },
            nodeCount,
        };
    }

    function saveExtractionFile() {
        if (!state.bundle) { return; }
        // Validate before asking for a name, so a selection that cannot be extracted
        // is refused without making the user name it first.
        const built = buildExtractionBundle();
        if (built.error) {
            setStatus(built.error);
            window.alert(built.error);
            return;
        }
        openNameDialog({
            title: 'Save extraction',
            label: 'Name for the extracted network',
            note: 'Names the new bundle and the file it is saved to. The extraction contains '
                + built.nodeCount.toLocaleString() + ' of '
                + state.bundle.graph.nodes.length.toLocaleString() + ' nodes.',
            value: built.bundle.name,
            confirmLabel: 'Save',
            onConfirm: name => {
                built.bundle.name = name;
                downloadBundleFile(
                    built.bundle,
                    name.replace(/\s+/g, '_'),
                    'Saved "' + name + '": an extraction of ' + built.nodeCount.toLocaleString() + ' nodes'
                );
            },
        });
    }
