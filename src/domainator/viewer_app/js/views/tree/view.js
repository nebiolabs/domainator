    // The tree view: the hooks the shared code (js/common/) calls while a .dvb bundle holding a
    // tree is open. See VIEWS in viewer.js.
    VIEWS.tree = {
        opens: bundle => Object.keys(bundle.views)[0] === 'tree',
        bundleNoun: 'tree',
        hasClusters: true,
        canvas: document.getElementById('tree-canvas'),
        nodeIds: () => state.bundle.nodes,
        repaintNow: renderTreeView,
        scheduleRepaint: scheduleTreeRender,
        // The clade each leaf is in, cutting the tree at the slider.
        exportColumns: () => {
            const clusterNumbers = treeClusterNumbers();
            return [{name: 'tree_cluster', valueFor: clusterNumbers.numberFor}];
        },
        clusterAssignments: () => {
            const clusterNumbers = treeClusterNumbers();
            return {
                defaultName: 'tree_cluster',
                clusterCount: clusterNumbers.clusterCount,
                numberFor: clusterNumbers.numberFor,
                description: () => {
                    const cutX = treeCutX();
                    return Number.isFinite(cutX) ? 'cutting the tree at ' + formatTreeLength(cutX) : 'with no cut (one per leaf)';
                },
            };
        },
        metadataNoteText: ({selectedCount, pagerActive, filterDescription, sortDescription}) => {
            const sorted = sortDescription ? ', sorted by ' + sortDescription : '';
            if (selectedCount === 0) {
                return (pagerActive ? 'No leaves selected. Browse the full table with the pager' : 'No leaves selected. Showing the full table') +
                    filterDescription + sorted + '.';
            }
            return pagerActive
                ? 'Showing one page of the selected rows in the table' + filterDescription + sorted + '. Use the pager to browse more rows. Export includes the full filtered selection.'
                : 'Click a leaf to toggle it, click a branch to toggle its clade, Shift-drag a box to add leaves, click table rows to stage them, shift-click to select or deselect row ranges, use the search box to filter rows, and click a column header to sort' +
                    (sortDescription ? ' by ' + sortDescription : '') + '.';
        },
        reset: bundle => {
            if (state.tree && state.tree.pendingFrame !== null) {
                cancelAnimationFrame(state.tree.pendingFrame);
            }
            state.tree = createTreeState(bundle);
        },
        enableControls: () => {
            const tree = state.tree;
            ['reset-view', 'export-png', 'export-png-scale', 'export-svg', 'tree-row-height', 'tree-cut',
                'tree-expand-all', 'export-newick', 'export-fasta'].forEach(id => {
                document.getElementById(id).disabled = false;
            });
            document.querySelectorAll('.tree-msa-control').forEach(element => { element.hidden = !tree.msa; });
            document.getElementById('msa-pane').hidden = !tree.msa;
            document.getElementById('tree-show-support').disabled = !tree.hasSupport;
            document.getElementById('tree-cut').value = String(TREE_CUT_DEFAULT);
            document.getElementById('metadata-add-cluster-column').title =
                'Number every leaf by the clade it is in when the tree is cut at the "Cut tree at" slider, largest clade first. Defaults to the name tree_cluster.';
        },
        afterSessionApplied: () => {},
        afterInstalled: () => {
            // A small tree gets a stage just tall enough for its rows at the default height.
            const rows = state.tree.leafCount[0];
            const wanted = rows * TREE_ROW_HEIGHT_DEFAULT + (state.tree.msa ? MSA_HEADER : TREE_HEADER_PLAIN) + TREE_PADDING_BOTTOM + 2;
            document.querySelector('.tree-stage').style.height = wanted < window.innerHeight * 0.78 ? Math.max(320, wanted) + 'px' : '';
            treeComputeLayout();
            updateTreeCutLabel();
            renderTreeView();
        },
        sessionFields: shared => [
            domValueField('view', 'tree_layout', 'tree-layout'),
            domCheckboxField('view', 'tree_ladderize', 'tree-ladderize'),
            domCheckboxField('view', 'tree_show_labels', 'tree-show-labels'),
            domCheckboxField('view', 'tree_show_color_bars', 'tree-show-color-bars'),
            domCheckboxField('view', 'tree_color_branches', 'tree-color-branches'),
            domCheckboxField('view', 'tree_show_support', 'tree-show-support'),
            domValueField('view', 'tree_cut', 'tree-cut'),
            domValueField('view', 'msa_scheme', 'msa-scheme'),
            domValueField('view', 'msa_column_width', 'msa-column-width'),
            shared.colorBy,
            shared.labelBy,
            shared.exportPngScale,
            shared.exportExtent,
            {
                section: 'view',
                key: 'tree_rooting',
                get: () => ({original: treeTopology(state.tree.rootOriginal), location: state.tree.rootLocation}),
                set: restoreTreeRootingState,
            },
            {
                section: 'view',
                key: 'tree_child_order',
                get: () => Array.from(state.tree.childOrder, ([node, children]) => ({node, children: children.slice()}))
                    .sort((left, right) => left.node - right.node),
                set: restoreTreeChildOrder,
            },
            {
                // Collapsed clades by their node's place in the tree's arrays, which a session
                // file carries along with the clades.
                section: 'view',
                key: 'tree_collapsed',
                get: () => Array.from(state.tree.collapsed).sort((left, right) => left - right),
                set: value => {
                    const nodes = (value || []).map(Number);
                    if (nodes.some(node => !Number.isInteger(node) || node <= 0 || node >= state.tree.count || treeIsLeaf(node))) {
                        throw new Error('not internal nodes of this tree');
                    }
                    state.tree.collapsed = new Set(nodes);
                },
            },
            {
                section: 'view',
                key: 'tree_scroll',
                get: () => ({
                    row_height: state.tree.rowHeight,
                    scroll_y: state.tree.scrollY,
                    msa_scroll_x: state.tree.msaScrollX,
                }),
                set: value => {
                    if (!value || !(value.row_height > 0) || !Number.isFinite(value.scroll_y) || !Number.isFinite(value.msa_scroll_x)) {
                        throw new Error('malformed tree scroll');
                    }
                    state.tree.rowHeight = Number(value.row_height);
                    state.tree.scrollY = Number(value.scroll_y);
                    state.tree.msaScrollX = Number(value.msa_scroll_x);
                },
            },
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
        exportPng: exportTreePNG,
        exportSvg: exportTreeSVG,
        focusSelection: treeFocusSelection,
        resetView: treeResetView,
        updateSelectionControls: selectedCount => {
            document.getElementById('focus-selection').disabled = selectedCount === 0;
            updateTreeRootControls();
            updateTreeSwapControl();
        },
    };
