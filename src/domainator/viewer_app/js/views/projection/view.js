    // The projection view: the hooks the shared code (js/common/) calls while a .dvb bundle
    // holding a projection is open. See VIEWS in viewer.js.
    VIEWS.projection = {
        opens: bundle => Object.keys(bundle.views)[0] === 'projection',
        bundleNoun: 'projection',
        hasClusters: false,
        canvas: document.getElementById('projection-canvas'),
        nodeIds: () => state.bundle.nodes,
        repaintNow: renderProjectionView,
        scheduleRepaint: scheduleProjectionRender,
        exportColumns: () => [],
        clusterAssignments: () => null,
        metadataNoteText: ({selectedCount, pagerActive, filterDescription, sortDescription}) => {
            const sorted = sortDescription ? ', sorted by ' + sortDescription : '';
            if (selectedCount === 0) {
                return (pagerActive ? 'No points selected. Browse the full table with the pager' : 'No points selected. Showing the full table') +
                    filterDescription + sorted + '.';
            }
            return pagerActive
                ? 'Showing one page of the selected rows in the table' + filterDescription + sorted + '. Use the pager to browse more rows. Export includes the full filtered selection.'
                : 'Click a point to toggle it, Shift-drag a box to add points, click table rows to stage them, shift-click to select or deselect row ranges, use the search box to filter rows, and click a column header to sort' +
                    (sortDescription ? ' by ' + sortDescription : '') + '.';
        },
        reset: bundle => {
            if (state.projection && state.projection.pendingFrame !== null) {
                cancelAnimationFrame(state.projection.pendingFrame);
            }
            state.projection = createProjectionState(bundle);
            document.getElementById('projection-point-size').value = String(defaultProjectionPointRadius(bundle.nodes.length));
        },
        enableControls: () => {
            ['reset-view', 'export-png', 'export-png-scale', 'export-svg'].forEach(id => {
                document.getElementById(id).disabled = false;
            });
            document.querySelector('.projection-3d-control').hidden = state.projection.view.axes.length !== 3;
        },
        afterSessionApplied: () => {},
        afterInstalled: renderProjectionView,
        sessionFields: shared => [
            domValueField('view', 'projection_axes', 'projection-axes'),
            domValueField('view', 'projection_point_size', 'projection-point-size'),
            shared.colorBy,
            shared.labelBy,
            shared.exportPngScale,
            {
                // The data at the plot's center and the data units per pixel.
                section: 'view',
                key: 'projection_view',
                get: () => state.projection.center === null ? null : {
                    center_x: state.projection.center.x,
                    center_y: state.projection.center.y,
                    units_per_pixel: state.projection.unitsPerPixel,
                },
                set: value => {
                    if (value === null || value === undefined) {
                        state.projection.center = null;
                        return;
                    }
                    if (!Number.isFinite(value.center_x) || !Number.isFinite(value.center_y) || !(value.units_per_pixel > 0)) {
                        throw new Error('malformed projection view');
                    }
                    state.projection.center = {x: value.center_x, y: value.center_y};
                    state.projection.unitsPerPixel = value.units_per_pixel;
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
        exportPng: exportProjectionPNG,
        exportSvg: exportProjectionSVG,
        focusSelection: projectionFocusSelection,
        resetView: projectionResetView,
        updateSelectionControls: selectedCount => {
            document.getElementById('focus-selection').disabled = selectedCount === 0;
        },
    };
