    // The listeners of the controls every view has: opening a file, coloring and labeling,
    // the color picker, the metadata table and the selection.
    function setupSharedViewerControls() {
        // The view toolbar's buttons act on whichever view is open.
        document.getElementById('export-png').addEventListener('click', () => activeView().exportPng());
        document.getElementById('export-svg').addEventListener('click', () => activeView().exportSvg());
        document.getElementById('focus-selection').addEventListener('click', () => activeView().focusSelection());
        document.getElementById('reset-view').addEventListener('click', () => activeView().resetView());
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
        document.getElementById('color-by').addEventListener('change', () => {
            rebuildNodeColorCache();
            repaintNow();
            if (colorPickerIsOpen()) {
                openColorPicker();
            }
        });
        document.getElementById('label-by').addEventListener('change', () => repaintNow());
        document.getElementById('export-selected').addEventListener('click', exportSelection);
        document.getElementById('customize-colors').addEventListener('click', openColorPicker);
        document.getElementById('save-session').addEventListener('click', saveSessionFile);
        document.getElementById('rename-network').addEventListener('click', openRenameNetworkDialog);
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
        document.getElementById('zoom-level').addEventListener('change', changeZoomLevel);
        document.getElementById('clear-selection').addEventListener('click', () => {
            state.selectedNodeIndices = new Set();
            resetMetadataPage();
            clearMetadataRowSelection();
            repaintNow();
            updateMetadataTable();
        });
    }
