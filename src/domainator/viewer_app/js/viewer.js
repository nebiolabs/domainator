    // Values the page is built with, from the viewer-config block (see viewer_html.py).
    const VIEWER_CONFIG = JSON.parse(document.getElementById('viewer-config').textContent);
    const state = {
        // The key in VIEWS of the view the open bundle is shown in.
        activeViewName: 'ssn',
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
    const EMBEDDED_BUNDLE_BASE64 = VIEWER_CONFIG.embeddedBundleBase64;
    const LAYOUT_WORKER_CODE = VIEWER_CONFIG.layoutWorkerCode;
    const DOMAINATOR_VERSION = VIEWER_CONFIG.domainatorVersion;
    const VIEWER_APP_NAME = VIEWER_CONFIG.viewerAppName;
    // The bundle format (.dvb) build_ssn, build_tree and build_projection write, and "Save
    // session" too (viewer_bundle.py).
    const VIEWER_BUNDLE_FORMAT = VIEWER_CONFIG.viewerBundleFormat;
    const VIEWER_BUNDLE_VERSION = VIEWER_CONFIG.viewerBundleVersion;
    const SUPPORTED_VIEWER_BUNDLE_VERSIONS = VIEWER_CONFIG.supportedViewerBundleVersions;
    // The network bundles (.dsnv) of earlier versions, which open as the .dvb bundle they
    // correspond to (ssn_bundle.py).
    const SSN_BUNDLE_FORMAT = VIEWER_CONFIG.ssnBundleFormat;
    const SUPPORTED_SSN_BUNDLE_VERSIONS = VIEWER_CONFIG.supportedSsnBundleVersions;

    // The views a bundle can open in, by name. Each js/views/<name>/view.js adds its own: which
    // bundles it opens and the hooks below. The shared code (js/common/) reaches the
    // open view only through these functions.
    const VIEWS = {};

    function activeView() {
        return VIEWS[state.activeViewName];
    }

    function viewerCanvas() {
        return activeView().canvas;
    }

    function viewerNodeIds() {
        return activeView().nodeIds();
    }

    function viewerNodeCount() {
        return viewerNodeIds().length;
    }

    // Redraws the view now, so that what follows sees the new pixels.
    function repaintNow() {
        activeView().repaintNow();
    }

    // Redraws the view on the next animation frame, once however often it is asked.
    function scheduleRepaint() {
        activeView().scheduleRepaint();
    }

    // [{name, valueFor(nodeIndex)}]: columns the view adds to the table export.
    function viewExportColumns() {
        return activeView().exportColumns();
    }

    // The view's clusters, for "Fill with cluster numbers": {defaultName, clusterCount,
    // numberFor(nodeIndex), description()}, or null when the view has none.
    function viewClusterAssignments() {
        return activeView().clusterAssignments();
    }

    // The .dvb bundle a .dsnv one corresponds to, as ssn_bundle.upgrade_ssn_viewer_bundle
    // makes it. Its graph's other keys are derived, and recomputed instead.
    function upgradeSsnBundle(bundle) {
        const upgraded = {
            format: VIEWER_BUNDLE_FORMAT,
            version: VIEWER_BUNDLE_VERSION,
            name: bundle.name,
            domainator_version: bundle.domainator_version,
            nodes: bundle.graph.nodes,
            metadata: bundle.metadata,
            defaults: bundle.defaults || {},
            views: {ssn: {
                mst_edges: bundle.graph.mst_edges,
                merge_impact_metric: bundle.graph.merge_impact_metric,
                hierarchy: bundle.graph.hierarchy,
            }},
        };
        if (bundle.app_state) {
            upgraded.app_state = bundle.app_state;
        }
        return upgraded;
    }

    function installBundle(bundle) {
        // A wrong `format` is fatal, but an unrecognized `version` is not: bundle
        // revisions so far have only added sections, so loading and reporting beats
        // refusing a file this viewer can very likely still render.
        const supportedVersions = bundle.format === SSN_BUNDLE_FORMAT ? SUPPORTED_SSN_BUNDLE_VERSIONS : SUPPORTED_VIEWER_BUNDLE_VERSIONS;
        const versionWarning = supportedVersions.includes(bundle.version)
            ? ''
            : ' (bundle version ' + bundle.version + ' is not one this viewer knows about ('
                + supportedVersions.join(', ') + '); loading anyway, some features may misbehave)';
        if (bundle.format === SSN_BUNDLE_FORMAT) {
            bundle = upgradeSsnBundle(bundle);
        }
        if (bundle.format !== VIEWER_BUNDLE_FORMAT) {
            throw new Error('Unsupported bundle format: ' + bundle.format);
        }
        const viewName = Object.keys(VIEWS).find(name => VIEWS[name].opens(bundle));
        if (viewName === undefined) {
            throw new Error('No view for this bundle\'s ' + Object.keys(bundle.views).join(', ') + ' view');
        }
        const view = VIEWS[viewName];
        state.bundleVersionWarning = versionWarning;
        state.activeViewName = viewName;
        document.body.dataset.activeView = viewName;
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
        view.reset(bundle);
        state.allNodeIndices = viewerNodeIds().map((_, i) => i);
        rebuildMetadataCaches();

        // Must precede rebuildNodeColorCache: the cache is keyed off the "Color by"
        // selection, which populateMetadataControls sets from bundle.defaults.
        populateMetadataControls();
        rebuildNodeColorCache();
        view.enableControls();
        document.getElementById('export-extent').disabled = false;
        document.getElementById('export-extent').value = 'viewport';
        // "Fill with cluster numbers" is for views that have clusters.
        document.getElementById('metadata-add-cluster-column').style.display = view.hasClusters ? '' : 'none';
        document.getElementById('customize-colors').disabled = false;
        document.getElementById('save-session').disabled = false;
        document.getElementById('rename-network').disabled = false;
        state.presetPreviewSlot = null;
        if (colorPickerIsOpen()) {
            openColorPicker();
        }
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
            'Loaded ' + (bundle.name || 'bundle') + ' with ' + viewerNodeCount().toLocaleString() + ' nodes.'
            + state.bundleVersionWarning + sessionNote
        );
        renderPresetSlots();
        view.afterInstalled();
        updateMetadataTable();
    }

// @include common/bundle_io.js
// @include common/controls.js
// @include common/metadata.js
// @include common/colors.js
// @include common/color_picker.js
// @include common/metadata_table.js
// @include common/export.js
// @include common/canvas_view.js
// @include views/ssn/slider.js
// @include views/ssn/layouts.js
// @include views/ssn/split_chart.js
// @include views/ssn/cluster_view.js
// @include views/ssn/controls.js
// @include views/ssn/view.js
// @include views/ssn/layout_core.js
// @include common/session_state.js
// @include common/selection_presets.js
// @include common/table_editing.js
// @include views/ssn/extraction.js
// @include common/gradient_stops.js
// @include common/column_charts.js
// @include views/ssn/split_chart_hover.js
// @include common/svg_context.js
// @include views/tree/tree.js
// @include views/tree/rooting.js
// @include views/tree/msa.js
// @include views/tree/render.js
// @include views/tree/controls.js
// @include views/tree/view.js
// @include views/projection/projection.js
// @include views/projection/view.js

    setupSharedViewerControls();
    setupSsnViewControls();
    setupTreeViewControls();
    setupProjectionViewControls();
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
