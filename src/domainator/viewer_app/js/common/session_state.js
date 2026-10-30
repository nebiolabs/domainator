/*
 * Session save/load: the `app_state` section of a v4 bundle.
 *
 * Included (`// @include`) into viewer.js, the viewer's one `<script>` scope, so
 * it can call `state`, `repaintNow()` and friends directly, and its
 * function declarations hoist above the listener wiring.
 */

    // ---------------------------------------------------------------------
    // Session state
    //
    // A saved session is an ordinary bundle whose `metadata` carries whatever
    // the user edited plus one extra top-level section, `app_state`. Because
    // the viewer's UI is expected to churn, state is *not* dumped ad hoc:
    // every persisted control is one entry in a view's sessionFields.
    // Serializing walks the registry; applying looks each saved key up in it.
    // That makes both directions tolerant -- a key this build dropped is
    // skipped and reported, and a key the file omits keeps its default -- so
    // old files open in new viewers and new files open in old ones.
    // ---------------------------------------------------------------------
    const SESSION_STATE_VERSION = 1;
    const PRESET_SLOT_COUNT = 10;

    function domValueField(section, key, elementId, options) {
        const parse = (options && options.parse) || (value => value);
        const format = (options && options.format) || (value => String(value));
        return {
            section,
            key,
            get: () => parse(document.getElementById(elementId).value),
            set: value => {
                const element = document.getElementById(elementId);
                const text = format(value);
                if (element.tagName === 'SELECT' &&
                        !Array.from(element.options).some(option => option.value === text)) {
                    // Restoring a column or option this bundle does not have would
                    // silently blank the control; report it as skipped instead.
                    throw new Error('no option "' + text + '"');
                }
                element.value = text;
            },
        };
    }

    function domCheckboxField(section, key, elementId) {
        return {
            section,
            key,
            get: () => document.getElementById(elementId).checked,
            set: value => { document.getElementById(elementId).checked = Boolean(value); },
        };
    }

    function nodeIndexByIdMap() {
        if (!state.nodeIndexById) {
            state.nodeIndexById = new Map(viewerNodeIds().map((id, index) => [id, index]));
        }
        return state.nodeIndexById;
    }

    // Selections persist as node ids, never indices: an index is a position in
    // the bundle's nodes and silently means a different sequence if the bundle is
    // rebuilt or subsetted. Unmatched ids are counted so the load can say so.
    function nodeIndicesFromIds(ids) {
        const map = nodeIndexByIdMap();
        const indices = [];
        let missing = 0;
        (ids || []).forEach(id => {
            const index = map.get(id);
            if (index === undefined) { missing += 1; } else { indices.push(index); }
        });
        state.sessionMissingNodeIds += missing;
        return indices;
    }

    function nodeIdsFromIndices(indices) {
        return Array.from(indices).sort((left, right) => left - right).map(index => nodeId(index));
    }

    function serializeSelectionPresets() {
        const presets = {};
        state.selectionPresets.forEach((preset, slot) => {
            presets[String(slot)] = {
                node_ids: nodeIdsFromIndices(preset.nodeIndices),
                saved_at: preset.savedAt || null,
            };
        });
        return presets;
    }

    function deserializeSelectionPresets(presets) {
        state.selectionPresets = new Map();
        Object.keys(presets || {}).forEach(slotKey => {
            const slot = Number.parseInt(slotKey, 10);
            if (!Number.isInteger(slot) || slot < 0 || slot >= PRESET_SLOT_COUNT) {
                throw new Error('bad preset slot "' + slotKey + '"');
            }
            const entry = presets[slotKey] || {};
            const indices = nodeIndicesFromIds(entry.node_ids);
            if (indices.length === 0) { return; }
            state.selectionPresets.set(slot, {
                nodeIndices: new Set(indices),
                savedAt: entry.saved_at || null,
            });
        });
    }

    // The fields every view saves, by key. A view's sessionFields lists them, in its own order,
    // among the fields only it has; the order is the order a session is applied in.
    function sharedSessionFields() {
        return {
            colorBy: {
                section: 'view', key: 'color_by', get: currentColorField,
                set: value => {
                    if (value !== '' && !metadataColumn(value)) { throw new Error('no option "' + value + '"'); }
                    state.nodeStyles[0].color_by = value;
                    renderAppearanceControls();
                },
            },
            labelBy: {
                section: 'view', key: 'label_by', get: currentLabelField,
                set: value => {
                    if (value !== '' && value !== '__node_id__' && !metadataColumn(value)) { throw new Error('no option "' + value + '"'); }
                    state.labelStyles[0].text_by = value || (state.activeViewName === 'tree' ? '__node_id__' : '');
                    renderAppearanceControls();
                },
            },
            exportPngScale: domValueField('view', 'export_png_scale', 'export-png-scale'),
            exportExtent: domValueField('view', 'export_extent', 'export-extent'),
            viewTransform: {
                section: 'view',
                key: 'view_transform',
                get: () => ({
                    scale: state.viewTransform.scale,
                    offsetX: state.viewTransform.offsetX,
                    offsetY: state.viewTransform.offsetY,
                }),
                set: value => {
                    if (!value || !Number.isFinite(value.scale) || value.scale <= 0 ||
                            !Number.isFinite(value.offsetX) || !Number.isFinite(value.offsetY)) {
                        throw new Error('malformed view transform');
                    }
                    // Consumed by applyComputedLayout once a layout exists, which is
                    // also what would otherwise auto-fit the view and discard this.
                    state.pendingRestoreViewTransform = {
                        scale: value.scale, offsetX: value.offsetX, offsetY: value.offsetY,
                    };
                },
            },
            filter: domValueField('table', 'filter', 'metadata-filter'),
            nullOrder: domValueField('table', 'null_order', 'metadata-null-order'),
            rowsPerPage: domValueField('table', 'rows_per_page', 'metadata-rows-per-page'),
            sort: {
                section: 'table',
                key: 'sort',
                get: () => ({
                    column_key: state.metadataSort.columnKey,
                    direction: state.metadataSort.direction,
                }),
                set: value => {
                    const columnKey = value && value.column_key;
                    if (columnKey && !metadataColumnKeys().includes(columnKey)) {
                        throw new Error('no column "' + columnKey + '"');
                    }
                    state.metadataSort = {
                        columnKey: columnKey || null,
                        direction: (value && value.direction === 'desc') ? 'desc' : 'asc',
                    };
                },
            },
            columnWidths: {
                section: 'table',
                key: 'column_widths',
                get: () => Object.fromEntries(state.metadataColumnWidths),
                set: value => {
                    const widths = new Map();
                    Object.keys(value || {}).forEach(columnKey => {
                        const width = Number(value[columnKey]);
                        if (!Number.isFinite(width)) { return; }
                        widths.set(columnKey, Math.max(96, Math.min(640, Math.round(width))));
                    });
                    state.metadataColumnWidths = widths;
                },
            },
            customPalettes: {
                section: 'colors',
                key: 'custom_palettes',
                get: () => state.customPalettes,
                set: value => {
                    if (value === null || typeof value !== 'object') { throw new Error('not an object'); }
                    // Sessions written before gradients became stop lists carry
                    // {low, mid, high}; rewrite them so nothing downstream has to
                    // know about that shape.
                    state.customPalettes = normalizeCustomPalettes(value);
                },
            },
            categoricalColumns: {
                section: 'colors',
                key: 'categorical_columns',
                get: () => Array.from(state.categoricalColumns),
                set: value => {
                    state.categoricalColumns = new Set(value || []);
                    // Whether a numeric column colors as categories is *derived* from this
                    // set, and rebuildMetadataCaches -- which derives it -- has already run
                    // by the time a session is applied. Without recomputing it here the
                    // checkbox would read categorical while the colors and the picker
                    // stayed in gradient mode until the box was toggled by hand.
                    state.metadataColorInfoByName.forEach((info, columnName) => {
                        if (info.baseType === 'numeric') {
                            info.type = state.categoricalColumns.has(columnName) ? 'categorical' : 'numeric';
                        }
                    });
                },
            },
            nodeIds: {
                section: 'selection',
                key: 'node_ids',
                get: () => nodeIdsFromIndices(state.selectedNodeIndices),
                set: value => { state.selectedNodeIndices = new Set(nodeIndicesFromIds(value)); },
            },
            presets: {
                section: 'selection',
                key: 'presets',
                get: () => serializeSelectionPresets(),
                set: value => { deserializeSelectionPresets(value); },
            },
        };
    }

    // The open view's session fields, built once per view.
    function viewStateFields() {
        const view = activeView();
        if (!view.stateFields) {
            view.stateFields = view.sessionFields(sharedSessionFields()).concat([
                {section: 'view', key: 'node_styles',
                    get: () => state.nodeStyles.map((style, index) => ({...style, color_by: index === 0 ? currentColorField() : style.color_by})),
                    set: value => restoreAppearanceStyles(value, 'node')},
                {section: 'view', key: 'label_styles',
                    get: () => state.labelStyles.map((style, index) => ({...style, text_by: index === 0 ? currentLabelField() : style.text_by})),
                    set: value => restoreAppearanceStyles(value, 'label')},
                {section: 'shapes', key: 'custom_mappings', get: () => state.customShapes,
                    set: value => { state.customShapes = validateShapeMappings(value); }},
                ...['node', 'label'].map(kind => ({
                    section: 'view', key: kind + '_settings_open',
                    get: () => document.getElementById(kind + '-settings').open,
                    set: value => {
                        if (typeof value !== 'boolean') { throw new Error('invalid section expansion'); }
                        document.getElementById(kind + '-settings').open = value;
                    },
                })),
            ]);
        }
        return view.stateFields;
    }

    function collectSessionState() {
        const sections = {};
        viewStateFields().forEach(field => {
            if (!sections[field.section]) { sections[field.section] = {}; }
            sections[field.section][field.key] = field.get();
        });
        return Object.assign({
            state_version: SESSION_STATE_VERSION,
            saved_by: DOMAINATOR_VERSION,
            saved_at: new Date().toISOString(),
        }, sections);
    }

    // Returns a human-readable summary of anything that could not be restored.
    function applySessionState(appState) {
        if (!appState || typeof appState !== 'object') { return ''; }
        state.sessionMissingNodeIds = 0;
        const fieldsByPath = new Map(
            viewStateFields().map(field => [field.section + '.' + field.key, field])
        );
        const skipped = [];
        Object.keys(appState).forEach(sectionName => {
            const section = appState[sectionName];
            // Scalars at the top level are metadata about the save itself
            // (state_version, saved_by, saved_at), not restorable fields.
            if (!section || typeof section !== 'object' || Array.isArray(section)) { return; }
            // Legacy aliases apply first, independent of the file's key order.
            const keys = Object.keys(section).sort((a, b) =>
                Number(['node_styles', 'label_styles'].includes(a)) - Number(['node_styles', 'label_styles'].includes(b)));
            keys.forEach(key => {
                const field = fieldsByPath.get(sectionName + '.' + key);
                if (!field) {
                    skipped.push(sectionName + '.' + key + ' (unknown)');
                    return;
                }
                try {
                    field.set(section[key]);
                } catch (error) {
                    skipped.push(sectionName + '.' + key + ' (' + error.message + ')');
                }
            });
        });

        activeView().afterSessionApplied();

        const notes = [];
        if (skipped.length > 0) {
            notes.push('skipped ' + skipped.length + ' saved setting' + (skipped.length === 1 ? '' : 's') +
                ': ' + skipped.join(', '));
        }
        if (state.sessionMissingNodeIds > 0) {
            notes.push(state.sessionMissingNodeIds.toLocaleString() +
                ' saved node id' + (state.sessionMissingNodeIds === 1 ? '' : 's') + ' are not in this bundle');
        }
        return notes.length > 0 ? ' (' + notes.join('; ') + ')' : '';
    }

    // One small dialog serves both "Rename network" and naming an extraction: the
    // network's name is a single field on the bundle, and it drives the heading, the
    // browser tab and the name of every file saved from here.
    function openNameDialog(options) {
        state.pendingNameAction = options.onConfirm;
        document.getElementById('name-dialog-title').textContent = options.title;
        document.getElementById('name-dialog-label').textContent = options.label;
        document.getElementById('name-dialog-note').textContent = options.note || '';
        document.getElementById('name-apply').textContent = options.confirmLabel || 'Apply';
        const input = document.getElementById('name-value');
        input.value = options.value || '';
        document.getElementById('name-overlay').hidden = false;
        input.focus();
        input.select();
    }

    function closeNameDialog() {
        document.getElementById('name-overlay').hidden = true;
        state.pendingNameAction = null;
    }

    function nameDialogIsOpen() {
        return !document.getElementById('name-overlay').hidden;
    }

    function confirmNameDialog() {
        const value = document.getElementById('name-value').value.trim();
        if (value === '') {
            setStatus('Enter a name.');
            return;
        }
        const action = state.pendingNameAction;
        closeNameDialog();
        if (action) {
            action(value);
        }
    }

    function renameNetwork(name) {
        if (!state.bundle) { return; }
        const previous = state.bundle.name || '';
        state.bundle.name = name;
        updateViewerTitle(name);
        setStatus('Renamed "' + previous + '" to "' + name + '". Files saved from here use the new name.');
    }

    function openRenameNetworkDialog() {
        if (!state.bundle) { return; }
        openNameDialog({
            title: 'Rename ' + activeView().bundleNoun,
            label: activeView().bundleNoun.charAt(0).toUpperCase() + activeView().bundleNoun.slice(1) + ' name',
            note: 'Shown in the heading and the browser tab, stored in the bundle, and used to name saved files.',
            value: state.bundle.name || '',
            confirmLabel: 'Rename',
            onConfirm: renameNetwork,
        });
    }

    function setupNameDialog() {
        document.getElementById('name-cancel').addEventListener('click', closeNameDialog);
        document.getElementById('name-apply').addEventListener('click', confirmNameDialog);
        document.getElementById('name-overlay').addEventListener('click', event => {
            if (event.target === event.currentTarget) { closeNameDialog(); }
        });
        document.getElementById('name-value').addEventListener('keydown', event => {
            if (event.key === 'Enter') {
                event.preventDefault();
                confirmNameDialog();
            }
        });
        document.addEventListener('keydown', event => {
            if (event.key === 'Escape' && nameDialogIsOpen()) { closeNameDialog(); }
        });
        // The heading shows the name, so it is the obvious thing to double-click.
        document.getElementById('viewer-title').addEventListener('dblclick', openRenameNetworkDialog);
    }

    function browserSupportsBundleSaving() {
        return 'CompressionStream' in window;
    }

    // Serialize a bundle and hand it to the browser as a download. Shared by
    // "Save session" and "Save extraction" so both produce the same file format.
    async function downloadBundleFile(payload, filenameBase, description) {
        try {
            setStatus(description + '...');
            const bytes = new TextEncoder().encode(JSON.stringify(payload));
            let blob;
            let extension;
            if (browserSupportsBundleSaving()) {
                const stream = new Blob([bytes]).stream().pipeThrough(new CompressionStream('gzip'));
                blob = await new Response(stream).blob();
                extension = '.dvb';
            } else {
                // decodeBundleFile falls back to plain JSON, so this still round-trips.
                blob = new Blob([bytes], {type: 'application/json'});
                extension = '.json';
            }
            triggerDownload(URL.createObjectURL(blob), filenameBase + extension, true);
            setStatus(description + ' (' + Math.round(blob.size / 1024).toLocaleString() + ' KB).');
            return true;
        } catch (error) {
            console.error(error);
            setStatus(description + ' failed: ' + error.message);
            return false;
        }
    }

    async function saveSessionFile() {
        if (!state.bundle) { return; }
        // state.metadataByNodeIndex is bundle.metadata.rows, so table edits are
        // already in here; the file stays a valid bundle that Python readers
        // and this viewer can both open.
        await downloadBundleFile(
            Object.assign({}, state.bundle, {version: VIEWER_BUNDLE_VERSION, app_state: collectSessionState()}),
            exportBaseName() + '_session',
            'Saved session'
        );
    }
