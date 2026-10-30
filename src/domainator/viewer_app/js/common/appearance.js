    function resetAppearance() {
        state.appearanceNextId = 2;
        state.nodeStyles = [{id: 'node-1', color_by: state.bundle.defaults.color_by || '', shape_by: '', show_node: true}];
        state.labelStyles = [{id: 'label-1', text_by: state.bundle.defaults.label_by ||
            (state.activeViewName === 'tree' ? '__node_id__' : '')}];
        state.customShapes = {};
        state.shapeCaches = new Map();
        state.styleColorCaches = new Map();
        state.colorPickerColumn = null;
        state.colorPickerStyleId = null;
        document.getElementById('tree-show-labels').checked = true;
        document.getElementById('tree-show-color-bars').checked = true;
        document.getElementById('tree-color-branches').checked = true;
        document.getElementById('render-nodes').checked = true;
        closeShapePicker();
    }

    function appearanceText(style, nodeIndex) {
        if (!style?.text_by) { return ''; }
        const value = style.text_by === '__node_id__' ? nodeId(nodeIndex) : metadataValue(nodeIndex, style.text_by);
        return value === null || value === undefined ? '' : String(value);
    }

    function appearanceColor(style, nodeIndex, fallback) {
        if (!style?.color_by) { return fallback; }
        if (!state.styleColorCaches.has(style.color_by)) {
            state.styleColorCaches.set(style.color_by, state.allNodeIndices.map(index => nodeColor(index, style.color_by)));
        }
        return state.styleColorCaches.get(style.color_by)[nodeIndex];
    }

    function colorPickerField() {
        return state.colorPickerColumn === null ? currentColorField() : state.colorPickerColumn;
    }

    function customizeAppearanceColors(columnName, styleId) {
        state.colorPickerColumn = columnName;
        state.colorPickerStyleId = styleId;
        openColorPicker();
    }

    function appearancePropertyMenu(style, key, title, id) {
        const container = document.createElement('div');
        container.className = 'control';
        const label = document.createElement('label');
        label.htmlFor = id;
        label.textContent = title;
        const select = document.createElement('select');
        select.id = id;
        select.dataset.appearanceKey = key;
        select.add(new Option('None', ''));
        if (key === 'text_by') { select.add(new Option('node_id', '__node_id__')); }
        state.metadataColumns.forEach(column => select.add(new Option(column.name, column.name)));
        select.value = style[key];
        container.append(label, select);
        return container;
    }

    function appearanceButton(text, action, disabled = false) {
        const button = document.createElement('button');
        button.type = 'button';
        button.textContent = text;
        button.dataset.appearanceAction = action;
        button.disabled = disabled;
        return button;
    }

    function renderAppearanceControls() {
        if (!state.bundle) { return; }
        const focusedId = document.activeElement?.id;
        [['node', state.nodeStyles], ['label', state.labelStyles]].forEach(([kind, styles]) => {
            const list = document.getElementById(kind + '-styles');
            list.replaceChildren();
            styles.forEach((style, index) => {
                const entry = document.createElement('div');
                entry.className = 'appearance-entry';
                entry.dataset.styleId = style.id;
                entry.dataset.styleKind = kind;
                const heading = document.createElement('strong');
                heading.textContent = (kind === 'node' ? 'Node ' : 'Label ') + (index + 1);
                entry.append(heading);
                if (kind === 'node') {
                    if (state.activeViewName === 'tree') {
                        const container = document.createElement('div');
                        container.className = 'control checkbox';
                        const checkbox = document.createElement('input');
                        checkbox.type = 'checkbox';
                        checkbox.id = style.id + '-show-node';
                        checkbox.dataset.appearanceKey = 'show_node';
                        checkbox.checked = style.show_node !== false;
                        const label = document.createElement('label');
                        label.htmlFor = checkbox.id;
                        label.textContent = 'Show node';
                        container.append(checkbox, label);
                        entry.append(container);
                    }
                    entry.append(appearancePropertyMenu(style, 'color_by', 'Color by', index === 0 ? 'color-by' : style.id + '-color'));
                    const colors = appearanceButton('Customize colors...', 'colors', !style.color_by);
                    if (index === 0) { colors.id = 'customize-colors'; }
                    entry.append(colors);
                    entry.append(appearancePropertyMenu(style, 'shape_by', 'Shape by', index === 0 ? 'shape-by' : style.id + '-shape'));
                    entry.append(appearanceButton('Customize shapes...', 'shapes', !style.shape_by));
                } else {
                    entry.append(appearancePropertyMenu(style, 'text_by', 'Text by', index === 0 ? 'label-by' : style.id + '-text'));
                }
                if (state.activeViewName === 'tree') {
                    const actions = document.createElement('div');
                    actions.className = 'appearance-actions';
                    [['Move up', 'up', index === 0], ['Move down', 'down', index === styles.length - 1], ['Remove', 'remove', false]]
                        .forEach(([text, action, disabled]) => {
                            const button = appearanceButton(text, action, disabled);
                            button.setAttribute('aria-label', text + ' ' + kind + ' ' + (index + 1));
                            actions.append(button);
                        });
                        if (focusedId) { document.getElementById(focusedId)?.focus({preventScroll: true}); }
                    entry.append(actions);
                }
                list.append(entry);
            });
            document.getElementById('add-' + kind + '-style').hidden = state.activeViewName !== 'tree';
        });
    }

    function appearanceChanged() {
        renderAppearanceControls();
        rebuildNodeColorCache();
        repaintNow();
        if (colorPickerIsOpen()) { openColorPicker(); }
    }

    function setupAppearanceControls() {
        ['node', 'label'].forEach(kind => {
            const list = document.getElementById(kind + '-styles');
            list.addEventListener('change', event => {
                const key = event.target.dataset.appearanceKey;
                if (!key) { return; }
                const styles = kind === 'node' ? state.nodeStyles : state.labelStyles;
                const style = styles.find(item => item.id === event.target.closest('[data-style-id]').dataset.styleId);
                style[key] = event.target.type === 'checkbox' ? event.target.checked : event.target.value;
                if (key === 'color_by' && state.colorPickerStyleId === style.id) {
                    state.colorPickerColumn = style.color_by;
                }
                if (key === 'shape_by' && state.shapePickerStyleId === style.id) {
                    if (style.shape_by) { state.shapePickerColumn = style.shape_by; renderShapePicker(); }
                    else { closeShapePicker(); }
                }
                appearanceChanged();
            });
            list.addEventListener('click', event => {
                const button = event.target.closest('[data-appearance-action]');
                if (!button) { return; }
                const styles = kind === 'node' ? state.nodeStyles : state.labelStyles;
                const index = styles.findIndex(item => item.id === button.closest('[data-style-id]').dataset.styleId);
                const style = styles[index];
                const action = button.dataset.appearanceAction;
                if (action === 'colors') { customizeAppearanceColors(style.color_by, style.id); return; }
                if (action === 'shapes') { openShapePicker(style.shape_by, style.id); return; }
                if (action === 'remove') {
                    if (state.colorPickerStyleId === style.id) { closeColorPicker(); }
                    if (state.shapePickerStyleId === style.id) { closeShapePicker(); }
                    styles.splice(index, 1);
                }
                else {
                    const other = index + (action === 'up' ? -1 : 1);
                    [styles[index], styles[other]] = [styles[other], styles[index]];
                }
                appearanceChanged();
            });
            document.getElementById('add-' + kind + '-style').addEventListener('click', () => {
                const styles = kind === 'node' ? state.nodeStyles : state.labelStyles;
                const id = kind + '-' + state.appearanceNextId++;
                styles.push(kind === 'node' ? {id, color_by: styles[0]?.color_by || '', shape_by: '', show_node: true}
                    : {id, text_by: '__node_id__'});
                appearanceChanged();
            });
        });
        document.getElementById('shape-picker-close').addEventListener('click', closeShapePicker);
        document.getElementById('shape-picker-overlay').addEventListener('click', event => {
            if (event.target === event.currentTarget) { closeShapePicker(); }
        });
        document.addEventListener('keydown', event => {
            if (event.key === 'Escape' && shapePickerIsOpen()) { closeShapePicker(); }
        });
        ['prev', 'next'].forEach(direction => {
            document.getElementById('shape-picker-' + direction).addEventListener('click', () => {
                state.shapePickerPage += direction === 'prev' ? -1 : 1;
                renderShapePicker();
            });
        });
        document.getElementById('shape-picker-reset').addEventListener('click', () => {
            delete state.customShapes[state.shapePickerColumn];
            rebuildNodeColorCache();
            repaintNow();
            renderShapePicker();
        });
    }

    function validateAppearanceStyles(value, kind) {
        if (!Array.isArray(value) || (state.activeViewName !== 'tree' && value.length !== 1)) {
            throw new Error('invalid ' + kind + ' styles');
        }
        const ids = new Set();
        const keys = kind === 'node' ? ['color_by', 'shape_by'] : ['text_by'];
        const styles = value.map(entry => {
            if (!entry || typeof entry.id !== 'string' || !entry.id.startsWith(kind + '-') ||
                    !/^(node|label)-[0-9]+$/.test(entry.id) ||
                    !Number.isSafeInteger(Number(entry.id.split('-')[1])) || ids.has(entry.id)) {
                throw new Error('invalid or duplicate style id');
            }
            ids.add(entry.id);
            const style = {id: entry.id};
            keys.forEach(key => {
                const field = entry[key];
                if (typeof field !== 'string' || (field !== '' && !(key === 'text_by' && field === '__node_id__') && !metadataColumn(field))) {
                    throw new Error('unknown property "' + field + '"');
                }
                style[key] = field;
            });
            if (kind === 'node') {
                if (entry.show_node !== undefined && typeof entry.show_node !== 'boolean') {
                    throw new Error('invalid node visibility');
                }
                style.show_node = entry.show_node ?? true;
            }
            return style;
        });
        return styles;
    }

    function restoreAppearanceStyles(value, kind) {
        const styles = validateAppearanceStyles(value, kind);
        if (kind === 'node') { state.nodeStyles = styles; }
        else { state.labelStyles = styles; }
        const ids = [...state.nodeStyles, ...state.labelStyles].map(style => Number(style.id.split('-')[1]));
        state.appearanceNextId = ids.reduce((max, id) => Math.max(max, id), 1) + 1;
        renderAppearanceControls();
    }

    function retargetAppearanceProperty(oldName, newName) {
        let changed = 0;
        [...state.nodeStyles, ...state.labelStyles].forEach(style => {
            ['color_by', 'shape_by', 'text_by'].forEach(key => {
                if (style[key] === oldName) { style[key] = newName; changed++; }
            });
        });
        if (Object.hasOwn(state.customShapes, oldName)) {
            if (newName) {
                Object.defineProperty(state.customShapes, newName, {
                    value: state.customShapes[oldName], writable: true, enumerable: true, configurable: true,
                });
            }
            delete state.customShapes[oldName];
        }
        if (state.colorPickerColumn === oldName) {
            if (newName) { state.colorPickerColumn = newName; } else { closeColorPicker(); state.colorPickerColumn = null; }
        }
        if (state.shapePickerColumn === oldName) {
            if (newName) { state.shapePickerColumn = newName; } else { closeShapePicker(); }
        }
        return changed;
    }
