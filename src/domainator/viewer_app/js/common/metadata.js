    function nodeId(nodeIndex) {
        return viewerNodeIds()[nodeIndex];
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

        state.metadataSearchTextByNodeIndex = viewerNodeIds().map((nodeName, nodeIndex) => {
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

    function populateMetadataControls() {
        renderAppearanceControls();
    }

    function currentColorField() {
        return document.getElementById('color-by')?.value || '';
    }

    function currentLabelField() {
        return document.getElementById('label-by')?.value || '';
    }

    function metadataColumn(name) {
        return state.metadataColumnByName.get(name) || null;
    }

    function metadataColumnKeys() {
        return ['node_id', ...state.metadataColumns.map(column => column.name)];
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
