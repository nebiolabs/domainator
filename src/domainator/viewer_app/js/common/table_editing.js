/*
 * Editable metadata cells, user-added columns, and bulk fills.
 *
 * Included into viewer.js's script scope; see session_state.js.
 */

    // ---------------------------------------------------------------------
    // Metadata editing
    //
    // Edits write straight into state.metadataByNodeIndex, which *is*
    // bundle.metadata.rows (installBundle assigns it by reference). A saved
    // session is therefore an ordinary bundle whose metadata simply reflects
    // the annotations -- no overlay to reconcile, and Python readers
    // (ssn_bundle, ssn_navigator) see the edits with no extra work.
    //
    // node_id is not editable: it is the graph key that node indices, the
    // hierarchy and every saved selection resolve through.
    // ---------------------------------------------------------------------

    function columnIsEditable(columnKey) {
        return Boolean(state.bundle) && columnKey !== 'node_id' && Boolean(metadataColumn(columnKey));
    }

    function columnIsViewerAdded(columnKey) {
        const column = metadataColumn(columnKey);
        return Boolean(column && column.origin === 'viewer');
    }

    // Parse user text into the column's declared type. Returns {ok, value, widenTo}
    // so callers can reject rather than silently coercing "abc" in a numeric
    // column to NaN (or quietly turning the whole column into text).
    //
    // Typing 12.5 into an int column widens the column to float rather than
    // truncating to 12: the type was only inferred from the values that happened
    // to be in the source TSV, so widening is what Python would have inferred had
    // the data looked like this to begin with, and it loses nothing.
    function parseMetadataInput(columnKey, rawText) {
        const text = String(rawText === null || rawText === undefined ? '' : rawText).trim();
        if (text === '' || text === '—') {
            return {ok: true, value: null};
        }
        const columnType = metadataColumnType(columnKey);
        if (columnType === 'int' || columnType === 'float') {
            // Reject "12abc": Number() is lenient about surrounding whitespace but
            // not about trailing garbage, which is the failure mode worth catching.
            const numeric = Number(text);
            if (!Number.isFinite(numeric)) {
                return {ok: false, reason: '"' + text + '" is not a number'};
            }
            if (columnType === 'int' && !Number.isInteger(numeric)) {
                return {ok: true, value: numeric, widenTo: 'float'};
            }
            return {ok: true, value: numeric};
        }
        return {ok: true, value: text};
    }

    // Applies a parse result's requested type widening, if any. Returns whether
    // the column list changed, so callers can refresh the column-driven menus.
    function applyColumnWidening(columnName, parsed) {
        if (!parsed.widenTo) { return false; }
        const column = metadataColumn(columnName);
        if (!column || column.type === parsed.widenTo) { return false; }
        column.type = parsed.widenTo;
        setStatus('Column "' + columnName + '" now holds fractional values; its type widened to float.');
        return true;
    }

    // Writes the value without touching caches; callers batch a single
    // refreshAfterMetadataEdit() so a bulk fill rebuilds caches once, not once
    // per row (rebuildMetadataCaches is O(rows x columns)).
    function writeMetadataValue(nodeIndex, columnName, value) {
        const columnIndex = state.metadataColumnIndexByName.get(columnName);
        if (columnIndex === undefined) { return false; }
        const row = state.metadataByNodeIndex[nodeIndex];
        if (!row) { return false; }
        row[columnIndex] = value;
        return true;
    }

    // Everything downstream of a metadata value: the search index and numeric
    // min/max live in rebuildMetadataCaches, node colors are cached separately.
    function refreshAfterMetadataEdit(options) {
        const settings = options || {};
        rebuildMetadataCaches();
        if (settings.columnsChanged) {
            populateMetadataControlsPreservingSelection();
            populateEditableColumnMenus();
        }
        rebuildNodeColorCache();
        updateMetadataTable();
        repaintNow();
        if (colorPickerIsOpen()) {
            openColorPicker();
        }
    }

    // populateMetadataControls resets color-by/label-by to the bundle defaults,
    // which would undo the user's choice every time a column is added.
    function populateMetadataControlsPreservingSelection() {
        const colorBy = document.getElementById('color-by');
        const labelBy = document.getElementById('label-by');
        const previousColor = colorBy.value;
        const previousLabel = labelBy.value;
        populateMetadataControls();
        if (Array.from(colorBy.options).some(option => option.value === previousColor)) {
            colorBy.value = previousColor;
        }
        if (Array.from(labelBy.options).some(option => option.value === previousLabel)) {
            labelBy.value = previousLabel;
        }
    }

    function setMetadataValue(nodeIndex, columnName, rawText) {
        if (!columnIsEditable(columnName)) { return false; }
        const parsed = parseMetadataInput(columnName, rawText);
        if (!parsed.ok) {
            setStatus('Cannot set "' + columnName + '": ' + parsed.reason + '.');
            return false;
        }
        if (!writeMetadataValue(nodeIndex, columnName, parsed.value)) { return false; }
        const widened = applyColumnWidening(columnName, parsed);
        refreshAfterMetadataEdit({columnsChanged: widened});
        return true;
    }

    // ---- in-place cell editing ----

    function beginMetadataCellEdit(cell) {
        if (state.metadataEditCell) { return; }
        const row = cell.closest('tr[data-node-index]');
        const columnKey = cell.dataset.columnKey;
        if (!row || !columnIsEditable(columnKey)) { return; }
        const nodeIndex = Number(row.dataset.nodeIndex);
        const value = metadataValue(nodeIndex, columnKey);

        const input = document.createElement('input');
        input.type = 'text';
        input.className = 'metadata-cell-input';
        input.value = isMissingMetadataValue(value) ? '' : String(value);
        input.setAttribute('aria-label', 'Edit ' + columnKey);
        cell.textContent = '';
        cell.appendChild(input);
        state.metadataEditCell = {cell, nodeIndex, columnKey};
        input.focus();
        input.select();

        const commit = () => finishMetadataCellEdit(true);
        input.addEventListener('blur', commit);
        input.addEventListener('keydown', event => {
            if (event.key === 'Enter') {
                event.preventDefault();
                input.removeEventListener('blur', commit);
                finishMetadataCellEdit(true);
            } else if (event.key === 'Escape') {
                event.preventDefault();
                input.removeEventListener('blur', commit);
                finishMetadataCellEdit(false);
            }
        });
    }

    function finishMetadataCellEdit(commit) {
        const editing = state.metadataEditCell;
        if (!editing) { return; }
        state.metadataEditCell = null;
        const input = editing.cell.querySelector('input');
        const rawText = input ? input.value : '';
        // Clicking away from an untouched cell is the common case, and it should
        // not pay for a cache rebuild -- it only needs the <input> swapped out.
        const parsed = commit ? parseMetadataInput(editing.columnKey, rawText) : null;
        const unchanged = parsed && parsed.ok && !parsed.widenTo &&
            parsed.value === metadataValue(editing.nodeIndex, editing.columnKey);
        if (!commit || unchanged ||
                !setMetadataValue(editing.nodeIndex, editing.columnKey, rawText)) {
            // Canceled, unchanged, or rejected: repaint to restore the cell text.
            updateMetadataTable();
        }
    }

    // ---- adding and removing columns ----

    // Shared by add and rename. `node_id` is reserved: it is the table's synthetic
    // first column, so a metadata column of that name would be unreachable.
    function validateColumnName(name) {
        const columnName = String(name || '').trim();
        if (columnName === '') {
            return {ok: false, reason: 'Enter a column name.'};
        }
        if (columnName === 'node_id' || state.metadataColumnIndexByName.has(columnName)) {
            return {ok: false, reason: 'A column named "' + columnName + '" already exists.'};
        }
        return {ok: true, name: columnName};
    }

    function addMetadataColumn(name, columnType) {
        if (!state.bundle) { return false; }
        const validated = validateColumnName(name);
        if (!validated.ok) {
            setStatus(validated.reason);
            return false;
        }
        const columnName = validated.name;
        // `origin` marks columns created here so they can be deleted again. It is
        // additive: the build tools write only {name, type} and no reader
        // inspects extra keys, so a bundle carrying it stays valid.
        state.metadataColumns.push({
            name: columnName,
            type: columnType === 'number' ? 'float' : 'str',
            origin: 'viewer',
        });
        state.metadataByNodeIndex.forEach(row => { row.push(null); });
        refreshAfterMetadataEdit({columnsChanged: true});
        setStatus('Added column "' + columnName + '".');
        return true;
    }

    // Menus whose option values are column names. Which column each one points at
    // is held in the DOM rather than in `state`, so it is one more thing keyed by
    // column name that a rename has to carry across.
    const COLUMN_NAME_MENU_IDS = [
        'color-by',
        'label-by',
        'metadata-fill-column',
        'metadata-paste-column',
        'metadata-rename-column',
        'metadata-delete-column',
        'metadata-select-column',
    ];

    // Retarget every menu currently pointing at `oldName`. Rewriting the option
    // updates its select's value, so the rebuild that follows -- which restores
    // whatever each menu was showing -- carries the choice over on its own.
    function renameColumnInMenus(oldName, newName) {
        COLUMN_NAME_MENU_IDS.forEach(menuId => {
            Array.from(document.getElementById(menuId).options).forEach(option => {
                if (option.value === oldName) {
                    option.value = newName;
                }
            });
        });
    }

    // Writes the view's cluster numbers (see viewClusterAssignments) into a column. An
    // existing column of the same name is overwritten.
    function addClusterColumn(name) {
        if (!state.bundle) { return false; }
        const clusters = viewClusterAssignments();
        if (!clusters) { return false; }
        const columnName = String(name || '').trim() || clusters.defaultName;
        if (columnName === 'node_id') {
            setStatus('"node_id" is reserved; choose another column name.');
            return false;
        }
        const existing = metadataColumn(columnName);
        if (existing && existing.origin !== 'viewer' &&
                !window.confirm('Overwrite the column "' + columnName + '" with cluster numbers?\n\n' +
                    'It came from this bundle\'s metadata, not from the viewer, and its current values will be lost.')) {
            return false;
        }

        if (!existing) {
            state.metadataColumns.push({name: columnName, type: 'int', origin: 'viewer'});
            state.metadataByNodeIndex.forEach(row => { row.push(null); });
        } else {
            existing.type = 'int';
        }
        rebuildMetadataCaches();  // so writeMetadataValue can find the new column index
        state.allNodeIndices.forEach(nodeIndex => {
            writeMetadataValue(nodeIndex, columnName, clusters.numberFor(nodeIndex));
        });
        // Cluster ids are labels, not magnitudes, so color them as discrete
        // categories -- the same thing the --categorical option of build_ssn, build_tree and build_projection asks for.
        state.categoricalColumns.add(columnName);

        refreshAfterMetadataEdit({columnsChanged: !existing});
        setStatus('Wrote ' + clusters.clusterCount.toLocaleString() + ' cluster numbers into "' +
            columnName + '" ' + clusters.description() + '.');
        return true;
    }

    // Moves any state that is keyed by column name along with the column. Values
    // themselves are positional, so the rows need no rewriting.
    function renameMetadataColumn(oldName, newName) {
        if (!state.bundle || !columnIsEditable(oldName)) { return false; }
        const trimmed = String(newName || '').trim();
        if (trimmed === oldName) { return true; }
        const validated = validateColumnName(trimmed);
        if (!validated.ok) {
            setStatus(validated.reason);
            return false;
        }
        const columnName = validated.name;

        // Both palette spellings, because paletteKey() suffixes a numeric column
        // that is currently colored categorically.
        ['', '\u0000categorical'].forEach(suffix => {
            const fromKey = oldName + suffix;
            if (Object.prototype.hasOwnProperty.call(state.customPalettes, fromKey)) {
                state.customPalettes[columnName + suffix] = state.customPalettes[fromKey];
                delete state.customPalettes[fromKey];
            }
        });
        if (state.categoricalColumns.delete(oldName)) {
            state.categoricalColumns.add(columnName);
        }
        if (state.metadataColumnWidths.has(oldName)) {
            state.metadataColumnWidths.set(columnName, state.metadataColumnWidths.get(oldName));
            state.metadataColumnWidths.delete(oldName);
        }
        if (state.metadataSort.columnKey === oldName) {
            state.metadataSort.columnKey = columnName;
        }
        renameColumnInMenus(oldName, columnName);
        metadataColumn(oldName).name = columnName;

        refreshAfterMetadataEdit({columnsChanged: true});
        setStatus('Renamed column "' + oldName + '" to "' + columnName + '".');
        return true;
    }

    function deleteMetadataColumn(columnName) {
        if (!columnIsEditable(columnName)) { return false; }
        // Dropping a column that came from the source data is the one destructive,
        // un-undoable edit in the viewer -- and it only really bites once the
        // session is saved over the original. Columns created here are cheap to
        // recreate, so they go without a prompt.
        if (!columnIsViewerAdded(columnName) &&
                !window.confirm('Delete the column "' + columnName + '"?\n\n' +
                    'It came from this bundle\'s metadata, not from the viewer. ' +
                    'Reload the bundle to get it back; a session saved after this will not contain it.')) {
            return false;
        }
        const columnIndex = state.metadataColumnIndexByName.get(columnName);
        state.metadataColumns.splice(columnIndex, 1);
        state.metadataByNodeIndex.forEach(row => { row.splice(columnIndex, 1); });
        state.metadataColumnWidths.delete(columnName);
        // Palettes are stored under paletteKey(), which suffixes a numeric column
        // that is being colored categorically; clear both spellings.
        delete state.customPalettes[columnName];
        delete state.customPalettes[columnName + '\u0000categorical'];
        state.categoricalColumns.delete(columnName);
        if (state.metadataSort.columnKey === columnName) {
            state.metadataSort = {columnKey: null, direction: 'asc'};
        }
        refreshAfterMetadataEdit({columnsChanged: true});
        setStatus('Deleted column "' + columnName + '".');
        return true;
    }

    // ---- bulk fill ----

    // The three fill targets. "selection" is the annotation workflow (select a
    // cluster on the canvas, label it); "rows" uses the table's staged rows;
    // "page" matches what the per-column Copy button copies.
    function bulkFillTargetNodeIndices(target) {
        if (target === 'rows') {
            return Array.from(state.selectedMetadataNodeIndices).sort((left, right) => left - right);
        }
        if (target === 'page') {
            return state.renderedNodeIndices.slice();
        }
        if (target === 'all') {
            return state.allNodeIndices.slice();
        }
        return Array.from(state.selectedNodeIndices).sort((left, right) => left - right);
    }

    function bulkFillTargetLabel(target) {
        if (target === 'rows') { return 'staged table rows'; }
        if (target === 'page') { return 'rows on this page'; }
        if (target === 'all') { return 'nodes in the network'; }
        return 'selected nodes';
    }

    function applyBulkFill() {
        if (!state.bundle) { return; }
        const columnName = document.getElementById('metadata-fill-column').value;
        const target = document.getElementById('metadata-fill-target').value;
        const rawText = document.getElementById('metadata-fill-value').value;
        if (!columnIsEditable(columnName)) {
            setStatus('Choose a column to fill.');
            return;
        }
        const nodeIndices = bulkFillTargetNodeIndices(target);
        if (nodeIndices.length === 0) {
            setStatus('No ' + bulkFillTargetLabel(target) + ' to fill.');
            return;
        }
        const parsed = parseMetadataInput(columnName, rawText);
        if (!parsed.ok) {
            setStatus('Cannot fill "' + columnName + '": ' + parsed.reason + '.');
            return;
        }
        nodeIndices.forEach(nodeIndex => writeMetadataValue(nodeIndex, columnName, parsed.value));
        const widened = applyColumnWidening(columnName, parsed);
        refreshAfterMetadataEdit({columnsChanged: widened});
        setStatus('Set "' + columnName + '" on ' + nodeIndices.length.toLocaleString() + ' ' +
            bulkFillTargetLabel(target) + '.');
    }

    // Paste applies to the rendered page in display order -- the exact inverse
    // of the per-column Copy button. Counts must match: a short or long paste
    // almost always means the page, filter or sort moved since the copy, and
    // silently filling a prefix would mislabel rows.
    function applyPastedColumn() {
        if (!state.bundle) { return; }
        const columnName = document.getElementById('metadata-paste-column').value;
        const text = document.getElementById('metadata-paste-values').value;
        if (!columnIsEditable(columnName)) {
            setStatus('Choose a column to paste into.');
            return;
        }
        const lines = text.split(/\r?\n/);
        while (lines.length > 0 && lines[lines.length - 1].trim() === '') {
            lines.pop();
        }
        const nodeIndices = state.renderedNodeIndices;
        if (lines.length !== nodeIndices.length) {
            setStatus('Paste has ' + lines.length.toLocaleString() + ' values but this page shows ' +
                nodeIndices.length.toLocaleString() + ' rows; they must match.');
            return;
        }
        const parsedValues = [];
        let widenTo = null;
        for (let i = 0; i < lines.length; i++) {
            // Tolerate a pasted multi-column block by taking the first field.
            const parsed = parseMetadataInput(columnName, lines[i].split('\t')[0]);
            if (!parsed.ok) {
                setStatus('Cannot paste into "' + columnName + '" at line ' + (i + 1) + ': ' + parsed.reason + '.');
                return;
            }
            widenTo = widenTo || parsed.widenTo || null;
            parsedValues.push(parsed.value);
        }
        nodeIndices.forEach((nodeIndex, position) => {
            writeMetadataValue(nodeIndex, columnName, parsedValues[position]);
        });
        const widened = applyColumnWidening(columnName, {widenTo});
        refreshAfterMetadataEdit({columnsChanged: widened});
        closeMetadataPasteDialog();
        setStatus('Pasted ' + parsedValues.length.toLocaleString() + ' values into "' + columnName + '".');
    }

    // ---- editing controls ----

    function populateEditableColumnMenus() {
        const editable = state.metadataColumns.map(column => column.name);
        [
            document.getElementById('metadata-fill-column'),
            document.getElementById('metadata-paste-column'),
            document.getElementById('metadata-rename-column'),
        ].forEach(select => {
            const previous = select.value;
            select.innerHTML = '';
            editable.forEach(name => {
                const option = document.createElement('option');
                option.value = name;
                option.textContent = name;
                select.appendChild(option);
            });
            if (editable.includes(previous)) { select.value = previous; }
            select.disabled = editable.length === 0;
        });

        const deleteSelect = document.getElementById('metadata-delete-column');
        const previousDelete = deleteSelect.value;
        deleteSelect.innerHTML = '';
        state.metadataColumns.forEach(column => {
            const option = document.createElement('option');
            option.value = column.name;
            // Say where each column came from, since deleting a source column
            // discards data the viewer cannot regenerate.
            option.textContent = column.name +
                (column.origin === 'viewer' ? ' (added here)' : ' (from bundle)');
            deleteSelect.appendChild(option);
        });
        if (editable.includes(previousDelete)) { deleteSelect.value = previousDelete; }
        deleteSelect.disabled = editable.length === 0;
        document.getElementById('metadata-delete-column-apply').disabled = editable.length === 0;
        document.getElementById('metadata-fill-apply').disabled = editable.length === 0;
        document.getElementById('metadata-paste-open').disabled = editable.length === 0;
        document.getElementById('metadata-rename-apply').disabled = editable.length === 0;
        updateDeleteColumnNote();
        populateMetadataSelectMenus();
        updateMetadataSelectOperands();
        updateMetadataSelectNote();
    }

    function updateDeleteColumnNote() {
        const note = document.getElementById('metadata-delete-note');
        const columnName = document.getElementById('metadata-delete-column').value;
        if (!columnIsEditable(columnName)) {
            note.textContent = 'There are no columns to delete.';
        } else if (columnIsViewerAdded(columnName)) {
            note.textContent = 'Created in the viewer.';
        } else {
            note.textContent = 'From the bundle — reload it to restore this column.';
        }
    }


    // ---------------------------------------------------------------------
    // Select by value
    //
    // A column/comparison/value query that edits the node selection. Almost
    // everything downstream -- the per-column charts, "Set column", Export table
    // TSV, Save extraction -- is driven by state.selectedNodeIndices, and until
    // now the only way to build that set was to find the nodes on the canvas by
    // eye or to sort the table and tick rows.
    //
    // Two conventions, both surfaced in the panel's own note:
    //   * Text comparisons run against the *stored* value, not the table's
    //     formatted rendering of it, so a thousands separator never has to be
    //     typed to match the cell that displays as "1,234".
    //   * Ordering comparisons go through compareMetadataValues(), the very
    //     comparator the table's column sort uses -- so "greater than" means what
    //     sorting that column already showed, on text columns as well as numeric
    //     ones. Cells with no value never match anything.
    // ---------------------------------------------------------------------

    // node_id is offered alongside the metadata columns: it is the one field every
    // node has, and an accession prefix is a common way into a network.
    const METADATA_SELECT_NODE_ID = '__node_id__';

    // `operands` drives which value fields the panel shows, so how many bounds a
    // comparison takes is stated once instead of re-derived at each use.
    const METADATA_SELECT_OPS = [
        {id: 'contains', label: 'contains', operands: 1},
        {id: 'exact', label: 'is exactly', operands: 1},
        {id: 'gt', label: 'greater than', operands: 1},
        {id: 'lt', label: 'less than', operands: 1},
        {id: 'between', label: 'between', operands: 2},
        {id: 'regex', label: 'matches regex', operands: 1},
    ];

    function metadataSelectOp(opId) {
        return METADATA_SELECT_OPS.find(op => op.id === opId) || METADATA_SELECT_OPS[0];
    }

    function metadataSelectFieldLabel(field) {
        return field === METADATA_SELECT_NODE_ID ? 'node_id' : field;
    }

    function metadataSelectFieldValue(nodeIndex, field) {
        return field === METADATA_SELECT_NODE_ID
            ? nodeId(nodeIndex)
            : metadataValue(nodeIndex, field);
    }

    // Returns {test} on success, or {error} carrying a message for the panel note.
    // `test` takes the stored value and its String() form, since the substring and
    // regex comparisons want text while the ordering ones want the value itself.
    function metadataSelectMatcher(opId, firstText, secondText) {
        const first = String(firstText ?? '').trim();
        const second = String(secondText ?? '').trim();
        if (first === '') {
            return {error: 'Type a value to compare against.'};
        }
        if (opId === 'regex') {
            let pattern;
            try {
                // Case-insensitive, and unanchored, so it reads like the search box
                // rather than like a fullmatch.
                pattern = new RegExp(first, 'i');
            } catch (error) {
                return {error: 'Not a valid regular expression: ' + error.message};
            }
            return {test: (raw, text) => pattern.test(text)};
        }
        if (opId === 'contains') {
            const needle = first.toLowerCase();
            return {test: (raw, text) => text.toLowerCase().includes(needle)};
        }
        if (opId === 'exact') {
            const wanted = first.toLowerCase();
            // Numeric equality as well as text, so "5.0" and "5" both match a float
            // column's 5 and the match does not depend on how the number was typed.
            const wantedNumber = Number(first);
            const numeric = Number.isFinite(wantedNumber);
            return {test: (raw, text) => text.toLowerCase() === wanted
                || (numeric && typeof raw === 'number' && raw === wantedNumber)};
        }
        if (opId === 'gt') {
            return {test: raw => compareMetadataValues(raw, first) > 0};
        }
        if (opId === 'lt') {
            return {test: raw => compareMetadataValues(raw, first) < 0};
        }
        if (opId === 'between') {
            if (second === '') {
                return {error: 'Type both bounds, or pick another comparison.'};
            }
            // Inclusive at both ends, and tolerant of bounds typed the wrong way
            // round -- ordering them by the same comparator that will test them.
            const reversed = compareMetadataValues(first, second) > 0;
            const low = reversed ? second : first;
            const high = reversed ? first : second;
            return {test: raw => compareMetadataValues(raw, low) >= 0
                && compareMetadataValues(raw, high) <= 0};
        }
        return {error: 'Unknown comparison.'};
    }

    function metadataSelectMatches(nodeIndices, field, matcher) {
        return nodeIndices.filter(nodeIndex => {
            const raw = metadataSelectFieldValue(nodeIndex, field);
            // A blank cell is an absence, not a value: it matches no comparison,
            // including "contains ''" (which the matcher rejects outright anyway).
            if (isMissingMetadataValue(raw)) {
                return false;
            }
            return matcher.test(raw, String(raw));
        });
    }

    function setMetadataSelectNote(text) {
        document.getElementById('metadata-select-note').textContent = text;
    }

    // Says how the chosen column will be compared, which is the one thing about
    // this panel that is not visible from its controls.
    function updateMetadataSelectNote() {
        const field = document.getElementById('metadata-select-column').value;
        if (!field) {
            setMetadataSelectNote('');
            return;
        }
        const columnType = field === METADATA_SELECT_NODE_ID ? 'string' : metadataColumnType(field);
        const kind = (columnType === 'int' || columnType === 'float') ? 'numbers' : 'text';
        setMetadataSelectNote('Compares ' + metadataSelectFieldLabel(field) + ' as ' + kind +
            '; blank cells never match.');
    }

    function updateMetadataSelectOperands() {
        const op = metadataSelectOp(document.getElementById('metadata-select-op').value);
        const twoOperands = op.operands === 2;
        document.getElementById('metadata-select-and').hidden = !twoOperands;
        document.getElementById('metadata-select-value2').hidden = !twoOperands;
        document.getElementById('metadata-select-value').placeholder = twoOperands ? 'lower bound' : 'value';
    }

    function populateMetadataSelectMenus() {
        const opSelect = document.getElementById('metadata-select-op');
        if (opSelect.options.length === 0) {
            METADATA_SELECT_OPS.forEach(op => {
                const option = document.createElement('option');
                option.value = op.id;
                option.textContent = op.label;
                opSelect.appendChild(option);
            });
        }
        const select = document.getElementById('metadata-select-column');
        const previous = select.value;
        select.innerHTML = '';
        const nodeIdOption = document.createElement('option');
        nodeIdOption.value = METADATA_SELECT_NODE_ID;
        nodeIdOption.textContent = 'node_id';
        select.appendChild(nodeIdOption);
        state.metadataColumns.forEach(column => {
            const option = document.createElement('option');
            option.value = column.name;
            option.textContent = column.name;
            select.appendChild(option);
        });
        const available = Array.from(select.options).map(option => option.value);
        select.value = available.includes(previous) ? previous : METADATA_SELECT_NODE_ID;
    }

    // `mode` is 'add' (search the whole network and union the matches in),
    // 'remove' (search the current selection and drop the matches) or 'subset'
    // (search the current selection and keep only the matches). Only 'add' can
    // grow the selection; that asymmetry is what makes the three buttons compose
    // into an intersection of queries.
    function applyMetadataSelectByValue(mode) {
        if (!state.bundle) {
            return false;
        }
        const field = document.getElementById('metadata-select-column').value;
        const opId = document.getElementById('metadata-select-op').value;
        const matcher = metadataSelectMatcher(
            opId,
            document.getElementById('metadata-select-value').value,
            document.getElementById('metadata-select-value2').value,
        );
        if (matcher.error) {
            setMetadataSelectNote(matcher.error);
            return false;
        }
        const scope = mode === 'add'
            ? state.allNodeIndices
            : Array.from(state.selectedNodeIndices).sort((left, right) => left - right);
        const matches = metadataSelectMatches(scope, field, matcher);
        const before = state.selectedNodeIndices.size;
        if (mode === 'add') {
            matches.forEach(nodeIndex => state.selectedNodeIndices.add(nodeIndex));
        } else if (mode === 'remove') {
            matches.forEach(nodeIndex => state.selectedNodeIndices.delete(nodeIndex));
        } else {
            state.selectedNodeIndices = new Set(matches);
        }
        const verb = {add: 'Matched', remove: 'Removed', subset: 'Kept'}[mode];
        const scopeLabel = mode === 'add' ? 'the network' : 'the selection';
        setMetadataSelectNote(verb + ' ' + matches.length.toLocaleString() + ' of ' +
            scope.length.toLocaleString() + ' node' + (scope.length === 1 ? '' : 's') +
            ' in ' + scopeLabel + '.');
        setStatus(metadataSelectFieldLabel(field) + ' ' + metadataSelectOp(opId).label +
            ': ' + verb.toLowerCase() + ' ' + matches.length.toLocaleString() + ' node' +
            (matches.length === 1 ? '' : 's') + '; ' +
            state.selectedNodeIndices.size.toLocaleString() + ' now selected (was ' +
            before.toLocaleString() + ').');
        resetMetadataPage();
        repaintNow();
        updateMetadataTable();
        return true;
    }

    // One panel at a time: showing every editing control at once made the panel
    // busy, and these are occasional actions rather than per-row ones.
    const METADATA_EDIT_PANELS = ['add', 'set', 'rename', 'delete', 'select'];

    function toggleMetadataEditPanel(name) {
        const wasOpen = document.getElementById('metadata-panel-' + name)
            .getAttribute('aria-expanded') === 'true';
        METADATA_EDIT_PANELS.forEach(panelName => {
            const open = !wasOpen && panelName === name;
            document.getElementById('metadata-panel-' + panelName)
                .setAttribute('aria-expanded', open ? 'true' : 'false');
            document.getElementById('metadata-' + panelName + '-panel').hidden = !open;
        });
        if (wasOpen) { return; }
        populateEditableColumnMenus();
        const focusTarget = {
            add: 'metadata-new-column-name',
            set: 'metadata-fill-value',
            rename: 'metadata-rename-value',
            delete: 'metadata-delete-column',
            select: 'metadata-select-value',
        }[name];
        document.getElementById(focusTarget).focus();
    }

    function closeMetadataEditPanels() {
        METADATA_EDIT_PANELS.forEach(panelName => {
            document.getElementById('metadata-panel-' + panelName).setAttribute('aria-expanded', 'false');
            document.getElementById('metadata-' + panelName + '-panel').hidden = true;
        });
    }

    function applyColumnRename() {
        const input = document.getElementById('metadata-rename-value');
        if (renameMetadataColumn(document.getElementById('metadata-rename-column').value, input.value)) {
            input.value = '';
        }
    }

    function openMetadataPasteDialog() {
        populateEditableColumnMenus();
        document.getElementById('metadata-paste-count').textContent =
            state.renderedNodeIndices.length.toLocaleString();
        document.getElementById('metadata-paste-values').value = '';
        document.getElementById('metadata-paste-overlay').hidden = false;
        document.getElementById('metadata-paste-values').focus();
    }

    function closeMetadataPasteDialog() {
        document.getElementById('metadata-paste-overlay').hidden = true;
    }

    function metadataPasteDialogIsOpen() {
        return !document.getElementById('metadata-paste-overlay').hidden;
    }

    function setupMetadataEditing() {
        METADATA_EDIT_PANELS.forEach(panelName => {
            document.getElementById('metadata-panel-' + panelName)
                .addEventListener('click', () => toggleMetadataEditPanel(panelName));
        });
        document.getElementById('metadata-delete-column')
            .addEventListener('change', updateDeleteColumnNote);
        document.getElementById('metadata-select-column')
            .addEventListener('change', updateMetadataSelectNote);
        document.getElementById('metadata-select-op').addEventListener('change', () => {
            updateMetadataSelectOperands();
            updateMetadataSelectNote();
        });
        ['add', 'remove', 'subset'].forEach(mode => {
            document.getElementById('metadata-select-' + mode)
                .addEventListener('click', () => applyMetadataSelectByValue(mode));
        });
        ['metadata-select-value', 'metadata-select-value2'].forEach(inputId => {
            document.getElementById(inputId).addEventListener('keydown', event => {
                if (event.key === 'Enter') {
                    event.preventDefault();
                    // Enter runs the widening action; the two narrowing ones are
                    // destructive enough to be worth an explicit click.
                    applyMetadataSelectByValue('add');
                }
            });
        });
        document.getElementById('metadata-rename-apply').addEventListener('click', applyColumnRename);
        document.getElementById('metadata-rename-value').addEventListener('keydown', event => {
            if (event.key === 'Enter') {
                event.preventDefault();
                applyColumnRename();
            }
        });

        const tbody = document.querySelector('#metadata-table tbody');
        tbody.addEventListener('dblclick', event => {
            const cell = event.target.closest('td[data-column-key]');
            if (cell) {
                event.preventDefault();
                beginMetadataCellEdit(cell);
            }
        });

        document.getElementById('metadata-add-column').addEventListener('click', () => {
            const nameInput = document.getElementById('metadata-new-column-name');
            if (addMetadataColumn(nameInput.value, document.getElementById('metadata-new-column-type').value)) {
                nameInput.value = '';
            }
        });
        document.getElementById('metadata-add-cluster-column').addEventListener('click', () => {
            const nameInput = document.getElementById('metadata-new-column-name');
            if (addClusterColumn(nameInput.value)) {
                nameInput.value = '';
            }
        });
        document.getElementById('metadata-new-column-name').addEventListener('keydown', event => {
            if (event.key === 'Enter') {
                event.preventDefault();
                document.getElementById('metadata-add-column').click();
            }
        });
        document.getElementById('metadata-delete-column-apply').addEventListener('click', () => {
            deleteMetadataColumn(document.getElementById('metadata-delete-column').value);
        });
        document.getElementById('metadata-fill-apply').addEventListener('click', applyBulkFill);
        document.getElementById('metadata-fill-value').addEventListener('keydown', event => {
            if (event.key === 'Enter') {
                event.preventDefault();
                applyBulkFill();
            }
        });
        document.getElementById('metadata-paste-open').addEventListener('click', openMetadataPasteDialog);
        document.getElementById('metadata-paste-cancel').addEventListener('click', closeMetadataPasteDialog);
        document.getElementById('metadata-paste-apply').addEventListener('click', applyPastedColumn);
        document.getElementById('metadata-paste-overlay').addEventListener('click', event => {
            if (event.target === event.currentTarget) { closeMetadataPasteDialog(); }
        });
        document.addEventListener('keydown', event => {
            if (event.key === 'Escape' && metadataPasteDialogIsOpen()) {
                closeMetadataPasteDialog();
            }
        });
    }
