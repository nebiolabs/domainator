/*
 * Ten keyboard-addressable selection presets with a hover preview.
 *
 * Included into viewer.js's script scope; see session_state.js.
 */

    // ---------------------------------------------------------------------
    // Selection presets
    //
    // Ten slots above the hierarchy view. Press 0-9 to recall a slot (or click
    // it), Shift+0-9 to store the current selection into one, and hover a slot
    // to outline its nodes in gray without disturbing the live selection.
    //
    // Storing is bound to Shift rather than Ctrl/Cmd because browsers reserve
    // Ctrl/Cmd+1-9 for tab switching and Ctrl/Cmd+0 for zoom reset; a page
    // cannot intercept those. Slots are keyed off event.code (Digit0..Digit9)
    // so neither the keyboard layout nor Shift's symbol mapping matters.
    // ---------------------------------------------------------------------

    function presetPreviewNodeSet() {
        if (state.presetPreviewSlot === null) { return null; }
        const preset = state.selectionPresets.get(state.presetPreviewSlot);
        return preset && preset.nodeIndices.size > 0 ? preset.nodeIndices : null;
    }

    function renderPresetSlots() {
        const container = document.getElementById('preset-slots');
        if (!container) { return; }
        container.querySelectorAll('[data-preset-slot]').forEach(button => {
            const slot = Number(button.dataset.presetSlot);
            const preset = state.selectionPresets.get(slot);
            const count = preset ? preset.nodeIndices.size : 0;
            button.classList.toggle('preset-slot-filled', count > 0);
            button.disabled = !state.bundle;
            button.title = count > 0
                ? 'Preset ' + slot + ': ' + count.toLocaleString() + ' nodes.' +
                  ' Click (or press ' + slot + ') to replace the selection with it;' +
                  ' Shift-click to add it; Alt+Shift-click to subtract it;' +
                  ' hover to preview; Shift+' + slot + ' to overwrite it.'
                : 'Preset ' + slot + ' is empty. Select nodes, then press Shift+' + slot + ' to store them here.';
        });
    }

    function storeSelectionPreset(slot) {
        if (!state.bundle) { return; }
        if (state.selectedNodeIndices.size === 0) {
            // Storing nothing is how you clear a slot.
            state.selectionPresets.delete(slot);
            setStatus('Cleared preset ' + slot + ' (nothing was selected).');
        } else {
            state.selectionPresets.set(slot, {
                nodeIndices: new Set(state.selectedNodeIndices),
                savedAt: new Date().toISOString(),
            });
            setStatus('Stored ' + state.selectedNodeIndices.size.toLocaleString() +
                ' nodes in preset ' + slot + '.');
        }
        renderPresetSlots();
        scheduleClusterRender();
    }

    // mode: 'replace' (the default), 'add' to union the preset into the current
    // selection, or 'subtract' to remove it. Shift/Alt+Shift on a slot button pick
    // the latter two, matching the canvas box-select gesture where Shift selects
    // and adding Alt deselects.
    function recallSelectionPreset(slot, mode = 'replace') {
        if (!state.bundle) { return; }
        const preset = state.selectionPresets.get(slot);
        if (!preset) {
            setStatus('Preset ' + slot + ' is empty.');
            return;
        }
        const before = state.selectedNodeIndices.size;
        if (mode === 'add') {
            preset.nodeIndices.forEach(nodeIndex => state.selectedNodeIndices.add(nodeIndex));
        } else if (mode === 'subtract') {
            preset.nodeIndices.forEach(nodeIndex => state.selectedNodeIndices.delete(nodeIndex));
        } else {
            state.selectedNodeIndices = new Set(preset.nodeIndices);
        }
        const after = state.selectedNodeIndices.size;
        resetMetadataPage();
        clearMetadataRowSelection();
        if (mode === 'add') {
            setStatus('Added preset ' + slot + ' to the selection (+' +
                (after - before).toLocaleString() + ' nodes, ' + after.toLocaleString() + ' total).');
        } else if (mode === 'subtract') {
            setStatus('Subtracted preset ' + slot + ' from the selection (-' +
                (before - after).toLocaleString() + ' nodes, ' + after.toLocaleString() + ' remaining).');
        } else {
            setStatus('Recalled preset ' + slot + ' (' + after.toLocaleString() + ' nodes).');
        }
        renderClusterView();
        updateMetadataTable();
    }

    function presetClickMode(event) {
        if (!event.shiftKey) { return 'replace'; }
        return event.altKey ? 'subtract' : 'add';
    }

    function setPresetPreviewSlot(slot) {
        if (state.presetPreviewSlot === slot) { return; }
        state.presetPreviewSlot = slot;
        scheduleClusterRender();
    }

    // True when a keystroke belongs to a control the user is typing into, in
    // which case the bare-digit shortcuts must not fire.
    function keyboardTargetIsTextEntry(target) {
        if (!target || !target.tagName) { return false; }
        if (target.isContentEditable) { return true; }
        return ['INPUT', 'SELECT', 'TEXTAREA', 'OPTION'].includes(target.tagName);
    }

    function handlePresetKeydown(event) {
        if (!state.bundle) { return; }
        if (event.ctrlKey || event.metaKey || event.altKey) { return; }
        if (keyboardTargetIsTextEntry(event.target)) { return; }
        if (colorPickerIsOpen() || metadataPasteDialogIsOpen() || nameDialogIsOpen()
            || columnChartIsOpen()) { return; }
        const match = /^Digit([0-9])$/.exec(event.code || '');
        if (!match) { return; }
        const slot = Number(match[1]);
        event.preventDefault();
        if (event.shiftKey) {
            storeSelectionPreset(slot);
        } else {
            recallSelectionPreset(slot);
        }
    }

    function setupSelectionPresets() {
        const container = document.getElementById('preset-slots');
        container.addEventListener('click', event => {
            const button = event.target.closest('[data-preset-slot]');
            if (button) {
                recallSelectionPreset(Number(button.dataset.presetSlot), presetClickMode(event));
            }
        });
        // mouseover/mouseout (not mouseenter/mouseleave) so one delegated pair
        // covers all ten buttons.
        container.addEventListener('mouseover', event => {
            const button = event.target.closest('[data-preset-slot]');
            if (button) { setPresetPreviewSlot(Number(button.dataset.presetSlot)); }
        });
        container.addEventListener('mouseout', event => {
            const button = event.target.closest('[data-preset-slot]');
            if (button) { setPresetPreviewSlot(null); }
        });
        // Keyboard focus should preview too, so the slots are usable without a mouse.
        container.addEventListener('focusin', event => {
            const button = event.target.closest('[data-preset-slot]');
            if (button) { setPresetPreviewSlot(Number(button.dataset.presetSlot)); }
        });
        container.addEventListener('focusout', () => setPresetPreviewSlot(null));
        document.addEventListener('keydown', handlePresetKeydown);
        renderPresetSlots();
    }
