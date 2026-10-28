    function defaultMetadataColumnWidth(columnKey) {
        if (columnKey === 'node_id') {
            return 230;
        }
        const column = metadataColumn(columnKey);
        if (!column) {
            return 180;
        }
        if (column.type === 'int' || column.type === 'float') {
            return 140;
        }
        if (column.type === 'bool' || column.type === 'boolean') {
            return 120;
        }
        return 190;
    }

    function metadataColumnWidth(columnKey) {
        return state.metadataColumnWidths.get(columnKey) || defaultMetadataColumnWidth(columnKey);
    }

    function applyMetadataColumnWidths() {
        const colgroup = document.querySelector('#metadata-table colgroup');
        if (!colgroup) {
            return;
        }
        colgroup.innerHTML = '';
        metadataColumnKeys().forEach(columnKey => {
            const col = document.createElement('col');
            col.style.width = metadataColumnWidth(columnKey) + 'px';
            colgroup.appendChild(col);
        });
    }

    function updateMetadataColumnWidth(columnKey, nextWidth) {
        const clampedWidth = Math.max(96, Math.min(640, Math.round(nextWidth)));
        state.metadataColumnWidths.set(columnKey, clampedWidth);
        applyMetadataColumnWidths();
    }

    function startMetadataColumnResize(columnKey, event) {
        event.preventDefault();
        event.stopPropagation();
        state.metadataResize = {
            columnKey,
            startX: event.clientX,
            startWidth: metadataColumnWidth(columnKey),
        };
        document.body.style.cursor = 'col-resize';
        document.body.style.userSelect = 'none';
    }

    function setupMetadataTableDelegation() {
        const thead = document.querySelector('#metadata-table thead');
        const tbody = document.querySelector('#metadata-table tbody');

        thead.addEventListener('click', event => {
            // Checked before the sort button: the chart glyph must never be
            // given a data-column-key, or clicking it would also sort.
            const chartBtn = event.target.closest('[data-chart-column]');
            if (chartBtn) {
                event.preventDefault();
                event.stopPropagation();
                openColumnChartMenu(chartBtn.dataset.chartColumn, chartBtn);
                return;
            }
            const sortBtn = event.target.closest('[data-column-key]');
            if (sortBtn) { toggleMetadataSort(sortBtn.dataset.columnKey); return; }
            const copyBtn = event.target.closest('[data-copy-column]');
            if (copyBtn) { event.preventDefault(); event.stopPropagation(); copyMetadataColumn(copyBtn.dataset.copyColumn); }
        });

        thead.addEventListener('pointerdown', event => {
            const resizeHandle = event.target.closest('[data-resize-column]');
            if (resizeHandle) { startMetadataColumnResize(resizeHandle.dataset.resizeColumn, event); }
        });

        tbody.addEventListener('mousedown', event => {
            if (event.shiftKey && !event.target.closest('button, a, input, select, textarea')) {
                event.preventDefault();
            }
        });

        tbody.addEventListener('click', event => {
            if (event.target.closest('button, a, input, select, textarea')) { return; }
            // The second click of a double-click opens the cell editor; letting it
            // through would toggle the row selection on the way in.
            if (event.detail > 1) { return; }
            if (!event.shiftKey && metadataTableHasActiveTextSelection()) { return; }
            const row = event.target.closest('tr[data-node-index]');
            if (!row) { return; }
            toggleMetadataRowSelection(Number(row.dataset.nodeIndex), state.renderedNodeIndices, {range: event.shiftKey});
        });

        tbody.addEventListener('keydown', event => {
            if (event.key !== 'Enter' && event.key !== ' ') { return; }
            if (event.target.closest('button, a, input, select, textarea')) { return; }
            const row = event.target.closest('tr[data-node-index]');
            if (!row) { return; }
            event.preventDefault();
            toggleMetadataRowSelection(Number(row.dataset.nodeIndex), state.renderedNodeIndices, {range: event.shiftKey});
        });
    }

    function updateMetadataTable() {
        // The table body is re-rendered wholesale, so any in-progress cell edit
        // is destroyed along with its <input>. The header goes the same way, so
        // an open chart menu would be left anchored to a detached button.
        state.metadataEditCell = null;
        closeColumnChartMenu();
        const selected = Array.from(state.selectedNodeIndices).sort((left, right) => left - right);
        const baseNodeIndices = metadataBaseNodeIndices();
        const filteredNodeIndices = filteredMetadataNodeIndices(baseNodeIndices);
        const sortedNodeIndices = sortedMetadataNodeIndices(filteredNodeIndices);
        const pagination = metadataPagination(sortedNodeIndices.length);
        state.metadataPage = pagination.pageIndex;
        const renderedNodeIndices = sortedNodeIndices.slice(pagination.start, pagination.end);
        state.renderedNodeIndices = renderedNodeIndices;
        pruneMetadataRowSelection(renderedNodeIndices);
        applyMetadataColumnWidths();
        const thead = document.querySelector('#metadata-table thead');
        const tbody = document.querySelector('#metadata-table tbody');

        const columnKeys = metadataColumnKeys();
        const headerCells = columnKeys.map(label => {
            const escapedLabel = htmlEscape(label);
            const ariasort = state.metadataSort.columnKey === label
                ? (state.metadataSort.direction === 'asc' ? 'ascending' : 'descending')
                : 'none';
            const indicator = htmlEscape(metadataSortIndicator(label));
            return '<th class="metadata-header" aria-sort="' + ariasort + '">'
                + '<div class="metadata-header-cell">'
                + '<button type="button" class="metadata-sort-button" data-column-key="' + escapedLabel + '" title="Sort by ' + escapedLabel + '">'
                + '<span>' + escapedLabel + '</span><span class="metadata-sort-indicator">' + indicator + '</span></button>'
                // No chart glyph on node_id: it is the graph key, unique by
                // construction, so every chart of it is degenerate (N rows of
                // count 1, N equal slices). A disabled button would only invite
                // a "why?" with no answer.
                + (label === 'node_id' ? ''
                    : '<button type="button" class="metadata-chart-button" data-chart-column="' + escapedLabel
                        + '" aria-haspopup="menu" aria-expanded="false" title="Charts and tables for ' + escapedLabel + '">'
                        + '<svg class="metadata-chart-glyph" viewBox="0 0 12 12" aria-hidden="true" focusable="false">'
                        + '<rect x="1" y="6" width="2.6" height="5"></rect>'
                        + '<rect x="4.7" y="3" width="2.6" height="8"></rect>'
                        + '<rect x="8.4" y="1" width="2.6" height="10"></rect></svg></button>')
                + '<button type="button" class="metadata-copy-button" data-copy-column="' + escapedLabel + '" title="Copy the currently displayed values from ' + escapedLabel + '">Copy</button>'
                + '</div>'
                + '<div class="metadata-resize-handle" role="separator" aria-orientation="vertical" title="Drag to resize column" data-resize-column="' + escapedLabel + '"></div>'
                + '</th>';
        });
        thead.innerHTML = '<tr>' + headerCells.join('') + '</tr>';

        let bodyHtml;
        if (sortedNodeIndices.length === 0) {
            const emptyMsg = htmlEscape(metadataFilterText()
                ? 'No metadata rows match the current filter.'
                : 'No metadata rows to display.');
            bodyHtml = '<tr class="metadata-empty-row"><td colspan="' + (state.metadataColumns.length + 1) + '">' + emptyMsg + '</td></tr>';
        } else {
            bodyHtml = renderedNodeIndices.map(nodeIndex => {
                const isSelected = state.selectedMetadataNodeIndices.has(nodeIndex);
                const rowClass = isSelected ? ' class="metadata-row-selected"' : '';
                const ariaSelected = isSelected ? 'true' : 'false';
                const idText = htmlEscape(nodeId(nodeIndex));
                let cellsHtml = '<td title="' + idText + '">' + idText + '</td>';
                for (const column of state.metadataColumns) {
                    const value = metadataValue(nodeIndex, column.name);
                    const cellClass = metadataCellClass(column.name, value);
                    const escapedName = htmlEscape(column.name);
                    const displayText = htmlEscape(formatMetadataDisplayValue(column.name, value));
                    cellsHtml += '<td data-column-key="' + escapedName + '"'
                        + ' class="metadata-cell-editable' + (cellClass ? ' ' + cellClass : '') + '"'
                        + ' title="' + displayText + ' (double-click to edit)">' + displayText + '</td>';
                }
                return '<tr' + rowClass + ' tabindex="0" aria-selected="' + ariaSelected + '" data-node-index="' + nodeIndex + '">' + cellsHtml + '</tr>';
            }).join('');
        }
        tbody.innerHTML = bodyHtml;

        document.getElementById('selection-summary').textContent = selected.length.toLocaleString() + ' nodes selected';
        const pageStatus = document.getElementById('metadata-page-status');
        if (sortedNodeIndices.length === 0) {
            pageStatus.textContent = metadataFilterText()
                ? 'No rows match the current filter.'
                : 'No metadata rows available.';
        } else if (pagination.showAll) {
            pageStatus.textContent = 'Showing all ' + sortedNodeIndices.length.toLocaleString() + ' rows on one page.';
        } else {
            pageStatus.textContent = 'Showing rows ' + (pagination.start + 1).toLocaleString() + ' to ' + pagination.end.toLocaleString() +
                ' of ' + sortedNodeIndices.length.toLocaleString() + ' (page ' + (pagination.pageIndex + 1).toLocaleString() +
                ' of ' + pagination.pageCount.toLocaleString() + ').';
        }
        document.getElementById('metadata-prev-page').disabled = !state.bundle || sortedNodeIndices.length === 0 || pagination.pageIndex === 0;
        document.getElementById('metadata-next-page').disabled = !state.bundle || sortedNodeIndices.length === 0 || pagination.pageIndex >= pagination.pageCount - 1;
        const sortDescription = metadataSortDescription();
        const filterText = metadataFilterText();
        const filterDescription = filterText ? ' matching filter "' + filterText + '"' : '';
        const pagerActive = pagination.pageCount > 1;
        state.metadataBaseNoteText = activeView().metadataNoteText({
            selectedCount: selected.length, pagerActive, filterDescription, sortDescription,
        });
        document.getElementById('export-selected').disabled = !state.bundle;
        document.getElementById('metadata-reset-sort').disabled = !state.bundle || !state.metadataSort.columnKey;
        document.getElementById('clear-selection').disabled = selected.length === 0;
        activeView().updateSelectionControls(selected.length);
        applyMetadataTableRowHighlights();
        // Charts summarize the rows this function just resolved, so a changed
        // selection, filter or edit has to reach an open chart.
        refreshColumnChartIfOpen();
    }

    function applyMetadataTableRowHighlights() {
        const tbody = document.querySelector('#metadata-table tbody');
        if (tbody) {
            tbody.querySelectorAll('tr[data-node-index]').forEach(row => {
                const idx = Number(row.dataset.nodeIndex);
                const sel = state.selectedMetadataNodeIndices.has(idx);
                row.className = sel ? 'metadata-row-selected' : '';
                row.setAttribute('aria-selected', sel ? 'true' : 'false');
            });
        }
        const metadataSelectionCount = state.selectedMetadataNodeIndices.size;
        const metadataSelectionDescription = metadataSelectionCount > 0
            ? metadataSelectionCount.toLocaleString() + ' table rows selected. Click Select nodes to promote them into the graph selection. '
            : '';
        document.getElementById('selection-note').textContent = metadataSelectionDescription + state.metadataBaseNoteText;
        document.getElementById('metadata-select-nodes').disabled = !state.bundle || metadataSelectionCount === 0;
        document.getElementById('metadata-deselect-rows').disabled = !state.bundle || metadataSelectionCount === 0;
        // "Add to selection" searches the whole network, so it works from a standing
        // start; the two narrowing actions have nothing to narrow without a selection.
        document.getElementById('metadata-select-add').disabled = !state.bundle;
        document.getElementById('metadata-select-remove').disabled = !state.bundle || state.selectedNodeIndices.size === 0;
        document.getElementById('metadata-select-subset').disabled = !state.bundle || state.selectedNodeIndices.size === 0;
    }

    function metadataSortIndicator(columnKey) {
        if (state.metadataSort.columnKey !== columnKey) {
            return '↕';
        }
        return state.metadataSort.direction === 'asc' ? '↑' : '↓';
    }

    function metadataSortDescription() {
        if (!state.metadataSort.columnKey) {
            return null;
        }
        return state.metadataSort.columnKey + ' (' + (state.metadataSort.direction === 'asc' ? 'ascending' : 'descending') + ')';
    }

    function metadataFilterText() {
        const input = document.getElementById('metadata-filter');
        return input ? input.value.trim().toLowerCase() : '';
    }

    function metadataNullPlacement() {
        const select = document.getElementById('metadata-null-order');
        return select ? (select.value || 'last') : 'last';
    }

    function isMissingMetadataValue(value) {
        return value === null || value === undefined || value === '';
    }

    function metadataSortValue(nodeIndex, columnKey) {
        if (columnKey === 'node_id') {
            return nodeId(nodeIndex);
        }
        const columnIndex = state.metadataColumnIndexByName.get(columnKey);
        if (columnIndex === undefined) {
            return null;
        }
        return state.metadataByNodeIndex[nodeIndex]?.[columnIndex] ?? null;
    }

    function metadataMatchesFilter(nodeIndex, filterText) {
        if (!filterText) {
            return true;
        }
        return (state.metadataSearchTextByNodeIndex[nodeIndex] || '').includes(filterText);
    }

    function scheduleMetadataFilterUpdate(delayMs = 120) {
        if (state.pendingMetadataFilterTimer !== null) {
            window.clearTimeout(state.pendingMetadataFilterTimer);
        }
        state.pendingMetadataFilterTimer = window.setTimeout(() => {
            state.pendingMetadataFilterTimer = null;
            resetMetadataPage();
            updateMetadataTable();
        }, delayMs);
    }

    function filteredMetadataNodeIndices(nodeIndices) {
        const filterText = metadataFilterText();
        if (!filterText) {
            return [...nodeIndices];
        }
        return nodeIndices.filter(nodeIndex => metadataMatchesFilter(nodeIndex, filterText));
    }

    function metadataDisplayNodeIndices(nodeIndices) {
        return sortedMetadataNodeIndices(filteredMetadataNodeIndices(nodeIndices));
    }

    function metadataRowsPerPageSetting() {
        const select = document.getElementById('metadata-rows-per-page');
        const value = select ? select.value : String(DEFAULT_METADATA_PAGE_SIZE);
        if (value === 'all') {
            return {showAll: true, pageSize: Number.POSITIVE_INFINITY};
        }
        const pageSize = Number.parseInt(value, 10);
        if (!Number.isFinite(pageSize) || pageSize <= 0) {
            return {showAll: false, pageSize: DEFAULT_METADATA_PAGE_SIZE};
        }
        return {showAll: false, pageSize};
    }

    function metadataPagination(totalRowCount) {
        const rowsPerPage = metadataRowsPerPageSetting();
        if (rowsPerPage.showAll) {
            return {pageCount: 1, pageIndex: 0, start: 0, end: totalRowCount, showAll: true};
        }
        const pageCount = Math.max(1, Math.ceil(totalRowCount / rowsPerPage.pageSize));
        const pageIndex = Math.max(0, Math.min(state.metadataPage, pageCount - 1));
        const start = totalRowCount === 0 ? 0 : pageIndex * rowsPerPage.pageSize;
        const end = Math.min(totalRowCount, start + rowsPerPage.pageSize);
        return {pageCount, pageIndex, start, end, showAll: false};
    }

    function resetMetadataPage() {
        state.metadataPage = 0;
    }

    function stepMetadataPage(delta) {
        if (!state.bundle) {
            return;
        }
        const totalRowCount = metadataDisplayNodeIndices(metadataBaseNodeIndices()).length;
        const pagination = metadataPagination(totalRowCount);
        const nextPage = Math.max(0, Math.min(pagination.pageCount - 1, pagination.pageIndex + delta));
        if (nextPage === pagination.pageIndex) {
            return;
        }
        state.metadataPage = nextPage;
        updateMetadataTable();
    }

    function metadataBaseNodeIndices() {
        const selected = Array.from(state.selectedNodeIndices).sort((left, right) => left - right);
        if (selected.length > 0) {
            return selected;
        }
        return state.allNodeIndices;
    }

    function metadataColumnValues(nodeIndices, columnKey) {
        if (columnKey === 'node_id') {
            return nodeIndices.map(nodeIndex => nodeId(nodeIndex));
        }
        return nodeIndices.map(nodeIndex => formatMetadataDisplayValue(columnKey, metadataSortValue(nodeIndex, columnKey)));
    }

    async function copyMetadataColumn(columnKey) {
        if (!state.bundle) {
            return;
        }
        const nodeIndices = metadataDisplayNodeIndices(metadataBaseNodeIndices());
        const values = metadataColumnValues(nodeIndices, columnKey);
        try {
            await writeTextToClipboard(values.join('\n'));
            setStatus('Copied column "' + columnKey + '" for ' + values.length.toLocaleString() + ' rows.');
        } catch (error) {
            console.error(error);
            setStatus('Failed to copy column "' + columnKey + '": ' + error.message);
        }
    }

    function pruneMetadataRowSelection(visibleNodeIndices) {
        const visibleSet = new Set(visibleNodeIndices);
        if (state.metadataRowSelectionAnchor !== null && !visibleSet.has(state.metadataRowSelectionAnchor)) {
            state.metadataRowSelectionAnchor = null;
        }
        if (state.selectedMetadataNodeIndices.size === 0) {
            return;
        }
        state.selectedMetadataNodeIndices = new Set(
            Array.from(state.selectedMetadataNodeIndices).filter(nodeIndex => visibleSet.has(nodeIndex))
        );
    }

    function clearMetadataRowSelection() {
        state.selectedMetadataNodeIndices = new Set();
        state.metadataRowSelectionAnchor = null;
    }

    function metadataTableHasActiveTextSelection() {
        const selection = window.getSelection ? window.getSelection() : null;
        if (!selection || selection.rangeCount === 0 || selection.isCollapsed) {
            return false;
        }
        const metadataTable = document.getElementById('metadata-table');
        if (!metadataTable) {
            return false;
        }
        for (let rangeIndex = 0; rangeIndex < selection.rangeCount; rangeIndex += 1) {
            const range = selection.getRangeAt(rangeIndex);
            if (metadataTable.contains(range.commonAncestorContainer)) {
                return true;
            }
        }
        return false;
    }

    function toggleMetadataRowSelection(nodeIndex, renderedNodeIndices, options = {}) {
        const nextSelection = new Set(state.selectedMetadataNodeIndices);
        if (options.range && state.metadataRowSelectionAnchor !== null) {
            const anchorPosition = renderedNodeIndices.indexOf(state.metadataRowSelectionAnchor);
            const currentPosition = renderedNodeIndices.indexOf(nodeIndex);
            if (anchorPosition >= 0 && currentPosition >= 0) {
                const start = Math.min(anchorPosition, currentPosition);
                const end = Math.max(anchorPosition, currentPosition);
                const shouldSelect = !nextSelection.has(nodeIndex);
                for (let position = start; position <= end; position += 1) {
                    const rangeNodeIndex = renderedNodeIndices[position];
                    if (shouldSelect) {
                        nextSelection.add(rangeNodeIndex);
                    } else {
                        nextSelection.delete(rangeNodeIndex);
                    }
                }
            } else if (nextSelection.has(nodeIndex)) {
                nextSelection.delete(nodeIndex);
            } else {
                nextSelection.add(nodeIndex);
            }
        } else if (nextSelection.has(nodeIndex)) {
            nextSelection.delete(nodeIndex);
        } else {
            nextSelection.add(nodeIndex);
        }
        state.selectedMetadataNodeIndices = nextSelection;
        state.metadataRowSelectionAnchor = nodeIndex;
        applyMetadataTableRowHighlights();
    }

    function selectNodesFromMetadataRows() {
        if (state.selectedMetadataNodeIndices.size === 0) {
            return;
        }
        state.selectedNodeIndices = new Set(state.selectedMetadataNodeIndices);
        resetMetadataPage();
        clearMetadataRowSelection();
        repaintNow();
        updateMetadataTable();
    }

    function metadataColumnType(columnKey) {
        if (columnKey === 'node_id') {
            return 'string';
        }
        return metadataColumn(columnKey)?.type || 'string';
    }

    function formatMetadataDisplayValue(columnKey, value) {
        if (isMissingMetadataValue(value)) {
            return '—';
        }
        const columnType = metadataColumnType(columnKey);
        if (columnType === 'int') {
            const numeric = typeof value === 'number' ? value : Number(value);
            return Number.isFinite(numeric) ? Math.trunc(numeric).toLocaleString() : String(value);
        }
        if (columnType === 'float') {
            const numeric = typeof value === 'number' ? value : Number(value);
            if (!Number.isFinite(numeric)) {
                return String(value);
            }
            const absValue = Math.abs(numeric);
            if (absValue >= 10000 || (absValue > 0 && absValue < 0.001)) {
                return numeric.toExponential(3);
            }
            return new Intl.NumberFormat(undefined, {maximumFractionDigits: 4}).format(numeric);
        }
        if (columnType === 'bool' || columnType === 'boolean') {
            return value ? 'True' : 'False';
        }
        if (Array.isArray(value)) {
            return value.join(', ');
        }
        return String(value);
    }

    function metadataCellClass(columnKey, value) {
        if (isMissingMetadataValue(value)) {
            return 'metadata-cell-null';
        }
        const columnType = metadataColumnType(columnKey);
        return columnType === 'int' || columnType === 'float' ? 'metadata-cell-number' : '';
    }

    function compareMetadataValues(leftValue, rightValue) {
        if (leftValue === rightValue) {
            return 0;
        }
        if (typeof leftValue === 'number' && typeof rightValue === 'number') {
            return leftValue - rightValue;
        }
        const leftText = String(leftValue);
        const rightText = String(rightValue);
        const leftNumeric = Number(leftText);
        const rightNumeric = Number(rightText);
        if (Number.isFinite(leftNumeric) && Number.isFinite(rightNumeric) && leftText.trim() !== '' && rightText.trim() !== '') {
            return leftNumeric - rightNumeric;
        }
        return leftText.localeCompare(rightText, undefined, {numeric: true, sensitivity: 'base'});
    }

    function sortedMetadataNodeIndices(nodeIndices) {
        if (!state.metadataSort.columnKey) {
            return [...nodeIndices];
        }
        const direction = state.metadataSort.direction === 'desc' ? -1 : 1;
        const nullPlacement = metadataNullPlacement();
        return [...nodeIndices].sort((leftIndex, rightIndex) => {
            const leftValue = metadataSortValue(leftIndex, state.metadataSort.columnKey);
            const rightValue = metadataSortValue(rightIndex, state.metadataSort.columnKey);
            const leftMissing = isMissingMetadataValue(leftValue);
            const rightMissing = isMissingMetadataValue(rightValue);
            if (leftMissing || rightMissing) {
                if (leftMissing && rightMissing) {
                    return nodeId(leftIndex).localeCompare(nodeId(rightIndex), undefined, {numeric: true, sensitivity: 'base'});
                }
                return leftMissing
                    ? (nullPlacement === 'first' ? -1 : 1)
                    : (nullPlacement === 'first' ? 1 : -1);
            }
            const comparison = compareMetadataValues(
                leftValue,
                rightValue,
            );
            if (comparison !== 0) {
                return comparison * direction;
            }
            return nodeId(leftIndex).localeCompare(nodeId(rightIndex), undefined, {numeric: true, sensitivity: 'base'});
        });
    }

    window.addEventListener('pointermove', event => {
        if (!state.metadataResize) {
            return;
        }
        updateMetadataColumnWidth(
            state.metadataResize.columnKey,
            state.metadataResize.startWidth + (event.clientX - state.metadataResize.startX),
        );
    });

    window.addEventListener('pointerup', () => {
        if (!state.metadataResize) {
            return;
        }
        state.metadataResize = null;
        document.body.style.cursor = '';
        document.body.style.userSelect = '';
    });

    window.addEventListener('pointercancel', () => {
        if (!state.metadataResize) {
            return;
        }
        state.metadataResize = null;
        document.body.style.cursor = '';
        document.body.style.userSelect = '';
    });

        document.getElementById('metadata-select-nodes').addEventListener('click', selectNodesFromMetadataRows);

    function toggleMetadataSort(columnKey) {
        if (state.metadataSort.columnKey === columnKey) {
            state.metadataSort.direction = state.metadataSort.direction === 'asc' ? 'desc' : 'asc';
        } else {
            state.metadataSort = {columnKey, direction: 'asc'};
        }
        resetMetadataPage();
        updateMetadataTable();
    }

    function resetMetadataSort() {
        state.metadataSort = {columnKey: null, direction: 'asc'};
        resetMetadataPage();
        updateMetadataTable();
    }

