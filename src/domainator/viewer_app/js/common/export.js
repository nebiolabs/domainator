    function exportSelection() {
        const selected = Array.from(state.selectedNodeIndices).sort((left, right) => left - right);
        if (!state.bundle) {
            return;
        }
        // Columns the view adds (the SSN view's cluster numbers), in place of metadata
        // columns of the same names.
        const extraColumns = viewExportColumns();
        const extraNames = new Set(extraColumns.map(column => column.name));
        const metadataColumns = state.metadataColumns.filter(column => !extraNames.has(column.name));
        const exportedNodeIndices = metadataDisplayNodeIndices(selected.length > 0
            ? selected
            : viewerNodeIds().map((_, nodeIndex) => nodeIndex));
        const header = ['node_id', ...extraNames, ...metadataColumns.map(column => column.name)];
        const lines = [header.join('\t')];
        exportedNodeIndices.forEach(nodeIndex => {
            const row = [
                nodeId(nodeIndex),
                ...extraColumns.map(column => column.valueFor(nodeIndex) ?? ''),
                ...metadataColumns.map(column => {
                    const value = metadataValue(nodeIndex, column.name);
                    return value === null || value === undefined ? '' : String(value);
                }),
            ];
            lines.push(row.join('\t'));
        });
        const blob = new Blob([lines.join('\n')], {type: 'text/tab-separated-values'});
        const link = document.createElement('a');
        link.href = URL.createObjectURL(blob);
        link.download = (
            selected.length > 0
                ? ((state.bundle?.name || 'selection') + '_selected.tsv')
                : ((state.bundle?.name || activeView().bundleNoun) + '_table.tsv')
        ).replace(/\s+/g, '_');
        link.click();
        URL.revokeObjectURL(link.href);
    }

