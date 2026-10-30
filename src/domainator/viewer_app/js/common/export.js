    function entireContentExport() {
        return document.getElementById('export-extent').value === 'entire';
    }

    function imageExportSuffix(extension, scaled = false) {
        const scale = selectedPngScale();
        return (entireContentExport() ? '_entire' : '_view') +
            (scaled && scale > 1 ? '@' + scale + 'x' : '') + '.' + extension;
    }

    function validateImageExportSize(width, height, raster = false) {
        if (![width, height].every(value => Number.isSafeInteger(value) && value > 0)) {
            throw new Error('Invalid image dimensions: ' + width + ' x ' + height + '.');
        }
        // Bound bitmap allocation before asking the browser for hundreds of MB.
        // Browsers with smaller canvas limits may still refuse encoding below this.
        if (raster && (width > 32767 || height > 32767 || width * height > 64000000)) {
            throw new Error('PNG is too large (' + width + ' x ' + height +
                ' px; maximum side 32767 px and area 64 million pixels). Use SVG or a lower PNG resolution.');
        }
    }

    async function runImageExport(action) {
        try {
            await action();
        } catch (error) {
            console.error(error);
            setStatus('Image export failed: ' + error.message);
        }
    }

    async function exportCanvasPNG(size, draw, suffix = imageExportSuffix('png')) {
        const scale = selectedPngScale();
        const width = Math.round(size.width * scale);
        const height = Math.round(size.height * scale);
        validateImageExportSize(width, height, true);
        const canvas = document.createElement('canvas');
        canvas.width = width;
        canvas.height = height;
        const ctx = canvas.getContext('2d');
        if (!ctx) { throw new Error('Could not allocate PNG canvas. Use SVG or a lower PNG resolution.'); }
        ctx.scale(scale, scale);
        draw(ctx);
        const blob = await new Promise((resolve, reject) => {
            canvas.toBlob(value => value ? resolve(value) :
                reject(new Error('Could not encode PNG (' + width + ' x ' + height +
                    ' px). Use SVG or a lower PNG resolution.')), 'image/png');
        });
        triggerDownload(URL.createObjectURL(blob), exportBaseName() + suffix, true);
    }

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
