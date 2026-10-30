    const NODE_SHAPES = [
        ['circle', 'Circle'], ['square', 'Square'], ['diamond', 'Diamond'],
        ['triangle_up', 'Triangle up'], ['triangle_down', 'Triangle down'], ['hexagon', 'Hexagon'],
    ];
    const SHAPE_PAGE_SIZE = 50;

    function shapeVertices(shape, x, y, radius) {
        if (shape === 'circle') { return null; }
        if (shape === 'cell') {
            return [[x - radius, y - radius], [x + radius, y - radius],
                [x + radius, y + radius], [x - radius, y + radius]];
        }
        const sides = shape === 'hexagon' ? 6 : shape.startsWith('triangle') ? 3 : 4;
        const angle = shape === 'square' ? Math.PI / 4
            : shape === 'triangle_down' ? Math.PI / 2 : -Math.PI / 2;
        return Array.from({length: sides}, (_, i) => [
            x + radius * Math.cos(angle + i * 2 * Math.PI / sides),
            y + radius * Math.sin(angle + i * 2 * Math.PI / sides),
        ]);
    }

    function appendNodeGlyph(ctx, shape, x, y, radius) {
        const vertices = shapeVertices(shape, x, y, radius);
        if (!vertices) {
            ctx.moveTo(x + radius, y);
            ctx.arc(x, y, radius, 0, 2 * Math.PI);
            return;
        }
        ctx.moveTo(...vertices[0]);
        vertices.slice(1).forEach(point => ctx.lineTo(...point));
        ctx.closePath();
    }

    function nodeGlyphSVG(shape, x, y, radius, fmt = String) {
        if (shape === 'circle') {
            return '<circle cx="' + fmt(x) + '" cy="' + fmt(y) + '" r="' + fmt(radius) + '"/>';
        }
        if (shape === 'cell') {
            return '<rect x="' + fmt(x - radius) + '" y="' + fmt(y - radius) +
                '" width="' + fmt(radius * 2) + '" height="' + fmt(radius * 2) + '"/>';
        }
        return '<polygon points="' + shapeVertices(shape, x, y, radius)
            .map(point => point.map(fmt).join(',')).join(' ') + '"/>';
    }

    function automaticShapes(columnName) {
        const values = {};
        sortPaletteKeys(distinctColumnValues(columnName, Infinity).values.map(entry => entry.key))
            .forEach((key, index) => {
                Object.defineProperty(values, key, {
                    value: NODE_SHAPES[index % NODE_SHAPES.length][0], enumerable: true,
                });
            });
        return values;
    }

    function shapeForNode(style, nodeIndex, fallback = 'circle') {
        if (!style || !style.shape_by) { return fallback; }
        const columnName = style.shape_by;
        if (!state.shapeCaches.has(columnName)) {
            const automatic = automaticShapes(columnName);
            const custom = Object.hasOwn(state.customShapes, columnName) ? state.customShapes[columnName] : null;
            state.shapeCaches.set(columnName, state.allNodeIndices.map(index => {
                const raw = metadataValue(index, columnName);
                if (raw === null || raw === undefined || raw === '') {
                    return custom?.null_shape || 'circle';
                }
                const key = String(raw);
                return (custom && Object.hasOwn(custom.values, key) ? custom.values[key] : automatic[key]);
            }));
        }
        return state.shapeCaches.get(columnName)[nodeIndex];
    }

    function shapePickerIsOpen() {
        return !document.getElementById('shape-picker-overlay').hidden;
    }

    function closeShapePicker() {
        document.getElementById('shape-picker-overlay').hidden = true;
        state.shapePickerColumn = null;
        state.shapePickerStyleId = null;
    }

    function openShapePicker(columnName, styleId) {
        if (!metadataColumn(columnName)) {
            setStatus('Choose a "Shape by" property first.');
            return;
        }
        state.shapePickerColumn = columnName;
        state.shapePickerStyleId = styleId;
        state.shapePickerPage = 0;
        document.getElementById('shape-picker-overlay').hidden = false;
        renderShapePicker();
        document.getElementById('shape-picker-close').focus();
    }

    function renderShapePicker() {
        const columnName = state.shapePickerColumn;
        if (!metadataColumn(columnName)) { closeShapePicker(); return; }
        const distinct = distinctColumnValues(columnName, Infinity);
        const keys = sortPaletteKeys(distinct.values.map(entry => entry.key));
        const automatic = automaticShapes(columnName);
        const custom = Object.hasOwn(state.customShapes, columnName) ? state.customShapes[columnName] : null;
        const pages = Math.max(1, Math.ceil(keys.length / SHAPE_PAGE_SIZE));
        state.shapePickerPage = Math.max(0, Math.min(pages - 1, state.shapePickerPage));
        document.getElementById('shape-picker-title').textContent = 'Shapes: ' + columnName;
        document.getElementById('shape-picker-note').textContent =
            keys.length + ' categories. ' + (keys.length > NODE_SHAPES.length
                ? 'Shapes repeat: there are ' + NODE_SHAPES.length + ' available shapes. ' : '') +
            'Page ' + (state.shapePickerPage + 1) + ' of ' + pages + '.';
        const list = document.getElementById('shape-picker-list');
        list.replaceChildren();
        function row(key, label, value, isNull) {
            const container = document.createElement('label');
            container.className = 'shape-assignment';
            const text = document.createElement('span');
            text.textContent = label;
            const preview = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
            preview.setAttribute('viewBox', '0 0 24 24');
            preview.setAttribute('aria-hidden', 'true');
            preview.innerHTML = nodeGlyphSVG(value, 12, 12, 9);
            const select = document.createElement('select');
            select.setAttribute('aria-label', 'Shape for ' + label);
            NODE_SHAPES.forEach(([name, title]) => select.add(new Option(title, name)));
            select.value = value;
            select.addEventListener('change', () => {
                if (!Object.hasOwn(state.customShapes, columnName)) {
                    Object.defineProperty(state.customShapes, columnName, {
                        value: {values: {}, null_shape: 'circle'}, writable: true, enumerable: true, configurable: true,
                    });
                }
                const mapping = state.customShapes[columnName];
                if (isNull) { mapping.null_shape = select.value; }
                else {
                    Object.defineProperty(mapping.values, key, {
                        value: select.value, writable: true, enumerable: true, configurable: true,
                    });
                }
                rebuildNodeColorCache();
                repaintNow();
                renderShapePicker();
            });
            container.append(text, preview, select);
            list.append(container);
        }
        keys.slice(state.shapePickerPage * SHAPE_PAGE_SIZE, (state.shapePickerPage + 1) * SHAPE_PAGE_SIZE)
            .forEach(key => row(key, key, custom && Object.hasOwn(custom.values, key) ? custom.values[key] : automatic[key], false));
        row('', 'No value', custom?.null_shape || 'circle', true);
        document.getElementById('shape-picker-prev').disabled = state.shapePickerPage === 0;
        document.getElementById('shape-picker-next').disabled = state.shapePickerPage === pages - 1;
    }

    function validateShapeMappings(value) {
        if (!value || typeof value !== 'object' || Array.isArray(value)) { throw new Error('not shape mappings'); }
        const valid = name => NODE_SHAPES.some(([id]) => id === name);
        Object.entries(value).forEach(([column, mapping]) => {
            if (!metadataColumn(column) || !mapping || !mapping.values ||
                    typeof mapping.values !== 'object' || Array.isArray(mapping.values) ||
                    !valid(mapping.null_shape) || !Object.values(mapping.values).every(valid)) {
                throw new Error('invalid shape mapping for "' + column + '"');
            }
        });
        return JSON.parse(JSON.stringify(value));
    }
