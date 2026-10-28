    function rebuildNodeColorCache() {
        if (!state.bundle) {
            state.nodeColorCache = [];
            return;
        }
        const nodeCount = viewerNodeCount();
        const cache = new Array(nodeCount);
        for (let i = 0; i < nodeCount; i++) {
            cache[i] = nodeColor(i);
        }
        state.nodeColorCache = cache;
        // A palette edit that recolors the nodes should recolor an open chart
        // of the same column too.
        refreshColumnChartIfOpen();
    }

    // Static menu of the named palettes. The disabled first entry is what a
    // hand-edited column shows; "Reset to defaults" is the way back to
    // DEFAULT_CATEGORICAL_PALETTE, so the menu carries no pseudo-entry for it.
    function populateColorPaletteMenu() {
        const select = document.getElementById('color-palette');
        select.innerHTML = '';
        const customOption = document.createElement('option');
        customOption.value = '';
        customOption.textContent = 'Custom colors';
        customOption.disabled = true;
        select.appendChild(customOption);
        CATEGORICAL_PALETTES.forEach(scheme => {
            const option = document.createElement('option');
            option.value = scheme.name;
            option.textContent = scheme.label;
            option.title = scheme.note;
            select.appendChild(option);
        });
    }

    // Normalize a color string to '#RRGGBB' (uppercase), mirroring Python's
    // normalize_color_hex. Returns null for anything that isn't a 6-digit hex.
    // Mirrors color_genbank.normalize_color_hex: an optional '#' and six hex digits,
    // canonicalized to '#RRGGBB'. Loosened by one case for hand-typed input -- CSS-style
    // three-digit shorthand, where each digit doubles. Returns null for anything else.
    function normalizeColorHex(text) {
        const trimmed = String(text).trim().toUpperCase();
        const short = /^#?([0-9A-F]{3})$/.exec(trimmed);
        if (short) {
            const digits = short[1];
            return '#' + digits.charAt(0) + digits.charAt(0)
                + digits.charAt(1) + digits.charAt(1)
                + digits.charAt(2) + digits.charAt(2);
        }
        const match = /^#?([0-9A-F]{6})$/.exec(trimmed);
        return match ? '#' + match[1] : null;
    }

    // Interpolate two '#rrggbb' colors channel-wise in RGB. t is clamped to [0, 1].
    function lerpHexColor(a, b, t) {
        const clampT = Math.max(0, Math.min(1, t));
        const parse = hex => {
            const value = parseInt(hex.slice(1), 16);
            return [(value >> 16) & 255, (value >> 8) & 255, value & 255];
        };
        const left = parse(a);
        const right = parse(b);
        const channel = index => Math.round(left[index] + ((right[index] - left[index]) * clampT));
        const toHex = n => n.toString(16).padStart(2, '0');
        return '#' + toHex(channel(0)) + toHex(channel(1)) + toHex(channel(2));
    }

    // Named qualitative palettes (from domainator.utils.NAMED_CATEGORICAL_PALETTES) and
    // the neutral gray get_palette gives values with no color of their own.
    const CATEGORICAL_PALETTES = VIEWER_CONFIG.categoricalPalettes;
    const PALETTE_NO_VALUE_COLOR = VIEWER_CONFIG.otherColor;

    function paletteSchemeByName(schemeName) {
        return CATEGORICAL_PALETTES.find(scheme => scheme.name === schemeName) || null;
    }

    // Python's str comparison is by code point; localeCompare is not, and the assignment
    // order has to match domainator.utils.sort_palette_values to reproduce its colors.
    function comparePaletteText(left, right) {
        if (left < right) {
            return -1;
        }
        return left > right ? 1 : 0;
    }

    // Mirror of domainator.utils.sort_palette_values: numeric-looking keys sort
    // numerically (ties broken by text), everything else sorts by text. Nulls never
    // reach here -- distinctColumnValues counts them separately.
    function sortPaletteKeys(keys) {
        const numeric = [];
        const other = [];
        keys.forEach(key => {
            const number = Number(key);
            if (key.trim() !== '' && Number.isFinite(number)) {
                numeric.push({number, key});
            } else {
                other.push(key);
            }
        });
        numeric.sort((a, b) => (a.number - b.number) || comparePaletteText(a.key, b.key));
        other.sort(comparePaletteText);
        return numeric.map(entry => entry.key).concat(other);
    }

    // Fill a column's palette from a named scheme: colors are handed out in
    // sort_palette_values order and cycle once the scheme runs out, exactly like
    // get_palette. Values keep their assigned color afterwards, so individual swatches,
    // the color-table TSV and the legend all stay editable/exportable. `offset` skips
    // that many of the scheme's colors first, like get_palette's offset, so a second
    // column can avoid the colors a first one with that many values used.
    function applyNamedPalette(columnName, schemeName, offset = 0) {
        const scheme = paletteSchemeByName(schemeName);
        if (!scheme) {
            return null;
        }
        const start = paletteOffsetValue(offset);
        const distinct = distinctColumnValues(columnName, Infinity);
        const colors = {};
        sortPaletteKeys(distinct.values.map(entry => entry.key)).forEach((key, index) => {
            colors[key] = scheme.colors[(start + index) % scheme.colors.length];
        });
        const palette = ensureCategoricalPalette(columnName);
        palette.colors = colors;
        palette.nullColor = PALETTE_NO_VALUE_COLOR;
        palette.scheme = scheme.name;
        palette.offset = start;
        return {scheme, assigned: distinct.values.length};
    }

    // A whole number of colors to skip; anything else (blank, negative, fractional,
    // or missing from a palette saved before offsets existed) is clamped into one.
    function paletteOffsetValue(value) {
        const number = Math.floor(Number(value));
        return Number.isFinite(number) && number > 0 ? number : 0;
    }

    // The palette a column has before anyone customizes it. Values used to be
    // colored by a hue hashed from their text -- stable, but arbitrary, and free to
    // put two adjacent categories on nearly the same color. They now get
    // get_palette()'s own 64 distinct colors in sort_palette_values order, so an
    // untouched viewer agrees with what domainator build_ssn would have drawn.
    const DEFAULT_CATEGORICAL_PALETTE = 'domainator';

    // Derived, not user data: never written into a session, and never reported as
    // "custom colors" by the picker. Keyed like state.customPalettes (by
    // paletteKey, so a numeric column flipped to discrete gets its own entry) and
    // dropped by rebuildMetadataCaches(), since the assignment depends on the whole
    // set of values in the column.
    function defaultCategoricalPalette(columnName) {
        const key = paletteKey(columnName);
        const cached = state.defaultPalettes.get(key);
        if (cached) {
            return cached;
        }
        const bundleTable = state.bundle.defaults.color_table;
        if (bundleTable && columnName === state.bundle.defaults.color_by) {
            const palette = bundleColorTablePalette(columnName, bundleTable);
            state.defaultPalettes.set(key, palette);
            return palette;
        }
        const scheme = paletteSchemeByName(DEFAULT_CATEGORICAL_PALETTE);
        const colors = {};
        if (scheme) {
            sortPaletteKeys(distinctColumnValues(columnName, Infinity).values.map(entry => entry.key))
                .forEach((valueKey, index) => {
                    colors[valueKey] = scheme.colors[index % scheme.colors.length];
                });
        }
        // nullColor is deliberately left unset: what color an empty cell gets is a
        // separate choice from the categorical assignment, and categoricalColor's
        // own fallback is already the color the rest of the viewer uses for one.
        const palette = {
            type: 'categorical',
            colors,
            nullColor: null,
            scheme: scheme ? scheme.name : null,
        };
        state.defaultPalettes.set(key, palette);
        return palette;
    }

    // The palette from build_ssn --color_table. Like build_ssn, values the table
    // leaves out are gray; the em-dash row is the no-value color, as in loadColorTableText.
    function bundleColorTablePalette(columnName, table) {
        const colors = {};
        distinctColumnValues(columnName, Infinity).values.forEach(entry => {
            colors[entry.key] = PALETTE_NO_VALUE_COLOR;
        });
        let nullColor = null;
        Object.keys(table).forEach(value => {
            if (value === '—') {
                nullColor = table[value];
            } else {
                colors[value] = table[value];
            }
        });
        return {type: 'categorical', colors, nullColor, scheme: null};
    }

    // `palette` is a categorical palette ({colors, nullColor}) -- normally the
    // column's own, from customPalette(), which now always supplies one for a
    // discrete column. The hue hashed from the value's text is the last resort for
    // a value the palette says nothing about: one typed in after a color table was
    // loaded, say, or after the palette was hand-edited.
    function categoricalColor(value, palette) {
        if (value === null || value === undefined || value === '') {
            return (palette && palette.nullColor) || '#b3a89d';
        }
        const text = String(value);
        if (palette && palette.colors && palette.colors[text]) {
            return palette.colors[text];
        }
        let hash = 0;
        for (let i = 0; i < text.length; i++) {
            hash = ((hash << 5) - hash) + text.charCodeAt(i);
            hash |= 0;
        }
        const hue = Math.abs(hash) % 360;
        return 'hsl(' + hue + ' 58% 54%)';
    }

    // `palette` is an optional numeric custom palette ({stops, nullColor}).
    // With low+high set, interpolate in RGB (low->mid->high split at 0.5 when mid is
    // present); otherwise fall back to the default cyan->orange hsl gradient.
    function numericColor(value, minValue, maxValue, palette) {
        if (value === null || value === undefined || Number.isNaN(value)) {
            return (palette && palette.nullColor) || '#b3a89d';
        }
        // A custom palette is a list of {value, color} stops, kept sorted by whoever
        // stored it. Interpolate between the pair that brackets the value and hold the
        // end colors flat outside the ends -- so two stops are a plain ramp, three
        // reproduce a low/mid/high midpoint, and more shape the ramp arbitrarily.
        const stops = palette && palette.stops;
        if (stops && stops.length >= 2) {
            if (value <= stops[0].value) {
                return stops[0].color;
            }
            const last = stops[stops.length - 1];
            if (value >= last.value) {
                return last.color;
            }
            for (let index = 1; index < stops.length; index++) {
                const upper = stops[index];
                if (value <= upper.value) {
                    const lower = stops[index - 1];
                    const span = upper.value - lower.value;
                    return span > 0
                        ? lerpHexColor(lower.color, upper.color, (value - lower.value) / span)
                        : upper.color;
                }
            }
            return last.color;
        }
        // No custom palette: the built-in hue ramp across the column's own range.
        const fraction = maxValue <= minValue ? 0.5 : (value - minValue) / (maxValue - minValue);
        const clamped = Math.max(0, Math.min(1, fraction));
        const hue = 200 - (160 * clamped);
        const light = 72 - (24 * clamped);
        return 'hsl(' + hue + ' 72% ' + light + '%)';
    }

    function colorInfo(columnName) {
        return state.metadataColorInfoByName.get(columnName) || null;
    }

    // True when `columnName` holds numbers, and so can be colored either as a gradient
    // or as discrete categories.
    function columnIsNumericType(columnName) {
        const info = colorInfo(columnName);
        return !!info && info.baseType === 'numeric';
    }

    function columnIsCategoricalNumeric(columnName) {
        return columnIsNumericType(columnName) && state.categoricalColumns.has(columnName);
    }

    // Switch a numeric column between gradient and discrete coloring. No-op for columns
    // that are not numeric (those are always categorical).
    function setColumnCategorical(columnName, categorical) {
        const info = colorInfo(columnName);
        if (!info || info.baseType !== 'numeric') {
            return;
        }
        if (categorical) {
            state.categoricalColumns.add(columnName);
        } else {
            state.categoricalColumns.delete(columnName);
        }
        info.type = categorical ? 'categorical' : 'numeric';
    }

    // Custom palettes are stored per column AND per coloring mode, so flipping a numeric
    // column to categorical and back does not discard the gradient the user set up.
    function paletteKey(columnName) {
        return columnIsCategoricalNumeric(columnName) ? columnName + '\u0000categorical' : columnName;
    }

    function customPalette(columnName) {
        const stored = state.customPalettes[paletteKey(columnName)];
        if (stored) {
            return stored;
        }
        // A discrete column with nothing stored is not uncolored -- it is on the
        // default named palette. Only gradient columns fall through to null, where
        // numericColor's own built-in ramp takes over.
        const info = colorInfo(columnName);
        return info && info.type === 'categorical' ? defaultCategoricalPalette(columnName) : null;
    }

    function nodeColor(nodeIndex) {
        const columnName = currentColorField();
        const info = colorInfo(columnName);
        const value = metadataValue(nodeIndex, columnName);
        if (!info) {
            return '#d88f3d';
        }
        const palette = customPalette(columnName);
        if (info.type === 'numeric') {
            return numericColor(value, info.min, info.max, palette);
        }
        return categoricalColor(value, palette);
    }

    function escapeXml(text) {
        return String(text)
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;');
    }

    // `css` as a canvas reads it: '#rrggbb', or rgba(...) when it has an alpha.
    // Black when the browser does not understand it.
    function canvasColor(css) {
        if (!state.colorParserContext) {
            state.colorParserContext = document.createElement('canvas').getContext('2d');
        }
        const context = state.colorParserContext;
        context.fillStyle = '#000000';
        context.fillStyle = css;
        return context.fillStyle;
    }

    // Normalize any CSS color (hex, rgb, rgba, hsl) into an SVG-friendly
    // {color, opacity} pair using the canvas's own color parser. SVG editors
    // such as Illustrator/Inkscape choke on rgba()/space-separated hsl(), so we
    // split out the alpha into a separate opacity attribute.
    function svgColorParts(css) {
        const normalized = canvasColor(css);
        const match = /^rgba\(\s*([\d.]+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)\s*\)$/.exec(normalized);
        if (match) {
            return {color: 'rgb(' + match[1] + ',' + match[2] + ',' + match[3] + ')', opacity: Number(match[4])};
        }
        return {color: normalized, opacity: 1};
    }

    // Normalize any CSS color (hex, rgb, hsl) to '#rrggbb' using the canvas parser.
    // <input type="color"> only accepts '#rrggbb', but the default schemes emit
    // hsl(...), so color inputs must be seeded through this. Alpha is dropped.
    function cssToHex(css) {
        const normalized = canvasColor(css);
        if (normalized.charAt(0) === '#') {
            return normalized.length === 7 ? normalized : normalized.slice(0, 7);
        }
        const match = /rgba?\(\s*([\d.]+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)/.exec(normalized);
        if (!match) {
            return '#000000';
        }
        const toHex = n => Math.round(Number(n)).toString(16).padStart(2, '0');
        return '#' + toHex(match[1]) + toHex(match[2]) + toHex(match[3]);
    }

