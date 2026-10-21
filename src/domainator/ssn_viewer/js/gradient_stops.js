/*
 * Editable gradient stops: two mandatory ends plus any number between.
 *
 * Included into viewer.js's script scope; see session_state.js.
 */

    // ---------------------------------------------------------------------
    // Gradient stops
    //
    // A numeric palette is an ordered list of {value, color} stops: two
    // mandatory ends plus any number of intermediates. Coloring interpolates
    // between the bracketing pair and holds the end colors flat outside the
    // ends, so two stops give a plain ramp and three reproduce the old
    // low/mid/high midpoint exactly.
    //
    // `state.gradientStops` is the working list the dialog edits; committing it
    // writes a sorted copy into the column's stored palette. The histogram, the
    // gradient bar and the slider all read that one list, so they cannot drift
    // apart, and every edit -- typing, dragging a knob, adding or removing a
    // stop -- funnels through commitGradientStops().
    // ---------------------------------------------------------------------

    const MIN_GRADIENT_STOPS = 2;
    // A cap, not a limitation of the model: past a dozen the rows stop fitting
    // in the dialog and the knobs start overlapping on the track.
    const MAX_GRADIENT_STOPS = 12;
    const DEFAULT_NO_VALUE_COLOR = '#b3a89d';

    function sortedGradientStops(stops) {
        return stops.slice().sort((left, right) => left.value - right.value);
    }

    function copyGradientStops(stops) {
        return stops.map(stop => ({value: stop.value, color: stop.color}));
    }

    // Read a stored palette as stops. Sessions written before stops existed
    // carry {low, mid, high} with separate lowValue/midValue/highValue bounds;
    // converting here means nothing downstream has to know about that shape.
    function numericPaletteStops(palette, minValue, maxValue) {
        if (!palette) {
            return null;
        }
        if (Array.isArray(palette.stops) && palette.stops.length >= MIN_GRADIENT_STOPS) {
            return sortedGradientStops(palette.stops
                .filter(stop => stop && Number.isFinite(Number(stop.value)) && stop.color)
                .map(stop => ({value: Number(stop.value), color: stop.color})));
        }
        if (!palette.low || !palette.high) {
            return null;
        }
        const low = (palette.lowValue !== null && palette.lowValue !== undefined)
            ? palette.lowValue : minValue;
        const high = (palette.highValue !== null && palette.highValue !== undefined)
            ? palette.highValue : maxValue;
        const stops = [{value: low, color: palette.low}, {value: high, color: palette.high}];
        if (palette.mid) {
            const mid = (palette.midValue !== null && palette.midValue !== undefined)
                ? palette.midValue : ((low + high) / 2);
            stops.push({value: mid, color: palette.mid});
        }
        return sortedGradientStops(stops);
    }

    // Rewrite any legacy numeric palettes in place, so a session saved by an
    // older viewer keeps its colors instead of silently reverting to defaults.
    function normalizeCustomPalettes(palettes) {
        Object.keys(palettes || {}).forEach(key => {
            const palette = palettes[key];
            if (!palette || palette.type !== 'numeric' || Array.isArray(palette.stops)) {
                return;
            }
            const stops = numericPaletteStops(palette, 0, 1);
            palettes[key] = {
                type: 'numeric',
                stops: stops || [],
                nullColor: palette.nullColor || null,
            };
        });
        return palettes;
    }

    function currentGradientStops() {
        return Array.isArray(state.gradientStops) ? state.gradientStops : [];
    }

    // ---------------------------------------------------------------------
    // The shared axis
    //
    // The histogram's binned extent, not the outermost stops: narrowing the
    // ramp then moves the knobs inward instead of rescaling the chart under
    // them.
    // ---------------------------------------------------------------------
    function gradientRangeDomain() {
        const histogram = state.colorHistogram;
        if (histogram && histogram.highEdge > histogram.lowEdge) {
            return {low: histogram.lowEdge, high: histogram.highEdge};
        }
        const info = colorInfo(currentColorField());
        if (info && Number.isFinite(info.min) && Number.isFinite(info.max) && info.max > info.min) {
            return {low: info.min, high: info.max};
        }
        return null;
    }

    // Where a value sits on the axis, as a percentage. Deliberately unclamped:
    // CSS gradient stops accept out-of-range percentages and extend the ramp
    // past the box, which is exactly how a stop outside the data range should
    // render. Knob placement clamps separately.
    function gradientRangePercent(value, domain) {
        return ((value - domain.low) / (domain.high - domain.low)) * 100;
    }

    function gradientValueAtPercent(percent, domain) {
        return domain.low + ((percent / 100) * (domain.high - domain.low));
    }

    // Integer columns snap to whole numbers; float columns round to a fixed
    // number of decimals chosen from the span, so dragging writes 3.42 rather
    // than 3.4166666666666665 into the number input.
    function roundGradientValue(value, domain) {
        const histogram = state.colorHistogram;
        if (histogram && histogram.integer) {
            return Math.round(value);
        }
        const span = domain.high - domain.low;
        const decimals = span > 0
            ? Math.max(0, Math.min(6, 3 - Math.floor(Math.log10(span))))
            : 2;
        return Number(value.toFixed(decimals));
    }

    function gradientKnobStep(domain) {
        const histogram = state.colorHistogram;
        if (histogram && histogram.integer) {
            return 1;
        }
        return (domain.high - domain.low) / 100;
    }

    // A stop is held between its neighbors so the ramp cannot fold over during
    // a drag. Typing into the number inputs is left unconstrained -- values are
    // re-sorted on commit instead, so a stop can be typed past its neighbor.
    function clampGradientStopValue(index, value, domain) {
        const stops = currentGradientStops();
        let lowLimit = domain.low;
        let highLimit = domain.high;
        if (index > 0) {
            lowLimit = Math.max(lowLimit, stops[index - 1].value);
        }
        if (index < stops.length - 1) {
            highLimit = Math.min(highLimit, stops[index + 1].value);
        }
        if (highLimit < lowLimit) {
            highLimit = lowLimit;
        }
        return Math.max(lowLimit, Math.min(highLimit, value));
    }

    // ---------------------------------------------------------------------
    // Committing an edit
    // ---------------------------------------------------------------------
    function commitGradientStops() {
        const columnName = currentColorField();
        if (!columnName) {
            return;
        }
        const stored = ensureNumericPalette(columnName);
        // Stored sorted so numericColor -- which runs once per node -- can walk
        // the stops without sorting them itself.
        stored.stops = copyGradientStops(sortedGradientStops(currentGradientStops()));
        stored.nullColor = document.getElementById('color-null').value;
        updateGradientPreview();
        drawColorHistogram(columnName);
        rebuildNodeColorCache();
        scheduleClusterRender();
    }

    function gradientStopRole(index, total) {
        if (index === 0) {
            return 'Min';
        }
        if (index === total - 1) {
            return 'Max';
        }
        return 'Stop ' + (index + 1);
    }

    function renderGradientStopRows() {
        const stops = currentGradientStops();
        const list = document.getElementById('color-stop-list');
        list.innerHTML = stops.map((stop, index) => {
            const role = gradientStopRole(index, stops.length);
            const removable = index > 0 && index < stops.length - 1;
            return '<div class="cp-stop-row" data-stop-index="' + index + '">'
                + '<span class="cp-stop-role">' + role + '</span>'
                + '<input type="color" data-stop-field="color" value="' + stop.color
                + '" aria-label="' + role + ' color" />'
                + '<input type="text" class="cp-stop-hex" data-stop-field="hex"'
                + ' value="' + normalizeColorHex(stop.color) + '" spellcheck="false"'
                + ' autocapitalize="off" autocomplete="off"'
                + ' aria-label="' + role + ' hex code" />'
                + '<input type="number" step="any" data-stop-field="value" value="' + stop.value
                + '" aria-label="' + role + ' value" />'
                + (removable
                    ? '<button type="button" class="cp-stop-remove" data-stop-remove'
                        + ' aria-label="Remove ' + role + '" title="Remove this stop">×</button>'
                    : '<span class="cp-stop-remove-gap"></span>')
                + '</div>';
        }).join('');
        const addButton = document.getElementById('color-add-stop');
        addButton.disabled = stops.length >= MAX_GRADIENT_STOPS;
        addButton.title = addButton.disabled
            ? ('At the limit of ' + MAX_GRADIENT_STOPS + ' stops.')
            : 'Add a stop in the widest gap, in the color the ramp already has there';
    }

    function addGradientStop() {
        const stops = currentGradientStops();
        const domain = gradientRangeDomain();
        if (stops.length < MIN_GRADIENT_STOPS || !domain) {
            return;
        }
        if (stops.length >= MAX_GRADIENT_STOPS) {
            setStatus('A gradient holds at most ' + MAX_GRADIENT_STOPS + ' stops.');
            return;
        }
        // Split the widest gap: with no other information that is where an extra
        // stop buys the most control.
        let gapIndex = 1;
        let widest = -Infinity;
        for (let index = 1; index < stops.length; index++) {
            const gap = stops[index].value - stops[index - 1].value;
            if (gap > widest) {
                widest = gap;
                gapIndex = index;
            }
        }
        const lower = stops[gapIndex - 1];
        const upper = stops[gapIndex];
        const midpoint = (lower.value + upper.value) / 2;
        const rounded = roundGradientValue(midpoint, domain);
        // Rounding an integer column would land the new stop on top of a
        // neighbor when the gap is a single unit; keep the exact midpoint there.
        const value = (rounded > lower.value && rounded < upper.value) ? rounded : midpoint;
        stops.splice(gapIndex, 0, {
            // Taking the color the ramp already has here means adding a stop
            // changes nothing until it is moved: you shape a ramp, not reset it.
            value,
            color: lerpHexColor(lower.color, upper.color, 0.5),
        });
        renderGradientStopRows();
        commitGradientStops();
        setStatus('Added a gradient stop at ' + formatValue(value) + '.');
    }

    function removeGradientStop(index) {
        const stops = currentGradientStops();
        if (index <= 0 || index >= stops.length - 1) {
            return;   // the two ends define the ramp and are not removable
        }
        const removed = stops.splice(index, 1)[0];
        renderGradientStopRows();
        commitGradientStops();
        setStatus('Removed the gradient stop at ' + formatValue(removed.value) + '.');
    }

    // "Reset values to data range": refit the ramp onto the data while keeping
    // the relative spacing of any stops the user has placed.
    function resetGradientStopValues() {
        const info = colorInfo(currentColorField());
        const stops = currentGradientStops();
        const domain = gradientRangeDomain();
        if (!info || stops.length < MIN_GRADIENT_STOPS) {
            return;
        }
        const first = stops[0].value;
        const span = stops[stops.length - 1].value - first;
        stops.forEach((stop, index) => {
            // Even spacing is the only sensible reading when every stop
            // currently sits on the same value.
            const fraction = span > 0
                ? (stop.value - first) / span
                : index / (stops.length - 1);
            const value = info.min + (fraction * (info.max - info.min));
            stop.value = domain ? roundGradientValue(value, domain) : value;
        });
        renderGradientStopRows();
        commitGradientStops();
    }

    function setGradientStopHexField(index, color) {
        const field = document.querySelector(
            '.cp-stop-row[data-stop-index="' + index + '"] [data-stop-field="hex"]');
        if (field) {
            field.value = normalizeColorHex(color) || color;
        }
    }

    function setupGradientStopEditing() {
        const list = document.getElementById('color-stop-list');
        const stopIndexOf = target => {
            const row = target.closest('[data-stop-index]');
            return row ? Number(row.dataset.stopIndex) : -1;
        };
        list.addEventListener('input', event => {
            const field = event.target.dataset ? event.target.dataset.stopField : null;
            const index = stopIndexOf(event.target);
            const stop = currentGradientStops()[index];
            if (!field || !stop) {
                return;
            }
            if (field === 'color') {
                stop.color = event.target.value;
                setGradientStopHexField(index, stop.color);
            } else if (field === 'hex') {
                const hex = normalizeColorHex(event.target.value);
                if (hex === null) {
                    return;   // mid-edit ("#", "#1a2"); the change handler has the say
                }
                // Stored lower case, which is the one form the color input and
                // lerpHexColor both produce; the field itself shows the canonical
                // upper-case spelling.
                stop.color = hex.toLowerCase();
                const swatch = document.querySelector(
                    '.cp-stop-row[data-stop-index="' + index + '"] [data-stop-field="color"]');
                if (swatch) { swatch.value = stop.color; }
            } else {
                const value = parseNumberOrNull(event.target.value);
                if (value === null) {
                    return;   // mid-edit ("", "-", "1e"); wait for something numeric
                }
                stop.value = value;
            }
            commitGradientStops();
        });
        // Leaving the hex field rewrites whatever was typed into the canonical
        // spelling; text that is not a color at all reverts to the stop's color
        // rather than being silently dropped or left sitting there looking valid.
        list.addEventListener('change', event => {
            if (!event.target.dataset || event.target.dataset.stopField !== 'hex') {
                return;
            }
            const index = stopIndexOf(event.target);
            const stop = currentGradientStops()[index];
            if (!stop) {
                return;
            }
            const typed = event.target.value.trim();
            // Clearing the field and tabbing away reads as canceling the edit, so
            // that reverts quietly; anything else that is not a color is worth saying.
            if (typed !== '' && normalizeColorHex(typed) === null) {
                setStatus('"' + typed + '" is not a hex color (expected #RGB or #RRGGBB).');
            }
            setGradientStopHexField(index, stop.color);
        });
        // Re-sorting mid-keystroke would yank the row out from under the cursor,
        // so a stop typed past its neighbor is only reordered on commit.
        list.addEventListener('change', event => {
            if (!event.target.dataset || event.target.dataset.stopField !== 'value') {
                return;
            }
            const stops = currentGradientStops();
            const sorted = sortedGradientStops(stops);
            if (sorted.some((stop, index) => stop !== stops[index])) {
                state.gradientStops = sorted;
                renderGradientStopRows();
                commitGradientStops();
            }
        });
        list.addEventListener('click', event => {
            const button = event.target.closest('[data-stop-remove]');
            if (button) {
                removeGradientStop(stopIndexOf(button));
            }
        });
        document.getElementById('color-add-stop').addEventListener('click', addGradientStop);
    }

    // ---------------------------------------------------------------------
    // The slider: one knob per stop
    // ---------------------------------------------------------------------
    function updateGradientSlider() {
        const slider = document.getElementById('color-range-slider');
        const domain = gradientRangeDomain();
        const stops = currentGradientStops();
        if (!domain || stops.length < MIN_GRADIENT_STOPS) {
            slider.hidden = true;
            return;
        }
        slider.hidden = false;
        let knobs = slider.querySelectorAll('.cp-range-knob');
        if (knobs.length !== stops.length) {
            // Rebuilt only when the count changes, so a drag is never
            // interrupted by its own knob being replaced underneath it.
            knobs.forEach(knob => knob.remove());
            stops.forEach((stop, index) => {
                const knob = document.createElement('div');
                knob.className = 'cp-range-knob';
                knob.dataset.stopIndex = String(index);
                knob.setAttribute('role', 'slider');
                knob.setAttribute('tabindex', '0');
                slider.appendChild(knob);
            });
            knobs = slider.querySelectorAll('.cp-range-knob');
        }
        knobs.forEach((knob, index) => {
            const stop = stops[index];
            const role = gradientStopRole(index, stops.length);
            const isEnd = index === 0 || index === stops.length - 1;
            knob.classList.toggle('cp-range-knob-mid', !isEnd);
            knob.style.left = Math.max(0, Math.min(100, gradientRangePercent(stop.value, domain))) + '%';
            knob.style.background = stop.color;
            knob.setAttribute('aria-label', role + ' value');
            knob.setAttribute('aria-valuemin', String(domain.low));
            knob.setAttribute('aria-valuemax', String(domain.high));
            knob.setAttribute('aria-valuenow', String(stop.value));
            knob.setAttribute('aria-valuetext', formatValue(stop.value));
        });
        // Shade the part of the axis the ramp actually spans.
        const span = document.getElementById('color-range-span');
        const lowPercent = Math.max(0, Math.min(100, gradientRangePercent(stops[0].value, domain)));
        const highPercent = Math.max(0, Math.min(100,
            gradientRangePercent(stops[stops.length - 1].value, domain)));
        span.hidden = false;
        span.style.left = Math.min(lowPercent, highPercent) + '%';
        span.style.width = Math.abs(highPercent - lowPercent) + '%';
    }

    // Move a stop from the slider and push the result through the same commit
    // path as typing, so there is one way for a gradient to change.
    function setGradientStopFromSlider(index, rawValue, domain) {
        const stop = currentGradientStops()[index];
        if (!stop) {
            return;
        }
        const value = roundGradientValue(clampGradientStopValue(index, rawValue, domain), domain);
        if (value === stop.value) {
            return;
        }
        stop.value = value;
        const input = document.querySelector(
            '.cp-stop-row[data-stop-index="' + index + '"] [data-stop-field="value"]');
        if (input) {
            input.value = value;
        }
        commitGradientStops();
    }

    function setupGradientRangeSlider() {
        const slider = document.getElementById('color-range-slider');
        slider.addEventListener('pointerdown', event => {
            const knob = event.target.closest('.cp-range-knob');
            const domain = gradientRangeDomain();
            if (!knob || !domain) {
                return;
            }
            const index = Number(knob.dataset.stopIndex);
            event.preventDefault();
            knob.focus();
            // Pointer capture keeps the drag alive past the knob's 18px, and
            // past the dialog edge, without a document-level move listener.
            knob.setPointerCapture(event.pointerId);
            const drag = moveEvent => {
                const rect = slider.getBoundingClientRect();
                if (rect.width <= 0) {
                    return;
                }
                const percent = ((moveEvent.clientX - rect.left) / rect.width) * 100;
                setGradientStopFromSlider(index, gradientValueAtPercent(percent, domain), domain);
            };
            const stop = () => {
                knob.removeEventListener('pointermove', drag);
                knob.removeEventListener('pointerup', stop);
                knob.removeEventListener('pointercancel', stop);
            };
            knob.addEventListener('pointermove', drag);
            knob.addEventListener('pointerup', stop);
            knob.addEventListener('pointercancel', stop);
            drag(event);
        });
        slider.addEventListener('keydown', event => {
            const knob = event.target.closest('.cp-range-knob');
            const domain = gradientRangeDomain();
            if (!knob || !domain) {
                return;
            }
            const index = Number(knob.dataset.stopIndex);
            const stop = currentGradientStops()[index];
            if (!stop) {
                return;
            }
            const step = gradientKnobStep(domain);
            let next = null;
            if (event.key === 'ArrowLeft' || event.key === 'ArrowDown') { next = stop.value - step; }
            else if (event.key === 'ArrowRight' || event.key === 'ArrowUp') { next = stop.value + step; }
            else if (event.key === 'PageDown') { next = stop.value - (step * 10); }
            else if (event.key === 'PageUp') { next = stop.value + (step * 10); }
            else if (event.key === 'Home') { next = domain.low; }
            else if (event.key === 'End') { next = domain.high; }
            else if (event.key === 'Delete' || event.key === 'Backspace') {
                event.preventDefault();
                removeGradientStop(index);
                return;
            }
            if (next === null) {
                return;
            }
            event.preventDefault();
            setGradientStopFromSlider(index, next, domain);
        });
    }
