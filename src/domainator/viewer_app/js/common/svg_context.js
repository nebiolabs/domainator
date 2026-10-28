    // A stand-in for a canvas 2D context that records what is drawn on it as SVG, so a view can
    // draw its screen, its PNG export and its SVG export with one function. It covers the part of
    // the canvas API those drawings use: paths of lines, rectangles and whole circles, filled or
    // stroked; rectangles; text; translate, scale, save and restore; dashes and globalAlpha.
    // Clipping is ignored, so a drawing must keep to its area itself.
    function createSvgContext(width, height) {
        const parts = [];
        let path = [];
        const stack = [];
        const transform = {tx: 0, ty: 0, sx: 1, sy: 1};
        // A real context measures text, so labels are laid out as they are on the screen.
        const measurer = document.createElement('canvas').getContext('2d');
        const round = value => Math.round(value * 100) / 100;
        const X = x => round(transform.tx + x * transform.sx);
        const Y = y => round(transform.ty + y * transform.sy);
        const W = w => round(w * transform.sx);
        const H = h => round(h * transform.sy);

        function paint(style, kind) {
            const {color, opacity} = svgColorParts(style);
            const alpha = opacity * ctx.globalAlpha;
            return kind + '="' + escapeXml(color) + '"' + (alpha < 1 ? ' ' + kind + '-opacity="' + round(alpha) + '"' : '');
        }

        function strokeAttributes() {
            let attributes = paint(ctx.strokeStyle, 'stroke') + ' stroke-width="' + round(ctx.lineWidth * transform.sx) + '"';
            if (ctx.dash.length > 0) {
                attributes += ' stroke-dasharray="' + ctx.dash.map(round).join(' ') + '"';
            }
            return attributes;
        }

        const ctx = {
            fillStyle: '#000000',
            strokeStyle: '#000000',
            lineWidth: 1,
            font: '10px sans-serif',
            textAlign: 'start',
            textBaseline: 'alphabetic',
            globalAlpha: 1,
            lineCap: 'butt',
            lineJoin: 'miter',
            dash: [],
            save() {
                stack.push({
                    fillStyle: ctx.fillStyle, strokeStyle: ctx.strokeStyle, lineWidth: ctx.lineWidth, font: ctx.font,
                    textAlign: ctx.textAlign, textBaseline: ctx.textBaseline, globalAlpha: ctx.globalAlpha,
                    dash: ctx.dash.slice(), transform: Object.assign({}, transform),
                });
            },
            restore() {
                const saved = stack.pop();
                if (!saved) { return; }
                Object.assign(transform, saved.transform);
                delete saved.transform;
                Object.assign(ctx, saved);
            },
            translate(x, y) {
                transform.tx += x * transform.sx;
                transform.ty += y * transform.sy;
            },
            scale(x, y) {
                transform.sx *= x;
                transform.sy *= y;
            },
            setTransform() {},
            setLineDash(segments) { ctx.dash = segments.slice(); },
            clip() {},
            clearRect() {},
            beginPath() { path = []; },
            closePath() { path.push('Z'); },
            moveTo(x, y) { path.push('M' + X(x) + ' ' + Y(y)); },
            lineTo(x, y) { path.push('L' + X(x) + ' ' + Y(y)); },
            rect(x, y, w, h) {
                path.push('M' + X(x) + ' ' + Y(y) + 'h' + W(w) + 'v' + H(h) + 'h' + W(-w) + 'Z');
            },
            // Whole circles only, which is all the views draw.
            arc(x, y, r) {
                const radius = W(r);
                path.push('M' + round(X(x) - radius) + ' ' + Y(y) +
                    'a' + radius + ' ' + radius + ' 0 1 0 ' + round(2 * radius) + ' 0' +
                    'a' + radius + ' ' + radius + ' 0 1 0 ' + round(-2 * radius) + ' 0');
            },
            fill() {
                if (path.length > 0) {
                    parts.push('<path d="' + path.join('') + '" ' + paint(ctx.fillStyle, 'fill') + '/>');
                }
            },
            stroke() {
                if (path.length > 0) {
                    parts.push('<path d="' + path.join('') + '" fill="none" ' + strokeAttributes() + '/>');
                }
            },
            fillRect(x, y, w, h) {
                parts.push('<rect x="' + X(x) + '" y="' + Y(y) + '" width="' + W(w) + '" height="' + H(h) + '" ' +
                    paint(ctx.fillStyle, 'fill') + '/>');
            },
            strokeRect(x, y, w, h) {
                parts.push('<rect x="' + X(x) + '" y="' + Y(y) + '" width="' + W(w) + '" height="' + H(h) + '" fill="none" ' +
                    strokeAttributes() + '/>');
            },
            measureText(text) {
                measurer.font = ctx.font;
                return measurer.measureText(text);
            },
            fillText(text, x, y) {
                const anchor = {center: 'middle', right: 'end', end: 'end'}[ctx.textAlign] || 'start';
                const baseline = {middle: 'central', top: 'hanging', hanging: 'hanging', bottom: 'text-after-edge'}[ctx.textBaseline];
                parts.push('<text x="' + X(x) + '" y="' + Y(y) + '" style="font: ' + escapeXml(ctx.font) + '" ' +
                    paint(ctx.fillStyle, 'fill') + (anchor !== 'start' ? ' text-anchor="' + anchor + '"' : '') +
                    (baseline ? ' dominant-baseline="' + baseline + '"' : '') + '>' + escapeXml(String(text)) + '</text>');
            },
            toSVG() {
                return '<svg xmlns="http://www.w3.org/2000/svg" width="' + width + '" height="' + height +
                    '" viewBox="0 0 ' + width + ' ' + height + '">' + parts.join('') + '</svg>';
            },
        };
        return ctx;
    }
