    function worldToScreenPoint(x, y) {
        return {
            x: (x * state.viewTransform.scale) + state.viewTransform.offsetX,
            y: (y * state.viewTransform.scale) + state.viewTransform.offsetY,
        };
    }

    function screenToWorldPoint(x, y) {
        return {
            x: (x - state.viewTransform.offsetX) / Math.max(state.viewTransform.scale, 1e-9),
            y: (y - state.viewTransform.offsetY) / Math.max(state.viewTransform.scale, 1e-9),
        };
    }

    function canvasCoordinatesFromEvent(event) {
        const rect = viewerCanvas().getBoundingClientRect();
        const scaleX = viewerCanvas().width / rect.width;
        const scaleY = viewerCanvas().height / rect.height;
        return {
            x: (event.clientX - rect.left) * scaleX,
            y: (event.clientY - rect.top) * scaleY,
        };
    }

    // The zoom box shows the scale (canvas pixels per layout unit) and sets it, so a
    // zoom can be copied from one viewer into another. A typed zoom keeps the center of
    // the canvas in place, as the wheel keeps the point under the pointer.
    function formatZoomLevel(scale) {
        return String(Number(scale.toPrecision(4)));
    }

    function updateZoomLevelControl() {
        const input = document.getElementById('zoom-level');
        // Not while the user is typing in it: the force layout repaints on its own.
        if (document.activeElement === input) {
            return;
        }
        input.value = formatZoomLevel(state.viewTransform.scale);
    }

    function changeZoomLevel() {
        const input = document.getElementById('zoom-level');
        const scale = Number(input.value);
        // Number('') is 0, so a cleared box is refused here too.
        if (state.bundle && Number.isFinite(scale) && scale > 0) {
            const transform = state.viewTransform;
            const nextScale = Math.max(transform.minScale, Math.min(transform.maxScale, scale));
            const centerX = viewerCanvas().width / 2;
            const centerY = viewerCanvas().height / 2;
            const anchor = screenToWorldPoint(centerX, centerY);
            transform.scale = nextScale;
            transform.offsetX = centerX - (anchor.x * nextScale);
            transform.offsetY = centerY - (anchor.y * nextScale);
            repaintNow();
        }
        // Shows the clamped scale, or puts back the current one after bad input.
        input.value = formatZoomLevel(state.viewTransform.scale);
    }

    function normalizedSelectionBox(box) {
        return {
            left: Math.min(box.startX, box.endX),
            top: Math.min(box.startY, box.endY),
            width: Math.abs(box.endX - box.startX),
            height: Math.abs(box.endY - box.startY),
        };
    }

    function circleIntersectsRect(cx, cy, radius, rect) {
        const nearestX = Math.max(rect.left, Math.min(cx, rect.left + rect.width));
        const nearestY = Math.max(rect.top, Math.min(cy, rect.top + rect.height));
        const dx = cx - nearestX;
        const dy = cy - nearestY;
        return (dx * dx) + (dy * dy) <= radius * radius;
    }

    // Axis-aligned overlap between a screen-space rectangle [ax0,ay0]-[ax1,ay1] (corners in
    // any order) and a selection box {left, top, width, height}.
    function rectIntersectsRect(ax0, ay0, ax1, ay1, rect) {
        const left = Math.min(ax0, ax1);
        const right = Math.max(ax0, ax1);
        const top = Math.min(ay0, ay1);
        const bottom = Math.max(ay0, ay1);
        return left <= (rect.left + rect.width) && right >= rect.left &&
            top <= (rect.top + rect.height) && bottom >= rect.top;
    }

    // The tooltip a view shows beside the pointer.
    function showViewTip(event, text) {
        const tip = document.getElementById('view-tip');
        tip.textContent = text;
        tip.hidden = false;
        tip.style.left = (event.clientX + 14) + 'px';
        tip.style.top = (event.clientY + 14) + 'px';
    }

    function hideViewTip() {
        document.getElementById('view-tip').hidden = true;
    }

    // ---- Drawing shared by the tree and projection views ----

    const VIEW_INK = '#1e2a2f';
    const VIEW_MUTED = '#5c6a70';
    const VIEW_ACCENT = '#c8553d';
    const VIEW_PREVIEW_OUTLINE = 'rgba(90, 90, 90, 0.7)';
    const VIEW_FONT = 'system-ui, -apple-system, "Segoe UI", sans-serif';

    function niceLength(target) {
        const power = Math.pow(10, Math.floor(Math.log10(target)));
        return [5, 2, 1].map(step => step * power).find(value => value <= target) || power;
    }

    // A canvas sized to its box in device pixels, with a context drawing in CSS pixels; null
    // while the canvas is not shown.
    function viewCanvasContext(canvas) {
        const rect = canvas.getBoundingClientRect();
        if (rect.width < 1 || rect.height < 1) {
            return null;
        }
        const ratio = window.devicePixelRatio || 1;
        const width = Math.round(rect.width * ratio);
        const height = Math.round(rect.height * ratio);
        if (canvas.width !== width || canvas.height !== height) {
            canvas.width = width;
            canvas.height = height;
        }
        const ctx = canvas.getContext('2d');
        ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
        return {ctx, width: rect.width, height: rect.height};
    }
