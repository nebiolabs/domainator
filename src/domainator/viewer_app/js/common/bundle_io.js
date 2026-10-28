    // The shell viewer is built before any bundle exists, and the file picker can
    // swap bundles at any time, so the heading and tab follow the loaded bundle.
    function updateViewerTitle(networkName) {
        const name = (networkName || '').trim();
        document.getElementById('viewer-title').textContent =
            name === '' ? VIEWER_APP_NAME : VIEWER_APP_NAME + ': ' + name;
        if (name !== '') {
            document.title = name;
        }
    }

    function setStatus(message) {
        document.getElementById('bundle-status').textContent = message;
    }

    async function writeTextToClipboard(text) {
        if (navigator.clipboard?.writeText) {
            await navigator.clipboard.writeText(text);
            return;
        }
        const helper = document.createElement('textarea');
        helper.value = text;
        helper.setAttribute('readonly', 'readonly');
        helper.style.position = 'fixed';
        helper.style.opacity = '0';
        document.body.appendChild(helper);
        helper.select();
        document.execCommand('copy');
        document.body.removeChild(helper);
    }


    function browserSupportsBundleLoading() {
        return 'DecompressionStream' in window;
    }

    function warnIfUnsupported() {
        if (!browserSupportsBundleLoading()) {
            document.getElementById('browser-warning').textContent = EMBEDDED_BUNDLE_BASE64
                ? 'This browser lacks DecompressionStream support; bundled data will still open, but loading external gzipped bundles is unavailable.'
                : 'This browser lacks DecompressionStream support; gzipped bundles cannot be loaded.';
        }
    }

    async function decodeBundleFile(file) {
        if (browserSupportsBundleLoading()) {
            try {
                const stream = file.stream().pipeThrough(new DecompressionStream('gzip'));
                const text = await new Response(stream).text();
                return JSON.parse(text);
            } catch (error) {
                console.warn('gzip decode failed, trying plain JSON', error);
            }
        }
        const text = await file.text();
        return JSON.parse(text);
    }

    function base64ToBytes(base64Value) {
        return Uint8Array.from(atob(base64Value), c => c.charCodeAt(0));
    }

    async function autoloadEmbeddedBundle() {
        if (!EMBEDDED_BUNDLE_BASE64) {
            return;
        }
        setStatus('Loading bundled data...');
        try {
            const bytes = base64ToBytes(EMBEDDED_BUNDLE_BASE64);
            let text;
            if (browserSupportsBundleLoading()) {
                const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));
                text = await new Response(stream).text();
            } else {
                text = new TextDecoder().decode(bytes);
            }
            installBundle(JSON.parse(text));
        } catch (error) {
            console.error(error);
            setStatus('Failed to load bundled data: ' + error.message);
        }
    }

    function formatValue(value) {
        if (value === null || value === undefined || value === '') {
            return '—';
        }
        if (typeof value === 'number') {
            return Number.isInteger(value) ? value.toLocaleString() : value.toFixed(2);
        }
        return String(value);
    }

    function htmlEscape(value) {
        return String(value).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
    }

    function exportBaseName() {
        return ((state.bundle && state.bundle.name) || activeView().bundleNoun).replace(/\s+/g, '_');
    }

    function triggerDownload(href, filename, revoke) {
        const link = document.createElement('a');
        link.href = href;
        link.download = filename;
        link.click();
        if (revoke) {
            URL.revokeObjectURL(href);
        }
    }

    function selectedPngScale() {
        const select = document.getElementById('export-png-scale');
        const value = select ? Number(select.value) : 1;
        return Number.isFinite(value) && value >= 1 ? value : 1;
    }

