// Saving data from a report: text files built in the page, handed to the browser as a download,
// and text copied to the clipboard.

function downloadText(fileName, mimeType, text) {
    downloadBlob(fileName, new Blob([text], {type: mimeType}));
}

function downloadBlob(fileName, blob) {
    const url = URL.createObjectURL(blob);
    const anchor = document.createElement("a");
    anchor.href = url;
    anchor.download = fileName;
    document.body.appendChild(anchor);
    anchor.click();
    anchor.remove();
    setTimeout(() => URL.revokeObjectURL(url), 0);
}

function csvField(value) {
    if (value === null || value === undefined) {
        return "";
    }
    const text = String(value);
    return /[",\r\n]/.test(text) ? `"${text.replace(/"/g, '""')}"` : text;
}

function tsvField(value) {
    if (value === null || value === undefined) {
        return "";
    }
    return String(value).replace(/[\t\r\n]/g, " ");
}

// rows: arrays of cell values, in the order of header.
function delimitedText(header, rows, delimiter) {
    const field = delimiter === "\t" ? tsvField : csvField;
    const lines = [header.map(field).join(delimiter)];
    for (const row of rows) {
        lines.push(row.map(field).join(delimiter));
    }
    return lines.join("\n") + "\n";
}

// Resolves to whether the text reached the clipboard. The Clipboard API needs a secure context
// (a file:// page is one in Chromium, not everywhere), so this falls back to a hidden textarea.
async function copyText(text) {
    if (navigator.clipboard && window.isSecureContext) {
        try {
            await navigator.clipboard.writeText(text);
            return true;
        } catch (error) {
            // Refused (no permission, or the page lost focus): try the fallback.
        }
    }
    const area = document.createElement("textarea");
    area.value = text;
    area.setAttribute("readonly", "");
    area.style.position = "fixed";
    area.style.opacity = "0";
    document.body.appendChild(area);
    area.select();
    let copied = false;
    try {
        copied = document.execCommand("copy");
    } catch (error) {
        copied = false;
    }
    area.remove();
    return copied;
}
