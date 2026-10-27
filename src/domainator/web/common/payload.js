// Reading the data a report was built with.
//
// report_html.py puts one envelope in <script type="application/json" id="report-data">:
// {"report": <name>, "payload_version": <n>, "encoding": "json" | "gzip+base64", "data": ...}.
// With "gzip+base64", data is a base64 string of gzipped JSON. A table report also has a
// <script type="application/x-ndjson" id="report-rows"> block, one JSON object per line.

async function gunzipBase64Json(base64Text) {
    if (typeof DecompressionStream === "undefined") {
        throw new Error("This browser cannot read the report's compressed data (DecompressionStream is unavailable).");
    }
    const binary = atob(base64Text);
    const bytes = new Uint8Array(binary.length);
    for (let i = 0; i < binary.length; i += 1) {
        bytes[i] = binary.charCodeAt(i);
    }
    const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"));
    return JSON.parse(await new Response(stream).text());
}

async function readReportData() {
    const element = document.getElementById("report-data");
    const envelope = JSON.parse(element.textContent);
    if (envelope.encoding === "gzip+base64") {
        envelope.data = await gunzipBase64Json(envelope.data);
    } else if (envelope.encoding !== "json") {
        throw new Error(`Unknown report data encoding: ${envelope.encoding}`);
    }
    return envelope;
}

function readNdjsonRows(elementId) {
    const element = document.getElementById(elementId);
    const rows = [];
    for (const line of element.textContent.split("\n")) {
        if (line.trim() !== "") {
            rows.push(JSON.parse(line));
        }
    }
    return rows;
}

function showReportError(container, error) {
    const message = document.createElement("div");
    message.className = "error";
    message.textContent = `This report could not be shown: ${error && error.message ? error.message : error}`;
    container.replaceChildren(message);
}
