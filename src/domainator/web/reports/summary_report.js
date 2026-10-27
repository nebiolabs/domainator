// The summary report (summary_report): contig counts and lengths, the domains found, how often
// other domains occur with the focus domains, and the taxonomy.
//
// report-data's data is {summary, length_histogram, taxonomy_tree}. summary is the document
// --json_out writes. length_histogram is null (no contigs) or {bin_edges, counts}, bin i
// counting the contig lengths in [bin_edges[i], bin_edges[i + 1]). taxonomy_tree is null
// (no --taxonomy) or a node {taxid, name, rank, count, children}, count being the contigs
// assigned to that taxon itself.

// @include common/payload.js
// @include common/download.js
// @include common/table.js
// @include common/charts.js
// @include common/hierarchy.js

function reportPanel(root, title, note = null) {
    const section = document.createElement("section");
    section.className = "panel";
    const heading = document.createElement("h2");
    heading.textContent = title;
    section.appendChild(heading);
    if (note) {
        const text = document.createElement("p");
        text.className = "note";
        text.textContent = note;
        section.appendChild(text);
    }
    root.appendChild(section);
    return section;
}

function panelPart(section, id) {
    const part = document.createElement("div");
    part.id = id;
    section.appendChild(part);
    return part;
}

// [label, value] pairs as a grid of figures.
function statGrid(entries) {
    const grid = document.createElement("dl");
    grid.className = "stat-grid";
    for (const [label, value] of entries) {
        const item = document.createElement("div");
        const term = document.createElement("dt");
        term.textContent = label;
        const figure = document.createElement("dd");
        figure.textContent = value;
        item.append(term, figure);
        grid.appendChild(item);
    }
    return grid;
}

function contigPanel(root, stats, histogram) {
    const section = reportPanel(root, "Contig Stats");
    const entries = [
        ["contigs", formatChartNumber(stats.contigs)],
        ["CDSs", formatChartNumber(stats.cdss)],
        ["partial CDSs", formatChartNumber(stats.partial_cdss)],
        ["protein fragments", formatChartNumber(stats.fragment_proteins)],
        ["CDSs per 10 kb", stats.cds_per_10kb.toFixed(2)],
    ];
    if (stats.length) {
        entries.push(
            ["shortest", formatChartNumber(stats.length.min)],
            ["median length", formatChartNumber(stats.length.median)],
            ["mean length", formatFixed1(stats.length.mean)],
            ["longest", formatChartNumber(stats.length.max)],
            ["total length", formatChartNumber(stats.length.total)],
        );
    }
    section.appendChild(statGrid(entries));
    const heading = document.createElement("h2");
    heading.textContent = "Contig Lengths";
    section.appendChild(heading);
    const chart = panelPart(section, "length-histogram");
    if (histogram) {
        drawHistogram(chart, {
            binEdges: histogram.bin_edges, counts: histogram.counts, xLabel: "contig length",
            yLabel: "contigs", itemLabel: "contigs", integerBins: true,
            fileStem: "summary_report_lengths", idPrefix: "lengths-",
        });
    } else {
        chart.innerHTML = '<div class="status">No contigs.</div>';
    }
}

function domainPanel(root, frequency) {
    const section = reportPanel(root, "Domain Frequency");
    createReportTable(panelPart(section, "domain-table"), {
        columns: [
            {name: "domain", type: "str"}, {name: "database", type: "str"}, {name: "count", type: "int"},
            {name: "description", type: "str"}, {name: "avg_score", type: "float"},
        ],
        rows: frequency,
        fileStem: "summary_report_domains",
        idPrefix: "domains-",
    });
}

// cooccurrence: {focus domain: {domain: percent}}, "" where the two never occur together.
function cooccurrencePanel(root, cooccurrence) {
    const section = reportPanel(root, "Domain Co-occurrence",
        "The percent of the contigs with the row's domain that also have the column's domain.");
    const names = [];
    const seen = new Set();
    for (const values of Object.values(cooccurrence)) {
        for (const name of Object.keys(values)) {
            if (!seen.has(name)) {
                seen.add(name);
                names.push(name);
            }
        }
    }
    const rows = Object.entries(cooccurrence).map(([focus, values]) => {
        const row = {domain: focus};
        for (const name of names) {
            const value = values[name];
            row[name] = value === "" || value === undefined ? null : value;
        }
        return row;
    });
    createReportTable(panelPart(section, "cooccurrence-table"), {
        columns: [{name: "domain", type: "str"}, ...names.map(name => ({name, type: "float"}))],
        rows,
        fileStem: "summary_report_cooccurrence",
        idPrefix: "cooccurrence-",
    });
}

function taxonomyPanel(root, tree, taxa) {
    const section = reportPanel(root, "Taxonomy",
        "Contigs by taxon. Click a taxon to zoom in to it, and click it again to zoom back out.");
    window.taxonomyChart = drawIcicle(panelPart(section, "taxonomy-chart"), tree, {
        itemLabel: "contigs", fileStem: "summary_report_taxonomy", idPrefix: "taxonomy-",
    });
    createReportTable(panelPart(section, "taxonomy-table"), {
        columns: [
            {name: "taxid", type: "int"}, {name: "name", type: "str"}, {name: "rank", type: "str"},
            {name: "count", type: "int"},
        ],
        rows: taxa || [],
        fileStem: "summary_report_taxonomy",
        idPrefix: "taxa-",
    });
}

(async () => {
    const root = document.getElementById("report-root");
    try {
        const {data} = await readReportData();
        const summary = data.summary;
        root.replaceChildren();
        contigPanel(root, summary.contig_stats, data.length_histogram);
        domainPanel(root, summary.domain_frequency);
        if (summary.domain_cooccurrence) {
            cooccurrencePanel(root, summary.domain_cooccurrence);
        }
        if (data.taxonomy_tree) {
            taxonomyPanel(root, data.taxonomy_tree, summary.taxonomy);
        }
    } catch (error) {
        showReportError(root, error);
    }
})();
