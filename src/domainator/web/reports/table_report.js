// The table reports (enum_report, hmmer_report): the rows in report-rows, one table.
// report-data's data is {columns: [{name, type}], max_height: <pixels> | null, file_stem}.

// @include common/payload.js
// @include common/download.js
// @include common/table.js

(async () => {
    const container = document.getElementById("report-table");
    try {
        const {data} = await readReportData();
        window.reportTable = createReportTable(container, {
            columns: data.columns,
            rows: readNdjsonRows("report-rows"),
            maxHeight: data.max_height,
            fileStem: data.file_stem,
        });
    } catch (error) {
        showReportError(container, error);
    }
})();
