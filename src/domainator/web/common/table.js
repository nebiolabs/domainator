// An interactive table: sort by a column, filter per column, page through the rows, drag
// columns to reorder them, drag a column's right edge to resize it (double-click the edge to
// fit its contents), right-click a cell to copy its value, and save the filtered rows as CSV,
// TSV or JSON.
//
// columns: [{name, type}], type "str", "int" or "float". rows: objects keyed by column name,
// null for a missing value. Floats are shown to one decimal place, as the TSV report writes
// them; the saved files keep the full values. idPrefix goes before the ids of the table's
// controls, so a page can hold more than one table.

const TABLE_PAGE_SIZES = [25, 50, 100, 1000];
const TABLE_COLUMN_WIDTH = 150;
const TABLE_MIN_COLUMN_WIDTH = 40;
const TABLE_MAX_FIT_WIDTH = 800;

// Python's "%0.1f": like toFixed(1), except that an exact tie rounds to the even digit
// (Python writes 1.25 as "1.2", toFixed as "1.3"). A double is exactly halfway between two
// tenths only when four times it is an odd integer (it ends in .25 or .75).
function formatFixed1(value) {
    if (!Number.isFinite(value)) {
        return Number.isNaN(value) ? "nan" : (value > 0 ? "inf" : "-inf");
    }
    if (Math.abs(value) >= 1e21) {
        return `${BigInt(value)}.0`;  // toFixed switches to exponent notation here
    }
    const quarters = value * 4;
    if (Number.isInteger(quarters) && quarters % 2 !== 0) {
        const tenths = Math.abs(value) * 10;  // exact: an odd multiple of 0.5
        const even = Math.floor(tenths) % 2 === 0 ? Math.floor(tenths) : Math.ceil(tenths);
        return `${value < 0 ? "-" : ""}${Math.floor(even / 10)}.${even % 10}`;
    }
    const text = value.toFixed(1);
    return Object.is(value, -0) ? "-0.0" : text;
}

// A float column's value in a saved CSV or TSV: every digit, and a whole number keeps its ".0",
// as Python writes a float.
function formatFullFloat(value) {
    if (value === null || value === undefined) {
        return value;
    }
    if (Object.is(value, -0)) {
        return "-0.0";
    }
    return Number.isInteger(value) && Math.abs(value) < 1e16 ? `${value}.0` : String(value);
}

function formatTableCell(value, type) {
    if (value === null || value === undefined) {
        return "";
    }
    if (type === "float") {
        return formatFixed1(Number(value));
    }
    return String(value);
}

// A number column's filter can be a comparison ("> 5", "<=0.1", "=3", "!= 0"); anything else,
// and every filter on a text column, matches a case-insensitive substring of the shown value.
function tableFilterPredicate(text, type) {
    const query = text.trim();
    if (query === "") {
        return null;
    }
    if (type !== "str") {
        const comparison = /^(>=|<=|!=|>|<|=)\s*(-?[0-9.]+(?:e[-+]?[0-9]+)?)$/i.exec(query);
        if (comparison) {
            const bound = Number(comparison[2]);
            const test = {
                ">": v => v > bound, ">=": v => v >= bound, "<": v => v < bound,
                "<=": v => v <= bound, "=": v => v === bound, "!=": v => v !== bound,
            }[comparison[1]];
            return value => value !== null && value !== undefined && test(Number(value));
        }
    }
    const needle = query.toLowerCase();
    return value => formatTableCell(value, type).toLowerCase().includes(needle);
}

function compareTableValues(a, b, type) {
    const aMissing = a === null || a === undefined;
    const bMissing = b === null || b === undefined;
    if (aMissing || bMissing) {
        return aMissing === bMissing ? 0 : (aMissing ? 1 : -1);
    }
    if (type === "str") {
        return String(a).localeCompare(String(b));
    }
    return Number(a) - Number(b);
}

function createReportTable(container, {columns, rows, maxHeight = null, fileStem = "data", idPrefix = ""}) {
    const state = {
        order: columns.map((_, index) => index),
        filters: columns.map(() => ""),
        sortColumn: null,
        sortDescending: false,
        pageSize: TABLE_PAGE_SIZES[0],
        page: 0,
        view: rows,
        widths: columns.map(() => TABLE_COLUMN_WIDTH),
    };

    const toolbar = document.createElement("div");
    toolbar.className = "toolbar";
    const downloads = [["CSV", "csv"], ["TSV", "tsv"], ["JSON", "json"]].map(([label, kind]) => {
        const button = document.createElement("button");
        button.type = "button";
        button.id = `${idPrefix}download-${kind}`;
        button.textContent = `Download ${label}`;
        button.addEventListener("click", () => download(kind));
        return button;
    });
    const spacer = document.createElement("div");
    spacer.className = "spacer";
    const pageSizeSelect = document.createElement("select");
    pageSizeSelect.setAttribute("aria-label", "Rows per page");
    for (const size of TABLE_PAGE_SIZES) {
        pageSizeSelect.add(new Option(`${size} rows`, String(size)));
    }
    const pagerButton = (id, label, text) => {
        const button = document.createElement("button");
        button.type = "button";
        button.id = idPrefix + id;
        button.textContent = text;
        button.setAttribute("aria-label", label);
        return button;
    };
    const first = pagerButton("page-first", "First page", "\u00ab First");
    const previous = pagerButton("page-prev", "Previous page", "\u2039 Prev");
    const next = pagerButton("page-next", "Next page", "Next \u203a");
    const last = pagerButton("page-last", "Last page", "Last \u00bb");
    // Type a page number and press Enter (or leave the box) to go to it: both fire "change".
    const pageInput = document.createElement("input");
    pageInput.type = "text";
    pageInput.inputMode = "numeric";
    pageInput.id = `${idPrefix}page-input`;
    pageInput.className = "page-input";
    pageInput.setAttribute("aria-label", "Page number");
    const pageCountLabel = document.createElement("span");
    pageCountLabel.className = "status";
    pageCountLabel.id = `${idPrefix}page-count`;
    const pager = document.createElement("span");
    pager.className = "pager";
    pager.append(first, previous, "Page", pageInput, pageCountLabel, next, last);
    const counter = document.createElement("span");
    counter.className = "status";
    counter.id = `${idPrefix}table-counter`;
    const copyStatus = document.createElement("span");
    copyStatus.className = "status";
    copyStatus.id = `${idPrefix}table-copy-status`;
    copyStatus.setAttribute("role", "status");
    toolbar.append(...downloads, copyStatus, spacer, counter, pager, pageSizeSelect);

    const wrap = document.createElement("div");
    wrap.className = "report-table-wrap";
    if (maxHeight !== null) {
        wrap.style.maxHeight = `${maxHeight}px`;
    }
    const table = document.createElement("table");
    table.className = "report-table";
    const colgroup = document.createElement("colgroup");
    table.appendChild(colgroup);
    const thead = table.createTHead();
    const tbody = table.createTBody();
    wrap.appendChild(table);
    container.replaceChildren(toolbar, wrap);

    pageSizeSelect.addEventListener("change", () => {
        state.pageSize = Number(pageSizeSelect.value);
        state.page = 0;
        renderBody();
    });
    // renderBody clamps the page, so "last" can ask for one past the end.
    first.addEventListener("click", () => { state.page = 0; renderBody(); });
    previous.addEventListener("click", () => { state.page -= 1; renderBody(); });
    next.addEventListener("click", () => { state.page += 1; renderBody(); });
    last.addEventListener("click", () => { state.page = Infinity; renderBody(); });
    pageInput.addEventListener("change", () => {
        const page = Number.parseInt(pageInput.value.trim(), 10);
        if (Number.isFinite(page)) {
            state.page = page - 1;
        }
        // An unreadable number puts the current page back in the box.
        renderBody();
    });

    function recompute() {
        const tests = state.filters
            .map((text, index) => [index, tableFilterPredicate(text, columns[index].type)])
            .filter(([, test]) => test !== null);
        let view = tests.length === 0
            ? rows.slice()
            : rows.filter(row => tests.every(([index, test]) => test(row[columns[index].name])));
        if (state.sortColumn !== null) {
            const {name, type} = columns[state.sortColumn];
            const sign = state.sortDescending ? -1 : 1;
            // Array.prototype.sort is stable, so ties keep the report's order.
            view.sort((a, b) => {
                const order = compareTableValues(a[name], b[name], type);
                // Missing values stay last in both directions.
                const aMissing = a[name] === null || a[name] === undefined;
                const bMissing = b[name] === null || b[name] === undefined;
                return aMissing || bMissing ? order : sign * order;
            });
        }
        state.view = view;
        state.page = 0;
    }

    function renderHead() {
        const row = document.createElement("tr");
        colgroup.replaceChildren(...state.order.map(index => {
            const col = document.createElement("col");
            col.dataset.column = String(index);
            return col;
        }));
        applyWidths();
        for (const index of state.order) {
            const column = columns[index];
            const th = document.createElement("th");
            th.draggable = true;
            th.dataset.column = String(index);
            th.title = column.name;
            const title = document.createElement("div");
            title.className = "column-title";
            const label = document.createElement("span");
            label.textContent = column.name;
            const mark = document.createElement("span");
            mark.className = "sort-mark";
            mark.textContent = state.sortColumn === index ? (state.sortDescending ? "▼" : "▲") : "";
            title.append(label, mark);
            title.addEventListener("click", () => {
                // Ascending, then descending, then the report's own order.
                if (state.sortColumn !== index) {
                    state.sortColumn = index;
                    state.sortDescending = false;
                } else if (!state.sortDescending) {
                    state.sortDescending = true;
                } else {
                    state.sortColumn = null;
                }
                recompute();
                renderHead();
                renderBody();
            });
            const filter = document.createElement("input");
            filter.type = "search";
            filter.placeholder = column.type === "str" ? "filter" : "filter, e.g. > 10";
            filter.setAttribute("aria-label", `Filter ${column.name}`);
            filter.value = state.filters[index];
            filter.addEventListener("input", () => {
                state.filters[index] = filter.value;
                recompute();
                renderBody();
            });
            filter.addEventListener("mousedown", event => event.stopPropagation());
            th.append(title, filter, resizeHandle(th, index, column.name));

            th.addEventListener("dragstart", event => {
                event.dataTransfer.setData("text/plain", String(index));
                event.dataTransfer.effectAllowed = "move";
            });
            th.addEventListener("dragover", event => {
                event.preventDefault();
                th.classList.add("drop-target");
            });
            th.addEventListener("dragleave", () => th.classList.remove("drop-target"));
            th.addEventListener("drop", event => {
                event.preventDefault();
                th.classList.remove("drop-target");
                moveColumn(Number(event.dataTransfer.getData("text/plain")), index);
            });
            row.appendChild(th);
        }
        thead.replaceChildren(row);
    }

    // The column widths are on <col> elements; with table-layout: fixed the table must be
    // exactly as wide as its columns, or the browser stretches them to fill it.
    function applyWidths() {
        let total = 0;
        for (const col of colgroup.children) {
            const width = state.widths[Number(col.dataset.column)];
            col.style.width = `${width}px`;
            total += width;
        }
        table.style.width = `${total}px`;
    }

    function setWidth(index, width) {
        state.widths[index] = Math.max(TABLE_MIN_COLUMN_WIDTH, Math.round(width));
        applyWidths();
    }

    // The width that shows the header and every value on the current page in full.
    function fitWidth(index) {
        const {name, type} = columns[index];
        const context = document.createElement("canvas").getContext("2d");
        const measure = (font, text) => {
            context.font = font;
            return context.measureText(text).width;
        };
        const th = thead.querySelector(`th[data-column="${index}"]`);
        const headStyle = getComputedStyle(th);
        const padding = parseFloat(headStyle.paddingLeft) + parseFloat(headStyle.paddingRight) + 2;
        // Room for the sort mark beside the name.
        let width = measure(headStyle.font, name) + measure(headStyle.font, "\u25b2") + 4;
        const cell = tbody.querySelector("td");
        if (cell) {
            const font = getComputedStyle(cell).font;
            const start = state.page * state.pageSize;
            for (const record of state.view.slice(start, start + state.pageSize)) {
                width = Math.max(width, measure(font, formatTableCell(record[name], type)));
            }
        }
        return Math.min(TABLE_MAX_FIT_WIDTH, Math.ceil(width + padding));
    }

    // A handle on a column's right edge: drag it to resize the column, double-click it to fit
    // the column to its contents. It takes focus, and the arrow keys resize too.
    function resizeHandle(th, index, name) {
        const handle = document.createElement("span");
        handle.className = "resize-handle";
        handle.tabIndex = 0;
        handle.setAttribute("role", "separator");
        handle.setAttribute("aria-orientation", "vertical");
        handle.setAttribute("aria-label", `Resize ${name}`);
        handle.title = "Drag to resize; double-click to fit";
        handle.addEventListener("pointerdown", event => {
            if (event.button !== 0) {
                return;
            }
            event.preventDefault();
            event.stopPropagation();
            // Otherwise the drag also starts moving the column.
            th.draggable = false;
            th.classList.add("resizing");
            document.body.classList.add("resizing-column");
            handle.setPointerCapture(event.pointerId);
            const startX = event.clientX;
            const startWidth = state.widths[index];
            const move = moveEvent => setWidth(index, startWidth + moveEvent.clientX - startX);
            const end = () => {
                handle.removeEventListener("pointermove", move);
                handle.removeEventListener("pointerup", end);
                handle.removeEventListener("pointercancel", end);
                th.draggable = true;
                th.classList.remove("resizing");
                document.body.classList.remove("resizing-column");
            };
            handle.addEventListener("pointermove", move);
            handle.addEventListener("pointerup", end);
            handle.addEventListener("pointercancel", end);
        });
        handle.addEventListener("dblclick", event => {
            event.stopPropagation();
            setWidth(index, fitWidth(index));
        });
        handle.addEventListener("keydown", event => {
            const step = event.shiftKey ? 50 : 10;
            if (event.key === "ArrowLeft" || event.key === "ArrowRight") {
                event.preventDefault();
                setWidth(index, state.widths[index] + (event.key === "ArrowRight" ? step : -step));
            }
        });
        handle.addEventListener("click", event => event.stopPropagation());
        return handle;
    }

    // Right-clicking a cell opens a menu with "Copy cell value", which copies the value itself:
    // a float with every digit, not rounded as it is shown. Shift + right-click opens the
    // browser's own menu.
    const menu = document.createElement("div");
    menu.className = "context-menu";
    menu.setAttribute("role", "menu");
    menu.hidden = true;
    const copyItem = document.createElement("button");
    copyItem.type = "button";
    copyItem.setAttribute("role", "menuitem");
    copyItem.textContent = "Copy cell value";
    menu.appendChild(copyItem);
    document.body.appendChild(menu);
    let menuCell = null;
    let copyStatusTimer = null;

    function hideMenu() {
        menu.hidden = true;
        menuCell = null;
    }

    function cellText(cell) {
        const value = state.view[cell.row][columns[cell.column].name];
        if (value === null || value === undefined) {
            return "";
        }
        return columns[cell.column].type === "float" ? formatFullFloat(value) : String(value);
    }

    tbody.addEventListener("contextmenu", event => {
        const td = event.target.closest("td");
        if (!td || event.shiftKey) {
            return;
        }
        event.preventDefault();
        menuCell = {row: Number(td.dataset.row), column: Number(td.dataset.column)};
        menu.hidden = false;
        // Kept inside the window.
        const x = Math.min(event.clientX, window.innerWidth - menu.offsetWidth - 4);
        const y = Math.min(event.clientY, window.innerHeight - menu.offsetHeight - 4);
        menu.style.left = `${Math.max(0, x)}px`;
        menu.style.top = `${Math.max(0, y)}px`;
        copyItem.focus();
    });
    copyItem.addEventListener("click", async () => {
        if (menuCell === null) {
            return;
        }
        const text = cellText(menuCell);
        hideMenu();
        const copied = await copyText(text);
        copyStatus.textContent = copied ? "Copied." : "Could not copy to the clipboard.";
        clearTimeout(copyStatusTimer);
        copyStatusTimer = setTimeout(() => { copyStatus.textContent = ""; }, 2000);
    });
    document.addEventListener("pointerdown", event => {
        if (!menu.hidden && !menu.contains(event.target)) {
            hideMenu();
        }
    });
    document.addEventListener("keydown", event => {
        if (event.key === "Escape" && !menu.hidden) {
            hideMenu();
        }
    });
    // Not closed on scroll: a scroll event can arrive just after the right-click that opened it.
    window.addEventListener("blur", hideMenu);
    window.addEventListener("resize", hideMenu);

    // Puts column `from` just before column `to`.
    function moveColumn(from, to) {
        if (from === to || !Number.isInteger(from)) {
            return;
        }
        const order = state.order.filter(index => index !== from);
        order.splice(order.indexOf(to), 0, from);
        state.order = order;
        renderHead();
        renderBody();
    }

    function renderBody() {
        const total = state.view.length;
        const pageCount = Math.max(1, Math.ceil(total / state.pageSize));
        state.page = Math.min(Math.max(state.page, 0), pageCount - 1);
        const start = state.page * state.pageSize;
        const end = Math.min(start + state.pageSize, total);
        const fragment = document.createDocumentFragment();
        for (let i = start; i < end; i += 1) {
            const record = state.view[i];
            const tr = document.createElement("tr");
            for (const index of state.order) {
                const {name, type} = columns[index];
                const td = document.createElement("td");
                td.dataset.row = String(i);
                td.dataset.column = String(index);
                const text = formatTableCell(record[name], type);
                td.textContent = text;
                td.title = text;
                if (type !== "str") {
                    td.className = "number";
                }
                tr.appendChild(td);
            }
            fragment.appendChild(tr);
        }
        tbody.replaceChildren(fragment);
        const filtered = total === rows.length ? "" : ` (filtered from ${rows.length})`;
        counter.textContent = total === 0
            ? `No rows${filtered}`
            : `Showing ${start + 1}-${end} of ${total} rows${filtered}`;
        first.disabled = previous.disabled = state.page === 0;
        next.disabled = last.disabled = state.page >= pageCount - 1;
        pageInput.value = String(state.page + 1);
        pageCountLabel.textContent = `of ${pageCount}`;
    }

    // Saves the rows the filters keep, in the shown order, with the columns as arranged.
    function download(kind) {
        const shown = state.order.map(index => columns[index].name);
        if (kind === "json") {
            const records = state.view.map(record => Object.fromEntries(shown.map(name => [name, record[name] ?? null])));
            downloadText(`${fileStem}.json`, "application/json", JSON.stringify(records, null, 1) + "\n");
            return;
        }
        const delimiter = kind === "tsv" ? "\t" : ",";
        const types = state.order.map(index => columns[index].type);
        const cells = state.view.map(record => shown.map((name, i) => types[i] === "float" ? formatFullFloat(record[name]) : record[name]));
        downloadText(`${fileStem}.${kind}`, kind === "tsv" ? "text/tab-separated-values" : "text/csv", delimitedText(shown, cells, delimiter));
    }

    recompute();
    renderHead();
    renderBody();
    return {
        filter(columnName, text) {
            const index = columns.findIndex(column => column.name === columnName);
            state.filters[index] = text;
            recompute();
            renderHead();
            renderBody();
        },
        get viewRows() { return state.view; },
    };
}
