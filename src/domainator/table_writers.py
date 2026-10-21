"""Writers for the tabular reports (enum_report, hmmer_report): TSV, NDJSON and HTML."""
import json


class TableWriter:
    """Writes rows of a table with named, typed ("str", "int" or "float") columns to a text handle."""

    def __init__(self, columns, column_types, out_handle):
        self.columns = columns
        self.column_types = column_types
        self.out_handle = out_handle

    def write_header(self):
        pass

    def write_row(self, values):
        raise NotImplementedError

    def write_footer(self):
        pass


class TSVTableWriter(TableWriter):
    def write_header(self):
        print("\t".join(self.columns), file=self.out_handle)

    def write_row(self, values):
        out_vals = list()
        try:
            for i, v in enumerate(values):
                if v is None:
                    out_vals.append("")
                    continue
                data_type = self.column_types[i]
                if data_type == "str":
                    out_vals.append(v.replace("\t", " "))
                if data_type == "int":
                    out_vals.append(str(int(v)))
                if data_type == "float":
                    out_vals.append(f"%0.1f" % float(v))
        except TypeError as exc:
            exc.add_note(f"while writing column {self.columns[i]!r} of {self.columns} (types {self.column_types}) from values {values}")
            raise
        print("\t".join(out_vals), file=self.out_handle)


class JSONTableWriter(TableWriter):
    """Streaming NDJSON writer: one compact JSON object per record, keyed by column name.

    NDJSON keeps the writer streaming (no need to buffer every row) while staying
    trivially parseable by agents. Values are coerced by column type; missing
    values are emitted as null so the schema is stable across rows.
    """
    def _coerce(self, value, data_type):
        if value is None:
            return None
        if data_type == "int":
            return int(value)
        if data_type == "float":
            return float(value)
        return str(value)

    def write_row(self, values):
        record = {
            column: self._coerce(value, self.column_types[i])
            for i, (column, value) in enumerate(zip(self.columns, values))
        }
        print(json.dumps(record, separators=(",", ":")), file=self.out_handle)


class HTMLTableWriter(TableWriter):
    """An interactive (tabulator) HTML table."""
    sorters = {"str": "string", "int": "number", "float": "number"}
    def __init__(self, columns, column_types, out_handle, max_height):
        super().__init__(columns, column_types, out_handle)
        self.max_height = max_height
        self.id_counter = 0

    def _format_cell(self, value, data_type):
        if data_type == "str":
            return json.dumps("" if value is None else value)
        if value is None:
            return "null"
        if data_type == "int":
            return str(int(value))
        if data_type == "float":
            return f"%0.1f" % float(value)
        raise ValueError(f"Unsupported column type: {data_type}")
    
    def write_header(self):
        print("""<!doctype html><html lang='en'>
<head>
<meta charset="UTF-8" />
<title>Enum Report</title>
<link href="https://unpkg.com/tabulator-tables@5.4.3/dist/css/tabulator.min.css" rel="stylesheet">
<script type="text/javascript" src="https://unpkg.com/tabulator-tables@5.4.3/dist/js/tabulator.min.js"></script>
</head>

<body>
<div>
    <button id="download-csv">Download CSV</button>
    <button id="download-html">Download HTML</button>
    <button id="download-json">Download JSON</button>
</div>
<div id="enum_report_table"></div>
<script>
var tabledata = [""",
        file=self.out_handle)
        self.id_counter = 0
        

    def write_row(self, values):
        out_vals = list()
        for i, v in enumerate(values):
            data_type = self.column_types[i]
            out_vals.append(self._format_cell(v, data_type))
        self.out_handle.write("{id:" + str(self.id_counter) + ",")
        for i,c in enumerate(out_vals):
            self.out_handle.write(f"""{json.dumps(self.columns[i])}:{c}, """)
        self.out_handle.write("},\n")
        self.id_counter += 1

    def write_footer(self):
        print("""];
    var table = new Tabulator("#enum_report_table", {
    data:tabledata,
    pagination:"local",
    paginationSize:25,
    paginationSizeSelector:[25, 50, 100, 1000],
    movableColumns:true,
    paginationCounter:"rows",""", file=self.out_handle)
        if self.max_height is not None:
            print(f"    maxHeight:{self.max_height},", file=self.out_handle)
        print("    columns:[", file=self.out_handle)
        for i, column in enumerate(self.columns):
            print(json.dumps({"title":column, "field":column, "width":150, "sorter": self.sorters[self.column_types[i]], "tooltip":True}) + ",", file=self.out_handle)
        print("    ],", file=self.out_handle)
        print("""});

//trigger download of data.csv file
document.getElementById("download-csv").addEventListener("click", function(){
    table.download("csv", "data.csv");
});

//trigger download of data.json file
document.getElementById("download-json").addEventListener("click", function(){
    table.download("json", "data.json");
});

//trigger download of data.html file
document.getElementById("download-html").addEventListener("click", function(){
    table.download("html", "data.html", {style:true});
});

</script>
</body></html>""", file=self.out_handle)
