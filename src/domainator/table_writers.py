"""Writers for the tabular reports (enum_report, hmmer_report): TSV, NDJSON and HTML."""
import json

from domainator import report_html


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


class HTMLTableWriter(JSONTableWriter):
    """An interactive HTML table (report_html.py's ``table_report`` page).

    The rows are the ones JSONTableWriter writes, streamed into the page; the page shows floats
    to one decimal place, as the TSV does, and its downloads keep the full values.
    """
    def __init__(self, columns, column_types, out_handle, max_height, title="Enum Report", report="enum_report"):
        super().__init__(columns, column_types, out_handle)
        self.max_height = max_height
        self.title = title
        self.report = report
        self._tail = None

    def write_header(self):
        data = {
            "columns": [{"name": column, "type": data_type} for column, data_type in zip(self.columns, self.column_types)],
            "max_height": self.max_height,
            "file_stem": self.report,
        }
        head, self._tail = report_html.report_parts("table_report", self.title, self.report, data)
        self.out_handle.write(head)

    def write_row(self, values):
        record = {
            column: self._coerce(value, self.column_types[i])
            for i, (column, value) in enumerate(zip(self.columns, values))
        }
        self.out_handle.write(report_html.report_row(record))

    def write_footer(self):
        self.out_handle.write(self._tail)
