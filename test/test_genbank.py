"""Tests for domainator.genbank, the GenBank layer on upstream Biopython."""
import glob
import io
import os

import pytest

from domainator.genbank import parse_genbank, write_genbank_records

_DATA = os.path.join(os.path.dirname(__file__), "data")
_UNPARSEABLE = {"simple_genpept_equals_second_line.gb"}  # deliberately malformed


def _signature(record):
    return (record.id, record.name, record.description, str(record.seq),
            sorted((key, str(value)) for key, value in record.annotations.items()),
            [(f.type, str(f.location), {k: list(v) for k, v in f.qualifiers.items()}) for f in record.features])


def _genbank_fixtures():
    paths = sorted(glob.glob(os.path.join(_DATA, "*.gb")) + glob.glob(os.path.join(_DATA, "*.gbk")))
    return [p for p in paths if os.path.getsize(p) > 0 and os.path.basename(p) not in _UNPARSEABLE]


@pytest.mark.parametrize("path", _genbank_fixtures(), ids=os.path.basename)
def test_every_fixture_round_trips(path):
    """Parsing, writing, and parsing again gives the same records, and writing again the same text.

    This guards the Biopython upper bound in pyproject.toml: the layer uses private parts of
    Biopython's GenBank parser and writer.
    """
    with open(path) as handle:
        records = list(parse_genbank(handle))
    assert records
    first = io.StringIO()
    write_genbank_records(records, first)
    again = list(parse_genbank(io.StringIO(first.getvalue())))
    assert [_signature(r) for r in again] == [_signature(r) for r in records]
    second = io.StringIO()
    write_genbank_records(again, second)
    assert second.getvalue() == first.getvalue()


def test_long_qualifier_values_are_not_wrapped():
    # Biopython wraps qualifier values at 80 columns, which puts spaces into a long value without any.
    with open(os.path.join(_DATA, "pDONR201_multi_genemark_domainator.gb")) as handle:
        record = next(parse_genbank(handle))
    value = "X" * 300
    record.features[0].qualifiers["note"] = [value]
    out = io.StringIO()
    write_genbank_records([record], out)
    back = next(parse_genbank(io.StringIO(out.getvalue())))
    assert back.features[0].qualifiers["note"] == [value]


def test_long_name_genpept_locus_line():
    # A protein LOCUS line whose name is too long for its column keeps its molecule type and date.
    with open(os.path.join(_DATA, "FeSOD_20.gb")) as handle:
        records = list(parse_genbank(handle))
    assert all(r.annotations["molecule_type"] == "protein" for r in records)
    assert all("date" in r.annotations for r in records)


def test_parsing_leaves_the_warning_filters_alone():
    import warnings
    before = list(warnings.filters)
    with open(os.path.join(_DATA, "pDONR201_multi_genemark_domainator.gb")) as handle:
        list(parse_genbank(handle))
    assert warnings.filters == before
