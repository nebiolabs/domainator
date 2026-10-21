"""Tests for the --columns option of the table reports."""
import re

import pytest

from domainator import enum_report, hmmer_report
from domainator.errors import InputError
from domainator.report_columns import parse_columns


def test_tokens_parse_in_order():
    tokens = ["domains", "taxid:genus:species", "qualifier:CDS:a:b", "append:x:int:1", "named_domain_count:CcdB"]
    assert parse_columns(tokens, enum_report.COLUMNS) == [
        "domains",
        ("taxid", ["genus", "species"]),
        ("qualifier", ["CDS", "a:b"]),  # the last argument keeps any ':'
        ("append", ["x", "int", "1"]),
        ("named_domain_count", ["CcdB"]),
    ]


@pytest.mark.parametrize("token,message", [
    ("domainz", "Unknown column 'domainz'"),
    ("domains:x", "takes no arguments"),
    ("taxid", "needs arguments: --columns taxid:rank[:rank...]"),
    ("append:x:int", "takes 3 arguments: --columns append:name:type:value"),
])
def test_bad_tokens(token, message):
    with pytest.raises(InputError, match=re.escape(message)):
        parse_columns([token], enum_report.COLUMNS)


def test_bad_token_is_a_usage_error(capsys):
    with pytest.raises(SystemExit) as exit_info:
        hmmer_report.main(["-i", "x.hmm", "-o", "x.tsv", "--columns", "taxid:genus"])
    assert exit_info.value.code == 2
    assert "Unknown column 'taxid'" in capsys.readouterr().err


def test_repeated_option_extends():
    params = enum_report.build_parser().parse_args(["-i", "x.gb", "-o", "y", "--columns", "domains", "--columns", "length"])
    assert params.columns == ["domains", "length"]
