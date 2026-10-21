"""The ``--columns`` option of the table reports (enum_report, hmmer_report).

Each report has a catalog of columns. ``--columns`` takes column tokens in the order they
should appear: a column's name, followed by its arguments if it has any, separated by ``:``,
for example ``--columns domains taxid:genus qualifier:CDS:product``. The last argument takes the
rest of the token, so it can contain ``:``. ``--columns`` can be given more than once.
"""
import argparse
from typing import Dict, List, NamedTuple, Sequence, Tuple, Union

from domainator.errors import InputError

VARIADIC = "..."
"""As the last argument name: the column takes one or more values for it."""


class Column(NamedTuple):
    name: str
    help: str
    args: Tuple[str, ...] = ()
    """Names of the column's arguments; a final VARIADIC means the one before it repeats."""


Analysis = Union[str, Tuple[str, List[str]]]
"""A parsed column: its name, or (name, argument values) for a column with arguments."""


def _token_format(column: Column) -> str:
    if not column.args:
        return column.name
    if column.args[-1] == VARIADIC:
        return f"{column.name}:{column.args[0]}[:{column.args[0]}...]"
    return ":".join((column.name,) + column.args)


def parse_columns(tokens: Sequence[str], catalog: Sequence[Column]) -> List[Analysis]:
    """Parse ``--columns`` tokens into the analyses a report computes, in order."""
    by_name: Dict[str, Column] = {column.name: column for column in catalog}
    analyses: List[Analysis] = []
    for token in tokens:
        name, _, rest = token.partition(":")
        if name not in by_name:
            raise InputError(f"Unknown column '{name}' in --columns. Columns: {', '.join(by_name)}.")
        column = by_name[name]
        if not column.args:
            if rest:
                raise InputError(f"Column '{name}' takes no arguments, but got '{token}'.")
            analyses.append(name)
            continue
        if not rest:
            raise InputError(f"Column '{name}' needs arguments: --columns {_token_format(column)}")
        if column.args[-1] == VARIADIC:
            values = rest.split(":")
        else:
            values = rest.split(":", len(column.args) - 1)
            if len(values) != len(column.args):
                raise InputError(f"Column '{name}' takes {len(column.args)} arguments: --columns {_token_format(column)}")
        analyses.append((name, values))
    return analyses


def columns_help(catalog: Sequence[Column]) -> str:
    lines = ["The columns to report, in order. A column with arguments takes them after ':', as in "
             "'qualifier:CDS:product'. Can be given more than once. Columns:"]
    lines += [f"{_token_format(column)}: {column.help}" for column in catalog]
    return "\n".join(lines)


class _ColumnToken:
    """The ``type`` of ``--columns``: checks each token when the command line is parsed."""
    __name__ = "column"

    def __init__(self, catalog: Sequence[Column]):
        self.catalog = catalog

    def __call__(self, token: str) -> str:
        try:
            parse_columns([token], self.catalog)
        except InputError as exc:
            raise argparse.ArgumentTypeError(str(exc)) from None
        return token


def add_columns_argument(parser, catalog: Sequence[Column]):
    """Add ``--columns`` for a report with the given column catalog."""
    return parser.add_argument("--columns", nargs="+", action="extend", type=_ColumnToken(catalog), default=None,
                               metavar="COLUMN", help=columns_help(catalog))
