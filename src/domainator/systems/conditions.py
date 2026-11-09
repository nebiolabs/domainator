"""The condition grammar of system models.

``&`` (and), ``|`` (or), ``~`` (not), parentheses, and ``<`` (order). ``<`` binds tightest,
then ``~``, then ``&``, then ``|``: ``~a < b`` is ``~(a < b)``. The operands of ``<`` are
names. Names are read as by ``domainator.expr``.

The functions ``cds()``, ``minimum()``, ``repeatable()``, ``minscore()`` and ``define``
aliases are format 1.1, and are refused here.

A condition is evaluated over the members of one window, in system orientation. Names that
are positive (under an even number of ``~``) are looked up among the window's members; names
under ``~`` are looked up among the members within its reach (see
docs/design/systems/detection.md), which in chain and anchor mode are the same members.
"""
from typing import FrozenSet, Optional, Sequence, Set, Tuple

from domainator.errors import InputError
from domainator.expr import tokenize

_OPERATORS = "&|~()<"
FORMAT_1_1_FUNCTIONS = ("cds", "minimum", "repeatable", "minscore")

# A parsed condition: ("name", n), ("order", (n1, n2, ...)), ("not", x), ("and", a, b), ("or", a, b).
_Node = Tuple


class Condition:
    """A parsed condition. Raises InputError when the text is malformed."""

    def __init__(self, text: str):
        self.text = str(text)
        self._tokens = tokenize(self.text, _OPERATORS)
        self._pos = 0
        if not self._tokens:
            raise InputError(f"Empty condition: {self.text!r}")
        self._tree = self._parse_or()
        if self._pos != len(self._tokens):
            self._error(f"unexpected {self._tokens[self._pos]!r}")
        positive: Set[str] = set()
        negative: Set[str] = set()
        self._collect(self._tree, True, positive, negative)
        self.positive_names: FrozenSet[str] = frozenset(positive)
        """Names used under an even number of ``~``."""
        self.negative_names: FrozenSet[str] = frozenset(negative)
        """Names used under an odd number of ``~``."""
        self.names: FrozenSet[str] = self.positive_names | self.negative_names
        if not self.positive_names:
            raise InputError(f"Condition {self.text!r} has no positive term: something must be present for a system to be found.")

    def _error(self, problem):
        raise InputError(f"Invalid condition {self.text!r}: {problem}.")

    def _peek(self, offset=0):
        pos = self._pos + offset
        return self._tokens[pos] if pos < len(self._tokens) else None

    def _parse_or(self) -> _Node:
        node = self._parse_and()
        while self._peek() == "|":
            self._pos += 1
            node = ("or", node, self._parse_and())
        return node

    def _parse_and(self) -> _Node:
        node = self._parse_not()
        while self._peek() == "&":
            self._pos += 1
            node = ("and", node, self._parse_not())
        return node

    def _parse_not(self) -> _Node:
        token = self._peek()
        if token is None:
            self._error("it ends where a name was expected")
        if token == "~":
            self._pos += 1
            return ("not", self._parse_not())
        if token == "(":
            self._pos += 1
            node = self._parse_or()
            if self._peek() != ")":
                self._error("a '(' is never closed")
            self._pos += 1
            return node
        names = [self._parse_name()]
        while self._peek() == "<":
            self._pos += 1
            names.append(self._parse_name())
        if len(names) == 1:
            return ("name", names[0])
        return ("order", tuple(names))

    def _parse_name(self) -> str:
        token = self._peek()
        if token is None:
            self._error("it ends where a name was expected")
        if token in _OPERATORS:
            self._error(f"{token!r} where a name was expected")
        if self._peek(1) == "(":
            if token in FORMAT_1_1_FUNCTIONS:
                raise InputError(f"Condition {self.text!r} uses {token}(), which needs system model format 1.1. "
                                 "This version of Domainator reads format 1.")
            self._error(f"{token!r} is followed by '(', but is not a function")
        self._pos += 1
        return token

    @classmethod
    def _collect(cls, node, positive, positive_names, negative_names):
        kind = node[0]
        if kind == "name":
            (positive_names if positive else negative_names).add(node[1])
        elif kind == "order":
            (positive_names if positive else negative_names).update(node[1])
        elif kind == "not":
            cls._collect(node[1], not positive, positive_names, negative_names)
        else:
            cls._collect(node[1], positive, positive_names, negative_names)
            cls._collect(node[2], positive, positive_names, negative_names)

    def evaluate(self, present: Sequence[FrozenSet[str]],
                 reach: Optional[Sequence[FrozenSet[str]]] = None) -> Tuple[bool, FrozenSet[str]]:
        """Whether the condition holds, and the names that made it true.

        present: the component names each member of the window fills, members in system order
            (upstream first).
        reach: the same for the members within reach, which negated names are looked up in.
            Default: ``present``.

        The names that made it true are the positive names whose presence the result relies
        on; negated terms contribute none.
        """
        if reach is None:
            reach = present
        present_names = frozenset().union(*present) if present else frozenset()
        reach_names = frozenset().union(*reach) if reach else frozenset()

        def ordered(names, members) -> bool:
            # greedy: the earliest member for each name in turn finds a chain if there is one
            position = 0
            for name in names:
                while position < len(members) and name not in members[position]:
                    position += 1
                if position == len(members):
                    return False
                position += 1
            return True

        def evaluate(node, positive) -> Tuple[bool, FrozenSet[str]]:
            kind = node[0]
            if kind == "name":
                value = node[1] in (present_names if positive else reach_names)
                return value, (frozenset((node[1],)) if value and positive else frozenset())
            if kind == "order":
                value = ordered(node[1], present if positive else reach)
                return value, (frozenset(node[1]) if value and positive else frozenset())
            if kind == "not":
                value, _ = evaluate(node[1], not positive)
                return not value, frozenset()
            left, left_names = evaluate(node[1], positive)
            right, right_names = evaluate(node[2], positive)
            if kind == "and":
                return (left and right), (left_names | right_names if left and right else frozenset())
            return (left or right), (left_names | right_names)

        return evaluate(self._tree, True)
