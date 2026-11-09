"""The condition grammar of system models.

``&`` (and), ``|`` (or), ``~`` (not), parentheses, and ``<`` (order). ``<`` binds tightest,
then ``~``, then ``&``, then ``|``: ``~a < b`` is ``~(a < b)``. Names are read as by
``domainator.expr``.

Format 1.1 adds antiSMASH's functions and aliases:

- ``cds(...)``: the expression inside holds on one CDS (a fusion). Inside it, ``a < b`` orders
  hits N- to C-terminal on that CDS.
- ``minimum(n, [a, b, c])``: at least n of the listed names are present.
- ``repeatable(n, [a, b, c])``: at least n (member, name) pairs among the listed names, so
  repeats count.
- ``minscore(a, s)``: a is present with a bit score of at least s.
- ``define`` aliases: a name that stands for a sub-expression. Where a list of names is
  expected (the operands of ``<`` and the lists of the functions), an alias must be a
  disjunction of names, and stands for all of them.

A condition is evaluated over the members of one window, in system orientation. Names that
are positive (under an even number of ``~``) are looked up among the window's members; names
under ``~`` are looked up among the members within its reach (see
docs/design/systems/detection.md), which in chain and anchor mode are the same members.
Inside ``cds()``, ``~`` is about that one CDS.
"""
import math
from typing import FrozenSet, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

from domainator.errors import InputError
from domainator.expr import tokenize

OPERATORS = "&|~()<,[]"
"""The characters that are tokens of a condition, and so can't be in a name."""
FUNCTIONS = ("cds", "minimum", "repeatable", "minscore")
FORMAT_1_1 = (1, 1)

# A parsed condition: ("name", n), ("order", (names1, names2, ...)), ("not", x), ("and", a, b),
# ("or", a, b), ("cds", x), ("minimum", n, names), ("repeatable", n, names), ("minscore", names, s).
# names are frozensets: an alias that is a disjunction of names stands for all of them.
_Node = Tuple


class _SetMember:
    """A window member given as a bare set of names: no scores, no hit positions."""
    __slots__ = ("names",)

    def __init__(self, names):
        self.names = frozenset(names)

    def score(self, name: str) -> float:
        return -math.inf

    def positions(self, name: str) -> Sequence[int]:
        return ()


def _member(value):
    return _SetMember(value) if isinstance(value, (set, frozenset)) else value


class _Parser:
    def __init__(self, text: str, define: Mapping[str, str], aliases: Tuple[str, ...], in_cds: bool, used: Set[str],
                 functions: Set[str]):
        self.text = text
        self.tokens = tokenize(text, OPERATORS)
        self.pos = 0
        self.define = define
        self.aliases = aliases  # the aliases being expanded, to catch cycles
        self.in_cds = in_cds
        self.used = used
        self.functions = functions

    def error(self, problem):
        where = f" (in alias {self.aliases[-1]!r})" if self.aliases else ""
        return InputError(f"Invalid condition {self.text!r}{where}: {problem}.")

    def peek(self, offset=0):
        pos = self.pos + offset
        return self.tokens[pos] if pos < len(self.tokens) else None

    def expect(self, token):
        if self.peek() != token:
            found = self.peek()
            raise self.error(f"expected {token!r}, found {'the end' if found is None else repr(found)}")
        self.pos += 1

    def parse(self) -> _Node:
        if not self.tokens:
            raise self.error("it is empty")
        node = self.parse_or()
        if self.pos != len(self.tokens):
            raise self.error(f"unexpected {self.tokens[self.pos]!r}")
        return node

    def parse_or(self) -> _Node:
        node = self.parse_and()
        while self.peek() == "|":
            self.pos += 1
            node = ("or", node, self.parse_and())
        return node

    def parse_and(self) -> _Node:
        node = self.parse_not()
        while self.peek() == "&":
            self.pos += 1
            node = ("and", node, self.parse_not())
        return node

    def parse_not(self) -> _Node:
        token = self.peek()
        if token is None:
            raise self.error("it ends where a name was expected")
        if token == "~":
            self.pos += 1
            return ("not", self.parse_not())
        if token == "(":
            self.pos += 1
            node = self.parse_or()
            if self.peek() != ")":
                raise self.error("a '(' is never closed")
            self.pos += 1
            return node
        if self.peek(1) == "(" and token not in OPERATORS:
            return self.parse_function()
        operands = [self.parse_operand()]
        while self.peek() == "<":
            self.pos += 1
            operands.append(self.parse_operand())
        if len(operands) == 1:
            return operands[0]
        return ("order", tuple(self.names_of(operand) for operand in operands))

    def parse_operand(self) -> _Node:
        """A name, or an alias as the node it stands for."""
        token = self.name_token()
        if token in self.define:
            return self.expand(token)
        return ("name", token)

    def names_of(self, node: _Node) -> FrozenSet[str]:
        if node[0] == "name":
            return frozenset((node[1],))
        if node[0] == "or":
            return self.names_of(node[1]) | self.names_of(node[2])
        raise self.error("the operands of '<' and the lists of minimum(), repeatable() and minscore() must be "
                         "names, or aliases that are disjunctions of names ('a | b')")

    def name_token(self) -> str:
        token = self.peek()
        if token is None:
            raise self.error("it ends where a name was expected")
        if token in OPERATORS:
            raise self.error(f"{token!r} where a name was expected")
        if self.peek(1) == "(":
            raise self.error(f"{token}() where a name was expected")
        self.pos += 1
        return token

    def expand(self, alias: str) -> _Node:
        if alias in self.aliases:
            raise self.error(f"alias {alias!r} is defined in terms of itself ({' -> '.join(self.aliases + (alias,))})")
        self.used.add(alias)
        parser = _Parser(str(self.define[alias]), self.define, self.aliases + (alias,), self.in_cds, self.used,
                         self.functions)
        return parser.parse()

    def number(self, what) -> float:
        token = self.peek()
        try:
            value = float(token)
        except (TypeError, ValueError):
            raise self.error(f"{what} must be a number, not {token!r}") from None
        self.pos += 1
        return value

    def name_list(self) -> FrozenSet[str]:
        self.expect("[")
        names: FrozenSet[str] = frozenset()
        while True:
            names |= self.names_of(self.parse_operand())
            if self.peek() == ",":
                self.pos += 1
                continue
            break
        self.expect("]")
        return names

    def parse_function(self) -> _Node:
        function = self.tokens[self.pos]
        if function not in FUNCTIONS:
            raise self.error(f"{function!r} is followed by '(', but is not a function (functions: {', '.join(FUNCTIONS)})")
        if self.in_cds:
            raise self.error(f"{function}() can't be used inside cds()")
        self.functions.add(function)
        self.pos += 2
        if function == "cds":
            self.in_cds = True
            node = ("cds", self.parse_or())
            self.in_cds = False
        elif function in ("minimum", "repeatable"):
            count = self.number(f"the count of {function}()")
            if count != int(count) or count < 0:
                raise self.error(f"the count of {function}() must be a whole number")
            self.expect(",")
            node = (function, int(count), self.name_list())
        else:
            names = self.names_of(self.parse_operand())
            self.expect(",")
            node = ("minscore", names, self.number("the score of minscore()"))
        self.expect(")")
        return node


class Condition:
    """A parsed condition. Raises InputError when the text is malformed.

    define: aliases, name -> expression (format 1.1).
    """

    def __init__(self, text: str, define: Optional[Mapping[str, str]] = None):
        self.text = str(text)
        define = dict(define or {})
        used: Set[str] = set()
        functions: Set[str] = set()
        if not self.text.strip():
            raise InputError(f"Empty condition: {self.text!r}")
        self._tree = _Parser(self.text, define, (), False, used, functions).parse()
        self.aliases: FrozenSet[str] = frozenset(used)
        """The aliases the condition uses, directly or through other aliases."""
        self.functions: FrozenSet[str] = frozenset(functions)
        """The functions the condition uses."""
        positive: Set[str] = set()
        negative: Set[str] = set()
        veto: Set[str] = set()
        order: List[str] = list()
        self._collect(self._tree, True, False, positive, negative, veto, order)
        self.ordered_names: Tuple[str, ...] = tuple(dict.fromkeys(order))
        """Every name, in order of first use (aliases expanded)."""
        self.positive_names: FrozenSet[str] = frozenset(positive)
        """Names whose presence can make the condition true (under an even number of ``~``)."""
        self.negative_names: FrozenSet[str] = frozenset(negative)
        """Names under an odd number of ``~``."""
        self.veto_names: FrozenSet[str] = frozenset(veto)
        """Names whose presence alone, anywhere in the window, makes the condition false where it was
        true: those negated as bare names, outside cds() and the functions."""
        self.names: FrozenSet[str] = self.positive_names | self.negative_names
        if not self.positive_names:
            raise InputError(f"Condition {self.text!r} has no positive term: something must be present for a system to be found.")

    @property
    def version(self) -> Tuple[int, int]:
        """The format version the condition needs: 1.1 when it uses functions or aliases."""
        return FORMAT_1_1 if self.functions or self.aliases else (1, 0)

    @classmethod
    def _collect(cls, node, positive, inside, positive_names, negative_names, veto_names, order):
        kind = node[0]
        if kind in ("and", "or"):
            cls._collect(node[1], positive, inside, positive_names, negative_names, veto_names, order)
            cls._collect(node[2], positive, inside, positive_names, negative_names, veto_names, order)
            return
        if kind == "not":
            cls._collect(node[1], not positive, inside, positive_names, negative_names, veto_names, order)
            return
        if kind == "cds":
            cls._collect(node[1], positive, True, positive_names, negative_names, veto_names, order)
            return
        if kind == "name":
            names = [node[1]]
            if not positive and not inside:
                veto_names.add(node[1])
        elif kind == "order":
            names = [name for operand in node[1] for name in sorted(operand)]
        elif kind == "minscore":
            names = sorted(node[1])
        else:  # minimum, repeatable
            names = sorted(node[2])
        (positive_names if positive else negative_names).update(names)
        order.extend(names)

    def evaluate(self, present: Sequence, reach: Optional[Sequence] = None) -> Tuple[bool, FrozenSet[str]]:
        """Whether the condition holds, and the names that made it true.

        present: the members of the window, in system order (upstream first). Each is a set of the
            component names it fills, or an object with ``names`` (that set), ``score(name)`` (the
            best bit score for a name) and ``positions(name)`` (the starts of its hits to a name,
            N- to C-terminal along the protein), which minscore() and cds(a < b) need.
        reach: the same for the members within reach, which negated terms are looked up in.
            Default: ``present``.

        The names that made it true are the positive names whose presence the result relies
        on; negated terms contribute none.
        """
        present = [_member(m) for m in present]
        reach = present if reach is None else [_member(m) for m in reach]
        present_names = frozenset().union(*(m.names for m in present))
        reach_names = frozenset().union(*(m.names for m in reach))

        def ordered(operands, members) -> bool:
            # greedy: the earliest member for each operand in turn finds a chain if there is one
            position = 0
            for names in operands:
                while position < len(members) and not (names & members[position].names):
                    position += 1
                if position == len(members):
                    return False
                position += 1
            return True

        def evaluate(node, positive) -> Tuple[bool, FrozenSet[str]]:
            kind = node[0]
            members = present if positive else reach
            names = present_names if positive else reach_names
            if kind == "name":
                value = node[1] in names
                return value, (frozenset((node[1],)) if value and positive else frozenset())
            if kind == "order":
                value = ordered(node[1], members)
                return value, (frozenset().union(*node[1]) & names if value and positive else frozenset())
            if kind == "not":
                value, _ = evaluate(node[1], not positive)
                return not value, frozenset()
            if kind in ("and", "or"):
                left, left_names = evaluate(node[1], positive)
                right, right_names = evaluate(node[2], positive)
                if kind == "and":
                    return (left and right), (left_names | right_names if left and right else frozenset())
                return (left or right), (left_names | right_names)
            if kind == "cds":
                value = False
                contributors: FrozenSet[str] = frozenset()
                for member in members:
                    one, one_names = _evaluate_cds(node[1], member, positive)
                    if one:
                        value = True
                        contributors |= one_names
                return value, (contributors if positive else frozenset())
            if kind == "minimum":
                found = node[2] & names
                value = len(found) >= node[1]
                return value, (found if value and positive else frozenset())
            if kind == "repeatable":
                value = sum(len(m.names & node[2]) for m in members) >= node[1]
                return value, (node[2] & names if value and positive else frozenset())
            # minscore
            found = frozenset(name for m in members for name in m.names & node[1] if m.score(name) >= node[2])
            return bool(found), (found if positive else frozenset())

        return evaluate(self._tree, True)


def _evaluate_cds(node, member, positive) -> Tuple[bool, FrozenSet[str]]:
    """A cds() expression on one member."""
    kind = node[0]
    if kind == "name":
        value = node[1] in member.names
        return value, (frozenset((node[1],)) if value and positive else frozenset())
    if kind == "order":
        last = -math.inf
        for names in node[1]:
            after = [p for name in names & member.names for p in member.positions(name) if p > last]
            if not after:
                return False, frozenset()
            last = min(after)
        return True, (frozenset().union(*node[1]) & member.names if positive else frozenset())
    if kind == "not":
        value, _ = _evaluate_cds(node[1], member, not positive)
        return not value, frozenset()
    left, left_names = _evaluate_cds(node[1], member, positive)
    right, right_names = _evaluate_cds(node[2], member, positive)
    if kind == "and":
        return (left and right), (left_names | right_names if left and right else frozenset())
    return (left or right), (left_names | right_names)


def alias_closure(define: Mapping[str, str], used: Iterable[str]) -> Tuple[Tuple[str, str], ...]:
    """The aliases in ``used``, as sorted (name, expression) pairs. ``used`` is a Condition's
    ``aliases``, which already holds the aliases those use in turn."""
    return tuple(sorted((name, str(define[name])) for name in used))
