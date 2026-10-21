"""Boolean expressions over domain names and taxonomy (``--domain_expr``, ``--taxonomy_expr``).
"""
from typing import Iterable, List, Tuple, Union

from domainator.errors import InputError

_OPERATORS = "&|~()"

# A parsed expression: a name, or ("not", node), ("and", left, right), ("or", left, right).
_Node = Union[str, Tuple]


class BooleanEvaluator():
    """A boolean expression over names: ``&`` (and), ``|`` (or), ``~`` (not) and parentheses.

    ``~`` binds tightest, then ``&``, then ``|``, as in Python. A name is any run of other
    non-space characters. The expression is parsed once, when the evaluator is made; a
    malformed expression raises InputError.
    """

    def __init__(self, expression):
        self.expression = str(expression)
        self._tokens = self._tokenize(self.expression)
        self._pos = 0
        if not self._tokens:
            raise InputError(f"Empty boolean expression: {self.expression!r}")
        self._tree = self._parse_or()
        if self._pos != len(self._tokens):
            self._error(f"unexpected {self._tokens[self._pos]!r}")
        self.names = frozenset(t for t in self._tokens if t not in _OPERATORS)
        """Every name the expression mentions."""

    @staticmethod
    def _tokenize(expression: str) -> List[str]:
        tokens = list()
        name = ""
        for char in expression:
            if char in _OPERATORS or char.isspace():
                if name:
                    tokens.append(name)
                    name = ""
                if not char.isspace():
                    tokens.append(char)
            else:
                name += char
        if name:
            tokens.append(name)
        return tokens

    def _error(self, problem):
        raise InputError(f"Invalid boolean expression {self.expression!r}: {problem}.")

    def _peek(self):
        return self._tokens[self._pos] if self._pos < len(self._tokens) else None

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
        self._pos += 1
        if token == "~":
            return ("not", self._parse_not())
        if token == "(":
            node = self._parse_or()
            if self._peek() != ")":
                self._error("a '(' is never closed")
            self._pos += 1
            return node
        if token in _OPERATORS:
            self._error(f"{token!r} where a name was expected")
        return token

    def check_expression(self, domains: Iterable[str]) -> bool:
        """True if the expression holds when exactly the names in ``domains`` are true."""
        present = domains if isinstance(domains, (set, frozenset)) else set(domains)

        def evaluate(node):
            if isinstance(node, str):
                return node in present
            if node[0] == "not":
                return not evaluate(node[1])
            if node[0] == "and":
                return evaluate(node[1]) and evaluate(node[2])
            return evaluate(node[1]) or evaluate(node[2])

        return evaluate(self._tree)

    @classmethod
    def sanitize_identifier(cls, identifier: str) -> str:
        """
        Convert an identifier string to a name safe for boolean expressions.
        
        The BooleanEvaluator uses these special characters:
            & - AND operator
            | - OR operator  
            ~ - NOT operator
            ( - grouping open
            ) - grouping close
            space - token separator
        
        This function replaces these characters to create valid identifiers.
        
        Args:
            identifier: The string to sanitize
        
        Returns:
            A sanitized string suitable for boolean expressions
        """
        # Characters that have special meaning in BooleanEvaluator
        # We use double/triple underscores to avoid collision with single underscores
        # that might already be in the identifier
        replacements = [
            ('(', '__'),    # grouping open
            (')', '___'),   # grouping close
            ('&', '_AND_'), # AND operator
            ('|', '_OR_'),  # OR operator
            ('~', '_NOT_'), # NOT operator
            (' ', '_'),     # space separator
            ('-', '_'),     # common separator in patterns (e.g., PROSITE)
            ('<', 'Nterm_'),  # N-terminal anchor (PROSITE)
            ('>', '_Cterm'),  # C-terminal anchor (PROSITE)
            ('.', ''),      # trailing period (PROSITE)
        ]
        
        result = identifier
        for old, new in replacements:
            result = result.replace(old, new)
        
        # Remove trailing underscores
        result = result.rstrip('_')
        
        return result
