"""The system condition grammar (domainator.systems.conditions)."""
import pytest

from domainator.errors import InputError
from domainator.systems.conditions import Condition


def members(*names):
    """Members in system order: each argument is the set of component names of one member, like "A" or "AB"."""
    return [frozenset(name) for name in names]


@pytest.mark.parametrize("text,present,expected", [
    ("A & B", members("A", "B"), True),
    ("A & B", members("A"), False),
    ("A | B", members("B"), True),
    ("A & (B | C)", members("A", "C"), True),
    ("A & ~B", members("A"), True),
    ("A & ~B", members("A", "B"), False),
    ("A < B", members("A", "B"), True),
    ("A < B", members("B", "A"), False),
    ("A < B", members("AB"), False),  # one member doesn't order against itself
    ("A < B < C", members("A", "B", "C"), True),
    ("A < B < C", members("A", "BC"), False),  # a chain needs distinct members
    ("A < B < C", members("A", "C", "B", "C"), True),
    ("~A < B & B", members("B", "A"), True),  # ~(A < B)
    ("~A < B & B", members("A", "B"), False),
    ("A < B & C", members("A", "B", "C"), True),  # (A < B) & C
])
def test_evaluate(text, present, expected):
    assert Condition(text).evaluate(present)[0] is expected


def test_negated_names_read_the_reach():
    condition = Condition("X & Y & ~W")
    assert condition.evaluate(members("X", "Y"), members("X", "Y", "W")) == (False, frozenset())
    assert condition.evaluate(members("X", "Y"), members("X", "Y")) == (True, frozenset("XY"))
    # positive names are read from the members only, not from the reach
    assert condition.evaluate(members("X"), members("X", "Y"))[0] is False


def test_contributors():
    condition = Condition("(A & B) | C")
    assert condition.evaluate(members("A", "C")) == (True, frozenset("C"))
    assert condition.evaluate(members("A", "B", "C")) == (True, frozenset("ABC"))


def test_names_and_polarity():
    condition = Condition("X & (Y | ~Z) & ~(W < V)")
    assert condition.positive_names == {"X", "Y"}
    assert condition.negative_names == {"Z", "W", "V"}


def test_names_read_like_expr():
    condition = Condition("db:PF00001.2 & Lamassu-Fam")
    assert condition.names == {"db:PF00001.2", "Lamassu-Fam"}


@pytest.mark.parametrize("text,message", [
    ("cds(A & B)", "format 1.1"),
    ("minimum(2, [A, B])", "format 1.1"),
    ("~A", "no positive term"),
    ("A &", "ends where a name was expected"),
    ("(A | B", "never closed"),
    ("A < (B)", "where a name was expected"),
    ("foo(A)", "not a function"),
    ("", "Empty"),
])
def test_refused(text, message):
    with pytest.raises(InputError, match=message):
        Condition(text)
