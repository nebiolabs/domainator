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


class Member:
    """A window member with hit scores and positions: hits is name -> [(position, score), ...]."""

    def __init__(self, **hits):
        self.hits = hits
        self.names = frozenset(hits)

    def score(self, name):
        return max(score for _, score in self.hits[name])

    def positions(self, name):
        return tuple(sorted(position for position, _ in self.hits[name]))


@pytest.mark.parametrize("text,present,expected", [
    ("cds(A & B)", members("A", "B"), False),
    ("cds(A & B)", members("A", "AB"), True),
    ("cds(A & ~B)", members("AB"), False),
    ("cds(A & ~B)", members("AB", "A"), True),
    ("~cds(A & B) & A", members("A", "B"), True),
    ("~cds(A & B) & A", members("AB"), False),
    ("minimum(2, [A, B, C])", members("A", "A"), False),
    ("minimum(2, [A, B, C])", members("A", "C"), True),
    ("minimum(0, [A]) & B", members("B"), True),
    ("repeatable(2, [A, B])", members("A", "A"), True),
    ("repeatable(3, [A, B])", members("A", "AB"), True),
    ("repeatable(3, [A, B])", members("A", "B"), False),
    ("D & ~minimum(2, [A, B])", members("D", "A"), True),
    ("D & ~minimum(2, [A, B])", members("D", "A", "B"), False),
])
def test_evaluate_functions(text, present, expected):
    assert Condition(text).evaluate(present)[0] is expected


def test_minscore():
    condition = Condition("minscore(A, 100) & B")
    assert condition.evaluate([Member(A=[(0, 99.5)]), Member(B=[(0, 1)])])[0] is False
    assert condition.evaluate([Member(A=[(0, 99.5)]), Member(A=[(0, 100)]), Member(B=[(0, 1)])]) == (True, frozenset("AB"))
    assert condition.evaluate([frozenset("A"), frozenset("B")])[0] is False  # bare sets have no scores


def test_cds_order():
    """Inside cds(), < orders the hits on one CDS, N- to C-terminal."""
    condition = Condition("cds(A < B < C)")
    assert condition.evaluate([Member(A=[(10, 1)], B=[(50, 1)], C=[(90, 1)])])[0] is True
    assert condition.evaluate([Member(A=[(10, 1)], B=[(95, 1)], C=[(90, 1)])])[0] is False
    assert condition.evaluate([Member(A=[(10, 1)], B=[(10, 1)], C=[(90, 1)])])[0] is False  # ties aren't ordered
    # B appears twice: the second copy is after C
    assert Condition("cds(C < B)").evaluate([Member(B=[(5, 1), (95, 1)], C=[(90, 1)])])[0] is True
    # across CDSs it is the members' order that counts, and one CDS doesn't order against itself
    assert Condition("A < B").evaluate([Member(A=[(10, 1)], B=[(50, 1)])])[0] is False


def test_define():
    define = {"ANY_KS": "PKS_KS | ene_KS", "KS_AT": "cds(PKS_AT & ANY_KS)"}
    condition = Condition("KS_AT | cds(PKS_AT & ~ANY_KS) & ANY_KS < PKS_AT", define)
    assert condition.aliases == {"ANY_KS", "KS_AT"}
    assert condition.names == {"PKS_KS", "ene_KS", "PKS_AT"}
    assert condition.evaluate(members({"PKS_AT", "ene_KS"}))[0] is True
    assert condition.evaluate(members({"PKS_AT"}))[0] is False
    assert condition.evaluate(members({"PKS_KS"}, {"PKS_AT"}))[0] is True
    assert condition.version == (1, 1)
    assert Condition("A & B").version == (1, 0)
    assert Condition("minimum(1, [ANY_KS])", define).evaluate(members({"ene_KS"}))[0] is True
    with pytest.raises(InputError, match="disjunctions of names"):
        Condition("ANY_KS < KS_AT", define)
    with pytest.raises(InputError, match="in terms of itself"):
        Condition("A & L1", {"L1": "B | L2", "L2": "L1"})


def test_veto_names():
    """Only a bare negated name vetoes on its own; the short form makes it forbidden."""
    condition = Condition("X & ~W & ~cds(X & V) & ~(Q < X) & ~minimum(2, [R, S])")
    assert condition.veto_names == {"W"}
    assert condition.positive_names == {"X"}
    assert condition.ordered_names == ("X", "W", "V", "Q", "R", "S")


@pytest.mark.parametrize("text,message", [
    ("~A", "no positive term"),
    ("cds(minimum(1, [A]))", "inside cds"),
    ("minimum(two, [A])", "must be a number"),
    ("minimum(1.5, [A])", "whole number"),
    ("minimum(1, [A & B])", "expected '\\]'"),
    ("minimum(1, A)", "expected '\\['"),
    ("minscore(A, 1", "expected '\\)'"),
    ("A < cds(B)", "where a name was expected"),
    ("A &", "ends where a name was expected"),
    ("(A | B", "never closed"),
    ("A < (B)", "where a name was expected"),
    ("foo(A)", "not a function"),
    ("", "Empty"),
])
def test_refused(text, message):
    with pytest.raises(InputError, match=message):
        Condition(text)
