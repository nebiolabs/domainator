"""Candidate windows for system detection: chain, span and anchor.

The functions here work on *items*: objects with ``start`` and ``end`` (bp, on an unwrapped
axis, so ``end`` may pass the contig length on a circular contig), ``rank`` (the number of
real CDSs before the item), ``is_cds`` (False for the pseudo-CDS of a contig-level nucleotide
hit, which never counts as a CDS), ``gene`` (an identity shared by the copies of one item on a
circular contig) and ``shifted(k)`` (a copy moved k times around the contig). Items are sorted
by ``start``. See docs/design/systems/detection.md, step 5.
"""
from typing import Callable, Dict, Hashable, List, Optional, Sequence, Tuple, TypeVar

Item = TypeVar("Item")


def cds_between(a, b) -> int:
    """Real CDSs strictly between item a and a later item b."""
    return b.rank - a.rank - (1 if a.is_cds else 0)


def cds_inclusive(items: Sequence) -> int:
    """Real CDSs from the first item to the last, inclusive."""
    first, last = items[0], items[-1]
    return last.rank + (1 if last.is_cds else 0) - first.rank


def span_bp(items: Sequence) -> int:
    """Length of the envelope of ``items``, in bp."""
    return max(item.end for item in items) - items[0].start


def gap_bp(a, b) -> int:
    """bp from the end of item a to the start of a later item b (negative if they overlap)."""
    return b.start - a.end


def fits_span(items: Sequence, max_cds: Optional[int], max_kb: Optional[float]) -> bool:
    """Whether ``items`` cover at most max_cds CDSs (inclusive) and at most max_kb kb."""
    if max_cds is not None and cds_inclusive(items) > max_cds:
        return False
    if max_kb is not None and span_bp(items) > max_kb * 1000:
        return False
    return True


def within_radius(a, b, max_cds: Optional[int], max_kb: Optional[float]) -> bool:
    """Whether items a and b are within max_cds intervening CDSs and max_kb kb (edge to edge) of each other."""
    if b.start < a.start:
        a, b = b, a
    if max_cds is not None and cds_between(a, b) > max_cds:
        return False
    if max_kb is not None and max(gap_bp(a, b), 0) > max_kb * 1000:
        return False
    return True


def unwrap(items: Sequence[Item], linked: Callable[[Item, Item], bool], circular: bool) -> List[Item]:
    """``items`` in the order to chain them in.

    On a linear contig that is ``items``. On a circular contig, the list is rotated to start
    after a break: a pair of neighbors (cyclically) that are not linked, or if every pair is
    linked, the pair with the largest gap. Items that come after the origin in the rotated order
    are shifted once around the contig, so the order is increasing.
    """
    items = list(items)
    if not circular or len(items) < 2:
        return items
    n = len(items)
    following = [items[i + 1] if i < n - 1 else items[0].shifted(1) for i in range(n)]
    unlinked = [i for i in range(n) if not linked(items[i], following[i])]
    if unlinked:
        cut = unlinked[0]
    else:
        cut = max(range(n), key=lambda i: gap_bp(items[i], following[i]))
    start = (cut + 1) % n
    return items[start:] + [item.shifted(1) for item in items[:start]]


def chain(items: Sequence[Item], linked: Callable[[Item, Item], bool]) -> List[List[Item]]:
    """Split ``items`` (in chaining order) into maximal runs of consecutively linked items."""
    clusters: List[List[Item]] = list()
    for item in items:
        if clusters and linked(clusters[-1][-1], item):
            clusters[-1].append(item)
        else:
            clusters.append([item])
    return clusters


def span_candidates(items: Sequence[Item], fits: Callable[[Sequence[Item]], bool],
                    evaluate: Callable[[Sequence[Item], Sequence[Item]], Optional[object]],
                    starts: Optional[range] = None) -> List[Tuple[int, int, object]]:
    """The minimal passing runs of a span window.

    Every contiguous run of ``items`` that ``fits`` is evaluated with ``evaluate(run, reach)``,
    where the reach is every item that could be added to the run with the result still fitting.
    ``evaluate`` returns None for a run that fails. A passing run is a candidate when neither of
    the runs one item shorter passes. ``starts`` limits the first item of a candidate (on a
    circular contig, to one copy of the items); runs never hold two copies of one gene.

    Returns (first index, last index, evaluate's result) of each candidate.
    """
    n = len(items)
    passing: Dict[Tuple[int, int], object] = dict()
    first_starts = range(n) if starts is None else range(max(starts.start - 1, 0), min(starts.stop + 1, n))
    for i in first_starts:
        genes = set()
        for j in range(i, n):
            if items[j].gene in genes or not fits(items[i:j + 1]):
                break
            genes.add(items[j].gene)
            left = i
            while left > 0 and items[left - 1].gene not in genes and fits(items[left - 1:j + 1]):
                left -= 1
            right = j
            while right + 1 < n and items[right + 1].gene not in genes and fits(items[i:right + 2]):
                right += 1
            reach = [item for item in items[left:right + 1] if item.gene not in genes] + list(items[i:j + 1])
            result = evaluate(items[i:j + 1], _dedup_genes(sorted(reach, key=lambda item: item.start)))
            if result is not None:
                passing[(i, j)] = result
    candidates = list()
    for (i, j), result in passing.items():
        if starts is not None and i not in starts:
            continue
        if j > i and ((i + 1, j) in passing or (i, j - 1) in passing):
            continue
        candidates.append((i, j, result))
    return candidates


def _dedup_genes(items: Sequence[Item]) -> List[Item]:
    seen = set()
    out = list()
    for item in items:
        if item.gene not in seen:
            seen.add(item.gene)
            out.append(item)
    return out


def select_greedily(candidates: Sequence[Tuple[Sequence[Item], object]],
                    key: Callable[[Sequence[Item], object], tuple]) -> List[Tuple[Sequence[Item], object]]:
    """Pick candidates in ``key`` order, skipping any that shares a gene with one already picked."""
    picked = list()
    used = set()
    for run, result in sorted(candidates, key=lambda candidate: key(*candidate)):
        genes = {item.gene for item in run}
        if genes & used:
            continue
        used |= genes
        picked.append((run, result))
    return picked


def anchor_window(items: Sequence[Item], anchor: Item, max_cds: Optional[int], max_kb: Optional[float]) -> List[Item]:
    """The items within the radius of ``anchor``, one copy per gene (the nearest)."""
    window = [item for item in items if within_radius(anchor, item, max_cds, max_kb)]
    nearest: Dict[Hashable, Item] = dict()
    for item in window:
        current = nearest.get(item.gene)
        if current is None or abs(item.start - anchor.start) < abs(current.start - anchor.start):
            nearest[item.gene] = item
    return sorted(nearest.values(), key=lambda item: item.start)
