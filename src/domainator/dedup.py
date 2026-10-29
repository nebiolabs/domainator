"""Exact output-record deduplication policies shared by search and selection.

Domainator assumes that records with the same ID have the same sequence, though
their annotations and hit scores may differ. Sequence uniqueness therefore also
enforces ID uniqueness; each policy needs only one deduplication key.
"""

import heapq
from typing import Dict, List, Optional, Set, Tuple, Union

from Bio.SeqRecord import SeqRecord

from domainator.errors import InputError
from domainator.model import get_hit_best_score, sequence_hashes


DEDUP_CHOICES = ("none", "id", "sequence", "sequence_rc")
_DedupKey = Union[str, bytes, int]


def check_dedup(dedup: str) -> str:
    if not isinstance(dedup, str) or dedup not in DEDUP_CHOICES:
        raise InputError(f"dedup must be one of {', '.join(DEDUP_CHOICES)}; got {dedup!r}.")
    return dedup


def _record_key(record: SeqRecord, dedup: str) -> Optional[_DedupKey]:
    both_strands = dedup == "sequence_rc"
    if both_strands and record.annotations.get("molecule_type") == "protein":
        raise InputError("--dedup sequence_rc is only supported for nucleotide output.")
    if dedup in ("sequence", "sequence_rc"):
        return min(sequence_hashes(record.seq, both_strands))
    if dedup == "id":
        return record.id
    return None


class RecordDeduplicator:
    """Keep the first record with an unseen key under the ID-to-sequence invariant."""

    def __init__(self, dedup: str = "none"):
        self.dedup = check_dedup(dedup)
        self.seen: Set[_DedupKey] = set()

    def keep(self, record: SeqRecord) -> bool:
        key = _record_key(record, self.dedup)
        if key is None:
            return True
        if key in self.seen:
            return False
        self.seen.add(key)
        return True


class TopHits:
    """Keep at most ``max_hits`` best-scoring records, unique by the requested policy.

    Better copies replace their representative; equal scores favor earlier arrivals.
    Only retained records and O(max_hits) heap metadata are held. IDs must identify
    sequences consistently, as they do for RecordDeduplicator.
    """

    def __init__(self, max_hits: int, dedup: str = "none"):
        if max_hits < 1:
            raise InputError("max_hits must be at least 1.")
        self.max_hits = max_hits
        self.dedup = check_dedup(dedup)
        self._records: Dict[_DedupKey, Tuple[float, int, SeqRecord]] = {}
        self._heap: List[Tuple[float, int, _DedupKey]] = []
        self._arrival = 0

    def __len__(self) -> int:
        return len(self._records)

    def add(self, record: SeqRecord) -> None:
        key = _record_key(record, self.dedup)
        arrival = self._arrival
        self._arrival += 1
        if key is None:
            key = arrival
        score = get_hit_best_score(record)
        rank = (score, -arrival)
        current = self._records.get(key)
        if current is not None:
            if rank <= current[:2]:
                return
        elif len(self._records) == self.max_hits:
            while self._heap:
                root_score, root_arrival, root_key = self._heap[0]
                root = self._records.get(root_key)
                if root is not None and root[:2] == (root_score, root_arrival):
                    break
                heapq.heappop(self._heap)
            if rank <= self._heap[0][:2]:
                return
            _, _, evicted_key = heapq.heappop(self._heap)
            del self._records[evicted_key]

        self._records[key] = (score, -arrival, record)
        heapq.heappush(self._heap, (score, -arrival, key))
        # Stale entries contain no records; compact them so repeated improvements
        # cannot accumulate unbounded metadata either.
        if len(self._heap) > 2 * self.max_hits:
            self._heap = [(score, arrival, key)
                          for key, (score, arrival, _) in self._records.items()]
            heapq.heapify(self._heap)

    def records(self) -> List[SeqRecord]:
        return [record for _, _, record in sorted(self._records.values(), reverse=True)]
