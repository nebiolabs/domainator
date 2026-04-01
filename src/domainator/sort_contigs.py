"""Sort contigs by Domain_Search score.

This matches the output ordering used by domain_search.py when --max_hits is set.
It is primarily intended for post-processing domain_search output written without
--max_hits, where the temporary hit_best_score annotation is not preserved in the
written GenBank file.
"""

import heapq
import sys
from typing import Iterable, Optional

from jsonargparse import ActionConfigFile, ArgumentParser

from domainator import DOMAIN_SEARCH_BEST_HIT_NAME, RawAndDefaultsFormatter, __version__
from domainator.utils import parse_seqfiles, write_genbank


def set_record_best_score(record):
    best_score = None
    for feature in record.features:
        if feature.type != DOMAIN_SEARCH_BEST_HIT_NAME:
            continue
        try:
            score = float(feature.qualifiers["score"][0])
        except (KeyError, IndexError, TypeError, ValueError) as exc:
            raise ValueError(
                f"Record {record.id} has a {DOMAIN_SEARCH_BEST_HIT_NAME} feature without a valid score qualifier."
            ) from exc
        if best_score is None or score > best_score:
            best_score = score

    if best_score is None:
        raise ValueError(
            f"Record {record.id} does not contain any {DOMAIN_SEARCH_BEST_HIT_NAME} features to sort by."
        )

    record.set_hit_best_score(best_score)


def sort_contigs(records: Iterable, max_hits: Optional[int] = None):
    out_heap = []

    for record in records:
        set_record_best_score(record)
        if max_hits is None or len(out_heap) < max_hits:
            heapq.heappush(out_heap, record)
        else:
            heapq.heappushpop(out_heap, record)

    out_heap.sort(reverse=True)
    for record in out_heap:
        yield record


def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument("-i", "--input", default=None, nargs="+", type=str, required=False,
                        help="Input genbank files. If not supplied, reads from stdin.")
    parser.add_argument("-o", "--output", default=None, type=str, required=False,
                        help="Output genbank filename. If not supplied, writes to stdout.")
    parser.add_argument("--max_hits", type=int, default=None,
                        help="Maximum number of contigs to return. Prioritized by bitscore of the best Domain_Search annotation. [default: return all contigs]")
    parser.add_argument("--config", action=ActionConfigFile)

    params = parser.parse_args(argv)

    if params.max_hits is not None and params.max_hits < 0:
        raise ValueError("max_hits must be >= 0.")

    if params.input is None:
        input_files = (sys.stdin,)
        filetype_override = "genbank"
    else:
        input_files = params.input
        filetype_override = None

    if params.output is None:
        output_handle = sys.stdout
    else:
        output_handle = open(params.output, "w")

    max_hits = None if params.max_hits == 0 else params.max_hits

    for record in sort_contigs(
        parse_seqfiles(input_files, filetype_override=filetype_override),
        max_hits=max_hits,
    ):
        write_genbank((record,), output_handle)

    if params.output is not None:
        output_handle.close()


def _entrypoint():
    main(sys.argv[1:])


if __name__ == '__main__':
    main(sys.argv[1:])