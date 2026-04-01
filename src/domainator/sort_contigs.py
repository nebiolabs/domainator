"""Sort contigs by Domain_Search score.

This matches the output ordering used by domain_search.py when --max_hits is set.
It is primarily intended for post-processing domain_search output written without
--max_hits, where the temporary hit_best_score annotation is not preserved in the
written GenBank file.
"""

import heapq
import sys
from typing import Iterable, Optional, Set

from jsonargparse import ActionConfigFile, ArgumentParser

from domainator import DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME, RawAndDefaultsFormatter, __version__
from domainator.utils import parse_seqfiles, write_genbank


def _get_feature_score(feature, record_id: str) -> float:
    try:
        return float(feature.qualifiers["score"][0])
    except (KeyError, IndexError, TypeError, ValueError) as exc:
        raise ValueError(
            f"Record {record_id} has a {feature.type} feature without a valid score qualifier."
        ) from exc


def set_record_best_score(record, domains: Optional[Set[str]] = None, databases: Optional[Set[str]] = None):
    best_score = None

    if domains is None:
        for feature in record.features:
            if feature.type != DOMAIN_SEARCH_BEST_HIT_NAME:
                continue
            score = _get_feature_score(feature, record.id)
            if best_score is None or score > best_score:
                best_score = score

        if best_score is None:
            raise ValueError(
                f"Record {record.id} does not contain any {DOMAIN_SEARCH_BEST_HIT_NAME} features to sort by."
            )
    else:
        for feature in record.features:
            if feature.type != DOMAIN_FEATURE_NAME:
                continue
            try:
                feature_name = feature.qualifiers["name"][0]
                feature_database = feature.qualifiers["database"][0]
            except (KeyError, IndexError, TypeError) as exc:
                raise ValueError(
                    f"Record {record.id} has a {DOMAIN_FEATURE_NAME} feature missing required name or database qualifiers."
                ) from exc

            if feature_name not in domains:
                continue
            if databases is not None and feature_database not in databases:
                continue

            score = _get_feature_score(feature, record.id)
            if best_score is None or score > best_score:
                best_score = score

        if best_score is None:
            if databases is None:
                raise ValueError(
                    f"Record {record.id} does not contain any {DOMAIN_FEATURE_NAME} features matching the requested domains."
                )
            raise ValueError(
                f"Record {record.id} does not contain any {DOMAIN_FEATURE_NAME} features matching the requested domains and databases."
            )

    record.set_hit_best_score(best_score)


def sort_contigs(records: Iterable, max_hits: Optional[int] = None, domains: Optional[Set[str]] = None, databases: Optional[Set[str]] = None):
    out_heap = []

    for record in records:
        set_record_best_score(record, domains=domains, databases=databases)
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
                        help="Maximum number of contigs to return. Prioritized by the best selected annotation score. [default: return all contigs]")
    parser.add_argument("--domains", type=str, default=None, nargs="+",
                        help="If supplied, sort by the best Domainator annotation score from domains with these names instead of Domain_Search annotation scores.")
    parser.add_argument("--databases", type=str, default=None, nargs="+",
                        help="When used with --domains, only consider Domainator annotations from these databases.")
    parser.add_argument("--config", action=ActionConfigFile)

    params = parser.parse_args(argv)

    if params.max_hits is not None and params.max_hits < 0:
        raise ValueError("max_hits must be >= 0.")

    domains = None if params.domains is None else set(params.domains)
    databases = None if params.databases is None else set(params.databases)
    if domains is None:
        databases = None

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
        domains=domains,
        databases=databases,
    ):
        write_genbank((record,), output_handle)

    if params.output is not None:
        output_handle.close()


def _entrypoint():
    main(sys.argv[1:])


if __name__ == '__main__':
    main(sys.argv[1:])