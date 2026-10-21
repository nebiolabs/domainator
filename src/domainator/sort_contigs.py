"""Sort contigs by Domain_Search score.

This matches the output ordering used by domainator domain_search when --max_hits is set.
It is primarily intended for post-processing domain_search output written without
--max_hits, where the temporary hit_best_score annotation is not preserved in the
written GenBank file.
"""

import heapq
from typing import Iterable, Optional, Set

from jsonargparse import ActionConfigFile

from domainator.cli.common import add_input, add_seq_outputs, check_outputs, tool_parser
from domainator.outputs import write_seq_outputs
from domainator import DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME
from domainator.utils import list_and_file_to_dict_keys, parse_seqfiles


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


def build_parser():
    parser = tool_parser(__doc__)

    add_input(parser, "GenBank files.")
    add_seq_outputs(parser)
    parser.add_argument("--max_hits", type=int, default=None,
                        help="Maximum number of contigs to return, keeping the best-scoring ones. If not supplied, return all contigs.")
    parser.add_argument("--domains", type=str, default=None, nargs="+",
                        help="If supplied, sort by the best Domainator annotation score from domains with these names instead of Domain_Search annotation scores. Additive with --domains_file.")
    parser.add_argument("--domains_file", type=str, default=None,
                        help="Like --domains, but read the domain names from this file, one per line. Additive with --domains.")
    parser.add_argument("--databases", type=str, default=None, nargs="+",
                        help="When used with --domains, only consider Domainator annotations from these databases.")
    parser.add_argument("--config", action=ActionConfigFile)
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)

    if params.max_hits is not None and params.max_hits < 1:
        parser.error("--max_hits must be at least 1. Leave it out to return all contigs.")

    domains = list_and_file_to_dict_keys(params.domains, params.domains_file, as_set=True)
    databases = None if params.databases is None else set(params.databases)
    if domains is None:
        databases = None

    check_outputs(parser, params, ["genbank_out", "fasta_out"])

    write_seq_outputs(sort_contigs(
        parse_seqfiles(params.input),
        max_hits=params.max_hits,
        domains=domains,
        databases=databases,
    ), genbank=params.genbank_out, fasta=params.fasta_out)
