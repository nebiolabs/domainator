"""Annotate systems: sets of domains found together in one genomic neighborhood.

Reads records that already carry domain annotations (from domainate), finds the systems that
the models describe, and adds a Domain_System feature for each one. A system is defined in a
system model file (--system_models), or on the command line with --domain_expr and a window:

    domainator find_systems -i annotated.gb --domain_expr "X & Y & (Q | R)" --window_cds 10 -o systems.gb

Detection uses only the annotations already in the records, so run domainate permissively
enough (a loose -e, no --max_hits_per_cds, --max_overlap 1) that the hits the models need are
there. See docs/file_formats.md for the model format and the Domain_System qualifiers.
"""
import os
from typing import Iterable, Iterator, Optional, Sequence, Set, Union

from Bio.SeqRecord import SeqRecord
from jsonargparse import ActionConfigFile
from jsonargparse.typing import PositiveFloat

from domainator import SYSTEM_FEATURE_NAME
from domainator.cli.common import (add_annotation_type, add_input, add_list_and_file, add_seq_outputs, add_switch,
                                   add_system_detection, add_system_models, add_window_model, check_outputs,
                                   declare_outputs, list_and_file_to_dict_keys, partial_help, tool_parser)
from domainator.errors import InputError
from domainator.model import ANNOTATION_TYPES
from domainator.outputs import write_seq_outputs
from domainator.partial import PARTIAL_CHOICES
from domainator.seqio import parse_seqfiles
from domainator.systems.detect import SystemDetector
from domainator.systems.features import system_feature
from domainator.systems.model import COMMAND_LINE_DATABASE, ModelCollection
from domainator.systems.native import read_system_models, select_models, window_model


def find_systems(records: Iterable[SeqRecord], *,
                 system_models: Optional[Sequence[Union[str, os.PathLike, ModelCollection]]] = None,
                 domain_expr: Optional[str] = None, window_cds: Optional[int] = None, window_kb: Optional[float] = None,
                 window_mode: str = "span", system_name: str = "domain_expr", same_strand: bool = False,
                 systems: Optional[Set[str]] = None, evalue: Optional[float] = None, annotation_type: str = "domainator",
                 databases: Optional[Set[str]] = None, min_score: Optional[float] = None,
                 min_completeness: Optional[float] = None, max_system_overlap: float = 1.0,
                 max_systems_per_contig: Optional[int] = None, truncated_systems: str = "exclude",
                 clear: bool = True, partial: str = "include") -> Iterator[SeqRecord]:
    """Annotate each record with the systems the models find, as Domain_System features.

    Models come from system_models (native model files, directories holding systems.yaml, or
    already-read ModelCollections) and from domain_expr with window_cds and/or window_kb (a
    model named system_name in the command_line database); at least one is required. systems
    restricts the models to these names (``name`` or ``database:name``).

    A hit counts as evidence when it passes its model's thresholds and evalue, comes from one
    of databases (None: any), and is one of the annotation_type annotations. partial selects
    which CDSs can be members, by fragment status. Systems scoring below min_score, or with a
    completeness below min_completeness, are dropped (models without a quorum have no
    completeness, and are not dropped by it). max_system_overlap and max_systems_per_contig
    resolve overlapping systems. truncated_systems="include" also reports systems that a linear
    contig end may have cut (see SystemDetector).

    With clear, existing Domain_System features from the models' databases are removed first
    (for the command-line model, only those with its name), so that re-running replaces systems.
    Records are yielded as they are annotated, so the input can be a stream.
    """
    collections = read_system_models(system_models or [])
    if domain_expr is not None:
        collections.append(window_model(domain_expr, window_cds=window_cds, window_kb=window_kb, window_mode=window_mode,
                                        system_name=system_name, same_strand=same_strand))
    if not collections:
        raise InputError("No system models: give system_models (--system_models), domain_expr (--domain_expr), or both.")
    if systems is not None:
        collections = select_models(collections, systems)
    yield from annotate_systems(records, collections, evalue=evalue, annotation_type=annotation_type, databases=databases,
                                min_score=min_score, min_completeness=min_completeness,
                                max_system_overlap=max_system_overlap, max_systems_per_contig=max_systems_per_contig,
                                truncated_systems=truncated_systems, clear=clear, partial=partial)


def annotate_systems(records: Iterable[SeqRecord], collections: Sequence[ModelCollection], *,
                     evalue: Optional[float] = None, annotation_type: str = "domainator",
                     databases: Optional[Set[str]] = None, min_score: Optional[float] = None,
                     min_completeness: Optional[float] = None, max_system_overlap: float = 1.0,
                     max_systems_per_contig: Optional[int] = None, truncated_systems: str = "exclude",
                     clear: bool = True, partial: str = "include") -> Iterator[SeqRecord]:
    """Add a Domain_System feature to each record for each system the collections find, as
    find_systems does, yielding the records as they are annotated. Protein records get no systems."""
    if partial not in PARTIAL_CHOICES:
        raise ValueError(f"partial must be one of {', '.join(PARTIAL_CHOICES)}, not {partial}")
    detector = SystemDetector(collections, evalue=evalue, databases=databases, annotation_type=annotation_type,
                              partial=partial, max_system_overlap=max_system_overlap,
                              max_systems_per_contig=max_systems_per_contig, truncated_systems=truncated_systems)
    cleared_databases = {c.database for c in collections if c.database != COMMAND_LINE_DATABASE}
    cleared_models = {(m.database, m.name) for c in collections if c.database == COMMAND_LINE_DATABASE for m in c.models}

    def cleared(feature) -> bool:
        if feature.type != SYSTEM_FEATURE_NAME:
            return False
        database = feature.qualifiers.get("database", [None])[0]
        return database in cleared_databases or (database, feature.qualifiers.get("name", [None])[0]) in cleared_models

    for record in records:
        if clear:
            record.features = [feature for feature in record.features if not cleared(feature)]
        if record.annotations.get("molecule_type") == "protein":
            yield record
            continue
        for hit in detector.detect(record):
            if min_score is not None and hit.score < min_score:
                continue
            if min_completeness is not None and hit.completeness is not None and hit.completeness < min_completeness:
                continue
            record.features.append(system_feature(hit))
        yield record


def build_parser():
    parser = tool_parser(__doc__)
    add_input(parser, "GenBank files with domain annotations.")
    add_seq_outputs(parser)
    add_system_models(parser)
    add_window_model(parser)
    add_list_and_file(parser, "systems", help="Only look for the models with these names (name or database:name).")
    parser.add_argument("-e", "--evalue", type=PositiveFloat, default=None,
                        help="Only hits with e-values at or below this count as evidence, in addition to each model's own "
                             "thresholds. Default: no cutoff.")
    add_annotation_type(parser, choices=ANNOTATION_TYPES, help="Which annotations count as evidence.")
    parser.add_argument("--databases", type=str, default=None, nargs="+",
                        help="Only hits from these databases count as evidence. Default: all databases.")
    parser.add_argument("--min_score", type=float, default=None,
                        help="Drop systems with a score below this (the sum of the best bit scores of their counted components).")
    parser.add_argument("--min_completeness", type=float, default=None,
                        help="Drop systems with a completeness below this (counted components found / model components). "
                             "Models without a quorum (condition only, like --domain_expr) have no completeness and are kept.")
    add_system_detection(parser)
    add_switch(parser, "clear", "Remove existing Domain_System features from the models' databases first "
                                "(for --domain_expr, those with the same --system_name).",
               default=True, negative_help="Keep existing Domain_System features.")
    parser.add_argument("--partial", type=str, default="include", choices=PARTIAL_CHOICES,
                        help=partial_help("Which CDSs can be system members"))
    parser.add_argument("--config", action=ActionConfigFile)
    declare_outputs(parser, ["genbank_out", "fasta_out"])
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)
    check_outputs(parser, params)
    if params.system_models is None and params.domain_expr is None:
        parser.error("Give --system_models, --domain_expr, or both.")
    if params.domain_expr is None and (params.window_cds is not None or params.window_kb is not None):
        parser.error("--window_cds and --window_kb define the --domain_expr system, so they need --domain_expr.")
    systems = list_and_file_to_dict_keys(params.systems, params.systems_file, as_set=True)
    write_seq_outputs(find_systems(
        parse_seqfiles(params.input),
        system_models=params.system_models,
        domain_expr=params.domain_expr,
        window_cds=params.window_cds,
        window_kb=params.window_kb,
        window_mode=params.window_mode,
        system_name=params.system_name,
        same_strand=params.same_strand,
        systems=systems,
        evalue=params.evalue,
        annotation_type=params.annotation_type,
        databases=None if params.databases is None else set(params.databases),
        min_score=params.min_score,
        min_completeness=params.min_completeness,
        max_system_overlap=params.max_system_overlap,
        max_systems_per_contig=params.max_systems_per_contig,
        truncated_systems=params.truncated_systems,
        clear=params.clear,
        partial=params.partial,
    ), genbank=params.genbank_out, fasta=params.fasta_out)
