"""Inspect system models, and write them as a model file, a table, or the profiles they need

Reads system models (--system_models, and/or one model defined on the command line with
--domain_expr and a window, as in find_systems) and writes:

- --yaml_out: a native system model file. Models that fit the short form are written in it, so
  this is how a command-line model becomes a file:

    domainator format_systems --domain_expr "X & Y & (Q | R)" --window_cds 10 --system_name XYZ_QR --yaml_out xyz.yaml

- --tsv_out: one row per model, component and reference, with its thresholds, to see what a
  model (for example a converted one) means;
- --hmm_out: the HMM profiles that the selected models name, from their reference_files, ready
  for domainate -r.

--systems/--systems_file and --system_categories select models, so --systems with --yaml_out
copies one model out of a collection into a file of its own, to edit.
"""
from dataclasses import replace
import os
from typing import Iterable, List, Optional, Sequence, Set, Union
import warnings

from jsonargparse import ActionConfigFile

from domainator.cli.common import (add_list_and_file, add_system_models, add_window_model, check_outputs,
                                   declare_outputs, list_and_file_to_dict_keys, tool_parser)
from domainator.errors import InputError
from domainator.hmmio import iter_hmms, open_writable_hmm_file
from domainator.seqio import db_name_from_path
from domainator.systems.model import COMMAND_LINE_DATABASE, ModelCollection
from domainator.systems.native import read_system_models, reference_paths, select_models, window_model, write_system_models

TSV_COLUMNS = ("database", "system", "category", "window", "min_core", "min_total", "condition", "component", "role",
               "reference", "reference_database", "evalue", "min_score", "min_reference_coverage",
               "min_target_coverage", "min_target_length", "max_target_length", "min_tmscore", "min_lddt", "min_prob",
               "exchangeable", "attach_to_cds")
_HMM_SUFFIXES = (".hmm", ".hmm.gz")


def _cell(value) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return f"{value:g}"
    return str(value).replace("\t", " ")


def _tsv_rows(collections: Sequence[ModelCollection]) -> Iterable[List[str]]:
    for collection in collections:
        for model in collection.models:
            for component in model.components:
                for rule in component.references:
                    yield [_cell(value) for value in (
                        model.database, model.name, model.category, model.window.describe(), model.min_core,
                        model.min_total, model.condition, component.name, component.role, rule.name, rule.database,
                        rule.evalue, rule.min_score, rule.min_reference_coverage, rule.min_target_coverage,
                        rule.min_target_length, rule.max_target_length, rule.min_tmscore, rule.min_lddt, rule.min_prob,
                        rule.exchangeable, rule.attach_to_cds)]


def _write_hmms(collections: Sequence[ModelCollection], hmm_out: str) -> None:
    wanted = {(rule.name, rule.database) for collection in collections for model in collection.models
              for component in model.components for rule in component.references}
    written = set()
    files = [file for collection in collections for file in reference_paths(collection)]
    if not files:
        raise InputError("The selected models list no reference_files, so there are no profiles to write.")
    out_database = db_name_from_path(hmm_out)
    for database in sorted({database for _, database in wanted if database is not None and database != out_database}):
        warnings.warn(f"Some components only accept hits from database {database!r}, but profiles annotated from "
                      f"{hmm_out} are labeled {out_database!r}. Name the file {database}.hmm to keep them matching.")
    handle = open_writable_hmm_file(hmm_out)
    try:
        for file in files:
            if not file.lower().endswith(_HMM_SUFFIXES):
                warnings.warn(f"Skipping reference file {file}: --hmm_out only copies HMM profiles (.hmm files).")
                continue
            if not os.path.exists(file):
                raise InputError(f"Reference file not found: {file}")
            database = db_name_from_path(file)
            for hmm in iter_hmms(file):
                if hmm.name not in written and ((hmm.name, None) in wanted or (hmm.name, database) in wanted):
                    hmm.write(handle)
                    written.add(hmm.name)
    finally:
        handle.close()
    missing = sorted({name for name, _ in wanted} - written)
    if missing:
        warnings.warn(f"No profile found in the reference files for: {', '.join(missing)}.")


def format_systems(system_models: Optional[Sequence[Union[str, os.PathLike, ModelCollection]]] = None, *,
                   domain_expr: Optional[str] = None, window_cds: Optional[int] = None,
                   window_kb: Optional[float] = None, window_mode: str = "span", system_name: str = "domain_expr",
                   same_strand: bool = False, systems: Optional[Set[str]] = None,
                   system_categories: Optional[Set[str]] = None, yaml_out: Optional[str] = None,
                   tsv_out: Optional[str] = None, hmm_out: Optional[str] = None) -> List[ModelCollection]:
    """Read system models, select some, and write them as YAML, a table and/or HMM profiles.

    Models come from system_models (native model files, directories holding systems.yaml, or
    ModelCollections) and from domain_expr with window_cds and/or window_kb, as in
    find_systems. systems (``name`` or ``database:name``) and system_categories select models.

    yaml_out writes the selected models as a native model file; they must come from one
    collection. A command-line model is written in the collection named after the file.
    tsv_out writes one row per model, component and reference. hmm_out writes the HMM
    profiles the selected models name, from their reference files.

    Returns:
        The selected collections.
    """
    collections = read_system_models(system_models or [])
    if domain_expr is not None:
        collections.append(window_model(domain_expr, window_cds=window_cds, window_kb=window_kb, window_mode=window_mode,
                                        system_name=system_name, same_strand=same_strand))
    if not collections:
        raise InputError("No system models: give system_models (--system_models), domain_expr (--domain_expr), or both.")
    collections = select_models(collections, systems, system_categories)
    if not collections:
        raise InputError("No models match --systems and --system_categories.")

    if yaml_out is not None:
        if len(collections) > 1:
            raise InputError("--yaml_out writes one collection, but the selected models come from "
                             f"{', '.join(c.database for c in collections)}. Select models from one with --systems.")
        collection, = collections
        if collection.database == COMMAND_LINE_DATABASE:
            database = db_name_from_path(yaml_out)
            collection = replace(collection, database=database,
                                 models=tuple(replace(model, database=database) for model in collection.models))
        write_system_models(collection, yaml_out)
    if tsv_out is not None:
        with open(tsv_out, "w") as handle:
            print("\t".join(TSV_COLUMNS), file=handle)
            for row in _tsv_rows(collections):
                print("\t".join(row), file=handle)
    if hmm_out is not None:
        _write_hmms(collections, hmm_out)
    return collections


def build_parser():
    parser = tool_parser(__doc__)
    add_system_models(parser)
    add_window_model(parser)
    add_list_and_file(parser, "systems", help="Only the models with these names (name or database:name).")
    parser.add_argument("--system_categories", type=str, nargs="+", default=None,
                        help="Only the models in these categories (e.g. defense or PKS). With --systems, a model must "
                             "match both.")
    parser.add_argument("--yaml_out", type=str, default=None,
                        help="Write the selected models as a native system model file. They must come from one "
                             "collection; a command-line model is written in a collection named after this file.")
    parser.add_argument("--tsv_out", type=str, default=None,
                        help="Write one row per model, component and reference, with its thresholds.")
    parser.add_argument("--hmm_out", type=str, default=None,
                        help="Write the HMM profiles the selected models name, from their reference files, for "
                             "domainate -r. Components that only accept hits from one database match hits labeled "
                             "with that database, which domainate takes from the file name, so name the file after it.")
    parser.add_argument("--config", action=ActionConfigFile)
    declare_outputs(parser, ["yaml_out", "tsv_out", "hmm_out"])
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)
    check_outputs(parser, params)
    format_systems(
        params.system_models,
        domain_expr=params.domain_expr,
        window_cds=params.window_cds,
        window_kb=params.window_kb,
        window_mode=params.window_mode,
        system_name=params.system_name,
        same_strand=params.same_strand,
        systems=list_and_file_to_dict_keys(params.systems, params.systems_file, as_set=True),
        system_categories=params.system_categories,
        yaml_out=params.yaml_out,
        tsv_out=params.tsv_out,
        hmm_out=params.hmm_out,
    )
