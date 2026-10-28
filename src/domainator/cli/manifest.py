"""``domainator describe``: a machine-readable description of every command (the tool manifest).

The manifest is built from the commands' parsers, so it can't drift from the command line. It
is meant for programs that drive Domainator, such as a GUI: they render a form per command
from it, write the form's values to a JSON config file (keys are the parameters' ``name``s),
and run ``domainator <command> --config <file>``. Its schema is in docs/tool_manifest.md.
"""
import argparse
import importlib
import inspect
import json
import os
import sys
from typing import Any, Dict, List, Optional

MANIFEST_VERSION = 1

# Options that hold a path, by name. Command options are named for what they hold, so the
# name decides: see docs/tool_manifest.md for the rules. test_manifest.py checks that every
# string option of every command is classified here or listed in TEXT_OPTIONS.
INPUT_PATHS = {"input", "reference", "metadata", "color_table", "coords", "foldseek", "newick", "msa",
               "esm2_3di_weights", "matrix", "ncbi_summary"}
DIRECTORIES = {"output_dir", "tmp_dir", "model_dir", "tmbed_model_dir", "datasets_workdir", "ncbi_taxonomy_path"}

TEXT_OPTIONS = {
    "acc", "after", "before", "categorical", "cdhit_args", "cds",
    "cluster_sep", "color_by", "column_names", "contains", "contigs",
    "contigs_regex", "database_name", "databases", "datasets_include", "decoys", "definition_regex",
    "desc", "device", "diamond_args", "domain_expr", "domains", "email", "esm2_3di_device", "exact",
    "exclude_accessions", "foldseek_args", "foldseek_device",
    "input_format", "keep_databases", "label_by", "motif", "name", "name_field",
    "ncbi_api_key", "node", "pca_args", "regex", "reseek_args", "search_hit_color",
    "sequence_regex", "subset", "taxonomy_expr", "threshold", "title", "tsne_args", "types", "umap_args",
    "usearch_args",
}
"""String options without choices that are not paths."""

_TYPES = {"str": "string", "lower": "string", "int": "integer", "float": "number", "bool": "boolean"}
_RESTRICTED = {  # jsonargparse restricted types and Domainator's own checking types
    "PositiveInt": ("integer", {"min": 1}), "NonNegativeInt": ("integer", {"min": 0}),
    "PositiveFloat": ("number", {"exclusive_min": 0}), "NonNegativeFloat": ("number", {"min": 0}),
    "ClosedUnitInterval": ("number", {"min": 0, "max": 1}), "knn_arg": ("integer", {"min": 1}),
    "_positive_int": ("integer", {"min": 1}), "column": ("string", {}),
}
_SKIPPED = {"help", "print_config", "config"}  # every command has them; see "invocation"


def path_role(dest: str) -> Optional[str]:
    """input, output, directory or program for an option that holds a path, else None."""
    if dest.endswith("_out"):
        return "output"
    if dest in DIRECTORIES or dest.endswith("_dir"):
        return "directory"
    if dest.endswith("_path"):
        return "program"
    if dest in INPUT_PATHS or dest.endswith("_file"):
        return "input"
    return None


def _type_name(action) -> str:
    typehint = getattr(action, "_typehint", None) or action.type
    return getattr(typehint, "__name__", "str") if typehint is not None else "str"


def _json_default(value):
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, (list, tuple)):
        return [_json_default(v) for v in value]
    return str(value)


def describe_parameter(action, columns=None) -> Dict[str, Any]:
    dest = action.dest
    flags = list(action.option_strings)
    parameter: Dict[str, Any] = {"name": dest, "flags": flags, "help": action.help or ""}
    if action.nargs == 0:
        parameter.update(kind="flag", type="boolean", default=bool(action.default),
                         sets=not isinstance(action, argparse._StoreFalseAction))
        return parameter
    type_name = _type_name(action)
    json_type, bounds = _RESTRICTED.get(type_name, (_TYPES.get(type_name, "string"), {}))
    role = path_role(dest)
    kind = "columns" if columns is not None else "path" if role else "value"
    parameter.update(kind=kind, type=json_type, **bounds)
    if role:
        parameter["role"] = role
    if action.nargs in ("+", "*") or isinstance(action, (argparse._AppendAction, argparse._ExtendAction)):
        parameter["multiple"] = True
    elif isinstance(action.nargs, int) and action.nargs > 1:
        parameter["count"] = action.nargs
    parameter["required"] = bool(action.required)
    default = _json_default(action.default)
    if isinstance(default, str) and default.startswith(os.path.expanduser("~")):
        parameter["default_dynamic"] = True  # depends on the user's environment, so not recorded
    else:
        parameter["default"] = default
    if action.choices is not None:
        choices = sorted(action.choices, key=str) if isinstance(action.choices, (set, frozenset)) else action.choices
        parameter["choices"] = [_json_default(c) for c in choices]
    if role == "input" and dest == "input":
        parameter["stdin"] = getattr(action, "stdin", False)
    if role == "output":
        # add_output records it; an output defined otherwise is taken to need a file
        parameter["stdout"] = getattr(action, "stdout", False)
    if columns is not None:
        parameter["columns"] = [{"name": c.name, "args": [a for a in c.args if a != "..."],
                                 "repeat_last": bool(c.args) and c.args[-1] == "...", "help": c.help} for c in columns]
    return parameter


def describe_tool(name: str) -> Dict[str, Any]:
    from domainator.cli.registry import TOOLS

    tool = TOOLS[name]
    module = importlib.import_module(f"domainator.{tool.module}")
    parser = module.build_parser()
    parameters = []
    by_name = {}
    for action in parser._actions:
        if not action.option_strings or action.dest in _SKIPPED or isinstance(action, argparse._HelpAction):
            continue
        if action.dest in by_name and action.nargs == 0:  # --no_<name> of a switch (cli.common.add_switch)
            by_name[action.dest]["negative_flags"] = list(action.option_strings)
            continue
        columns = getattr(module, "COLUMNS", None) if action.dest == "columns" else None
        parameter = describe_parameter(action, columns)
        by_name[action.dest] = parameter
        parameters.append(parameter)
    exclusive = [[a.dest for a in group._group_actions] for group in parser._mutually_exclusive_groups]
    return {
        "name": name,
        "category": tool.category,
        "summary": tool.summary,
        "description": inspect.cleandoc(module.__doc__ or ""),  # dedented as on Python 3.13+
        "parameters": parameters,
        "required_outputs": getattr(parser, "required_outputs", None),
        "mutually_exclusive": exclusive,
    }


def manifest(tools: Optional[List[str]] = None) -> Dict[str, Any]:
    from domainator import __version__
    from domainator.cli.registry import CATEGORIES, TOOLS

    names = list(TOOLS) if tools is None else tools
    return {
        "manifest_version": MANIFEST_VERSION,
        "domainator_version": __version__,
        "invocation": {"command": ["domainator", "{tool}"], "config_option": "--config",
                       "config_keys": "parameter name", "stdio_path": "-"},
        "categories": [{"name": key, "title": title} for key, title in CATEGORIES],
        "tools": [describe_tool(name) for name in names],
    }


def main(argv: Optional[List[str]] = None) -> int:
    from domainator.cli.registry import TOOLS

    parser = argparse.ArgumentParser(prog="domainator describe", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tool", nargs="+", choices=list(TOOLS), default=None,
                        help="Describe only these commands.")
    parser.add_argument("--indent", type=int, default=None, help="Indent the JSON by this many spaces.")
    params = parser.parse_args(argv)
    json.dump(manifest(params.tool), sys.stdout, indent=params.indent)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
