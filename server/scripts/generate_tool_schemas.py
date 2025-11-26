#!/usr/bin/env python3
"""Generate Domainator tool schema JSON files from CLI definitions."""
from __future__ import annotations

import argparse
import importlib
import json
import re
import sys
import typing
from argparse import SUPPRESS, _AppendConstAction, _StoreFalseAction, _StoreTrueAction
from pathlib import Path
from types import ModuleType
from typing import Any, Dict, List, Optional, get_args, get_origin

try:  # Python 3.10+
    from types import UnionType as TypesUnionType  # type: ignore[attr-defined]
except ImportError:  # Python < 3.10
    TypesUnionType = None  # type: ignore[assignment]

TYPING_ANNOTATED = getattr(typing, "Annotated", None)
TYPING_LITERAL = getattr(typing, "Literal", None)
from unittest import mock

try:  # Python 3.11+
    import tomllib  # type: ignore[attr-defined]
except ModuleNotFoundError:  # Python 3.10 fallback
    import tomli as tomllib  # type: ignore[no-redef]

from jsonargparse import ArgumentParser as JsonArgParser  # type: ignore


class ParserCaptured(Exception):
    """Raised to short-circuit CLI execution once the parser is available."""


captured_parser: Optional[JsonArgParser] = None


class CaptureArgumentParser(JsonArgParser):
    """ArgumentParser subclass that captures the constructed parser."""

    def parse_args(self, *args: Any, **kwargs: Any):  # type: ignore[override]
        global captured_parser
        captured_parser = self
        raise ParserCaptured


def parse_pyproject(path: Path) -> Dict[str, str]:
    data = tomllib.loads(path.read_text(encoding="utf-8"))
    scripts = data.get("project", {}).get("scripts", {})
    if not isinstance(scripts, dict):
        raise ValueError("[project.scripts] section is missing or invalid")
    return {name: target for name, target in scripts.items() if isinstance(target, str)}


def capture_parser(cli_module_name: str) -> tuple[JsonArgParser, ModuleType]:
    """Import a CLI module and capture its ArgumentParser instance."""

    global captured_parser
    captured_parser = None

    with mock.patch("jsonargparse.ArgumentParser", CaptureArgumentParser):
        cli_module = importlib.import_module(cli_module_name)
        cli_module = importlib.reload(cli_module)
        runner = getattr(cli_module, "main", None)
        if runner is None:
            raise RuntimeError(f"Module '{cli_module_name}' does not define a main() function")
        try:
            runner([])
        except ParserCaptured:
            pass

    parser = captured_parser
    captured_parser = None

    if parser is None:
        raise RuntimeError(f"Failed to capture parser for '{cli_module_name}'")

    return parser, cli_module


def friendly_tool_name(tool_id: str) -> str:
    return tool_id.replace(".py", "").replace("_", " ").title()


def _slugify(value: Any) -> str:
    if isinstance(value, (list, tuple)):
        text = "_".join(str(v) for v in value)
    else:
        text = str(value)
    text = re.sub(r"[^0-9a-zA-Z]+", "_", text).strip("_")
    return text.lower() or "option"


def _resolve_type_hint(annotation: Any) -> Optional[type]:
    if annotation is None or annotation is Any:
        return None
    if isinstance(annotation, type):
        return annotation

    origin = get_origin(annotation)
    if origin is None:
        return annotation if isinstance(annotation, type) else None

    if origin in (list, tuple, set, frozenset):
        args = get_args(annotation)
        if args:
            return _resolve_type_hint(args[0])
        return None

    if origin is dict:
        args = get_args(annotation)
        if len(args) == 2:
            return _resolve_type_hint(args[1])
        return None

    if origin is typing.Union or (TypesUnionType is not None and origin is TypesUnionType):
        args = [arg for arg in get_args(annotation) if arg is not type(None)]  # noqa: E721
        if len(args) == 1:
            return _resolve_type_hint(args[0])
        return None

    if TYPING_ANNOTATED is not None and origin is TYPING_ANNOTATED:
        args = get_args(annotation)
        if args:
            return _resolve_type_hint(args[0])
        return None

    if TYPING_LITERAL is not None and origin is TYPING_LITERAL:
        args = get_args(annotation)
        if args:
            return type(args[0])
        return None

    return origin if isinstance(origin, type) else None


def _infer_action_base_type(action) -> Optional[type]:
    if isinstance(action, (_StoreTrueAction, _StoreFalseAction)):
        return bool

    hint = _resolve_type_hint(getattr(action, "_typehint", None))
    if hint is None:
        hint = _resolve_type_hint(getattr(action, "type", None))

    return hint if isinstance(hint, type) else None


def guess_parameter_type(action) -> str:
    dest = action.dest or ""
    dest_lower = dest.lower()
    option_lower = " ".join(action.option_strings or ()).lower()

    base_type = _infer_action_base_type(action)

    if base_type is bool:
        return "boolean"
    if base_type is not None:
        # bool is a subclass of int, so check explicitly before integers
        if issubclass(base_type, bool):
            return "boolean"
        if issubclass(base_type, int):
            return "integer"
        if issubclass(base_type, float):
            return "number"
        if issubclass(base_type, Path):
            return "file"
    if "output" in dest_lower or dest_lower.endswith("_out"):
        return "output"
    if any(keyword in dest_lower for keyword in ("file", "path", "dir", "input", "database")):
        return "file"
    if any(keyword in option_lower for keyword in ("--output", "--out", "--sparse", "--dense", "--dense_text", "--html", "--xgmml")):
        return "output"
    if any(keyword in option_lower for keyword in ("--file", "--path", "--dir", "--input")):
        return "file"
    return "string"


def friendly_parameter_name(action) -> str:
    if action.option_strings:
        long_opt = next((opt for opt in action.option_strings if opt.startswith("--")), action.option_strings[0])
        base = long_opt.lstrip("-")
    else:
        base = action.dest or "parameter"
    return base.replace("_", " ").strip().title()


def coerce_json(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, (list, tuple)):
        return [coerce_json(v) for v in value]
    if isinstance(value, dict):
        return {str(k): coerce_json(v) for k, v in value.items()}
    return str(value)


def action_to_parameter(action) -> Optional[Dict[str, Any]]:
    if action.dest in {"help", "config", "print_config"}:
        return None
    if any(opt in {"-h", "--help"} for opt in action.option_strings or () ):
        return None

    # jsonargparse appends special flags multiple times; generate unique parameter ids when needed
    def _base_parameter_dict(name: str, parameter: str, help_text: str, param_type: str) -> Dict[str, Any]:
        data: Dict[str, Any] = {
            "name": name,
            "parameter": parameter,
            "help": help_text,
            "type": param_type,
        }
        if getattr(action, "option_strings", None):
            data["flags"] = list(action.option_strings)
        return data

    if isinstance(action, _AppendConstAction):
        target = action.dest or friendly_parameter_name(action)
        const_value = coerce_json(action.const)
        parameter_id = f"{target}__{_slugify(const_value)}"
        help_text = (action.help or "").strip() or "Include this column in the report."
        parameter = _base_parameter_dict(friendly_parameter_name(action), parameter_id, help_text, "boolean")
        parameter["behavior"] = "append_const"
        parameter["append_const"] = {"target": action.dest, "value": const_value}
        parameter["target_parameter"] = action.dest
        if action.default is not None and action.default is not SUPPRESS:
            parameter["default"] = bool(action.default)
        return parameter

    if action.__class__.__name__.lower() == "dynamicarg":
        target = action.dest or friendly_parameter_name(action)
        parameter_id = f"{target}__{_slugify(action.const)}"
        help_text = (action.help or "").strip()
        guidance = "Enter one entry per line; separate multiple values with spaces or commas."
        help_combined = f"{help_text} {guidance}" if help_text else guidance
        parameter = _base_parameter_dict(friendly_parameter_name(action), parameter_id, help_combined, "dynamic")
        parameter["behavior"] = "dynamic"
        parameter["dynamic_action"] = {
            "target": action.dest,
            "const": coerce_json(action.const),
            "nargs": action.nargs,
        }
        parameter["target_parameter"] = action.dest
        return parameter

    parameter: Dict[str, Any] = {
        "name": friendly_parameter_name(action),
        "parameter": action.dest,
        "help": (action.help or "").strip() or "No description provided.",
        "type": guess_parameter_type(action),
    }

    if action.required:
        parameter["required"] = True

    multiple = False
    if action.nargs in ("+", "*"):
        multiple = True
    elif isinstance(action.nargs, int):
        multiple = action.nargs > 1
    if multiple:
        parameter["multiple"] = True

    if action.choices is not None:
        parameter["choices"] = list(action.choices)

    if action.default is not None and action.default is not SUPPRESS:
        parameter["default"] = coerce_json(action.default)

    if getattr(action, "option_strings", None):
        parameter["flags"] = list(action.option_strings)

    return parameter


def module_description(module: ModuleType) -> str:
    doc = module.__doc__ or ""
    doc = doc.strip()
    return doc


def resolve_cli_module(target: str) -> tuple[str, str]:
    module_name, sep, attr_path = target.partition(":")
    if not sep:
        raise ValueError(f"Target '{target}' does not contain ':'")
    parts = [segment for segment in attr_path.split(".") if segment]
    if not parts:
        raise ValueError(f"Target '{target}' has no attribute path")
    if len(parts) == 1:
        cli_module = f"{module_name}.{parts[0]}"
    else:
        cli_module = ".".join([module_name] + parts[:-1])
    return module_name, cli_module


def build_schema(tool_id: str, target: str) -> Dict[str, Any]:
    module_name, cli_module_name = resolve_cli_module(target)
    parser, cli_module = capture_parser(cli_module_name)

    parameters: List[Dict[str, Any]] = []
    for action in parser._actions:
        param = action_to_parameter(action)
        if param:
            parameters.append(param)

    schema: Dict[str, Any] = {
        "id": tool_id.replace(".py", ""),
        "display_name": friendly_tool_name(tool_id),
        "category": "Domainator",
        "description": module_description(cli_module),
        "runner": "module",
        "entry_point": cli_module_name,
        "parameters": parameters,
        "advanced_parameters": [],
    }
    return schema


def write_schema(schema: Dict[str, Any], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    tool_id = schema["id"]
    output_path = output_dir / f"{tool_id}.json"
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(schema, handle, indent=2)
        handle.write("\n")


def main(argv: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(description="Generate tool schema files from Domainator CLIs")
    parser.add_argument("--pyproject", type=Path, default=Path("pyproject.toml"), help="Path to the Domainator pyproject.toml")
    parser.add_argument("--output-dir", type=Path, required=True, help="Destination directory for generated schemas")
    parser.add_argument("--include", nargs="*", help="Optional subset of script names to generate (match project.scripts keys)")
    args = parser.parse_args(argv)

    scripts = parse_pyproject(args.pyproject)
    if args.include:
        missing = [name for name in args.include if name not in scripts]
        if missing:
            raise SystemExit(f"Requested script(s) not found: {', '.join(missing)}")
        scripts = {name: scripts[name] for name in args.include}

    for tool_id, target in sorted(scripts.items()):
        schema = build_schema(tool_id, target)
        write_schema(schema, args.output_dir)
        print(f"Generated schema for {tool_id}")


if __name__ == "__main__":
    main()
