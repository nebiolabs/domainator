#!/usr/bin/env python3
"""Generate Domainator tool schema JSON files from CLI definitions."""
from __future__ import annotations

import argparse
import importlib
import json
import sys
from argparse import SUPPRESS, _StoreFalseAction, _StoreTrueAction
from pathlib import Path
from types import ModuleType
from typing import Any, Dict, List, Optional
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


def guess_parameter_type(action) -> str:
    dest = action.dest or ""
    dest_lower = dest.lower()
    option_lower = " ".join(action.option_strings or ()).lower()

    if isinstance(action, (_StoreTrueAction, _StoreFalseAction)):
        return "boolean"
    if action.type is int:
        return "integer"
    if action.type is float:
        return "number"
    if "output" in dest_lower or dest_lower.endswith("_out"):
        return "output"
    if any(keyword in dest_lower for keyword in ("file", "path", "dir", "input", "database")):
        return "file"
    if any(keyword in option_lower for keyword in ("--output", "--out")):
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
    if action.dest in {"help", "config"}:
        return None
    if any(opt in {"-h", "--help"} for opt in action.option_strings or () ):
        return None

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
    elif isinstance(action.nargs, int) and action.nargs != 1:
        multiple = True
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
