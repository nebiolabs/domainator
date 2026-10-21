"""Tests for the tool manifest (domainator describe).

The manifest is derived from the parsers; these tests keep the parsers describable: every
option has a role, output options follow the naming rule, and a GUI's config file (keys are the
parameter names) parses back to the same values.
"""
import importlib
import json
from pathlib import Path

import pytest

from domainator.cli import manifest as manifest_module
from domainator.cli.registry import TOOLS

SNAPSHOT = Path(__file__).parent / "data" / "manifest_snapshot.json"


@pytest.fixture(scope="module")
def manifest():
    return manifest_module.manifest()


def _interface(manifest):
    """The parts of the manifest that are the command-line interface (not the help text)."""
    keys = ("name", "flags", "negative_flags", "kind", "type", "role", "multiple", "count", "required", "choices")
    return {tool["name"]: [{k: p[k] for k in keys if k in p} for p in tool["parameters"]] for tool in manifest["tools"]}


def test_every_command_is_described_and_serializes(manifest):
    assert [tool["name"] for tool in manifest["tools"]] == list(TOOLS)
    assert json.loads(json.dumps(manifest)) == manifest
    assert manifest_module.manifest() == manifest  # deterministic


def test_string_options_are_classified(manifest):
    used = set()
    for tool in manifest["tools"]:
        for parameter in tool["parameters"]:
            if parameter["kind"] == "value" and parameter["type"] == "string" and "choices" not in parameter:
                used.add(parameter["name"])
                assert parameter["name"] in manifest_module.TEXT_OPTIONS, (
                    f"{tool['name']} --{parameter['name']}: say in cli/manifest.py whether it is a path "
                    "(INPUT_PATHS, DIRECTORIES, or a name ending in _out, _file, _dir or _path) or text (TEXT_OPTIONS)")
    assert not manifest_module.TEXT_OPTIONS - used, "unused TEXT_OPTIONS entries"


def test_output_naming_rule(manifest):
    # An option ends in _out exactly when it names an output.
    for tool in manifest["tools"]:
        names = {p["name"]: p for p in tool["parameters"]}
        for parameter in tool["parameters"]:
            assert (parameter.get("role") == "output") == parameter["name"].endswith("_out"), (tool["name"], parameter["name"])
        for dest in tool["required_outputs"] or []:
            assert names[dest]["role"] == "output", (tool["name"], dest)


def test_parsers_reject_abbreviations():
    for tool in TOOLS.values():
        assert importlib.import_module(f"domainator.{tool.module}").build_parser().allow_abbrev is False


@pytest.mark.parametrize("tool_name", list(TOOLS))
def test_config_with_manifest_defaults_parses(manifest, tool_name, tmp_path):
    # What a GUI does: write every parameter's value to a config file keyed by name.
    tool = next(t for t in manifest["tools"] if t["name"] == tool_name)
    config = {}
    for parameter in tool["parameters"]:
        if parameter.get("default_dynamic"):
            continue
        value = parameter.get("default")
        if parameter.get("required") and value is None:
            value = parameter["choices"][0] if "choices" in parameter else {"integer": 1, "number": 1.0}.get(parameter["type"], "x")
            value = [value] if parameter.get("multiple") else value
        config[parameter["name"]] = value
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    parser = importlib.import_module(f"domainator.{TOOLS[tool_name].module}").build_parser()
    params = parser.parse_args(["--config", str(path)])
    for name, value in config.items():
        assert getattr(params, name) == value, (tool_name, name)


def test_interface_snapshot(manifest):
    """The command-line interface is part of the 1.x stability contract (docs/stability.md).

    If you changed an option on purpose, regenerate the snapshot:
        python -c "import json; from test.test_manifest import _interface; from domainator.cli.manifest import manifest; print(json.dumps(_interface(manifest()), indent=1))" > test/data/manifest_snapshot.json
    """
    assert _interface(manifest) == json.loads(SNAPSHOT.read_text())


def test_describe_command(capsys):
    from domainator import cli
    assert cli.main(["describe", "--tool", "sort_contigs"]) == 0
    described = json.loads(capsys.readouterr().out)
    assert described["manifest_version"] == manifest_module.MANIFEST_VERSION
    (tool,) = described["tools"]
    parameters = {p["name"]: p for p in tool["parameters"]}
    assert tool["required_outputs"] == ["genbank_out", "fasta_out"]
    assert parameters["input"]["role"] == "input" and parameters["input"]["stdin"] is True
    assert parameters["genbank_out"]["flags"] == ["-o", "--genbank_out"]
    assert parameters["max_hits"]["type"] == "integer"
