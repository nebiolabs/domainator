"""The ``domainator`` dispatcher, its command registry, and shell completion."""

import importlib
import subprocess
import sys

import pytest

from domainator import __version__, cli
from domainator.cli.registry import CATEGORIES, TOOLS


@pytest.mark.parametrize("command", sorted(TOOLS))
def test_every_command_builds_its_parser(command):
    module = importlib.import_module(f"domainator.{TOOLS[command].module}")
    assert callable(module.main)
    parser = module.build_parser()
    assert parser.allow_abbrev is False
    assert parser.description


def test_registry_categories_are_known():
    keys = {key for key, _ in CATEGORIES}
    assert {tool.category for tool in TOOLS.values()} <= keys
    assert not set(TOOLS) & set(cli.BUILTIN_COMMANDS)


def test_help_lists_every_command(capsys):
    assert cli.main(["--help"]) == 0
    out = capsys.readouterr().out
    for command in list(TOOLS) + list(cli.BUILTIN_COMMANDS):
        assert f"  {command} " in out


def test_no_arguments_prints_help(capsys):
    assert cli.main([]) == 0
    assert "usage: domainator <command>" in capsys.readouterr().out


def test_version(capsys):
    assert cli.main(["--version"]) == 0
    assert capsys.readouterr().out.strip() == f"domainator {__version__}"


def test_unknown_command_suggests_close_matches(capsys):
    assert cli.main(["domainate.py"]) == 2
    err = capsys.readouterr().err
    assert "unknown command 'domainate.py'" in err
    assert "domainate" in err


def test_tool_usage_names_the_subcommand(capsys):
    original_prog = sys.argv[0]
    with pytest.raises(SystemExit) as exc:
        cli.main(["genbank_to_fasta", "--help"])
    assert exc.value.code == 0
    assert "usage: domainator genbank_to_fasta" in capsys.readouterr().out
    assert sys.argv[0] == original_prog


def test_abbreviated_options_are_rejected(capsys):
    with pytest.raises(SystemExit) as exc:
        cli.main(["enum_report", "--name"])
    assert exc.value.code == 2
    assert "unrecognized arguments: --name" in capsys.readouterr().err


def test_dispatch_runs_the_tool(shared_datadir, tmp_path):
    out = tmp_path / "out.fasta"
    assert cli.main(["genbank_to_fasta", "-i", str(shared_datadir / "pDONR201.gb"), "-o", str(out)]) == 0
    assert out.read_text().startswith(">")


def test_python_dash_m_runs_the_dispatcher():
    result = subprocess.run([sys.executable, "-m", "domainator", "--version"], capture_output=True, text=True)
    assert result.returncode == 0
    assert result.stdout.strip() == f"domainator {__version__}"


@pytest.mark.parametrize("shell", ["bash", "zsh", "tcsh"])
def test_completion_script_is_deterministic(shell):
    pytest.importorskip("shtab")
    from domainator.cli.completion import completion_script

    script = completion_script(shell)
    assert script == completion_script(shell)
    for command in TOOLS:
        assert command in script


def test_bash_completion_completes_commands_and_choices(tmp_path):
    pytest.importorskip("shtab")
    from domainator.cli.completion import completion_script

    script = tmp_path / "domainator.bash"
    script.write_text(completion_script("bash"))
    driver = f"""
source {script}
f=$(complete -p domainator | grep -o "_shtab[a-z_]*")
try() {{ COMP_WORDS=("$@"); COMP_CWORD=$((${{#COMP_WORDS[@]}}-1)); COMP_LINE="${{COMP_WORDS[*]}}"; COMP_POINT=${{#COMP_LINE}}; COMPREPLY=(); $f; echo "${{COMPREPLY[*]}}"; }}
try domainator domain_s
try domainator seq_dist --symmetrize ""
"""
    result = subprocess.run(["bash", "-c", driver], capture_output=True, text=True, check=True)
    commands, choices = result.stdout.strip().split("\n")
    assert commands.split() == ["domain_search"]
    assert choices.split() == ["max", "mean", "min", "none"]
