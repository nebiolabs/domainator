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


@pytest.mark.parametrize("command", sorted(TOOLS))
def test_evalue_options_share_one_validator(command):
    """-e/--evalue must be > 0 and --min_evalue >= 0 everywhere, checked at parse time."""
    from jsonargparse.typing import NonNegativeFloat, PositiveFloat

    parser = importlib.import_module(f"domainator.{TOOLS[command].module}").build_parser()
    for action in parser._actions:
        if "--evalue" in action.option_strings:
            assert getattr(action, "_typehint", action.type) is PositiveFloat
        if "--min_evalue" in action.option_strings:
            assert getattr(action, "_typehint", action.type) is NonNegativeFloat


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


def test_domainator_errors_are_one_line_messages(shared_datadir, tmp_path, capsys):
    """A DomainatorError is printed without a traceback and exits with status 1."""
    bad = tmp_path / "bad.hdf5"
    import h5py
    with h5py.File(bad, "w") as f:
        f.attrs["MATRIX_FILE_VERSION"] = "99.0"
    assert cli.main(["build_tree", "-i", str(bad), "-o", str(tmp_path / "out.nwk")]) == 1
    err = capsys.readouterr().err
    assert err.startswith("domainator build_tree: error: Unsupported matrix file version '99.0'")
    assert "Traceback" not in err


def test_verbose_keeps_the_traceback(tmp_path):
    from domainator.errors import FormatVersionError
    bad = tmp_path / "bad.hdf5"
    import h5py
    with h5py.File(bad, "w") as f:
        f.attrs["MATRIX_FILE_VERSION"] = "99.0"
    with pytest.raises(FormatVersionError):
        cli.main(["build_tree", "-v", "-i", str(bad), "-o", str(tmp_path / "out.nwk")])


def test_quiet_and_verbose_set_the_log_level(shared_datadir, tmp_path, capsys):
    import logging
    out = str(tmp_path / "out.gb")
    args = ["-i", str(shared_datadir / "pDONR201_multi_genemark_domainator.gb"), "-o", out, "--contigs", "pDONR201_1", "--cds_range", "3"]
    assert cli.main(["trim_contigs", *args]) == 0
    assert logging.getLogger("domainator").level == logging.INFO
    assert cli.main(["trim_contigs", "-q", *args]) == 0
    assert logging.getLogger("domainator").level == logging.ERROR
    assert cli.main(["trim_contigs", "--verbose", *args]) == 0
    assert logging.getLogger("domainator").level == logging.DEBUG
    with pytest.raises(SystemExit):
        cli.main(["trim_contigs", "-q", "-v", *args])
    assert "not allowed with" in capsys.readouterr().err


def test_warnings_are_logged_with_a_prefix(tmp_path, capsys):
    from domainator.Bio.Seq import Seq
    from domainator.Bio.SeqRecord import SeqRecord
    from domainator.Bio.SeqFeature import FeatureLocation, SeqFeature
    from domainator.utils import write_genbank
    # Every CDS is unannotated, so --unannotated would strip them all: trim_contigs warns and skips the contig.
    record = SeqRecord(Seq("ATG" * 200), id="rec", annotations={"molecule_type": "DNA"})
    record.features = [SeqFeature(FeatureLocation(0, 300, 1), type="CDS"), SeqFeature(FeatureLocation(400, 600, 1), type="CDS")]
    in_gb = str(tmp_path / "in.gb")
    write_genbank([record], in_gb)
    args = ["trim_contigs", "-i", in_gb, "-o", str(tmp_path / "out.gb"), "--unannotated"]
    assert cli.main(args) == 0
    assert "warning: All CDSs in contig rec are composed of domains to be stripped" in capsys.readouterr().err
    assert cli.main(args[:1] + ["-q"] + args[1:]) == 0
    assert capsys.readouterr().err == ""
