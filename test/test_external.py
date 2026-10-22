"""Tests for external.py: finding and running external programs, and optional imports."""
import os
import sys

import pytest

from domainator.errors import ExternalToolError
from domainator.external import append_extra_args, import_optional, resolve_binary, run_command


def test_resolve_binary_uses_an_explicit_path():
    assert resolve_binary("python", sys.executable) == sys.executable


def test_resolve_binary_rejects_a_path_that_is_not_executable(tmp_path):
    not_executable = tmp_path / "diamond"
    not_executable.write_text("")
    with pytest.raises(ExternalToolError, match="is not an executable file"):
        resolve_binary("diamond", str(not_executable))


def test_resolve_binary_searches_the_path(monkeypatch):
    monkeypatch.setenv("PATH", os.path.dirname(sys.executable))
    assert resolve_binary(os.path.basename(sys.executable)) == sys.executable


def test_resolve_binary_names_the_option_and_install_hint(monkeypatch, tmp_path):
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(ExternalToolError, match=r"Could not find 'cd-hit' on PATH\. Install it\. .*--cdhit_path"):
        resolve_binary("cd-hit", option="--cdhit_path", install_hint="Install it.")
    with pytest.raises(ExternalToolError, match="--diamond_path"):
        resolve_binary("diamond")


def test_run_command_success():
    assert run_command([sys.executable, "-c", "pass"], "noop") is None
    completed = run_command([sys.executable, "-c", "print('hi')"], "echo", capture=True)
    assert completed.stdout.decode().strip() == "hi"


def test_run_command_failure_reports_the_end_of_stderr():
    script = "import sys\nfor i in range(40): print(f'line {i}', file=sys.stderr)\nsys.exit(3)"
    with pytest.raises(ExternalToolError) as error:
        run_command([sys.executable, "-c", script], "the program")
    message = str(error.value)
    assert message.startswith("the program failed with exit code 3.")
    assert "line 39" in message and "line 15" in message and "line 14\n" not in message


def test_run_command_failure_falls_back_to_stdout():
    with pytest.raises(ExternalToolError, match="only on stdout"):
        run_command([sys.executable, "-c", "print('only on stdout'); raise SystemExit(1)"], "the program")


def test_run_command_accepts_a_failure_that_produced_the_expected_output(tmp_path):
    expected = tmp_path / "result"
    expected.write_text("done")
    assert run_command([sys.executable, "-c", "raise SystemExit(1)"], "foldseek", expect=str(expected)) is None
    with pytest.raises(ExternalToolError):
        run_command([sys.executable, "-c", "raise SystemExit(1)"], "foldseek", expect=str(tmp_path / "missing"))


def test_run_command_missing_program(tmp_path):
    with pytest.raises(ExternalToolError, match="could not run"):
        run_command([str(tmp_path / "no_such_program")], "the program")


def test_append_extra_args():
    assert append_extra_args(["diamond", "-p", "4"], ["--sensitive", "-k", "5"], "--diamond_args") == \
        ["diamond", "-p", "4", "--sensitive", "-k", "5"]
    assert append_extra_args(["diamond"], None, "--diamond_args") == ["diamond"]
    with pytest.raises(ValueError, match="--diamond_args cannot set -k, -p: Domainator sets them itself"):
        append_extra_args(["diamond", "-p", "4", "-k", "1"], ["-p", "2", "-k", "3"], "--diamond_args")


def test_import_optional():
    assert import_optional("json", "unused").dumps(1) == "1"
    with pytest.raises(ExternalToolError, match="optional package 'no_such_package'.*pip install it"):
        import_optional("no_such_package.sub", "pip install it")
