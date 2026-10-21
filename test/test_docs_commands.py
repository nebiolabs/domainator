"""Every option shown with a Domainator command in the documentation must exist on that command.

Options must be spelled out in full: abbreviations are rejected (allow_abbrev=False). Complete
commands in code blocks must also name an output, since no command writes to stdout by default.
"""

import contextlib
import importlib
import io
import re
import shlex
from pathlib import Path

import pytest

from domainator.cli.registry import TOOLS

REPO = Path(__file__).resolve().parents[1]
DOCS = [REPO / "README.md", REPO / "AGENTS.md", *sorted((REPO / "docs").rglob("*.md"))]
# Documents that discuss options as history or plans rather than as current usage.
EXCLUDED_DOCS = {REPO / "docs" / "v1_release_plan.md"}

_MODULE_TO_COMMAND = {tool.module: name for name, tool in TOOLS.items()}
_COMMAND_RE = re.compile(r"(?:^|[\s|;(])(?:domainator\s+([a-z_]+)|([a-z_]+)\.py)(?=\s)")


def _code_snippets(text):
    """Yield (line_number, code, fenced) for each line of a fenced code block (shell comments removed) and each inline code span."""
    in_fence = False
    for line_number, line in enumerate(text.split("\n"), start=1):
        if line.lstrip().startswith("```"):
            in_fence = not in_fence
            continue
        if in_fence:
            yield line_number, re.sub(r"(^|\s)#.*", "", line), True
        else:
            for span in re.findall(r"`([^`]+)`", line):
                yield line_number, span, False


def _documented_invocations():
    for doc in DOCS:
        if doc in EXCLUDED_DOCS:
            continue
        text = doc.read_text(encoding="utf-8")
        text = re.sub(r"\\\n\s*", " ", text)  # join shell line continuations
        for line_number, code, fenced in _code_snippets(text):
            for match in _COMMAND_RE.finditer(" " + code):
                command = match.group(1) or _MODULE_TO_COMMAND.get(match.group(2), match.group(2))
                if command not in TOOLS:
                    continue
                rest = (" " + code)[match.end():]
                rest = re.split(r"\s[|;>]\s|\s&&\s", rest)[0]  # stop at the end of this command
                try:
                    tokens = shlex.split(rest)
                except ValueError:
                    tokens = rest.split()
                options = [t.split("=", 1)[0] for t in tokens if re.fullmatch(r"--?[A-Za-z][\w-]*(=.*)?", t)]
                if options:
                    yield pytest.param(command, options, tokens if fenced else None,
                                       id=f"{doc.relative_to(REPO)}:{line_number}:{command}")


_PARSERS = {}


def _parser(command):
    if command not in _PARSERS:
        _PARSERS[command] = importlib.import_module(f"domainator.{TOOLS[command].module}").build_parser()
    return _PARSERS[command]


def _option_strings(command):
    return {s for action in _parser(command)._actions for s in action.option_strings}


@pytest.mark.parametrize("command,options,tokens", list(_documented_invocations()))
def test_documented_options_exist(command, options, tokens):
    unknown = [option for option in options if option not in _option_strings(command)]
    assert not unknown, f"{command} has no option(s) {unknown}"

    if tokens is None or "--help" in tokens:
        return
    parser = _parser(command)
    try:
        with contextlib.redirect_stderr(io.StringIO()):
            params = parser.parse_args(tokens)
    except SystemExit:  # placeholders such as "..." or "<threshold>"
        return
    outputs = [action.dest for action in parser._actions if hasattr(action, "stdout")]
    if outputs:
        assert any(getattr(params, dest, None) is not None for dest in outputs), \
            f"{command} names none of its outputs ({', '.join(outputs)}); use - to write to stdout"
