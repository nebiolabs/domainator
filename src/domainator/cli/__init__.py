"""The ``domainator`` command: ``domainator <command> [options]``.

Each command is a tool module listed in :mod:`domainator.cli.registry`. Only the module
of the requested command is imported, so ``domainator --help`` and dispatching stay fast
even though some tools import heavy dependencies.
"""

import difflib
import importlib
import logging
import sys
from typing import List, Optional

from domainator import __version__
from domainator.cli.registry import CATEGORIES, TOOLS
from domainator.errors import DomainatorError

BUILTIN_COMMANDS = {
    "completion": "Print a shell completion script: domainator completion {bash,zsh,tcsh}.",
    "doctor": "List the optional programs and packages, and which of them are installed.",
}
"""Commands of the dispatcher itself, rather than tools."""


def _help_text() -> str:
    width = max(len(name) for name in list(TOOLS) + list(BUILTIN_COMMANDS))
    lines = [
        "usage: domainator <command> [options]",
        "       domainator <command> --help",
        "",
        "Domain-based gene neighborhood and protein search, extraction, and clustering.",
    ]
    for key, heading in CATEGORIES:
        lines += ["", f"{heading}:"]
        lines += [f"  {name:<{width}}  {tool.summary}" for name, tool in TOOLS.items() if tool.category == key]
    lines += ["", "Other commands:"]
    lines += [f"  {name:<{width}}  {summary}" for name, summary in BUILTIN_COMMANDS.items()]
    lines += ["", "options:", "  -h, --help     show this message and exit", "  --version      print the Domainator version and exit",
              "", "Every command also takes -v/--verbose (debug output, and tracebacks for errors) and",
              "-q/--quiet (errors only).",
              "", f"Domainator {__version__}"]
    return "\n".join(lines)


def run_tool(command: str, argv: List[str]) -> int:
    """Run one registered command with an explicit argument list and return its exit code.

    A :class:`~domainator.errors.DomainatorError` is reported as a one-line message with exit
    code 1, unless -v/--verbose asked for debug output, which keeps the traceback.
    """
    module = importlib.import_module(f"domainator.{TOOLS[command].module}")
    original_prog = sys.argv[0] if sys.argv else None
    # argparse names the program in usage and error messages after sys.argv[0].
    sys.argv[0:1] = [f"domainator {command}"]
    try:
        status = module.main(argv)
    except DomainatorError as exc:
        if logging.getLogger("domainator").isEnabledFor(logging.DEBUG):
            raise
        print(f"domainator {command}: error: {exc}", file=sys.stderr)
        return 1
    finally:
        sys.argv[0:1] = [] if original_prog is None else [original_prog]
    return 0 if status is None else status


def main(argv: Optional[List[str]] = None) -> int:
    argv = sys.argv[1:] if argv is None else list(argv)
    if not argv or argv[0] in ("-h", "--help"):
        print(_help_text())
        return 0
    if argv[0] == "--version":
        print(f"domainator {__version__}")
        return 0

    command, args = argv[0], argv[1:]
    if command == "completion":
        from domainator.cli.completion import main as completion_main
        return completion_main(args)
    if command == "doctor":
        from domainator.cli.doctor import main as doctor_main
        return doctor_main(args)
    if command not in TOOLS:
        message = f"domainator: unknown command '{command}'."
        close = difflib.get_close_matches(command.removesuffix(".py"), list(TOOLS) + list(BUILTIN_COMMANDS), n=3)
        if close:
            message += " Did you mean: " + ", ".join(close) + "?"
        print(message + " Run 'domainator --help' for the list of commands.", file=sys.stderr)
        return 2
    return run_tool(command, args)
