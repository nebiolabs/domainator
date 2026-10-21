"""``domainator completion {bash,zsh,tcsh}``: print a static shell completion script.

The script is generated from every command's parser, so completing is instant at the
shell (no Python runs on <Tab>). Regenerate it after upgrading Domainator.
"""

import argparse
import importlib
import sys
from typing import List

from domainator import __version__
from domainator.cli.common import tool_parser
from domainator.cli.registry import TOOLS

SHELLS = ("bash", "zsh", "tcsh")


def build_completion_parser() -> argparse.ArgumentParser:
    """A parser for the whole ``domainator`` command, with one subcommand per tool."""
    import shtab

    root = tool_parser("Domainator")
    root.prog = "domainator"
    root.add_argument("--version", action="version", version=f"domainator {__version__}")
    commands = root.add_subcommands(required=True, dest="command")
    completion = tool_parser(__doc__)
    completion.add_argument("shell", choices=SHELLS)
    commands.add_subcommand("completion", completion, help="Print a shell completion script.")
    commands.add_subcommand("doctor", tool_parser("List the optional programs and packages."), help="List the optional programs and packages.")
    for name, tool in TOOLS.items():
        parser = importlib.import_module(f"domainator.{tool.module}").build_parser()
        for action in parser._actions:
            if action.choices is not None:
                # Sets would make the script differ from run to run.
                action.choices = sorted(action.choices, key=str)
            elif action.option_strings and action.nargs != 0 and not isinstance(action, argparse._HelpAction):
                action.complete = shtab.FILE
        commands.add_subcommand(name, parser, help=tool.summary)
    return root


def completion_script(shell: str) -> str:
    import shtab

    return shtab.complete(build_completion_parser(), shell)


def main(argv: List[str]) -> int:
    parser = tool_parser(__doc__ + """
To install for bash:
    domainator completion bash > ~/.local/share/bash-completion/completions/domainator
For zsh, write it to a directory on your $fpath, named _domainator:
    domainator completion zsh > ~/.zfunc/_domainator
""")
    parser.prog = "domainator completion"
    parser.add_argument("shell", choices=SHELLS)
    params = parser.parse_args(argv)
    try:
        import shtab  # noqa: F401
    except ImportError:
        print("domainator completion needs the 'shtab' package: pip install shtab", file=sys.stderr)
        return 1
    sys.stdout.write(completion_script(params.shell))
    return 0
