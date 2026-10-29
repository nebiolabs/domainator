"""Parser construction and option conventions shared by every Domainator command.

Conventions (see MIGRATION.md):

- ``-i/--input`` is the file being edited or read. Standard input is read only when asked for,
  with ``-i -``, and only by commands that read their input once.
- Every output file has its own option ending in ``_out`` (``--genbank_out``, ``--fasta_out``,
  ``--tsv_out``, ...), and the main one also has ``-o``. Standard output is written only when
  asked for, with ``-``; at most one output can go there.
- Every command takes ``-v/--verbose`` and ``-q/--quiet``. Progress messages go to the
  ``domainator`` logger, which a command line run shows on stderr at INFO level.
"""

import argparse
import inspect
import logging
import os
import sys
import warnings
from typing import Iterable, Optional, Sequence

from jsonargparse import ArgumentParser

from domainator import __version__
from domainator.hmmio import iter_hmm_names, read_hmms
from domainator.model import ANNOTATION_TYPES
from domainator.seqio import STDIN_PATHS, STDOUT_PATHS, get_file_type, parse_seqfiles


LOGGER_NAME = "domainator"


class RawAndDefaultsFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    """Keeps the line breaks of descriptions and help strings, and appends each option's default.

    Empty defaults (None, False, an empty list) are left out: the help text says what leaving
    the option out means.
    """

    def _get_help_string(self, action):
        default = action.default
        if default is None or default is False or (isinstance(default, (list, tuple, str)) and len(default) == 0):
            return action.help
        return super()._get_help_string(action)


class _StderrHandler(logging.StreamHandler):
    """Writes to whatever ``sys.stderr`` is when a message is logged, so redirection and capture work."""

    def __init__(self):
        super().__init__(sys.stderr)

    @property
    def stream(self):
        return sys.stderr

    @stream.setter
    def stream(self, value):
        pass


class _MessageFormatter(logging.Formatter):
    """Plain messages for progress; ``warning: ...`` and ``error: ...`` for problems."""

    def format(self, record):
        message = super().format(record)
        if record.levelno >= logging.WARNING:
            return f"{record.levelname.lower()}: {message}"
        return message


def configure_logging(verbose: bool = False, quiet: bool = False) -> None:
    """Show the ``domainator`` logger on stderr: INFO by default, DEBUG with verbose, errors only with quiet.

    quiet also hides Python warnings raised by Domainator modules.
    """
    logger = logging.getLogger(LOGGER_NAME)
    if not any(isinstance(handler, _StderrHandler) for handler in logger.handlers):
        handler = _StderrHandler()
        handler.setFormatter(_MessageFormatter("%(message)s"))
        logger.addHandler(handler)
    logger.setLevel(logging.DEBUG if verbose else logging.ERROR if quiet else logging.INFO)
    if quiet:
        warnings.filterwarnings("ignore", module=r"domainator(\.|$)")


class ToolParser(ArgumentParser):
    """ArgumentParser that applies -v/--verbose and -q/--quiet once the arguments are parsed."""

    def parse_args(self, *args, **kwargs):
        params = super().parse_args(*args, **kwargs)
        configure_logging(getattr(params, "verbose", False), getattr(params, "quiet", False))
        return params


def tool_parser(description: str) -> ArgumentParser:
    """Return the ArgumentParser every Domainator command builds its options on.

    Abbreviated options are rejected, so a typo or a removed option is an error rather than a
    silent match to a longer option (``--name`` used to match ``--named_domain_count``).
    The description is dedented, as Python 3.13+ does to docstrings, so help is the same on
    every Python version.
    """
    parser = ToolParser(
        description=inspect.cleandoc(description),
        epilog=f"Domainator {__version__}",
        formatter_class=RawAndDefaultsFormatter,
        allow_abbrev=False,
    )
    verbosity = parser.add_mutually_exclusive_group()
    verbosity.add_argument("-v", "--verbose", action="store_true",
                           help="show debug output, and a traceback for errors.")
    verbosity.add_argument("-q", "--quiet", action="store_true",
                           help="show errors only: no progress messages or warnings.")
    return parser


def add_input(parser, help: str, *, multiple: bool = True, required: bool = True, stdin: bool = True):
    """Add ``-i/--input``. ``stdin=False`` for commands that read their input more than once."""
    if stdin:
        help += " Use - to read standard input."
    action = parser.add_argument("-i", "--input", type=str, default=None, required=required,
                                 nargs="+" if multiple else None, help=help)
    action.stdin = stdin
    return action


def add_output(parser, content: str, help: str, *, primary: bool = False, required: bool = False,
               dest: Optional[str] = None, stdout: bool = True):
    """Add ``--<content>_out`` (and ``-o`` for the command's main output).

    ``stdout=False`` for an output that cannot be written to standard output; ``check_outputs``
    then refuses ``-`` for it.
    """
    flags = (["-o"] if primary else []) + [f"--{content}_out"]
    kwargs = {"dest": dest} if dest else {}
    action = parser.add_argument(*flags, type=str, default=None, required=required,
                                 help=help + (" Use - for standard output." if stdout else ""), **kwargs)
    action.stdout = stdout
    return action


def add_seq_outputs(parser, *, genbank_help: str = "Write the records to this GenBank file.",
                    fasta_help: str = "Write the records to this FASTA file."):
    """Add the ``-o/--genbank_out`` and ``--fasta_out`` outputs of a command that writes sequence records."""
    add_output(parser, "genbank", genbank_help, primary=True)
    add_output(parser, "fasta", fasta_help)


def add_switch(parser, name: str, help: str, *, default: bool, negative_help: str):
    """Add ``--<name>`` and ``--no_<name>``, which set the same boolean ``<name>``."""
    parser.add_argument(f"--{name}", dest=name, action="store_true", default=default, help=help)
    parser.add_argument(f"--no_{name}", dest=name, action="store_false", default=default, help=negative_help)



def add_annotation_type(parser, *, choices: Sequence[str] = ANNOTATION_TYPES, default: str = "domainator", help: str = ""):
    """Add ``--annotation_type``: which annotations a command reads.

    ``domainator`` is the domain annotations written by domainate (``Domainator`` features),
    ``domain_search`` the best hits written by domain_search (``Domain_Search`` features).
    """
    return parser.add_argument("--annotation_type", type=str, default=default, choices=list(choices),
                               help=(help or "Which annotations to use.") + " domainator: domain annotations (from domainate). "
                                    "domain_search: best-hit annotations (from domain_search)."
                                    + (" all: both." if "all" in choices else ""))


def add_has_search_hit(parser):
    """Add the independent, additive search-hit selection criterion."""
    return parser.add_argument("--has_search_hit", action="store_true", default=False,
                               help="Also select existing domain_search hits, regardless of annotation_type, "
                                    "domain names, domain expressions, or databases. Hits must pass the e-value "
                                    "and score cutoffs where applicable. Participates in --invert where available.")


def add_alphabet(parser):
    """Add ``--alphabet {nucleotide,protein}``: the molecule type of FASTA input. GenBank records declare their own."""
    return parser.add_argument("--alphabet", type=str, default="protein", choices=["nucleotide", "protein"],
                               help="Whether the sequences in FASTA inputs are protein or nucleotide. "
                                    "GenBank records declare their own molecule type.")


def declare_outputs(parser, dests: Sequence[str]):
    """Record that the command writes at least one of the outputs ``dests``. check_outputs enforces
    it, and ``domainator describe`` reports it."""
    parser.required_outputs = list(dests)


def check_outputs(parser, params, dests: Optional[Sequence[str]] = None, *, required: bool = True):
    """Require at least one of the outputs ``dests`` (by default the ones declare_outputs recorded),
    and at most one of them on standard output."""
    if dests is None:
        dests = parser.required_outputs
    given = [dest for dest in dests if getattr(params, dest, None) is not None]
    if required and not given:
        options = ", ".join(_option_for(parser, dest) for dest in dests)
        parser.error(f"No output specified. Give at least one of: {options}")
    on_stdout = [dest for dest in given if getattr(params, dest) in STDOUT_PATHS]
    for dest in on_stdout:
        if not getattr(_action_for(parser, dest), "stdout", True):
            parser.error(f"{_option_for(parser, dest)} cannot be written to standard output. Give a file name.")
    if len(on_stdout) > 1:
        parser.error("Only one output can be written to standard output, but got: "
                     + ", ".join(_option_for(parser, dest) for dest in on_stdout))


def reject_stdin(parser, paths: Optional[Iterable[str]], reason: str = "it reads its input more than once"):
    """Refuse ``-i -`` for a command that cannot stream its input."""
    if paths is None:
        return
    if isinstance(paths, str):
        paths = [paths]
    if any(path in STDIN_PATHS for path in paths):
        parser.error(f"This command cannot read standard input, because {reason}. Give a file instead.")


def parse_keyword_args(text: Optional[str], option: str) -> dict:
    """Parse ``"key=value key2=value2"`` (split like a shell command line) into a dict.

    Values are Python literals where they parse as one (``30``, ``0.1``, ``True``, ``'x'``),
    otherwise strings. ``option`` names the command-line option in error messages.
    """
    import ast
    import shlex

    out = {}
    for token in shlex.split(text or ""):
        key, sep, value = token.partition("=")
        if not sep or not key.isidentifier():
            raise ValueError(f"{option} expects key=value pairs, such as {option}=\"n_neighbors=30 min_dist=0.1\", but got {token!r}.")
        try:
            out[key] = ast.literal_eval(value)
        except (ValueError, SyntaxError):
            out[key] = value
    return out


def _action_for(parser, dest: str):
    for action in parser._actions:
        if action.dest == dest and action.option_strings:
            return action
    return None


def _option_for(parser, dest: str) -> str:
    action = _action_for(parser, dest)
    return max(action.option_strings, key=len) if action is not None else dest


OUTPUT_COMPRESSION_HELP = (
    " Writes gzip-compressed output if the path ends in '.gz'; '.bgz' is refused because "
    "HMMER cannot read BGZF. Output to stdout is never compressed."
)

def list_and_file_to_dict_keys(input_list=None, input_file=None, as_set=False):
    """
        args:
            input_list: a list of items to add, in order, as keys to the output dict.

            input_file: The path to a text file, either simple list of strings, a fasta file, a genbank file, or an hmm file. 
                        If a sequence file is passed, then return the sequence ids,
                        If an hmm file is passed then return the hmm profile names
                        If any other kind of file is passed, then just read the lines and pass the lines individually

        output:
            None if input_list and input_file are both None.
            else, a dict where keys are values read from input_list and input_file, and values are None.
    """

    out = None  # dict() we're treating this as an ordered set where values = None
    if input_list is not None and len(input_list) > 0:
        out = dict()
        for c in input_list:
            out[c] = None
    if input_file is not None:
        if out is None:
            out = dict()
        input_file_type = get_file_type(input_file)
        if input_file_type not in {"hmm", "fasta", "genbank"}: #if it's not an hmm file, a fasta file, or a genbank file, we assume it is a text file with one value per line
            with open(input_file, "r") as input_file_handle:
                for line in input_file_handle:
                    line = line.strip()
                    if len(line) != 0:
                        out[line] = None
        elif input_file_type == "hmm":
            # Only the names are wanted, so scan the NAME lines rather than
            # building (and discarding) every profile in the file.
            names = list(iter_hmm_names(input_file))
            if not names and os.path.getsize(input_file) > 0:
                # Not recognizable as hmm text after all; fall back to parsing so
                # that any real format error surfaces from pyhmmer, not silently
                # as an empty selection.
                for hmmdb in read_hmms([input_file]).values():
                    names = list(hmmdb)
            for hmm_name in names:
                out[hmm_name] = None
        elif input_file_type in {"fasta","genbank"}:
            for rec in parse_seqfiles([input_file]):
                out[rec.id] = None

    if as_set and out is not None:
        return set(out.keys())

    return out

def partial_help(what:str) -> str:
    """
        Help text for a --partial argument. what: what is selected, e.g. "Which CDSs/proteins to search".
    """
    return (f"{what}, by fragment status. include: all. exclude: skip fragments. only: only fragments. "
            "A CDS is a fragment if its location has '<' or '>', or it has a /partial qualifier. A protein record is a fragment if it has a UniProt ' (Fragment)' or MGnify 'FL=0' header, "
            "a 'Flags: Fragment' definition or NON_TER features (UniProt text entries), or a CDS/Protein feature with '<' or '>' spanning it.")
