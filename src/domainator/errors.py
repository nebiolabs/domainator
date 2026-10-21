"""Exceptions Domainator raises for problems a user can fix.

The ``domainator`` command prints a ``DomainatorError`` as a one-line message, without a
traceback (``-v`` shows the traceback). Each subclass also derives from the built-in exception
it refines, so ``except ValueError`` and ``except RuntimeError`` keep catching them.
"""


class DomainatorError(Exception):
    """Base class for errors caused by inputs, options, or the environment, not by a bug."""


class InputError(DomainatorError, ValueError):
    """An input file or option value can't be used."""


class ExternalToolError(DomainatorError, RuntimeError):
    """An external program (diamond, cd-hit, foldseek, ...) is missing or failed."""


class OutputSizeLimitExceeded(DomainatorError, RuntimeError):
    """An output would exceed the --max_output_gb guardrail."""


class FormatVersionError(DomainatorError, ValueError):
    """A file was written in a format version this Domainator can't read."""
