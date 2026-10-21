"""Running external programs (diamond, cd-hit, foldseek, ...).

Every tool finds a program with :func:`resolve_binary` and runs it with :func:`run_command`,
so a missing or failing program is always an ``ExternalToolError`` that names the program,
the ``--<program>_path`` option, and the end of the program's error output.
"""
import os
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional, Sequence

from domainator.errors import ExternalToolError


def resolve_binary(name: str, path: Optional[str] = None, *, option: Optional[str] = None, install_hint: str = "") -> str:
    """Return the executable to run for the program ``name``.

    ``path`` is the user's explicit path (the ``--<name>_path`` option, or ``option`` if that's
    spelled differently); without one, ``name`` is looked up on the PATH.
    """
    if path is not None:
        if not (os.path.isfile(path) and os.access(path, os.X_OK)) and shutil.which(path) is None:
            raise ExternalToolError(f"'{path}' is not an executable file.")
        return path
    resolved = shutil.which(name)
    if resolved is None:
        hint = f" {install_hint}" if install_hint else ""
        raise ExternalToolError(
            f"Could not find '{name}' on PATH.{hint} "
            f"Alternatively pass an explicit path with {option or f'--{name}_path'}."
        )
    return resolved


def run_command(argv: Sequence[str], description: str, expect: Optional[str] = None, env=None,
                capture: bool = False):
    """Run an external command, raising ExternalToolError with the end of its stderr on failure.

    A nonzero exit code alone is not treated as failure when ``expect`` names an artifact
    that was produced: foldseek's multi-step workflows exit nonzero during their cleanup
    phase even when the run succeeded.

    subprocess.run rather than Popen followed by stderr.read(), which deadlocks once the
    child writes more than the pipe buffer.

    Returns the CompletedProcess when capture is set, for callers that need to inspect the
    command's own output (reseek reports how many chains it skipped, for example).
    """
    argv = [str(a) for a in argv]
    try:
        completed = subprocess.run(argv, capture_output=True, env=env)
    except FileNotFoundError:
        raise ExternalToolError(f"{description}: could not run '{argv[0]}'; it is not an executable file.") from None
    if completed.returncode == 0:
        return completed if capture else None
    if expect is not None and Path(expect).exists():
        return completed if capture else None
    output = completed.stderr.decode("utf-8", errors="replace").strip() or completed.stdout.decode("utf-8", errors="replace").strip()
    tail = "\n".join(output.split("\n")[-25:])
    raise ExternalToolError(
        f"{description} failed with exit code {completed.returncode}.\n"
        f"Command: {' '.join(argv)}\n{tail}"
    )


def append_extra_args(argv: List[str], extra_args: Optional[Sequence[str]], option: str) -> List[str]:
    """Return ``argv`` followed by a user's extra arguments for an external program.

    Raises ValueError if the extra arguments set an option that ``argv`` already sets, naming
    ``option`` (the Domainator option they came from) in the message.
    """
    extra_args = list(extra_args or ())
    conflicts = sorted({token for token in argv if isinstance(token, str) and token.startswith("-")}
                       .intersection(token for token in extra_args if token.startswith("-")))
    if conflicts:
        raise ValueError(f"{option} cannot set {', '.join(conflicts)}: Domainator sets "
                         f"{'it' if len(conflicts) == 1 else 'them'} itself.")
    return list(argv) + extra_args
