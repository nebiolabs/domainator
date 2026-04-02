"""Build a covariance model (CM) from a multiple sequence alignment (MSA).

NOTE: EXPERIMENTAL!

Allows the user to specify the ACC, NAME, and DESC fields of the CM profile.
Some of these options are not available in the standard cmbuild tool from Infernal.
The input alignment is passed directly to cmbuild, so any aligned MSA format it
accepts can be used here. Formats without structure annotation, such as aligned
FASTA, generally require --noss.
"""
from jsonargparse import ArgumentParser, ActionConfigFile
from domainator import __version__, RawAndDefaultsFormatter
import os
import re
import shutil
import subprocess
import sys
import tempfile
from typing import BinaryIO, Optional, Union


def sanitize_string(s: str) -> str:
    return re.sub(r"[^ \w\d_\-\.;:]", "_", s)


def _count_summary_models(summary_text: str) -> int:
    count = 0
    for line in summary_text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        count += 1
    return count


def _set_header_field(header_lines, key: str, value: Optional[str]):
    if value is None:
        return

    formatted = f"{key:<9}{sanitize_string(value)}\n"
    for index, line in enumerate(header_lines):
        if line.startswith(f"{key:<9}"):
            header_lines[index] = formatted
            return

    insert_at = len(header_lines)
    preferred_order = ["NAME", "ACC", "DESC"]
    key_position = preferred_order.index(key)
    for previous_key in reversed(preferred_order[:key_position]):
        for index, line in enumerate(header_lines):
            if line.startswith(f"{previous_key:<9}"):
                insert_at = index + 1
                header_lines.insert(insert_at, formatted)
                return

    for index, line in enumerate(header_lines):
        if line.startswith("STATES"):
            insert_at = index
            break
    header_lines.insert(insert_at, formatted)


def _patch_cm_header(cm_text: str, name: Optional[str], acc: Optional[str], desc: Optional[str]) -> str:
    header_lines = []
    body_lines = []
    in_body = False

    for line in cm_text.splitlines(keepends=True):
        if not in_body and line.strip() == "CM":
            in_body = True
        if in_body:
            body_lines.append(line)
        else:
            header_lines.append(line)

    if not body_lines:
        raise ValueError("Invalid CM output: could not find CM body")

    _set_header_field(header_lines, "NAME", name)
    _set_header_field(header_lines, "ACC", acc)
    _set_header_field(header_lines, "DESC", desc)
    return "".join(header_lines + body_lines)


def _resolve_cmbuild_path(cmbuild_path: Optional[str]) -> str:
    if cmbuild_path is not None:
        return cmbuild_path

    resolved = shutil.which("cmbuild")
    if resolved is None:
        raise RuntimeError("Could not find 'cmbuild' on PATH. Install Infernal or pass --cmbuild_path.")
    return resolved


def infernal_build(
    file: Union[str, os.PathLike, BinaryIO],
    name: str,
    acc: Optional[str] = None,
    desc: Optional[str] = None,
    informat: Optional[str] = None,
    cmbuild_path: Optional[str] = None,
    noss: bool = False,
    hand: bool = False,
    symfrac: Optional[float] = None,
) -> str:
    cmbuild = _resolve_cmbuild_path(cmbuild_path)

    with tempfile.TemporaryDirectory() as tmpdir:
        raw_output_path = os.path.join(tmpdir, "output.cm")
        summary_path = os.path.join(tmpdir, "summary.txt")

        if isinstance(file, (str, os.PathLike)):
            input_path = os.fspath(file)
        else:
            input_path = os.path.join(tmpdir, "input.msa")
            input_content = file.read()
            if isinstance(input_content, str):
                input_content = input_content.encode()
            with open(input_path, "wb") as input_handle:
                input_handle.write(input_content)

        args = [cmbuild, "-F", "-o", summary_path, "-n", sanitize_string(name)]
        if informat is not None:
            args.extend(["--informat", informat])
        if noss:
            args.append("--noss")
        if hand:
            args.append("--hand")
        if symfrac is not None:
            args.extend(["--symfrac", str(symfrac)])
        args.extend([raw_output_path, input_path])

        result = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8")
        if result.returncode != 0:
            raise RuntimeError(
                f"cmbuild failed for input {input_path}:\n{result.stdout}\n{result.stderr}".strip()
            )

        with open(summary_path, "r", encoding="utf-8") as summary_handle:
            summary_text = summary_handle.read()
        if _count_summary_models(summary_text) != 1:
            raise ValueError("infernal_build.py currently supports exactly one alignment per invocation")

        with open(raw_output_path, "r", encoding="utf-8") as cm_handle:
            cm_text = cm_handle.read()
        return _patch_cm_header(cm_text, name=name, acc=acc, desc=desc)


def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument(
        "-i",
        "--input",
        default=None,
        required=False,
        type=str,
        help="Path of input aligned MSA. Any format accepted by cmbuild can be used. Formats without structure annotation, such as aligned FASTA, generally require --noss. If not supplied, reads from stdin.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        required=False,
        type=str,
        help="CM output file path. If not supplied writes to stdout.",
    )
    parser.add_argument("--name", default=None, required=True, type=str, help="Name of the CM profile.")
    parser.add_argument("--acc", default=None, required=False, type=str, help="Accession of the CM profile.")
    parser.add_argument("--desc", default=None, required=False, type=str, help="Description of the CM profile.")
    parser.add_argument(
        "--informat",
        default=None,
        required=False,
        type=str,
        help="Optional alignment format to pass through to cmbuild, such as stockholm, selex, or afa.",
    )
    parser.add_argument("--cmbuild_path", default=None, required=False, type=str, help="Path to the cmbuild executable.")
    parser.add_argument("--noss", action="store_true", help="Ignore secondary-structure annotation and build a zero-basepair model.")
    parser.add_argument("--hand", action="store_true", help="Use Stockholm RF annotation to define consensus columns.")
    parser.add_argument("--symfrac", default=None, required=False, type=float, help="Consensus-column threshold to pass through to cmbuild.")
    parser.add_argument("--config", action=ActionConfigFile)

    params = parser.parse_args(argv)

    if params.input is None:
        input_file = sys.stdin.buffer
    else:
        input_file = params.input

    cm_text = infernal_build(
        file=input_file,
        name=params.name,
        acc=params.acc,
        desc=params.desc,
        informat=params.informat,
        cmbuild_path=params.cmbuild_path,
        noss=params.noss,
        hand=params.hand,
        symfrac=params.symfrac,
    )

    if params.output is None:
        sys.stdout.write(cm_text)
    else:
        with open(params.output, "w", encoding="utf-8") as output_handle:
            output_handle.write(cm_text)


def _entrypoint():
    main(sys.argv[1:])


if __name__ == "__main__":
    main(sys.argv[1:])