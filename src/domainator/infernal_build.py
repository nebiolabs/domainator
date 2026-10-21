"""Build a covariance model (CM) from a multiple sequence alignment (MSA).

NOTE: EXPERIMENTAL!

Allows the user to specify the ACC, NAME, and DESC fields of the CM profile.
Some of these options are not available in the standard cmbuild tool from Infernal.
The input alignment is passed directly to cmbuild, so any aligned MSA format it
accepts can be used here. Formats without structure annotation, such as aligned
FASTA, generally require --noss.
"""
from jsonargparse import ActionConfigFile
from domainator.cli.common import add_input, add_output, check_outputs, declare_outputs, STDIN_PATHS, tool_parser
from domainator.seqio import open_writable_seqfile
from domainator.external import resolve_binary, run_command
import os
import re
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


def infernal_build(
    file: Union[str, os.PathLike, BinaryIO],
    *,
    name: str,
    acc: Optional[str] = None,
    desc: Optional[str] = None,
    informat: Optional[str] = None,
    cmbuild_path: Optional[str] = None,
    noss: bool = False,
    hand: bool = False,
    symfrac: Optional[float] = None,
) -> str:
    cmbuild = resolve_binary("cmbuild", cmbuild_path, install_hint="Install Infernal.")

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

        run_command(args, f"cmbuild on {input_path}")

        with open(summary_path, "r", encoding="utf-8") as summary_handle:
            summary_text = summary_handle.read()
        if _count_summary_models(summary_text) != 1:
            raise ValueError("domainator infernal_build currently supports exactly one alignment per invocation")

        with open(raw_output_path, "r", encoding="utf-8") as cm_handle:
            cm_text = cm_handle.read()
        return _patch_cm_header(cm_text, name=name, acc=acc, desc=desc)


def build_parser():
    parser = tool_parser(__doc__)

    add_input(parser, "Path of input aligned MSA. Any format accepted by cmbuild can be used. Formats without structure annotation, such as aligned FASTA, generally require --noss.", multiple=False)
    add_output(parser, "cm", "Write the CM to this file.", primary=True)
    parser.add_argument("--name", default=None, required=True, type=str, help="Name of the CM profile.")
    parser.add_argument("--acc", default=None, required=False, type=str, help="Accession of the CM profile.")
    parser.add_argument("--desc", default=None, required=False, type=str, help="Description of the CM profile.")
    parser.add_argument(
        "--input_format",
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
    declare_outputs(parser, ["cm_out"])
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)
    check_outputs(parser, params)
    input_file = sys.stdin.buffer if params.input in STDIN_PATHS else params.input

    cm_text = infernal_build(
        file=input_file,
        name=params.name,
        acc=params.acc,
        desc=params.desc,
        informat=params.input_format,
        cmbuild_path=params.cmbuild_path,
        noss=params.noss,
        hand=params.hand,
        symfrac=params.symfrac,
    )

    with open_writable_seqfile(params.cm_out) as output_handle:
        output_handle.write(cm_text)
