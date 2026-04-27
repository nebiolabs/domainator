"""Shared helpers for output-size guardrails."""

from __future__ import annotations

import os
from pathlib import Path
import tempfile


DEFAULT_MAX_OUTPUT_GB = 25.0
BYTES_PER_GB = 1_000_000_000
DEFAULT_MAX_OUTPUT_BYTES = int(DEFAULT_MAX_OUTPUT_GB * BYTES_PER_GB)


class OutputSizeLimitExceeded(RuntimeError):
    """Raised when a projected output exceeds the configured size limit."""


def max_output_gb_to_bytes(max_output_gb: float) -> int | None:
    """Convert a CLI size limit in GB to bytes, where 0 disables the guardrail."""
    if max_output_gb < 0:
        raise ValueError("--max_output_gb must be >= 0.")
    if max_output_gb == 0:
        return None
    return int(max_output_gb * BYTES_PER_GB)


def add_max_output_gb_argument(parser) -> None:
    """Add the shared output-size limit argument to a CLI parser."""
    parser.add_argument(
        "--max_output_gb",
        type=float,
        default=DEFAULT_MAX_OUTPUT_GB,
        help=(
            "Fail if a projected output file would exceed this many GB. "
            "Set to 0 to disable the guardrail."
        ),
    )


def format_bytes(num_bytes: int) -> str:
    """Return a human-readable decimal byte count."""
    if num_bytes < 1000:
        return f"{num_bytes} B"

    value = float(num_bytes)
    for unit in ("KB", "MB", "GB", "TB", "PB"):
        value /= 1000.0
        if value < 1000 or unit == "PB":
            return f"{value:.2f} {unit}"

    return f"{num_bytes} B"


def build_output_limit_message(
    *,
    output_description: str,
    projected_bytes: int,
    max_output_bytes: int,
    mitigation_options: list[str] | tuple[str, ...] | None = None,
    extra_guidance: str | None = None,
) -> str:
    """Build a standardized size-limit failure message."""
    message = (
        f"Projected {output_description} size is {format_bytes(projected_bytes)}, "
        f"which exceeds the configured limit of {format_bytes(max_output_bytes)}."
    )
    message += " Reduce the projected output size or raise the limit with --max_output_gb."

    if mitigation_options:
        message += " Relevant options: " + ", ".join(mitigation_options) + "."

    if extra_guidance:
        message += f" {extra_guidance}"

    return message


def enforce_output_limit(
    *,
    projected_bytes: int,
    max_output_bytes: int | None,
    output_description: str,
    mitigation_options: list[str] | tuple[str, ...] | None = None,
    extra_guidance: str | None = None,
) -> None:
    """Raise when a projected output exceeds the configured size limit."""
    if max_output_bytes is None:
        return
    if projected_bytes <= max_output_bytes:
        return

    raise OutputSizeLimitExceeded(
        build_output_limit_message(
            output_description=output_description,
            projected_bytes=projected_bytes,
            max_output_bytes=max_output_bytes,
            mitigation_options=mitigation_options,
            extra_guidance=extra_guidance,
        )
    )


def enforce_matrix_output_limit(
    *,
    output_type: str,
    matrix,
    row_names: list[str],
    col_names: list[str],
    max_output_bytes: int | None,
    output_path: str,
    row_lengths=None,
    col_lengths=None,
    data_type: str = "",
    mitigation_options: list[str] | tuple[str, ...] | None = None,
    extra_guidance: str | None = None,
) -> int:
    """Estimate a matrix output and enforce the configured size limit."""
    from domainator.data_matrix import DataMatrix

    projected_bytes = DataMatrix.estimate_write_size(
        output_type,
        matrix,
        row_names,
        col_names,
        row_lengths,
        col_lengths,
        data_type,
    )
    enforce_output_limit(
        projected_bytes=projected_bytes,
        max_output_bytes=max_output_bytes,
        output_description=f"{output_type} matrix output '{output_path}'",
        mitigation_options=mitigation_options,
        extra_guidance=extra_guidance,
    )
    return projected_bytes


def make_temporary_output_path(target_path: str) -> str:
    """Create a temporary output path in the target directory."""
    target = Path(target_path)
    file_descriptor, temp_path = tempfile.mkstemp(
        prefix=f".{target.name}.",
        suffix=".tmp",
        dir=str(target.parent),
    )
    os.close(file_descriptor)
    return temp_path