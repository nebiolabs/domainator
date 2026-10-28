"""Plain-text plots and tables for the text reports."""
from typing import Sequence

import numpy as np


def text_histogram(values: Sequence[float], *, title: str, bins: int = 50, height: int = 20) -> str:
    """A histogram of ``values`` drawn with characters, one column per bin, ``height`` rows tall."""
    counts, edges = np.histogram(np.asarray(values, dtype=float), bins=bins)
    top = int(counts.max()) if len(counts) else 0
    if top == 0:
        return f" {title}\n (no values)\n"
    label_width = len(str(top))
    lines = [f" {title}"]
    for row in range(height, 0, -1):
        threshold = top * (row - 0.5) / height
        label = str(round(top * row / height)) if row in (height, 1) or row % 5 == 0 else ""
        lines.append(f" {label:>{label_width}} | " + "".join("o" if count >= threshold else " " for count in counts))
    lines.append(" " + " " * label_width + " +-" + "-" * len(counts))
    lines.append(" " + " " * label_width + f"   {edges[0]:g} to {edges[-1]:g}, {len(counts)} bins")
    return "\n".join(lines) + "\n"


def text_table(header: Sequence[str], rows: Sequence[Sequence[str]], right_aligned: Sequence[bool]) -> str:
    """Columns of text two spaces apart, each as wide as its widest cell, a header line first.

    The cells are already formatted. A column is right-aligned (numbers) or left-aligned
    (text), its header with it, and every line loses its trailing spaces, so that the layout
    depends only on the cells' lengths in characters."""
    widths = [len(name) for name in header]
    for row in rows:
        widths = [max(width, len(cell)) for width, cell in zip(widths, row)]

    def line(cells):
        return "  ".join(cell.rjust(width) if right else cell.ljust(width)
                         for cell, width, right in zip(cells, widths, right_aligned)).rstrip()

    return "".join(line(cells) + "\n" for cells in [header, *rows])
