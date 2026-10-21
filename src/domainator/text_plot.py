"""Plain-text plots for the text reports."""
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
