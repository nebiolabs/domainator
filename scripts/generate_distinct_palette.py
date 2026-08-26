"""Regenerate the DISTINCT_COLORS palette used by domainator.utils.get_palette.

Greedy maximin ("glasbey") selection in CIELAB: repeatedly pick the candidate color
whose nearest already-chosen color is as far away as possible. Seeded with
matplotlib's tab10 so the first ten entries stay familiar, and with white/black so
nothing lands too close to a plot background.

Candidates are restricted to 25 < L* < 85: darker colors are hard to tell apart and
lighter ones wash out on a white background.

Usage:
    python scripts/generate_distinct_palette.py [N]

Prints the hex list ready to paste into src/domainator/utils.py, plus the minimum
pairwise deltaE, which is the number that matters: as long as it stays comfortably
above ~20, every pair of colors in the palette is distinguishable, so consecutively
numbered groups are guaranteed to look different.
"""

import itertools
import sys

import matplotlib
import numpy as np

# Candidate grid resolution per channel, and the L* window colors must fall in.
GRID_STEPS = 21
MIN_LIGHTNESS = 25
MAX_LIGHTNESS = 85


def srgb_to_lab(rgb):
    """rgb: (..., 3) array of floats in [0, 1]. Returns (..., 3) CIELAB (D65)."""
    rgb = np.asarray(rgb, dtype=float)
    linear = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    to_xyz = np.array(
        [
            [0.4124564, 0.3575761, 0.1804375],
            [0.2126729, 0.7151522, 0.0721750],
            [0.0193339, 0.1191920, 0.9503041],
        ]
    )
    xyz = (linear @ to_xyz.T) / np.array([0.95047, 1.0, 1.08883])
    delta = 6 / 29
    f = np.where(xyz > delta ** 3, np.cbrt(xyz), xyz / (3 * delta * delta) + 4 / 29)
    return np.stack(
        [116 * f[..., 1] - 16, 500 * (f[..., 0] - f[..., 1]), 200 * (f[..., 1] - f[..., 2])],
        axis=-1,
    )


def generate(n_colors):
    grid = np.linspace(0, 1, GRID_STEPS)
    candidates = np.array(list(itertools.product(grid, grid, grid)))
    candidate_lab = srgb_to_lab(candidates)

    in_range = (candidate_lab[..., 0] > MIN_LIGHTNESS) & (candidate_lab[..., 0] < MAX_LIGHTNESS)
    candidates, candidate_lab = candidates[in_range], candidate_lab[in_range]

    anchors = [matplotlib.colors.to_rgb(c) for c in matplotlib.colormaps["tab10"].colors]
    palette = list(anchors)
    palette_lab = list(srgb_to_lab(np.array(anchors)))

    # Distance from each candidate to the nearest color we must stay away from.
    backgrounds = srgb_to_lab(np.array([[1.0, 1.0, 1.0], [0.0, 0.0, 0.0]]))
    nearest = np.linalg.norm(candidate_lab[:, None, :] - backgrounds[None, :, :], axis=-1).min(1)
    for lab in palette_lab:
        nearest = np.minimum(nearest, np.linalg.norm(candidate_lab - lab, axis=-1))

    while len(palette) < n_colors:
        best = int(np.argmax(nearest))
        palette.append(candidates[best])
        palette_lab.append(candidate_lab[best])
        nearest = np.minimum(nearest, np.linalg.norm(candidate_lab - candidate_lab[best], axis=-1))

    return palette[:n_colors], np.array(palette_lab[:n_colors])


def main(argv):
    n_colors = int(argv[0]) if argv else 64
    palette, palette_lab = generate(n_colors)

    distances = np.linalg.norm(palette_lab[:, None, :] - palette_lab[None, :, :], axis=-1)
    np.fill_diagonal(distances, np.inf)
    print(f"# {n_colors} colors, minimum pairwise deltaE = {distances.min():.1f}")

    hexes = ["#%02X%02X%02X" % tuple(int(round(c * 255)) for c in rgb) for rgb in palette]
    print("DISTINCT_COLORS = (")
    for i in range(0, n_colors, 6):
        print("    " + " ".join('"%s",' % h for h in hexes[i : i + 6]))
    print(")")


if __name__ == "__main__":
    main(sys.argv[1:])
