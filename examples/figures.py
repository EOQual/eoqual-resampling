#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figures de la documentation (PNG dans ``resources/``).

    python examples/figures.py

Nécessite matplotlib (non requis par le paquet lui-même).
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

from eoqual.resampling.geometry import target_positions  # noqa: E402

RESOURCES = Path(__file__).resolve().parents[1] / "resources"


def figure_grille(k=3, n_source=9, chemin=RESOURCES / "grille_cible.png"):
    """Conventions de grille "center" et "origin" pour un facteur k."""
    fig, ax = plt.subplots(figsize=(9, 3.6), layout="constrained")
    lignes = [("pixels source", 2.0), ('grid="center"', 0.8), ('grid="origin"', -0.6)]
    n_out = int(n_source // k)          # pixels cibles produits : floor(n / k)

    # Pixels source : intervalles [j - 1/2, j + 1/2].
    for j in range(n_source):
        ax.add_patch(Rectangle((j - 0.5, 1.75), 1, 0.5, fc="0.93", ec="0.4"))
        ax.text(j, 2.0, str(j), ha="center", va="center", fontsize=10)

    # Pixels cibles : centres x_i, empreinte [x_i - k/2, x_i + k/2].
    couleurs = {"center": "#1f77b4", "origin": "#d62728"}
    for grid, y in (("center", 0.8), ("origin", -0.6)):
        x = target_positions(n_out, k, grid)
        for i, xi in enumerate(x):
            lo, hi = xi - k / 2, xi + k / 2
            lo_vis, hi_vis = max(lo, -1.0), min(hi, n_source - 0.5)
            ax.add_patch(Rectangle((lo_vis, y - 0.25), hi_vis - lo_vis, 0.5,
                                   fc=couleurs[grid], alpha=0.18, ec=couleurs[grid], lw=1.5))
            ax.plot(xi, y, "o", color=couleurs[grid], ms=6)
            ax.text(xi, y - 0.42, f"cible {i}\n$x_{i}$ = {xi:g}", ha="center", va="top", fontsize=8.5)
            ax.plot([xi, xi], [y + 0.25, 1.75], ":", color=couleurs[grid], lw=1)

    for nom, y in lignes:
        ax.text(-1.3, y, nom, ha="right", va="center", fontsize=10,
                family="monospace" if "grid" in nom else None)
    ax.set_xlim(-3.6, n_source - 0.3)
    ax.set_ylim(-1.5, 2.4)
    ax.axis("off")
    ax.set_title(f"Position des pixels cibles, facteur k = {k} (coordonnées en pixels source)", fontsize=11)
    fig.savefig(chemin, dpi=150)
    plt.close(fig)
    return chemin


if __name__ == "__main__":
    print(figure_grille())
