"""
Figures Plotly (via eoqual-viz) : fonctions de transfert des méthodes.
"""

from __future__ import annotations

from typing import Sequence

from eoqual.viz import plot_curves

from .analysis import transfer_function

__all__ = ["plot_transfer_functions"]


def plot_transfer_functions(
    methods: Sequence[str],
    factor: float,
    params: dict[str, dict] | None = None,
    title: str = "",
):
    """
    Superpose le gain (valeur absolue) des méthodes ``methods`` pour un
    facteur ``factor``, avec un repère au Nyquist cible. ``params`` :
    paramètres par méthode, ``{nom: {param: valeur}}``.
    """
    params = params or {}
    curves = []
    for m in methods:
        f, g = transfer_function(m, factor, **params.get(m, {}))
        curves.append((m, f, abs(g)))
    return plot_curves(
        curves,
        vlines=[(0.5, "Nyquist cible", "black")],
        xlabel="fréquence (cycles par pixel cible)",
        ylabel="|gain|",
        title=title or f"Fonctions de transfert, facteur {factor:g}",
    )
