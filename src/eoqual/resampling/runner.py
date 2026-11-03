"""
Point d'entrée unique : réduire une image par une méthode choisie par son
nom, et lister les méthodes disponibles.

>>> from eoqual.resampling import resample, list_methods
>>> small = resample(image, 3)                        # sinc fenêtré par Kaiser, facteur 3
>>> small = resample(image, shape=(400, 600), method="sinc_kaiser")
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .config import DEFAULT_METHOD, METHODS_CONFIGS
from .geometry import resolve_factors
from .methods import check_method, resample_2d

__all__ = ["resample", "list_methods"]


def list_methods() -> None:
    """Affiche les méthodes, leur famille, leurs paramètres par défaut et leur référence."""
    from rich.console import Console
    from rich.table import Table

    table = Table(title="Méthodes de eoqual.resampling", header_style="bold magenta")
    table.add_column("Méthode", style="bold cyan", no_wrap=True)
    table.add_column("Famille")
    table.add_column("Paramètres (défaut)")
    table.add_column("Référence")
    for name, cfg in METHODS_CONFIGS.items():
        params = ", ".join(f"{k}={v}" for k, v in cfg["params"].items()) or "—"
        label = name + (" (défaut)" if name == DEFAULT_METHOD else "")
        table.add_row(label, cfg["family"], params, cfg["reference"])
    Console().print(table)


def resample(
    image: ArrayLike,
    factor: float | Sequence[float] | None = None,
    *,
    shape: Sequence[int] | None = None,
    method: str = DEFAULT_METHOD,
    grid: str = "center",
    boundary: str = "symmetric",
    channel_axis: int | None = None,
    **params,
) -> NDArray[np.float64]:
    """
    Réduit la résolution d'une image.

    Parameters
    ----------
    image : array_like
        Image 2D, ou 3D avec ``channel_axis`` (chaque bande est réduite de
        la même façon). Valeurs finies (pas de NaN).
    factor : float ou (float, float), optional
        Facteur de réduction ``>= 1`` (taille d'un pixel cible en pixels
        source), commun ou par axe ``(lignes, colonnes)``, entier ou non.
        La taille de sortie est ``floor(n / factor)``.
    shape : (int, int), optional
        Taille de sortie, à la place de ``factor`` : le facteur vaut alors
        ``n / shape`` sur chaque axe (emprise conservée).
    method : str
        Une des méthodes de :func:`list_methods` ; ``"sinc_kaiser"`` par
        défaut.
    grid : {"center", "origin"}
        Convention de la grille cible (voir
        :mod:`eoqual.resampling.geometry`) : ``"center"`` (défaut) aligne
        chaque pixel cible sur le bloc de pixels source qu'il recouvre.
    boundary : {"symmetric", "reflect", "edge", "wrap"}
        Prolongement de l'image au-delà de ses bords (méthodes
        séparables). Les méthodes spectrales traitent les bords par la
        décomposition périodique + lisse, sauf ``"wrap"`` (image
        périodique).
    channel_axis : int, optional
        Axe des bandes d'une image multi-bandes.
    **params
        Paramètres de la méthode (voir ``METHODS_CONFIGS[method]["params"]``
        et ``REFERENCE_TECHNIQUE.md``).

    Returns
    -------
    np.ndarray
        Image réduite, ``float64``. Une image uniforme est conservée
        exactement ; la moyenne d'une image quelconque l'est au repliement
        près.

    Raises
    ------
    NotImplementedError
        Si un facteur est inférieur à 1 (agrandissement non implémenté).
    """
    check_method(method)
    arr = np.asarray(image, dtype=np.float64)
    if channel_axis is None:
        if arr.ndim != 2:
            raise ValueError("image 2D attendue (ou préciser channel_axis pour une image multi-bandes)")
        bands = None
    else:
        if arr.ndim != 3:
            raise ValueError("channel_axis suppose une image 3D")
        bands = np.moveaxis(arr, channel_axis, 0)
    if not np.all(np.isfinite(arr)):
        raise ValueError("l'image contient des valeurs non finies (NaN, inf)")

    plane_shape = arr.shape if bands is None else bands.shape[1:]
    factors, out_shape = resolve_factors(plane_shape, factor, shape)

    if bands is None:
        return resample_2d(arr, method, factors, out_shape, grid, boundary, **params)
    out = np.stack([resample_2d(b, method, factors, out_shape, grid, boundary, **params) for b in bands])
    return np.moveaxis(out, 0, channel_axis)
