"""
Géométrie de la grille cible : facteur, taille de sortie, positions des
pixels cibles exprimées en coordonnées de pixels source.

Coordonnées : le centre du pixel source ``j`` est à la position ``j``
(pixel ``j`` = intervalle ``[j - 1/2, j + 1/2]``). Une réduction d'un
facteur ``k`` (``k >= 1``) produit des pixels cibles de taille ``k``.

Deux conventions de grille (``grid``) :

- ``"center"`` (défaut) : le pixel cible ``i`` couvre exactement le bloc
  source ``[i k - 1/2, (i + 1) k - 1/2]`` ; son centre est en
  ``x_i = (i + 1/2) k - 1/2``. C'est la convention de l'agrégation de
  pixels (binning) : pour ``k`` entier, le pixel cible 0 est la moyenne
  des pixels source ``0 .. k-1``, et l'image réduite couvre la même
  emprise que l'image source.
- ``"origin"`` : le pixel cible ``i`` est centré en ``x_i = i k`` ; le
  premier pixel cible est centré sur le premier pixel source (convention
  d'un sous-échantillonnage « un pixel sur k »). L'emprise est décalée de
  ``(k - 1) / 2`` pixel source par rapport à ``"center"``.
"""

from __future__ import annotations

import math
from typing import Sequence

import numpy as np
from numpy.typing import NDArray

__all__ = ["GRIDS", "resolve_factors", "output_size", "target_positions"]

#: Conventions de grille reconnues.
GRIDS = ("center", "origin")

# Tolérance pour qu'un facteur « presque entier » (ex. 512 / 3 * 3) donne
# la taille attendue malgré l'arrondi flottant.
_EPS = 1e-9


def _check_grid(grid: str) -> None:
    if grid not in GRIDS:
        raise ValueError(f"grid='{grid}' inconnu ; valeurs possibles : {GRIDS}")


def output_size(n_in: int, factor: float) -> int:
    """
    Nombre de pixels cibles entiers contenus dans ``n_in`` pixels source
    pour un facteur ``factor`` : ``floor(n_in / factor)``.

    Les pixels cibles incomplets en fin de ligne sont écartés : chaque
    pixel cible est entièrement couvert par l'image source.
    """
    n = int(math.floor(n_in / factor + _EPS))
    if n < 1:
        raise ValueError(
            f"facteur {factor} trop grand pour une dimension de {n_in} pixels "
            "(aucun pixel cible complet)"
        )
    return n


def resolve_factors(
    in_shape: tuple[int, int],
    factor: float | Sequence[float] | None = None,
    shape: Sequence[int] | None = None,
) -> tuple[tuple[float, float], tuple[int, int]]:
    """
    Facteurs ``(k_lignes, k_colonnes)`` et taille de sortie, à partir de
    ``factor`` **ou** de ``shape``.

    - ``factor`` : scalaire ou couple, ``>= 1`` ; la taille de sortie est
      ``floor(n / k)`` (voir :func:`output_size`).
    - ``shape`` : taille de sortie ``(lignes, colonnes)`` ; le facteur
      vaut alors exactement ``n_in / n_out`` (l'emprise est conservée).

    Raises
    ------
    ValueError
        Si aucun ou les deux paramètres sont donnés, ou si un facteur est
        invalide.
    NotImplementedError
        Si un facteur est inférieur à 1 (agrandissement, prévu mais non
        implémenté dans cette version).
    """
    if (factor is None) == (shape is None):
        raise ValueError("donner soit 'factor', soit 'shape' (un seul des deux)")

    if shape is not None:
        if len(shape) != 2:
            raise ValueError("'shape' doit être un couple (lignes, colonnes)")
        out = (int(shape[0]), int(shape[1]))
        if min(out) < 1:
            raise ValueError("'shape' doit être strictement positif")
        factors = (in_shape[0] / out[0], in_shape[1] / out[1])
    else:
        if np.ndim(factor) == 0:
            factors = (float(factor), float(factor))
        else:
            if len(factor) != 2:
                raise ValueError("'factor' doit être un scalaire ou un couple")
            factors = (float(factor[0]), float(factor[1]))
        if not all(np.isfinite(factors)) or min(factors) <= 0:
            raise ValueError(f"facteur invalide : {factor}")
        out = (output_size(in_shape[0], factors[0]), output_size(in_shape[1], factors[1]))

    if min(factors) < 1.0 - _EPS:
        raise NotImplementedError(
            "agrandissement (facteur < 1) non implémenté dans cette version : "
            "eoqual-resampling ne traite que la réduction de résolution "
            "(voir REFERENCE_TECHNIQUE.md, « Évolutions envisagées »)"
        )
    return factors, out


def target_positions(n_out: int, factor: float, grid: str = "center") -> NDArray[np.float64]:
    """
    Positions des centres des ``n_out`` pixels cibles, en coordonnées de
    pixels source (voir l'en-tête du module pour les conventions).
    """
    _check_grid(grid)
    i = np.arange(n_out, dtype=np.float64)
    if grid == "center":
        return (i + 0.5) * factor - 0.5
    return i * factor
