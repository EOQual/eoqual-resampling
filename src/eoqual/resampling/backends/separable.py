"""
Rééchantillonnage séparable par noyau évalué **directement aux positions
cibles**.

Pour chaque axe, une matrice creuse ``W`` (``n_out × n_in``) contient,
ligne ``i``, les poids ``h(x_i - j)`` des pixels source ``j`` voisins de
la position cible ``x_i`` (voir :mod:`eoqual.resampling.geometry`).
L'image réduite vaut ``W_lignes · image · W_colonnesᵀ``.

Cette évaluation directe est exacte pour un facteur quelconque (non
entier, différent par axe) : il n'y a ni filtrage sur la grille source
suivi d'une interpolation, ni approximation polyphasée. Chaque ligne de
``W`` est normalisée à une somme de 1 : une image uniforme est conservée
exactement, quel que soit le décalage sous-pixel de la position cible.
"""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import NDArray
from scipy import sparse

from ..geometry import target_positions
from ..kernels import Kernel

__all__ = ["BOUNDARIES", "fold_indices", "resampling_matrix", "box_matrix", "apply_separable"]

#: Traitements des bords reconnus.
BOUNDARIES = ("symmetric", "reflect", "edge", "wrap")


def fold_indices(j: NDArray[np.int64], n: int, boundary: str) -> NDArray[np.int64]:
    """
    Ramène des indices quelconques dans ``[0, n[`` selon ``boundary`` :

    - ``"symmetric"`` : miroir demi-pixel (``-1 → 0``, ``n → n-1``),
      cohérent avec la convention « centre de pixel » ; défaut ;
    - ``"reflect"`` : miroir sur le pixel de bord (``-1 → 1``) ;
    - ``"edge"`` : répétition du pixel de bord ;
    - ``"wrap"`` : image périodique.
    """
    if boundary == "symmetric":
        period = 2 * n
        m = np.mod(j, period)
        return np.where(m < n, m, period - 1 - m)
    if boundary == "reflect":
        if n == 1:
            return np.zeros_like(j)
        period = 2 * n - 2
        m = np.mod(j, period)
        return np.where(m < n, m, period - m)
    if boundary == "edge":
        return np.clip(j, 0, n - 1)
    if boundary == "wrap":
        return np.mod(j, n)
    raise ValueError(f"boundary='{boundary}' inconnu ; valeurs possibles : {BOUNDARIES}")


def _assemble(weights, cols, n_out, n_in, boundary) -> sparse.csr_matrix:
    """Normalise les lignes, replie les indices et construit la matrice creuse."""
    sums = weights.sum(axis=1, keepdims=True)
    if np.any(np.abs(sums) < 1e-12):
        raise ValueError("noyau de somme nulle à une position cible (support trop court ?)")
    weights = weights / sums
    cols = fold_indices(cols, n_in, boundary)
    rows = np.repeat(np.arange(n_out), weights.shape[1])
    # Les doublons (indices repliés sur le même pixel) sont additionnés.
    return sparse.csr_matrix((weights.ravel(), (rows, cols.ravel())), shape=(n_out, n_in))


def resampling_matrix(
    n_in: int,
    n_out: int,
    factor: float,
    kernel: Kernel,
    grid: str = "center",
    boundary: str = "symmetric",
) -> sparse.csr_matrix:
    """
    Matrice ``W`` (``n_out × n_in``) d'un axe pour le noyau ``kernel``.

    Le noyau est dilaté du facteur si ``kernel.antialias`` (filtre
    anti-repliement), sinon évalué tel quel (interpolation sans
    préfiltrage).
    """
    x = target_positions(n_out, factor, grid)
    scale = max(factor, 1.0) if kernel.antialias else 1.0
    r = kernel.radius * scale
    span = int(math.ceil(2 * r)) + 2
    start = np.floor(x - r).astype(np.int64)
    cols = start[:, None] + np.arange(span)[None, :]
    u = (x[:, None] - cols) / scale
    # Troncature stricte au support déclaré : sinon un poids résiduel d'un
    # seul côté (noyau non nul au-delà de son rayon) décalerait l'image.
    weights = np.where(np.abs(u) < kernel.radius, kernel.phi(u), 0.0)
    return _assemble(weights, cols, n_out, n_in, boundary)


def box_matrix(
    n_in: int,
    n_out: int,
    factor: float,
    grid: str = "center",
    boundary: str = "symmetric",
) -> sparse.csr_matrix:
    """
    Agrégation exacte (« boîte », moyenne par surface) : le poids du pixel
    source ``j`` est la longueur de son recouvrement avec l'empreinte
    ``[x_i - k/2, x_i + k/2]`` du pixel cible. Pour ``k`` entier et la
    grille ``"center"``, c'est exactement la moyenne par blocs (binning).
    """
    x = target_positions(n_out, factor, grid)
    half = 0.5 * factor
    span = int(math.ceil(factor)) + 3
    start = np.floor(x - half - 0.5).astype(np.int64)
    cols = start[:, None] + np.arange(span)[None, :]
    lo = np.maximum(cols - 0.5, (x - half)[:, None])
    hi = np.minimum(cols + 0.5, (x + half)[:, None])
    weights = np.clip(hi - lo, 0.0, None)
    return _assemble(weights, cols, n_out, n_in, boundary)


def apply_separable(
    image: NDArray[np.float64],
    w_rows: sparse.csr_matrix,
    w_cols: sparse.csr_matrix,
) -> NDArray[np.float64]:
    """``W_lignes · image · W_colonnesᵀ`` (image 2D float64)."""
    tmp = w_rows @ image            # (n_out_lignes, n_in_colonnes)
    return np.asarray((w_cols @ tmp.T).T)
