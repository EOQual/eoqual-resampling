"""
Modèles de MTF pour la méthode ``mtf_target``.

Une MTF est une fonction ``mtf(fy, fx)`` de fréquences exprimées en
**cycles par pixel de la grille où elle s'applique** : la grille cible
pour la MTF visée (``target``), la grille source pour la MTF de l'image
d'entrée (``source``). Le Nyquist de chaque grille vaut 0,5.

>>> from eoqual.resampling.mtf import gaussian_mtf, pixel_mtf, product
>>> target = product(gaussian_mtf(0.3), pixel_mtf())   # optique + détecteur
"""

from __future__ import annotations

import math
from typing import Callable

import numpy as np
from numpy.typing import NDArray

__all__ = ["MTF", "gaussian_mtf", "pixel_mtf", "diffraction_mtf", "product"]

#: MTF : ``mtf(fy, fx)`` → tableau de même forme, valeurs dans [0, 1].
MTF = Callable[[NDArray[np.float64], NDArray[np.float64]], NDArray[np.float64]]


def gaussian_mtf(value_at_nyquist: float) -> MTF:
    """
    MTF gaussienne isotrope ``exp(-a f²)`` valant ``value_at_nyquist`` au
    Nyquist (``f = 0,5``) : modèle simple d'une optique et de ses flous.
    """
    if not 0 < value_at_nyquist <= 1:
        raise ValueError("value_at_nyquist doit être dans ]0, 1]")
    a = -math.log(value_at_nyquist) / 0.25

    def mtf(fy, fx):
        return np.exp(-a * (np.asarray(fy) ** 2 + np.asarray(fx) ** 2))

    return mtf


def pixel_mtf(fill_factor: float = 1.0) -> MTF:
    """
    MTF d'intégration sur un pixel carré de côté ``fill_factor`` (en pas
    d'échantillonnage) : ``|sinc(c fy) sinc(c fx)|`` ; 0,64 au Nyquist pour
    un facteur de remplissage de 1.
    """
    if not 0 < fill_factor <= 1:
        raise ValueError("fill_factor doit être dans ]0, 1]")

    def mtf(fy, fx):
        return np.abs(np.sinc(fill_factor * np.asarray(fy)) * np.sinc(fill_factor * np.asarray(fx)))

    return mtf


def diffraction_mtf(cutoff: float) -> MTF:
    """
    MTF de diffraction d'une pupille circulaire non obstruée, de fréquence
    de coupure ``cutoff`` (cycles par pixel, ``D / (λ F)`` multiplié par
    le pas) : ``(2/π)(acos ν - ν √(1 - ν²))``, ``ν = f / cutoff``.
    """
    if cutoff <= 0:
        raise ValueError("cutoff doit être > 0")

    def mtf(fy, fx):
        nu = np.clip(np.hypot(fy, fx) / cutoff, 0.0, 1.0)
        return (2.0 / np.pi) * (np.arccos(nu) - nu * np.sqrt(1.0 - nu * nu))

    return mtf


def product(*mtfs: MTF) -> MTF:
    """Produit de MTF (chaîne optique × détecteur × ...)."""
    if not mtfs:
        raise ValueError("au moins une MTF")

    def mtf(fy, fx):
        out = np.ones(np.broadcast(np.asarray(fy), np.asarray(fx)).shape)
        for m in mtfs:
            out = out * m(fy, fx)
        return out

    return mtf
