"""
Noyaux continus 1D des méthodes séparables.

Un noyau est décrit par :

- ``phi(u)`` : fonction vectorisée de la variable sans dimension ``u`` ;
- ``radius`` : demi-support (``phi(u) = 0`` pour ``|u| >= radius``) ;
- ``antialias`` : si vrai, le noyau est **dilaté du facteur de réduction**
  ``k`` : appliqué en pixels source, il vaut ``h(x) = phi(x / k) / k``
  (``u`` est alors en pixels **cibles**, et la fréquence de coupure suit
  le nouveau Nyquist). Sinon, ``h(x) = phi(x)`` : simple interpolation
  aux positions cibles, **sans préfiltrage** (méthodes de référence qui
  replient le spectre).

Les noyaux ne sont pas normalisés ici : chaque ligne de la matrice de
rééchantillonnage est ramenée à une somme de 1 (voir
:mod:`eoqual.resampling.backends.separable`), ce qui conserve exactement
les niveaux d'une image uniforme quel que soit le décalage sous-pixel.

Références : Keys (1981) pour la cubique, Mitchell & Netravali (1988),
Duchon (1979) pour Lanczos, Kaiser (1974) pour la fenêtre de Kaiser —
voir ``BIBLIOGRAPHY.md``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable

import numpy as np
from numpy.typing import NDArray
from scipy.special import i0

__all__ = [
    "Kernel",
    "nearest_kernel",
    "linear_kernel",
    "keys_kernel",
    "mitchell_kernel",
    "gaussian_kernel",
    "lanczos_kernel",
    "sinc_gaussian_kernel",
    "kaiser_beta",
    "kaiser_radius",
    "sinc_kaiser_kernel",
]

FWHM_TO_SIGMA = 1.0 / (2.0 * math.sqrt(2.0 * math.log(2.0)))


@dataclass(frozen=True)
class Kernel:
    """Noyau continu 1D (voir l'en-tête du module)."""

    phi: Callable[[NDArray[np.float64]], NDArray[np.float64]]
    radius: float
    antialias: bool


# --- Méthodes de référence, sans préfiltrage --------------------------------

def nearest_kernel() -> Kernel:
    """
    Plus proche voisin : 1 sur ``[-1/2, 1/2[``. Avec la grille
    ``"center"`` et un facteur pair, le centre du pixel cible tombe entre
    deux pixels source : le second est retenu (décalage d'un demi-pixel
    source, inhérent à la méthode).
    """
    return Kernel(lambda u: ((u >= -0.5) & (u < 0.5)).astype(np.float64), 0.5 + 1e-12, False)


def linear_kernel(antialias: bool = False) -> Kernel:
    """Triangle (interpolation bilinéaire en séparable)."""
    return Kernel(lambda u: np.clip(1.0 - np.abs(u), 0.0, None), 1.0, antialias)


def _cubic(u: NDArray[np.float64], b: float, c: float) -> NDArray[np.float64]:
    """Famille de cubiques par morceaux de Mitchell & Netravali (1988)."""
    a = np.abs(u)
    a2, a3 = a * a, a * a * a
    inner = ((12 - 9 * b - 6 * c) * a3 + (-18 + 12 * b + 6 * c) * a2 + (6 - 2 * b)) / 6.0
    outer = ((-b - 6 * c) * a3 + (6 * b + 30 * c) * a2 + (-12 * b - 48 * c) * a + (8 * b + 24 * c)) / 6.0
    return np.where(a < 1.0, inner, np.where(a < 2.0, outer, 0.0))


def keys_kernel(a: float = -0.5, antialias: bool = False) -> Kernel:
    """
    Cubique de convolution de Keys (1981), paramètre ``a`` (``-0.5`` :
    ordre 3 d'approximation). C'est le cas ``B = 0, C = -a`` de
    Mitchell-Netravali.
    """
    return Kernel(lambda u: _cubic(u, 0.0, -a), 2.0, antialias)


def mitchell_kernel(b: float = 1.0 / 3.0, c: float = 1.0 / 3.0) -> Kernel:
    """Filtre de Mitchell & Netravali (1988), ``B = C = 1/3`` par défaut ; dilaté."""
    return Kernel(lambda u: _cubic(u, b, c), 2.0, True)


# --- Filtres anti-repliement --------------------------------------------------

def gaussian_kernel(sigma: float = 0.5) -> Kernel:
    """
    Gaussienne d'écart-type ``sigma`` **en pixels cibles**, tronquée à
    4 σ. MTF : ``exp(-2 π² σ² f²)``, soit ``exp(-π² σ² / 2)`` au Nyquist
    cible (0,29 pour σ = 0,5).
    """
    if sigma <= 0:
        raise ValueError("sigma doit être > 0")
    return Kernel(lambda u: np.exp(-0.5 * (u / sigma) ** 2), max(4.0 * sigma, 1.0), True)


def lanczos_kernel(a: int = 3) -> Kernel:
    """Lanczos-``a`` : ``sinc(u) sinc(u / a)`` sur ``|u| < a`` (Duchon, 1979)."""
    if a < 1:
        raise ValueError("a doit être >= 1")
    return Kernel(lambda u: np.where(np.abs(u) < a, np.sinc(u) * np.sinc(u / a), 0.0), float(a), True)


def sinc_gaussian_kernel(fwhm: float = 6.0, radius: float = 6.0, cutoff: float = 0.5) -> Kernel:
    """
    Sinus cardinal de coupure ``cutoff`` (cycles par pixel cible, 0,5 =
    Nyquist cible) apodisé par une gaussienne de largeur à mi-hauteur
    ``fwhm`` (pixels cibles), tronqué à ``radius`` pixels cibles.
    """
    if fwhm <= 0 or radius <= 0 or not 0 < cutoff <= 0.5:
        raise ValueError("fwhm et radius doivent être > 0, cutoff dans ]0, 0.5]")
    sigma = fwhm * FWHM_TO_SIGMA

    def phi(u):
        return 2 * cutoff * np.sinc(2 * cutoff * u) * np.exp(-0.5 * (u / sigma) ** 2)

    return Kernel(phi, float(radius), True)


def kaiser_beta(attenuation: float) -> float:
    """
    Paramètre β de la fenêtre de Kaiser pour une atténuation en bande
    coupée de ``attenuation`` dB (formule empirique de Kaiser, 1974).
    """
    a = float(attenuation)
    if a > 50:
        return 0.1102 * (a - 8.7)
    if a >= 21:
        return 0.5842 * (a - 21) ** 0.4 + 0.07886 * (a - 21)
    return 0.0


def kaiser_radius(attenuation: float, transition: float) -> float:
    """
    Demi-support (pixels cibles) d'un filtre de Kaiser atteignant
    ``attenuation`` dB avec une bande de transition de largeur
    ``transition`` (cycles par pixel cible) : longueur
    ``L = (A - 7,95) / (2,285 · 2π Δf)`` (Kaiser, 1974).
    """
    if transition <= 0:
        raise ValueError("transition doit être > 0")
    length = max(attenuation - 7.95, 1.0) / (2.285 * 2.0 * math.pi * transition)
    return max(length / 2.0, 1.0)


def sinc_kaiser_kernel(
    cutoff: float = 0.5,
    attenuation: float = 60.0,
    transition: float = 0.2,
    radius: float | None = None,
) -> Kernel:
    """
    Sinus cardinal fenêtré par Kaiser. ``cutoff`` : fréquence de coupure
    (cycles par pixel cible, gain 1/2) ; ``attenuation`` (dB) et
    ``transition`` (largeur de la bande de transition, cycles par pixel
    cible) fixent β et le support. ``radius`` (pixels cibles) impose le
    support : la transition obtenue s'en déduit.
    """
    if not 0 < cutoff <= 0.5:
        raise ValueError("cutoff doit être dans ]0, 0.5]")
    beta = kaiser_beta(attenuation)
    r = float(radius) if radius is not None else kaiser_radius(attenuation, transition)
    norm = float(i0(beta))

    def phi(u):
        t = np.clip(1.0 - (u / r) ** 2, 0.0, None)
        w = np.where(np.abs(u) < r, i0(beta * np.sqrt(t)) / norm, 0.0)
        return 2 * cutoff * np.sinc(2 * cutoff * u) * w

    return Kernel(phi, r, True)
