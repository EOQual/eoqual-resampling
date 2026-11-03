"""
Diagnostic des méthodes : fonction de transfert et indicateurs du
compromis repliement / flou / oscillations / bruit.

Toutes les fréquences sont en cycles par pixel **cible** : le Nyquist
cible vaut 0,5, le Nyquist source ``k / 2``.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from .backends.separable import box_matrix, resampling_matrix
from .geometry import target_positions
from .methods import build_kernel, build_response, check_method
from .runner import resample

__all__ = ["transfer_function", "kernel_metrics", "step_response"]


def _separable_taps(method: str, factor: float, **params):
    """Lignes de la matrice d'un axe, loin des bords, et positions cibles."""
    kernel = build_kernel(method, **params)
    n_out = 64
    n_in = int(np.ceil(n_out * factor))
    if kernel == "box":
        w = box_matrix(n_in, n_out, factor, "center", "symmetric")
    else:
        w = resampling_matrix(n_in, n_out, factor, kernel, "center", "symmetric")
    x = target_positions(n_out, factor, "center")
    rows = range(n_out // 4, 3 * n_out // 4)
    return w, x, rows


def transfer_function(
    method: str,
    factor: float,
    freqs: NDArray[np.float64] | None = None,
    **params,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Réponse en fréquence 1D (le long d'un axe) de ``method`` pour un
    facteur ``factor``.

    Méthodes séparables : transformée de Fourier à temps discret des
    poids réellement appliqués, autour de chaque position cible, moyennée
    sur les positions (pour un facteur non entier, les poids varient d'un
    pixel cible à l'autre). Méthodes spectrales : la réponse imposée,
    ``G(0, f)``.

    Returns
    -------
    freqs, gain
        Fréquences (cycles par pixel cible, de 0 au Nyquist source
        ``factor / 2`` par défaut) et gain réel (signé : un noyau à lobes
        négatifs peut avoir un gain négatif). Pour une méthode qui décale
        l'image (plus proche voisin, facteur pair, grille ``"center"``),
        le module du gain.
    """
    cfg = check_method(method)
    if freqs is None:
        freqs = np.linspace(0.0, factor / 2.0, 801)
    freqs = np.asarray(freqs, dtype=np.float64)
    if cfg["backend"] == "spectral":
        response, _ = build_response(method, (factor, factor), **params)
        return freqs, np.asarray(response(np.zeros_like(freqs), freqs), dtype=np.float64)

    w, x, rows = _separable_taps(method, factor, **params)
    f_src = freqs / factor
    acc = np.zeros(freqs.size, dtype=np.complex128)
    for i in rows:
        start, end = w.indptr[i], w.indptr[i + 1]
        j, wij = w.indices[start:end], w.data[start:end]
        acc += np.exp(-2j * np.pi * np.outer(f_src, j - x[i])) @ wij
    acc /= len(rows)
    # Noyau symétrique : réponse réelle (signée). Sinon (plus proche voisin,
    # facteur pair, grille "center") la méthode décale d'un demi-pixel
    # source : on renvoie le module.
    if np.max(np.abs(acc.imag)) < 1e-9:
        return freqs, acc.real
    return freqs, np.abs(acc)


def step_response(method: str, factor: float, n_phases: int = 8, **params) -> NDArray[np.float64]:
    """
    Profils réduits d'un échelon 0 → 1 placé à ``n_phases`` positions
    sous-pixel (une ligne par position), pour mesurer les oscillations.
    """
    n = int(np.ceil(48 * factor))
    rows = int(np.ceil(2 * factor))
    profiles = []
    for ph in np.linspace(0.0, 1.0, n_phases, endpoint=False):
        edge = n / 2 + ph * factor
        img = np.broadcast_to((np.arange(n) >= edge).astype(np.float64), (rows, n))
        out = resample(img, (1.0, factor), method=method, boundary="edge", **params)
        profiles.append(out[out.shape[0] // 2])
    return np.array(profiles)


def kernel_metrics(method: str, factor: float, **params) -> dict:
    """
    Indicateurs de ``method`` pour un facteur ``factor`` :

    - ``mtf_nyquist`` : gain au Nyquist cible ;
    - ``f90`` / ``f10`` : fréquence où le gain passe sous 0,9 / 0,1 ;
      ``transition`` = ``f10 - f90`` ;
    - ``aliasing`` : énergie au-delà du Nyquist cible rapportée à
      l'énergie en deçà (spectre d'entrée blanc), soit la part de signal
      qui se replie ;
    - ``noise_variance_ratio`` : variance d'un bruit blanc après réduction
      divisée par sa variance avant (image 2D) ;
    - ``overshoot`` : dépassement maximal d'un échelon unité (pire
      position sous-pixel), mesure des oscillations près des bords.
    """
    f, g = transfer_function(method, factor, **params)
    a = np.abs(g)
    df = f[1] - f[0]
    below = f <= 0.5
    above = f > 0.5

    def first_below(level):
        idx = np.nonzero(a < level)[0]
        return float(f[idx[0]]) if idx.size else float("nan")

    f90, f10 = first_below(0.9), first_below(0.1)
    e_in = np.sum(g[below] ** 2) * df
    e_out = np.sum(g[above] ** 2) * df
    # Variance d'un bruit blanc : ∫ G² sur la bande source [-1/2, 1/2]
    # (cycles par pixel source), soit (2/k) ∫_0^{k/2} G² en fréquences cibles.
    noise_1d = 2.0 / factor * np.sum(g ** 2) * df
    steps = step_response(method, factor, **params)
    return {
        "mtf_nyquist": float(np.interp(0.5, f, a)),
        "f90": f90,
        "f10": f10,
        "transition": f10 - f90,
        "aliasing": float(e_out / e_in) if e_in > 0 else float("nan"),
        "noise_variance_ratio": float(noise_1d ** 2),
        "overshoot": max(float(steps.max() - 1.0), float(-steps.min()), 0.0) + 0.0,
    }
