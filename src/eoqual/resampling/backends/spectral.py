"""
Rééchantillonnage spectral : filtrage par une réponse en fréquence
choisie, puis évaluation **exacte** de l'image filtrée aux positions
cibles.

1. L'image est décomposée en ``u = p + s`` (composante périodique +
   composante lisse, Moisan 2011, ``eoqual.core.fourier.perdecomp``) :
   la composante périodique ne porte plus le saut de bord qui donnerait la
   « croix » spectrale et des oscillations sur les bords de l'image.
2. Le spectre de ``p`` est multiplié par la réponse ``G(fy, fx)``
   (fréquences en cycles par pixel **cible**).
3. Le résultat est l'interpolant trigonométrique de ``p`` filtré, évalué
   aux positions cibles ``x_i`` (voir :mod:`eoqual.resampling.geometry`) :

   - si le pas de la grille cible vaut exactement ``N / M`` (facteur issu
     de ``shape=``, ou ``M·k`` entier — l'image est alors recadrée à
     l'emprise ``M·k`` des pixels cibles), l'échantillonnage replie les
     fréquences modulo ``M`` : spectre replié puis FFT inverse de taille
     ``M``, en ``O(N² log N)`` ;
   - sinon, l'interpolant est évalué par une transformée de Fourier non
     uniforme séparable (deux produits matriciels, exacts, plus lents).

   La fréquence de Nyquist d'une taille paire est partagée par moitié
   entre ``+N/2`` et ``-N/2`` (interpolant réel).
4. La composante lisse ``s`` (très basses fréquences) est réduite par un
   Lanczos-3 séparable, de réponse ≈ 1 à ces fréquences.

Ce schéma réalise **toute** réponse ``G`` : un gabarit coupant au
Nyquist cible (aucun repliement), ou une MTF d'instrument (qui replie
comme le ferait le vrai capteur).
"""

from __future__ import annotations

import math
from typing import Callable

import numpy as np
from numpy.typing import NDArray
from scipy import sparse

from eoqual.core.fourier import perdecomp

from ..geometry import target_positions
from ..kernels import lanczos_kernel
from .separable import apply_separable, resampling_matrix

__all__ = ["Response", "spectral_resample", "raised_cosine"]

#: Réponse en fréquence : ``G(fy, fx)`` en cycles par pixel cible
#: (tableaux de même forme), à valeurs réelles.
Response = Callable[[NDArray[np.float64], NDArray[np.float64]], NDArray[np.float64]]

# Taille des blocs de lignes de la transformée non uniforme (mémoire).
_CHUNK = 256


def raised_cosine(f: NDArray[np.float64], passband: float, stopband: float) -> NDArray[np.float64]:
    """1 jusqu'à ``passband``, transition en cosinus surélevé, 0 à partir de ``stopband``."""
    a = np.abs(f)
    if stopband <= passband:
        return (a <= passband).astype(np.float64)
    t = np.clip((a - passband) / (stopband - passband), 0.0, 1.0)
    return 0.5 * (1.0 + np.cos(np.pi * t))


def _extended_indices(n: int) -> tuple[NDArray[np.int64], NDArray[np.int64], NDArray[np.float64]]:
    """
    Fréquences entières ``m`` de l'interpolant réel de taille ``n``,
    l'indice FFT correspondant et le poids (1/2 pour les deux moitiés du
    Nyquist d'une taille paire).
    """
    if n % 2 == 0:
        m = np.arange(-n // 2, n // 2 + 1)
        w = np.ones(m.size)
        w[0] = w[-1] = 0.5
    else:
        m = np.arange(-(n // 2), n // 2 + 1)
        w = np.ones(m.size)
    return m, np.mod(m, n), w


def _fold_matrix(m: NDArray[np.int64], n_out: int) -> sparse.csr_matrix:
    """Matrice 0/1 qui replie les fréquences ``m`` modulo ``n_out``."""
    rows = np.mod(m, n_out)
    return sparse.csr_matrix((np.ones(m.size), (rows, np.arange(m.size))), shape=(n_out, m.size))


def _nudft_rows(x: NDArray[np.float64], m: NDArray[np.int64], n: int) -> NDArray[np.complex128]:
    """Matrice ``exp(2iπ m x / n) / n`` (positions × fréquences)."""
    return np.exp(2j * np.pi * np.outer(x, m) / n) / n


def spectral_resample(
    image: NDArray[np.float64],
    factors: tuple[float, float],
    out_shape: tuple[int, int],
    grid: str,
    response: Response,
    periodic: bool = False,
) -> NDArray[np.float64]:
    """
    Réduit ``image`` (2D, float64) par la réponse ``response`` (voir
    l'en-tête du module).

    ``periodic=True`` traite l'image comme périodique (pas de
    décomposition ``p + s``) : utile pour une image réellement périodique
    (tests, simulation sur un tore).
    """
    ky, kx = factors
    my, mx = out_shape

    # Recadrage à l'emprise des pixels cibles quand M·k est entier : le pas
    # cible vaut alors exactement N'/M et le repliement est exact et rapide.
    # Jamais pour une image périodique : recadrer casserait sa période.
    crop = []
    for n, k, mo in ((image.shape[0], ky, my), (image.shape[1], kx, mx)):
        span = mo * k
        exact = abs(span - round(span)) < 1e-9 and round(span) <= n
        crop.append(int(round(span)) if exact and not periodic else n)
    u = image[: crop[0], : crop[1]]

    if periodic:
        p, s = u, None
    else:
        p, s = perdecomp(u)

    ny, nx = p.shape
    spec = np.fft.fft2(p)
    m_y, idx_y, w_y = _extended_indices(ny)
    m_x, idx_x, w_x = _extended_indices(nx)
    ext = spec[np.ix_(idx_y, idx_x)] * np.outer(w_y, w_x)

    # Réponse sur la grille des fréquences, en cycles par pixel cible.
    fy = (m_y / ny * ky)[:, None]
    fx = (m_x / nx * kx)[None, :]
    gain = np.asarray(response(np.broadcast_to(fy, ext.shape), np.broadcast_to(fx, ext.shape)),
                      dtype=np.float64)
    ext = ext * gain

    y = target_positions(my, ky, grid)
    x = target_positions(mx, kx, grid)
    fast_y = abs(ny / my - ky) < 1e-9
    fast_x = abs(nx / mx - kx) < 1e-9

    if fast_y and fast_x:
        # Valeur en x_i = i·N/M + δ : phase exp(2iπ m δ / N), repli modulo M.
        ext = ext * np.outer(np.exp(2j * np.pi * m_y * y[0] / ny), np.exp(2j * np.pi * m_x * x[0] / nx))
        folded = _fold_matrix(m_y, my) @ ext
        folded = (_fold_matrix(m_x, mx) @ folded.T).T
        out = np.real(np.fft.ifft2(folded)) * (my / ny) * (mx / nx)
    else:
        a_x = _nudft_rows(x, m_x, nx)                      # (mx, nx+)
        tmp = ext @ a_x.T                                  # (ny+, mx)
        out = np.empty((my, mx))
        for start in range(0, my, _CHUNK):
            a_y = _nudft_rows(y[start:start + _CHUNK], m_y, ny)
            out[start:start + _CHUNK] = np.real(a_y @ tmp)

    if s is not None:
        kern = lanczos_kernel(3)
        w_r = resampling_matrix(s.shape[0], my, ky, kern, grid, "symmetric")
        w_c = resampling_matrix(s.shape[1], mx, kx, kern, grid, "symmetric")
        out = out + apply_separable(s, w_r, w_c)
    return out


def check_response_args(passband: float, stopband: float) -> None:
    """Contrôle des bornes d'un gabarit en cosinus surélevé."""
    if not (0 < passband <= stopband) or not math.isfinite(stopband):
        raise ValueError("il faut 0 < passband <= stopband (cycles par pixel cible)")
