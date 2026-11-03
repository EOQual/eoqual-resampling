"""
Construction de chaque méthode à partir de ses paramètres : noyau 1D
(méthodes séparables) ou réponse en fréquence (méthodes spectrales).

Les paramètres sont en unités **cibles** (pixels cibles, cycles par
pixel cible) : un même réglage garde le même sens quel que soit le
facteur de réduction.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import NDArray

from .backends.separable import BOUNDARIES, apply_separable, box_matrix, resampling_matrix
from .backends.spectral import check_response_args, raised_cosine, spectral_resample
from .config import METHODS_CONFIGS
from .kernels import (
    Kernel,
    gaussian_kernel,
    keys_kernel,
    lanczos_kernel,
    linear_kernel,
    mitchell_kernel,
    nearest_kernel,
    sinc_gaussian_kernel,
    sinc_kaiser_kernel,
)

__all__ = ["check_method", "build_kernel", "build_response", "resample_2d"]


def check_method(method: str) -> dict:
    """Configuration de ``method`` ; ``ValueError`` si elle est inconnue."""
    if method not in METHODS_CONFIGS:
        raise ValueError(
            f"Méthode '{method}' inconnue. Méthodes disponibles : {list(METHODS_CONFIGS)}"
        )
    return METHODS_CONFIGS[method]


def _check_params(method: str, params: dict[str, Any]) -> None:
    allowed = set(METHODS_CONFIGS[method]["params"])
    unknown = set(params) - allowed
    if unknown:
        raise TypeError(
            f"paramètre(s) inconnu(s) pour '{method}' : {sorted(unknown)} ; "
            f"paramètres acceptés : {sorted(allowed) or 'aucun'}"
        )


def build_kernel(method: str, **params) -> Kernel | str:
    """Noyau d'une méthode séparable (``"box"`` pour l'agrégation exacte)."""
    check_method(method)
    _check_params(method, params)
    if method == "nearest":
        return nearest_kernel()
    if method == "bilinear":
        return linear_kernel(antialias=False)
    if method == "bicubic":
        return keys_kernel(params.get("a", -0.5), antialias=False)
    if method == "box":
        return "box"
    if method == "gaussian":
        return gaussian_kernel(params.get("sigma", 0.5))
    if method == "bicubic_aa":
        return keys_kernel(params.get("a", -0.5), antialias=True)
    if method == "mitchell":
        return mitchell_kernel(params.get("b", 1.0 / 3.0), params.get("c", 1.0 / 3.0))
    if method == "lanczos":
        return lanczos_kernel(int(params.get("a", 3)))
    if method == "sinc_gaussian":
        return sinc_gaussian_kernel(params.get("fwhm", 6.0), params.get("radius", 6.0),
                                    params.get("cutoff", 0.5))
    if method == "sinc_kaiser":
        return sinc_kaiser_kernel(params.get("cutoff", 0.5), params.get("attenuation", 60.0),
                                  params.get("transition", 0.2), params.get("radius"))
    raise ValueError(f"'{method}' n'est pas une méthode séparable")


def build_response(method: str, factors: tuple[float, float], **params):
    """
    Réponse ``G(fy, fx)`` (cycles par pixel cible) d'une méthode
    spectrale, et l'option ``periodic``.
    """
    check_method(method)
    _check_params(method, params)
    periodic = bool(params.get("periodic", False))
    if method == "fourier":
        fp, fs = params.get("passband", 0.4), params.get("stopband", 0.5)
        check_response_args(fp, fs)

        def response(fy, fx):
            return raised_cosine(fy, fp, fs) * raised_cosine(fx, fp, fs)

        return response, periodic

    if method == "mtf_target":
        target = params.get("target")
        if target is None or not callable(target):
            raise ValueError(
                "mtf_target : 'target' (MTF visée, fonction (fy, fx) en cycles par pixel "
                "cible) est requis — voir eoqual.resampling.mtf"
            )
        source = params.get("source")
        floor = float(params.get("floor", 0.05))
        if not 0 < floor <= 1:
            raise ValueError("floor doit être dans ]0, 1]")
        ky, kx = factors

        def response(fy, fx):
            t = np.asarray(target(fy, fx), dtype=np.float64)
            if source is None:
                return t
            s = np.asarray(source(fy / ky, fx / kx), dtype=np.float64)
            return t / np.maximum(s, floor)

        return response, periodic
    raise ValueError(f"'{method}' n'est pas une méthode spectrale")


def resample_2d(
    image: NDArray[np.float64],
    method: str,
    factors: tuple[float, float],
    out_shape: tuple[int, int],
    grid: str = "center",
    boundary: str = "symmetric",
    **params,
) -> NDArray[np.float64]:
    """Réduit une image 2D float64 (pas de contrôle des entrées ici)."""
    cfg = check_method(method)
    if boundary not in BOUNDARIES:
        raise ValueError(f"boundary='{boundary}' inconnu ; valeurs possibles : {BOUNDARIES}")
    if cfg["backend"] == "spectral":
        response, periodic = build_response(method, factors, **params)
        return spectral_resample(image, factors, out_shape, grid, response,
                                 periodic=periodic or boundary == "wrap")
    kernel = build_kernel(method, **params)
    (ny, nx), (my, mx), (ky, kx) = image.shape, out_shape, factors
    if kernel == "box":
        w_r = box_matrix(ny, my, ky, grid, boundary)
        w_c = box_matrix(nx, mx, kx, grid, boundary)
    else:
        w_r = resampling_matrix(ny, my, ky, kernel, grid, boundary)
        w_c = resampling_matrix(nx, mx, kx, kernel, grid, boundary)
    return apply_separable(image, w_r, w_c)
