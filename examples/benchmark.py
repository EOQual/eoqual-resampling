#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Benchmark des méthodes de réduction : génère les tableaux de BENCHMARK.md.

    python examples/benchmark.py > /tmp/benchmark_tables.md

1. Indicateurs de la fonction de transfert (``kernel_metrics``) pour les
   facteurs 2, 3 et 4.
2. Mire de Fresnel (« zone plate », fréquence locale croissant du centre
   vers les coins jusqu'au Nyquist source) réduite d'un facteur 3 :
   contraste conservé dans l'anneau où la fréquence cible est < 0,35, et
   contraste résiduel (repliement) dans l'anneau où elle est > 0,6 — une
   réduction idéale y donnerait une image uniforme.
3. Temps de calcul sur une image 2048 × 2048 (facteur 3 et 2,5).

Tout est synthétique et déterministe.
"""

import time

import numpy as np

from eoqual.resampling import METHODS_CONFIGS, kernel_metrics, resample
from eoqual.resampling.geometry import target_positions
from eoqual.resampling.mtf import gaussian_mtf, pixel_mtf, product

#: Paramètres utilisés pour mtf_target dans le benchmark : optique
#: gaussienne (0,5 au Nyquist) × pixel, soit 0,32 au Nyquist cible.
MTF_TARGET = {"target": product(gaussian_mtf(0.5), pixel_mtf())}

VARIANTS = [(m, MTF_TARGET if m == "mtf_target" else {}) for m in METHODS_CONFIGS]
VARIANTS += [
    ("sinc_kaiser", {"transition": 0.1}),
    ("sinc_kaiser", {"radius": 6.0}),
    ("fourier", {"passband": 0.45, "stopband": 0.55}),
]


def label(method, params):
    if not params or method == "mtf_target":
        return f"`{method}`"
    return f"`{method}` ({', '.join(f'{k}={v}' for k, v in params.items())})"


def fmt(x, nd=3):
    return "—" if not np.isfinite(x) else f"{x:.{nd}f}"


def table_metrics(factor):
    print(f"\n### Facteur {factor}\n")
    print("| Méthode | MTF au Nyquist | f90 | f10 | transition | énergie repliée / conservée | variance du bruit | dépassement |")
    print("|---|---|---|---|---|---|---|---|")
    for m, p in VARIANTS:
        r = kernel_metrics(m, factor, **p)
        print(f"| {label(m, p)} | {fmt(r['mtf_nyquist'])} | {fmt(r['f90'])} | {fmt(r['f10'])} | "
              f"{fmt(r['transition'])} | {fmt(r['aliasing'], 4)} | {fmt(r['noise_variance_ratio'])} | "
              f"{fmt(100 * r['overshoot'], 1)} % |")


def zone_plate(n):
    """Mire de Fresnel : fréquence locale 0,5 · r / r_max cycle par pixel source."""
    c = (n - 1) / 2.0
    y, x = np.mgrid[0:n, 0:n].astype(np.float64)
    r_max = c * np.sqrt(2.0)
    phase = np.pi * ((x - c) ** 2 + (y - c) ** 2) / (2.0 * r_max)
    return 0.5 + 0.5 * np.cos(phase), c, r_max


def table_zone_plate(factor=3, n=1536):
    img, c, r_max = zone_plate(n)
    print(f"\n### Mire de Fresnel {n} × {n}, facteur {factor}\n")
    print("Contraste RMS rapporté à celui de la mire (0,354) : conservé en bande passante "
          "(idéal 1), résiduel au-delà du Nyquist cible (idéal 0).\n")
    print("f : plus grande composante de la fréquence locale, en cycles par pixel cible.\n")
    print("| Méthode | contraste conservé (f < 0,35) | repliement résiduel (f > 0,6) |")
    print("|---|---|---|")
    for m, p in VARIANTS:
        out = resample(img, factor, method=m, **p)
        yt = target_positions(out.shape[0], factor)
        xt = target_positions(out.shape[1], factor)
        dy, dx = yt[:, None] - c, xt[None, :] - c
        # Fréquence locale : vecteur radial de module 0,5 r / r_max (cycles par
        # pixel source). La grille cible est carrée : une composante est
        # représentable si elle est sous 0,5 par axe — on classe donc par la
        # plus grande des deux composantes, pas par le module.
        f_t = factor * 0.5 * np.maximum(np.abs(dy), np.abs(dx)) / r_max
        ref = 0.5 / np.sqrt(2.0)
        inner = (f_t < 0.35) & (f_t > 0.05)
        outer = (f_t > 0.6) & (np.hypot(dy, dx) < c - 2 * factor)      # loin des bords
        keep = np.std(out[inner] - 0.5) / ref
        alias = np.std(out[outer] - 0.5) / ref
        print(f"| {label(m, p)} | {keep:.3f} | {alias:.3f} |")


def table_timing():
    rng = np.random.default_rng(0)
    img = rng.random((2048, 2048))
    print("\n### Temps de calcul, image 2048 × 2048 (meilleur de 3 essais)\n")
    print("| Méthode | facteur 3 (s) | facteur 2,5 (s) |")
    print("|---|---|---|")
    for m, p in VARIANTS[: len(METHODS_CONFIGS)]:
        times = []
        for k in (3, 2.5):
            best = np.inf
            for _ in range(3):
                t0 = time.perf_counter()
                resample(img, k, method=m, **p)
                best = min(best, time.perf_counter() - t0)
            times.append(best)
        print(f"| {label(m, p)} | {times[0]:.2f} | {times[1]:.2f} |")


if __name__ == "__main__":
    print("## Indicateurs de la fonction de transfert")
    for k in (2, 3, 4):
        table_metrics(k)
    print("\n## Mire de Fresnel")
    table_zone_plate()
    print("\n## Temps de calcul")
    table_timing()
