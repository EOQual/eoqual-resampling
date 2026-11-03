#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Démarrage rapide : réduire une image, choisir une méthode, simuler un
capteur moins résolu, et comparer les fonctions de transfert.

    python examples/quickstart.py            # affiche les résultats
    python examples/quickstart.py --plot     # + figure Plotly des fonctions de transfert
"""

import sys

import numpy as np

from eoqual.resampling import kernel_metrics, list_methods, resample
from eoqual.resampling.mtf import gaussian_mtf, pixel_mtf, product


def scene(n=600, seed=0):
    """Scène synthétique : texture + bords nets + mire radiale."""
    rng = np.random.default_rng(seed)
    y, x = np.mgrid[0:n, 0:n].astype(np.float64)
    from scipy.ndimage import gaussian_filter

    img = 0.3 * gaussian_filter(rng.random((n, n)), 1.5)
    img += 0.4 * ((x > n / 3) & (y > n / 3) & (x < 2 * n / 3) & (y < 2 * n / 3))
    img += 0.2 * np.cos(np.arctan2(y - n / 2, x - n / 2) * 60)
    return img


def main(plot=False):
    list_methods()
    img = scene()

    # 1. Réduction d'un facteur 3 par la méthode par défaut (sinc fenêtré par Kaiser).
    small = resample(img, 3)
    print(f"facteur 3 : {img.shape} -> {small.shape}, moyenne {img.mean():.4f} -> {small.mean():.4f}")

    # 2. Facteur non entier, ou taille imposée, autre méthode.
    print("facteur 2,5 :", resample(img, 2.5, method="sinc_gaussian").shape)
    print("taille imposée :", resample(img, shape=(150, 200), method="fourier").shape)

    # 3. Image d'un capteur simulé : MTF visée = optique (0,5 au Nyquist) × pixel.
    target = product(gaussian_mtf(0.5), pixel_mtf())
    sim = resample(img, 4, method="mtf_target", target=target)
    print("capteur simulé, facteur 4 :", sim.shape)

    # 4. Compromis de quelques méthodes au facteur 3.
    for m in ("box", "lanczos", "sinc_gaussian", "sinc_kaiser", "fourier"):
        r = kernel_metrics(m, 3)
        print(f"{m:14s} MTF@Nyquist {r['mtf_nyquist']:.2f}  repliement {r['aliasing']:.4f}  "
              f"dépassement {100 * r['overshoot']:.1f} %")

    if plot:
        from eoqual.resampling.plotting import plot_transfer_functions

        plot_transfer_functions(["box", "gaussian", "lanczos", "sinc_gaussian", "sinc_kaiser", "fourier"], 3).show()


if __name__ == "__main__":
    main(plot="--plot" in sys.argv)
