#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Scènes synthétiques communes aux tests (vérité connue)."""
import os

import numpy as np
import pytest

from eoqual.resampling import METHODS_CONFIGS
from eoqual.resampling.mtf import gaussian_mtf

#: Graine des tirages aléatoires, modifiable par EOQUAL_TEST_SEED.
SEED = int(os.environ.get("EOQUAL_TEST_SEED", "0"))

#: Toutes les méthodes du registre.
ALL_METHODS = list(METHODS_CONFIGS)

#: Méthodes sans préfiltrage (replient le spectre).
REFERENCE_METHODS = [m for m, c in METHODS_CONFIGS.items() if c["family"] == "référence"]


def method_params(method: str) -> dict:
    """Paramètres minimaux pour appeler ``method`` (MTF requise pour mtf_target)."""
    return {"target": gaussian_mtf(0.3)} if method == "mtf_target" else {}


@pytest.fixture(scope="session")
def rng():
    return np.random.default_rng(SEED)


@pytest.fixture(scope="session")
def texture(rng):
    """Texture aléatoire lissée 192×192 (spectre large, sans structure)."""
    from scipy.ndimage import gaussian_filter

    return gaussian_filter(rng.random((192, 192)), 0.8)


def blob(shape, cy, cx, sigma):
    """Tache gaussienne de centre (cy, cx) en pixels source."""
    y, x = np.mgrid[0:shape[0], 0:shape[1]].astype(np.float64)
    return np.exp(-((x - cx) ** 2 + (y - cy) ** 2) / (2.0 * sigma ** 2))


def centroid(image):
    y, x = np.mgrid[0:image.shape[0], 0:image.shape[1]]
    s = image.sum()
    return (image * y).sum() / s, (image * x).sum() / s


def sinusoid_gain(method, factor, f_target, n=600, **params):
    """
    Gain mesuré sur une sinusoïde 1D de fréquence ``f_target`` (cycles
    par pixel cible) : amplitude de la sortie à la fréquence (éventuellement
    repliée) attendue, rapportée à 1, avec son signe.
    """
    from eoqual.resampling import resample
    from eoqual.resampling.geometry import target_positions

    x = np.arange(n, dtype=np.float64)
    f_src = f_target / factor
    img = np.broadcast_to(np.cos(2 * np.pi * f_src * x), (8, n))
    out = resample(img, (1.0, factor), method=method, **params)[4]
    xt = target_positions(out.size, factor)
    keep = slice(25, -25)
    a = np.column_stack([np.cos(2 * np.pi * f_src * xt), np.sin(2 * np.pi * f_src * xt)])[keep]
    coef, *_ = np.linalg.lstsq(a, out[keep], rcond=None)
    return float(coef[0]), float(coef[1])
