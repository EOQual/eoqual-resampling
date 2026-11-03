#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Comportement fréquentiel : fonction de transfert, repliement, exactitude spectrale, bruit."""
import numpy as np
import pytest

from conftest import ALL_METHODS, REFERENCE_METHODS, method_params, sinusoid_gain
from eoqual.resampling import kernel_metrics, resample, transfer_function
from eoqual.resampling.geometry import target_positions
from eoqual.resampling.mtf import diffraction_mtf, gaussian_mtf, pixel_mtf, product


@pytest.mark.parametrize("method", ALL_METHODS)
@pytest.mark.parametrize("factor", [2, 3])
def test_transfer_function_matches_measured_gain(method, factor):
    """Le gain annoncé par transfer_function est celui mesuré sur une sinusoïde (repliée ou non)."""
    params = method_params(method)
    for f in (0.1, 0.3, 0.45, 0.7):
        c, s = sinusoid_gain(method, factor, f, **params)
        _, g = transfer_function(method, factor, freqs=np.array([f]), **params)
        if method == "nearest" and factor % 2 == 0:
            # Seul le plus proche voisin n'est pas symétrique : pour un facteur
            # pair, le centre cible tombe entre deux pixels (décalage d'un
            # demi-pixel source) ; transfer_function donne alors le module.
            assert abs(np.hypot(c, s) - g[0]) < 2e-3, (f, c, s, g[0])
        else:
            assert abs(c - g[0]) < 2e-3 and abs(s) < 2e-3, (f, c, s, g[0])


@pytest.mark.parametrize("method", REFERENCE_METHODS)
def test_reference_methods_alias_fully(method):
    """Sans préfiltrage, une fréquence au-delà du Nyquist cible se replie intégralement."""
    c, _ = sinusoid_gain(method, 3, 0.7)
    assert c > 0.95


# fourier : le résidu de 1e-6 vient de la composante lisse (sinusoïde non
# périodique sur l'image), réduite par un Lanczos-3.
@pytest.mark.parametrize("method,max_alias", [("sinc_kaiser", 2e-3), ("fourier", 1e-5)])
def test_antialias_rejects_above_nyquist(method, max_alias):
    c, s = sinusoid_gain(method, 3, 0.7)
    assert np.hypot(c, s) < max_alias


# Lanczos-3 a une ondulation de +1,2 % en bande passante (propriété connue).
@pytest.mark.parametrize("method,tol", [("sinc_gaussian", 0.005), ("sinc_kaiser", 0.005),
                                        ("lanczos", 0.02), ("fourier", 1e-6)])
def test_passband_is_preserved(method, tol):
    c, _ = sinusoid_gain(method, 3, 0.2)
    assert abs(c - 1.0) < tol


@pytest.mark.parametrize("n,factor", [(240, 3), (240, 2.5), (241, 3.7), (250, 3)])
def test_fourier_is_exact_on_band_limited_periodic_image(n, factor, rng):
    """Image périodique à bande limitée : la sortie est l'échantillonnage exact aux positions cibles."""
    fmax = int(0.35 * n / factor)
    comps = [(*rng.integers(-fmax, fmax + 1, 2), rng.normal(), rng.uniform(0, 2 * np.pi)) for _ in range(10)]

    def scene(y, x):
        return sum(a * np.cos(2 * np.pi * (my * y + mx * x) / n + ph) for my, mx, a, ph in comps)

    yy, xx = np.mgrid[0:n, 0:n].astype(np.float64)
    out = resample(scene(yy, xx), factor, method="fourier", boundary="wrap")
    yt, xt = target_positions(out.shape[0], factor), target_positions(out.shape[1], factor)
    np.testing.assert_allclose(out, scene(yt[:, None], xt[None, :]), atol=1e-10)


def test_fourier_non_periodic_image_has_no_edge_artefact():
    """Grâce à la décomposition p + s, une rampe (non périodique) reste une rampe."""
    yy, xx = np.mgrid[0:240, 0:240].astype(np.float64)
    ramp = 0.01 * xx + 0.02 * yy
    out = resample(ramp, 3, method="fourier")
    yt, xt = target_positions(80, 3), target_positions(80, 3)
    expected = 0.01 * xt[None, :] + 0.02 * yt[:, None]
    assert np.abs(out - expected)[3:-3, 3:-3].max() < 2e-3


def test_mtf_target_reproduces_target_mtf_and_its_aliasing():
    """La MTF visée est appliquée telle quelle, y compris au-delà du Nyquist (repliement physique)."""
    target = product(gaussian_mtf(0.3), pixel_mtf())
    for f in (0.2, 0.45, 0.7):
        c, _ = sinusoid_gain("mtf_target", 3, f, target=target)
        assert abs(c - float(target(np.array(0.0), np.array(f)))) < 2e-3


def test_mtf_target_compensates_source_mtf():
    """Avec une MTF source déclarée, la sortie a la MTF cible, pas le produit des deux."""
    source = gaussian_mtf(0.6)                    # MTF de l'image d'entrée (cycles / pixel source)
    target = gaussian_mtf(0.3)
    x = np.arange(600, dtype=np.float64)
    f_t = 0.3
    f_s = f_t / 3
    attenuation = float(source(np.array(0.0), np.array(f_s)))
    img = np.broadcast_to(attenuation * np.cos(2 * np.pi * f_s * x), (8, 600))
    out = resample(img, (1, 3), method="mtf_target", target=target, source=source)[4]
    xt = target_positions(out.size, 3)
    a = np.column_stack([np.cos(2 * np.pi * f_s * xt), np.sin(2 * np.pi * f_s * xt)])[25:-25]
    coef, *_ = np.linalg.lstsq(a, out[25:-25], rcond=None)
    assert abs(coef[0] - float(target(np.array(0.0), np.array(f_t)))) < 2e-3


def test_mtf_models():
    zero = np.array(0.0)
    assert abs(float(gaussian_mtf(0.3)(zero, np.array(0.5))) - 0.3) < 1e-12
    assert abs(float(pixel_mtf()(zero, np.array(0.5))) - 2 / np.pi) < 1e-12
    d = diffraction_mtf(0.8)
    assert float(d(zero, zero)) == pytest.approx(1.0)
    assert float(d(zero, np.array(0.8))) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("method", ["box", "gaussian", "sinc_kaiser", "fourier"])
def test_noise_variance_ratio_is_predicted(method, rng):
    noise = rng.normal(size=(600, 600))
    out = resample(noise, 3, method=method, boundary="wrap")
    predicted = kernel_metrics(method, 3)["noise_variance_ratio"]
    assert abs(out.var() / predicted - 1.0) < 0.08


def test_kernel_metrics_ordering():
    """Ordre attendu du compromis : repliement, transition, oscillations."""
    m = {name: kernel_metrics(name, 3) for name in ("box", "gaussian", "lanczos", "sinc_kaiser", "fourier")}
    assert m["fourier"]["aliasing"] < 1e-12
    assert m["sinc_kaiser"]["aliasing"] < m["lanczos"]["aliasing"] < m["box"]["aliasing"]
    assert m["sinc_kaiser"]["transition"] < m["lanczos"]["transition"] < m["gaussian"]["transition"]
    assert m["gaussian"]["overshoot"] < 1e-6 and m["box"]["overshoot"] < 1e-6
    assert m["sinc_kaiser"]["overshoot"] > m["lanczos"]["overshoot"]
    assert abs(m["box"]["noise_variance_ratio"] - 1 / 9) < 1e-3


def test_kaiser_radius_controls_transition():
    wide = kernel_metrics("sinc_kaiser", 3, transition=0.3)["transition"]
    narrow = kernel_metrics("sinc_kaiser", 3, transition=0.1)["transition"]
    assert narrow < wide
