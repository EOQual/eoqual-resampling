#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Contrat d'appel, géométrie et conservation des niveaux, pour toutes les méthodes."""
import numpy as np
import pytest

from conftest import ALL_METHODS, REFERENCE_METHODS, blob, centroid, method_params
from eoqual.resampling import list_methods, resample


@pytest.mark.parametrize("method", ALL_METHODS)
@pytest.mark.parametrize("factor", [2, 2.5, 3, (2, 3.5)])
def test_contract(method, factor, texture):
    out = resample(texture, factor, method=method, **method_params(method))
    k = np.broadcast_to(factor, 2)
    assert out.shape == (int(192 / k[0]), int(192 / k[1]))
    assert out.dtype == np.float64
    assert np.all(np.isfinite(out))


@pytest.mark.parametrize("method", ALL_METHODS)
@pytest.mark.parametrize("factor", [2, 2.5, 3.7])
def test_uniform_image_is_preserved(method, factor):
    out = resample(np.full((150, 170), 3.25), factor, method=method, **method_params(method))
    np.testing.assert_allclose(out, 3.25, atol=1e-12)


@pytest.mark.parametrize("method", [m for m in ALL_METHODS if m not in REFERENCE_METHODS])
@pytest.mark.parametrize("grid", ["center", "origin"])
@pytest.mark.parametrize("factor", [2, 2.5, 3, 4])
def test_no_geometric_shift(method, grid, factor):
    """Le centroïde d'une tache tombe exactement à la position attendue sur la grille cible."""
    cy, cx = 118.7, 121.3
    img = blob((241, 241), cy, cx, 6.0)
    out = resample(img, factor, method=method, grid=grid, **method_params(method))
    oy, ox = centroid(out)
    if grid == "center":
        ey, ex = (cy + 0.5) / factor - 0.5, (cx + 0.5) / factor - 0.5
    else:
        ey, ex = cy / factor, cx / factor
    assert abs(oy - ey) < 2e-3 and abs(ox - ex) < 2e-3


def test_shape_argument_sets_exact_factor(texture):
    out = resample(texture, shape=(70, 50), method="sinc_kaiser")
    assert out.shape == (70, 50)


@pytest.mark.parametrize("method", ["sinc_gaussian", "fourier"])
def test_channel_axis(method, texture):
    rgb = np.stack([texture, 2 * texture, texture + 1], axis=-1)
    out = resample(rgb, 3, method=method, channel_axis=-1)
    assert out.shape == (64, 64, 3)
    np.testing.assert_allclose(out[..., 1], 2 * out[..., 0], atol=1e-12)
    np.testing.assert_allclose(out[..., 2], out[..., 0] + 1, atol=1e-12)


def test_box_is_block_mean_for_integer_factor(texture):
    k = 4
    out = resample(texture, k, method="box")
    blocks = texture.reshape(48, k, 48, k).mean(axis=(1, 3))
    np.testing.assert_allclose(out, blocks, atol=1e-12)


def test_box_non_integer_factor_preserves_mean(rng):
    img = rng.random((150, 150))
    out = resample(img, 2.5, method="box")      # 150 / 2.5 = 60 : emprise exacte
    assert abs(out.mean() - img.mean()) < 1e-12


def test_origin_grid_nearest_is_plain_decimation(texture):
    out = resample(texture, 3, method="nearest", grid="origin")
    np.testing.assert_array_equal(out, texture[::3, ::3])


@pytest.mark.parametrize("boundary", ["symmetric", "reflect", "edge", "wrap"])
def test_boundaries(boundary, texture):
    out = resample(texture, 3, method="lanczos", boundary=boundary)
    assert out.shape == (64, 64) and np.all(np.isfinite(out))


def test_errors(texture):
    with pytest.raises(ValueError, match="inconnue"):
        resample(texture, 2, method="nope")
    with pytest.raises(NotImplementedError, match="agrandissement"):
        resample(texture, 0.5)
    with pytest.raises(ValueError, match="un seul"):
        resample(texture, 2, shape=(10, 10))
    with pytest.raises(ValueError, match="un seul"):
        resample(texture)
    with pytest.raises(TypeError, match="inconnu"):
        resample(texture, 2, method="lanczos", sigma=1.0)
    with pytest.raises(ValueError, match="target"):
        resample(texture, 2, method="mtf_target")
    with pytest.raises(ValueError, match="grid"):
        resample(texture, 2, grid="corner")
    with pytest.raises(ValueError, match="boundary"):
        resample(texture, 2, boundary="zero")
    bad = texture.copy()
    bad[3, 3] = np.nan
    with pytest.raises(ValueError, match="non finies"):
        resample(bad, 2)
    with pytest.raises(ValueError, match="2D"):
        resample(np.zeros((4, 4, 3)), 2)
    with pytest.raises(ValueError, match="trop grand"):
        resample(np.zeros((4, 4)), 5)


def test_list_methods_runs(capsys):
    list_methods()
    assert "sinc_gaussian" in capsys.readouterr().out
