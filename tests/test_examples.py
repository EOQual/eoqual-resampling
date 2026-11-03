#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Les exemples s'exécutent ; la figure des fonctions de transfert se construit."""
import importlib.util
from pathlib import Path

import pytest

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, EXAMPLES / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_quickstart_runs(capsys):
    _load("quickstart").main(plot=False)
    assert "capteur simulé" in capsys.readouterr().out


@pytest.mark.slow
def test_benchmark_tables_run(capsys):
    bench = _load("benchmark")
    bench.table_metrics(3)
    bench.table_zone_plate(factor=3, n=384)
    assert "Mire de Fresnel" in capsys.readouterr().out


def test_plot_transfer_functions():
    from eoqual.resampling.plotting import plot_transfer_functions

    fig = plot_transfer_functions(["box", "sinc_kaiser", "fourier"], 3)
    assert len(fig.data) == 3


def test_documentation_figure(tmp_path):
    pytest.importorskip("matplotlib")
    chemin = _load("figures").figure_grille(chemin=tmp_path / "grille.png")
    assert chemin.stat().st_size > 0


def test_orthophoto_pipeline_on_synthetic_rgb():
    """Chaîne d'évaluation de l'exemple orthophoto sur une petite image RGB synthétique."""
    pytest.importorskip("skimage")
    import numpy as np
    from scipy.ndimage import gaussian_filter

    rng = np.random.default_rng(0)
    rgb = np.stack([gaussian_filter(rng.random((256, 256)), 1.0) for _ in range(3)], -1)
    rgb = np.clip(255 * (rgb - rgb.min()) / np.ptp(rgb), 0, 255)
    sorties, lignes = _load("orthophoto").evaluer(rgb)
    par_nom = {l["libelle"]: l for l in lignes}
    assert par_nom["box"]["ecart_box"] == 0.0
    assert par_nom["fourier"]["repli"] < 1e-12
    assert par_nom["sinc_kaiser"]["repli"] < par_nom["box"]["repli"] < par_nom["nearest"]["repli"]
    # Grille "origin" : décalage de (k - 1) / 2 pixel source = 0,375 pixel cible.
    assert np.allclose(par_nom['sinc_kaiser, grid="origin"']["decalage"], 0.375, atol=0.03)
    assert all(np.allclose(l["decalage"], 0, atol=0.01) for l in lignes
               if l["methode"] != "nearest" and "origin" not in l["libelle"])
