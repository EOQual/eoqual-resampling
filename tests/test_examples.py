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
