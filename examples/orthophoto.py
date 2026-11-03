#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cas réel : orthophotographie à 5 cm réduite à 20 cm (facteur 4).

    python examples/orthophoto.py capitole_10240.tif              # tableau + figures
    python examples/orthophoto.py capitole_10240.tif --save sorties/   # + images 20 cm (GeoTIFF)

Données : extrait de 10240 × 10240 pixels (512 m × 512 m, place du
Capitole) de l'**Orthophotoplan 2024 de Toulouse Métropole**, publié sous
Licence Ouverte v2.0 (Etalab) :
https://data.toulouse-metropole.fr/explore/dataset/orthophotoplan-2024/
Emprise de l'extrait (EPSG:3943, RGF93 / CC43) :
``gdal_translate -projwin 1574165 2268568.5 1574677 2268056.5 <orthophotoplan> capitole_10240.tif``.

Indicateurs, calculés sur la luminance pour chaque méthode :

- **temps** de réduction de l'image RGB complète ;
- **contraste des détails fins** et **repliement / détails fins** : à
  partir du spectre réel de la scène à 5 cm et de la réponse de chaque
  méthode. Les détails fins sont les fréquences de 0,25 à 0,5 cycle par
  pixel de 20 cm (périodes de 40 à 80 cm) : le contraste y est la racine
  du rapport des énergies après/avant, le repliement l'énergie repliée
  rapportée à l'énergie de ces détails. C'est le comportement attendu sur
  *cette* scène, pas sur un bruit blanc ;
- **décalage** par rapport à l'agrégation exacte (`box`), mesuré par
  corrélation de phase sous-pixel (pixels de 20 cm) ;
- **écart à la moyenne par blocs** (écart-type de la différence avec
  `box`, en niveaux numériques) : ce qui sépare la méthode d'une simple
  agrégation de pixels ;
- **pixels hors [0, 255]** : dépassements dus aux oscillations, qui
  seront écrêtés à l'enregistrement en 8 bits.

Figures (``resources/orthophoto_*.png``) : aperçu, spectre de la scène,
extraits « toits » et « parking » à 20 cm pour huit méthodes.

Dépendances de l'exemple (non requises par le paquet) : tifffile,
matplotlib, scikit-image ; rasterio pour ``--save``. Mémoire : ≈ 8 Go
pour les méthodes spectrales (traitement global, ≈ 20 s par bande).
"""

import argparse
import time
from pathlib import Path

import numpy as np

from eoqual.core.fourier import perdecomp
from eoqual.resampling import METHODS_CONFIGS, resample, transfer_function
from eoqual.resampling.methods import build_response
from eoqual.resampling.mtf import gaussian_mtf, pixel_mtf, product

RESOURCES = Path(__file__).resolve().parents[1] / "resources"
K = 4                                   # 5 cm → 20 cm
ATTRIBUTION = "Orthophotoplan 2024 © Toulouse Métropole, Licence Ouverte v2.0 (Etalab)"

#: Capteur 20 cm simulé pour mtf_target : optique (0,5 au Nyquist) × pixel.
MTF_CAPTEUR = product(gaussian_mtf(0.5), pixel_mtf())

#: (libellé, méthode, paramètres) ; la dernière ligne illustre l'erreur de
#: convention de grille (origine sur le premier pixel au lieu du bloc).
CAS = [(m, m, {"target": MTF_CAPTEUR} if m == "mtf_target" else {}) for m in METHODS_CONFIGS]
CAS.append(('sinc_kaiser, grid="origin"', "sinc_kaiser", {"grid": "origin"}))

#: Extraits illustrés (ligne, colonne, côté en pixels de 5 cm).
EXTRAITS = {"toits": (2700, 5300, 512), "parking": (600, 8500, 512)}
#: Méthodes montrées dans les figures d'extraits.
MONTREES = ["nearest", "bilinear", "box", "gaussian", "lanczos", "sinc_kaiser", "fourier"]


def luminance(rgb):
    return rgb @ np.array([0.299, 0.587, 0.114])


def lire(chemin):
    import tifffile

    img = tifffile.imread(chemin)
    if img.ndim != 3 or img.shape[2] < 3:
        raise ValueError("image RGB attendue")
    return img[..., :3]


def reponse_2d(methode, params, fy, fx):
    """Réponse 2D (cycles par pixel cible) d'une méthode : séparable ou spectrale."""
    params = {k: v for k, v in params.items() if k != "grid"}
    if METHODS_CONFIGS[methode]["backend"] == "spectral":
        response, _ = build_response(methode, (K, K), **params)
        return response(fy, fx)
    axe = np.unique(np.abs(np.concatenate([fy.ravel(), fx.ravel()])))
    _, g = transfer_function(methode, K, freqs=axe, **params)
    return np.interp(np.abs(fy), axe, g) * np.interp(np.abs(fx), axe, g)


def spectre_scene(y):
    """Densité spectrale de la luminance à 5 cm (composante périodique, sans moyenne)."""
    p, _ = perdecomp(y)
    p = p - p.mean()
    psd = np.abs(np.fft.rfft2(p)) ** 2
    fy = np.fft.fftfreq(y.shape[0])[:, None] * K          # cycles par pixel cible
    fx = np.fft.rfftfreq(y.shape[1])[None, :] * K
    poids = np.where((np.arange(psd.shape[1]) == 0) | (np.arange(psd.shape[1]) == psd.shape[1] - 1), 1.0, 2.0)
    return psd * poids[None, :], fy, fx


def decalage(ref, img, cote=1024):
    """Décalage sous-pixel (lignes, colonnes) de img par rapport à ref, au centre."""
    from skimage.registration import phase_cross_correlation

    cote = min(cote, *ref.shape)
    c0, c1 = (ref.shape[0] - cote) // 2, (ref.shape[1] - cote) // 2
    sl = (slice(c0, c0 + cote), slice(c1, c1 + cote))
    shift, _, _ = phase_cross_correlation(ref[sl], img[sl], upsample_factor=200, normalization=None)
    return -shift


def evaluer(img):
    """Réduit l'image par chaque méthode ; renvoie résultats et indicateurs."""
    rgb = img.astype(np.float64)
    y_in = luminance(rgb)
    psd, fy, fx = spectre_scene(y_in)
    fmax = np.maximum(np.abs(fy), np.abs(fx))
    dehors = fmax > 0.5                              # se replie à 20 cm
    fin = (fmax > 0.25) & ~dehors                    # détails fins représentables
    e_fin = np.sum(psd * fin)

    sorties, lignes = {}, []
    for libelle, methode, params in CAS:
        t0 = time.perf_counter()
        out = resample(rgb, K, method=methode, channel_axis=-1, **params)
        duree = time.perf_counter() - t0
        sorties[libelle] = out
        g2 = reponse_2d(methode, params, fy, fx) ** 2
        lignes.append({"libelle": libelle, "methode": methode, "temps": duree,
                       "contraste": np.sqrt(np.sum(psd * g2 * fin) / e_fin),
                       "repli": np.sum(psd * g2 * dehors) / np.sum(psd * g2 * fin),
                       "hors": 100.0 * np.mean((out < -0.5) | (out > 255.5)),
                       "moyenne": abs(luminance(out).mean() - y_in.mean())})

    y_box = luminance(sorties["box"])
    for ligne in lignes:
        y_out = luminance(sorties[ligne["libelle"]])
        ligne["decalage"] = decalage(y_box, y_out)
        ligne["ecart_box"] = float(np.std(y_out - y_box))
    return sorties, lignes


def tableau_markdown(lignes):
    fr = lambda v, nd: f"{v:.{nd}f}".replace(".", ",")
    out = ["| Méthode | temps (s) | contraste des détails fins | repliement / détails fins | décalage (px 20 cm) | écart à la moyenne par blocs (NN) | pixels hors [0, 255] |",
           "|---|---|---|---|---|---|---|"]
    for l in lignes:
        dy, dx = l["decalage"]
        dec = "0" if max(abs(dy), abs(dx)) < 0.005 else f"{fr(dy, 3)} ; {fr(dx, 3)}"
        out.append(f"| `{l['libelle']}` | {fr(l['temps'], 1)} | {fr(l['contraste'], 3)} | "
                   f"{fr(100 * l['repli'], 2)} % | {dec} | {fr(l['ecart_box'], 2)} | {fr(l['hors'], 3)} % |")
    return "\n".join(out)


# --- Figures --------------------------------------------------------------------

def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def figure_apercu(img, dossier):
    plt = _plt()
    from matplotlib.patches import Rectangle

    pas = 16
    fig, ax = plt.subplots(figsize=(5.6, 6.0), layout="constrained")
    ax.imshow(img[::pas, ::pas])
    for nom, (r, c, n) in EXTRAITS.items():
        ax.add_patch(Rectangle((c / pas, r / pas), n / pas, n / pas, fill=False, ec="yellow", lw=2))
        ax.text(c / pas, r / pas - 6, nom, color="yellow", fontsize=11, weight="bold")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("Place du Capitole, 512 m × 512 m à 5 cm (10240 × 10240 px)", fontsize=10)
    ax.set_xlabel(ATTRIBUTION, fontsize=7)
    chemin = dossier / "orthophoto_apercu.png"
    fig.savefig(chemin, dpi=110)
    plt.close(fig)
    # Vue de situation seulement : palette de 256 couleurs pour alléger le PNG.
    from PIL import Image

    Image.open(chemin).convert("RGB").quantize(256, dither=Image.Dither.NONE).save(chemin, optimize=True)
    return chemin


def figure_spectre(img, dossier):
    """Spectre de la luminance à 5 cm et réponses de quelques méthodes."""
    plt = _plt()
    psd, fy, fx = spectre_scene(luminance(img.astype(np.float64)))
    fmax = np.maximum(np.abs(fy), np.abs(fx))
    bords = np.linspace(0, 2, 81)
    idx = np.digitize(fmax.ravel(), bords) - 1
    somme = np.bincount(idx, psd.ravel(), minlength=bords.size)[:bords.size - 1]
    nombre = np.bincount(idx, minlength=bords.size)[:bords.size - 1]
    centres = 0.5 * (bords[1:] + bords[:-1])
    dsp = somme / np.maximum(nombre, 1)

    fig, ax = plt.subplots(figsize=(8, 4), layout="constrained")
    ax.axvspan(0.5, 2, color="0.92", zorder=0)
    ax.semilogy(centres, dsp / dsp[0], "k", lw=2, label="spectre de la scène à 5 cm")
    ax.set_xlabel("fréquence (cycles par pixel de 20 cm ; max des deux composantes)")
    ax.set_ylabel("densité spectrale (relative)")
    ax.text(0.53, 0.5, "repliement à 20 cm", transform=ax.get_xaxis_transform(), fontsize=9)
    ax2 = ax.twinx()
    f = np.linspace(0, 2, 401)
    for m in ("box", "gaussian", "sinc_kaiser", "fourier"):
        ax2.plot(f, np.abs(transfer_function(m, K, freqs=f)[1]), lw=1.2, ls="--", label=m)
    ax2.set_ylabel("gain de la méthode")
    ax2.set_ylim(0, 1.05)
    ax.set_xlim(0, 2)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper right")
    ax.axvline(0.5, color="k", lw=0.8)
    ax.set_title(ATTRIBUTION, fontsize=7, loc="right")
    chemin = dossier / "orthophoto_spectre.png"
    fig.savefig(chemin, dpi=130)
    plt.close(fig)
    return chemin


def figure_extrait(img, sorties, nom, dossier):
    plt = _plt()
    r, c, n = EXTRAITS[nom]
    fig, axs = plt.subplots(2, 4, figsize=(12, 6.6), layout="constrained")
    axs = axs.ravel()
    axs[0].imshow(img[r:r + n, c:c + n], interpolation="nearest")
    axs[0].set_title("origine, 5 cm", fontsize=10)
    for ax, m in zip(axs[1:], MONTREES):
        bloc = sorties[m][r // K:(r + n) // K, c // K:(c + n) // K]
        ax.imshow(np.clip(np.rint(bloc), 0, 255).astype(np.uint8), interpolation="nearest")
        ax.set_title(f"{m}, 20 cm", fontsize=10)
    for ax in axs:
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle(f"Extrait « {nom} » ({n * 0.05:g} m × {n * 0.05:g} m) — {ATTRIBUTION}", fontsize=9)
    chemin = dossier / f"orthophoto_{nom}.png"
    fig.savefig(chemin, dpi=110)
    plt.close(fig)
    return chemin


def enregistrer(source, sorties, dossier):
    """Écrit chaque résultat en GeoTIFF 8 bits géoréférencé à 20 cm (nécessite rasterio)."""
    import rasterio
    from rasterio.transform import Affine

    dossier.mkdir(parents=True, exist_ok=True)
    with rasterio.open(source) as src:
        profil = src.profile
        t = src.transform
    for libelle, out in sorties.items():
        profil.update(width=out.shape[1], height=out.shape[0], count=3, dtype="uint8",
                      transform=t * Affine.scale(K), compress="deflate", predictor=2,
                      tiled=True, blockxsize=256, blockysize=256)
        if "origin" in libelle:      # le premier pixel est centré sur le premier pixel source
            profil.update(transform=t * Affine.translation(-(K - 1) / 2, -(K - 1) / 2) * Affine.scale(K))
        nom = libelle.split(",")[0] + ("_origin" if "origin" in libelle else "")
        with rasterio.open(dossier / f"capitole_20cm_{nom}.tif", "w", **profil) as dst:
            dst.write(np.clip(np.rint(out), 0, 255).astype(np.uint8).transpose(2, 0, 1))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("image", type=Path, help="orthophotographie RGB à 5 cm (GeoTIFF)")
    parser.add_argument("--figures", type=Path, default=RESOURCES, help="dossier des figures PNG")
    parser.add_argument("--save", type=Path, default=None, help="dossier où écrire les images à 20 cm")
    args = parser.parse_args()

    img = lire(args.image)
    sorties, lignes = evaluer(img)
    print(tableau_markdown(lignes))
    args.figures.mkdir(parents=True, exist_ok=True)
    print(figure_apercu(img, args.figures))
    print(figure_spectre(img, args.figures))
    for nom in EXTRAITS:
        print(figure_extrait(img, sorties, nom, args.figures))
    if args.save is not None:
        enregistrer(args.image, sorties, args.save)


if __name__ == "__main__":
    main()
