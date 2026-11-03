<p align="center">
  <img src="resources/eoqual_rounded.png" alt="EOQual" width="280">
</p>

# eoqual-resampling

Réduction de résolution d'images pensée pour l'observation de la Terre :
à partir d'une image haute résolution, produire une image plus basse
résolution « idéale », sans repliement, sans décalage géométrique et sans
perte des niveaux, pour un facteur quelconque (entier ou non, différent
par axe). Voir **[REFERENCE_TECHNIQUE.md](REFERENCE_TECHNIQUE.md)** pour
le principe, l'implémentation et les limites de chaque méthode,
**[BENCHMARK.md](BENCHMARK.md)** pour la comparaison chiffrée,
**[THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md)** pour l'audit des
licences et **[BIBLIOGRAPHY.md](BIBLIOGRAPHY.md)** pour les publications
de référence.

## Sommaire

1. [Installation](#installation)
2. [Démarrage rapide](#démarrage-rapide)
3. [Catalogue des méthodes](#catalogue-des-méthodes)
4. [Quelle méthode choisir ?](#quelle-méthode-choisir-)
5. [Conventions](#conventions)
6. [Modules EOQual utilisés](#modules-eoqual-utilisés)
7. [Licence](#licence)
8. [Contribuer](#contribuer)

---

## Installation

```bash
pip install git+https://github.com/EOQual/eoqual-resampling.git
```

Dépendances : numpy, scipy, plotly, rich, et les modules EOQual
eoqual-core et eoqual-viz, tous sous licence permissive.

## Démarrage rapide

```python
from eoqual.resampling import resample, list_methods, kernel_metrics

list_methods()                                   # méthodes, famille, paramètres, référence

small = resample(image, 3)                       # sinc fenêtré par Kaiser (défaut), facteur 3
small = resample(image, 2.5, method="sinc_gaussian")        # facteur non entier
small = resample(image, shape=(400, 600), method="fourier") # taille imposée
small = resample(rgb, 4, channel_axis=-1)                   # image multi-bandes

kernel_metrics("sinc_kaiser", 3)   # MTF au Nyquist, repliement, oscillations, bruit
```

Simuler l'image d'un capteur moins résolu (MTF visée imposée) :

```python
from eoqual.resampling.mtf import gaussian_mtf, pixel_mtf, product

target = product(gaussian_mtf(0.5), pixel_mtf())   # optique × intégration sur le pixel
sim = resample(image, 4, method="mtf_target", target=target)
```

Exemples complets : `examples/quickstart.py` (avec `--plot` pour la
figure des fonctions de transfert), `examples/benchmark.py` (mires
synthétiques) et `examples/orthophoto.py` — cas réel : une
orthophotographie de Toulouse réduite de 5 cm à 20 cm par les douze
méthodes, avec indicateurs et extraits commentés dans
[REFERENCE_TECHNIQUE.md](REFERENCE_TECHNIQUE.md) §8.

## Catalogue des méthodes

| Méthode | Famille | Principe | Paramètres (défaut) |
|---|---|---|---|
| `nearest` | référence | plus proche voisin, sans préfiltrage | — |
| `bilinear` | référence | interpolation linéaire, sans préfiltrage | — |
| `bicubic` | référence | cubique de Keys, sans préfiltrage | `a=-0.5` |
| `box` | anti-repliement | agrégation exacte par surface (binning) | — |
| `gaussian` | anti-repliement | gaussienne puis échantillonnage | `sigma=0.5` |
| `bicubic_aa` | anti-repliement | cubique de Keys dilatée du facteur | `a=-0.5` |
| `mitchell` | anti-repliement | Mitchell-Netravali dilaté du facteur | `b=c=1/3` |
| `lanczos` | anti-repliement | Lanczos (sinc fenêtré par un sinc) | `a=3` |
| `sinc_gaussian` | anti-repliement | sinc apodisé par une gaussienne | `fwhm=6`, `radius=6`, `cutoff=0.5` |
| **`sinc_kaiser`** | anti-repliement | sinc fenêtré par Kaiser, réglé par l'atténuation et la transition visées (**défaut**) | `attenuation=60`, `transition=0.2`, `cutoff=0.5` |
| `fourier` | spectral | gabarit en cosinus surélevé, aucun repliement | `passband=0.4`, `stopband=0.5` |
| `mtf_target` | spectral | MTF visée imposée (capteur simulé), MTF source compensée | `target` (requis), `source`, `floor=0.05` |

Les paramètres sont en unités **cibles** (pixels cibles, cycles par pixel
cible ; Nyquist cible = 0,5) : un réglage garde le même sens quel que
soit le facteur. Les méthodes « référence » ne filtrent pas avant
d'échantillonner : elles replient tout le spectre et ne servent qu'à la
comparaison. Toutes les méthodes sont des implémentations originales
sous MIT, écrites à partir des publications citées.

Résultats au facteur 3 (détail dans [BENCHMARK.md](BENCHMARK.md)) :

| Méthode | MTF au Nyquist | énergie repliée / conservée | repliement mesuré sur mire | dépassement d'un bord |
|---|---|---|---|---|
| `box` | 0,67 | 26 % | 0,26 | 0 % |
| `lanczos` | 0,50 | 3,0 % | 0,039 | 6,5 % |
| `sinc_gaussian` | 0,50 | 1,6 % | 0,009 | 7,4 % |
| `sinc_kaiser` | 0,50 | 1,1 % | 0,003 | 7,7 % |
| `fourier` | 0 | 0 % | 0,002 | 9,0 % |

## Quelle méthode choisir ?

Aucune méthode n'est la meilleure sur tout : réduire, c'est arbitrer
entre **repliement** (fréquences au-delà du nouveau Nyquist qui
reviennent en fausses basses fréquences), **flou** (contraste perdu en
deçà) et **oscillations** près des bords nets (phénomène de Gibbs).

- **Image réduite « idéale »** (référence, vérité terrain d'une
  super-résolution, comparaison de produits) : `sinc_kaiser` (défaut) —
  même netteté que `sinc_gaussian`, trois fois moins de repliement
  mesuré, et le compromis se règle directement (`attenuation`,
  `transition`). `sinc_gaussian`, à support 1,5 fois plus court, reste un
  très bon choix.
- **Aucun repliement toléré** : `fourier` (gabarit qui s'annule avant le
  Nyquist cible ; exact pour une image à bande limitée), au prix d'un
  contraste nul au Nyquist et d'un calcul global (≈ 20 fois plus lent).
- **Pas d'oscillation tolérée** (masques, bords nets, images
  thématiques) : `gaussian` ou `box`, plus flous ou plus repliés.
- **Simuler un capteur** : `mtf_target`, avec la MTF de l'instrument
  visé ; elle replie comme le ferait le vrai capteur.
- **Agréger des pixels** (binning) : `box`, seule méthode à conserver
  exactement la moyenne — donc le flux — de chaque bloc de pixels source.

## Conventions

- **Appel unifié** : `resample(image, factor | shape=, method=, **params)`
  renvoie une image `float64`.
- **Image** : 2D ; une image multi-bandes se traite avec
  `channel_axis=`. Les valeurs non finies (NaN) sont refusées.
- **Facteur** : `>= 1` (réduction), scalaire ou par axe ; taille de sortie
  `floor(n / k)`. Avec `shape=`, le facteur vaut exactement `n / shape`.
  L'agrandissement (facteur < 1) n'est pas implémenté dans cette version
  (erreur explicite) ; l'architecture le permet (REFERENCE_TECHNIQUE.md
  §9).
- **Grille** (`grid=`) : `"center"` (défaut) — chaque pixel cible couvre
  exactement un bloc de `k` pixels source, l'image réduite a la même
  emprise que l'image source (convention du binning) ; `"origin"` — le
  premier pixel cible est centré sur le premier pixel source. Aucune
  méthode ne décale l'image (vérifié par test), sauf `nearest` à facteur
  pair sur la grille `"center"` (demi-pixel source, inhérent à la
  méthode).
- **Bords** (`boundary=`) : `"symmetric"` (miroir demi-pixel, défaut),
  `"reflect"`, `"edge"`, `"wrap"`. Les méthodes spectrales traitent les
  bords par la décomposition périodique + lisse (sauf `"wrap"`).
- **Niveaux** : une image uniforme est conservée exactement par toutes
  les méthodes (vérifié par test).

## Modules EOQual utilisés

| Module | Usage |
|---|---|
| [eoqual-core](https://github.com/EOQual/eoqual-core) | décomposition périodique + lisse (méthodes spectrales) |
| [eoqual-viz](https://github.com/EOQual/eoqual-viz) | figure des fonctions de transfert |

---

## Licence

MIT — voir `license.txt`. Licences des dépendances :
[THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md).

## Contribuer

Voir [CONTRIBUTING.md](CONTRIBUTING.md).
