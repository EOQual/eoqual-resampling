# Benchmark — eoqual-resampling 0.1.0

Comparaison chiffrée des méthodes de réduction. Tableaux générés par
`python examples/benchmark.py` (données synthétiques, déterministes) ;
définitions des indicateurs : [REFERENCE_TECHNIQUE.md](REFERENCE_TECHNIQUE.md) §7.

Pour un cas réel (orthophotographie de 5 cm à 20 cm, extraits commentés),
voir [REFERENCE_TECHNIQUE.md](REFERENCE_TECHNIQUE.md) §8.

## 1. Protocole

1. **Fonction de transfert** (`kernel_metrics`), facteurs 2, 3 et 4 :
   gain au Nyquist cible, fréquences de gain 0,9 et 0,1 (cycles par pixel
   cible), énergie repliée rapportée à l'énergie conservée (spectre blanc), rapport des
   variances d'un bruit blanc (2D), dépassement d'un échelon unité (pire
   phase sous-pixel).
2. **Mire de Fresnel** 1536 × 1536 (fréquence locale croissant du centre
   vers les bords), réduite d'un facteur 3 : contraste conservé là où la
   fréquence cible est inférieure à 0,35 (idéal 1) et contraste résiduel
   là où elle dépasse 0,6 (idéal 0 : ce contraste est du repliement).
   La fréquence considérée est la plus grande des deux composantes de la
   fréquence locale : la grille cible est carrée, une fréquence
   diagonale de composantes inférieures à 0,5 est représentable.
3. **Temps de calcul**, image 2048 × 2048, Apple M1 Max, un cœur
   (numpy/scipy), meilleur de 3 essais.

`mtf_target` est évalué avec une MTF visée optique gaussienne (0,5 au
Nyquist) × pixel, soit 0,32 au Nyquist cible : ce n'est pas un filtre
anti-repliement, il replie comme le capteur simulé.

## 2. Fonction de transfert

### Facteur 2

| Méthode | MTF au Nyquist | f90 | f10 | transition | énergie repliée / conservée | variance du bruit | dépassement |
|---|---|---|---|---|---|---|---|
| `nearest` | 1.000 | — | — | — | 0.9975 | 1.003 | 0.0 % |
| `bilinear` | 0.707 | 0.288 | 0.936 | 0.649 | 0.2208 | 0.251 | 0.0 % |
| `bicubic` | 0.884 | 0.480 | 0.958 | 0.478 | 0.3436 | 0.411 | 6.2 % |
| `box` | 0.707 | 0.288 | 0.936 | 0.649 | 0.2208 | 0.251 | 0.0 % |
| `gaussian` | 0.291 | 0.146 | 0.684 | 0.537 | 0.0267 | 0.080 | 0.0 % |
| `bicubic_aa` | 0.486 | 0.285 | 0.715 | 0.430 | 0.0483 | 0.164 | 4.7 % |
| `mitchell` | 0.378 | 0.204 | 0.682 | 0.479 | 0.0324 | 0.115 | 1.9 % |
| `lanczos` | 0.502 | 0.366 | 0.636 | 0.270 | 0.0292 | 0.200 | 8.2 % |
| `sinc_gaussian` | 0.500 | 0.421 | 0.581 | 0.160 | 0.0157 | 0.218 | 9.3 % |
| `sinc_kaiser` | 0.500 | 0.448 | 0.554 | 0.106 | 0.0100 | 0.228 | 9.6 % |
| `fourier` | 0.000 | 0.421 | 0.480 | 0.059 | 0.0000 | 0.192 | 9.3 % |
| `mtf_target` | 0.318 | 0.155 | 0.694 | 0.539 | 0.0304 | 0.087 | 0.1 % |
| `sinc_kaiser` (transition=0.1) | 0.500 | 0.474 | 0.527 | 0.054 | 0.0047 | 0.240 | 9.9 % |
| `sinc_kaiser` (radius=6.0) | 0.500 | 0.421 | 0.580 | 0.159 | 0.0158 | 0.218 | 9.3 % |
| `fourier` (passband=0.45, stopband=0.55) | 0.500 | 0.471 | 0.530 | 0.059 | 0.0055 | 0.238 | 9.9 % |

### Facteur 3

| Méthode | MTF au Nyquist | f90 | f10 | transition | énergie repliée / conservée | variance du bruit | dépassement |
|---|---|---|---|---|---|---|---|
| `nearest` | 1.000 | — | — | — | 2.0000 | 1.003 | 0.0 % |
| `bilinear` | 1.000 | — | — | — | 2.0000 | 1.003 | 0.0 % |
| `bicubic` | 1.000 | — | — | — | 2.0000 | 1.003 | 0.0 % |
| `box` | 0.667 | 0.266 | 0.922 | 0.656 | 0.2615 | 0.112 | 0.0 % |
| `gaussian` | 0.291 | 0.148 | 0.684 | 0.536 | 0.0270 | 0.036 | 0.0 % |
| `bicubic_aa` | 0.494 | 0.291 | 0.713 | 0.422 | 0.0493 | 0.074 | 3.7 % |
| `mitchell` | 0.384 | 0.204 | 0.682 | 0.478 | 0.0333 | 0.052 | 1.9 % |
| `lanczos` | 0.502 | 0.367 | 0.636 | 0.268 | 0.0297 | 0.089 | 6.5 % |
| `sinc_gaussian` | 0.501 | 0.422 | 0.581 | 0.159 | 0.0162 | 0.097 | 7.4 % |
| `sinc_kaiser` | 0.500 | 0.448 | 0.553 | 0.105 | 0.0105 | 0.102 | 7.7 % |
| `fourier` | 0.000 | 0.422 | 0.480 | 0.058 | 0.0000 | 0.085 | 9.0 % |
| `mtf_target` | 0.318 | 0.156 | 0.694 | 0.538 | 0.0308 | 0.039 | 0.1 % |
| `sinc_kaiser` (transition=0.1) | 0.500 | 0.474 | 0.527 | 0.052 | 0.0052 | 0.107 | 8.0 % |
| `sinc_kaiser` (radius=6.0) | 0.500 | 0.422 | 0.581 | 0.159 | 0.0163 | 0.097 | 7.4 % |
| `fourier` (passband=0.45, stopband=0.55) | 0.500 | 0.471 | 0.531 | 0.060 | 0.0060 | 0.106 | 7.9 % |

### Facteur 4

| Méthode | MTF au Nyquist | f90 | f10 | transition | énergie repliée / conservée | variance du bruit | dépassement |
|---|---|---|---|---|---|---|---|
| `nearest` | 1.000 | — | — | — | 2.9851 | 1.003 | 0.0 % |
| `bilinear` | 0.924 | 0.575 | 1.873 | 1.297 | 1.0973 | 0.251 | 0.0 % |
| `bicubic` | 0.992 | 0.960 | 1.915 | 0.955 | 1.5611 | 0.411 | 6.2 % |
| `box` | 0.653 | 0.260 | 0.917 | 0.657 | 0.2724 | 0.063 | 0.0 % |
| `gaussian` | 0.291 | 0.147 | 0.685 | 0.538 | 0.0265 | 0.020 | 0.0 % |
| `bicubic_aa` | 0.493 | 0.290 | 0.713 | 0.423 | 0.0483 | 0.042 | 4.3 % |
| `mitchell` | 0.383 | 0.205 | 0.682 | 0.477 | 0.0325 | 0.029 | 1.7 % |
| `lanczos` | 0.502 | 0.367 | 0.638 | 0.270 | 0.0287 | 0.050 | 7.5 % |
| `sinc_gaussian` | 0.500 | 0.422 | 0.583 | 0.160 | 0.0154 | 0.055 | 8.5 % |
| `sinc_kaiser` | 0.500 | 0.448 | 0.555 | 0.108 | 0.0097 | 0.057 | 8.9 % |
| `fourier` | 0.000 | 0.422 | 0.480 | 0.057 | 0.0000 | 0.048 | 8.6 % |
| `mtf_target` | 0.318 | 0.155 | 0.695 | 0.540 | 0.0302 | 0.022 | 0.1 % |
| `sinc_kaiser` (transition=0.1) | 0.500 | 0.475 | 0.527 | 0.052 | 0.0044 | 0.060 | 9.1 % |
| `sinc_kaiser` (radius=6.0) | 0.500 | 0.422 | 0.580 | 0.157 | 0.0154 | 0.055 | 8.5 % |
| `fourier` (passband=0.45, stopband=0.55) | 0.500 | 0.473 | 0.530 | 0.057 | 0.0052 | 0.060 | 9.1 % |

## 3. Mire de Fresnel

### Mire de Fresnel 1536 × 1536, facteur 3

Contraste RMS rapporté à celui de la mire (0,354) : conservé en bande passante (idéal 1), résiduel au-delà du Nyquist cible (idéal 0).

f : plus grande composante de la fréquence locale, en cycles par pixel cible.

| Méthode | contraste conservé (f < 0,35) | repliement résiduel (f > 0,6) |
|---|---|---|
| `nearest` | 1.002 | 1.002 |
| `bilinear` | 1.002 | 1.002 |
| `bicubic` | 1.002 | 1.002 |
| `box` | 0.890 | 0.258 |
| `gaussian` | 0.707 | 0.049 |
| `bicubic_aa` | 0.920 | 0.081 |
| `mitchell` | 0.816 | 0.055 |
| `lanczos` | 1.000 | 0.039 |
| `sinc_gaussian` | 1.004 | 0.009 |
| `sinc_kaiser` | 1.001 | 0.003 |
| `fourier` | 1.001 | 0.002 |
| `mtf_target` | 0.727 | 0.053 |
| `sinc_kaiser` (transition=0.1) | 1.002 | 0.003 |
| `sinc_kaiser` (radius=6.0) | 1.003 | 0.008 |
| `fourier` (passband=0.45, stopband=0.55) | 1.002 | 0.003 |

## 4. Temps de calcul

### Temps de calcul, image 2048 × 2048 (meilleur de 3 essais)

| Méthode | facteur 3 (s) | facteur 2,5 (s) |
|---|---|---|
| `nearest` | 0.01 | 0.01 |
| `bilinear` | 0.01 | 0.01 |
| `bicubic` | 0.01 | 0.01 |
| `box` | 0.01 | 0.01 |
| `gaussian` | 0.01 | 0.01 |
| `bicubic_aa` | 0.01 | 0.01 |
| `mitchell` | 0.01 | 0.01 |
| `lanczos` | 0.01 | 0.01 |
| `sinc_gaussian` | 0.02 | 0.02 |
| `sinc_kaiser` | 0.03 | 0.03 |
| `fourier` | 0.43 | 0.51 |
| `mtf_target` | 0.48 | 0.52 |


## 5. Lecture

- **Méthodes de référence** (`nearest`, `bilinear`, `bicubic`) : sans
  préfiltrage, tout se replie (contraste résiduel 1 sur la mire). Pour
  un facteur pair et la grille `"center"`, `bilinear` moyenne deux pixels
  et se comporte comme `box`.
- **`box`** : aucune oscillation, mais une énergie repliée égale à 26 % de
  l'énergie conservée, et le
  contraste baisse dès les basses fréquences. À réserver au binning.
- **`gaussian`, `mitchell`, `bicubic_aa`** : peu de repliement, peu ou pas
  d'oscillation, mais flou net (contraste conservé 0,71 à 0,92).
- **Sinus cardinaux fenêtrés** (`lanczos`, `sinc_gaussian`,
  `sinc_kaiser`) : contraste intégralement conservé en bande passante,
  repliement résiduel 0,039 / 0,009 / 0,003 sur la mire, au prix de 6,5 à
  8 % de dépassement près des bords nets. `sinc_kaiser` (méthode par
  défaut) divise le repliement mesuré par trois par rapport à
  `sinc_gaussian` ; à support égal (`radius=6`), les deux sont
  équivalents.
- **`fourier`** : repliement nul par construction (énergie au-delà du
  Nyquist nulle), contraste conservé, transition la plus étroite ; plus
  d'oscillations (9 %) et environ 20 fois plus lent que les méthodes
  séparables. Le résidu de 0,002 sur la mire vient de la composante
  lisse et des bords.
- **Bruit** : plus la méthode est nette, plus le bruit réduit est fort
  (et corrélé) ; `box` atteint la borne `1 / k²`.
