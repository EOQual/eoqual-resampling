# Méthodes de réduction

## 3. Méthodes séparables

Toutes les méthodes sauf `fourier` et `mtf_target` appliquent un noyau
1D `h` séparable : l'image réduite vaut `W_lignes · image · W_colonnesᵀ`,
où la ligne `i` de la matrice creuse `W` (`n_out × n_in`) contient les
poids `h(x_i - j)`.

- **Évaluation directe aux positions cibles.** Le noyau est évalué
  exactement en `x_i - j` pour chaque pixel cible : pas de filtrage sur
  la grille source suivi d'une interpolation, donc pas de flou
  supplémentaire variable d'un pixel à l'autre quand `k` n'est pas
  entier.
- **Dilatation.** Un noyau anti-repliement est défini en pixels cibles,
  `phi(u)`, et appliqué en pixels source sous la forme
  `h(x) = phi(x / k) / k` : sa fréquence de coupure suit le nouveau
  Nyquist. Les méthodes « référence » (`nearest`, `bilinear`, `bicubic`)
  ne sont **pas** dilatées : simple interpolation aux positions cibles,
  sans préfiltrage — elles replient tout le spectre.
- **Normalisation.** Chaque ligne de `W` est ramenée à une somme de 1 :
  une image uniforme est conservée exactement, quelle que soit la phase
  sous-pixel de `x_i`.
- **Troncature stricte.** Les poids sont nuls au-delà du rayon déclaré du
  noyau : un poids résiduel d'un seul côté décalerait l'image (voir
  [§11](conception.md#11-choix-de-conception-et-variantes-écartées)).
- **Bords** (`boundary=`) : les indices hors image sont repliés —
  `"symmetric"` (miroir demi-pixel, cohérent avec la convention centre de
  pixel ; défaut), `"reflect"` (miroir sur le pixel de bord), `"edge"`
  (pixel de bord répété), `"wrap"` (périodique).
- **Coût** : `O(n · support)` par axe, matrices creuses (scipy). Une
  image 2048 × 2048 se réduit en 10 à 30 ms.

`box` utilise la même mécanique avec un poids exact : la longueur du
recouvrement entre le pixel source `[j - 1/2, j + 1/2]` et l'empreinte
`[x_i - k/2, x_i + k/2]` du pixel cible (moyenne par surface, « area »).
Pour `k` entier et la grille `"center"`, c'est exactement la moyenne par
blocs (vérifié par test).

## 4. Les noyaux

`u` est en pixels cibles pour les noyaux dilatés.

### 4.1 Méthodes de référence (sans préfiltrage)

| Méthode | `phi(u)` | Remarque |
|---|---|---|
| `nearest` | 1 sur `[-1/2, 1/2[` | pour `k` pair et `grid="center"`, `x_i` tombe entre deux pixels : le second est retenu (décalage d'un demi-pixel source, inhérent) |
| `bilinear` | `max(0, 1 - abs(u))` | |
| `bicubic` | cubique de Keys, `a = -0,5` | |

Pour `k` impair et `grid="center"`, les positions cibles tombent sur des
pixels source : ces trois méthodes reviennent à prendre un pixel sur `k`
(gain 1 à toutes les fréquences, repliement total). Elles servent de
borne basse au benchmark.

### 4.2 `box` — agrégation

MTF : `sinc(f)` (pixel cible carré), 0,64 au Nyquist pour `k` entier
grand (0,67 à `k = 3` sur la grille discrète). Aucune oscillation, mais
repliement important (lobes secondaires du sinc : énergie repliée égale
à 26 % de l'énergie conservée à `k = 3`, voir [§7](diagnostic.md#7-diagnostic)). C'est le modèle d'intégration d'un détecteur aux pixels `k`
fois plus grands, sans optique.

### 4.3 `gaussian`

`phi(u) = exp(-u² / 2σ²)`, `σ = 0,5` pixel cible par défaut, tronquée à
`4σ`. MTF `exp(-2π² σ² f²)`, soit 0,29 au Nyquist. Aucune oscillation,
repliement faible, mais flou net dès les basses fréquences (gain 0,9 dès
`f = 0,15`). C'est le filtre des pyramides d'images (Burt & Adelson,
1983).

### 4.4 `bicubic_aa`, `mitchell`

Famille de cubiques par morceaux de Mitchell & Netravali (1988),
paramètres `(B, C)` : `bicubic_aa` est la cubique de Keys (`B = 0`,
`C = -a = 0,5`), `mitchell` le compromis `B = C = 1/3` recommandé par ses
auteurs. Dilatées du facteur, comme le fait Pillow : bons compromis
pratiques, à support court (2 pixels cibles), mais plus floues et plus
repliées que les sinus cardinaux fenêtrés.

### 4.5 `lanczos`

`phi(u) = sinc(u) · sinc(u / a)` sur `|u| < a`, `a = 3` par défaut
(Duchon, 1979). Sinus cardinal fenêtré par le lobe central d'un sinus
cardinal plus large. Ondulation de +1,2 % en bande passante,
transition 0,37 → 0,64.

### 4.6 `sinc_gaussian` — sinc apodisé par une gaussienne

`phi(u) = 2 f_c sinc(2 f_c u) · exp(-u² / 2σ²)` sur `|u| < radius`, avec
`f_c = 0,5` (coupure au Nyquist cible : gain 1/2 au Nyquist),
`σ = fwhm / 2,355`, `fwhm = 6` et `radius = 6` pixels cibles. C'est le
sinus cardinal apodisé des notes de travail de l'auteur (voir [§11](conception.md#11-choix-de-conception-et-variantes-écartées)). Transition 0,42 → 0,58, énergie repliée 1,6 % de l'énergie conservée.

### 4.7 `sinc_kaiser` — sinc fenêtré réglé par l'objectif (méthode par défaut)

`phi(u) = 2 f_c sinc(2 f_c u) · I0(β √(1 - (u/R)²)) / I0(β)` sur
`|u| < R` (fenêtre de Kaiser, 1974). Au lieu de régler la fenêtre, on
fixe l'objectif :

- `attenuation` (dB, défaut 60) → `β` par la formule de Kaiser
  (`β = 0,1102 (A - 8,7)` pour `A > 50`) ;
- `transition` (largeur de la bande de transition, cycles par pixel
  cible, défaut 0,2) → longueur `L = (A - 7,95) / (2,285 · 2π Δf)`,
  `R = L / 2` (9,1 pixels cibles par défaut) ;
- `radius` impose le support à la place de `transition`.

**À support égal, Kaiser et gaussienne se valent** (`radius=6` : mêmes
indicateurs que `sinc_gaussian`, voir [benchmark](benchmark.md)) : le gain vient du
support plus long, choisi pour atteindre l'objectif. Réglage par défaut :
transition 0,45 → 0,55, énergie repliée 1,1 %, repliement trois fois moindre
que `sinc_gaussian` sur la mire de Fresnel. `transition=0.1` divise
encore le repliement par deux (support 18 pixels cibles).

## 5. Méthodes spectrales

`fourier` et `mtf_target` imposent une réponse en fréquence `G(fy, fx)`
quelconque, puis échantillonnent exactement.

1. **Décomposition périodique + lisse** (Moisan, 2011 ;
   `eoqual.core.fourier.perdecomp`) : `u = p + s`, où `p` est périodique
   (pas de saut de bord, donc ni « croix » spectrale ni oscillations de
   bord) et `s` est lisse, de moyenne nulle.
2. **Filtrage de `p`** : le spectre de `p` est multiplié par `G`.
3. **Échantillonnage exact** : la sortie est l'interpolant
   trigonométrique de `p` filtré, évalué aux positions cibles `x_i`.
   - Si le pas cible vaut exactement `N / M` (facteur issu de `shape=`,
     ou `M·k` entier — l'image est alors recadrée à l'emprise `M·k` des
     pixels cibles avant la décomposition), l'évaluation en
     `x_i = i N / M + δ` **replie les fréquences modulo `M`** : le
     spectre (multiplié par la phase `exp(2iπ m δ / N)`) est replié puis
     inversé par une FFT de taille `M`. C'est exact et en `O(N² log N)`.
   - Sinon, l'interpolant est évalué par une transformée de Fourier non
     uniforme séparable (deux produits matriciels par blocs) : exact,
     plus lent.
   - Pour une taille paire, le coefficient de Nyquist est partagé par
     moitié entre `+N/2` et `-N/2` (interpolant réel).
4. **Composante lisse** : `s` (très basses fréquences) est réduite par un
   Lanczos-3 séparable, de gain ≈ 1 à ces fréquences.

Une image périodique (`boundary="wrap"` ou `periodic=True`) n'est ni
décomposée ni recadrée (le recadrage casserait sa période) : pour une
image à bande limitée, la sortie est alors l'échantillonnage exact aux
positions cibles, à 10⁻¹³ près (vérifié par test pour des facteurs
entiers et non entiers, tailles paires et impaires).

### 5.1 `fourier` — gabarit sans repliement

`G(fy, fx) = g(fy) g(fx)`, `g` valant 1 jusqu'à `passband` (0,4), 0 à
partir de `stopband` (0,5), avec une transition en cosinus surélevé. Le
gabarit s'annule **avant** le Nyquist : aucun repliement (l'énergie
au-delà du Nyquist est nulle par construction), au prix d'un gain nul au
Nyquist. `passband=0.45, stopband=0.55` donne un gain de 1/2 au Nyquist,
comme les sinus cardinaux, avec une transition deux fois plus étroite
que `sinc_kaiser`. Coût : ≈ 0,4 s pour une image 2048 × 2048, soit
environ 20 fois une méthode séparable.

Limite : traitement global de l'image. Un traitement par tuiles
demanderait un recouvrement (voir [§9](conception.md#9-pistes)).

## 6. Simuler un capteur : `mtf_target`

But : produire l'image qu'aurait donnée un instrument de MTF connue, à la
résolution cible — c'est la méthode « physique » de réduction.

`G(fy, fx) = T(fy, fx) / max(S(fy/ky, fx/kx), floor)` :

- `target` (`T`, requis) : MTF visée, en cycles par pixel **cible** ;
- `source` (`S`, optionnel) : MTF de l'image d'entrée, en cycles par
  pixel **source** — elle est compensée, pour que la sortie ait la MTF
  `T` et non le produit `S · T` ;
- `floor` (0,05) : plancher de `S`, qui borne l'amplification là où la
  MTF source est faible (l'image d'entrée n'y contient presque que du
  bruit).

Contrairement aux autres méthodes, `T` n'a aucune raison de s'annuler au
Nyquist cible : **la sortie replie comme le ferait le vrai capteur**
(vérifié par test : une sinusoïde à `f = 0,7` ressort, repliée en 0,3,
avec le gain `T(0,7)`). C'est voulu.

Modèles fournis (`eoqual.resampling.mtf`), à combiner par `product` :

| Modèle | MTF | Au Nyquist |
|---|---|---|
| `gaussian_mtf(v)` | `exp(-a f²)`, isotrope | `v` |
| `pixel_mtf(c)` | `abs(sinc(c fy) · sinc(c fx))`, pixel carré de côté `c` | 0,64 pour `c = 1` |
| `diffraction_mtf(fc)` | pupille circulaire non obstruée, coupure `fc` | — |

Toute fonction `mtf(fy, fx)` convient (MTF mesurée interpolée, par
exemple).
