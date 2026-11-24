# Pistes, tests et choix de conception

## 9. Pistes

- **Agrandissement** (facteur < 1) : les noyaux et le moteur séparable le
  permettent déjà (dilatation par `max(k, 1)`, coupure au Nyquist
  source) ; il reste à l'exposer, à le tester et à documenter les noyaux
  d'interpolation : leur MTF et la variance du bruit dépendent de la
  phase d'interpolation (à un demi-pixel, tout noyau symétrique annule le
  Nyquist et corrèle le bruit). Il est refusé explicitement.
- **Pixels invalides** (NaN, masque de validité) : normalisation des
  poids sur les seuls pixels valides (méthodes séparables).
- **Grandes images** : traitement par tuiles avec recouvrement du support
  (`eoqual.core.tiling`), lecture et écriture par blocs.
- **Splines aux moindres carrés** (Unser, Aldroubi & Eden, 1995) :
  réduction optimale au sens L2 dans l'espace des splines.
- **Noyau BCO** d'Orfeo ToolBox et autres noyaux des chaînes de production,
  pour comparer à ce qu'elles produisent.
- **Brique commune** : si un deuxième module EOQual a besoin de réduire
  une image (simulateur d'eoqual-measure, pyramide de la déconvolution
  aveugle d'eoqual-deconvolution), le noyau séparable et la géométrie
  pourront passer dans eoqual-core.

## 10. Tests

Données synthétiques uniquement (`pytest`, 220 tests, < 2 s) :

- **contrat** : forme, type, valeurs finies, facteurs entiers, non
  entiers et différents par axe, `shape=`, multi-bandes, bords, erreurs
  explicites ;
- **niveaux** : image uniforme conservée exactement par toutes les
  méthodes ; `box` = moyenne par blocs pour `k` entier ;
- **géométrie** : centroïde d'une tache à moins de 0,002 pixel cible de
  la position attendue, pour les deux grilles et `k` = 2 ; 2,5 ; 3 ; 4 ;
- **spectre** : gain annoncé = gain mesuré sur sinusoïdes (repliées ou
  non) ; repliement total des méthodes de référence, rejet au-delà du
  Nyquist des méthodes anti-repliement ; exactitude de `fourier` sur une
  image à bande limitée ; `mtf_target` reproduit la MTF visée, son
  repliement et compense la MTF source ;
- **bruit** : variance mesurée = variance prédite à 8 % près ;
- **exemples** : `quickstart.py`, la figure de la documentation
  (`figures.py`, si matplotlib est installé), la chaîne d'évaluation de
  `orthophoto.py` sur une image RGB synthétique (si scikit-image est
  installé) et (marqueur `slow`) les tableaux du benchmark.

## 11. Choix de conception et variantes écartées

Le noyau `sinc_gaussian` reproduit le sinus cardinal apodisé des notes de
travail de l'auteur (sinc × gaussienne de 12 pixels cibles, largeur à
mi-hauteur 6 ; écart < 10⁻¹³ sur une image de test, à la géométrie près).
Le module n'en reprend pas le code : il est écrit à partir des publications
de la [bibliographie](bibliographie.md), et évite les pièges suivants.

- **Noyau de taille paire et convolution « valid »** : un noyau de `12 k`
  échantillons combiné à une convolution en mode « valid » échantillonne
  en `i k - 1` au lieu de `i k`, soit un décalage d'**un pixel source**.
- **Convention de grille implicite** (`"origin"`) : écart de `(k - 1)/2`
  pixel source avec l'emprise du binning. Les deux conventions sont
  explicites ([§2](principe.md#2-géométrie-de-la-grille-cible)), `"center"` par défaut.
- **Interpolation linéaire après filtrage** pour un facteur non entier :
  flou variable selon la phase. Le noyau est évalué directement aux
  positions cibles.
- **Convolution 2D non séparable** (`(12k)²` poids par pixel) : deux
  passes 1D la remplacent.
- **Poids résiduel au-delà du rayon** : `sinc_gaussian`, non nul au-delà
  de son rayon, recevrait un poids résiduel d'un seul côté du support
  (décalage de 0,004 pixel cible à `k = 3`) : troncature stricte au rayon
  pour tous les noyaux.
- **Recadrage à l'emprise `M·k`** (chemin rapide des méthodes
  spectrales) : il casse la période d'une image périodique ; il n'est fait
  que sur une image décomposée en `p + s`.
- **Méthode par défaut.** `sinc_kaiser` plutôt que `sinc_gaussian` : même
  netteté, repliement mesuré divisé par trois sur la mire de Fresnel, au
  prix d'un support 1,5 fois plus long (9 pixels cibles au lieu de 6,
  temps de calcul négligeable dans les deux cas).
