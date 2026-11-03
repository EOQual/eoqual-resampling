# Bibliographie

Publications à partir desquelles les méthodes ont été implémentées, et
références des choix de conception. Aucun code tiers n'a été consulté
pour écrire les méthodes de ce paquet (voir
[THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md)).

## Échantillonnage et repliement

- C. E. Shannon, "Communication in the presence of noise", *Proc. IRE*
  37(1), 10-21, 1949. — théorème d'échantillonnage, Nyquist.
- A. V. Oppenheim, R. W. Schafer, *Discrete-Time Signal Processing*,
  3ᵉ éd., Pearson, 2010, chap. 4 et 7. — changement de cadence,
  conception de filtres RIF par fenêtrage, fenêtre de Kaiser.
- K. Turkowski, "Filters for common resampling tasks", in A. S. Glassner
  (éd.), *Graphics Gems*, Academic Press, 1990, 147-165. — comparaison
  des noyaux de réduction (boîte, gaussienne, Lanczos).
- P. Thévenaz, T. Blu, M. Unser, "Interpolation revisited", *IEEE Trans.
  Medical Imaging* 19(7), 739-758, 2000. — noyaux d'interpolation,
  conventions de grille.

## Noyaux

- R. G. Keys, "Cubic convolution interpolation for digital image
  processing", *IEEE Trans. Acoustics, Speech, and Signal Processing*
  29(6), 1153-1160, 1981. — `bicubic`, `bicubic_aa`.
- D. P. Mitchell, A. N. Netravali, "Reconstruction filters in
  computer graphics", *SIGGRAPH Computer Graphics* 22(4), 221-228,
  1988. — `mitchell`, famille de cubiques `(B, C)`.
- C. E. Duchon, "Lanczos filtering in one and two dimensions", *Journal
  of Applied Meteorology* 18(8), 1016-1022, 1979. — `lanczos`.
- J. F. Kaiser, "Nonrecursive digital filter design using the I0-sinh
  window function", *Proc. IEEE International Symposium on Circuits and
  Systems*, 20-23, 1974. — `sinc_kaiser` (β et longueur du filtre en
  fonction de l'atténuation et de la transition).
- P. J. Burt, E. H. Adelson, "The Laplacian pyramid as a compact image
  code", *IEEE Trans. Communications* 31(4), 532-540, 1983. — réduction
  par filtre gaussien (`gaussian`).

## Méthodes spectrales

- L. Moisan, "Periodic plus smooth image decomposition", *Journal of
  Mathematical Imaging and Vision* 39(2), 161-179, 2011. — décomposition
  `p + s` (eoqual-core), traitement des bords de `fourier` et
  `mtf_target`.

## Évolutions envisagées

- M. Unser, A. Aldroubi, M. Eden, "Enlargement or reduction of digital
  images with minimum loss of information", *IEEE Trans. Image
  Processing* 4(3), 247-258, 1995. — réduction aux moindres carrés dans
  l'espace des splines.
