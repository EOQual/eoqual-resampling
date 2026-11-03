# Licences tierces — inventaire et statut

> **Avertissement.** Cet inventaire est tenu au meilleur effort, à partir
> des fichiers LICENSE et des métadonnées PyPI des paquets, vérifiés le
> 2026-11-03. Si une erreur y est identifiée, elle sera corrigée dès que
> signalée — et si nécessaire, l'intégration concernée sera retirée
> plutôt que sa licence « arrangée » a posteriori. Pour signaler un
> problème, voir `CONTRIBUTING.md`.

Ce document recense le code qui fait réellement le calcul et sa licence,
au regard de la licence MIT d'`eoqual-resampling` (`license.txt`).

**Principe** : **aucun code tiers n'est copié dans ce dépôt.** Toutes les
méthodes sont des implémentations originales, écrites à partir des
publications de [BIBLIOGRAPHY.md](BIBLIOGRAPHY.md).

## 1. Méthodes

Les 12 méthodes (`nearest`, `bilinear`, `bicubic`, `box`, `gaussian`,
`bicubic_aa`, `mitchell`, `lanczos`, `sinc_gaussian`, `sinc_kaiser`,
`fourier`, `mtf_target`) : 🟢 MIT (ce paquet).

## 2. Dépendances

| Paquet | Usage | Licence |
|---|---|---|
| numpy | calcul | 🟢 BSD-3-Clause |
| scipy | matrices creuses, fonction de Bessel `i0` | 🟢 BSD-3-Clause |
| plotly | figure des fonctions de transfert | 🟢 MIT |
| rich | `list_methods` | 🟢 MIT |
| eoqual-core | décomposition périodique + lisse | 🟢 MIT |
| eoqual-viz | `plot_curves` | 🟢 MIT |

Aucune dépendance sous licence restrictive, aucun extra optionnel.
