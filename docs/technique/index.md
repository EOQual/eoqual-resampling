# Documents techniques

Principe, implémentation et limites de chaque méthode de réduction de
résolution. Les chapitres gardent une numérotation unique d'un document à
l'autre : un renvoi « §n » désigne toujours le même chapitre. Ces documents
décrivent l'état présent du code ; ce qui change d'une version à l'autre
est dans les [notes de version](../changelog.md).

| Document | Chapitres | Contenu |
|---|---|---|
| [Le problème et la géométrie de la grille cible](principe.md) | 1, 2 | Repliement, flou, décalage ; conventions de grille |
| [Méthodes de réduction](methodes.md) | 3, 4, 5, 6 | Moteur séparable, noyaux (dont sinc_kaiser par défaut), méthodes spectrales, simulation d'un capteur par sa MTF |
| [Diagnostic et cas réel](diagnostic.md) | 7, 8 | Indicateurs du compromis repliement / flou / oscillations ; orthophotographie de 5 cm à 20 cm |
| [Pistes, tests et choix de conception](conception.md) | 9, 10, 11 | Pistes, tests, choix de conception et variantes écartées |
| [Benchmark](benchmark.md) | | Comparaison chiffrée des méthodes |
| [Bibliographie](bibliographie.md) | | Références des méthodes |

```{toctree}
:hidden:
:maxdepth: 2

principe
methodes
diagnostic
conception
benchmark
bibliographie
```
