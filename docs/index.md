# eoqual-resampling

Réduction de résolution sans repliement ni décalage géométrique, pour un facteur quelconque (non entier, différent par axe) : 12 méthodes, des interpolations sans préfiltrage (références) à l'agrégation de pixels, aux sinus cardinaux fenêtrés (gaussienne, Kaiser réglé par l'atténuation visée) et aux méthodes spectrales (gabarit sans repliement, simulation d'un capteur par sa MTF) ; conventions de grille explicites, conservation exacte des niveaux, diagnostic du compromis repliement / flou / oscillations et benchmark.

- [Documents techniques](technique/index.md) : principe, mise en œuvre,
  validation et limites des méthodes.
- [Notes de version](changelog.md) : ce qui change d'une version à l'autre.
- [Référence de l'API](api/index) : générée à partir du code source.

```{toctree}
:hidden:
:maxdepth: 2

technique/index
changelog
api/index
```
