# Change Log

## Non publié

- Référence technique découpée en documents spécialisés dans
  `docs/technique/` (principe, méthodes, diagnostic et cas réel,
  conception) ; les chapitres gardent leur numéro.
- Documentation réorganisée : documents techniques dans
  `docs/technique/`, publiés avec l'API et les notes de version par
  Sphinx et MyST (`docs/`). Les documents décrivent l'état présent ; ce
  qui change d'une version à l'autre est dans ce fichier, et les
  variantes de méthode écartées sont décrites au présent dans les
  documents.
- Exemple sur données réelles (`examples/orthophoto.py`) :
  orthophotographie de Toulouse (Licence Ouverte v2.0) réduite de 5 cm à
  20 cm par les douze méthodes ; indicateurs (contraste des détails fins,
  repliement sur le spectre réel, décalage, écart à la moyenne par blocs)
  et figures commentées dans `docs/technique/diagnostic.md` §8.
- Figure PNG des conventions de grille (`examples/figures.py`).

## 0.1.0 (2026-11-03)

- Première version : réduction de résolution d'images par 12 méthodes
  (3 de référence sans préfiltrage ; 7 noyaux anti-repliement séparables,
  dont le sinus cardinal apodisé par une gaussienne et le sinc fenêtré
  par Kaiser, méthode par défaut ; 2 méthodes spectrales : gabarit sans
  repliement et MTF cible),
  facteur quelconque par axe, conventions de grille explicites,
  conservation exacte des niveaux.
- Diagnostic : fonction de transfert, repliement, transition,
  dépassement, variance du bruit ; figure Plotly.
- Benchmark synthétique (`examples/benchmark.py`, `BENCHMARK.md`),
  documentation Sphinx (GitHub Pages) et tests GitHub Actions.
