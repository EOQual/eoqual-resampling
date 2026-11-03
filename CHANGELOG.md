# Change Log

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
