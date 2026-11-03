# Contribuer à eoqual-resampling

Les contributions sont bienvenues, sous toutes leurs formes : signalement
de bug, correction, proposition de nouvelle méthode, implémentation d'une
méthode déjà proposée, amélioration de la documentation.

## Signaler un bug

Ouvrez une issue avec :
- la méthode concernée (`list_methods()` pour voir le registre) ;
- un exemple minimal reproduisant le problème (image ou scène
  synthétique, facteur, paramètres, code) ;
- le résultat obtenu vs. attendu.

## Proposer une méthode

Ouvrez une issue décrivant la méthode (nom, référence bibliographique,
ce qu'elle apporte par rapport au catalogue — voir
`REFERENCE_TECHNIQUE.md`). Pas besoin de savoir l'implémenter pour la
proposer.

## Implémenter une méthode

Convention du projet :

- une méthode **séparable** est un noyau 1D `Kernel(phi, radius,
  antialias)` dans `src/eoqual/resampling/kernels.py`, ses paramètres en
  unités **cibles** (pixels cibles, cycles par pixel cible) ;
- une méthode **spectrale** est une réponse `G(fy, fx)` construite dans
  `methods.py::build_response` ;
- elle est déclarée dans `config.py::METHODS_CONFIGS` (famille, moteur,
  paramètres par défaut, référence) et construite dans `methods.py` ;
- elle ne décale pas l'image et conserve une image uniforme (tests
  communs de `tests/test_methods.py`, qui la couvrent automatiquement) ;
- sa référence est ajoutée à `BIBLIOGRAPHY.md`.

Une brique qui sert à deux modules EOQual au moins va dans eoqual-core
(numpy et scipy seulement) ; une figure générique va dans eoqual-viz
(Plotly).

**Licence — condition de fusion, pas une formalité.**

- **Implémentation depuis une publication** : bienvenue sous MIT.
  L'en-tête du module cite les publications utilisées ; aucun code tiers
  n'est consulté ni repris.
- **Reprise d'un code publié ailleurs** : seulement si sa licence
  vérifiée est permissive, avec conservation de sa notice. Un code sans
  licence n'est pas acceptable.
- **Bibliothèque sous licence restrictive ou dépendance lourde** :
  seulement en extra pip optionnel, jamais importée à l'import du
  paquet, avec un avertissement dans `THIRD_PARTY_LICENSES.md`.

## Tests

Les tests ne s'appuient que sur des données synthétiques (géométrie et
spectre connus). Les tests longs portent le marqueur `slow`
(`pytest -m "not slow"` pour une passe rapide).

Mettre à jour `BENCHMARK.md` avec `examples/benchmark.py`.

## Documentation

Toute méthode nouvelle ou modifiée doit rester reflétée dans `README.md`
(catalogue), `REFERENCE_TECHNIQUE.md` (principe, implémentation,
limites), `THIRD_PARTY_LICENSES.md` (licence) et `BIBLIOGRAPHY.md`
(référence).
