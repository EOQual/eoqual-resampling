"""
Registre central des méthodes de réduction de résolution de
``eoqual.resampling``.

Chaque entrée décrit :

- ``family`` : ``"référence"`` (interpolation sans préfiltrage : replie
  le spectre, pour comparaison), ``"anti-repliement"`` (noyau local
  dilaté du facteur) ou ``"spectral"`` (réponse en fréquence imposée) ;
- ``backend`` : ``"separable"`` ou ``"spectral"`` ;
- ``params`` : paramètres et valeurs par défaut (unités : pixels ou
  cycles par pixel **cibles**) ;
- ``license`` : licence du code qui fait le calcul ;
- ``reference`` : référence courte (complète dans ``BIBLIOGRAPHY.md``).
"""

from collections import OrderedDict

_OWN = "MIT (ce paquet)"

METHODS_CONFIGS: OrderedDict = OrderedDict(
    {
        "nearest": {
            "family": "référence",
            "backend": "separable",
            "params": {},
            "license": _OWN,
            "reference": "plus proche voisin, sans préfiltrage",
        },
        "bilinear": {
            "family": "référence",
            "backend": "separable",
            "params": {},
            "license": _OWN,
            "reference": "interpolation linéaire, sans préfiltrage",
        },
        "bicubic": {
            "family": "référence",
            "backend": "separable",
            "params": {"a": -0.5},
            "license": _OWN,
            "reference": "Keys (1981), sans préfiltrage",
        },
        "box": {
            "family": "anti-repliement",
            "backend": "separable",
            "params": {},
            "license": _OWN,
            "reference": "agrégation par surface (binning, « INTER_AREA »)",
        },
        "gaussian": {
            "family": "anti-repliement",
            "backend": "separable",
            "params": {"sigma": 0.5},
            "license": _OWN,
            "reference": "gaussienne puis échantillonnage (Burt & Adelson, 1983)",
        },
        "bicubic_aa": {
            "family": "anti-repliement",
            "backend": "separable",
            "params": {"a": -0.5},
            "license": _OWN,
            "reference": "cubique de Keys (1981) dilatée du facteur",
        },
        "mitchell": {
            "family": "anti-repliement",
            "backend": "separable",
            "params": {"b": 1.0 / 3.0, "c": 1.0 / 3.0},
            "license": _OWN,
            "reference": "Mitchell & Netravali (1988)",
        },
        "lanczos": {
            "family": "anti-repliement",
            "backend": "separable",
            "params": {"a": 3},
            "license": _OWN,
            "reference": "Duchon (1979) ; Turkowski (1990)",
        },
        "sinc_gaussian": {
            "family": "anti-repliement",
            "backend": "separable",
            "params": {"fwhm": 6.0, "radius": 6.0, "cutoff": 0.5},
            "license": _OWN,
            "reference": "sinus cardinal apodisé par une gaussienne",
        },
        "sinc_kaiser": {
            "family": "anti-repliement",
            "backend": "separable",
            "params": {"cutoff": 0.5, "attenuation": 60.0, "transition": 0.2, "radius": None},
            "license": _OWN,
            "reference": "fenêtre de Kaiser (Kaiser, 1974 ; Oppenheim & Schafer)",
        },
        "fourier": {
            "family": "spectral",
            "backend": "spectral",
            "params": {"passband": 0.4, "stopband": 0.5, "periodic": False},
            "license": _OWN,
            "reference": "gabarit en cosinus surélevé ; décomposition p + s (Moisan, 2011)",
        },
        "mtf_target": {
            "family": "spectral",
            "backend": "spectral",
            "params": {"target": "(requis)", "source": None, "floor": 0.05, "periodic": False},
            "license": _OWN,
            "reference": "rapport MTF cible / MTF source puis échantillonnage",
        },
    }
)

#: Méthode utilisée par défaut par :func:`eoqual.resampling.resample`.
DEFAULT_METHOD = "sinc_kaiser"
