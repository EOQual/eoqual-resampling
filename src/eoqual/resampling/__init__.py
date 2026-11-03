"""
eoqual.resampling — réduction de résolution d'images sans repliement ni
déformation géométrique, pensée pour l'observation de la Terre.

>>> from eoqual.resampling import resample, list_methods
>>> small = resample(image, 3)                         # sinc fenêtré par Kaiser, facteur 3
>>> small = resample(image, 2.5)                       # facteur non entier
>>> from eoqual.resampling.mtf import gaussian_mtf, pixel_mtf, product
>>> sim = resample(image, 4, method="mtf_target",
...                target=product(gaussian_mtf(0.3), pixel_mtf()))

Voir ``README.md`` et ``REFERENCE_TECHNIQUE.md``.
"""

__version__ = "0.1.0"
__title__ = "eoqual.resampling"
__description__ = "Image downsampling methods for Earth observation"
__url__ = "https://github.com/EOQual/eoqual-resampling"
__uri__ = __url__
__doc__ = __description__ + " <" + __uri__ + ">"
__author__ = "Olivier Amram"
__email__ = "olivier.amram@free.fr"
__license__ = "MIT"
__copyright__ = "Copyright (c) 2026 Olivier Amram"

from .analysis import kernel_metrics, transfer_function  # noqa: F401
from .config import DEFAULT_METHOD, METHODS_CONFIGS  # noqa: F401
from .runner import list_methods, resample  # noqa: F401

__all__ = [
    "resample",
    "list_methods",
    "transfer_function",
    "kernel_metrics",
    "METHODS_CONFIGS",
    "DEFAULT_METHOD",
]
