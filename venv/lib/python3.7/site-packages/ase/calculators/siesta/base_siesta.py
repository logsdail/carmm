# File solely for gently deprecating this BaseSiesta class.
import warnings
import numpy as np
from ase.calculators.siesta.siesta import Siesta


class BaseSiesta(Siesta):
    def __init__(self, *args, **kwargs):
        warnings.warn(
            "The BaseSiesta calculator class will no longer be supported. "
            "Use 'ase.calculators.siesta.Siesta in stead.",
            np.VisibleDeprecationWarning)
        Siesta.__init__(self, *args, **kwargs)
