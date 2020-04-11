# Copyright 2008, 2009 CAMd
# (see accompanying license files for details).

"""Atomic Simulation Environment."""

import sys
from distutils.version import LooseVersion
import numpy as np


if sys.version_info[0] == 2:
    raise ImportError('ASE requires Python3. This is Python2.')


if LooseVersion(np.__version__) < '1.9':
    raise ImportError(
        'ASE needs NumPy-1.9.0 or later. You have: %s' % np.__version__)


__all__ = ['Atoms', 'Atom']
__version__ = '3.19.1'


from ase.atom import Atom
from ase.atoms import Atoms

# import ase.parallel early to avoid circular import problems when
# ase.parallel does "from gpaw.mpi import world":
import ase.parallel  # noqa
ase.parallel  # silence pyflakes
