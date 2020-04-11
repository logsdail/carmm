import numpy as np

from ase import Atoms, Atom
from ase.build import fcc111, fcc211, add_adsorbate, bulk, surface
import math

atoms = fcc211('Au', (3, 5, 8), vacuum=10.)
assert len(atoms) == 120

atoms = atoms.repeat((2, 1, 1))
assert np.allclose(atoms.get_distance(0, 130), 2.88499566724)

atoms = fcc111('Ni', (2, 2, 4), orthogonal=True)
add_adsorbate(atoms, 'H', 1, 'bridge')
add_adsorbate(atoms, Atom('O'), 1, 'fcc')
add_adsorbate(atoms, Atoms('F'), 1, 'hcp')

# The next test ensures that a simple string of multiple atoms cannot be used,
# which should fail with a KeyError that reports the name of the molecule due
# to the string failing to work with Atom().
failed = False
try:
    add_adsorbate(atoms, 'CN', 1, 'ontop')
except KeyError as e:
    failed = True
    assert e.args[0] == 'CN'
assert failed

# This test ensures that the default periodic behavior remains unchanged
cubic_fcc = bulk("Al", a=4.05, cubic=True)
surface_fcc = surface(cubic_fcc, (1,1,1), 3)

assert list(surface_fcc.pbc) == [True, True, False]
assert surface_fcc.cell[2][2] == 0

# This test checks the new periodic option
cubic_fcc = bulk("Al", a=4.05, cubic=True)
surface_fcc = surface(cubic_fcc, (1,1,1), 3, periodic=True)

assert (list(surface_fcc.pbc) == [True, True, True])
expected_length = 4.05*3**0.5 #for FCC with a=4.05
assert math.isclose(surface_fcc.cell[2][2], expected_length)

