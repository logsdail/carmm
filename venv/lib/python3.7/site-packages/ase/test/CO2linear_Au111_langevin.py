from math import pi, cos, sin
import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixLinearTriatomic
from ase.md import Langevin
from ase.build import fcc111, add_adsorbate
import ase.units as units
"""Test Langevin with constraints for rigid linear
triatomic molecules""" 

rng = np.random.RandomState(0)
eref = 3.1356 

zpos = cos(134.3 / 2.0 * pi / 180.0) * 1.197
xpos = sin(134.3 / 2.0 * pi / 180.0) * 1.19
co2 = Atoms('COO', positions=[(-xpos + 1.2, 0, -zpos),
                              (-xpos + 1.2, -1.1, -zpos),
                              (-xpos + 1.2, 1.1, -zpos)])

slab = fcc111('Au', size=(2, 2, 4), vacuum=2 * 5, orthogonal=True)
slab.center()
add_adsorbate(slab, co2, 1.5, 'bridge')
slab.set_pbc((True, True, False))
d0 = co2.get_distance(-3, -2)
d1 = co2.get_distance(-3, -1)
d2 = co2.get_distance(-2, -1)

calc = EMT()
slab.set_calculator(calc)
constraint = FixLinearTriatomic(triples=[(-2, -3, -1)])
slab.set_constraint(constraint)

fr = 0.1
dyn = Langevin(slab, 2.0 * units.fs,
               300 * units.kB, fr, 
               trajectory='langevin_%.1f.traj' % fr,
               logfile='langevin_%.1f.log' % fr,
               loginterval=20, rng=rng)
dyn.run(100)

# Check that the temperature is within a reasonable range 
T = slab.get_temperature()
assert T > 100
assert T < 500

# Check that the constraints work
assert abs(slab.get_distance(-3, -2, mic=1) - d0) < 1e-9
assert abs(slab.get_distance(-3, -1, mic=1) - d1) < 1e-9
assert abs(slab.get_distance(-2, -1, mic=1) - d2) < 1e-9

# If the energy differs from the reference energy
# it is most probable that the redistribution of 
# random forces in Langevin is not working properly
assert abs(slab.get_potential_energy() - eref) < 1e-4
