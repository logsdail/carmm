import numpy as np

from ase.build import bulk
from ase.constraints import FixScaled
from ase.calculators.emt import EMT

a = bulk("Ni", cubic=True)
a.set_calculator(EMT())

pos = a.get_positions()

a.set_constraint(
    FixScaled(
        a.cell,
        0,
    )
)

a.set_positions(
    pos * 1.01
)

assert np.sum(np.abs(a.get_forces()[0])) < 1e-12
assert np.sum(np.abs(a.get_positions() - pos)[0]) < 1e-12
assert np.sum(np.abs(a.get_positions() - pos*1.01)[1:].flatten()) < 1e-12

