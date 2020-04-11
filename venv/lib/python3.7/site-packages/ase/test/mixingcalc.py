"""This test checks the basic functionality of the MixingCalculators.
The example system is based on the SinglePointCalculator test case.
"""
import numpy as np

from ase.build import fcc111
from ase.calculators.emt import EMT
from ase.calculators.mixing import SumCalculator, LinearCombinationCalculator, AverageCalculator, MixedCalculator
from ase.constraints import FixAtoms

# Calculate reference values:
atoms = fcc111('Cu', (2, 2, 1), vacuum=10.)
atoms[0].x += 0.2

# First run the test with EMT similarly to the test of the single point calculator.
calc = EMT()
atoms.set_calculator(calc)
forces = atoms.get_forces()

# SumCalculator: Alternative ways to associate a calculator with an atoms object.
atoms1 = atoms.copy()
calc1 = SumCalculator([EMT(), EMT()])
atoms1.set_calculator(calc1)

atoms2 = atoms.copy()
calc2 = SumCalculator(calcs=[EMT(), EMT()], atoms=atoms2)

# Check the results.
assert np.isclose(2 * forces, atoms1.get_forces()).all()
assert np.isclose(2 * forces, atoms2.get_forces()).all()

# testing  step
atoms1[0].x += 0.2
assert not np.isclose(2 * forces, atoms1.get_forces()).all()

# Check constraints
atoms1.set_constraint(FixAtoms(indices=[atom.index for atom in atoms]))
assert np.isclose(0, atoms1.get_forces()).all()

# AverageCalculator:

atoms1 = atoms.copy()
calc1 = AverageCalculator([EMT(), EMT()])
atoms1.set_calculator(calc1)

# LinearCombinationCalculator:

atoms2 = atoms.copy()
calc2 = LinearCombinationCalculator([EMT(), EMT()], weights=[.5, .5], atoms=atoms2)

# Check the results (it should be the same because it is tha average of the same values).
assert np.isclose(forces, atoms1.get_forces()).all()
assert np.isclose(forces, atoms2.get_forces()).all()

# testing  step
atoms1[0].x += 0.2
assert not np.isclose(2 * forces, atoms1.get_forces()).all()

try:
    calc1 = LinearCombinationCalculator([], [])
except ValueError:
    assert True

try:
    calc1 = AverageCalculator([])
except ValueError:
    assert True


# test  MixedCalculator and energy contributions
w1, w2 = 0.78, 0.22
atoms1 = atoms.copy()
atoms1.set_calculator(EMT())
E_tot = atoms1.get_potential_energy()

calc1 = MixedCalculator(EMT(), EMT(), w1, w2)
E1, E2 = calc1.get_energy_contributions(atoms1)
assert np.isclose(E1, E_tot)
assert np.isclose(E2, E_tot)
