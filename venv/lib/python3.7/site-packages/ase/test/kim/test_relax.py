"""
Test that a static relaxation that requires multiple neighbor list
rebuilds can be carried out successfully.  This is verified by relaxing
an icosahedral cluster of atoms and checking that the relaxed energy
matches a known precomputed value for an example model.
"""
import numpy as np
from ase.cluster import Icosahedron
from ase.calculators.kim import KIM
from ase.optimize import BFGS

energy_ref = -0.5420939378624228  # eV

# Create structure and calculator
atoms = Icosahedron("Ar", latticeconstant=3.0, noshells=2)
calc = KIM("ex_model_Ar_P_Morse_07C")
atoms.set_calculator(calc)

opt = BFGS(atoms, logfile=None)
opt.run(fmax=0.05)

assert np.isclose(atoms.get_potential_energy(), energy_ref)
