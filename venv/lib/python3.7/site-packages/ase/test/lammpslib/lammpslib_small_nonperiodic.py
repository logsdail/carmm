# test that lammpslib handle nonperiodic cases where the cell size
# in some directions is small (for example for a dimer).
import numpy as np
from ase.calculators.lammpslib import LAMMPSlib
from ase import Atoms

cmds = ["pair_style eam/alloy",
        "pair_coeff * * NiAlH_jea.eam.alloy Ni H"]
lammps = LAMMPSlib(lmpcmds=cmds,
                   atom_types={'Ni': 1, 'H': 2},
                   log_file='test.log', keep_alive=True)
a = 2.0
dimer = Atoms("NiNi", positions=[(0, 0, 0), (a, 0, 0)],
              cell=(1000*a, 1000*a, 1000*a), pbc=(0, 0, 0))
dimer.set_calculator(lammps)

energy_ref = -1.10756669119
energy = dimer.get_potential_energy()
print("Computed energy: {}".format(energy))
np.testing.assert_allclose(energy, energy_ref, atol=1e-4, rtol=1e-4)

np.set_printoptions(precision=16)
forces_ref = np.array([[-0.9420162329811532, 0., 0.],
                       [+0.9420162329811532, 0., 0.]])
forces = dimer.get_forces()
print(np.array2string(forces))
np.testing.assert_allclose(forces, forces_ref, atol=1e-4, rtol=1e-4)
