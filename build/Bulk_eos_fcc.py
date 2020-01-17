import numpy as np
from ase import Atoms
from ase.io.trajectory import Trajectory
from ase.calculators.emt import EMT
from ase.calculators.aims import Aims
from ase.visualize import view

a = 4.0 # approximate lattice constant
b = a / 2

calc = Aims(xc='pbe',
           spin='none',
           k_grid=(9,9,9),
           vdw_correction_hirshfeld="True",
           relativistic=('atomic_zora','scalar'),
           #use_dipole_correction='True',
           compute_forces="true",
           output=['mulliken'],
          # elsi_restart=("write",1)
           )

bulk = Atoms('Pd',
           cell=[(0, b, b), (b, 0, b), (b, b, 0)],
           pbc=1,
           calculator=calc)  # use EMT potential

"""
## bulk structure can also be defined manually

bulk = Atoms('Pd2Cu2',
              scaled_positions=[(0, 0, 0),
                                (0.5, 0.5, 0),
                                (0.5, 0, 0.5),
                                (0, 0.5, 0.5),
                                ],

              cell=[a, a, a],
              pbc=True,
              calculator=EMT())

view(bulk)
"""

cell = bulk.get_cell()
traj = Trajectory('bulk.traj', 'w')
for x in np.linspace(0.90, 1.10, 5):
    bulk.set_cell(cell * x, scale_atoms=True)
    bulk.get_potential_energy()
    traj.write(bulk)

########## EOS #############

from ase.io import read
from ase.units import kJ
from ase.eos import EquationOfState
configs = read('bulk.traj@0:5')  # read 10 configurations
# Extract volumes and energies:
volumes = [bulk.get_volume() for bulk in configs]
energies = [bulk.get_potential_energy() for bulk in configs]
eos = EquationOfState(volumes, energies, eos='birchmurnaghan')
v0, e0, B = eos.fit()

print("Volume: ", v0)
print("FCC Lattice parameter a0: ", (4*v0)**(1/3), "Angstrom")
#print("FCC Lattice parameter a0: ", (4*v0/NUMBER_OF_ATOMS_IN_UNIT_CELL)**(1/3), "Angstrom")
print("Bulk modulus: ", B / kJ * 1.0e24, 'GPa')
eos.plot('bulk-eos.png')

#view('bulk.traj')
