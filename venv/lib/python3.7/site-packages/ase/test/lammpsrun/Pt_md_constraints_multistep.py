from ase.calculators.lammpsrun import LAMMPS
from numpy.testing import assert_allclose
from ase.test.eam_pot import Pt_u3
from ase.build import fcc111
import os


pot_fn = 'Pt_u3.eam'
f = open(pot_fn, 'w')
f.write(Pt_u3)
f.close()

slab = fcc111('Pt', size=(2, 2, 5), vacuum=30.0)
# We use fully periodic boundary conditions because the Lammpsrun
# calculator does not know if it can convert the cell correctly with
# mixed ones and will give a warning.
slab.pbc = 1

params = {}
params['pair_style'] = 'eam'
params['pair_coeff'] = ['1 1 {}'.format(pot_fn)]

calc = LAMMPS(specorder=['Pt'], files=[pot_fn], **params)
slab.set_calculator(calc)

assert_allclose(slab.get_potential_energy(), -110.3455014595596,
                atol=1e-4, rtol=1e-4)

params['group'] = ['lower_atoms id '
                   + ' '.join([str(i+1) for i,
                              tag in enumerate(slab.get_tags()) if tag >= 4])]
params['fix'] = ['freeze_lower_atoms lower_atoms setforce 0.0 0.0 0.0']
params['run'] = 100
params['timestep'] = 0.0005
params['dump_period'] = 10
params['write_velocities'] = True
calc.parameters = params
# set_atoms=True to read final coordinates and velocities after NVE simulation
calc.run(set_atoms=True)

new_slab = calc.atoms.copy()

Ek = calc.atoms.copy().get_kinetic_energy()
assert_allclose(Ek, 0.1014556059885532, atol=1e-4, rtol=1e-4)
assert_allclose(Ek, calc.thermo_content[-1]['ke'],
                atol=1e-4, rtol=1e-4)
assert_allclose(slab.get_potential_energy(), -110.4469605087525,
                atol=1e-4, rtol=1e-4)

os.remove(pot_fn)
