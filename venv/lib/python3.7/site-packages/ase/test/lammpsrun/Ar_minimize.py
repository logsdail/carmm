from ase.calculators.lammpsrun import LAMMPS
from ase.cluster.icosahedron import Icosahedron
from ase.data import atomic_numbers,  atomic_masses
from numpy.testing import assert_allclose
from ase.optimize import LBFGS


ar_nc = Icosahedron('Ar', noshells=2)
ar_nc.cell = [[300, 0, 0], [0, 300, 0], [0, 0, 300]]
ar_nc.pbc = True

params = {}
params['pair_style'] = 'lj/cut 8.0'
params['pair_coeff'] = ['1 1 0.0108102 3.345']
params['masses'] = ['1 {}'.format(atomic_masses[atomic_numbers['Ar']])]

calc = LAMMPS(specorder=['Ar'], **params)

ar_nc.set_calculator(calc)

assert_allclose(ar_nc.get_potential_energy(), -0.468147667942117,
                atol=1e-4, rtol=1e-4)
assert_allclose(ar_nc.get_forces(), calc.calculate_numerical_forces(ar_nc),
                atol=1e-4, rtol=1e-4)

dyn = LBFGS(ar_nc, force_consistent=False)
dyn.run(fmax=1E-6)

assert_allclose(ar_nc.get_potential_energy(), -0.4791815886953914,
                atol=1e-4, rtol=1e-4)
assert_allclose(ar_nc.get_forces(), calc.calculate_numerical_forces(ar_nc),
                atol=1e-4, rtol=1e-4)
