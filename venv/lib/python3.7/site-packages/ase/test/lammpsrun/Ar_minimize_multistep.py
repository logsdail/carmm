from ase.calculators.lammpsrun import LAMMPS
from ase.cluster.icosahedron import Icosahedron
from ase.data import atomic_numbers,  atomic_masses
from numpy.testing import assert_allclose

ar_nc = Icosahedron('Ar', noshells=2)
ar_nc.cell = [[300, 0, 0], [0, 300, 0], [0, 0, 300]]
ar_nc.pbc = True

params = {}
params['pair_style'] = 'lj/cut 8.0'
params['pair_coeff'] = ['1 1 0.0108102 3.345']
params['masses'] = ['1 {}'.format(atomic_masses[atomic_numbers['Ar']])]

calc = LAMMPS(specorder=['Ar'], **params)

ar_nc.set_calculator(calc)
F1_numer = calc.calculate_numerical_forces(ar_nc)

assert_allclose(ar_nc.get_potential_energy(), -0.468147667942117,
                atol=1e-4, rtol=1e-4)
assert_allclose(ar_nc.get_forces(), calc.calculate_numerical_forces(ar_nc),
                atol=1e-4, rtol=1e-4)

params['minimize'] = '1.0e-15 1.0e-6 2000 4000'   # add minimize
calc.parameters = params

# set_atoms=True to read final coordinates after minimization
calc.run(set_atoms=True)

# get final coordinates after minimization
ar_nc.set_positions(calc.atoms.positions)

assert_allclose(ar_nc.get_potential_energy(), -0.4791815887032201,
                atol=1e-4, rtol=1e-4)
assert_allclose(ar_nc.get_forces(), calc.calculate_numerical_forces(ar_nc),
                atol=1e-4, rtol=1e-4)
