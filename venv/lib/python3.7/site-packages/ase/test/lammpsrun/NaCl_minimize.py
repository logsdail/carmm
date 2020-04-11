from ase.calculators.lammpsrun import LAMMPS
from ase.spacegroup import crystal
from ase.data import atomic_numbers,  atomic_masses
from ase.optimize import QuasiNewton
from ase.constraints import UnitCellFilter
from numpy.testing import assert_allclose


a = 6.15
n = 4
nacl = crystal(['Na', 'Cl'], [(0, 0, 0), (0.5, 0.5, 0.5)], spacegroup=225,
               cellpar=[a, a, a, 90, 90, 90]).repeat((n, n, n))

# Buckingham parameters from
# https://physics.stackexchange.com/questions/250018


pair_style = 'buck/coul/long 12.0'
pair_coeff = ['1 1 3796.9 0.2603 124.90']
pair_coeff += ['2 2 1227.2 0.3214 124.90']
pair_coeff += ['1 2 4117.9 0.3048 0.0']
masses = ['1 {}'.format(atomic_masses[atomic_numbers['Na']]),
          '2 {}'.format(atomic_masses[atomic_numbers['Cl']])]

calc = LAMMPS(specorder=['Na', 'Cl'],
              pair_style=pair_style,
              pair_coeff=pair_coeff,
              masses=masses,
              atom_style='charge',
              kspace_style='pppm 1.0e-5',
              keep_tmp_files=True,
              )

for a in nacl:
    if a.symbol == 'Na':
        a.charge = +1.
    else:
        a.charge = -1.

nacl.set_calculator(calc)

assert_allclose(nacl.get_potential_energy(), -1896.216737561538,
                atol=1e-4, rtol=1e-4)

E = nacl.get_potential_energy()

ucf = UnitCellFilter(nacl)
dyn = QuasiNewton(ucf, force_consistent=False)
dyn.run(fmax=1.0E-2)

assert_allclose(nacl.get_potential_energy(), -1897.208861729178,
                atol=1e-4, rtol=1e-4)
