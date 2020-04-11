import numpy as np
from numpy.testing import assert_allclose
from ase.calculators.lammpsrun import LAMMPS
from ase.build import bulk
from ase.test.eam_pot import Pt_u3
from ase.constraints import ExpCellFilter
from ase.optimize import BFGS

# (For now) reuse eam file stuff from other lammps test:
pot_fn = 'Pt_u3.eam'
f = open(pot_fn, 'w')
f.write(Pt_u3)
f.close()
params = {}
params['pair_style'] = 'eam'
params['pair_coeff'] = ['1 1 {}'.format(pot_fn)]
calc = LAMMPS(specorder=['Pt'], files=[pot_fn], **params)

rng = np.random.RandomState(17)

atoms = bulk('Pt') * (2, 2, 2)
atoms.rattle(stdev=0.1)
atoms.cell += 2 * rng.rand(3, 3)
atoms.calc = calc

assert_allclose(atoms.get_stress(), calc.calculate_numerical_stress(atoms),
                atol=1e-4, rtol=1e-4)

opt = BFGS(ExpCellFilter(atoms), trajectory='opt.traj')
for i, _ in enumerate(opt.irun(fmax=0.05)):
    pass

cell1_ref = np.array([
    [0.16298762, 3.89912471, 3.92825365],
    [4.21007577, 0.63362427, 5.04668170],
    [4.42895706, 3.29171414, 0.44623618]])


assert_allclose(np.asarray(atoms.cell), cell1_ref, atol=1e-4, rtol=1e-4)
assert_allclose(atoms.get_stress(), calc.calculate_numerical_stress(atoms),
                atol=1e-4, rtol=1e-4)

assert i < 80, 'Expected 59 iterations, got many more: {}'.format(i)
