import numpy as np

from ase.test.vasp import installed2 as installed
from ase.io import read
from ase.optimize import BFGS
from ase.build import fcc111
from ase.constraints import FixAtoms
from ase.calculators.vasp import Vasp2 as Vasp

assert installed()


def create_slab_with_constraints():
    slab = fcc111('Al', size=(1, 1, 3), periodic=True)
    slab.center(vacuum=4, axis=2)
    con = FixAtoms(indices=[0, 1])
    slab.set_constraint(con)
    return slab


def ase_relax():
    slab = create_slab_with_constraints()
    calc = Vasp(xc='LDA', ediffg=-1e-3, lwave=False, lcharg=False)
    slab.set_calculator(calc)
    opt = BFGS(slab, logfile=None)
    opt.run(fmax=0.1, steps=3)

    init_slab = create_slab_with_constraints()
    res = read('OUTCAR')
    assert np.allclose(res.positions[0], init_slab.positions[0])
    assert not np.allclose(res.positions[2], init_slab.positions[2])


def vasp_relax():
    slab = create_slab_with_constraints()
    calc = Vasp(xc='LDA', isif=0, nsw=3, ibrion=1,
                ediffg=-1e-3, lwave=False, lcharg=False)
    calc.calculate(slab)

    init_slab = create_slab_with_constraints()
    res = read('OUTCAR')
    assert np.allclose(res.positions[0], init_slab.positions[0])
    assert not np.allclose(res.positions[2], init_slab.positions[2])


def run_tests():
    ase_relax()
    vasp_relax()


run_tests()
