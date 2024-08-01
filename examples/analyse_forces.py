#!/usr/bin/env python3

"""
Created on Fri 21/04/21

@author: Igor Kowalec

This script checks if forces data stored in the calculator attached to an Atoms
object meets the convergence criterion used in geometry optimisation.

Also works if system includes constraints, i.e. knows which atoms to check forces for.
"""

def test_analyse_forces():
    from ase.build import molecule
    from ase.calculators.emt import EMT
    from ase.optimize import BFGSLineSearch
    from ase.build import bulk
    from carmm.analyse.forces import is_converged
    from ase.constraints import FixAtoms

    # Ensures correct import
    from carmm.utils.python_env_check import ase_env_check
    if ase_env_check('3.23.0'):
        from ase.filters import FrechetCellFilter
    else:
        from ase.constraints import ExpCellFilter

    fmax = 0.01 # eV/Angstrom
    atoms = molecule("CO2")
    atoms.calc = EMT()

    # Returns False prior to optimisation (no forces in calculator)
    opt = BFGSLineSearch(atoms)
    # print(is_converged(atoms, fmax))
    assert not is_converged(atoms, fmax)

    # Returns False if optimisation has not reached to/below fmax
    opt.run(fmax=fmax*30)
    # print(is_converged(atoms, fmax))
    assert not is_converged(atoms, fmax)

    # Returns True if optimised to or below desired fmax with constraints
    c = FixAtoms(indices=[0,1])
    atoms.set_constraint(c)
    opt.run(fmax=fmax)
    # print(is_converged(atoms, fmax))
    assert is_converged(atoms, fmax)

    # Returns True if optimised to or below desired fmax without constraints
    atoms.set_constraint()
    opt.run(fmax=fmax)
    # print(is_converged(atoms, fmax))
    assert is_converged(atoms, fmax)

    # Returns True if optimised to or below desired fmax without constraints
    crystal = bulk("Cu")
    crystal.calc = EMT()
    if ase_env_check('3.23.0'):
        cell_relaxation = FrechetCellFilter(crystal)
    else:
        cell_relaxation = ExpCellFilter(crystal)
    opt = BFGSLineSearch(cell_relaxation)
    opt.run(fmax=fmax)
    # Explicitly remove the forces from the results - some calculators (e.g. MACE) only store stress for bulk
    crystal.calc.results.pop("forces")
    assert is_converged(crystal, fmax)

test_analyse_forces()
