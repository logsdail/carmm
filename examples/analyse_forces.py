#!/usr/bin/env python3

"""
Created on Fri 21/04/21

@author: Igor Kowalec

This script checks if forces data stored in the calculator attached to an Atoms
object meets the convergence criterion used in geometry optimisation.
"""

def test_analyse_forces():
    from ase.build import molecule
    from ase.calculators.emt import EMT
    from ase.optimize import BFGSLineSearch
    from carmm.analyse.forces import is_converged

    fmax = 0.05 # eV/Angstrom
    atoms = molecule("H2")
    atoms.calc = EMT()

    # Returns False prior to optimisation (no forces in calculator)
    opt = BFGSLineSearch(atoms)
    print(is_converged(atoms, fmax))
    assert(is_converged(atoms, fmax) == False)

    # Returns False if optimisation has not reached to/below fmax
    opt.run(fmax=fmax*30)
    print(is_converged(atoms, fmax))
    assert(is_converged(atoms, fmax) == False)

    # Returns True if optimised to or below desired fmax
    opt.run(fmax=fmax)
    print(is_converged(atoms, fmax))
    assert(is_converged(atoms, fmax) == True)

test_analyse_forces()
