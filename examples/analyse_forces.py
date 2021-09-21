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
    from ase.constraints import FixAtoms

    fmax = 0.01 # eV/Angstrom
    atoms = molecule("CO2")
    atoms.calc = EMT()

    # Returns False prior to optimisation (no forces in calculator)
    opt = BFGSLineSearch(atoms)
    print(is_converged(atoms, fmax))
    assert(is_converged(atoms, fmax) == False)
    print(atoms.calc.forces)

    # Returns False if optimisation has not reached to/below fmax
    opt.run(fmax=fmax*30)
    print(is_converged(atoms, fmax))
    assert(is_converged(atoms, fmax) == False)

    # Returns True if optimised to or below desired fmax with constraints
    c = FixAtoms(indices=[0,1])
    atoms.set_constraint(c)
    opt.run(fmax=fmax)
    print(is_converged(atoms, fmax))
    assert (is_converged(atoms, fmax) == True)

    # Returns True if optimised to or below desired fmax without constraints
    atoms.set_constraint()
    opt.run(fmax=fmax)
    print(is_converged(atoms, fmax))
    assert(is_converged(atoms, fmax) == True)


test_analyse_forces()
