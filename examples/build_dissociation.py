#!/usr/bin/env python3

'''
Created on Fri 19/06/2020

@author: Igor Kowalec

Example use of a tool for investigating bond dissociation.
Bond length of interest is fixed and is increased by step_size in each iteration.
Returns a list of Atoms objects with changed positions and constraints applied,
which can be optimised by the user - toy model of 2 Cu ad-atoms on Au surface.

'''

def test_dissociation():
    from carmm.build.neb.bond_length_scan import dissociation
    from examples.data.model_gen import get_example_slab
    
    from ase import Atoms
    from ase.build import add_adsorbate
    from ase.calculators.emt import EMT
    from ase.optimize import BFGS

    # Generate and optimise slab prior to dissociation
    slab = get_example_slab()
    adatoms = Atoms("2Cu", positions=[(0,0,0), (2,0,0)])
    add_adsorbate(slab, adatoms, height=2, position=(slab[0].x, slab[0].y))
    opt = BFGS(slab)
    opt.run(fmax=0.05)

    # move adsorbed Cu atoms (indices 18,19) away from each other
    atoms_list, distance_list = dissociation(slab, 18, 19, step_size=0.2, n_steps=10)

    # Do we need the optimisation for all species as part of the test? 
    # No idea how long this takes - could be super quick!
    for atoms in atoms_list:
        atoms.set_calculator(EMT())
        opt = BFGS(atoms)
        opt.run(fmax=0.05)

    #TODO: Add in an assertion test so that functionality can actually be verified.
    
    #from ase.visualize import view
    #view(atoms_list)
    #print(distance_list)

test_dissociation()
