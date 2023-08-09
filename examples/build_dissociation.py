#!/usr/bin/env python3

"""
Created on Fri 19/06/2020

@author: Igor Kowalec

Example use of a tool for investigating bond dissociation.
Bond length of interest is fixed and is increased by step_size in each iteration.

Returns a list of Atoms objects with changed positions and constraints applied,
which can be optimised by the user - toy model of 2 Cu ad-atoms on Au surface.

Note: Optimisation is disabled by default in test. Enable when working with this
personally so you can see how everything would work in a practical application

"""


def test_dissociation():
    from carmm.build.neb.bond_length_scan import dissociation
    from examples.data.model_gen import get_example_slab
    from ase.constraints import FixAtoms

    # Generate and optimise slab prior to dissociation
    slab = get_example_slab(adsorbate=True, type="2Cu")
    # move one of the Cu atoms down, so that the bond is not parallel to the surface
    slab[-2].z -= 0.7

    # move adsorbed Cu atoms (indices 18,19) away from each other
    # if z_bias is True, the distance increment is not exactly step_size value
    z_bias = True
    atoms_list, distance_list = dissociation(slab, 18, 19, step_size=0.2, n_steps=10, z_bias=z_bias, group_move=[19])

    # Assertion tests - checking no one has broken the code.
    assert len(distance_list) == len(atoms_list) == 10
    assert(distance_list[0] - 3.111653 < 1e-5)
    assert(distance_list[9] - 4.939653 < 1e-5)

    # if z_bias is a float, the atom or group of atoms will be biased to stay around the z_bias value in z-coordinate
    z_bias = 16
    atoms_list, distance_list = dissociation(slab, 18, 19, step_size=0.2, n_steps=10, z_bias=z_bias, group_move=[19])

    # Assertion tests - checking no one has broken the code.
    assert len(distance_list) == len(atoms_list) == 10
    assert (distance_list[0] - 3.2525338 < 1e-5)
    assert (distance_list[9] - 4.881640 < 1e-5)



    # Set constraints and remove z_bias to test functionality
    slab.set_constraint([FixAtoms([atom.index for atom in slab if atom.tag > 1])])
    z_bias = False
    # group_move set to False since only one atom is moving
    atoms_list, distance_list = dissociation(slab, 18, 19, step_size=0.2, n_steps=10, z_bias=z_bias, group_move=False)
    assert(distance_list[0] - 3.168703 < 1e-5)
    assert(distance_list[9] - 4.968703 < 1e-5)


    """
    # Optimisation is disabled to ensure quick testing
    # If you are running this to see the functionality, enable this code block
    from ase.calculators.emt import EMT
    from ase.optimize import BFGS
    opt = BFGS(slab)
    opt.run(fmax=0.05)
    for atoms in atoms_list:
        atoms.calc = EMT()
        opt = BFGS(atoms)
        opt.run(fmax=0.05)

    from ase.visualize import view
    view(atoms_list)
    """

test_dissociation()
