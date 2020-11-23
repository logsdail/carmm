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

    from ase import Atoms
    from ase.build import add_adsorbate

    # Optimisation is disabled to ensure quick testing
    # If you are running this to see the functionality, enable the optimisation
    optimisation = False
    if optimisation:
        from ase.calculators.emt import EMT
        from ase.optimize import BFGS

    # Generate and optimise slab prior to dissociation
    # @ikowalec: Are you sure this is correct? The models when plotted seem a bit ... crazy ...
    slab = get_example_slab()
    adatoms = Atoms("2Cu", positions=[(2, 0, 0), (0, 0, 0)])
    add_adsorbate(slab, adatoms, height=2, position=(slab[0].x, slab[0].y))

    if optimisation:
        opt = BFGS(slab)
        opt.run(fmax=0.05)

    # move adsorbed Cu atoms (indices 18,19) away from each other
    atoms_list, distance_list = dissociation(slab, 18, 19, step_size=0.2, n_steps=10, z_bias=True, group_move=[19])

    # Assertion tests - checking no one has broken the code.
    assert (len(distance_list) == len(atoms_list) == 10)
    if optimisation:
        assert(distance_list[0] - 3.145244 < 1e-5)
        assert(distance_list[9] - 6.063336 < 1e-5)
    else:
        assert(distance_list[0] - 2.685424 < 1e-5)
        assert(distance_list[9] - 5.539200 < 1e-5)

    if optimisation:
        for atoms in atoms_list:
            atoms.calc = EMT()
            opt = BFGS(atoms)
            opt.run(fmax=0.05)

        from ase.visualize import view
        view(atoms_list)


test_dissociation()
