#!/usr/bin/env python3

def get_example_slab_with_adsorbate():
    '''Example model generation to show tool functionality'''
    from math import sqrt
    from ase.build import molecule, add_adsorbate, fcc111
    from ase.optimize import BFGS
    from ase.calculators.emt import EMT

    #### Traditional ASE functionality #####

    element='Au'
    lattice_parameter=2.939
    width=3
    depth=2
    vacuum=10.0

    # Create surface
    slab = fcc111(element, a=lattice_parameter*sqrt(2),
                  size=(width, width, depth), vacuum=vacuum)


    CO2 = molecule("CO2")
    CO2.rotate(90, 'x')
    CO2.rotate(50, 'z')

    add_adsorbate(slab, CO2, 3.0, position=(
        slab[17].position[0], slab[17].position[1]))

    return slab
