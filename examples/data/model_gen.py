#!/usr/bin/env python3

def get_example_slab(adsorbate=False):
    '''
    Example model generation to show tool functionality
    Parameters:

    adsorbate: Boolean
        Whether or not to include an adsorbate CO2 on the surface
    '''
    from math import sqrt
    from ase.build import molecule, add_adsorbate, fcc111

    #### Traditional ASE functionality #####
    element='Au'
    lattice_parameter=2.939
    width=3
    depth=2
    vacuum=10.0

    # Create surface
    slab = fcc111(element, a=lattice_parameter*sqrt(2),
                  size=(width, width, depth), vacuum=vacuum)

    # Put CO2 on the surface
    if adsorbate:
        CO2 = get_example_adsorbate()
        CO2.rotate(90, 'x')
        CO2.rotate(50, 'z')

        add_adsorbate(slab, CO2, 3.0, position=(
            slab[17].position[0], slab[17].position[1]))

    return slab

def get_example_adsorbate():
    '''
    Small function to return CO2 molecule for testing
    '''

    from ase.build import molecule

    return molecule("CO2")
