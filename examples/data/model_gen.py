#!/usr/bin/env python3

def get_example_slab(adsorbate=False, type="CO2"):
    '''
    Example model generation to show tool functionality
    Parameters:

    adsorbate: Boolean
        Whether or not to include an adsorbate CO2 on the surface
    type: str
        Two adsorbate types - "CO2" and "2Cu". CO2 behaviour during optimisation
        is unrealistic in EMT.
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

    # Put adsorbate on the surface
    if adsorbate:
        if type == "CO2":
            position = (slab[17].position[0], slab[17].position[1])
        elif type == "2Cu":
            position = (slab[5].x, slab[5].y)

        species = get_example_adsorbate(type)
        add_adsorbate(slab, species, 3.0, position=position)

    # Make the model a bit more technically complete - include a calculator.
    from ase.calculators.emt import EMT
    slab.set_calculator(EMT())

    return slab

def get_example_adsorbate(type="CO2"):
    '''
    Small function to return CO2 or Cu atoms
    Parameters:

    type: str
        Two adsorbate types - "CO2" and "2Cu". CO2 behaviour during optimisation
        is unrealistic in EMT.
    '''

    from ase.build import molecule
    from ase import Atoms
    if type == "CO2":
        atoms = molecule("CO2")
        atoms.rotate(90, 'x')
        atoms.rotate(50, 'z')
    elif type == "2Cu":
        atoms = Atoms("2Cu", positions=[(0,0,0), (-4.08/(2**(1/2)),0,0)])


    # Make the model a bit more technically complete - include a calculator.
    from ase.calculators.emt import EMT
    atoms.set_calculator(EMT())

    return atoms
