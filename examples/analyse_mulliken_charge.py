#!/usr/bin/env python3

'''
This script presents an example and QA test for analysing the Mulliken charge data
in an FHI-aims output, by imposing this charge data on an Atoms object and then
allowing you to visualise the data in the ASE GUI.

This is useful when trying to understand the electronic structure of your system
'''

def test_analyse_mulliken_charge():

    from carmm.analyse.mulliken import extract_mulliken_charge

    # The filename needed is the aims.out file not the fhiaims.hawk.output.log.XXXXXXX as said previously

    filename = "data/CO/co_light.log"

    from ase.io import read
    # Read atoms in
    atoms = read(filename)

    # Read in the Mulliken data
    mulliken_charge = extract_mulliken_charge(filename, len(atoms))
    # If you wanted to compare charges, you could subtract one value for Mulliken_charge from another to give a difference
    # and then set this as a charge on your model e.g. when comparing surfaces with and without adsorbates.
    atoms.set_initial_charges(mulliken_charge)

    # This opens the visualiser to see colours
    # To change colours, go to: View -> Colurs -> By Initial Charge (or just press "C")
    # The colour map options are available here: https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html
    # The classical colour map would be bwr = Blue - White - Red; to get this working you must
    # reverse the min/max aspects of the colour scale, as the default (unusual) colours are -ve = blue, +ve = red.
    ### ENABLE TO VIEW
    from ase.visualize import view
    view(atoms)

    # Confirm we are reading the Mulliken Charge correctly
    assert(mulliken_charge == ['0.088700', '-0.088700'])

# Run the example
test_analyse_mulliken_charge()