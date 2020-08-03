#!/usr/bin/env python3

'''
This script presents an example and QA test for analysing the Mulliken charge data
in an FHI-aims output, by imposing this charge data on an Atoms object and then
allowing you to visualise the data in the ASE GUI.

This is useful when trying to understand the electronic structure of your system
'''

def test_analyse_mulliken_spin():

    from carmm.analyse.mulliken import extract_mulliken_spin

    # This is the filename that would need changing for each FHI-aims input e.g. fhiaims.hawk.output.log.XXXXXXX
    filename = "data/Fe-CO/fe-co_light.log"

    from ase.io import read
    # Read atoms in
    atoms = read(filename)

    # Read in the Mulliken data
    mulliken_spin = extract_mulliken_spin(filename, len(atoms))
    # If you wanted to compare charges, you could subtract one value for Mulliken_charge from another to give a difference
    # and then set this as a charge on your model e.g. when comparing surfaces with and without adsorbates.
    atoms.set_initial_charges(mulliken_spin)

    # This opens the visualiser to see colours
    # To change colours, go to: View -> Colurs -> By Initial Charge (or just press "C")
    # The colour map options are available here: https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html
    # The classical colour map would be bwr = Blue - White - Red; to get this working you must
    # reverse the min/max aspects of the colour scale, as the default (unusual) colours are -ve = blue, +ve = red.
    ### ENABLE TO VIEW
    #from ase.visualize import view
    #view(atoms)

    # Assertion test to confirm the code works correctly.
    assert(mulliken_spin == ['-2.242867', '0.164614', '0.078252'])

# Run the example
test_analyse_mulliken_spin()