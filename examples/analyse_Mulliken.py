#!/usr/bin/env python3

def write_example_output_for_testing(fn):
    '''
    This is just to write an example output for testing ...
    '''
    with open(fn, 'a') as f:
        f.write('  ------------------------------------------------------------\n')
        f.write('                   Starting Mulliken Analysis\n')
        f.write('  ------------------------------------------------------------\n')
        f.write('                \n')
        f.write('  Performing Mulliken charge analysis on all atoms.\n')
        f.write(" Full analysis will be written to separate file 'Mulliken.out'.\n")
        f.write(' Summary of the per-atom charge analysis:\n')
        f.write('  |\n')
        f.write('  |  atom       electrons          charge             l=0             l=1             l=2\n')
        f.write('  |     1        5.911300        0.088700        3.781721        2.024840        0.104739\n')
        f.write('  |     2        8.088700       -0.088700        3.784537        4.266947        0.037216\n')

def delete_example_output(fn):
    '''
    ... and this deletes the example input to tidy up at the end!
    '''
    import os
    if os.path.exists(fn):
        os.remove(fn)

from software.analyse.mulliken import extract_mulliken_charge

# This is the filename that would need changing for each FHI-aims input e.g. fhiaims.hawk.output.log.XXXXXXX
filename = "test"
# This is just to write a temporary file to QA test the Mulliken extraction. Delete in personal scripts
if "test" in filename:
    from ase.build import molecule
    # Create test data
    write_example_output_for_testing(filename)
    atoms = molecule('CO')
else:
    from ase.io import read
    # Read atoms, more like real situation
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
#from ase.visualize import view
#view(atoms)

if "test" in filename:
    delete_example_output(filename)

#TODO: Add in assertion so this is a rigorous QA test.