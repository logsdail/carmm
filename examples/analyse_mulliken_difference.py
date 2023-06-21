

'''
This function takes the aims.out file from after adsorption and before adsorption and 
returns the difference between the mulliken charge of the two. If you want to compare the bare slab with adsorbed
species, this code intelligently assigns zero to adsorbate and calculates the difference. Alternatively, you can also
provide adsorbate mulliken charges and it will account for it in the calculation. 

syntax: call this file followed by <post adsorption file>, <pre-adsorption file>, <adsorbate file (optional)>
'''


from carmm.analyse.mulliken import extract_mulliken_charge
from ase.io import read
import numpy as np
import sys
from ase.visualize import view


def set_initial_charges(final_file, initial_file, adsorbate_file=None):
    final_r = read(final_file)
    initial_r = read(initial_file)

    mulliken_charge_1 = extract_mulliken_charge(final_file, len(final_r))
    mulliken_charge_1 = [float(x) for x in mulliken_charge_1]
    mulliken_charge_2 = extract_mulliken_charge(initial_file, len(initial_r))
    mulliken_charge_2 = [float(x) for x in mulliken_charge_2]

    if adsorbate_file is not None:
        adsorbate_r = read(adsorbate_file)
        adsorbate_length = len(adsorbate_r)

        adsorbate_mulliken_charges = extract_mulliken_charge(adsorbate_file, adsorbate_length)
        adsorbate_mulliken_charges = [float(x) for x in adsorbate_mulliken_charges]

        if len(final_r) == len(initial_r) + len(adsorbate_r):
            mulliken_charge_2.extend(adsorbate_mulliken_charges)
        else:
            print("Mismatch in the total number of atoms. Check all geometries for consistency.")
            mulliken_charge_2.extend([0] * adsorbate_length)

    else:
        if len(final_r) != len(initial_r):
            length_diff = abs(len(final_r) - len(initial_r))
            mulliken_charge_2.extend([0] * length_diff)

    mulliken_charge = np.subtract(mulliken_charge_1, mulliken_charge_2)

    final_r.set_initial_charges(mulliken_charge)
    view(final_r)    

args = sys.argv[1:]
final_file = args[0]
initial_file = args[1]
adsorbate_file = args[2] if len(args) >= 3 else None

# Call the function with command-line arguments
#set_initial_charges(final_file, initial_file, adsorbate_file)   #uncomment this line for use. 
set_initial_charges ('data/mulliken/aims.out', 'data/mulliken/aims_1.out', 'data/mulliken/pt.in')