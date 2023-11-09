from ase.io import read
from carmm.analyse.counterpoise_onepot import create_counterpoise_input
import subprocess


molecular_complex = read("data/NH3-H3O_traj/nh3-h3o.traj")
H2O_index = [0, 1, 2]
NH4_index = [atom.index for atom in molecular_complex if atom.index not in H2O_index]
dir_list = ['A_only', 'A_plus_ghost', 'B_only', 'B_plus_ghost']
create_counterpoise_input(complex_struc=molecular_complex, a_id=H2O_index, b_id=NH4_index, symbol_not_index=False)
third_H = "empty -0.0000000053828583 0.0220527604462460 1.6426686462733051 H\n"
file = open("A_plus_ghost/geometry.in", "r")
assert third_H in file.readlines()
file.close()
[subprocess.run(['rm', '-r', directory]) for directory in dir_list]

adsorption_complex = read("data/CO2_Cu/Cu_110_CO2c.traj")
CO2_ad = ['C', 'O']
Cu_slab = ['Cu']
create_counterpoise_input(complex_struc=adsorption_complex, a_id=CO2_ad, b_id=Cu_slab)
the_C = "atom 4.8402691096940629 3.3101933764080429 28.7941590443751956 C\n"
file = open("B_only/geometry.in", "r")
assert the_C not in file.readlines()
file.close()
[subprocess.run(['rm', '-r', directory]) for directory in dir_list]