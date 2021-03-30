def Calculate_planes(atoms, print_output=False):

from typing import List
from ase.io import read
from ase.neighborlist import neighbor_list
from ase import neighborlist
import numpy as np
from scipy import sparse
from ase.visualize import view
from ase.geometry.analysis import Analysis

atoms = read("aims.out")
analysis = Analysis(atoms)
cutOff = neighborlist.natural_cutoffs(atoms)
neighborList = neighborlist.NeighborList(cutOff, self_interaction=False, bothways=True)
neighborList.update(atoms)
matrix = neighborList.get_connectivity_matrix()
n_components, component_list = sparse.csgraph.connected_components(matrix)

molecules = []
for n in range(n_components):
    atomsIdxs: List[int] = [i for i in range(len(component_list)) if component_list[i] == n]
    print("The following atoms are part of molecule {}: {}".format(n, atomsIdxs))
    molecules.append(atomsIdxs)

distances = atoms.get_all_distances(mic=True, vector=False)

A = atomsIdxs[
   [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 102, 104, 106, 108, 110, 112, 114, 116, 118, 138, 140, 142, 144, 146, 148, 150, 152, 154, 156, 158, 160, 162, 164, 166, 168, 170, 172, 174, 176, 178, 180, 182, 184, 186, 188, 190, 192, 194, 196, 198, 200, 202, 204, 206, 208, 210, 212, 214, 216, 218, 220, 237, 239, 241, 243, 245, 247],

B = atomsIdxs[
   [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97, 99, 101, 103, 105, 107, 109, 111, 113, 115, 117, 119, 139, 141, 143, 145, 147, 149, 151, 153, 155, 157, 159, 161, 163, 165, 167, 169, 171, 173, 175, 177, 179, 181, 183, 185, 187, 189, 191, 193, 195, 197, 199, 201, 203, 205, 207, 209, 211, 213, 215, 217, 219, 221, 238, 240, 242, 244, 246, 248],

print_plane = A + "-" + B

AB_Bonds = analysis.get_bonds(A, B)
if AB_Bonds == [[]]:
   AB_BondsValues = None
else:
   AB_BondsValues = analysis.get_values(AB_Bonds)

if verbose and AB_BondsValues is not None:
   if not multirow:
       print_bond_table_header()
       # Table contents
#    import numpy as np
 #   print('{:<8.8s}{:<6.0f}{:>4.6f}{:^12.6f}{:>4.6f}'.format(
 #       print_AB, len(AB_BondsValues[0]), np.average(AB_BondsValues),
 #       np.amin(AB_BondsValues), np.amax(AB_BondsValues)))

#return print_AB, AB_Bonds, AB_BondsValues

print(atoms.symbols.get_chemical_formula())
print(molecules[0])
print(type(molecules[0]))
print(molecules.symbols.get_chemical_formula('hill', 'empirical'))

# idx = 0
# while idx <= 251:
#    atomsIdx = component_list[idx]
#    print("There are {} molecules in the system".format(n_components))
#    print("Atom {} is part of molecule {}".format(idx, atomsIdx))
#    atomsIdxs = [i for i in range(len(component_list)) if component_list[i] == atomsIdx]
#   print("The following atoms are part of molecule {}: {}".format(atomsIdx, atomsIdxs))
#    idx = idx + 1


view(atoms)
# for atom_index in range(0, len(atoms)):
#    assert isinstance(offsets, object)
#    print(atom_index, indices, offsets)
