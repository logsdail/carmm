from ase.io import read
from ase.neighborlist import neighbor_list
from ase import neighborlist
import numpy as np
from scipy import sparse
from ase.visualize import view

atoms = read("HNBBC.traj", )
cutOff = neighborlist.natural_cutoffs(atoms)
neighborList = neighborlist.NeighborList(cutOff, self_interaction=False, bothways=True)
neighborList.update(atoms)
matrix = neighborList.get_connectivity_matrix()
n_components, component_list = sparse.csgraph.connected_components(matrix)

molecules = []
for n in range(n_components):
    atomsIdxs = [i for i in range(len(component_list)) if component_list[i] == n]
    print("The following atoms are part of molecule {}: {}".format(n, atomsIdxs))
    molecules.append(atomsIdxs)

print(molecules[0])
print(type(molecules[0]))


#idx = 0
#while idx <= 251:
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
