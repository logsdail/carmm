#def Calculate_planes(atoms, print_output=False);

from typing import List
from ase.io import read
from ase.neighborlist import neighbor_list
from ase import neighborlist
import numpy as np
from scipy import sparse
from ase.visualize import view
from ase.geometry.analysis import Analysis
#from carmm.analyse.bonds import



atoms = read("aims.out")
analysis = Analysis(atoms)

from carmm.analyse.molecules import calculate_molecules


molecules = calculate_molecules(atoms)

distances = atoms.get_all_distances(mic=True, vector=False)
A = molecules[0]
B = molecules[1]


A_mol = atoms[A]
B_mol = atoms[B]

def get_lowest_distances(A_mol, B_mol):
    '''Definition'''
    measured = []
    for a in A_mol :
        bond_list_atom = []
        for b in B_mol:
            pos_diff = np.linalg.norm(a.position - b.position)
            bond_list_atom+= [pos_diff]
        measured+= [bond_list_atom]
    lowest_distances = [np.amin(i) for i in measured]

    return lowest_distances


print(get_lowest_distances(A_mol, B_mol))

#return print_AB, AB_Bonds, AB_BondsValues

#print(atoms.symbols.get_chemical_formula())
#print(molecules[0])
#print(type(molecules[0]))
#print(molecules.symbols.get_chemical_formula('hill', 'empirical'))

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
