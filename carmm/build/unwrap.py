import numpy as np
from scipy.sparse.csgraph import minimum_spanning_tree, breadth_first_order

def unwrap(mol, anchor_atom_index = 0):
    ''' Unwraps molecule in periodic cell. Makes minimum spanning tree of the molecule space continuous.

    Parameters:
    mol: Atoms object
            Input periodic model of interest.
    anchor_atom_index: index of the Atom that will be used as an anchor for the rest of molecule. 
            Position of this atom will not be changed. 

    Returns: Unwrapped molecule
    '''
    assert 0 <= anchor_atom_index < len(mol)
    V = mol.get_all_distances(True, True)
    D = np.linalg.norm(V, axis=-1)
    st = minimum_spanning_tree(D).toarray()
    order, preds = breadth_first_order(st, anchor_atom_index, directed=False, return_predecessors=True)
    for i in order[1:]:
        mol.positions[i] = mol.positions[preds[i]] + V[preds[i], i, :]
    return mol

