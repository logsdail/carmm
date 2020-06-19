'''
Created on Fri 19/06/2020

@author: Igor Kowalec, David Willock
'''

def dissociation(atoms, i1, i2, step_size=0.05, n_steps=20):

    '''
    This function is a tool for investigating bond dissociation.
    Bond length of interest is fixed and is increased by step_size in each iteration.
    This aims to help with obtaining activation energies for surface calculations, e.g. hydrogenation,
    where metastability of optimal starting position is often low and thus hard to obtain.
    Returns a list of Atoms objects with changed positions and constraints applied,
    which can be optimised by the user.

    Args:
        atoms: Atoms object
        i1: int
            Index of atom remaining as part of a molecule
        i2: int
            Index of atom dissociating from molecule
        step_size: float
            Distance moved during dissociation in Angstrom per iteration
        n_steps: int
            Total number of steps, including initial
    '''

    from ase.constraints import FixBondLength
    import copy
    import numpy as np

    # retrieve initial atom - atom distance
    pos_diff = atoms[i1].position - atoms[i2].position
    initial_dist = np.linalg.norm(pos_diff)

    # operate on deepcopy to avoid changes to original object in memory
    atoms = copy.deepcopy(atoms)
    atoms_list = []
    distance_list = []

    # Retrieve previous constraints info
    if not atoms.constraints:
        initial_constraint = None
    else:
        initial_constraint = copy.deepcopy(atoms.constraints)

    for i in range(0, n_steps):
        # operate on a deepcopy for intended functionality
        atoms = copy.deepcopy(atoms)
        #remove previous contraints and set up new ones
        atoms.set_constraint()
        atoms.set_distance(i1, i2, (initial_dist + i * step_size), fix=0)

        new_constraint = FixBondLength(i1, i2)

        if initial_constraint is not None:
            atoms.set_constraint([initial_constraint, new_constraint])
        else:
            atoms.set_constraint(new_constraint)

        atoms_list += [copy.deepcopy(atoms)]
        distance_list += [initial_dist + i * step_size]

    return atoms_list, distance_list