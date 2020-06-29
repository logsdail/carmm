'''
Created on Fri 19/06/2020

@author: Igor Kowalec, David Willock
'''

def dissociation(atoms, i1, i2, step_size=0.05, n_steps=20, reverse=False, final_distance=None):

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
            Total number of steps
        reverse: boolean
            if True allows to examine association instead
        final_distance: float
            User can specify the final distance, the increments will be then based
            on a fraction of n_steps/final_distance instead of step_size
    '''

    from ase.constraints import FixBondLength
    import copy
    import numpy as np

    # retrieve initial atom - atom distance
    pos_diff = atoms[i1].position - atoms[i2].position
    initial_dist = np.linalg.norm(pos_diff)

    atoms_list = []
    distance_list = []

    # Retrieve previous constraints info
    if not atoms.constraints:
        initial_constraint = None
    else:
        initial_constraint = copy.deepcopy(atoms.constraints)

    for i in range(1, n_steps+1):
        # operate on a deepcopy for intended functionality
        atoms = copy.deepcopy(atoms)
        # remove previous contraints and set up new ones
        atoms.set_constraint()

        # move atoms and fix bond length in fixed increments or fraction of final_distance
        # For reverse, reduce the distance rather than increase (association)
        if reverse is False:
            if final_distance:
                measured_distance = initial_dist + (
                    i/n_steps * (final_distance - initial_dist))
            else:
                measured_distance = (initial_dist + i * step_size)
        elif reverse if True:
            if final_distance:
                measured_distance = (initial_dist - (
                    i/n_steps * (final_distance - initial_dist)))
            else:
                measured_distance = (initial_dist - i * step_size)

        atoms.set_distance(i1, i2, measured_distance, fix=0)

        # adjust contraints
        if initial_constraint is not None:
            new_constraint = initial_constraint + [FixBondLength(i1, i2)]
            atoms.set_constraint(new_constraint)
        else:
            new_constraint = FixBondLength(i1, i2)
            atoms.set_constraint(new_constraint)

        # Record the size of fixed bond
        atoms_list += [copy.deepcopy(atoms)]
        if reverse is False:
            distance_list += measured_distance
        elif reverse is True:
            distance_list += measured_distance

    return atoms_list, distance_list
