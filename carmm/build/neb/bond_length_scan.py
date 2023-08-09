'''
Created on Fri 19/06/2020

@author: Igor Kowalec, David Willock
'''

def dissociation(atoms, i1, i2, step_size=0.05, n_steps=20, final_distance=None, group_move=None, z_bias=False):

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
            Distance moved during dissociation in Angstrom per iteration, not used
            when final_distance is specified.
            If negative value - Association is examined instead
        n_steps: int
            Total number of steps
        final_distance: None/float
            User can specify the final distance, the increments will be then based
            on a fraction of n_steps/final_distance instead of step_size
        group_move: list of integers
            User can specify a list of indices of atoms that need to be moved
            together with atom with index i2, e.g. OH group etc.
        z_bias: boolean of float
            If float - bias z-coord of moving atom to approach set value
            WARNING - If TRUE this will make steps vary from defined step_size!
            Bias to adjust the Z-coordinate of atom/group moving to approach
            the surface in periodic calculations rather than just elongate the
            bond.
    '''

    from ase.constraints import FixBondLength
    import copy
    import numpy as np
    from ase.neighborlist import natural_cutoffs

    # retrieve natural cutoffs for atomic radii
    cutOff = natural_cutoffs(atoms, mult=2)

    # retrieve initial atom - atom distance and atom[i2] position
    pos_diff = atoms[i1].position - atoms[i2].position
    initial_dist = np.linalg.norm(pos_diff)
    initial_i2_pos = copy.deepcopy(atoms[i2].position)

    # retrieve z-coordinate for z_bias
    if z_bias:
        # Define a list of atom chemical symbols and their count
        # Take the z coordinate of the most abundant element in the surface slab
        chem_symbol_count = [[x, atoms.get_chemical_symbols().count(x)] for x in set(atoms.get_chemical_symbols())]
        # Sort the list by the count
        chem_symbol_count.sort(key=lambda k: k[1])
        surf_z_list = [atom.z for atom in atoms if atom.symbol == chem_symbol_count[-1][0] and atom.tag > 0]
        # The maximum z-coordinate of the surface atoms is retrieved
        surf_z = np.amax(surf_z_list)
        surf_atom_top = [atom for atom in atoms if atom.z == surf_z][0]

        if isinstance(z_bias, (int, float)) and not isinstance(z_bias, bool):
            z_threshold_min = z_bias
        else:
            z_threshold_min = surf_z + cutOff[surf_atom_top.index]

        assert z_threshold_min >= surf_z + cutOff[surf_atom_top.index], \
            f"Your z-coordinate bias is too close to the surface, Surface z: {surf_z}, with a radius of: " \
            f"{cutOff[surf_atom_top.index]} z-bias: {z_bias}."


    atoms_list = []
    distance_list = []

    # Retrieve previous constraints info
    if not atoms.constraints:
        initial_constraint = None
    else:
        initial_constraint = copy.deepcopy(atoms.constraints)

    for i in range(1, n_steps+1):
        # operate on a deepcopy for intended functionality
        atoms = atoms.copy()
        # remove previous constraints and set up new ones
        atoms.set_constraint()

        # initial moving atom position
        imap = atoms.copy()[i2].position

        # move atoms and fix bond length in fixed increments or fraction of final_distance
        measured_distance = (initial_dist + i * step_size)

        if final_distance:
            measured_distance = initial_dist + (
                i/n_steps * (final_distance - initial_dist))

        atoms.set_distance(i1, i2, measured_distance, fix=0)

        # apply bias in z-coordinate towards the surface atoms
        if z_bias:
            # min distance in Angstrom from surf atoms
            if group_move:
                pass
            else:
                z_threshold_min = surf_z
            # make sure atoms from a group do not clash into surface atoms
            # move towards the surface or away if necessary
            if group_move:
                z_threshold_max = np.amin([atom.z for atom in atoms[group_move]])
            else:
                z_threshold_max = atoms[i2].z


            if z_threshold_max > z_threshold_min:
                atoms[i2].z -= (initial_i2_pos[2] - surf_z)/n_steps
            elif z_threshold_max < z_threshold_min:
                atoms[i2].z += (initial_i2_pos[2] - surf_z)/n_steps

        # move other specified atoms as part of a molecule, e.g. OH group
        if group_move:
            for m in group_move:
                if not m == i2:
                    atoms[m].position = atoms[m].position + (atoms[i2].position - imap)

        # adjust contraints
        if initial_constraint is not None:
            new_constraint = initial_constraint + [FixBondLength(i1, i2)]
            atoms.set_constraint(new_constraint)
        else:
            new_constraint = FixBondLength(i1, i2)
            atoms.set_constraint(new_constraint)

        # Record the size of fixed bond
        atoms_list += [atoms.copy()]
        distance_list += [np.linalg.norm(atoms[i1].position - atoms[i2].position)]

    return atoms_list, distance_list
