def cutout_sphere(atoms, centre, distance_cutoff=5.0):
    '''
    Returns a spherical cutout of a structure
    TODO: This doesn't work with periodic boundary conditions.
    It'd be nice to get that working and update the regression to test this.

    Parameters:

    atoms: Atoms object or String
        Input structure to cutout from, or the file containing this information
    centre: Integer
        Index of central atom in cutout
    distance_cutoff: Float
        Distance outside of which atoms are removed
    '''

    import numpy as np

    if isinstance(atoms, str):
        from ase.io import read
        atoms = read(atoms)

    atoms_to_delete = []
    for i in range(len(atoms)):

        #get distances between atom of interest and others - then constrain
        #all atoms beyond a certain radius

        #################### Edit atom tag ###############################
        distance_ab = np.linalg.norm((atoms.positions[centre] - atoms.positions[i]))

        ################ Edit distance in Angstrom here ###################
        if distance_ab > distance_cutoff:
        ###################################################################
            atoms_to_delete.append(i)

    del atoms[atoms_to_delete]

    return atoms