def neighbour_cutout_sphere(atoms, centre, distance_cutoff=5.0):
    '''
    Returns a spherical cutout of a structure
    TODO: This doesn't work with periodic boundary conditions.
    TODO:Could probably integrate this with the neighbor function bellow

    Parameters:

    atoms: Atoms object
        Input structure to cutout from
    centre: Integer
        Index of central atom in cutout
    distance_cutoff: Float
        Distance which inside atoms are counted as neighbours
    '''

    import numpy as np

    # Prevents unexpected editing of parent object in place
    # Now ensures returned object is different to incoming atoms
    atoms = atoms.copy()
    selection = []
    for i in range(len(atoms)):

        # get distances between atom of interest and others


        #################### Edit atom tag ###############################
        distance_ab = np.linalg.norm((atoms.positions[centre] - atoms.positions[i]))

        ################ Edit distance in Angstrom here ###################
        if distance_ab < distance_cutoff:
            selection.append(i)


    return selection


def neighbours(atoms, centre, nearest_neighbors):
    ''' Returns a list of indices of atomic neighbors from a central atom
    TODO: Needs simplifying
    TODO: made into a function but now it's broken, will need to have a look

    Parameters:
    atoms : Atoms object
        Input structure to count neighbours
    centre : Integer
        Index of atom to start counting from
    nearest_neighbors : Integer
        Size of the neighbour shell, 1st neighbours, 2nd neighbours ....

    Returns:
        List of all atoms within specified neighbour distances
    '''

    # Object storing all neighbour information
    all_neighbours = [centre]
    # Object to store neighbours in "current" shell - start with 0th shell
    current_neighbors = [centre]

    # Creates an emtpy list an appends atom indices whose distances are
    # x amount nearest neighbours away from centre
    for neighbors in range(nearest_neighbors):
        new_neighbors = []
        for atom in current_neighbors:
            # TODO: Make the cutoff dynamic!
            neighbour_list = neighbour_cutout_sphere(atoms, centre=atom, distance_cutoff=3)
            for neighbour_to_atom in neighbour_list:
                if neighbour_to_atom not in all_neighbours:
                    # What is this?
                    # atoms.numbers[index] = neighbors
                    new_neighbors.append(neighbour_to_atom)
                    all_neighbours.append(neighbour_to_atom)

        # Copy the list of new neighbours to be our current "shell"
        current_neighbors = new_neighbors.copy()

    return all_neighbours
