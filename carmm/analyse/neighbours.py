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
        centre : Index of atom to start counting from

        nearest_neighbors : Size of the neihgbour shell, 1st neighbours, 2nd neighbours ....
     '''
    all_neighbours = [0]
    #Choose atom index to start counting neighbours from
    current_neighbors = [centre]
    # Creates an emtpy list an appends atom indices whose distances are x amount nearest neighbours away from centre
    for neighbors in range(nearest_neighbors):
        new_neighbors = []
        for i in current_neighbors:
            list = neighbour_cutout_sphere(atoms, centre=i, distance_cutoff=2.5)
            for index in list:
                if index not in all_neighbours:
                    atoms.numbers[index] = neighbors
                    new_neighbors.append(index)
                    all_neighbours.append(index)
        current_neighbors = new_neighbors.copy()

