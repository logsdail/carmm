def cutout_sphere(atoms, centre, distance_cutoff=5.0):
    '''
    Returns a spherical cutout of a structure
    TODO: This doesn't work with periodic boundary conditions.
    It'd be nice to get that working and update the regression to test this.

    Parameters:

    atoms: Atoms object
        Input structure to cutout from
    centre: Integer
        Index of central atom in cutout
    distance_cutoff: Float
        Distance outside of which atoms are removed
    '''

    import numpy as np

    # Prevents unexpected editing of parent object in place
    # Now ensures returned object is different to incoming atoms
    atoms = atoms.copy()
    atoms_to_delete = []
    for i in range(len(atoms)):

        # get distances between atom of interest and others - then removes
        # all atoms beyond a certain radius

        #################### Edit atom tag ###############################
        distance_ab = np.linalg.norm((atoms.positions[centre] - atoms.positions[i]))

        ################ Edit distance in Angstrom here ###################
        if distance_ab > distance_cutoff:
            atoms_to_delete.append(i)

    del atoms[atoms_to_delete]

    return atoms

def transpose(periodic,cluster, start, stop, centre_periodic, centre_cluster, file_name)
    ''' Returns a ChemShell cluster representation of a periodic model 


    Parameters:
    periodic: Atoms object
             Input periodic model of interest
    cluster: Atoms object
             Input cluster model of interest 
    start: Integer 
             Starting index of atoms to copy into cluster
    stop:  Integer
              End index of atoms to copy into cluster
    centre_periodic: Integer
              Centre of periodic model
    centre_cluster: Integer 
              Centre of cluster model
    name: String 
              Name of new model to save 

    '''

    from ase.io import write

    # Start by specifying which atoms in periodic cluster you want to move to a cluster
    to_move = periodic[start:stop]
    #Define central atom in periodic for projection in cluster
    ref_old = periodic[centre_periodic].position
    for i in to_move:
        i.position -= ref_old
    #Define central atom in cluster for reference with periodic
    ref_new = cluster[centre_cluster].position
    for i in to_move:
        i.position += ref_new

    cluster = cluster + to_move
    scaled = str(cluster.get_scaled_positions())
    write = write(file_name, cluster)

