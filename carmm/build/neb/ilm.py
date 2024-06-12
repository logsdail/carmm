def multiple_local_extrema(filename="last_predicted_path.traj"):
    '''
    This function will detect local maxima in the trajectory file of a NEB
    calculation.

    Parameters:
    filename: str
        Name of the file containing latest Minimum Energy Path, default
        is "last_predicted_path.traj" as in MLNEB by CatLearn.
    Return:
        False for single maximum; True for multiple maxima.
    '''

    import numpy as np
    from scipy.signal import argrelextrema
    from ase.io.trajectory import Trajectory

    # for local maxima
    #argrelextrema(x, np.greater)

    # for local minima
    #argrelextrema(x, np.less)

    # Read all structures
    traj = Trajectory(filename, 'r')

    # Define variable
    list_to_1d_array = []
    mep = []
    minima = []

    # Retrieve energies and list them
    for atoms in traj:
        e = atoms.get_potential_energy()
        list_to_1d_array.append(e)
        mep += [atoms]

    # Find and identify local maxima in array of energies
    list_to_1d_array = np.array(list_to_1d_array)
    indices = argrelextrema(list_to_1d_array, np.greater)
    indices = [x for x in indices[0]]

    if len(indices) > 1:
        return True
    else:
        return False




def neb_identify_local_minima(filename="last_predicted_path.traj"):

    '''
    This function will detect local minima in the trajectory file of a NEB
    calculation and will output the associated structures and their
    index. The user can then optimise these and use as an alternative starting
    structures for the NEB calculations.

    Parameters:
    filename: str
        Name of the file containing latest Minimum Energy Path, default
        is "last_predicted_path.traj" as in MLNEB by CatLearn.
    '''

    import numpy as np
    from scipy.signal import argrelextrema
    from ase.io.trajectory import Trajectory

    # for local maxima
    #argrelextrema(x, np.greater)

    # for local minima
    #argrelextrema(x, np.less)

    # Read all structures
    traj = Trajectory(filename, 'r')

    # Define variable
    list_to_1d_array = []
    mep = []
    minima = []

    # Retrieve energies and list them
    for atoms in traj:
        e = atoms.get_potential_energy()
        list_to_1d_array.append(e)
        mep += [atoms]

    # Find and identify local minima in array of energies
    list_to_1d_array = np.array(list_to_1d_array)
    indices = argrelextrema(list_to_1d_array, np.less)
    indices = [x for x in indices[0]]

    # Images closest to input structure to avoid unnecessary
    # optimisation resulting in the same structure as input.
    # Remove minima of little significance.
    for i in indices:
        if i < int(len(mep)/10):
            indices = indices - [i]
        elif i > len(mep) - int(len(mep)/10):
            indices = indices - [i]


    for i in indices:
        minima += [mep[i]]

    return minima, indices
