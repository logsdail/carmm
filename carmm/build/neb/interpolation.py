def check_interpolation(initial, final, n_max, interpolation="linear", verbose=True, save=True):
    '''
    Interpolates the provided geometries with n_max total images
    and checks whether any bond lengths are below sane defaults.
    Saves the interpolation in interpolation.traj

    Parameters:

    initial: Atoms object or string
        Starting geometry for interpolation.
    final: Atoms object or string
        End point geometry for interpolation
    n_max: integer
        Desired total number of images for the interpolation
        including start and end point.
    interpolation: string
        "linear" or "idpp". First better for error identification, latter for
        use in NEB calculation
    verbose: boolean
        If verbose output of information is required
    save: boolean
        Whether to save the trajectory for transfer on to an NEB calculation
    '''

    from ase.neb import NEB
    from carmm.analyse.bonds import search_abnormal_bonds
    from ase.io.trajectory import Trajectory
    from ase.io import read

    # Pre-requirements
    if not isinstance(n_max, int):
        raise ValueError
        print("Max number of images must be an integer.")

    # Make a band consisting of 10 images:
    images = [initial]
    images += [initial.copy() for i in range(n_max-2)]
    images += [final]
    neb = NEB(images)
    # Interpolate linearly the potisions of the middle images:
    neb.interpolate(interpolation, apply_constraint=True)

    #TODO: Tidy up this horrible mix of if statements.
    if save:
        t = Trajectory('interpolation.traj', 'w')

    flag = True
    for i in range(0, n_max):
        if verbose:
            print("Assessing image", str(i+1) + '.')
        updated_flag = search_abnormal_bonds(images[i], verbose)
        if save:
            t.write(images[i])
        if (not updated_flag):
            flag = updated_flag

    if save:
        t.close()

    return flag
