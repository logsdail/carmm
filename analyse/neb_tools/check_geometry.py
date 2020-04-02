def switch_indices(model, A, B):
    '''
    Function for rearranging atomic indices in the structure.
    Returns a changed atoms object with retained calculator
    and ensures array of forces is rearranged accordingly.

    Parameters:

    model: Atoms object or string
        Structure requiring atom index rearrangement.
        If a string, e.g. 'name.traj', a file of this name will
        be read.
    A, B: integer
        Indices of atoms that need to be switched

    '''

    from ase import Atoms
    from ase.io import read

    if isinstance(model, str) is True:
        model = read(model)

    # Index manipulation
    if not isinstance(A, int) and isinstance(B, int):
        raise ValueError
        print("Indices must be integers.")

    # Defining list of manipulated atoms
    list_of_atoms = []

    # Retrieve calculator information as forces array
    # needs to be adjusted.
    prev_calc = model.get_calculator()

    prev_calc_results = prev_calc.results
    f = prev_calc_results["forces"]

    # Other properties require rearrangement too.
    t = model.get_tags()
    # TODO: magmoms + others

    # Atom count used later
    no_atoms = len(t)

    # Generate list of atoms that will be rearranged
    for i in range(0, no_atoms):
        list_of_atoms += [model[i]]

    # Rearrange indices and Atoms object properties
    list_of_atoms[A], list_of_atoms[B] = list_of_atoms[B], list_of_atoms[A]
    t[A], t[B] = t[B], t[A]
    # Ensure function works if force information empty
    if not f == []:
        f[A], f[B] = f[B], f[A]

    # Generate a new model based on switched indices
    new_model = Atoms(list_of_atoms,
                      pbc=model.get_pbc(),
                      cell=model.get_cell(),
                      tags=t,
                      constraint=model._get_constraints(),
                      calculator=prev_calc)

    # Trick calculator check_state by replacing atoms information
    # Can now use energy and forces as no changes in geometry detected
    prev_calc.atoms = new_model
    

    # User can interact with the new model
    return new_model


def check_interpolation(initial, final, n_max):
    '''
    Interpolates the provided geometries with n_max total images
    and checks whether any bond lengths below 0.74 Angstrom exist
    saves the interpolation in interpolation.traj

    # TODO: incorporate ase.neighborlist.natural_cutoff
    # for abnormal bond lengths based on typical A-B bonds

    Parameters:

    initial: Atoms object or string
        If a string, e.g. 'initial.traj', a file of this name will
        be read. Starting geometry for interpolation.
    final: Atoms object or string
        If a string, e.g. 'final.traj', a file of this name will
        be read. End point geometry for interpolation
    n_max: integer
        Desired total number of images for the interpolation
        including start and end point.
    '''

    from ase.neb import NEB
    from software.analyse.Interatomic_distances.analyse_bonds import search_abnormal_bonds
    from ase.io.trajectory import Trajectory
    from ase.io import read

    # Pre-requirements
    if isinstance(initial, str) is True:
        initial = read(initial)
    if isinstance(final, str) is True:
        final = read(final)
    if not isinstance(n_max, int):
        raise ValueError
        print("Max number of images must be an integer.")

    # Make a band consisting of 10 images:
    images = [initial]
    images += [initial.copy() for i in range(n_max-2)]
    images += [final]
    neb = NEB(images, climb=True)
    # Interpolate linearly the potisions of the middle images:
    neb.interpolate()

    t = Trajectory('interpolation.traj', 'w')

    for i in range(0, n_max):
        print("Assessing image", str(i+1) + '.')
        search_abnormal_bonds(images[i])
        t.write(images[i])

    t.close()


'''
def switch_indices_old(model, A, B, output_file):
    ## Takes 4 arguments: filename/Atoms object,
    ## indices A,B to switch around
    ## and output_file 'name.traj'
    from ase.io.trajectory import Trajectory

    if isinstance(model, str) is True:
       model = read(model)

    e=model.get_potential_energy()
    f=model.get_forces()
    new_model=Atoms()

    # index manipulation
    if not isinstance(A, int) and isinstance(B, int):
        raise ValueError
        print("Indices must be integers.")

    elif A ==B:
        print("Chosen idices are identical.")

    elif B < A:
            A,B = B,A # make sure indices are in the right order

    for i in range(0,A):
        new_model += model[i]
    new_model +=model[B]
    for i in range(A+1,B):
        new_model += model[i]
    new_model += model[A]
    for i in range(B+1,len(model.get_tags())):
        new_model += model[i]
    # Copy parameters into the new model
    new_model.set_cell(model.get_cell())
    new_model.set_pbc(model.get_pbc())
    new_model.set_constraint(model._get_constraints())

    # writer section
    #write('Model_corrected.traj', new_model)
    t1=Trajectory(str(output_file), 'w')
    t1.write(new_model, energy=e, forces=f)
    t1.close()
'''
