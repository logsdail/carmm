def switch_indices(model, A, B):
    '''
    Function for rearranging atomic indices in the structure.
    Returns a changed atoms object with retained calculator
    and ensures array of forces is rearranged accordingly.

    Parameters:

    model: Atoms object
        Structure requiring atom index rearrangement.
    A, B: integer
        Indices of atoms that need to be switched

    '''

    # Index manipulation
    if not isinstance(A, int) and isinstance(B, int):
        raise ValueError
        print("Indices must be integers.")

    # Defining list of manipulated atoms
    t = [atom.index for atom in model]

    # Retrieve calculator information
    # TODO: Move this after the creation of the atoms object, so we reduce if statements.
    if model.get_calculator() is not None:
        # If it exists, forces array needs to be adjusted.
        prev_calc = model.get_calculator()
        prev_calc_results = prev_calc.results
        if "forces" in prev_calc_results:
            f = prev_calc_results["forces"]
        else:
            f = []
    else:
        f = []

    t[A], t[B] = t[B], t[A]
    # Ensure function works if force information empty
    # TODO: Where is this information subsequently used?
    if not f == []:
        f[A], f[B] = f[B], f[A]

    # Generate a new model based on switched indices
    new_model = model[t]

    if model.get_calculator() is not None:
        new_model.set_calculator(prev_calc)
        # Trick calculator check_state by replacing atoms information
        # Can now use energy and forces as no changes in geometry detected
        prev_calc.atoms = new_model

    # User can interact with the new model
    return new_model

def switch_all_indices(model, new_indices):
    '''
    Method to update all indices in a model all at once,
    assuming the new indexing order is available from e.g. compare_structures

    Parameters:

    model: ASE atoms object
        Input structure to be reordered
    new_indices: List of Integers
        New indices for each atom in model

    TODO: Add an example to the QA tests.
    '''

    # List of indices that have been swapped
    swapped = []

    # TODO: Add a sanity check to make sure the length of new_indices matches the size
    #       of the model!
    for i in range(len(new_indices)):
        if new_indices[i] is not i and i not in swapped:
            model = switch_indices(model, i, new_indices[i])
            swapped.append(new_indices[i])

    return model

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
