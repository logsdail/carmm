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
    t[A], t[B] = t[B], t[A]

    # Generate a new model based on switched indices
    new_model = model[t]

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

    # Prevents unexpected editing of parent object in place
    # Now ensures returned object is different to incoming atoms
    new_model = model.copy()

    # List of indices that have been swapped
    swapped = []

    assert len(new_indices) == len(model), "Wrong number of indices provided"

    for i in range(len(new_indices)):
        if new_indices[i] is not i and i not in swapped:
            new_model = switch_indices(new_model, i, new_indices[i])
            swapped.append(new_indices[i])

    return new_model

def sort_by_symbols(model):
    '''
    Method for arranging the indices in the Atoms object based on their chemical
    symbols, enables interpolation if 'Atomic ordering' flags were raised.

    Parameters:
        model - Atoms object
    '''

    import copy
    blank = copy.deepcopy(model)
    del blank[range(len(model))]
    symbols = list(set(model.get_chemical_symbols()))
    symbols.sort()
    
    for i in symbols:
        blank += model[[atom.index for atom in model if atom.symbol == i]]

    return blank
