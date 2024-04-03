def analyse_all_angles(model, verbose=True):
    '''
    Returns a table of bond angle analysis for the supplied model.

    Parameters:

    model: Atoms object
        The structure that needs to be interrogated
    verbose: Boolean
        Whether to print information to screen
    Returns:
        - list of all elemental combinations
        - list of indices for each elemental combination
        - list of all angle values for each combination of indices
    '''

    # set() to ensure unique chemical symbols list
    list_of_symbols = list(set(model.get_chemical_symbols()))
    # Product to get all possible arrangements
    from itertools import product
    angles_element_combinations = product(list_of_symbols, repeat=3)

    # Table heading
    if verbose:
        print_angles_table_header()

    angles_elements = []
    angles_indices = []
    angles_values = []
    # Iterate over all arrangements of chemical symbols
    for angles in angles_element_combinations:
        ABC_indices, ABC_values = analyse_angles(model, angles[0], angles[1], angles[2], verbose=verbose, multirow=True)
        if len(ABC_indices[0]) > 0:
            angles_elements.append(angles)
            angles_indices.append(ABC_indices[0])
            angles_values.append(ABC_values[0])

    return angles_elements, angles_indices, angles_values

def analyse_angles(model, A, B, C, verbose=True, multirow=False):
    '''
    Check A-B-C angles present in the model.

    Parameters:
    model: Atoms object
        The ASE atoms object that needs to be considered 
    A: string, chemical symbol, e.g. "O"
    B: string, chemical symbol, e.g. "C"
    C: string, chemical symbol, e.g. "O"
    verbose: Boolean
        Whether to print information to screen
    multirow: Boolean
        Whether we are returning multiple sets of results in a Table
    '''

    from ase.geometry.analysis import Analysis
    analysis = Analysis(model)

    print_ABC = A + "-" + B + "-" + C
    # Retrieve bonds and values
    ABC_indices = analysis.get_angles(A, B, C)
    if len(ABC_indices[0]) == 0:
        ABC_values = None
    else:
        ABC_values = analysis.get_values(ABC_indices)

    if verbose and ABC_values is not None:
        # Table header
        if not multirow:
            print_angles_table_header()
        # Table contents
        import numpy as np
        print('{:<9.8s}{:<6.0f}{:>4.4f}{:^12.4f}{:>4.4f}'.format(
            print_ABC, len(ABC_indices[0]), np.average(ABC_values),
            np.amin(ABC_values), np.amax(ABC_values)))

    return ABC_indices, ABC_values

def print_angles_table_header():
    print("-" * 40)
    print('{:<9.8s}{:<6.5s}{:>4.10s}{:^13.10s}{:>4.10s}'.format(
        "Angle", "Count", "Average", "Minimum", "Maximum"))
    print("-" * 40)

'''
## not working as intended as specific indices are needed
def analyse_dihedrals(model):
    from itertools import product
    if isinstance(model, str) is True: #read file or Atoms object
        from ase.io import read
        model = read(model)

    analysis = Analysis(model)
    dash = "-"*40
    list_of_symbols = list(set(model.get_chemical_symbols())) #set to ensure unique chemical symbols list
    all_angles = product(list_of_symbols, repeat=4)

    print(dash)     #Table heading
    print('{:<6.5s}{:<6.5s}{:>4.10s}{:^13.10s}{:>4.10s}'.format("Angle","Count","Average", "Minimum", "Maximum"))
    print(dash)

    for angles in all_angles:     #iterate over all arrangements of chemical symbols
        A = angles[0]
        B = angles[1]
        C = angles[2]
        D = angles[3]

        print_ABC = A+'-'+B+'-'+C+'-'+D
        ABC_Angle = analysis.get_dihedrals(A,B,C,D)

        if not ABC_Angle == [[]]:   #make sure angles exist before retrieving values
            ABC_AngleValues = analysis.get_values(ABC_Angle)
            print('{:<6.8s}{:<6.0f}{:>4.4f}{:^12.4f}{:>4.4f}'.format(
            print_ABC,len(ABC_Angle[0]), np.average(ABC_AngleValues),np.amin(ABC_AngleValues),np.amax(ABC_AngleValues)))
'''
