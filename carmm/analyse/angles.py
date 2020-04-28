def analyse_all_angles(model, verbose=True):
    '''
    Returns a table of bond angle analysis for the supplied model.
    TODO: - Setup method to return information
          - Setup routine functionality to use analyse_angless

    Parameters:

    model: Atoms object or string. If string it will read a file
    in the same folder, e.g. "name.traj"
    verbose: Boolean
        Whether to print information to screen
    '''
    import numpy as np
    # Product to get all possible arrangements
    from itertools import product
    # Read file or Atoms object
    if isinstance(model, str) is True:
        from ase.io import read
        model = read(model)

    from ase.geometry.analysis import Analysis
    analysis = Analysis(model)
    dash = "-" * 40
    # set() to ensure unique chemical symbols list
    list_of_symbols = list(set(model.get_chemical_symbols()))
    all_angles = product(list_of_symbols, repeat=3)

    # Table heading
    if verbose:
        print(dash)
        print('{:<9.8s}{:<6.5s}{:>4.10s}{:^13.10s}{:>4.10s}'.format(
            "Angle", "Count", "Average", "Minimum", "Maximum"))
        print(dash)

    # Iterate over all arrangements of chemical symbols
    for angles in all_angles:
        A = angles[0]
        B = angles[1]
        C = angles[2]

        print_ABC = A+'-'+B+'-'+C
        ABC_Angle = analysis.get_angles(A, B, C)

        # Make sure angles exist before retrieving values, print table contents
        if not ABC_Angle == [[]]:
            ABC_AngleValues = analysis.get_values(ABC_Angle)
            if verbose:
                print('{:<9.8s}{:<6.0f}{:>4.4f}{:^12.4f}{:>4.4f}'.format(
                    print_ABC, len(ABC_Angle[0]), np.average(ABC_AngleValues),
                    np.amin(ABC_AngleValues), np.amax(ABC_AngleValues)))

def analyse_angles(model, A, B, C, verbose=True):
    '''
    Check A-B distances present in the model.
    TODO: Setup method to return information

    Parameters:
    model: Atoms object or string. If string it will read a file
        in the same folder, e.g. "name.traj"
    A: string, chemical symbol, e.g. "O"
    B: string, chemical symbol, e.g. "C"
    C: string, chemical symbol, e.g. "O"
    verbose: Boolean
        Whether to print information to screen
    '''
    import numpy as np
    # Read file or Atoms object
    if isinstance(model, str) is True:
        from ase.io import read
        model = read(model)

    from ase.geometry.analysis import Analysis
    analysis = Analysis(model)
    dash = "-"*40
    print_ABC = A + "-" + B + "-" + C
    # Retrieve bonds and values
    ABC_Angle = analysis.get_angles(A, B, C)
    ABC_AngleValues = analysis.get_values(ABC_Angle)
    if verbose:
        # Table header
        print(dash)
        print(print_ABC+"       Angle / Degrees")
        print(dash)
        print('{:<6.5s}{:>4.10s}{:^13.10s}{:>4.10s}'.format(
            "count", "average", "minimum", "maximum"))
        # Table contents
        print('{:<6.0f}{:>4.4f}{:^12.4f}{:>4.4f}'.format(
            len(ABC_Angle[0]), np.average(ABC_AngleValues),
            np.amin(ABC_AngleValues), np.amax(ABC_AngleValues)))

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
