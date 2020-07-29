def analyse_all_bonds(model, verbose=True, abnormal=True):
    '''
    Returns all abnormal bond types and list of these
    TODO: Make this more bullet proof - what happens if abnormal bonds aren't requested.
    A table of bond distance analysis for the supplied model is also possible

    Parameters:

    model: Atoms object
        Structure for which the analysis is to be conducted
    verbose: Boolean
        Determines whether the output should be printed to screen
    abnormal: Boolean
        Collect information about rogue looking bond lengths.
        (Does enabling this by default add a large time overhead?)
    '''

    # set() to ensure unique chemical symbols list
    list_of_symbols = list(set(model.get_chemical_symbols()))
    # Combination as AB = BA for bonds, avoiding redundancy
    from itertools import combinations_with_replacement
    all_bonds = combinations_with_replacement(list_of_symbols, 2)

    # Define lists to collect abnormal observations
    abnormal_bonds = []
    list_of_abnormal_bonds = []

    # Table heading
    if verbose:
        print_bond_table_header()

    from ase.data import chemical_symbols, covalent_radii
    # Iterate over all arrangements of chemical symbols
    for bonds in all_bonds:
        print_AB, AB_Bonds, AB_BondsValues = analyse_bonds(model, bonds[0], bonds[1], verbose=verbose, multirow=True)

        if abnormal and AB_BondsValues is not None:
            sum_of_covalent_radii = covalent_radii[chemical_symbols.index(bonds[0])] + covalent_radii[chemical_symbols.index(bonds[1])]
            abnormal_cutoff = max(0.4, sum_of_covalent_radii*0.75)

            for values in AB_BondsValues:
                abnormal_values = [ i for i in values if i < abnormal_cutoff ]
                if len(abnormal_values):
                    # Why do we add this value of 1? unclear and not tested in regression.
                    # @Igor: Is this a counter? I can't tell, and the QA test isn't thorough enough to be clear
                    # If it is a counter, should it be changed to += len(abnormal_values)
                    abnormal_bonds.append(len(abnormal_values))
                    list_of_abnormal_bonds.append(print_AB)

    # This now returns empty arrays if no abnormal bond checks are done,
    # or if genuinely there are no abnormal bonds.
    return abnormal_bonds, list_of_abnormal_bonds

def analyse_bonds(model, A, B, verbose=True, multirow=False):
    '''
    Check A-B distances present in the model.

    Parameters:
    model: Atoms object
        XXX
    A: string, chemical symbol, e.g. "H"
    B: string, chemical symbol, e.g. "H"
    verbose: Boolean
        Whether to print information to screen
    multirow: Boolean
        Whether we are working with analyse_all_bonds, so the output is multirow,
        or just one specific analysis of a bond, in which case the table header is needed.
    '''

    from ase.geometry.analysis import Analysis
    analysis = Analysis(model)

    print_AB = A + "-" + B
    # Retrieve bonds and values
    AB_Bonds = analysis.get_bonds(A, B)
    if AB_Bonds == [[]]:
        AB_BondsValues = None
    else:
        AB_BondsValues = analysis.get_values(AB_Bonds)

    if verbose and AB_BondsValues is not None:
        if not multirow:
            print_bond_table_header()
        # Table contents
        import numpy as np
        print('{:<8.8s}{:<6.0f}{:>4.6f}{:^12.6f}{:>4.6f}'.format(
            print_AB, len(AB_BondsValues[0]), np.average(AB_BondsValues),
            np.amin(AB_BondsValues), np.amax(AB_BondsValues)))

    return print_AB, AB_Bonds, AB_BondsValues

def print_bond_table_header():
    print("-" * 40)
    print('{:<6.5s}{:<6.5s}{:>4.10s}{:^13.10s}{:>4.10s}'.format(
        "Bond", "Count", "Average", "Minimum", "Maximum"))
    print("-" * 40)

def search_abnormal_bonds(model, verbose=True):
    '''
    Check all bond lengths in the model for abnormally
    short ones, ie. less than 0.74 Angstrom.

    Parameters:
    model: Atoms object or string. If string it will read a file
        in the same folder, e.g. "name.traj"
    '''

    # Abnormality check
    abnormal_bonds, list_of_abnormal_bonds = analyse_all_bonds(model, verbose=verbose, abnormal=True)

    # is it possible to make a loop with different possible values instead of 0.75 and takes the average
    if len(abnormal_bonds) > 0:
        if verbose:
            print("-"*40)
            print("A total of", len(abnormal_bonds),
            "abnormal bond lengths observed (<" + str(max(0.4, sum_of_covalent_radii*0.75)) + " A).")
            print("Identities:", list_of_abnormal_bonds)
            print("-"*40)
        return False
    else:
        return True

def compare_structures(atoms1, atoms2, label=None):
    '''

    Comparison of two input structures to identify equivalent atoms but incorrect index ordering

    Parameters:

    atoms1: Atoms object or trajectory of individual atoms
        An atoms object
    atoms2: Atoms object or trajectory of individual atoms
        Another atoms object
    label: String of elemental character
        Only necessary to limit search to specific atomic species
    '''
    from math import sqrt

    if len(atoms1) != len(atoms2):
        print("The inputs don't contain the same number of atoms.")
        exit()

    # Configure arrays
    differences = []
    atoms2_indices = []

    # Iterate over indices of all atoms in structure 1 and compare to structure 2.
    for i in range(len(atoms1.positions)):
        xyz = atoms1.positions[i]
        distance_sq = 999999.9
        temp_index = 0
        for j in range(len(atoms2.positions)):
            if atoms1.symbols[i] == atoms2.symbols[j] and (atoms1.symbols[i] == label or label == None):
                temp_distance_sq = ((atoms2.positions[j][0] - xyz[0]) * (atoms2.positions[j][0] - xyz[0])
                                    + (atoms2.positions[j][1] - xyz[1]) * (atoms2.positions[j][1] - xyz[1])
                                    + (atoms2.positions[j][2] - xyz[2]) * (atoms2.positions[j][2] - xyz[2]))

                if distance_sq > temp_distance_sq:
                    distance_sq = temp_distance_sq
                    temp_index = j

        atoms2_indices.append(temp_index)
        differences.append(sqrt(distance_sq))

    return atoms2_indices, differences

def get_indices_of_elements(list_of_symbols, symbol):
    '''

    Check an atoms object for occurences of symbols given and return indices

    Parameters:

    list_of_symbols: List of strings
        Symbols from an atoms object in structural order
    symbol:
        Symbol to search for
    '''
    return [i for i, x in enumerate(list_of_symbols) if x == symbol.capitalize()]
