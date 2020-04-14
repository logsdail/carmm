import numpy as np
from ase.io import read
from ase.geometry.analysis import Analysis


def analyse_all_bonds(model):
    '''
    Returns a table of bond distance analysis for the supplied model.

    Parameters:

    model: Atoms object or string. If string it will read a file
    in the same folder, e.g. "name.traj"
    '''
    # Combination as AB = BA for bonds, avoiding redundancy
    from itertools import combinations_with_replacement

    # Read file or Atoms object
    if isinstance(model, str) is True:
        model = read(model)

    analysis = Analysis(model)
    dash = "-" * 40

    # set() to ensure unique chemical symbols list
    list_of_symbols = list(set(model.get_chemical_symbols()))
    all_bonds = combinations_with_replacement(list_of_symbols, 2)

    # Table heading
    print(dash)
    print('{:<6.5s}{:<6.5s}{:>4.10s}{:^13.10s}{:>4.10s}'.format(
        "Bond", "Count", "Average", "Minimum", "Maximum"))
    print(dash)

    # Iterate over all arrangements of chemical symbols
    for bonds in all_bonds:
        A = bonds[0]
        B = bonds[1]

        print_AB = A+'-'+B
        AB_Bonds = analysis.get_bonds(A, B)

        # Make sure bond exist before retrieving values, then print contents
        if not AB_Bonds == [[]]:
            AB_BondsValues = analysis.get_values(AB_Bonds)
            print('{:<8.8s}{:<6.0f}{:>4.6f}{:^12.6f}{:>4.6f}'.format(
                print_AB, len(AB_BondsValues[0]), np.average(AB_BondsValues),
                np.amin(AB_BondsValues), np.amax(AB_BondsValues)))


def analyse_all_angles(model):
    '''
    Returns a table of bond angle analysis for the supplied model.

    Parameters:

    model: Atoms object or string. If string it will read a file
    in the same folder, e.g. "name.traj"

    '''

    # Product to get all possible arrangements
    from itertools import product
    # Read file or Atoms object
    if isinstance(model, str) is True:
        model = read(model)

    analysis = Analysis(model)
    dash = "-" * 40
    # set() to ensure unique chemical symbols list
    list_of_symbols = list(set(model.get_chemical_symbols()))
    all_angles = product(list_of_symbols, repeat=3)

    # Table heading
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
            print('{:<9.8s}{:<6.0f}{:>4.4f}{:^12.4f}{:>4.4f}'.format(
                print_ABC, len(ABC_Angle[0]), np.average(ABC_AngleValues),
                np.amin(ABC_AngleValues), np.amax(ABC_AngleValues)))


def analyse_bonds(model, A, B):
    '''
    Check A-B distances present in the model.
        model: Atoms object or string. If string it will read a file
        in the same folder, e.g. "name.traj"

    Parameters:
    A: string, chemical symbol, e.g. "H"
    B: string, chemical symbol, e.g. "H"
    '''
    # Read file or Atoms object
    if isinstance(model, str) is True:
        model = read(model)

    analysis = Analysis(model)
    dash = "-" * 40
    print_AB = A + "-" + B
    # Retrieve bonds and values
    AB_Bonds = analysis.get_bonds(A, B)
    AB_BondsValues = analysis.get_values(AB_Bonds)
    # Table header
    print(dash)
    print(print_AB+"       Distance / Angstrom")
    print(dash)
    print('{:<6.5s}{:>4.10s}{:^13.10s}{:>4.10s}'.format(
        "count", "average", "minimum", "maximum"))
    # Table contents
    print('{:<6.0f}{:>4.6f}{:^12.6f}{:>4.6f}'.format(
        len(AB_BondsValues[0]), np.average(AB_BondsValues),
        np.amin(AB_BondsValues), np.amax(AB_BondsValues)))


def analyse_angles(model, A, B, C):
    '''
    Check A-B distances present in the model.

    Parameters:
    model: Atoms object or string. If string it will read a file
        in the same folder, e.g. "name.traj"
    A: string, chemical symbol, e.g. "O"
    B: string, chemical symbol, e.g. "C"
    C: string, chemical symbol, e.g. "O"
    '''
    # Read file or Atoms object
    if isinstance(model, str) is True:
        model = read(model)

    analysis = Analysis(model)
    dash = "-"*40
    print_ABC = A + "-" + B + "-" + C
    # Retrieve bonds and values
    ABC_Angle = analysis.get_angles(A, B, C)
    ABC_AngleValues = analysis.get_values(ABC_Angle)
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


def search_abnormal_bonds(model, verbose=True):
    '''
    Check all bond lengths in the model for abnormally
    short ones, ie. less than 0.74 Angstrom.

    Parameters:
    model: Atoms object or string. If string it will read a file
        in the same folder, e.g. "name.traj"
    '''

    # Combination as AB = BA for bonds, avoiding redundancy
    from itertools import combinations_with_replacement
    # Imports necessary to work out accurate minimum bond distances
    from ase.data import chemical_symbols, covalent_radii

    # Read file or Atoms object
    if isinstance(model, str) is True:
        model = read(model)

    # Define lists of variables
    abnormal_bonds = []
    list_of_abnormal_bonds = []

    analysis = Analysis(model)
    # set() to ensure unique chemical symbols list
    list_of_symbols = list(set(model.get_chemical_symbols()))
    all_bonds = combinations_with_replacement(list_of_symbols, 2)

    # Iterate over all arrangements of chemical symbols

    for bonds in all_bonds:
        A = bonds[0]
        B = bonds[1]
        # For softcoded bond cutoff
        sum_of_covalent_radii = covalent_radii[chemical_symbols.index(A)]+covalent_radii[chemical_symbols.index(B)]

        print_AB = A+'-'+B
        AB_Bonds = analysis.get_bonds(A, B)

        # Make sure bond exist before retrieving values
        if not AB_Bonds == [[]]:
            AB_BondsValues = analysis.get_values(AB_Bonds)

            for i in range(0, len(AB_BondsValues)):
                for values in AB_BondsValues[i]:
                    # TODO: move the 75% of sum_of_covalent_radii before the loops
                    if values < max(0.4, sum_of_covalent_radii*0.75):
                        abnormal_bonds += [1]
                        list_of_abnormal_bonds = list_of_abnormal_bonds + [print_AB]

    # Abnormality check
    # is it possible to make a loop with different possible values instead of 0.75 and takes the average
    if len(abnormal_bonds) > 0:
        if verbose:
            print("A total of", len(abnormal_bonds),
            "abnormal bond lengths observed (<" + str(max(0.4, sum_of_covalent_radii*0.75)) + " A).")
            print("Identities:", list_of_abnormal_bonds)
        return False
    else:
        if verbose:
            print("OK")
        return True

def compare_structures(atoms1, atoms2):
    '''

    Description

    Parameters:

    atoms1: Atoms object or trajectory of individual atoms

    atoms2: Atoms object or trajectory of individual atoms

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
            temp_distance_sq = ((atoms2.positions[j][0] - xyz[0]) * (atoms2.positions[j][0] - xyz[0])
                                + (atoms2.positions[j][1] - xyz[1]) * (atoms2.positions[j][1] - xyz[1])
                                + (atoms2.positions[j][2] - xyz[2]) * (atoms2.positions[j][2] - xyz[2]))

            if distance_sq > temp_distance_sq and atoms1.symbols[i] == atoms2.symbols[j]:
                distance_sq = temp_distance_sq
                temp_index = j

        atoms2_indices.append(temp_index)
        differences.append(sqrt(distance_sq))

    return atoms2_indices, differences

def get_indices_of_elements(list_of_symbols, symbol):
    '''

    Description

    Parameters:

    list_of_symbols:

    symbol:
    '''
    return [i for i, x in enumerate(list_of_symbols) if x == symbol.capitalize()]



'''
## not working as intended as specific indices are needed
def analyse_dihedrals(model):
    from itertools import product
    if isinstance(model, str) is True: #read file or Atoms object
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
