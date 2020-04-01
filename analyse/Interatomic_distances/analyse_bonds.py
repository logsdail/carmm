import numpy as np
from ase.io import read
from ase.geometry.analysis import Analysis


def analyse_all_bonds(model):
    '''
    Returns a table of bond distance analysis for the supplied model.
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


def search_abnormal_bonds(model):
    '''
    Check all bond lengths in the model for abnormally
    short ones, ie. less than 0.74 Angstrom.
        model: Atoms object or string. If string it will read a file
        in the same folder, e.g. "name.traj"
    '''

    # Combination as AB = BA for bonds, avoiding redundancy
    from itertools import combinations_with_replacement
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

        print_AB = A+'-'+B
        AB_Bonds = analysis.get_bonds(A, B)

        # Make sure bond exist before retrieving values
        if not AB_Bonds == [[]]:
            AB_BondsValues = analysis.get_values(AB_Bonds)

            for i in range(0, len(AB_BondsValues)):
                for values in AB_BondsValues[i]:
                    if values < 0.74:
                        abnormal_bonds += [1]
                        list_of_abnormal_bonds = list_of_abnormal_bonds + [print_AB]

    # Abnormality check
    if not len(abnormal_bonds) == 0:
        print("A total of", len(abnormal_bonds),
        "abnormal bond lengths observed (<0.74 A).")
        print("Identities:", list_of_abnormal_bonds)
    else:
        print("OK")




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
