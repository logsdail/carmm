import numpy as np
import os
from ase.io.trajectory import Trajectory
from ase.io import read
from ase.geometry.analysis import Analysis
from ase import Atoms

## Use analyse_all_bonds and analyse_all_angles unless you are interested
## in specific values that are hard to get using visualizer.
## The functions take either a filename or a Atoms object.


def analyse_all_bonds(model):
    from itertools import combinations_with_replacement #combination as AB = BA for bonds, avoiding redundancy
    if isinstance(model, str) is True: #read file or Atoms object
        model = read(model)

    analysis = Analysis(model)
    dash = "-"*40
    list_of_symbols = list(set(model.get_chemical_symbols())) #set to ensure unique chemical symbols list
    all_bonds = combinations_with_replacement(list_of_symbols, 2)

    print(dash) #Table heading
    print('{:<6.5s}{:<6.5s}{:>4.10s}{:^13.10s}{:>4.10s}'.format("Bond","Count","Average", "Minimum", "Maximum"))
    print(dash)

    for bonds in all_bonds:     #iterate over all arrangements of chemical symbols
        A = bonds[0]
        B = bonds[1]

        print_AB = A+'-'+B
        AB_Bonds = analysis.get_bonds(A,B)

        if not AB_Bonds == [[]]:   #make sure bond exist before retrieving values
            AB_BondsValues = analysis.get_values(AB_Bonds)
            print('{:<8.8s}{:<6.0f}{:>4.6f}{:^12.6f}{:>4.6f}'.format(
            print_AB,len(AB_BondsValues[0]),np.average(AB_BondsValues),np.amin(AB_BondsValues),np.amax(AB_BondsValues)))

def analyse_all_angles(model):
    from itertools import product #product to get all possible arrangements
    if isinstance(model, str) is True: #read file or Atoms object
        model = read(model)

    analysis = Analysis(model)
    dash = "-"*40
    list_of_symbols = list(set(model.get_chemical_symbols())) #set to ensure unique chemical symbols list
    all_angles = product(list_of_symbols, repeat=3)

    print(dash)     #Table heading
    print('{:<8.8s}{:<6.5s}{:>4.10s}{:^13.10s}{:>4.10s}'.format("Angle","Count","Average", "Minimum", "Maximum"))
    print(dash)

    for angles in all_angles:     #iterate over all arrangements of chemical symbols
        A = angles[0]
        B = angles[1]
        C = angles[2]

        print_ABC = A+'-'+B+'-'+C
        ABC_Angle = analysis.get_angles(A,B,C)

        if not ABC_Angle == [[]]:   #make sure angles exist before retrieving values
            ABC_AngleValues = analysis.get_values(ABC_Angle)
            print('{:<8.8s}{:<6.0f}{:>4.4f}{:^12.4f}{:>4.4f}'.format(
            print_ABC,len(ABC_Angle[0]), np.average(ABC_AngleValues),np.amin(ABC_AngleValues),np.amax(ABC_AngleValues)))


# check A-B distance and A-B-C angles
# A,B,C are chemical symbols
def analyse_bonds(model, A, B):
    if isinstance(model, str) is True: #read file or Atoms object
        model = read(model)
    analysis = Analysis(model)
    dash = "-"*40

    print_AB = A +"-"+B
    AB_Bonds = analysis.get_bonds(A,B)
    AB_BondsValues = analysis.get_values(AB_Bonds)

    print(dash)
    print(print_AB+"       Distance / Angstrom")
    print(dash)
    print('{:<6.5s}{:>4.10s}{:^13.10s}{:>4.10s}'.format("count","average", "minimum", "maximum"))
    print('{:<6.0f}{:>4.6f}{:^12.6f}{:>4.6f}'.format(len(AB_BondsValues[0]),np.average(AB_BondsValues),np.amin(AB_BondsValues),np.amax(AB_BondsValues)))

def analyse_angles(model, A, B, C):
    if isinstance(model, str) is True: #read file or Atoms object
        model = read(model)
    analysis = Analysis(model)
    dash = "-"*40

    print_ABC = A +"-"+B+"-"+C

    ABC_Angle = analysis.get_angles(A,B,C)
    ABC_AngleValues = analysis.get_values(ABC_Angle)
    print(dash)
    print(print_ABC+"       Angle / Degrees")
    print(dash)
    print('{:<6.5s}{:>4.10s}{:^13.10s}{:>4.10s}'.format("count","average", "minimum", "maximum"))
    print('{:<6.0f}{:>4.4f}{:^12.4f}{:>4.4f}'.format(len(ABC_Angle[0]), np.average(ABC_AngleValues),np.amin(ABC_AngleValues),np.amax(ABC_AngleValues)))

def search_abnormal_bonds(model):
    from itertools import combinations_with_replacement #combination as AB = BA for bonds, avoiding redundancy
    if isinstance(model, str) is True: #read file or Atoms object
        model = read(model)

    abnormal_bonds = []         #define variables
    list_of_abnormal_bonds = []

    analysis = Analysis(model)
    list_of_symbols = list(set(model.get_chemical_symbols())) #set to ensure unique chemical symbols list
    all_bonds = combinations_with_replacement(list_of_symbols, 2)

    for bonds in all_bonds:     #iterate over all arrangements of chemical symbols
        A = bonds[0]
        B = bonds[1]

        print_AB = A+'-'+B
        AB_Bonds = analysis.get_bonds(A,B)

        if not AB_Bonds == [[]]:   #make sure bond exist before retrieving values
            AB_BondsValues = analysis.get_values(AB_Bonds)

            for i in range(0,len(AB_BondsValues)):
                for values in AB_BondsValues[i]:
                    #print(values)
                    if values < 0.74:
                        abnormal_bonds += [1]
                        list_of_abnormal_bonds = list_of_abnormal_bonds + [print_AB]

    if not len(abnormal_bonds)==0:
        print("A total of", len(abnormal_bonds), "abnormal bond lengths observed (<0.74 A).")
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
