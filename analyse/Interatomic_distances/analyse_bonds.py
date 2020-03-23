import numpy as np
import os
from ase.io.trajectory import Trajectory
from ase.io import read
from ase.geometry.analysis import Analysis
from ase import Atoms

# check A-B distance and A-B-C angles
# A,B,C are chemical symbols
def analyse_bonds(file, A, B, C):
    file = read(file)
    analysis = Analysis(file)

    print_AB = A +"-"+B
    print_ABC = A +"-"+B+"-"+C

    AB_Bonds = analysis.get_bonds(A,B)
    ABC_Angle = analysis.get_angles(A,B,C)

    AB_BondsValues = analysis.get_values(AB_Bonds)
    ABC_AngleValues = analysis.get_values(ABC_Angle)

    dash = "-"*40

    print(dash)
    print(print_AB+"       Distance / Angstrom")
    print(dash)
    print('{:<6.5s}{:>4.10s}{:^13.10s}{:>4.10s}'.format("count","average", "minimum", "maximum"))
    print('{:<6.0f}{:>4.6f}{:^12.6f}{:>4.6f}'.format(len(AB_BondsValues[0]),np.average(AB_BondsValues),np.amin(AB_BondsValues),np.amax(AB_BondsValues)))

    print(dash)
    print(print_ABC+"       Angle / Degrees")
    print(dash)
    print('{:<6.5s}{:>4.10s}{:^13.10s}{:>4.10s}'.format("count","average", "minimum", "maximum"))
    print('{:<6.0f}{:>4.4f}{:^12.4f}{:>4.4f}'.format(len(ABC_Angle[0]), np.average(ABC_AngleValues),np.amin(ABC_AngleValues),np.amax(ABC_AngleValues)))
