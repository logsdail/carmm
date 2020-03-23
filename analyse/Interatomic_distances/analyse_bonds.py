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

    print("there are {} "+print_AB+" bonds in BETA".format(len(AB_Bonds[0])))
    print("there are {} "+print_ABC+" angles in BETA".format(len(ABC_Angle[0])))

    AB_BondsValues = analysis.get_values(AB_Bonds)
    ABC_AngleValues = analysis.get_values(ABC_Angle)

    print("bond length data:")
    print("the average "+print_AB+" bond length is {}.".format(np.average(AB_BondsValues)))
    print("the minimum "+print_AB+" Distance is:", np.amin(AB_BondsValues))
    print("the maximum "+print_AB+" Distance is:", np.amax(AB_BondsValues))

    print("bond angle data:")
    print("the average "+print_ABC+" angle is {}.".format(np.average(ABC_AngleValues)))
    print("the maximum "+print_ABC+" angle is:", np.amax(ABC_AngleValues))
    print("the minimum "+print_ABC+" angle is:", np.amin(ABC_AngleValues))
