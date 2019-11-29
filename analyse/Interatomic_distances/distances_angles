from ase.io.trajectory import Trajectory
import numpy as np
from ase.io import read
from ase.geometry.analysis import Analysis

BEA = read('BEA.cif')

traj = Trajectory('BEA.traj', 'w')
traj.write(BEA)
BEA = read('BEA.traj')

ana = Analysis(BEA)

SiOBonds = ana.get_bonds('Si', 'O')
SiOSiAngles = ana.get_angles('Si', 'O', 'Si')

print("there are {} Si-O bonds in BETA".format(len(SiOBonds[0])))
print("there are {} Si-O-Si angles in BETA".format(len(SiOSiAngles[0])))

SiOBondsValues = ana.get_values(SiOBonds)
SiOSiAngleValues = ana.get_values(SiOSiAngles)

print("bond length data:")
print("the average Si-O bond length is {}.".format(np.average(SiOBondsValues)))
print("the minimum Si-O Distance is:", np.amin(SiOBondsValues))
print("the maximum Si-O Distance is:", np.amax(SiOBondsValues))

print("bond angle data:")
print("the average Si-O-Si angle is {}.".format(np.average(SiOSiAngleValues)))
print("the maximum Si-O-Si angle is:", np.amax(SiOSiAngleValues))
print("the minimum Si-O-Si angle is:", np.amin(SiOSiAngleValues))
