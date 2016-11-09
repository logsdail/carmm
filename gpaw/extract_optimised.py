from ase.io import write, PickleTrajectory
import os

#Read atoms
for file in os.listdir("."):
    if file.endswith(".traj"):
        traj=PickleTrajectory(file)
        atoms=traj[-1]
	write(file[:-5]+"_opt.xyz",atoms)
