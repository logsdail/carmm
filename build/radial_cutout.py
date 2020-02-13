import numpy as np
from ase.io.trajectory import Trajectory
from ase.io import read
from ase.build import molecule
from ase import Atoms
from ase.visualize import view
from ase.constraints import FixAtoms

# Read your trajectory/geometry file 
atoms = read("geometry_supercell.in")

#count number of atoms
tags = atoms.get_tags()
n = len(tags)

#Read atomic positions
coordinates = atoms.get_positions()

#generate a mask based on distance from neighbours
mask =[]
for i in range(n):

    #get distances between atom of interest and others - then constrain
    #all atoms beyond a certain radius
    
    #################### Edit atom tag ###############################
    ABdistance=np.linalg.norm((coordinates[237] - coordinates[i]))
    
    ################ Edit distance in Angstrom here ###################
    if ABdistance > 8.0:
    ###################################################################
        select = 1
        mask.append(select)
    else:
        select = 0
        mask.append(select)

atoms.set_constraint()
atoms.set_constraint(FixAtoms(mask = mask))
view(atoms)

## in ase gui select all constrained atoms and delete!
