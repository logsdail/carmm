import numpy as np
from ase.io.trajectory import Trajectory
from ase.io import read
from ase.build import molecule
from ase import Atoms

# Read your trajectory file (last config)
model = read("name.traj")

#Read atomic positions
coordinates = model.get_positions()

#unhash the following to find the number of X atoms in model:
# need to figure this out

#unhash the following to identify atoms
#print(coordinates, model.get_atomic_numbers())

#write Distances in a .txt file
w = open("Distances.txt", "w+")

#iterate over n Pd atoms in a slab, n is an integer
for i in range(n):

    #get XYZ coordinates relative to atom of interest
    #get a norm of these coordinates = distances between A and atoms in a slab
    ABdistance=np.linalg.norm((coordinates[n] - coordinates[i]))
    ABdistance = str(ABdistance)

    #write the distance into a file, separated by a space
    w.write(ABdistance+" ")

#stop writing into a file once loop has finished
w.close()


r = open("Distances.txt", "r")          #read the created file containing distances
contents=r.read()

                                        #The Distances.txt file can be read in Excel and transposed

                                        #to make it usable it needs to be coverted into a list


contents = contents.split()             #split into a list
contents=map(float, contents)           #map the list for floats
list = [float(s) for s in contents]     #output a list of floats

print("The minimum C-Pd Distance is: ", np.amin(list))  #find the minimum distance
