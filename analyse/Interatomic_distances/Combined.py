import fnmatch
import os
import numpy as np
from ase.io.trajectory import Trajectory
from ase.io import read
from ase.build import molecule
from ase import Atoms

for layer1 in os.listdir('.'):
    if fnmatch.fnmatch(layer1, 'FCC*'):
        os.chdir(layer1)
        for sites in ["bridge", "fcc", "hcp", "ontop"]:
            os.chdir(sites)
            for layer3 in os.listdir('.'):
              if fnmatch.fnmatch(layer3, "CO2onPd*"):
                  os.chdir(layer3)

                  # Read trajectory file (last config)
                  model = read("CO2Pd.traj")

                  #Read atomic positions
                  coordinates = model.get_positions()

                  #write Distances in a .txt file
                  w = open("Distances.txt", "w+")

                  #iterate over 36 Pd atoms
                  for i in range(35):
                      #get XYZ coordinates relative to 37th ie. C atom
                      #get a norm of these coordinates = distance between C and all Pd atoms
                      CXdistance=np.linalg.norm((coordinates[36] - coordinates[i]))
                      CXdistance = str(CXdistance)

                      #write the distance into a file, separated by a space
                      w.write(CXdistance+" ")

                  #stop writing into a file once loop has finished
                  w.close()


                  r = open("Distances.txt", "r")          #read the created file containing distances
                  contents=r.read()


                  contents = contents.split()             #split into a list
                  contents=map(float, contents)           #map the list for floats
                  list = [float(s) for s in contents]     #output a string of floats

                  #printing current directory to identify my values when importing to excel
                  print(os.getcwd()+": The_minimum_C-Pd_Distance_is: ", np.amin(list))
                
                  #if your folder structure is relevant to your data table turn the path into a list instead, e.g.:
                  #current=os.getcwd()
                  #split = current.split("/")
                  #del split[0:3]                  
                  #path = split
                  #print(path, "The_minimum_C-Pd_Distance_is: ", np.amin(list))


                  os.chdir("..")

            os.chdir("..")
        os.chdir("..")
