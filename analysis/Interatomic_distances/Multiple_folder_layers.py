import fnmatch
import os

for layer1 in os.listdir('.'):                              #list of files in the current directory
    if fnmatch.fnmatch(layer1, 'FCC*'):                     #match files to the desired name, * works as a wildcard
        os.chdir(layer1)                                    #Enter the 1st layer of directories
        for sites in ["bridge", "fcc", "hcp", "ontop"]:     #loop over the next layers, here expressed as a list
            os.chdir(sites)
                                                            #(or do same as before - listdir+fnmatch)

            for layer3 in os.listdir('.'):                  #repeat step 1 for 3rd layer
              if fnmatch.fnmatch(layer3, "CO2onPd*"):       #Find a match for a filename* (wildcard useful if you like keeping the slurm ID on your folders)
                  os.chdir(layer3)
                                                            #add more layers as necessary
                  #print(os.getcwd())                       #unhash print to see if the script is entering the directories correctly

                                                            #insert any script, e.g. bond analysis, file conversion etc.

                  os.chdir("..")                            #change to previous directory

            os.chdir("..")
        os.chdir("..")
