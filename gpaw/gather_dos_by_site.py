#/usr/bin/python

import sys
import os
import numpy
import pickle
import matplotlib.mlab as mlab

# This all assumes we are doing everything spin-paired
angular_contributions = 'spdf'
# This is the variance if we want to expand our data with gaussian
# 0.02 gives us a nice plot if we want it
gaussian_variance = 0.0

def write_dos_to_file(ef,energy,dos,output):
    pickle.dump((ef, energy, dos), open(output, 'w'))
    
    # Included for Shiny - writing everything as dat files
    with open(output+".dat", "w") as dat_file:
        dat_file.write("Ef: {}\n".format(ef))
        for i in range(len(energy)):
            dat_file.write('{} , {}\n'.format(energy[i],dos[i]))

def write_dos_to_files(ef,energy,dos_data,dos_filename):

    if gaussian_variance > 0.0:
        new_energy, new_dos_data = expand_with_gaussians(energy,dos_data)
    else:
        new_energy = energy
        new_dos_data = dos_data

    write_dos_to_file(ef,new_energy,new_dos_data[0],dos_filename)

    # Perhaps here we could break this down to import the angular contributions....
    for angular in range(len(angular_contributions)):
        angular_dos_filename = dos_filename+"."+angular_contributions[angular]
        write_dos_to_file(ef,new_energy,new_dos_data[angular+1],angular_dos_filename) 

def expand_with_gaussians(energy,dos_data):

    # We need to make space on the energy scale
    # Otherwise the gaussian expansion comes to an abrupt finish
    new_size = int(len(energy)*1.1)
    new_energy = [0]*new_size
    new_dos_data = [[0 for i in range(new_size)] for j in range(len(dos_data))]
    offset = (new_size-len(energy))/2
    stepsize = abs(energy[0]-energy[1])
    start_energy = min(energy) - (offset*stepsize)

    # Set the new energy scale
    for i in range(len(new_energy)):
        new_energy[i] = start_energy + (i*stepsize)

    for i in range(len(dos_data)):
        # OK so now prepare the gaussian expansion
        new_values = [0]*new_size
        # Values for normal distribution
        sigma = numpy.sqrt(gaussian_variance)
        x = numpy.linspace(min(new_energy),max(new_energy),new_size)

        # Expand the previous energy data, and scale in amplitude by dos contributions
        for j in range(len(energy)):
            new_values += mlab.normpdf(x,energy[j],sigma)*dos_data[i][j] #This scales the data to the current DOS content

        # Scale all values to reasonable size
        if sum(dos_data[i]) > 0:
            scalar = sum(dos_data[i])/sum(new_values)
        else:
            scalar = 0
        new_values *= scalar
        new_dos_data[i] = numpy.copy(new_values)

    # Return the newly scaled data
    return new_energy, new_dos_data

def read_dos_from_file(filename):
    ef = 0
    energy = []
    dos = []

    try:
        ef, energy, dos = pickle.load(open(filename))
    except:
        try:
            # See if it is a spin-polarised data set
            ef, energy, dos_spin1 = pickle.load(open(filename+".spin1"))
            ef, energy, dos_spin2 = pickle.load(open(filename+".spin2"))
            # Combine them and return
            dos = dos_spin1 + dos_spin2
        except:
            print "A problem occurred reading the dos from ", filename
            sys.exit(1)

    return ef, energy, dos

def read_dos_from_files(dos_filename,dos_data):
    # Read in information for dos
    ef, energy, dos = read_dos_from_file(dos_filename)

    # Add total dos to all other total dos quantites
    dos_data[0] += dos

    # Perhaps here we could break this down to import the angular contributions....
    for angular in range(len(angular_contributions)):
        angular_dos_filename = dos_filename+"."+angular_contributions[angular]
        ef, energy, dos = read_dos_from_file(angular_dos_filename)
        dos_data[angular+1] += dos # +1 to offest from zero, which is the total dos

    # At this point we have read in all the dos information.    
    return dos_data

def check_dos_data_is_correct(dos_data):
    sum_residual=abs((sum(dos_data[0]) - sum(dos_data[1] + dos_data[2] + dos_data[3] + dos_data[4])))
    if sum_residual > 1e-8:
        print "A problem occured during checks. Residual of differences is ", sum_residual
        sys.exit(1)

# Partition dos into regions: core/face/edge/vertex

# Read input, make sure it is valid
# Find out number of atoms for comparison with template

inputs = sys.argv[1:]
template = inputs[0]
dos_prefix = inputs[1]
core_shell = False

# Allow definition of template species
if len(inputs) > 2:
    species = inputs[2:]
    output_labels = inputs[2:]
# Otherwise assume we are working with core/shell system
else:
    species = ['H','He','Li','Be']
    output_labels = ['core','face','edge','vertex']
    core_shell = True

# Determine if a template is provided, if not use the one suggested
# Assume H/He/Li/Be labelling is consistent with original version
# Condense input file template into list containing just species labels if present

# File to hold all the information from the template
working_template = []

if os.path.exists(template):
    with open(template, 'r') as temp_template:
        template_lines = temp_template.readlines()  
        for item in template_lines[2:]:
            working_template.append(item.split()[0])
else:
    print "Error: No template ", template, " available"
    sys.exit(1)

# Break up the list of charges into regions.
# Define our common variables for reading in data
ef = 0
energy = []
dos = []
# We will make 5 arrays for each site type. 
# The first array to hold the total dos, and the next four to hold the angular contributions
# We can check the the angular contributions add up to match the total dos at the end
# This information will be held in a 3D array (hopefully)

# Assumptions:
# Ef and energy do not change as we read in the data.
# This will be correct for each file

for atom_number in range(len(working_template)):

    # Organise dos filename
    dos_filename = dos_prefix+"."+str(atom_number)

# We need to initialise our arrays, so lets do this for the first information
    if atom_number == 0:
        ef, energy, dos_temp = read_dos_from_file(dos_filename)
        # We need to quickly resize the arrays so they are suitable for loading information
        dos = [[[0 for i in range(len(dos_temp))] for j in range(5)] for k in range(len(species))]
 
    identified = False
    for i in range(len(species)):
        if working_template[atom_number] == species[i]:
            identified = True
            # Read in information for dos and add to appropriate values
            dos[i] = read_dos_from_files(dos_filename,dos[i])

    if not identified: 
        print "A problem occurred during identification for site type ", working_template[atom_number]
        print "Provided species are ", species
        sys.exit(1)

# All information is loaded
# Checks and Outputs

for i in range(len(species)):
    check_dos_data_is_correct(dos[i])
    write_dos_to_files(ef,energy,dos[i],dos_prefix+"."+output_labels[i])

# We'll write the information to file (pickles), and then we can use our old scripts to load in the data we want!

# Combine shell information

if core_shell:
    shell_dos = numpy.copy(dos[1]) + numpy.copy(dos[2]) + numpy.copy(dos[3])
    check_dos_data_is_correct(shell_dos)
    write_dos_to_files(ef,energy,shell_dos,dos_prefix+".shell")

# Write total dos in angular decompositions
total_dos = [[0 for i in range(len(dos[0][0]))] for j in range(5)]
for i in range(len(species)):
    total_dos += numpy.copy(dos[i])

check_dos_data_is_correct(total_dos)
write_dos_to_files(ef,energy,total_dos,dos_prefix+".decomp")
