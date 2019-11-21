#/usr/bin/python

import sys
import math
import os
import numpy
from ase import Atoms
from ase.io import read
from ase.visualize import view

import random

solid_file = sys.argv[1]
new_species = sys.argv[2]
new_species_total = int(sys.argv[3])
if len(sys.argv) > 4:
    regions_file = sys.argv[4]
else:
    regions_file = ''

# Exit if the specified file is not present
if not os.path.exists(solid_file):
    print "Error: Input file not found"
    sys.exit(1)

atoms = read(solid_file)

#view(atoms)

count = 0
labels = atoms.get_chemical_symbols()

#Using Random Number Generator
#random.seed(1)

while new_species_total > count:
    i = int(math.floor(random.random()*len(atoms)))
    if labels[i] != new_species:
        labels[i] = new_species
        count = count + 1
        
#Double random
count_all = 0
new_labels = ['']*len(atoms)

while len(atoms) > count_all:
    j = int(math.floor(random.random()*len(atoms)))
    if new_labels[j] == '':
        new_labels[j] = labels[count_all]
        count_all += 1

#Selected every other atom in order of file
#for i in range(0,len(atoms),2):
#    labels[i] = new_species
#    count = count + 1

print count, " of species ", new_species, " added"
atoms.set_chemical_symbols(labels)

#view(atoms)

#print new_labels
atoms.set_chemical_symbols(new_labels)

#Reording positions so we have Pd/Au alternating from Z-axis up
#positions = atoms.get_positions()
#col = 2 # Z-axis
#positions = positions[numpy.argsort(positions[:,col])]
#col = 1 # Y-axis
#positions = positions[numpy.argsort(positions[:,col])]
#col = 0 # X-axis
#positions = positions[numpy.argsort(positions[:,col])]

#positions = numpy.sort(positions.view('f8,f8,f8'), order=['f0'], axis=0).view(numpy.float)
#positions = numpy.sort(positions.view('f8,f8,f8'), order=['f1'], axis=0).view(numpy.float)
#positions = numpy.sort(positions.view('f8,f8,f8'), order=['f2'], axis=0).view(numpy.float)

#new_positions = numpy.zeros((len(atoms),3),dtype=('f8'))

#count = 0
#while (len(atoms)-1) > count:
#    j = int(math.floor(random.random()*len(atoms)))
#    if new_positions[j][0] == 0.0 and new_positions[j][1] == 0.0 and new_positions[j][2] == 0.0:
#        new_positions[j] = positions[count]
#        count += 1

#print positions

#atoms.set_positions(positions)
#atoms.set_positions(new_positions)

view(atoms)

# Determine if a template is provided, if not use the one suggested
# Assume H/He/Li/Be labelling is consistent with original version
# Condense input file template into list containing just species labels if present

working_template = []

if len(regions_file) > 0:
    if os.path.exists(regions_file):
        with open(regions_file, 'r') as temp_template:
            template_lines = temp_template.readlines()  
            for item in template_lines[2:]:
                working_template.append(item.split()[0])
    else:
        print "Error: No template", template_option, "available"
        print "Available built-in templates are:", templates.keys()
        sys.exit(1)

    # Break up the list of charges into regions    
    core_species = []
    face_species = []
    edge_species = []
    vertex_species = []

    for i in range(len(new_labels)):
        if working_template[i] == 'H':
            core_species.append(new_labels[i])
        elif working_template[i] == 'He':
            face_species.append(new_labels[i])
        elif working_template[i] == 'Li':
            edge_species.append(new_labels[i])
        elif working_template[i] == 'Be':
            vertex_species.append(new_labels[i])
        else:
            print "A problem occurred during species breakdown."
            sys.exit(1)

    print "Vertices: ", len(vertex_species), "; Pd/Au: ", vertex_species.count('Pd'), "/", vertex_species.count('Au')
    print "Edges   : ", len(edge_species), "; Pd/Au: ", edge_species.count('Pd'), "/", edge_species.count('Au')
    print "Faces   : ", len(face_species), "; Pd/Au: ", face_species.count('Pd'), "/", face_species.count('Au')
    print "Core    : ", len(core_species), "; Pd/Au: ", core_species.count('Pd'), "/", core_species.count('Au')

