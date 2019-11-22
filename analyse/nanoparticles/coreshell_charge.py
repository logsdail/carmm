#/usr/bin/python

import sys
import os

# Read 2 xyz files with charge data
# Assume the atoms are in the same order to avoid comparing coordinates
# Input format: cs_charge.py file1 file2
# Where file1 is a core/shell cluster and file2 is a solid cluster
input_files = sys.argv[1:]

coreshell_file = input_files[0]
solid_file = input_files[1]

# Exit if the specified file is not present
if not (os.path.exists(coreshell_file) and os.path.exists(solid_file)):
    print "Error: At least one input file not found"
    sys.exit(1)

with open(coreshell_file, 'r') as coreshell_input:
    coreshell_lines = coreshell_input.readlines()

with open(solid_file, 'r') as solid_input:
    solid_lines = solid_input.readlines()

# Check that the files have the same number of atoms
if not coreshell_lines[0] == solid_lines[0]:
    print "Error: Input files do not contain the same number of atoms"
    sys.exit(1)
else:
    file_length = int(coreshell_lines[0])

# 2 lines at the top of an xyz file not to be included
# Create an array specifyin wether the atoms match between the two files
match_array = [0]*file_length
count = 0
spare_lines = 2
for atom in range(spare_lines, file_length + spare_lines):
    if coreshell_lines[atom][0] == solid_lines[atom][0]:
        match_array[count] = 1
    else:
        match_array[count] = 0
    count += 1

# Divide the xyzq data into sets depending on if the atoms in the core/shell file
# are of the same element as those in the solid cluster file
coreshell_matching_atoms = []
coreshell_other_atoms = []
solid_matching_atoms = []
solid_other_atoms = []
match_array_count = 0
for item in match_array:
    if item == 0:
        coreshell_other_atoms.append(coreshell_lines[match_array_count + spare_lines])
        solid_other_atoms.append(solid_lines[match_array_count + spare_lines])
        match_array_count += 1
    elif item == 1:
        coreshell_matching_atoms.append(coreshell_lines[match_array_count + spare_lines])
        solid_matching_atoms.append(solid_lines[match_array_count + spare_lines])
        match_array_count += 1

# Sum the charges for each set of atoms
total_dq_coreshell_matching = 0.0
total_dq_coreshell_other = 0.0
total_dq_solid_matching = 0.0
total_dq_solid_other = 0.0

for data in coreshell_matching_atoms:
    line = data.split()
    total_dq_coreshell_matching += float(line[4])

for data in coreshell_other_atoms:
    line = data.split()
    total_dq_coreshell_other += float(line[4])

for data in solid_matching_atoms:
    line = data.split()
    total_dq_solid_matching += float(line[4])

for data in solid_other_atoms:
    line = data.split()
    total_dq_solid_other += float(line[4])

# Account for scaling of the charges for visualization
scaling = 20
scaled_dq_coreshell_matching = total_dq_coreshell_matching / scaling
scaled_dq_coreshell_other = total_dq_coreshell_other / scaling
scaled_dq_solid_matching = total_dq_solid_matching / scaling
scaled_dq_solid_other = total_dq_solid_other / scaling

total_matching =  len(coreshell_matching_atoms)
total_other = len(coreshell_other_atoms)

# Calculate per-atom dQ
peratom_dq_cs_matching = scaled_dq_coreshell_matching / total_matching
peratom_dq_cs_other = scaled_dq_coreshell_other / total_other
peratom_dq_s_matching = scaled_dq_solid_matching / total_matching
peratom_dq_s_other = scaled_dq_solid_other / total_other

# Divider
print ''

# Print outputs (wrt whole cluster}
print "Charge segregation - whole cluster:"
print coreshell_file, "- dQ on species", (coreshell_matching_atoms[0].split()[0]), scaled_dq_coreshell_matching

print coreshell_file, "- dQ on species", (coreshell_other_atoms[0].split()[0]), scaled_dq_coreshell_other

print "For the second species, with the same core/shell ordering..."

print solid_file, "- dQ on species", (solid_matching_atoms[0].split()[0]), scaled_dq_solid_matching

print solid_file, "- dQ on species", (solid_other_atoms[0].split()[0]), scaled_dq_solid_other

# Divider
print ''

# Print outputs (wrt single atoms}
print "Charge segregation - per atom:"
print coreshell_file, "- dQ on species", (coreshell_matching_atoms[0].split()[0]), peratom_dq_cs_matching

print coreshell_file, "- dQ on species", (coreshell_other_atoms[0].split()[0]), peratom_dq_cs_other

print "For the second species, with the same core/shell ordering..."

print solid_file, "- dQ on species", (solid_matching_atoms[0].split()[0]), peratom_dq_s_matching

print solid_file, "- dQ on species", (solid_other_atoms[0].split()[0]), peratom_dq_s_other

