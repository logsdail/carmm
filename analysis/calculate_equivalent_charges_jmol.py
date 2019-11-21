import sys
import os

input_files = sys.argv[1:]
charge_array = [[0 for x in xrange(0)] for x in xrange(len(input_files))]

if len(input_files) < 1:
    print "Error: Input file not found"
    sys.exit(1)

max_charge = 0
min_charge = 0

for i in range(len(input_files)):
    file = input_files[i]
    with open(file, 'r') as file_input:
        input_lines = file_input.readlines()
    for atom_data in input_lines[2:]:
        atom_data_list = atom_data.split()
        charge_array[i].append(float(atom_data_list[-1]))

    if min_charge > min(charge_array[i]):
        min_charge = min(charge_array[i])
    if max_charge < max(charge_array[i]):
        max_charge = max(charge_array[i])

print 'Overall minimum:', min_charge
print 'Overall_maximum:', max_charge
print 'Add 2 extra atoms with same coords as atom1 in each file with the above charges.'
