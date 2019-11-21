import sys
import os

user_input = sys.argv[1:]
if not len(user_input) == 3:
    print "Error: One input file and one displacement value and one output file required as inputs"
    sys.exit(1)

input_file = user_input[0]
if not os.path.exists(input_file):
    print "Error: Input file not found"
    sys.exit(1)

input_value = user_input[1]
try:
    floated_value = float(input_value)
except ValueError:
    print("Error: second argument must be a numerical value for displacement")
    sys.exit(1)

output_file = user_input[2]

# User inputs are now input_file, output_file and floated_value

with open(input_file, 'r') as structure_input:
    structure_lines = structure_input.readlines()

new_file_data = []

for atom in structure_lines[2:]:
    atom_data = atom.split()
    label = atom_data[0]
    label_components = list(label)
    region_label = label_components[-1]
# To get proper behaviour for bulk clusters, need to had a condition here related to z-component
# e.g. any with z larger than a specified origin have alternative displacement rules
    if int(region_label) == 3:
        disp = floated_value
    elif int(region_label) == 2:
        disp = 1.5*floated_value
    elif int(region_label) == 1:
        disp = 2*floated_value
    else:
        disp = 0
# Remove numerical labelling - better for jmol (alternative: str(atom_data[0]) for labelled)
    jmol_label = ''.join(label_components[:-1])
    new_label = str(jmol_label)
    new_x = float(atom_data[1])
    new_y = float(atom_data[2])
    new_z = float(atom_data[3]) + disp
    new_atom_data = [new_label, new_x, new_y, new_z]
    if not int(region_label) == 5:
        new_file_data.append(new_atom_data)

#if not len(new_file_data) == int(structure_lines[0].split()[0]):
#    print "Error - number of atoms in old file header and new file don't match."
#    print "Header of old file:", len(new_file_data)
#    print "Number of atoms displaced:", structure_lines[0].split()[0]
#    sys.exit(1)
#else:  
#    print "Header of old file:", len(new_file_data)
#    print "Number of atoms displaced:", structure_lines[0].split()[0]

with open(output_file, 'w') as outfile:
    print >> outfile, len(new_file_data)
    print >> outfile, "Displaced xyz for visualization of ChemShell clusters"
    for output_data in new_file_data:
        print >> outfile, '  '.join([str(i) for i in output_data])


        



