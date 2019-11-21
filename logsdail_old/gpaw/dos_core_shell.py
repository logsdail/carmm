from ase import *
from gpaw import *
import gpaw.mpi as mpi
import sys
import os
import pickle

def write_dos_to_file(ef,energy,dos,output):
    if mpi.rank == 0:
        pickle.dump((ef, energy, dos), open(output, 'w'))

input_files = sys.argv[1:]
template = input_files[0]
cluster = input_files[1]

atoms, calc = restart(cluster+'.gpw')
system_energy = atoms.get_potential_energy()

# Determine if a template is provided, if not use the one suggested
# Assume H/He/Li/Be labelling is consistent with original version
# Condense input file template into list containing just species labels if present

working_template = []

if os.path.exists(template):
    with open(template, 'r') as temp_template:
        template_lines = temp_template.readlines()  
        for item in template_lines[2:]:
            working_template.append(item.split()[0])
else:
    print "Error: No template", template, "available"
    sys.exit(1)

# Break up the list of charges into regions    
core_charge = []
face_charge = []
edge_charge = []
vertex_charge = []

for i in range(len(atoms)):
    if working_template[i] == 'H':
        core_charge.append(i)
    elif working_template[i] == 'He':
        face_charge.append(i)
    elif working_template[i] == 'Li':
        edge_charge.append(i)
    elif working_template[i] == 'Be':
        vertex_charge.append(i)
    else:
        print "A problem occurred during charge breakdown."
        sys.exit(1)

# For now lets just collate all shell charges
shell_charge = face_charge + edge_charge + vertex_charge
        
# Get Fermi Level

try:
    ef = calc.get_fermi_level()
except ValueError:
    ef = 0

width = 0.2 
npts = 2001

# Get Dos for Core Atoms

if len(core_charge) > 0:

    total_dos = [0]*npts
            
    for i in range(len(core_charge)):
        print "Core: ", i, core_charge[i]
        lenergy, ldos = calc.get_orbital_ldos(a=core_charge[i], spin=0, angular='spdf', npts=npts, width=width)
        total_dos += ldos
        # Get Spin down as well if necessary and sum together
        if calc.get_number_of_spins() == 2:
            lenergy, ldos = calc.get_orbital_ldos(a=core_charge[i], spin=1, angular='spdf', npts=npts, width=width)
            total_dos += ldos
    
    write_dos_to_file(ef,lenergy,total_dos,'dos.core')

    print "Core done"

# Save Decomposed Dos of Core

    for angular in 'spdf':
        print "Orbital: ", angular
        total_dos = [0]*npts
        for i in range(len(core_charge)):
            print "Core: ", i, angular, core_charge[i]
            lenergy, ldos = calc.get_orbital_ldos(a=core_charge[i], angular=angular, npts=npts, width=width)
            total_dos += ldos
            # Get Spin down as well if necessary and sum together
            if calc.get_number_of_spins() == 2:
                lenergy, ldos = calc.get_orbital_ldos(a=core_charge[i], spin=1, angular=angular, npts=npts, width=width)
                total_dos += ldos

        write_dos_to_file(ef,lenergy,total_dos,'dos.core'+angular)

    print "Core decomposed done"

if len(shell_charge) > 0:

    total_dos = [0]*npts

    for i in range(len(shell_charge)):
        print "Shell: ", i, shell_charge[i]
        lenergy, ldos = calc.get_orbital_ldos(a=shell_charge[i], spin=0, angular='spdf', npts=npts, width=width)
        total_dos += ldos
        # Get Spin down as well if necessary and sum together
        if calc.get_number_of_spins() == 2:
            lenergy, ldos = calc.get_orbital_ldos(a=shell_charge[i], spin=1, angular='spdf', npts=npts, width=width)
            total_dos += ldos

    write_dos_to_file(ef,lenergy,total_dos,'dos.shell')

    print "Shell done"

# Save Decomposed Dos of Core

    for angular in 'spdf':
        print "Orbital: ", angular
        total_dos = [0]*npts
        for i in range(len(shell_charge)):
            print "Shell: ", i, angular, shell_charge[i]
            lenergy, ldos = calc.get_orbital_ldos(a=shell_charge[i], angular=angular, npts=npts, width=width)
            total_dos += ldos
            # Get Spin down as well if necessary and sum together
            if calc.get_number_of_spins() == 2:
                lenergy, ldos = calc.get_orbital_ldos(a=shell_charge[i], spin=1, angular=angular, npts=npts, width=width)
                total_dos += ldos

        write_dos_to_file(ef,lenergy,total_dos,'dos.shell'+angular)

    print "Shell decomposed done"
