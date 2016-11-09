from ase.io import read, write, PickleTrajectory
from gpaw import GPAW, Mixer, MixerSum, FermiDirac
from gpaw.poisson import PoissonSolver
from ase.optimize import QuasiNewton
from gpaw.mpi import world
import gpaw.mpi as mpi
import os

def isodd(num):
    return num & 1 and True or False

file="Cub-Pd2Au3"
mode="w"

#Read atoms
if os.path.exists(file+'.restart.gpw'):
        calc = GPAW(file+'.restart.gpw')
        atoms = read(file+'.traj')
        mode="a"
else:
    if os.path.exists(file+'.traj'):
        atoms=read(file+'.traj')
        mode="a"
    else:
        atoms = read(file+'.xyz')
        atoms.center(vacuum=6)
        atoms.set_pbc((False,False,False))

    #Assign array for magnetic moments
    magmoms=[0]*len(atoms)
    #Get centre of mass
    com=atoms.get_center_of_mass()

    #Variables
    min_dij=999999
    min_i=0
    odd_count=0
    spinpol=False
    mixer=Mixer(0.05,5)

    #Check number of Au atoms and distance from COM
    for i in range(len(atoms)):
        # Count number of Au
        if atoms.numbers[i] == 79:
            odd_count=odd_count+1
        # Distance from this to middle of cluster
            dij=(((atoms.positions[i][0]-com[0])**2)+((atoms.positions[i][1]-com[1])**2)+((atoms.positions[i][2]-com[2])**2))**(0.5)
            if dij < min_dij:
                min_i=i

    # Assign spin if we have an odd number of Au atoms
    if isodd(odd_count):
        magmoms[min_i] = 1
        # We only need to actually assign this for an initial calculation
        if mode == "w":
            atoms.set_initial_magnetic_moments(magmoms=magmoms)
        #Ensure spinpol calculation even if this is a restart
        spinpol=True
        mixer=MixerSum(0.05,5)

    #Assign calculator
    calc = GPAW(h=0.18,
                xc = 'PBE',
                spinpol=spinpol,
                charge=0.0,
                nbands=-80,
                convergence={'energy':0.0001,'density':1.0e-5,'eigenstates':1.0e-8},
                mixer=mixer,
                poissonsolver=PoissonSolver(eps=1e-12),
                eigensolver='rmm-diis',
                occupations=FermiDirac(0.1),
                txt=file+'.txt',
                maxiter=200,
                )

    world.barrier()

# Attach a calculator
atoms.set_calculator(calc)

world.barrier()

# Save wavefunctions for more efficient restarts
# calc.attach(calc.write, 400, file+'.restart.gpw')
# Set up geometry optimisation
traj = trajectory=PickleTrajectory(file+'.traj', mode=mode, atoms=atoms)
geom = QuasiNewton(atoms, trajectory=traj)
geom.run(fmax=0.01)

############# COPY AND PASTE OF CUBE_DOS, TO SAVE RE-RUNNING

system_energy = atoms.get_potential_energy()

# Get Fermi Level

try:
     ef = calc.get_fermi_level()
except ValueError:
     ef = 0

# Get Total Dos

width = 0.2
npts = 2001
energy, dos = calc.get_dos(spin=0, npts=npts, width=width)
output = "dos.spin1"
if mpi.rank == 0:
    pickle.dump((ef,energy, dos), open(output, 'w'))

# For Unpaired system get spin down also

if calc.get_number_of_spins() == 2:
     energy, dos = calc.get_dos(spin=1, npts=npts, width=width)
     output = "dos.spin2"
     if mpi.rank == 0:
         pickle.dump((ef,energy, dos), open(output, 'w'))

print "All done"

# Get Dos for Au Atom

gold = 0
palladium = 0
total = len(atoms)
indices = [0]*total

# Mark Au as 1, Pd as 0
for i in range(total):
     if atoms.numbers[i] == 79:
         indices[i] = 1
         gold += 1
     else:
         indices[i] = 0
         palladium +=1

print "Gold:    ", gold
print "Palladium:  ", palladium

# Gold atoms

if gold > 0:

     total_dos = [0]*npts

     for i in range(total):
         if indices[i] == 1:
             print "Gold: ", i, atoms.numbers[i]
             lenergy, ldos = calc.get_orbital_ldos(a=i, spin=0,
angular='spdf', npts=npts, width=width)
             total_dos += ldos
             # Get Spin down as well if necessary and sum together
             if calc.get_number_of_spins() == 2:
                 lenergy, ldos = calc.get_orbital_ldos(a=i, spin=0,
angular='spdf', npts=npts, width=width)
                 total_dos += ldos

     output = "dos.Au"
     if mpi.rank == 0:
         pickle.dump((ef,lenergy, total_dos), open(output, 'w'))

     print "Au done"

# Save Decomposed Dos of Au

     for angular in 'spdf':
         print "Orbital: ", angular
         total_dos = [0]*npts
         for i in range(total):
             if indices[i] == 1:
                 print "Gold: ", i, angular, atoms.numbers[i]
                 lenergy, ldos = calc.get_orbital_ldos(a=i,
angular=angular, npts=npts, width=width)
                 total_dos += ldos
                 # Get Spin down as well if necessary and sum together
                 if calc.get_number_of_spins() == 2:
                     lenergy, ldos = calc.get_orbital_ldos(a=i, spin=0,
angular=angular, npts=npts, width=width)
                     total_dos += ldos

             output = "dos.Au." + angular
             if mpi.rank == 0:
                 pickle.dump((ef, lenergy, total_dos), open(output, 'w'))

     print "Au decomposed done"

# Palladium atoms

if palladium > 0:

     total_dos = [0]*npts

     for i in range(total):
         if indices[i] == 0:
             print "Palladium: ", i, atoms.numbers[i]
             lenergy, ldos = calc.get_orbital_ldos(a=i, spin=0,
angular='spdf', npts=npts, width=width)
             total_dos += ldos
             # Get Spin down as well if necessary and sum together
             if calc.get_number_of_spins() == 2:
                 lenergy, ldos = calc.get_orbital_ldos(a=i, spin=0,
angular='spdf', npts=npts, width=width)
                 total_dos += ldos

     output = "dos.Pd"
     if mpi.rank == 0:
         pickle.dump((ef,lenergy, total_dos), open(output, 'w'))

     print "Pd done"

# Save Decomposed Dos of Pd

     for angular in 'spdf':
         print "Orbital: ", angular
         total_dos = [0]*npts
         for i in range(total):
             if indices[i] == 0:
                 print "Palladium: ", i, angular, atoms.numbers[i]
                 lenergy, ldos = calc.get_orbital_ldos(a=i,
angular=angular, npts=npts, width=width)
                 total_dos += ldos
                 # Get Spin down as well if necessary and sum together
                 if calc.get_number_of_spins() == 2:
                     lenergy, ldos = calc.get_orbital_ldos(a=i, spin=0,
angular=angular, npts=npts, width=width)
                     total_dos += ldos

             output = "dos.Pd." + angular
             if mpi.rank == 0:
                 pickle.dump((ef, lenergy, total_dos), open(output, 'w'))

     print "Pd decomposed done"

# Write Cube File
rho = atoms.calc.get_all_electron_density(gridrefinement=2, broadcast=False)
if mpi.rank == 0:
        rho_scaled = rho * Bohr**3
        write('bader_analysis.cube', atoms, data=rho_scaled)

#########################

template = 'labelled.xyz'

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

for i in range(total):
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
    
    output = "dos.core"
    if mpi.rank == 0:
        pickle.dump((ef,lenergy, total_dos), open(output, 'w'))

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

        output = "dos.core." + angular
        if mpi.rank == 0:
            pickle.dump((ef, lenergy, total_dos), open(output, 'w'))

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

    output = "dos.shell"
    if mpi.rank == 0:
        pickle.dump((ef,lenergy, total_dos), open(output, 'w'))

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

        output = "dos.shell." + angular
        if mpi.rank == 0:
            pickle.dump((ef, lenergy, total_dos), open(output, 'w'))

    print "Shell decomposed done"

############################################################
# Save final output
calc.write(file+'.gpw', mode ='all')

