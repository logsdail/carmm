from ase import Atoms
from ase.calculators.lammpsrun import LAMMPS
from ase.lattice.cubic import BodyCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase import units
from ase.io import read

# Set executable names for LAMMPS
import os
os.environ['LAMMPS_COMMAND'] = "/apps/compilers/intel/2019.5/compilers_and_libraries_2019/linux/mpi/intel64/bin/mpirun -np 40 /home/scw1057/software/LAMMPS/src/lmp_mpi" # ASE 3.17
#os.environ['ASE_LAMMPSRUN_COMMAND'] = "mpirun -np 40 /home/c.sacal6/shared/software/LAMMPS/src/lmp_mpi" # ASE 3.18

# Variables for settings
#lattice_constant = 2.85 #angstrom #2.85102
#supercell = (5, 5, 5)
timestep = 2 #fs
start_temperature = 273 #K
oven_temperature = 1173 #K
temperature_steps = 100 #K

# Variables for MD calculation length
time_at_each_temperature = int(100000/timestep) #1ps as fs, 1ps at each temperature on ramp
time_between_outputs = int(1000/timestep) #fs

parameters = {'pair_style': 'meam/c',
              'pair_coeff': ['* * Jelinek_2012_meamf Al Si Mg Cu Fe Jelinek_2012_meam.alsimgcufe Fe']}

files = [ "Jelinek_2012_meamf", "Jelinek_2012_meam.alsimgcufe" ]

lammps = LAMMPS(#command="mpirun -np 4 lmp_openmpi", 
                parameters=parameters, files=files,
                keep_tmp_files=False, tmp_dir="tmp",
                always_triclinic=False, no_data_file=False )

fe = read(filename='Fe_0.15_doped.traj')

#defining fe
#fe = BodyCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
 #                         size=supercell, symbol='Fe', pbc = True, latticeconstant = lattice_constant)

#lammps calc and cell
fe.set_calculator(lammps)

# Only need for initialisation! So won't include in temperature ramp
# momenta and tot energy from ase MD tutorial
MaxwellBoltzmannDistribution(fe, start_temperature * units.kB)

from ase.md.npt import NPT
# NPT
ptime = 8.66
dyn = NPT(fe, timestep=timestep*units.fs, temperature=start_temperature*units.kB, externalstress=0.0000006324209, ttime=20*units.fs, pfactor= (ptime**2) *units.fs * 1.129469733, trajectory='mdFe.traj', logfile='mdFe.log')
#external stress = 1 atm converted into ev/A
#pfactor ptime x B(in eV/A). Converted 180.961 gPa from Fe energies graph

# Function to print energies and geometric observables.
def printenergy(a):
    """Function to print the potential, kinetic and total energy"""
    epot = a.get_potential_energy()
    ekin = a.get_kinetic_energy()
    length_angles = a.get_cell_lengths_and_angles()
    print('E per Unit Cell: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
              'Etot = %.3f Lengths A, B, C and Angles Alpha, Beta, Gamma = %.3f %.3f %.3f %.3f %.3f %.3f'
              % (epot, ekin, ekin / (1.5 * units.kB*len(a)), epot + ekin, length_angles[0], length_angles[1], length_angles[2], length_angles[3], length_angles[4], length_angles[5] ))

# dynamics run from ase MD tutorial
# Print initial conditions
printenergy(fe)

# heat the system up
for current_temperature in range(start_temperature, oven_temperature, temperature_steps):
    dyn.set_temperature(current_temperature*units.kB)

  
    # perform constant T and P dynamics for 100 ps
    for j in range(int(time_at_each_temperature/time_between_outputs)):
        dyn.run(time_between_outputs) # 1000fs
        printenergy(fe)

