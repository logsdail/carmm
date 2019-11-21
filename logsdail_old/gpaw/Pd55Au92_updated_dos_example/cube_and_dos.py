from ase.io import read, write, PickleTrajectory
from gpaw import GPAW, MixerSum, FermiDirac
from ase import *
from ase.units import Bohr
from gpaw.poisson import PoissonSolver
import sys
import os
import pickle
import gpaw.mpi as mpi

def write_dos_to_file(ef,energy,dos,output):
    if mpi.rank == 0:
        pickle.dump((ef, energy, dos), open(output, 'w'))

file="Cub-Pd2Au3"
width = 0.01
npts = 2001

# Set all GPAW options identical to those in the original calculation
# atoms, calc = restart(file+'.gpw')

#Read atoms
if os.path.exists(file+'.gpw'):
    atoms, calc = GPAW(file+'.gpw')
else:
    atoms=read(file+'.traj')

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

    mpi.world.barrier()

    # Attach a calculator
    atoms.set_calculator(calc)

system_energy = atoms.get_potential_energy()

# Get Fermi Level

try:
    ef = calc.get_fermi_level()
except ValueError:
    ef = 0

# Get Total Dos

energy, dos = calc.get_dos(spin=0, npts=npts, width=width)
write_dos_to_file(ef, energy, dos, 'dos.spin1')

# For Unpaired system get spin down also

if calc.get_number_of_spins() == 2:
    energy, dos = calc.get_dos(spin=1, npts=npts, width=width)
    write_dos_to_file(ef, energy, dos, 'dos.spin2')

# Get Dos for All Atoms. Disect species later.

for i in range(len(atoms)):
    total_dos_up = []
    total_dos_down = []

    lenergy, ldos = calc.get_orbital_ldos(a=i, spin=0, angular='spdf', npts=npts, width=width)
    write_dos_to_file(ef, lenergy, ldos, 'dos.'+str(i)+'.spin1')

    # Get Spin down as well if necessary and sum together
    if calc.get_number_of_spins() == 2:
        lenergy, ldos = calc.get_orbital_ldos(a=i, spin=0, angular='spdf', npts=npts, width=width)
        write_dos_to_file(ef, lenergy, ldos, 'dos.'+str(i)+'.spin2')

    for angular in 'spdf':
        lenergy, ldos = calc.get_orbital_ldos(a=i, angular=angular, npts=npts, width=width)
        write_dos_to_file(ef, lenergy, ldos, 'dos.'+str(i)+'.'+angular+'.spin1')

        # Get Spin down as well if necessary and sum together
        if calc.get_number_of_spins() == 2:
             lenergy, ldos = calc.get_orbital_ldos(a=i, spin=0, angular=angular, npts=npts, width=width)
             write_dos_to_file(ef, lenergy, ldos, 'dos.'+str(i)+'.'+angular+'.spin2')


# Write Cube File
rho = atoms.calc.get_all_electron_density(gridrefinement=2, broadcast=False)
if mpi.rank == 0:
    rho_scaled = rho * Bohr**3
    write('bader_analysis.cube', atoms, data=rho_scaled)

#if not os.path.exists(file+'.gpw'):
#    calc.write(file+'.gpw', mode ='all')
