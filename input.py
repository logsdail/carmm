from ase import Atoms
from ase.calculators.aims import Aims, AimsCube
from ase.optimize import QuasiNewton
import os
from ase.io import read
import numpy as np
from ase.io.read_vasp import (filename=POSCAR)
os.environ['ASE_AIMS_COMMAND']="mpirun -np "+os.environ['SLURM_NTASKS']+" /home/scw1057/software/fhi-aims/bin/aims."+os.environ['VERSION']+".scalapack.mpi.x"
os.environ['AIMS_SPECIES_DIR']="/home/scw1057/software/fhi-aims/species_defaults/light"   # Light settings
#os.environ['AIMS_SPECIES_DIR']="/home/c.sacal6/software/fhi-aims-species-defaults/tight" # Tight settings

atoms = read('POSCAR')


calc = Aims(xc='PBE',
            output=['dipole'],
            sc_accuracy_etot=1e-6,
            sc_accuracy_eev=1e-3,
            sc_accuracy_rho=1e-6,
            sc_accuracy_forces=1e-4,
            cubes=HBBC_cube)


HNBBC.set_calculator(calc)
dynamics = QuasiNewton(atoms, trajectory='HNBBC.traj')
dynamics.run(fmax=0.01)
