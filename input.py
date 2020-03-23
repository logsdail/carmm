from ase import Atoms
from ase.calculators.aims import Aims, AimsCube
from ase.optimize import QuasiNewton
import os

os.environ['ASE_AIMS_COMMAND']="mpirun -np "+os.environ['SLURM_NTASKS']+" /home/scw1057/software/fhi-aims/bin/aims."+os.environ['VERSION']+".scalapack.mpi.x"
os.environ['AIMS_SPECIES_DIR']="/home/scw1057/software/fhi-aims/species_defaults/light"   # Light settings
#os.environ['AIMS_SPECIES_DIR']="/home/c.sacal6/software/fhi-aims-species-defaults/tight" # Tight settings

atom = Atoms('N', calculator=Aims(xc='pbe'))
e_atom = atom.get_potential_energy()

d=1.1
molecule = Atoms('2N', [(0., 0., 0.), (0., 0., d)])
molecule.set_calculator(Aims(xc='pbe'))
e_molecule = molecule.get_potential_energy()

e_atomization = e_molecule - 2* e_atom 

print('Nitrogen atom energy: %5.2f eV' % e_atom)
print('Nitrogen molecule energy: %5.2f eV' % e_molecule)
print('Atomization energy: %5.2f eV' % -e_atomization)



