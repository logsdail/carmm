from ase.io import read, write, PickleTrajectory
from gpaw import GPAW, Mixer, MixerSum, FermiDirac
from gpaw.poisson import PoissonSolver
from ase.optimize import QuasiNewton
from gpaw.mpi import world
import os

file="147-Cub-AuAg-Alloy"
mode="w"

#Read atoms
if os.path.exists(file+'.traj'):
    atoms=read(file+'.traj')
    mode="a"
else:
    atoms = read(file+'.xyz')
    atoms.center(vacuum=6)
    atoms.set_pbc((False,False,False))

# We only need to assign magnetic moments for an initial calculation
if mode == "w":
    magmoms=[0]*len(atoms)
    magmoms[-1]=1
    atoms.set_initial_magnetic_moments(magmoms=magmoms)

#Assign calculator
calc = GPAW(h=0.18,
            xc = 'PBE',
            spinpol=True,
            charge=0.0,
            nbands=-80,
            convergence={'energy':0.0001,'density':1.0e-5,'eigenstates':1.0e-8},
            mixer=MixerSum(0.05,5),
            poissonsolver=PoissonSolver(eps=1e-12),
            eigensolver='rmm-diis',
            occupations=FermiDirac(0.1),
            txt=file+'.txt',
            maxiter=200,
            )

world.barrier()

# Attach a calculator
atoms.set_calculator(calc)
# Set up geometry optimisation
traj = trajectory=PickleTrajectory(file+'.traj', mode=mode, atoms=atoms)
geom = QuasiNewton(atoms, trajectory=traj)
geom.run(fmax=0.01)
calc.write(file+'.gpw', mode ='all')
