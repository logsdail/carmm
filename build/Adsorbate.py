from ase.io import read
from ase import Atoms
from ase.calculators.emt import EMT
from ase.calculators.aims import Aims
from ase.constraints import FixAtoms, FixBondLengths
from ase.optimize import BFGS
from ase.build import fcc111, add_adsorbate, surface, molecule, bulk
import numpy as np
from ase.visualize import view

########### ISAMBARD/FHI-AIMS SPECIFIC ###############
#import os
#command="aprun -n "+os.environ["NPROCS"]+" /home/ca-alogsdail/fhi-aims-gnu/bin/aims."+os.environ["VERSION"]+".scalapack.mpi.x"
#os.environ["ASE_AIMS_COMMAND"]=command
#os.environ["AIMS_SPECIES_DIR"]="/home/ca-alogsdail/fhi-aims-gnu/species_defaults/light"
######################################################

## Set separate calculators for each species:
## Adsorbate, surface and both
## Calculator in this case -  FHI-AIMS, other choices available - see ase.calculators
## Check configuration/bugtest with EMT calculator

def get_aims_calculator(n):
    #"gas" for gas-phase reactants and "periodic" for a periodic systems
    if(n=="gas"):
        return Aims(xc='pbe',
           spin='none',
           vdw_correction_hirshfeld="True",
           relativistic=('atomic_zora','scalar'),
           #use_dipole_correction='True',
           #default_initial_moment=2.0,
           compute_forces="true"
           #output=['mulliken']
           )
    else:
        if(n=="periodic"):
            return Aims(xc='pbe',
                spin='none',
                k_grid=(3,3,1),
                #vdw_correction_hirshfeld="True",
                relativistic=('atomic_zora','scalar'),
                use_dipole_correction='True',
                #default_initial_moment=2.0,
                compute_forces="true",
                #output=['mulliken']
                )

## Make a new directory for all generated files
#import os
#if not os.path.exists("data_folder"): os.mkdir("data_folder")
#os.chdir("data_folder")

######ADSORBATE#########
## Define adsorbate, rotate and pre-optimize
molecule=molecule("H2")
#molecule.rotate(90, 'x')
#molecule.set_calculator(get_aims_calculator("gas"))
molecule.set_calculator(EMT())

## Optimize
molecule_opt = BFGS(molecule, trajectory="adsorbate.traj", restart="adsorbate.pckl")
molecule_opt.run(fmax=0.01)

## You can also read the file geometry e.g. .traj file instead of calculating every time
#molecule = read("adsorbate.traj")

############SURFACE################################
from math import sqrt

## Edit these
atomic_species='Au'
shortest_M_M_distance=2.939
unit_cell_depth=3
unit_cell_width=3
slab_depth=4
vacuum_region_size=10.0

## Create surface
slab = fcc111(atomic_species, a=shortest_M_M_distance*sqrt(2), size=(unit_cell_width,unit_cell_depth,slab_depth))
#slab = fcc100(atomic_species, a=lattice_parameter*sqrt(2), size=(unit_cell_width,unit_cell_depth,slab_depth))

## Add vacuum
slab.center(vacuum=vacuum_region_size, axis=2)

## In short Pd surface
#slab2 = fcc111('Pd', a=3.914, size=(3,3,4), vacuum=20)

## For non-common(e.g. polymetallic) surfaces need to
## Define bulk and then cut it
#a = 3.78
#bulk = Atoms('Pd2Cu2',
#              scaled_positions=[(0, 0, 0),
#                                (0.5, 0.5, 0),
#                                (0.5, 0, 0.5),
#                                (0, 0.5, 0.5),
#                                ],
#              cell=[a, a, a],
#              pbc=True)
#slab3 = surface(bulk, (1, 1, 1), 4)
#slab3.center(vacuum=10, axis=2)

## Set constraint for surface/bulk characteristics
top_layers = 2      #leave x top layers relaxed
mask0 = [atom.tag > top_layers for atom in slab]
constraint0 = FixAtoms(mask=mask0)
slab.set_constraint([constraint0])

#slab.set_calculator(get_aims_calculator("periodic"))
slab.set_calculator(EMT())
#print(slab.get_potential_energy())

## Preview your generated surface prior to running calculations
## Don't use this on HAWK or ISAMBARD! Job will end prematurely.
#view(slab)
#view(molecule)

## Optimize the surface
surface_opt = BFGS(slab, trajectory="surface.traj", restart='surface.pckl')
surface_opt.run(fmax=0.01)

## If happy with your structures, remove triple quotes to proceed
"""
##########ADDING ADSORBATE ONTO THE SURFACE######
## Remove constraints and add adsorbate
slab.set_constraint()

## If build using ase.build can specify adsorption sites
add_adsorbate(slab, molecule, 2.0, 'ontop')
## Oherwise x-y coordinates - offset is specified in Angstrom
#add_adsorbate(slab, molecule, 2.0, position=(4.0, 2.4))

## Generate a new mask based on the changed number of atoms and constrain last two layers
mask0 = [atom.tag > top_layers for atom in slab]
constraint0 = FixAtoms(mask=mask0)
slab.set_constraint([constraint0])
#slab.set_calculator(get_aims_calculator("periodic"))
slab.set_calculator(EMT())

## Optimize
dyn = BFGS(slab, trajectory='ads_slab.traj', restart="ads_slab.pckl")
dyn.run(fmax=0.01)   #tighten to min 0.01 eV/A for actual calculation

## Extract Binding energy from optimised structures
optimised_molecule=read("adsorbate.traj")
optimised_surface=read("surface.traj")
optimised_slab=read("ads_slab.traj")
e_opt_molecule = optimised_molecule.get_potential_energy()
e_opt_surface = optimised_surface.get_potential_energy()
e_opt_slab = optimised_slab.get_potential_energy()
print("Energy of adsorbate: ", e_opt_molecule)
print("Energy of a clean surface: ", e_opt_surface)
print("Energy of an optimised slab: ", e_opt_slab)
print("Binding Energy: ", (e_opt_slab-(e_opt_molecule+e_opt_surface)))
"""
