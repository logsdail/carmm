#!/usr/bin/env python3
'''
TODO: Description

'''

from ase.build import fcc111, add_adsorbate, molecule
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from software.run.aims_calculator import get_aims_calculator
from software.run.aims_path import set_aims_command

emt=True

######ADSORBATE#########
## Define adsorbate, rotate and pre-optimize
molecule=molecule("H2")
#molecule.rotate(90, 'x')
if emt:
    molecule.set_calculator(EMT())
else:
    set_aims_command()
    molecule.set_calculator(get_aims_calculator("gas"))

## Optimize
molecule_opt = BFGS(molecule) #, trajectory="adsorbate.traj", restart="adsorbate.pckl")
molecule_opt.run(fmax=0.01)
e_opt_molecule = molecule.get_potential_energy()

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
slab = fcc111(atomic_species, a=shortest_M_M_distance*sqrt(2),
              size=(unit_cell_width,unit_cell_depth,slab_depth),
              vacuum=vacuum_region_size)


## Set constraint for surface/bulk characteristics
top_layers = 2      #leave x top layers relaxed
mask0 = [atom.tag > top_layers for atom in slab]
constraint0 = FixAtoms(mask=mask0)
slab.set_constraint([constraint0])

if emt:
    slab.set_calculator(EMT())
else:
    slab.set_calculator(get_aims_calculator("periodic"))

## Preview your generated surface prior to running calculations
## Don't use this on HAWK or ISAMBARD! Job will end prematurely.
#from ase.visualize import view
#view(slab)
#view(molecule)

## Optimize the surface
surface_opt = BFGS(slab) #, trajectory="surface.traj", restart='surface.pckl')
surface_opt.run(fmax=0.01)
e_opt_surface = slab.get_potential_energy()

## If happy with your structures, remove triple quotes to proceed

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

if emt:
    slab.set_calculator(EMT())
else:
    slab.set_calculator(get_aims_calculator("periodic"))

## Optimize
dyn = BFGS(slab) #, trajectory='ads_slab.traj', restart="ads_slab.pckl")
dyn.run(fmax=0.01)   #tighten to min 0.01 eV/A for actual calculation
e_opt_slab = slab.get_potential_energy()

#print("Energy of adsorbate: ", e_opt_molecule)
#print("Energy of a clean surface: ", e_opt_surface)
#print("Energy of an optimised slab: ", e_opt_slab)
Eb = (e_opt_slab-(e_opt_molecule+e_opt_surface))
#print("Binding Energy: ", Eb)

#assert(abs(Eb - -0.173372) < 1e-5)