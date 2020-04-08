#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 17:40:17 2020

@author: larakabalan
"""
# work in progress
from ase.build import fcc111
from ase.calculators.emt import EMT
from ase.data import atomic_numbers, reference_states
from ase.ga.data import PrepareDB
from ase.ga import set_raw_score
from ase.visualize import view
from ase.io import write
import random

# defining the new lattice parameters which is the average of the lattice parameters of Cu and Pt
def get_avg_lattice_constant(syms):
    a = 0.
    for m in set(syms):
        a += syms.count(m) * lattice_constants[m]
    return a / len(syms)


metals = ['Cu', 'Pd']
# Use experimental lattice constants
lattice_constants = dict((m, reference_states[atomic_numbers[m]]['a'])
                         for m in metals)

# Create the references (pure slabs) manually
pure_slabs = []
refs = {}
print('Reference energies:')
for m in metals:
    slab = fcc111(m, size=(2, 4, 3), a=lattice_constants[m],
                  vacuum=5, orthogonal=True)
    slab.set_calculator(EMT())

    # We save the reference energy as E_A / N
    e = slab.get_potential_energy()
    e_per_atom = e / len(slab)
    refs[m] = e_per_atom
    print('{0} = {1:.3f} eV/atom'.format(m, e_per_atom))

    # The mixing energy for the pure slab is 0 by definition
    set_raw_score(slab, 0.0)
    pure_slabs.append(slab)

# The population size should be at least the number of different compositions, slab lenght is x*y*z
pop_size = 2 * len(slab) #which is here 2*24=48

# We prepare the db and write a few constants that we are going to use later, hull.db should be removed or renamed everytime we are rerunning the script
db = PrepareDB('hull.db', population_size=pop_size,
               reference_energies=refs, metals=metals,
               lattice_constants=lattice_constants)

# We add the pure slabs to the database as relaxed because we have already
# set the raw_score
for slab in pure_slabs:
    db.add_relaxed_candidate(slab,
                             atoms_string=''.join(slab.get_chemical_symbols()))


# Now we create the rest of the candidates for the initial population, we are asking here Cu and Pd to be choosen randomley
for i in range(pop_size - 2):  # here we are asking for pos_size-2 of possibilities = 48-2 in this example
    # How many of each metal is picked at random, making sure that
    # we do not pick pure slabs, here for example we are asking that 2 of Cu to be picked randomley and spread with the other elements of Cu and Pd
    nA = random.randint(0, len(slab) - 2) # we are definnig nA to be from a value 0 to a maximium value of len(slab)-2
    nB = len(slab) - 2 - nA
    symbols = [metals[0]] * nA + [metals[1]] * nB + metals
    print ("chemical symbols = ", symbols)

    # Making a generic slab with the correct lattice constant
    slab = fcc111('X', size=(2, 4, 3),
                  a=get_avg_lattice_constant(symbols),
                  vacuum=5, orthogonal=True)

    # Setting the symbols and randomizing the order
    slab.set_chemical_symbols(symbols)
    random.shuffle(slab.numbers)

    # Add these candidates as unrelaxed, we will relax them later
    atoms_string = ''.join(slab.get_chemical_symbols())
    db.add_unrelaxed_candidate(slab, atoms_string=atoms_string)
    write("slab.xyz", slab)
    print(atoms_string)
    print("average lattice parameters = ", a)
    view(slab)    
