#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 17:40:17 2020

@author: larakabalan
"""

from ase.build import fcc111
from ase.calculators.emt import EMT
from ase.data import atomic_numbers, reference_states
from ase.ga.data import PrepareDB
from ase.ga import set_raw_score
from ase.visualize import view
import random


def get_avg_lattice_constant(syms):
    a = 0.
    for m in set(syms):
        a += syms.count(m) * lattice_constants[m]
    return a / len(syms)


metals = ['Cu', 'Pt']
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

# The population size should be at least the number of different compositions
pop_size = 2 * len(slab)

# We prepare the db and write a few constants that we are going to use later
db = PrepareDB('hull.db', population_size=pop_size,
               reference_energies=refs, metals=metals,
               lattice_constants=lattice_constants)

# We add the pure slabs to the database as relaxed because we have already
# set the raw_score
for slab in pure_slabs:
    db.add_relaxed_candidate(slab,
                             atoms_string=''.join(slab.get_chemical_symbols()))


# Now we create the rest of the candidates for the initial population
for i in range(pop_size - 2):
    # How many of each metal is picked at random, making sure that
    # we do not pick pure slabs
    nA = random.randint(0, len(slab) - 2)
    nB = len(slab) - 2 - nA
    symbols = [metals[0]] * nA + [metals[1]] * nB + metals

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
    
    view(slab)    