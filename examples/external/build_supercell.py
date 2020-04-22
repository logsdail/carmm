#!/usr/bin/env python3

'''
Simple example to show functionality to create a supercell of any given model
'''

from ase.build import bulk
from software.build.facets import create_supercell

atoms = bulk('Pd', crystalstructure='fcc', a=3.909)
supercell = create_supercell(atoms, 3, 3, 3)

#### Assertion test ####
assert(len(supercell) == 27)
#########