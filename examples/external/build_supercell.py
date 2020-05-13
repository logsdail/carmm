#!/usr/bin/env python3

'''
Simple example to show functionality to create a supercell of any given model
'''

from ase.build import bulk

atoms = bulk('Pd', crystalstructure='fcc', a=3.909)
supercell = atoms*(3,3,3)

from ase.visualize import view
view(supercell)
