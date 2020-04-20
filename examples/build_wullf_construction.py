#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 21:59:44 2020

@author: larakabalan

TODO: Description Needed

"""

from wulffpack import SingleCrystal
from wulffpack import Decahedron
from wulffpack import Icosahedron
from ase.build import bulk
from ase.io import write
prim = bulk('Pd', crystalstructure='fcc', a=3.909)
surface_energies = {(1, 1, 1): 1.1,
                    (1, 0, 1): 1.0,
                    (2, 1, 0): 1.0}
particle = SingleCrystal(surface_energies,
                         primitive_structure=prim,
                         natoms=5000)
particle.view()
write('fcc_single_crystal.xyz', particle.atoms)
