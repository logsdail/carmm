#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 21:59:44 2020

@author: larakabalan

TODO: Description Needed

"""

from wulffpack import SingleCrystal
from ase.build import bulk
from ase.visualize import view

prim = bulk('Pd')

#### Create an Octahedron ####
surface_energies = {(1, 1, 1): 1,
                    (1, 0, 0): 10,}
particle = SingleCrystal(surface_energies,
                         primitive_structure=prim,
                         natoms=309)

#particle.view()
#view(particle.atoms)

#### Create a Truncated Octahedron ####
surface_energies = {(1, 1, 1): 1,
                      (1, 0, 0): 1}
                      #(1, 0, 0): 2 / 3**0.5}

particle = SingleCrystal(surface_energies,
                         primitive_structure=prim,
                         natoms=309)

#view(particle.atoms)

#### Create Decahedron ####
from wulffpack import Decahedron
particle = Decahedron(surface_energies,
                      twin_energy=0.03,
                      primitive_structure=prim,
                      natoms=309)

#view(particle.atoms)

#### Create an Icosahedron ####
from wulffpack import Icosahedron
particle = Icosahedron(surface_energies,
                       twin_energy=0.03,
                       primitive_structure=prim,
                       natoms=309)

#view(particle.atoms)

