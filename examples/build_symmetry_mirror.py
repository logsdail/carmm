#!/usr/bin/env python3

'''
TODO: Description Needed

'''

def get_example_model():
    '''Example model generation to show tool functionality'''
    from math import sqrt
    from ase.build import molecule, add_adsorbate, fcc111
    from ase.optimize import BFGS
    from ase.calculators.emt import EMT

    #### Traditional ASE functionality #####

    element='Au'
    lattice_parameter=2.939
    width=3
    depth=2
    vacuum=10.0

    # Create surface
    slab = fcc111(element, a=lattice_parameter*sqrt(2),
                  size=(width, width, depth), vacuum=vacuum)


    CO2 = molecule("CO2")
    CO2.rotate(90, 'x')
    CO2.rotate(50, 'z')

    add_adsorbate(slab, CO2, 3.0, "ontop")

    return slab

###### EXAMPLE OF USE ########
from ase.visualize import view
from math import sqrt
from software.analyse.neb_tools.symmetry import mirror, translation, rotate_fcc

# Toy model of CO2 on top of Au FCC(111)
model = get_example_model()
#view(model)

# Retrieve index of the C atom
index = [atom.index for atom in model if atom.symbol == "C"]
# Mirror lattice in the x plane with respect to C atom
# C atom remains in place, the rest of unit cell is shifted
model = mirror(model, center_index=index[0], plane='x', surf="111")
#view(model)


kwargs_x={"axis":0, "surface":"111"}
kwargs_y={"axis":1, "surface":"111"}

# Translate one row of atoms at a time in x and y
model = translation(model, **kwargs_x)
model = translation(model, **kwargs_y)
model = translation(model, **kwargs_x)
model = translation(model, **kwargs_y)
#view(model)

### ASSERTION ###
# Check if Oxygen moves as expected
eps = 1e-8
assert((model[19].position - [4.01974776, 2.0844678, 15.39968345] < [eps, eps, eps]).all())

model = rotate_fcc(model, center_index=index[0], surf=surf)
#view(model)
