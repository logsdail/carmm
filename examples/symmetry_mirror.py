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
    width_a=3
    width_b=3
    depth=2
    vacuum=10.0

    # Create surface
    slab = fcc111(element, a=lattice_parameter*sqrt(2), size=(width_a, width_b, depth))

    # Add vacuum
    slab.center(vacuum=vacuum, axis=2)
    slab.set_calculator(EMT())
    opt = BFGS(slab)
    opt.run(fmax=0.05)

    CO2 = molecule("CO2")
    CO2.rotate(90, 'x')
    CO2.rotate(50, 'z')

    add_adsorbate(slab, CO2, 3.0, "ontop")

    return slab

###### EXAMPLE OF USE ########
from ase.visualize import view
from math import sqrt
from software.analyse.neb_tools.symmetry import mirror, translation

# Toy model of CO2 on top of Au FCC(111)
model = get_example_model()
#view(model)

# Retrieve index of the C atom
index = [atom.index for atom in model if atom.symbol == "C"]
# Mirror lattice in the x plane with respect to C atom
# C atom remains in place, the rest of unit cell is shifted
model = mirror(model, center_index=index[0], plane='x', surf="111")
#view(model)

surf = "111"
lattice_parameter = 2.939/(sqrt(2)/2)
kwargs_x={"a":3.914, "axis":0, "surface":surf}
kwargs_y={"a":3.914, "axis":1, "surface":surf}

# Translate one row of atoms at a time in x and y
model = translation(model, **kwargs_x)
model = translation(model, **kwargs_y)
model = translation(model, **kwargs_x)
model = translation(model, **kwargs_y)
#view(model)

### ASSERTION ###
# Check if Oxygen moves as expected
eps = 1e-8
assert((model[19].position - [4.01974776, 2.0844678, 15.26547259] < [eps, eps, eps]).all())
