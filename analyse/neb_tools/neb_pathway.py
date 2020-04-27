'''
This is a work case for automatic NEB setup for a toy model
of CH3O + H --> CH3OH on Au(111)
TODO: Design a process where optimal initial/final adsorbate
    arrangements are established automatically. Predict outcomes
    of symmetry operations and form a sequence before attempting
    the operations to reduce computational time.
'''
from carmm.examples.data.model_gen import get_example_slab
from ase.build import molecule, add_adsorbate
from ase import Atom
from carmm.analyse.bonds import compare_structures
from carmm.build.neb import switch_indices, switch_all_indices

# Work case for automatic neb setup
# initial geometry
initial = get_example_slab()
ad1 = molecule("CH3O")
ad1.rotate(-90, 'x')
add_adsorbate(initial, ad1, height=2.0)
bridge = initial[12].position + ((initial[13].position - initial[12].position)/2)
ad2 = Atom("H", position=bridge)
initial += ad2

# final geometry
final = get_example_slab()
bridge2 = initial[13].position + ((initial[16].position - initial[13].position)/2)
ad3 = molecule("CH3OH")
ad3.rotate(60, 'x')
add_adsorbate(final, ad3, 3.0, position=(bridge2[0], bridge[1]))

# scramble indices to ensure compare_structures works
final = switch_indices(final, 18, 19)
final = switch_indices(final, 21, 18)

# CH3O + H --> CH3OH on Au surface
#from ase.visualize import view
#view(initial)
#view(final)

# start the process with index correction

final = switch_all_indices(final, compare_structures(final, initial)[0])

for i in [atom.index for atom in final]:
    if initial[i].symbol is not final[i].symbol:
        print("Chemical symbols different on atoms with the same indices")
