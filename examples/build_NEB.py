#!/usr/bin/env python3

from software.build.neb import switch_indices, check_interpolation
from software.analyse.bonds import compare_structures
from ase.build import molecule

initial = molecule("CO2")
final = molecule("CO2")

#check_interpolation('initial.traj','final.traj',10)
#### Assertion tests ####
assert(check_interpolation(initial, final, 10, verbose=False, save=False))
########

atom_to_swap = 1
other_atom_to_swap = 2

#updated_final = switch_indices('final.traj', atom_to_swap, other_atom_to_swap):
updated_final = switch_indices(final, atom_to_swap, other_atom_to_swap)

#### Assertion tests ####
assert(not check_interpolation(initial, updated_final, 10, verbose=False, save=False))
########

#### Functionality to identify indices that would need swapping automatically
indices_to_swap, distances = compare_structures(initial, updated_final)
for i in range(len(indices_to_swap)):
    if indices_to_swap[i] is not i:
        reupdated_final = switch_indices(updated_final, i, indices_to_swap[i])

#### Assertion tests ####
assert(check_interpolation(initial, reupdated_final, 10, verbose=False, save=False))