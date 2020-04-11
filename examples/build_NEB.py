#!/usr/bin/env python3

from software.build.neb import switch_indices, check_interpolation
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