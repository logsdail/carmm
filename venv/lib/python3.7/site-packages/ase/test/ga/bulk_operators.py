import os
import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.ga.utilities import closest_distances_generator, atoms_too_close
from ase.ga.bulk_utilities import CellBounds
from ase.ga.bulk_startgenerator import StartGenerator
from ase.ga.bulk_crossovers import CutAndSplicePairing
from ase.ga.bulk_mutations import (SoftMutation, RotationalMutation,
                                   RattleRotationalMutation, StrainMutation)
from ase.ga.ofp_comparator import OFPComparator
from ase.ga.offspring_creator import CombinationMutation
from ase.ga.standardmutations import RattleMutation, PermutationMutation

h2 = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.75]])
blocks = [('H', 4), ('H2O', 3), (h2, 2)]  # the building blocks
volume = 40. * sum([x[1] for x in blocks])  # cell volume in angstrom^3
splits = {(2,): 1, (1,): 1}  # cell splitting scheme

stoichiometry = []
for block, count in blocks:
    if type(block) == str:
        stoichiometry += list(Atoms(block).numbers) * count
    else:
        stoichiometry += list(block.numbers) * count

atom_numbers = list(set(stoichiometry))
blmin = closest_distances_generator(atom_numbers=atom_numbers,
                                    ratio_of_covalent_radii=1.3)

cellbounds = CellBounds(bounds={'phi': [30, 150], 'chi': [30, 150],
                                'psi': [30, 150], 'a': [3, 50],
                                'b': [3, 50], 'c': [3, 50]})

sg = StartGenerator(blocks, blmin, volume, cellbounds=cellbounds,
                    splits=splits)

# Generate 2 candidates
a1 = sg.get_new_candidate()
a1.info['confid'] = 1
a2 = sg.get_new_candidate()
a2.info['confid'] = 2

# Define and test genetic operators
pairing = CutAndSplicePairing(blmin, p1=1., p2=0., minfrac=0.15,
                              cellbounds=cellbounds, use_tags=True)

a3, desc = pairing.get_new_individual([a1, a2])
cell = a3.get_cell()
assert cellbounds.is_within_bounds(cell)
assert not atoms_too_close(a3, blmin, use_tags=True)

n_top = len(a1)
strainmut = StrainMutation(blmin, stddev=0.7, cellbounds=cellbounds,
                           use_tags=True)
softmut = SoftMutation(blmin, bounds=[2., 5.], used_modes_file=None,
                       use_tags=True)
rotmut = RotationalMutation(blmin, fraction=0.3, min_angle=0.5 * np.pi)
rattlemut = RattleMutation(blmin, n_top, rattle_prop=0.3, rattle_strength=0.5,
                           use_tags=True, test_dist_to_slab=False)
rattlerotmut = RattleRotationalMutation(rattlemut, rotmut)
permut = PermutationMutation(n_top, probability=0.33, test_dist_to_slab=False,
                             use_tags=True, blmin=blmin)
combmut = CombinationMutation(rattlemut, rotmut, verbose=True)
mutations = [strainmut, softmut, rotmut,
             rattlemut, rattlerotmut, permut, combmut]

for i, mut in enumerate(mutations):
    a = [a1, a2][i % 2]
    a3 = None
    while a3 is None:
        a3, desc = mut.get_new_individual([a])

    cell = a3.get_cell()
    assert cellbounds.is_within_bounds(cell)
    assert np.all(a3.numbers == a.numbers)
    assert not atoms_too_close(a3, blmin, use_tags=True)

modes_file = 'modes.txt'
softmut_with = SoftMutation(blmin, bounds=[2., 5.], use_tags=True,
                            used_modes_file=modes_file)
no_muts = 3
for _ in range(no_muts):
    softmut_with.get_new_individual([a1])
softmut_with.read_used_modes(modes_file)
assert len(list(softmut_with.used_modes.values())[0]) == no_muts
os.remove(modes_file)

comparator = OFPComparator(recalculate=True)
gold = bulk('Au') * (2, 2, 2)
assert comparator.looks_like(gold, gold)

# This move should not exceed the default threshold
gc = gold.copy()
gc[0].x += .1
assert comparator.looks_like(gold, gc)

# An additional step will exceed the threshold
gc[0].x += .2
assert not comparator.looks_like(gold, gc)
