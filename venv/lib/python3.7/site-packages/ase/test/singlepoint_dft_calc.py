import numpy as np
from ase.calculators.singlepoint import (SinglePointDFTCalculator,
                                         arrays_to_kpoints)
from ase.build import bulk


rng = np.random.RandomState(17)
nspins, nkpts, nbands = shape = 2, 4, 5
eps = 2 * rng.rand(*shape)
occ = rng.rand(*shape)
weights = rng.rand(nkpts)

kpts = arrays_to_kpoints(eps, occ, weights)



atoms = bulk('Au')

calc = SinglePointDFTCalculator(atoms)

calc.kpts = kpts

assert calc.get_number_of_spins() == nspins
assert calc.get_spin_polarized()
assert np.allclose(calc.get_k_point_weights(), weights)

for s in range(nspins):
    for k in range(nkpts):
        eps1 = calc.get_eigenvalues(kpt=k, spin=s)
        occ1 = calc.get_occupation_numbers(kpt=k, spin=s)
        assert np.allclose(eps1, eps[s, k])
        assert np.allclose(occ1, occ[s, k])


# XXX Should check more stuff.
