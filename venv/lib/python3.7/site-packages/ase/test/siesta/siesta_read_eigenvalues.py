import ase.build
from ase.calculators.siesta.siesta import Siesta

# Test real calculation of the lithium bulk which produced a gapped .EIG file
atoms = ase.build.bulk('Li', cubic=True)
calc = Siesta(kpts=[2,1,1])
atoms.set_calculator(calc)
atoms.get_potential_energy()

assert calc.results['eigenvalues'].shape[:2] == (1, 2)  # spins x bands
assert calc.get_k_point_weights().shape == (2,)
assert calc.get_ibz_k_points().shape == (2, 3)
