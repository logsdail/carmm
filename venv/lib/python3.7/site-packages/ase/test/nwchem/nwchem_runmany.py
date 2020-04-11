import os
from ase.build import molecule
from ase.calculators.nwchem import NWChem
from numpy.testing import assert_allclose


def _try_delete(theory, prefix, suffix, sep='.'):
    fname = os.path.join(theory, sep.join([prefix, suffix]))
    try:
        os.remove(fname)
    except FileNotFoundError:
        pass


def _run_calc(atoms_in, theory, eref, forces=False, **kwargs):
    atoms = atoms_in.copy()
    calc = NWChem(label=theory, theory=theory, **kwargs)
    atoms.set_calculator(calc)
    assert_allclose(atoms.get_potential_energy(), eref, atol=1e-4, rtol=1e-4)
    if forces:
        assert_allclose(atoms.get_forces(),
                        calc.calculate_numerical_forces(atoms),
                        atol=1e-4, rtol=1e-4)

    # Delete all perm/scratch files to ensure tests are idempotent
    for suffix in ['db', 'movecs', 'cfock', 'mp2nos', 't2']:
        _try_delete(theory, theory, suffix)

    for element in ['H', 'O']:
        for suffix in ['psp', 'vpp', 'cpp', 'jpp']:
            _try_delete(theory, element, suffix)
        _try_delete(theory, element, 'basis', sep='_')

    _try_delete(theory, 'junk', 'inp')


def main():
    atoms = molecule('H2O')
    # GTO calculations
    _run_calc(atoms, 'dft', -2051.9802410863354, basis='3-21G', forces=True)
    _run_calc(atoms, 'scf', -2056.7877421222634, basis='3-21G', forces=True)
    _run_calc(atoms, 'mp2', -2060.1413846247333, basis='3-21G', forces=True)
    _run_calc(atoms, 'ccsd', -2060.3418911515882, basis='3-21G')
    _run_calc(atoms, 'tce', -2060.319141863451, basis='3-21G', tce={'ccd': ''})

    # Plane wave calculations
    atoms.center(vacuum=2)
    atoms.pbc = True
    _run_calc(atoms, 'pspw', -465.1290581383751)
    _run_calc(atoms, 'band', -465.1290611316276)
    _run_calc(atoms, 'paw', -2065.6600649367365)


main()
