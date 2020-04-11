import os
from ase.build import bulk
from ase.calculators.nwchem import NWChem
from numpy.testing import assert_allclose


def main():
    atoms = bulk('C')

    testname = 'stress_test'

    calc = NWChem(theory='pspw',
                  label=testname,
                  nwpw={'lmbfgs': None,
                        'tolerances': '1e-9 1e-9'},
                  )
    atoms.set_calculator(calc)

    assert_allclose(atoms.get_stress(), calc.calculate_numerical_stress(atoms),
                    atol=1e-3, rtol=1e-3)

    # remove scratch files created by NWChem
    os.remove(os.path.join(testname, 'junk.inp'))
    os.remove(os.path.join(testname, testname + '.movecs'))
    os.remove(os.path.join(testname, testname + '.db'))
    for suffix in ['psp', 'vpp', 'vpp2']:
        os.remove(os.path.join(testname, 'C.' + suffix))


main()
