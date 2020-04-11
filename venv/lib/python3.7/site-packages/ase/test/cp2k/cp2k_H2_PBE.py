# -*- coding: utf-8 -*-

"""Test suit for the CP2K ASE calulator.

http://www.cp2k.org
Author: Ole Schuett <ole.schuett@mat.ethz.ch>
"""

from ase.build import molecule
from ase.calculators.cp2k import CP2K


def main():
    calc = CP2K(xc='PBE', label='test_H2_PBE')
    h2 = molecule('H2', calculator=calc)
    h2.center(vacuum=2.0)
    energy = h2.get_potential_energy()
    energy_ref = -31.5917284949
    diff = abs((energy - energy_ref) / energy_ref)
    assert diff < 1e-10
    print('passed test "H2_PBE"')


main()
