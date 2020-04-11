"""Bravais lattice type check.

1) For each Bravais variant, check that we recognize the
standard cell correctly.

2) For those Bravais lattices that we can recognize in non-standard form,
   Niggli-reduce them and recognize them as well."""
import numpy as np
from ase.lattice import (get_lattice_from_canonical_cell, all_variants,
                         identify_lattice)

for lat in all_variants():
    if lat.ndim == 2:
        break

    cell = lat.tocell()

    def check(lat1):
        print('check', repr(lat), '-->', repr(lat1))
        err = np.abs(cell.cellpar() - lat1.cellpar()).max()
        assert err < 1e-5, err

    check(get_lattice_from_canonical_cell(cell))

    if lat.name == 'TRI':
        # The TRI lattices generally permute (the ones we produce as
        # all_variants() are reduced to a form with smaller
        # orthogonality defect) which might be desirable but would
        # trigger an error in this test.
        continue

    stdcell, op = identify_lattice(cell, 1e-4)
    check(stdcell)
    rcell, op = cell.niggli_reduce()
    stdcell, op = identify_lattice(rcell, 1e-4)
    check(stdcell)
