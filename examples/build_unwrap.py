#!/usr/bin/env python3

def test_build_unwrap():

    from carmm.build.unwrap import unwrap
    import numpy as np
    from ase.build import molecule

    a = molecule('CH4')
    orig_positions = a.positions.copy()
    a.cell = np.eye(3) * 10
    a.pbc = True

    a.wrap()
    unwrap(a)

    np.testing.assert_allclose(orig_positions, a.positions)

    
test_build_unwrap()
