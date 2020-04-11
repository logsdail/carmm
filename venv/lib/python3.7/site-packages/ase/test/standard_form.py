import numpy as np
from numpy.testing import assert_allclose
from ase.cell import Cell


TOL = 1E-10
rng = np.random.RandomState(0)

for i in range(20):
    cell0 = rng.uniform(-1, 1, (3, 3))
    for sign in [-1, 1]:
        cell = Cell(sign * cell0)
        rcell, Q = cell.standard_form()
        assert_allclose(rcell @ Q, cell, atol=TOL)
        assert_allclose(np.linalg.det(rcell), np.linalg.det(cell))
        assert_allclose(rcell.ravel()[[1, 2, 5]], 0, atol=TOL)
