import numpy as np
from ase.geometry import minkowski_reduce
from ase.cell import Cell


tol = 1E-14
rng = np.random.RandomState(0)
np.seterr(all='raise')


cell = Cell([[8.972058879514716, 0.0009788104586639142, 0.0005932485724084841],
             [4.485181755775297, 7.770520334862034, 0.00043663339838788054],
             [4.484671994095723, 2.5902066679984634, 16.25695615743613]])
cell.minkowski_reduce()

for i in range(40):
    B = rng.uniform(-1, 1, (3, 3))
    R, H = minkowski_reduce(B)
    assert np.allclose(H @ B, R, atol=tol)
    assert np.sign(np.linalg.det(B)) == np.sign(np.linalg.det(R))

    norms = np.linalg.norm(R, axis=1)
    assert (np.argsort(norms) == range(3)).all()

    # Test idempotency
    _, _H = minkowski_reduce(R)
    assert (_H == np.eye(3).astype(np.int)).all()

    rcell, _ = Cell(B).minkowski_reduce()
    assert np.allclose(rcell, R, atol=tol)


cell = np.array([[1, 1, 2], [0, 1, 4], [0, 0, 1]])
unimodular = np.array([[1, 2, 2], [0, 1, 2], [0, 0, 1]])
assert np.linalg.det(unimodular) == 1
lcell = unimodular.T @ cell

# test 3D
rcell, op = minkowski_reduce(lcell)
assert np.linalg.det(rcell) == 1

for pbc in [1, True, (1, 1, 1)]:
    rcell, op = minkowski_reduce(lcell, pbc=pbc)
    assert np.linalg.det(rcell) == 1
    assert np.sign(np.linalg.det(rcell)) == np.sign(np.linalg.det(lcell))

# test 0D
rcell, op = minkowski_reduce(lcell, pbc=[0, 0, 0])
assert (rcell == lcell).all()    # 0D reduction does nothing

# test 1D
for i in range(3):
    rcell, op = minkowski_reduce(lcell, pbc=np.roll([1, 0, 0], i))
    assert (rcell == lcell).all()    # 1D reduction does nothing

zcell = np.zeros((3, 3))
zcell[0] = lcell[0]
rcell, _ = Cell(zcell).minkowski_reduce()
assert np.allclose(rcell, zcell, atol=tol)

# test 2D
for i in range(3):
    pbc = np.roll([0, 1, 1], i)
    rcell, op = minkowski_reduce(lcell.astype(np.float), pbc=pbc)
    assert (rcell[i] == lcell[i]).all()

    zcell = np.copy(lcell.astype(np.float))
    zcell[i] = 0
    rzcell, _ = Cell(zcell).minkowski_reduce()
    rcell[i] = 0
    assert np.allclose(rzcell, rcell, atol=tol)
