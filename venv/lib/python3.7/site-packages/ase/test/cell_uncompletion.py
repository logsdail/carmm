import itertools
from ase.cell import Cell


def all_pbcs():
    values = [False, True]
    yield from itertools.product(values, values, values)


def test(cell):
    for pbc in all_pbcs():
        ucell = cell.uncomplete(pbc)
        assert ucell.rank == sum(pbc & cell.any(1))

    assert all(cell.uncomplete(True).any(1) == cell.any(1)), (cell.uncomplete(True), cell)
    assert all(cell.uncomplete(1).any(1) == cell.any(1))
    assert cell.uncomplete(False).rank == 0
    assert cell.uncomplete(0).rank == 0

test(Cell.new([3, 4, 5]))
test(Cell.new([2, 0, 3]))
