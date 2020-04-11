"""Verify 2D Bravais lattice and band path versus pbc information."""

from ase.cell import Cell

cell = Cell([[1.,0.,0.],
             [.1,1.,0.],
             [0.,0.,0.]])
lat = cell.get_bravais_lattice()
print(cell.cellpar())
print(lat)
assert lat.name == 'OBL'

cell[2, 2] = 7
lat3d = cell.get_bravais_lattice()
print(lat3d)
assert lat3d.name == 'MCL'
lat2d_pbc = cell.get_bravais_lattice(pbc=[1, 1, 0])
print(lat2d_pbc)
assert lat2d_pbc.name == 'OBL'

path = cell.bandpath()
print(path)

path2d = cell.bandpath(pbc=[1, 1, 0])
print(path2d)
assert path2d.cell.rank == 2
assert path2d.cell.get_bravais_lattice().name == 'OBL'
