import numpy as np
from ase.cell import Cell

# This tests a BCT cell which would be mischaracterized as MCLC
# depending on comparson's precision (fix: c432fd52ecfdca).
# The cell should actually be MCLC for small tolerances,
# and BCT with larger ones.  But it would always come out MCLC.
#
# The solution is that the Niggli reduction must run with a more
# coarse precision than the lattice recognition algorithm.
#
# Danger: Since the two mechanisms (Niggli, lattice recognition)
# define their precisions differently, it is not certain whether this
# problem is entirely gone.
cellpar = np.array([3.42864, 3.42864, 3.42864, 125.788, 125.788, 80.236])
cell = Cell.new(cellpar)
mclc = cell.get_bravais_lattice(eps=1e-4)
bct = cell.get_bravais_lattice(eps=1e-3)

print(mclc)
print(bct)
assert mclc.name == 'MCLC'
assert bct.name == 'BCT'

# Original cell is not perfect (rounding).
perfect_bct_cell = bct.tocell()
perfect_bct_cellpar = bct.cellpar()
assert perfect_bct_cell.get_bravais_lattice().name == 'BCT'
