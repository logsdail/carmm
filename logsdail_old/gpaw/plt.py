import sys
from gpaw import GPAW
from gpaw import restart
from ase.io import write

# Requires Orthorhombic Unit Cells
filename = sys.argv[1]
atoms, calc = restart(filename)

density = calc.get_pseudo_density()
fname = filename + ".plt"
## Write to file
print "Writing ELF to file ", fname
write(fname, atoms, data=density)
