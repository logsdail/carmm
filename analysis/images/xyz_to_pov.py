import sys
from gpaw import GPAW
from ase.io import write, read

filename = sys.argv[1]
atoms = read(filename)

fname = filename[:-4] + ".pov"
## Write to file
print "Writing POV to file ", fname
write(fname, atoms)
