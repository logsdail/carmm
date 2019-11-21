import ase.io.vasp
import sys
if len(sys.argv) != 2:
    print 'Please give just an xyz filename'
else:
    cell = ase.io.read(sys.argv[1])
    ase.io.vasp.write_vasp("POSCAR",cell,direct=False)
