import ase.io.vasp
import sys
if len(sys.argv) != 2:
    print 'Please give just a POSCAR/CONTCAR'
else:
    cell = ase.io.read(sys.argv[1])
    ase.io.vasp.write_vasp(sys.argv[1]+'.DIRECT',cell,direct=True)
