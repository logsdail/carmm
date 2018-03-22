#!/usr/bin/python
import ase.io.vasp
import sys
if len(sys.argv) != 5:
    print 'Please give just a POSCAR, and expansion of the super cell in x, y and z'
else:
    cell = ase.io.read(sys.argv[1])
    ase.io.vasp.write_vasp(sys.argv[1]+'.SUPER',cell*(int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])),direct=True)
