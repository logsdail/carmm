import ase.io.vasp
import sys
if len(sys.argv) != 11:
    print 'Please give a POSCAR/CONTCAR, and unit cell vectors x1, x2, x3, y1, y2, y3, z1, z2, z3'
else:
    cell = ase.io.read(sys.argv[1])
   
    x1 = sys.argv[2]
    x2 = sys.argv[3]
    x3 = sys.argv[4]
    y1 = sys.argv[5]
    y2 = sys.argv[6]
    y3 = sys.argv[7]
    z1 = sys.argv[8]
    z2 = sys.argv[9]
    z3 = sys.argv[10]

    cell.set_cell([[x1,x2,x3],
                   [y1,y2,y3],
                   [z1,z2,z3]])
    ase.io.vasp.write_vasp(sys.argv[1]+'.UNIT_CELL',cell)
