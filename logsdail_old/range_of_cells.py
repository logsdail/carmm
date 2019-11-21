#!/usr/bin/python3

import ase.io.vasp
import sys
import numpy as np

if len(sys.argv) != 6:
    print 'Please give a POSCAR/CONTCAR, and lattice constant/angles min and max'
else:
    cell = ase.io.read(sys.argv[1])
    positions =   np.array([[0.0000000000000000,  0.0000000000000000,  0.0000000000000000],
  [0.0000000000000000,  0.5000000000000000,  0.5000000000000000],
  [0.5000000000000000,  0.0000000000000000,  0.5000000000000000],
  [0.5000000000000000,  0.5000000000000000,  0.0000000000000000],
  [0.0000000000000000,  0.0000000000000000,  0.5000000000000000],
  [0.0000000000000000,  0.5000000000000000,  0.0000000000000000],
  [0.5000000000000000,  0.5000000000000000,  0.5000000000000000],
  [0.5000000000000000,  0.0000000000000000,  0.0000000000000000],
  [0.2500000000000000,  0.2500000000000000,  0.2500000000000000],
  [0.2500000000000000,  0.2500000000000000,  0.7500000000000000],
  [0.2500000000000000,  0.7500000000000000,  0.2500000000000000],
  [0.2500000000000000,  0.7500000000000000,  0.7500000000000000],
  [0.7500000000000000,  0.2500000000000000,  0.2500000000000000],
  [0.7500000000000000,  0.2500000000000000,  0.7500000000000000],
  [0.7500000000000000,  0.7500000000000000,  0.2500000000000000],
  [0.7500000000000000,  0.7500000000000000,  0.7500000000000000]], dtype=float)
   
    x_min = sys.argv[2]
    x_max = sys.argv[3]
    alpha_min = sys.argv[4]
    alpha_max = sys.argv[5] 

    x_range = np.linspace(float(x_min),float(x_max),num=1+(float(x_max)-float(x_min))/0.01,endpoint=True)
    alpha_range = np.linspace(float(alpha_min),float(alpha_max),num=1+(float(alpha_max)-float(alpha_min))/0.1,endpoint=True)

    for x in x_range:
        initial_cell = np.array([[ 0, x/np.sqrt(2), x/np.sqrt(2) ],
                                 [ x/np.sqrt(2), 0, x/np.sqrt(2) ],
                                 [ x/np.sqrt(2), x/np.sqrt(2), 0 ]],
                                 dtype = float)
        #print initial_cell
        for alpha in alpha_range:
            
            #work out current angle
            a_dot_b = initial_cell[0][0]*initial_cell[1][0]+initial_cell[0][1]*initial_cell[1][1]+initial_cell[0][2]*initial_cell[1][2]
            #print a_dot_b, x
            cos_theta = a_dot_b / (x*x)
            #print cos_theta
            theta = np.arccos(cos_theta)*(180/np.pi)*3/2
            change = (theta-alpha)/22.325
            #print change

            new_sides = np.sqrt((x*x)-(change*change))/np.sqrt(2)

            new_cell = np.array([[ change, new_sides, new_sides ],
                                 [ new_sides, change, new_sides ],
                                 [ new_sides, new_sides, change ]],
                                 dtype = float)

            new_a_dot_b = new_cell[0][0]*new_cell[1][0]+new_cell[0][1]*new_cell[1][1]+new_cell[0][2]*new_cell[1][2]
            #print a_dot_b, x
            new_cos_theta = new_a_dot_b / (x*x)
            #print cos_theta
            new_theta = np.arccos(new_cos_theta)*(180/np.pi)*3/2
            #print theta, new_theta, x, np.sqrt(2*new_sides*new_sides+change*change)  
            cell.set_cell(new_cell)
            cell.set_scaled_positions(positions)
            #print cell.get_cell_lengths_and_angles()
            ase.io.vasp.write_vasp(sys.argv[1]+'.'+str(x)+'.'+"{0:.2f}".format(new_theta),cell,direct=True)

