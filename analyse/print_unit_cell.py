#!/usr/bin/python3

import ase.io.vasp
import sys
from numpy import pi

max_distance = 3

cell = ase.io.read(sys.argv[1])

print(cell.get_cell_lengths_and_angles())

for i in range(len(cell)):
    if cell.numbers[i] != 8:
        for j in range(len(cell)):
            distance = cell.get_distance(i,j)
            if cell.numbers[i] != cell.numbers[j] and distance < max_distance:
                for k in range(len(cell)):
                    distance = cell.get_distance(j,k)
                    if cell.numbers[i] == cell.numbers[k] and distance < max_distance:
                        angle = cell.get_angle([i,j,k])*180. / pi
                        if angle > 80 and angle < 100:
                            print(i, j, k, cell.numbers[i], cell.numbers[j], cell.numbers[k], distance, angle)
