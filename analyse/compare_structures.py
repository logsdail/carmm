#!/usr/bin/python3

from ase.io import read
import sys

def print_non_zeros(square_distances, indices, filter_distance=0.0):
    from math import sqrt

    distances = []
    print("Filtering print-out for distances greater than ", sqrt(filter_distance), " \AA")
    for i in range(len(square_distances)):
        if square_distances[i] > 0.0:
            distances.append(sqrt(square_distances[i]))
            if square_distances[i] > filter_distance:
                print("Index structure 1: ", i, "; Index structure 2: ", indices[i], "; Difference: ", distances[-1])
    print("")
    print("RMSD (non-zero differences): ", sum(distances)/len(distances))
            

if len(sys.argv) < 3:
    print("Please provide at least two sturcture files")
    exit()

structure_one = read(sys.argv[1])
structure_two = read(sys.argv[2])
try:
    filter_cutoff = float(sys.argv[3])*float(sys.argv[3])
except:
    filter_cutoff = 0.01
     
if len(structure_one) != len(structure_two):
    print("The inputs don't contain the same number of atoms. How do you expect to compare them sensibly?")
    exit()

## Sadly, this doesn't work if the inputs are structured differently
#diff_positions = structure_two.positions - structure_one.positions
#
##print(diff_positions)
square_deviations = []
#for xyz in diff_positions:
#    square_deviations.append((xyz[0]*xyz[0])+(xyz[1]*xyz[1])+(xyz[2]*xyz[2]))

structure_two_indices = []
for i in range(len(structure_one.positions)):
    xyz = structure_one.positions[i]
    distance_square = 999999
    temp_index = 0
    for j in range(len(structure_two.positions)):
        temp_distance_square = ((structure_two.positions[j][0]-xyz[0])*(structure_two.positions[j][0]-xyz[0])
                               +(structure_two.positions[j][1]-xyz[1])*(structure_two.positions[j][1]-xyz[1])
                               +(structure_two.positions[j][2]-xyz[2])*(structure_two.positions[j][2]-xyz[2]))

        if distance_square > temp_distance_square and structure_one.symbols[i] == structure_two.symbols[j]:
            distance_square = temp_distance_square
            temp_index = j

    structure_two_indices.append(temp_index)
    square_deviations.append(distance_square)

print("Atoms that differ in position:")
print_non_zeros(square_deviations, structure_two_indices, filter_cutoff)
