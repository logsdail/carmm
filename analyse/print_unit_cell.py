#!/usr/bin/python3

from ase.io import read
import sys

if len(sys.argv) < 2):
    print("Runtime arguments needed: ./script filename")
    exit()

for filename in sys.argv[1:]:
    cell = read(filename)
    print(filename, cell.get_cell_lengths_and_angles())

