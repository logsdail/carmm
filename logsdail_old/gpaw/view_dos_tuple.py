#/usr/bin/python

import sys
import pickle
from array import array

def read_dos_from_file(filename):
    ef = 0
    energy = [] 
    dos = []

    try:
        ef, energy, dos = pickle.load(open(filename))
    except:
        print "A problem occurred reading the dos from ", filename
        sys.exit(1)        

    return ef, energy, dos

# Read input, make sure it is valid

inputs = sys.argv[1:]
filename = inputs[0]

ef = 0
energy = [] 
dos = []

ef, energy, dos = read_dos_from_file(filename)

print "Ef: ", ef
for i in range(len(energy)):
    print energy[i], ",", dos[i]
