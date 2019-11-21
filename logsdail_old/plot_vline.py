import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
import sys
import pylab

pylab.rcParams.update({'font.size': 32})
pylab.rcParams.update({'axes.linewidth': 2})

input_files = sys.argv[1:]
z_potential = []
x = []

if len(input_files) != 1:
    print "Error: Input file error"
    sys.exit(1)
else:
    file = input_files[0]

with open(file, 'r') as file_input:
    input_lines = file_input.readlines()
    del input_lines[0]
for word in input_lines:
    x.append(float(word.split()[0]))
    z_potential.append(float(word.split()[1]))

plt.plot(z_potential,x,lw ='3', color='red')

pylab.rcParams.update({'mathtext.fontset': 'stix'})

plt.xlabel('V(z)')
plt.ylabel(r'Z coordinate ($\AA$)')
plt.axis((-14,+10,x[0],x[-1]))
plt.gca().axes.get_yaxis().set_ticks([])
plt.tight_layout(pad=0.5)

plt.show()

